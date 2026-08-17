#include "xpclr.hpp"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>

static int g_fails = 0;

static void expect_true(bool cond, const char* msg) {
    if (!cond) {
        std::cerr << "FAIL: " << msg << "\n";
        ++g_fails;
    }
}

static void expect_near(double a, double b, double eps, const char* msg) {
    if (!(std::isfinite(a) && std::isfinite(b) && std::fabs(a - b) <= eps)) {
        std::cerr << "FAIL: " << msg << " got=" << a << " want=" << b << "\n";
        ++g_fails;
    }
}

static xpclr::SnpData snp(int64_t pos, int x, int na, double q2) {
    xpclr::SnpData s;
    s.pos = pos;
    s.x_alt = x;
    s.n_a = na;
    s.q2 = q2;
    return s;
}

int main() {
    // determine_c
    expect_near(xpclr::determine_c(0.01, 0.0, 20000.0), 1.0, 0.0, "c(s=0)");
    expect_near(xpclr::determine_c(0.01, -1.0, 20000.0), 1.0, 0.0, "c(s<0)");
    {
        const double c = xpclr::determine_c(1e-8, 0.1, 20000.0);
        expect_true(c > 0.0 && c <= 1.0, "c in (0,1]");
        // rounded to 1e-5 grid
        expect_near(c, std::round(c * 1e5) / 1e5, 0.0, "c sf=5");
    }

    // estimate_omega raw mean
    std::vector<xpclr::SnpData> snps;
    // q1=0.5,q2=0.5 -> r=0
    snps.push_back(snp(100, 10, 20, 0.5));
    // q1=0.2,q2=0.5 -> ( -0.3)^2 / 0.25 = 0.36
    snps.push_back(snp(200, 4, 20, 0.5));
    // q1=0.8,q2=0.5 -> 0.36
    snps.push_back(snp(300, 16, 20, 0.5));
    // one huge outlier r: q1=1,q2=0.1 -> 0.9^2/(0.1*0.9)=9
    snps.push_back(snp(400, 20, 20, 0.1));

    const double raw = xpclr::estimate_omega(snps, 0.0);
    expect_near(raw, (0.0 + 0.36 + 0.36 + 9.0) / 4.0, 1e-12, "omega raw");

    // trim 0.25 drops 1 of 4 highest -> drop 9, mean of 0,0.36,0.36
    const double tr = xpclr::estimate_omega(snps, 0.25);
    expect_near(tr, (0.0 + 0.36 + 0.36) / 3.0, 1e-12, "omega trim0.25");

    // chen_likelihood finite at neutral-ish point
    {
        const double var = raw * 0.5 * 0.5;
        const double ll = xpclr::chen_likelihood(10, 20, 1.0, 0.5, var);
        expect_true(std::isfinite(ll), "chen_likelihood finite");
        expect_true(ll > -1800.0, "chen_likelihood not floor");
    }

    // determine_weights: all four LdMode paths on synthetic SnpSet.
    {
        using xpclr::LdMode;
        using xpclr::detail::determine_weights;
        auto mk_set = [](std::vector<std::vector<int8_t>> cols, int n_b) {
            xpclr::SnpSet s;
            s.n_b = n_b;
            for (size_t i = 0; i < cols.size(); ++i) {
                xpclr::SnpData d;
                d.pos = static_cast<int64_t>(i + 1);
                d.n_a = 20; d.n_b = 2 * n_b; d.q2 = 0.5;
                s.snps.push_back(d);
                s.dosage_b.insert(s.dosage_b.end(), cols[i].begin(), cols[i].end());
            }
            return s;
        };
        std::vector<int> all_ix = {0, 1, 2};

        // dosage_fill (mode 0): rows 0121 / 0121 / 2101, missing filled as-is
        {
            auto s = mk_set({{0,1,2,1}, {0,1,2,1}, {2,1,0,1}}, 4);
            auto w = determine_weights(s, all_ix, 0.95, LdMode::dosage_fill);
            // r01=+1, r02=-1, r12=-1; r^2=1 > 0.95 for all pairs -> each SNP
            // above with 2 partners -> w = 1/3
            expect_near(w[0], 1.0/3.0, 1e-9, "weights dosage_fill all-linked");
            expect_near(w[1], 1.0/3.0, 1e-9, "weights dosage_fill s1");
            expect_near(w[2], 1.0/3.0, 1e-9, "weights dosage_fill s2");
        }
        // uncorrelated: rows 0121 / 1030-ish -> r^2 low -> w=1
        {
            auto s = mk_set({{0,1,2,1}, {1,1,1,1}}, 4);
            auto w = determine_weights(s, {0, 1}, 0.95, LdMode::dosage_fill);
            // snp1 constant -> ss=0 -> r NaN -> counted above by isnan rule
            expect_near(w[0], 0.5, 1e-9, "weights NaN-partner counted");
            expect_near(w[1], 0.5, 1e-9, "weights NaN-partner counted s1");
        }
        // pairwise (mode 3): per-pair complete samples only
        {
            // snp0 = 0 1 -9, snp1 = 0 -9 1: only k=0 shared valid -> cnt<2 -> NaN
            auto s = mk_set({{0,1,-9}, {0,-9,1}}, 3);
            auto w = determine_weights(s, {0, 1}, 0.95, LdMode::pairwise);
            expect_near(w[0], 0.5, 1e-9, "weights pairwise too-few-shared -> NaN");
            expect_near(w[1], 0.5, 1e-9, "weights pairwise s1");
        }
        // pairwise with enough shared samples and perfect correlation
        {
            auto s = mk_set({{0,1,2,-9}, {0,1,2,1}}, 4);
            auto w = determine_weights(s, {0, 1}, 0.95, LdMode::pairwise);
            // shared k=0..2: r=1 over 3 samples -> above -> w=1/2
            expect_near(w[0], 0.5, 1e-9, "weights pairwise perfect r");
        }
        // phased (mode 1): haplotype rows 0/1; identical rows -> r=1
        {
            auto s = mk_set({{0,1,1,0}, {0,1,1,0}}, 4);
            auto w = determine_weights(s, {0, 1}, 0.95, LdMode::phased);
            expect_near(w[0], 0.5, 1e-9, "weights phased identical haplotypes");
        }
        // em (mode 2): dosages; two samples snp0=0/1 snp1=0/1
        // haplotypes: ab, AB -> pAB=0.5, D=0.25, r=1
        {
            auto s = mk_set({{0,1}, {0,1}}, 2);
            auto w = determine_weights(s, {0, 1}, 0.95, LdMode::em);
            expect_near(w[0], 0.5, 1e-6, "weights EM perfect LD");
        }
        // em: uncorrelated double-het only -> D=0 -> r=0 -> w=1
        {
            auto s = mk_set({{0,1,1,0}, {1,0,0,1}}, 4);
            auto w = determine_weights(s, {0, 1}, 0.95, LdMode::em);
            // all four samples double het; EM on symmetric counts -> D=0, r=0
            expect_near(w[0], 1.0, 1e-6, "weights EM repulsion symmetric");
        }
    }

    if (g_fails) {
        std::cerr << g_fails << " failure(s)\n";
        return 1;
    }
    std::cout << "test_math ok\n";
    return 0;
}
