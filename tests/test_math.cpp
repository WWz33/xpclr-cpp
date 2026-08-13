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

    if (g_fails) {
        std::cerr << g_fails << " failure(s)\n";
        return 1;
    }
    std::cout << "test_math ok\n";
    return 0;
}
