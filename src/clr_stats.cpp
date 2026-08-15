#include "xpclr.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <mutex>
#include <vector>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_gamma.h>

// glibc hides M_PI under strict -std=c++17 (__STRICT_ANSI__).
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace xpclr {

// ---- math helpers (Python methods.py) ----
// Values fixed for hardingnj/python parity — do not change.
namespace {
constexpr double kQuadLo = 0.001;
constexpr double kQuadHi = 0.999;
constexpr double kQuadEpsRel = 0.001;
constexpr size_t kQuadLimit = 1000;
constexpr double kLikeFloor = -1800.0;  // when integral vanishes
}  // namespace

double determine_c(double r, double s, double ne, double min_rd, int sf) {
    if (s <= 0.0) return 1.0;
    double c = 1.0 - std::exp(-std::log(2.0 * ne) * std::max(r, min_rd) / s);
    // default sf=5 → scale 1e5 (match hardingnj round)
    const double scale = (sf == 5) ? 1e5 : std::pow(10.0, sf);
    return std::round(c * scale) / scale;
}

// methods.pdf for a scalar p1.
// IMPORTANT: when c > 0.5, for p1 in (1-c, c) BOTH left and right terms are added
// (Python vectorized path uses independent bisect slices, not if/else).
static double pdf_scalar(double p1, double c, double p2, double var) {
    if (!(var > 0.0) || !(c > 0.0) || !std::isfinite(p1)) return 0.0;
    double a_term = 1.0 / std::sqrt(2.0 * M_PI * var);
    double r = 0.0;
    if (p1 < c) {
        double b = (c - p1) / (c * c);
        double ct = (p1 - c * p2) * (p1 - c * p2) / (2.0 * c * c * var);
        r += a_term * b * std::exp(-ct);
    }
    if (p1 > 1.0 - c) {
        double b = (p1 + c - 1.0) / (c * c);
        double ct =
            (p1 + c - 1.0 - c * p2) * (p1 + c - 1.0 - c * p2) / (2.0 * c * c * var);
        r += a_term * b * std::exp(-ct);
    }
    return r;
}

static double binom_pmf_pre(int x, int n, double p, double logc) {
    if (x < 0 || x > n) return 0.0;
    if (p <= 0.0) return (x == 0) ? 1.0 : 0.0;
    if (p >= 1.0) return (x == n) ? 1.0 : 0.0;
    double lp = logc + x * std::log(p) + (n - x) * std::log1p(-p);
    if (!std::isfinite(lp)) return 0.0;
    return std::exp(lp);
}

struct QuadParams {
    int xj = 0;
    int nj = 0;
    double c = 0;
    double p2 = 0;
    double var = 0;
    double logc = 0;  // gsl_sf_lnchoose(nj, xj), once per likelihood
    bool with_binom = false;
};

static double integrand_gsl(double p1, void* p) {
    auto* st = static_cast<QuadParams*>(p);
    double dens = pdf_scalar(p1, st->c, st->p2, st->var);
    if (!st->with_binom) return dens;
    return dens * binom_pmf_pre(st->xj, st->nj, p1, st->logc);
}

// scipy.integrate.quad(..., epsrel=0.001, epsabs=0)
// Suppress GSL abort on roundoff/singularity; return best estimate / 0.
static void gsl_error_off(const char*, const char*, int, int) {}

static double gsl_quad(QuadParams* st, double a, double b) {
    // The error handler is process-wide; set it exactly once across all threads.
    static std::once_flag gsl_handler_once;
    std::call_once(gsl_handler_once,
                   [] { gsl_set_error_handler(&gsl_error_off); });
    gsl_function F;
    F.function = &integrand_gsl;
    F.params = st;
    // thread-local workspace (OpenMP safe)
    static thread_local gsl_integration_workspace* ws = nullptr;
    if (!ws) ws = gsl_integration_workspace_alloc(kQuadLimit);
    double result = 0.0, abserr = 0.0;
    // Primary: match scipy.integrate.quad(epsrel=0.001, epsabs=0).
    int status = gsl_integration_qags(&F, a, b, 0.0, kQuadEpsRel, kQuadLimit, ws,
                                      &result, &abserr);
    // Fallback only when integrator fails or returns non-finite/negative mass.
    // Do not coarsen a successful primary integral (preserves normal results).
    if (status != GSL_SUCCESS || !std::isfinite(result) || result < 0.0) {
        double fb = 0.0, fb_err = 0.0;
        const int st2 = gsl_integration_qags(&F, a, b, 1e-8, 1e-2, kQuadLimit, ws, &fb,
                                             &fb_err);
        if (st2 == GSL_SUCCESS && std::isfinite(fb) && fb >= 0.0) {
            result = fb;
        } else {
            result = 0.0;  // chen_likelihood maps 0 -> kLikeFloor
        }
    }
    return result;
}

double chen_likelihood(int xj, int nj, double c, double p2, double var) {
    QuadParams st;
    st.xj = xj;
    st.nj = nj;
    st.c = c;
    st.p2 = p2;
    st.var = var;
    st.logc = (xj >= 0 && xj <= nj) ? gsl_sf_lnchoose(nj, xj) : 0.0;
    st.with_binom = true;
    double like_i = gsl_quad(&st, kQuadLo, kQuadHi);
    st.with_binom = false;
    double like_b = gsl_quad(&st, kQuadLo, kQuadHi);
    if (like_i == 0.0 || like_b == 0.0) return kLikeFloor;
    return std::log(like_i) - std::log(like_b);
}

double estimate_omega(const SnpSet& snps, double trim) {
    return estimate_omega(snps.snps, trim);
}

double estimate_omega(const std::vector<SnpData>& snps, double trim) {
    if (!(trim >= 0.0 && trim < 1.0) || !std::isfinite(trim))
        die("estimate_omega: trim must be in [0,1)");
    std::vector<double> r;
    r.reserve(snps.size());
    for (const auto& s : snps) {
        if (!(s.q2 > 0.0 && s.q2 < 1.0) || s.n_a <= 0) continue;
        const double q1 = static_cast<double>(s.x_alt) / static_cast<double>(s.n_a);
        const double den = s.q2 * (1.0 - s.q2);
        if (!(den > 0.0) || !std::isfinite(den)) continue;
        const double ri = (q1 - s.q2) * (q1 - s.q2) / den;
        if (!std::isfinite(ri)) continue;
        r.push_back(ri);
    }
    if (r.empty()) die("estimate_omega: no valid SNPs");
    if (trim == 0.0 || r.size() < 2) {
        double sum = 0.0;
        for (double v : r) sum += v;
        return sum / static_cast<double>(r.size());
    }
    // Drop highest floor(n*trim) values; keep at least 1.
    std::sort(r.begin(), r.end());
    size_t drop = static_cast<size_t>(std::floor(static_cast<double>(r.size()) * trim));
    if (drop >= r.size()) drop = r.size() - 1;
    const size_t keep = r.size() - drop;
    double sum = 0.0;
    for (size_t i = 0; i < keep; ++i) sum += r[i];
    return sum / static_cast<double>(keep);
}

// Pearson corr pairwise (Rogers-Huff r on dosage), condensed upper triangle.
// If skip_missing, values < 0 (sentinel -9) are excluded per SNP pair.
static std::vector<float> rogers_huff_r(const std::vector<const int8_t*>& rows,
                                        int nvar, int nsamp, bool skip_missing) {
    if (skip_missing) {
        // Pairwise-complete: recompute mean/ss/cov per pair using only valid samples.
        std::vector<float> out;
        out.reserve(static_cast<size_t>(nvar) * (nvar - 1) / 2);
        for (int i = 0; i < nvar; ++i) {
            for (int j = i + 1; j < nvar; ++j) {
                double si = 0, sj = 0, si2 = 0, sj2 = 0, cov = 0;
                int cnt = 0;
                for (int k = 0; k < nsamp; ++k) {
                    int8_t vi = rows[i][k], vj = rows[j][k];
                    if (vi < 0 || vj < 0) continue;
                    si += vi; sj += vj;
                    si2 += (double)vi * vi;
                    sj2 += (double)vj * vj;
                    cov += (double)vi * vj;
                    ++cnt;
                }
                if (cnt < 2) { out.push_back(std::numeric_limits<float>::quiet_NaN()); continue; }
                double mi = si / cnt, mj = sj / cnt;
                double vi = si2 / cnt - mi * mi;
                double vj = sj2 / cnt - mj * mj;
                if (vi <= 0.0 || vj <= 0.0) { out.push_back(std::numeric_limits<float>::quiet_NaN()); continue; }
                double cv = cov / cnt - mi * mj;
                out.push_back(static_cast<float>(cv / std::sqrt(vi * vj)));
            }
        }
        return out;
    }
    // Original: use all samples (missing was filled 0).
    std::vector<double> mean(nvar, 0.0), ss(nvar, 0.0);
    for (int i = 0; i < nvar; ++i) {
        double s = 0.0;
        for (int j = 0; j < nsamp; ++j) s += rows[i][j];
        mean[i] = s / nsamp;
        for (int j = 0; j < nsamp; ++j) {
            double d = rows[i][j] - mean[i];
            ss[i] += d * d;
        }
    }
    std::vector<float> out;
    out.reserve(static_cast<size_t>(nvar) * (nvar - 1) / 2);
    for (int i = 0; i < nvar; ++i) {
        for (int j = i + 1; j < nvar; ++j) {
            if (ss[i] == 0.0 || ss[j] == 0.0) {
                out.push_back(std::numeric_limits<float>::quiet_NaN());
                continue;
            }
            double cov = 0.0;
            for (int k = 0; k < nsamp; ++k) {
                cov += (rows[i][k] - mean[i]) * (rows[j][k] - mean[j]);
            }
            out.push_back(static_cast<float>(cov / std::sqrt(ss[i] * ss[j])));
        }
    }
    return out;
}

// Haplotype-based Pearson r (raw -p1): rows are haplotypes (0/1, -9 missing).
// r = |pAB - pA*pB| / sqrt(pA(1-pA) pB(1-pB))
static std::vector<float> haplotype_r(const std::vector<const int8_t*>& rows,
                                      int nvar, int nsamp) {
    std::vector<float> out;
    out.reserve(static_cast<size_t>(nvar) * (nvar - 1) / 2);
    for (int i = 0; i < nvar; ++i) {
        for (int j = i + 1; j < nvar; ++j) {
            int cA = 0, cB = 0, cAB = 0, cnt = 0;
            for (int k = 0; k < nsamp; ++k) {
                int8_t vi = rows[i][k], vj = rows[j][k];
                if (vi < 0 || vj < 0) continue;
                cA += vi;
                cB += vj;
                cAB += vi * vj;
                ++cnt;
            }
            if (cnt == 0) { out.push_back(std::numeric_limits<float>::quiet_NaN()); continue; }
            double pA = (double)cA / cnt, pB = (double)cB / cnt, pAB = (double)cAB / cnt;
            double da = pA * (1.0 - pA), db = pB * (1.0 - pB);
            if (da <= 0.0 || db <= 0.0) { out.push_back(std::numeric_limits<float>::quiet_NaN()); continue; }
            double r = (pAB - pA * pB) / std::sqrt(da * db);
            out.push_back(static_cast<float>(std::fabs(r)));
        }
    }
    return out;
}

// EM two-locus phase inference (raw -p0): rows are dosages (0/1/2, -9 missing).
// Infers haplotype freq pAB from genotype counts via EM, then Pearson r.
static std::vector<float> em_haplotype_r(const std::vector<const int8_t*>& rows,
                                          int nvar, int nsamp) {
    std::vector<float> out;
    out.reserve(static_cast<size_t>(nvar) * (nvar - 1) / 2);
    for (int i = 0; i < nvar; ++i) {
        for (int j = i + 1; j < nvar; ++j) {
            // Count diploid genotypes: 0/1/2 at locus i and j.
            // Genotype categories for two biallelic loci:
            //   double heterozygote (1,1): phase unknown -> EM
            //   others: phase known
            int n00=0, n01=0, n02=0, n10=0, n11=0, n12=0, n20=0, n21=0, n22=0;
            int ndip = 0;
            for (int k = 0; k < nsamp; ++k) {
                int8_t di = rows[i][k], dj = rows[j][k];
                if (di < 0 || dj < 0) continue;
                int gi = di, gj = dj;  // 0/1/2
                ++ndip;
                if      (gi == 0 && gj == 0) ++n00;
                else if (gi == 0 && gj == 1) ++n01;
                else if (gi == 0 && gj == 2) ++n02;
                else if (gi == 1 && gj == 0) ++n10;
                else if (gi == 1 && gj == 1) ++n11;
                else if (gi == 1 && gj == 2) ++n12;
                else if (gi == 2 && gj == 0) ++n20;
                else if (gi == 2 && gj == 1) ++n21;
                else if (gi == 2 && gj == 2) ++n22;
            }
            if (ndip == 0) { out.push_back(std::numeric_limits<float>::quiet_NaN()); continue; }
            // Haplotype freqs: pAB, pAb, paB, pab (A=alt at i, B=alt at j)
            // Genotype -> haplotype contributions (diploid, 2 haplotypes each):
            //   (0,0) aa/bb -> 2 pab        (0,1) aa/AB -> 1 paB + 1 pab
            //   (0,2) aa/BB -> 2 paB        (1,0) AB/bb -> 1 pAb + 1 pab
            //   (1,2) AB/BB -> 1 pAB + 1 paB  (2,0) AA/bb -> 2 pAb
            //   (2,1) AA/AB -> 1 pAB + 1 pAb  (2,2) AA/BB -> 2 pAB
            //   (1,1) AB/AB -> double het, phase unknown -> EM
            // Known-phase haplotype counts (each diploid contributes 2 haplotypes):
            double kAB = 2.0 * n22 + 1.0 * (n12 + n21);
            double kAb = 2.0 * n20 + 1.0 * (n10 + n21);  // (2,0)->2Ab, (1,0)->1Ab+1ab, (2,1)->1AB+1Ab
            double kaB = 2.0 * n02 + 1.0 * (n01 + n12);  // (0,2)->2aB, (0,1)->1aB+1ab, (1,2)->1AB+1aB
            double kab = 2.0 * n00 + 1.0 * (n01 + n10);  // (0,0)->2ab, (0,1)->1aB+1ab, (1,0)->1Ab+1ab
            double ktot = kAB + kAb + kaB + kab;
            double pAB = 0.25, pAb = 0.25, paB = 0.25, pab = 0.25;
            if (ktot > 0) { pAB = kAB/ktot; pAb = kAb/ktot; paB = kaB/ktot; pab = kab/ktot; }
            // EM iterations for the n11 double hets.
            // Double het (1,1) has two possible phases:
            //   AB/ab -> 1 pAB + 1 pab   (coupling)
            //   Ab/aB -> 1 pAb + 1 paB   (repulsion)
            for (int iter = 0; iter < 100; ++iter) {
                double denom = pAB * pab + pAb * paB;
                if (denom <= 0.0) break;
                double fABab = pAB * pab / denom;   // P(AB/ab | double het)
                double fAbaB = pAb * paB / denom;   // P(Ab/aB | double het)
                // Total haplotype counts = known-phase + expected from double hets
                double nAB = kAB + fABab * n11;
                double nAb = kAb + fAbaB * n11;
                double naB = kaB + fAbaB * n11;
                double nab = kab + fABab * n11;
                double tot = nAB + nAb + naB + nab;
                if (tot <= 0.0) break;
                double nABn = nAB / tot, nAbn = nAb / tot, naBn = naB / tot, nabn = nab / tot;
                double d1 = std::fabs(nABn - pAB), d2 = std::fabs(nAbn - pAb),
                       d3 = std::fabs(naBn - paB), d4 = std::fabs(nabn - pab);
                pAB = nABn; pAb = nAbn; paB = naBn; pab = nabn;
                if (d1 < 1e-6 && d2 < 1e-6 && d3 < 1e-6 && d4 < 1e-6) break;
            }
            double pA = pAB + pAb, pB = pAB + paB;
            double da = pA * (1.0 - pA), db = pB * (1.0 - pB);
            if (da <= 0.0 || db <= 0.0) { out.push_back(std::numeric_limits<float>::quiet_NaN()); continue; }
            double D = pAB - pA * pB;
            double r = D / std::sqrt(da * db);
            out.push_back(static_cast<float>(std::fabs(r)));
        }
    }
    return out;
}

namespace detail {

std::vector<double> determine_weights(
    const SnpSet& snps, const std::vector<int>& ix, double ldcutoff, LdMode mode) {
    int n = static_cast<int>(ix.size());
    int nsamp = snps.n_b;
    std::vector<const int8_t*> rows(n);
    for (int i = 0; i < n; ++i) rows[i] = snps.dosage_row(static_cast<size_t>(ix[i]));

    std::vector<float> r;
    switch (mode) {
        case LdMode::pairwise:
            r = rogers_huff_r(rows, n, nsamp, /*skip_missing=*/true);
            break;
        case LdMode::phased:
            r = haplotype_r(rows, n, nsamp);
            break;
        case LdMode::em:
            r = em_haplotype_r(rows, n, nsamp);
            break;
        default:  // dosage_fill
            r = rogers_huff_r(rows, n, nsamp, /*skip_missing=*/false);
            break;
    }
    std::vector<double> w(n, 0.0);
    for (int i = 0; i < n; ++i) {
        int above = 0;
        for (int j = 0; j < n; ++j) {
            if (i == j) continue;
            int a = std::min(i, j), b = std::max(i, j);
            size_t k = static_cast<size_t>(a) * (2 * n - a - 1) / 2 + (b - a - 1);
            float rij = r[k];
            float r2 = rij * rij;
            if (std::isnan(rij) || r2 > static_cast<float>(ldcutoff)) ++above;
        }
        w[i] = 1.0 / (1.0 + above);
    }
    return w;
}

const double kSelCoefs[] = {
    0.0,    0.00001, 0.00005, 0.0001, 0.0002, 0.0004, 0.0006, 0.0008,
    0.001,  0.003,   0.005,   0.01,   0.05,   0.08,   0.1,    0.15};
const int kNSel = sizeof(kSelCoefs) / sizeof(kSelCoefs[0]);

double calculate_cl(double sc, const std::vector<WinRow>& dat, double ne) {
    if (!(sc >= 0.0 && sc < 1.0)) return std::numeric_limits<double>::infinity();
    double ml = 0.0;
    for (const auto& row : dat) {
        const double var = row.omega * row.p2 * (1.0 - row.p2);
        const double c = determine_c(row.rd, sc, ne);
        const double cl = chen_likelihood(row.xj, row.nj, c, row.p2, var);
        ml += row.weight * cl;
    }
    return -ml;
}

void compute_xpclr(const std::vector<WinRow>& dat, double ne,
                   bool unimodal_s, double& modelL, double& nullL, double& sel) {
    double maximum_li = std::numeric_limits<double>::infinity();
    double maxli_sc = 0.0;
    std::vector<double> lliks;
    lliks.reserve(kNSel);
    for (int i = 0; i < kNSel; ++i) {
        double ll = calculate_cl(kSelCoefs[i], dat, ne);
        lliks.push_back(ll);
        if (ll < maximum_li) {
            maximum_li = ll;
            maxli_sc = kSelCoefs[i];
        } else if (unimodal_s) {
            break;
        }
    }
    nullL = -lliks[0];
    modelL = -maximum_li;
    sel = maxli_sc;
}

}  // namespace detail

}  // namespace xpclr
