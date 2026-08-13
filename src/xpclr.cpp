#include "xpclr.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <limits>
#include <numeric>
#include <random>
#include <sstream>
#include <vector>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_gamma.h>

#ifdef _OPENMP
#include <omp.h>
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
constexpr size_t kIdBufSize = 160;
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
    static thread_local bool gsl_handler_set = false;
    if (!gsl_handler_set) {
        gsl_set_error_handler(&gsl_error_off);
        gsl_handler_set = true;
    }
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

static std::vector<double> determine_weights(
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

static const double kSelCoefs[] = {
    0.0,    0.00001, 0.00005, 0.0001, 0.0002, 0.0004, 0.0006, 0.0008,
    0.001,  0.003,   0.005,   0.01,   0.05,   0.08,   0.1,    0.15};
static const int kNSel = sizeof(kSelCoefs) / sizeof(kSelCoefs[0]);

struct WinRow {
    int xj;
    int nj;
    double rd;
    double p2;
    double omega;
    double weight;
};

static double calculate_cl(double sc, const std::vector<WinRow>& dat, double ne) {
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

static void compute_xpclr(const std::vector<WinRow>& dat, double ne,
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

static std::vector<int> choose_indices(int start_ix, int stop_ix, int maximum_size,
                                       uint64_t seed, int win_i) {
    int n = stop_ix - start_ix;
    std::vector<int> all(n);
    std::iota(all.begin(), all.end(), start_ix);
    if (n <= maximum_size) return all;
    std::seed_seq seq{static_cast<uint32_t>(seed & 0xffffffffu),
                      static_cast<uint32_t>(seed >> 32),
                      static_cast<uint32_t>(win_i),
                      static_cast<uint32_t>(start_ix),
                      static_cast<uint32_t>(stop_ix)};
    std::mt19937 rng(seq);
    for (int i = 0; i < maximum_size; ++i) {
        std::uniform_int_distribution<int> dist(i, n - 1);
        int j = dist(rng);
        std::swap(all[i], all[j]);
    }
    all.resize(maximum_size);
    std::sort(all.begin(), all.end());
    return all;
}


// Optional map: whitespace CHROM POS GDIST (header optional). Same chrom, sorted POS.
// Genetic distance at SNP POS via linear interpolation; POS outside range clamps.
static std::vector<double> genetic_distance_for_snps(const SnpSet& snps,
                                                    const Options& opt,
                                                    const std::string& chrom) {
    std::vector<double> g(snps.size(), 0.0);
    if (opt.gmap_path.empty()) {
        for (size_t i = 0; i < snps.size(); ++i)
            g[i] = static_cast<double>(snps.snps[i].pos) * opt.rrate;
        return g;
    }
    std::ifstream in(opt.gmap_path);
    if (!in) die("cannot open --gmap: " + opt.gmap_path);
    std::vector<int64_t> mpos;
    std::vector<double> mg;
    std::string line;
    int nline = 0;
    while (std::getline(in, line)) {
        ++nline;
        if (line.empty() || line[0] == '#') continue;
        std::istringstream iss(line);
        std::string c; int64_t pos = 0; double gd = 0.0;
        if (!(iss >> c >> pos >> gd)) {
            // skip header-ish lines
            continue;
        }
        if (c != chrom) continue;
        if (pos < 1 || !std::isfinite(gd)) continue;
        mpos.push_back(pos);
        mg.push_back(gd);
    }
    if (mpos.size() < 2)
        die("--gmap has <2 points for chrom " + chrom + ": " + opt.gmap_path);
    // require non-decreasing POS
    for (size_t i = 1; i < mpos.size(); ++i) {
        if (mpos[i] < mpos[i - 1])
            die("--gmap POS not sorted ascending for chrom " + chrom);
    }
    for (size_t i = 0; i < snps.size(); ++i) {
        const int64_t pos = snps.snps[i].pos;
        if (pos <= mpos.front()) { g[i] = mg.front(); continue; }
        if (pos >= mpos.back()) { g[i] = mg.back(); continue; }
        auto it = std::lower_bound(mpos.begin(), mpos.end(), pos);
        size_t hi = static_cast<size_t>(it - mpos.begin());
        size_t lo = hi - 1;
        if (mpos[hi] == pos) { g[i] = mg[hi]; continue; }
        const double tspan = static_cast<double>(mpos[hi] - mpos[lo]);
        const double w = tspan > 0.0 ? static_cast<double>(pos - mpos[lo]) / tspan : 0.0;
        g[i] = mg[lo] * (1.0 - w) + mg[hi] * w;
    }
    return g;
}

std::vector<WindowResult> xpclr_scan(const SnpSet& snps,
                                     const Options& opt,
                                     const std::string& chrom,
                                     int64_t win_start, int64_t win_stop) {
    const int nsnps = static_cast<int>(snps.size());
    std::vector<int64_t> pos(nsnps);
    for (int i = 0; i < nsnps; ++i) pos[i] = snps.snps[i].pos;

    double omega = estimate_omega(snps, opt.omega_trim);
    {
        std::ostringstream oss;
        oss << "Omega estimated as : " << std::fixed << std::setprecision(6) << omega
            << " (trim=" << opt.omega_trim << ")";
        log_info(opt, oss.str());
    }
    if (snps.size() < 50) {
        log_warn(opt, "omega estimated from only " + std::to_string(snps.size()) +
                         " SNPs on " + chrom + "; background scale may be unstable");
    }
    const std::vector<double> gdist = genetic_distance_for_snps(snps, opt, chrom);
    if (opt.gmap_path.empty()) {
        std::ostringstream oss;
        oss << "No genetic map; using rrate of " << std::scientific << opt.rrate << "/bp";
        log_info(opt, oss.str());
    } else {
        log_info(opt, "Genetic map: " + opt.gmap_path + " (chrom " + chrom + ")");
    }

    int64_t stop = win_stop;
    if (stop <= 0) stop = pos.back();
    int64_t start = win_start > 0 ? win_start : 1;
    if (start > stop) {
        log_warn(opt, "window start > stop on " + chrom + "; no windows");
        return {};
    }

    std::vector<std::pair<int64_t, int64_t>> windows;
    for (int64_t s = start; s < stop; s += opt.step) {
        windows.emplace_back(s, s - 1 + opt.size);
    }
    log_info(opt, "Windows on " + chrom + ": " + std::to_string(windows.size()) +
                      " (size=" + std::to_string(opt.size) +
                      ", step=" + std::to_string(opt.step) +
                      ", grid=[" + std::to_string(start) + "," +
                      std::to_string(stop) + "))");
    if (opt.unimodal_s) {
        log_info(opt, "Selection grid: --unimodal-s (stop at first LL decline; "
                      "hardingnj/python-like)");
    } else {
        log_info(opt, "Selection grid: full max over s (default)");
    }

    std::vector<WindowResult> out(windows.size());

#pragma omp parallel for schedule(dynamic) if (opt.threads > 1)
    for (int i = 0; i < static_cast<int>(windows.size()); ++i) {
        int64_t wstart = windows[i].first;
        int64_t wstop = windows[i].second;
        WindowResult wr;
        wr.chrom = chrom;
        wr.start = wstart;
        wr.stop = wstop;

        auto lo = std::lower_bound(pos.begin(), pos.end(), wstart);
        auto hi = std::lower_bound(pos.begin(), pos.end(), wstop);
        int start_ix = static_cast<int>(lo - pos.begin());
        int stop_ix = static_cast<int>(hi - pos.begin());
        int n_avail = stop_ix - start_ix;
        wr.nSNPs_avail = n_avail;
        wr.nSNPs = std::max(0, n_avail);

        if (n_avail < opt.minsnps) {
            wr.valid = false;
            out[i] = std::move(wr);
            continue;
        }

        auto ix = choose_indices(start_ix, stop_ix, opt.maxsnps, opt.seed, i);
        wr.nSNPs = static_cast<int>(ix.size());
        if (wr.nSNPs < opt.minsnps) {
            wr.valid = false;
            out[i] = std::move(wr);
            continue;
        }
        wr.pos_start = snps.snps[ix.front()].pos;
        wr.pos_stop = snps.snps[ix.back()].pos;

        auto weights = determine_weights(snps, ix, opt.ldcutoff, opt.ld_mode);

        // rd: distance from each SNP to the window center (genetic position).
        double center_gd;
        if (opt.gmap_path.empty()) {
            center_gd = static_cast<double>(wstart + wstop) / 2.0 * opt.rrate;
        } else {
            // interpolate gdist at window center between flanking SNPs
            int64_t center_bp = static_cast<int64_t>(
                static_cast<double>(wstart + wstop) / 2.0);
            auto it = std::lower_bound(pos.begin(), pos.end(), center_bp);
            size_t hi = static_cast<size_t>(it - pos.begin());
            if (hi == 0) {
                center_gd = gdist[0];
            } else if (hi >= pos.size()) {
                center_gd = gdist[pos.size() - 1];
            } else if (pos[hi] == center_bp) {
                center_gd = gdist[hi];
            } else {
                size_t lo = hi - 1;
                double span = static_cast<double>(pos[hi] - pos[lo]);
                double w = span > 0.0
                    ? static_cast<double>(center_bp - pos[lo]) / span : 0.0;
                center_gd = gdist[lo] * (1.0 - w) + gdist[hi] * w;
            }
        }

        std::vector<WinRow> dat(ix.size());
        for (size_t k = 0; k < ix.size(); ++k) {
            const auto& s = snps.snps[ix[k]];
            dat[k].xj = s.x_alt;
            dat[k].nj = s.n_a;
            dat[k].rd = std::fabs(gdist[static_cast<size_t>(ix[k])] - center_gd);
            dat[k].p2 = s.q2;
            dat[k].omega = omega;
            dat[k].weight = weights[k];
        }

        compute_xpclr(dat, opt.ne, opt.unimodal_s, wr.modelL, wr.nullL, wr.sel_coef);
        wr.valid = true;
        out[i] = std::move(wr);
    }

    return out;
}

void write_results(const std::string& path, const std::vector<WindowResult>& rows) {
    std::vector<double> xp(rows.size(), std::numeric_limits<double>::quiet_NaN());
    for (size_t i = 0; i < rows.size(); ++i) {
        if (!rows[i].valid) continue;
        xp[i] = 2.0 * (rows[i].modelL - rows[i].nullL);
    }
    // Population mean/sd of finite xpclr (same as two-pass mean then var/cnt).
    double sum = 0.0, sumsq = 0.0, cnt = 0.0;
    for (double v : xp) {
        if (!std::isfinite(v)) continue;
        sum += v;
        sumsq += v * v;
        cnt += 1.0;
    }
    const double mean = cnt > 0 ? sum / cnt : 0.0;
    const double var = cnt > 0 ? (sumsq - sum * mean) : 0.0;  // = sum (v-mean)^2
    const double sd =
        cnt > 0 ? std::sqrt(var / cnt) : std::numeric_limits<double>::quiet_NaN();

    std::ofstream out(path);
    if (!out) die("cannot write output: " + path);
    out << "id\tchrom\tstart\tstop\tpos_start\tpos_stop\tmodelL\tnullL\tsel_coef\t"
           "nSNPs\tnSNPs_avail\txpclr\txpclr_norm\n";
    out << std::setprecision(12);
    for (size_t i = 0; i < rows.size(); ++i) {
        const auto& r = rows[i];
        char idbuf[kIdBufSize];
        std::snprintf(idbuf, sizeof(idbuf), "%s_%08lld_%08lld", r.chrom.c_str(),
                      static_cast<long long>(r.start),
                      static_cast<long long>(r.stop));
        out << idbuf << "\t" << r.chrom << "\t" << r.start << "\t" << r.stop << "\t";
        if (!r.valid) {
            out << "nan\tnan\tnan\tnan\tnan\t" << r.nSNPs << "\t" << r.nSNPs_avail
                << "\tnan\tnan\n";
            continue;
        }
        double x = xp[i];
        double xn =
            (std::isfinite(x) && sd > 0) ? (x - mean) / sd : std::numeric_limits<double>::quiet_NaN();
        out << r.pos_start << "\t" << r.pos_stop << "\t" << r.modelL << "\t" << r.nullL
            << "\t" << r.sel_coef << "\t" << r.nSNPs << "\t" << r.nSNPs_avail << "\t"
            << x << "\t" << xn << "\n";
    }
}

int run_xpclr(const Options& opt) {
#ifdef _OPENMP
    omp_set_num_threads(opt.threads);
#endif
    // One open + header + index for the whole run (all contigs/regions).
    VcfSession* vcf = vcf_session_open(opt);
    const VcfHeaderInfo& hdr_info = vcf_session_info(vcf);

    auto pop = load_pop_file(opt.pop_file, opt);
    auto plan = resolve_samples(hdr_info.samples, pop, opt);

    std::vector<RegionTarget> targets;
    if (opt.region.empty()) {
        if (hdr_info.contigs.empty()) {
            vcf_session_close(vcf);
            die("VCF has no contigs in header: " + opt.vcf);
        }
        log_info(opt, "No -r/--regions: scanning all " +
                          std::to_string(hdr_info.contigs.size()) + " contigs");
        for (const auto& c : hdr_info.contigs) {
            RegionTarget t;
            t.chrom = c;
            targets.push_back(std::move(t));
        }
    } else {
        targets.push_back(parse_region_string(opt.region));
        log_info(opt, "Region: " + opt.region);
    }

    std::vector<WindowResult> all_rows;
    int n_ok = 0;
    for (const auto& t : targets) {
        auto snps = load_snps(vcf, opt, plan, t);
        if (snps.empty()) {
            log_warn(opt, "skip " + t.chrom + ": no SNPs after filters");
            continue;
        }
        int64_t win_start = t.has_beg ? t.beg : 1;
        int64_t win_stop = t.has_end ? t.end : 0;
        auto rows = xpclr_scan(snps, opt, t.chrom, win_start, win_stop);
        all_rows.insert(all_rows.end(), std::make_move_iterator(rows.begin()),
                        std::make_move_iterator(rows.end()));
        ++n_ok;
    }

    vcf_session_close(vcf);

    if (all_rows.empty())
        die("no windows produced (no usable SNPs in selected region(s))");
    write_results(opt.out, all_rows);
    log_info(opt, "Analysis complete. Output file " + opt.out + " (" +
                      std::to_string(all_rows.size()) + " windows, " +
                      std::to_string(n_ok) + " contigs)");
    return 0;
}

}  // namespace xpclr
