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
    // epsabs=0, epsrel as Python
    int status = gsl_integration_qags(&F, a, b, 0.0, kQuadEpsRel, kQuadLimit, ws,
                                      &result, &abserr);
    if (status != GSL_SUCCESS) {
        // fallback: coarser tolerance once
        status = gsl_integration_qags(&F, a, b, 1e-8, 1e-2, kQuadLimit, ws, &result,
                                      &abserr);
    }
    if (status != GSL_SUCCESS || !std::isfinite(result) || result < 0.0) {
        if (!std::isfinite(result) || result < 0.0) result = 0.0;
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

double estimate_omega(const std::vector<SnpData>& snps) {
    double sum = 0.0;
    size_t n = 0;
    for (const auto& s : snps) {
        if (!(s.q2 > 0.0 && s.q2 < 1.0) || s.n_a <= 0) continue;
        const double q1 = static_cast<double>(s.x_alt) / static_cast<double>(s.n_a);
        sum += (q1 - s.q2) * (q1 - s.q2) / (s.q2 * (1.0 - s.q2));
        ++n;
    }
    if (n == 0) die("estimate_omega: no valid SNPs");
    return sum / static_cast<double>(n);
}

// Pearson corr pairwise (Rogers-Huff r on dosage), condensed upper triangle
static std::vector<float> rogers_huff_r(const std::vector<const int8_t*>& rows,
                                        int nvar, int nsamp) {
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

static std::vector<double> determine_weights(
    const std::vector<SnpData>& snps, const std::vector<int>& ix, double ldcutoff) {
    int n = static_cast<int>(ix.size());
    int nsamp = snps.empty() ? 0 : static_cast<int>(snps[ix[0]].dosage_b.size());
    std::vector<const int8_t*> rows(n);
    for (int i = 0; i < n; ++i) rows[i] = snps[ix[i]].dosage_b.data();

    auto r = rogers_huff_r(rows, n, nsamp);
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

static double calculate_cl(double sc, const std::vector<WinRow>& dat) {
    if (!(sc >= 0.0 && sc < 1.0)) return std::numeric_limits<double>::infinity();
    double ml = 0.0;
    for (const auto& row : dat) {
        const double var = row.omega * row.p2 * (1.0 - row.p2);
        const double c = determine_c(row.rd, sc);
        const double cl = chen_likelihood(row.xj, row.nj, c, row.p2, var);
        ml += row.weight * cl;
    }
    return -ml;
}

static void compute_xpclr(const std::vector<WinRow>& dat, bool unimodal_s,
                          double& modelL, double& nullL, double& sel) {
    double maximum_li = std::numeric_limits<double>::infinity();
    double maxli_sc = 0.0;
    std::vector<double> lliks;
    lliks.reserve(kNSel);
    for (int i = 0; i < kNSel; ++i) {
        double ll = calculate_cl(kSelCoefs[i], dat);
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

std::vector<WindowResult> xpclr_scan(const std::vector<SnpData>& snps,
                                     const Options& opt,
                                     const std::string& chrom,
                                     int64_t win_start, int64_t win_stop) {
    const int nsnps = static_cast<int>(snps.size());
    std::vector<int64_t> pos(nsnps);
    for (int i = 0; i < nsnps; ++i) pos[i] = snps[i].pos;

    double omega = estimate_omega(snps);
    {
        std::ostringstream oss;
        oss << "Omega estimated as : " << std::fixed << std::setprecision(6) << omega;
        log_info(opt, oss.str());
    }
    {
        std::ostringstream oss;
        oss << "No genetic distance provided; using rrate of " << std::scientific
            << opt.rrate << "/bp";
        log_info(opt, oss.str());
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
            out[i] = wr;
            continue;
        }

        auto ix = choose_indices(start_ix, stop_ix, opt.maxsnps, opt.seed, i);
        wr.nSNPs = static_cast<int>(ix.size());
        if (wr.nSNPs < opt.minsnps) {
            wr.valid = false;
            out[i] = wr;
            continue;
        }
        wr.pos_start = snps[ix.front()].pos;
        wr.pos_stop = snps[ix.back()].pos;

        auto weights = determine_weights(snps, ix, opt.ldcutoff);

        std::vector<double> dq(ix.size());
        double mean_dq = 0.0;
        for (size_t k = 0; k < ix.size(); ++k) {
            dq[k] = static_cast<double>(snps[ix[k]].pos) * opt.rrate;
            mean_dq += dq[k];
        }
        mean_dq /= static_cast<double>(ix.size());

        std::vector<WinRow> dat(ix.size());
        for (size_t k = 0; k < ix.size(); ++k) {
            const auto& s = snps[ix[k]];
            dat[k].xj = s.x_alt;
            dat[k].nj = s.n_a;
            dat[k].rd = std::fabs(dq[k] - mean_dq);
            dat[k].p2 = s.q2;
            dat[k].omega = omega;
            dat[k].weight = weights[k];
        }

        compute_xpclr(dat, opt.unimodal_s, wr.modelL, wr.nullL, wr.sel_coef);
        wr.valid = true;
        out[i] = wr;
    }

    return out;
}

void write_results(const std::string& path, const std::vector<WindowResult>& rows) {
    std::vector<double> xp(rows.size(), std::numeric_limits<double>::quiet_NaN());
    for (size_t i = 0; i < rows.size(); ++i) {
        if (!rows[i].valid) continue;
        xp[i] = 2.0 * (rows[i].modelL - rows[i].nullL);
    }
    double mean = 0.0, cnt = 0.0;
    for (double v : xp) {
        if (std::isfinite(v)) {
            mean += v;
            cnt += 1.0;
        }
    }
    if (cnt > 0) mean /= cnt;
    double var = 0.0;
    for (double v : xp) {
        if (std::isfinite(v)) {
            double d = v - mean;
            var += d * d;
        }
    }
    double sd = cnt > 0 ? std::sqrt(var / cnt) : std::numeric_limits<double>::quiet_NaN();

    std::ofstream out(path);
    if (!out) die("cannot write output: " + path);
    out << "id\tchrom\tstart\tstop\tpos_start\tpos_stop\tmodelL\tnullL\tsel_coef\t"
           "nSNPs\tnSNPs_avail\txpclr\txpclr_norm\n";
    out << std::setprecision(12);
    for (size_t i = 0; i < rows.size(); ++i) {
        const auto& r = rows[i];
        char idbuf[160];
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
    auto hdr_info = read_vcf_header_info(opt.vcf);
    auto pop = load_pop_file(opt.pop_file, opt);
    auto plan = resolve_samples(hdr_info.samples, pop, opt);

    std::vector<RegionTarget> targets;
    if (opt.region.empty()) {
        if (hdr_info.contigs.empty()) die("VCF has no contigs in header: " + opt.vcf);
        log_info(opt, "No -r/--regions: scanning all " +
                          std::to_string(hdr_info.contigs.size()) + " contigs");
        for (auto& c : hdr_info.contigs) {
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
        auto snps = load_snps(opt, plan, t);
        if (snps.empty()) {
            log_warn(opt, "skip " + t.chrom + ": no SNPs after filters");
            continue;
        }
        int64_t win_start = t.has_beg ? t.beg : 1;
        int64_t win_stop = t.has_end ? t.end : 0;
        auto rows = xpclr_scan(snps, opt, t.chrom, win_start, win_stop);
        all_rows.insert(all_rows.end(), rows.begin(), rows.end());
        ++n_ok;
    }

    if (all_rows.empty())
        die("no windows produced (no usable SNPs in selected region(s))");
    write_results(opt.out, all_rows);
    log_info(opt, "Analysis complete. Output file " + opt.out + " (" +
                      std::to_string(all_rows.size()) + " windows, " +
                      std::to_string(n_ok) + " contigs)");
    return 0;
}

}  // namespace xpclr
