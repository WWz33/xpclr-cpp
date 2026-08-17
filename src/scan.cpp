#include "xpclr.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <numeric>
#include <random>
#include <sstream>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace xpclr {

namespace {
constexpr size_t kIdBufSize = 160;
}  // namespace

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
    while (std::getline(in, line)) {
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

        auto weights = detail::determine_weights(snps, ix, opt.ldcutoff, opt.ld_mode);

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

        std::vector<detail::WinRow> dat(ix.size());
        for (size_t k = 0; k < ix.size(); ++k) {
            const auto& s = snps.snps[ix[k]];
            dat[k].x_alt = s.x_alt;
            dat[k].n_a = s.n_a;
            dat[k].rd = std::fabs(gdist[static_cast<size_t>(ix[k])] - center_gd);
            dat[k].q2 = s.q2;
            dat[k].omega = omega;
            dat[k].weight = weights[k];
        }

        detail::compute_xpclr(dat, opt.ne, opt.unimodal_s, wr.modelL, wr.nullL, wr.sel_coef);
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
    out.flush();
    if (!out.good()) die("write failed (disk full or stream error): " + path);
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

    {
        const char* mode_name = "unknown";
        switch (opt.ld_mode) {
            case LdMode::dosage_fill: mode_name = "dosage-fill (unphased, missing->0)"; break;
            case LdMode::phased:      mode_name = "phased (haplotype r)"; break;
            case LdMode::em:           mode_name = "EM (two-locus phase inference)"; break;
            case LdMode::pairwise:    mode_name = "pairwise-complete (skip missing per pair)"; break;
        }
        log_info(opt, std::string("LD weight mode: ") + mode_name);
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
