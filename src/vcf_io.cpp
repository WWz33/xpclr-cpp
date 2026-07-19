#include "xpclr.hpp"

#include <htslib/bgzf.h>
#include <htslib/hts.h>
#include <htslib/tbx.h>
#include <htslib/vcf.h>

#include <cmath>
#include <cstring>
#include <sstream>
#include <string>
#include <vector>

namespace xpclr {

// bgzf/htslib decompression threads (IO phase; separate in time from OpenMP scan).
static void set_io_threads(htsFile* fp, int n) {
    if (!fp || n <= 1) return;
    if (hts_set_threads(fp, n) < 0) {
        // non-fatal: older files / plain VCF still work single-threaded
    }
}

// Genomic interval for indexed fetch (1-based inclusive POS on both ends when closed).
// When --stop is set, extend by --size so the last sliding window can still pull SNPs.
// When --stop is 0, open-ended on the right (whole contig from --start).
static void load_region_bounds(const Options& opt, int64_t& pos_lo, int64_t& pos_hi,
                               bool& open_end) {
    pos_lo = opt.start > 0 ? opt.start : 1;
    open_end = (opt.stop <= 0);
    if (open_end) {
        pos_hi = 0;
    } else {
        // last window [s, s-1+size] with s < stop => SNP upper exclusive ~ stop+size-1
        int64_t hi = opt.stop + opt.size;
        if (hi < pos_lo) hi = pos_lo;
        pos_hi = hi;
    }
}

static std::string region_string(const Options& opt) {
    int64_t lo = 0, hi = 0;
    bool open_end = false;
    load_region_bounds(opt, lo, hi, open_end);
    std::ostringstream oss;
    if (open_end) {
        // Whole contig query; process_rec still drops POS < --start.
        // Avoid "chr:start-" which some htslib builds reject.
        if (lo <= 1)
            oss << opt.chrom;
        else
            oss << opt.chrom << ":" << lo << "-";
    } else {
        oss << opt.chrom << ":" << lo << "-" << hi;
    }
    return oss.str();
}

std::vector<std::string> read_vcf_samples(const std::string& path) {
    htsFile* fp = bcf_open(path.c_str(), "r");
    if (!fp) die("cannot open VCF/BCF: " + path);
    bcf_hdr_t* hdr = bcf_hdr_read(fp);
    if (!hdr) {
        bcf_close(fp);
        die("cannot read VCF header: " + path);
    }
    int n = bcf_hdr_nsamples(hdr);
    std::vector<std::string> samples;
    samples.reserve(n);
    for (int i = 0; i < n; ++i) samples.emplace_back(hdr->samples[i]);
    bcf_hdr_destroy(hdr);
    bcf_close(fp);
    return samples;
}

static bool is_snp_biallelic(const bcf_hdr_t* hdr, bcf1_t* rec) {
    if (rec->n_allele != 2) return false;
    const char* ref = rec->d.allele[0];
    const char* alt = rec->d.allele[1];
    if (!ref || !alt) return false;
    // single-base alleles only (SNP)
    if (std::strlen(ref) != 1 || std::strlen(alt) != 1) return false;
    if (ref[0] == '.' || alt[0] == '.') return false;
    (void)hdr;
    return true;
}

// Count alleles for a subset of samples. Returns alt count and total called alleles.
// Multiallelic GT (allele >1) marks multi=true.
static void count_pop(const int32_t* gt, int nsmpl, const std::vector<int>& idx,
                      int& alt, int& ncall, bool& multi, bool& any_allele_gt1) {
    alt = 0;
    ncall = 0;
    multi = false;
    any_allele_gt1 = false;
    for (int si : idx) {
        if (si < 0 || si >= nsmpl) continue;
        int32_t a0 = gt[si * 2];
        int32_t a1 = gt[si * 2 + 1];
        if (bcf_gt_is_missing(a0) || bcf_gt_is_missing(a1)) continue;
        int aa = bcf_gt_allele(a0);
        int bb = bcf_gt_allele(a1);
        if (aa < 0 || bb < 0) continue;
        if (aa > 1 || bb > 1) {
            any_allele_gt1 = true;
            multi = true;
        }
        // still count only 0/1 for biallelic path; multi sites dropped later
        if (aa <= 1 && bb <= 1) {
            alt += (aa == 1) + (bb == 1);
            ncall += 2;
        }
    }
}

// dosage: n_alt for diploid; missing -> 0 (Python fill=0)
static void dosage_pop(const int32_t* gt, int nsmpl, const std::vector<int>& idx,
                       std::vector<int8_t>& out) {
    out.resize(idx.size());
    for (size_t i = 0; i < idx.size(); ++i) {
        int si = idx[i];
        int32_t a0 = gt[si * 2];
        int32_t a1 = gt[si * 2 + 1];
        if (bcf_gt_is_missing(a0) || bcf_gt_is_missing(a1)) {
            out[i] = 0;
            continue;
        }
        int aa = bcf_gt_allele(a0);
        int bb = bcf_gt_allele(a1);
        if (aa < 0 || bb < 0) {
            out[i] = 0;
            continue;
        }
        int d = 0;
        if (aa == 1) ++d;
        if (bb == 1) ++d;
        out[i] = static_cast<int8_t>(d);
    }
}

std::vector<SnpData> load_snps(const Options& opt, const SamplePlan& plan) {
    htsFile* fp = bcf_open(opt.vcf.c_str(), "r");
    if (!fp) die("cannot open VCF/BCF: " + opt.vcf);
    set_io_threads(fp, opt.threads);

    bcf_hdr_t* hdr = bcf_hdr_read(fp);
    if (!hdr) {
        bcf_close(fp);
        die("cannot read VCF header: " + opt.vcf);
    }

    // optional region via index
    tbx_t* tbx = nullptr;
    hts_idx_t* idx = nullptr;
    hts_itr_t* itr = nullptr;
    int tid = bcf_hdr_name2id(hdr, opt.chrom.c_str());
    if (tid < 0) {
        bcf_hdr_destroy(hdr);
        bcf_close(fp);
        die("chromosome not found in VCF header: " + opt.chrom);
    }

    int64_t pos_lo = 1, pos_hi = 0;
    bool open_end = true;
    load_region_bounds(opt, pos_lo, pos_hi, open_end);
    const std::string region = region_string(opt);

    // Prefer CSI/TBI; query chrom[:start-stop(+size)] to avoid loading whole contig.
    if (fp->format.format == bcf) {
        idx = bcf_index_load(opt.vcf.c_str());
        if (idx) {
            itr = bcf_itr_querys(idx, hdr, region.c_str());
            if (!itr) {
                // fallback: whole contig (e.g. empty region)
                itr = bcf_itr_querys(idx, hdr, opt.chrom.c_str());
            }
        }
    } else {
        tbx = tbx_index_load(opt.vcf.c_str());
        if (tbx) {
            itr = tbx_itr_querys(tbx, region.c_str());
            if (!itr) itr = tbx_itr_querys(tbx, opt.chrom.c_str());
        }
    }

    if (itr) {
        log_info(opt, "Indexed VCF load region: " + region +
                          " (IO threads=" + std::to_string(opt.threads) + ")");
    } else {
        log_warn(opt, "no index found for " + opt.vcf +
                          "; scanning file for chrom=" + opt.chrom +
                          " (consider: bcftools index -t " + opt.vcf + ")");
        if (!open_end) {
            log_info(opt, "Will keep only POS in [" + std::to_string(pos_lo) +
                              "," + std::to_string(pos_hi) + "] while scanning");
        }
    }

    bcf1_t* rec = bcf_init();
    int32_t* gt_arr = nullptr;
    int ngt_arr = 0;

    std::vector<SnpData> snps;
    snps.reserve(100000);

    int64_t n_total = 0, n_multi = 0, n_missing_pop = 0, n_fixed_p2 = 0,
            n_kept = 0, n_out_of_region = 0;

    auto process_rec = [&](bcf1_t* r) {
        if (r->rid != tid) return;
        // 1-based POS; drop sites outside load interval (safety for full-chrom query
        // and for unindexed scans).
        int64_t pos1 = static_cast<int64_t>(r->pos) + 1;
        if (pos1 < pos_lo) return;
        if (!open_end && pos1 > pos_hi) {
            ++n_out_of_region;
            return;
        }
        if (bcf_unpack(r, BCF_UN_STR | BCF_UN_FMT) < 0) return;
        ++n_total;

        // multi-allelic / non-SNP
        if (!is_snp_biallelic(hdr, r)) {
            ++n_multi;
            return;
        }

        int ngt = bcf_get_genotypes(hdr, r, &gt_arr, &ngt_arr);
        if (ngt <= 0) {
            ++n_missing_pop;
            return;
        }
        int nsmpl = bcf_hdr_nsamples(hdr);
        if (ngt < nsmpl * 2) {
            // not diploid layout; skip
            ++n_missing_pop;
            return;
        }

        int alt_a = 0, n_a = 0, alt_b = 0, n_b = 0;
        bool multi_a = false, multi_b = false, any1 = false, any2 = false;
        count_pop(gt_arr, nsmpl, plan.idx_a, alt_a, n_a, multi_a, any1);
        count_pop(gt_arr, nsmpl, plan.idx_b, alt_b, n_b, multi_b, any2);
        if (multi_a || multi_b || any1 || any2) {
            // biallelic site shouldn't hit; safety
            ++n_multi;
            return;
        }
        if (n_a == 0 || n_b == 0) {
            ++n_missing_pop;
            return;
        }

        // fixed / singleton in pop2 (hardingnj bin/xpclr)
        // ac2.is_non_segregating() | is_singleton(0) | is_singleton(1)
        // non-segregating: only one allele present among called
        // singleton(0): exactly one REF allele among all alleles? scikit-allel:
        // is_singleton(allele) true if count of that allele == 1
        int ref_b = n_b - alt_b;
        bool fixed_p2 = false;
        if (alt_b == 0 || ref_b == 0) fixed_p2 = true;           // non-segregating
        if (alt_b == 1 || ref_b == 1) fixed_p2 = true;           // singleton allele
        if (fixed_p2) {
            ++n_fixed_p2;
            return;
        }

        SnpData s;
        s.pos = r->pos + 1;  // htslib 0-based -> 1-based POS
        s.x_alt = alt_a;
        s.n_a = n_a;
        s.n_b = n_b;
        s.q2 = static_cast<double>(alt_b) / static_cast<double>(n_b);
        dosage_pop(gt_arr, nsmpl, plan.idx_b, s.dosage_b);
        snps.push_back(std::move(s));
        ++n_kept;
    };

    if (itr) {
        if (tbx) {
            kstring_t str = {0, 0, nullptr};
            while (tbx_itr_next(fp, tbx, itr, &str) >= 0) {
                if (vcf_parse(&str, hdr, rec) < 0) continue;
                process_rec(rec);
            }
            free(str.s);
        } else {
            while (bcf_itr_next(fp, itr, rec) >= 0) process_rec(rec);
        }
        hts_itr_destroy(itr);
    } else {
        while (bcf_read(fp, hdr, rec) == 0) process_rec(rec);
    }

    free(gt_arr);
    bcf_destroy(rec);
    if (tbx) tbx_destroy(tbx);
    if (idx) hts_idx_destroy(idx);
    bcf_hdr_destroy(hdr);
    bcf_close(fp);

    log_info(opt, std::to_string(n_total) + " records considered on " + opt.chrom +
                      " (load " + region + ")");
    if (n_out_of_region > 0) {
        log_info(opt, std::to_string(n_out_of_region) +
                          " records skipped outside load interval");
    }
    log_info(opt, std::to_string(n_multi) + " SNPs excluded as multiallelic/non-SNP");
    log_info(opt, std::to_string(n_missing_pop) +
                      " SNPs excluded as missing in all samples in a population");
    log_info(opt, std::to_string(n_fixed_p2) +
                      " SNPs excluded as invariant or singleton in population 2");
    log_info(opt, std::to_string(n_kept) + "/" + std::to_string(n_total) +
                      " SNPs included in the analysis");

    if (snps.empty()) die("no SNPs left after filters on chrom " + opt.chrom);
    return snps;
}

}  // namespace xpclr
