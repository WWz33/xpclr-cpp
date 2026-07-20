#include "xpclr.hpp"

#include <htslib/hts.h>
#include <htslib/tbx.h>
#include <htslib/vcf.h>

#include <cmath>
#include <cstring>
#include <sstream>
#include <string>
#include <vector>

namespace xpclr {

static void set_io_threads(htsFile* fp, int n) {
    if (!fp || n <= 1) return;
    (void)hts_set_threads(fp, n);
}

// Load interval for indexed fetch / POS filter (1-based).
// When region has end, extend by --size so last sliding window still sees SNPs.
static void load_bounds(const RegionTarget& t, int64_t size, int64_t& pos_lo,
                        int64_t& pos_hi, bool& open_end) {
    pos_lo = t.has_beg ? t.beg : 1;
    open_end = !t.has_end;
    if (open_end) {
        pos_hi = 0;
    } else {
        pos_hi = t.end + size;
        if (pos_hi < pos_lo) pos_hi = pos_lo;
    }
}

static std::string hts_region_query(const RegionTarget& t, int64_t size) {
    int64_t lo = 0, hi = 0;
    bool open_end = true;
    load_bounds(t, size, lo, hi, open_end);
    std::ostringstream oss;
    if (!t.has_beg && !t.has_end) {
        oss << t.chrom;
    } else if (open_end) {
        if (lo <= 1)
            oss << t.chrom;
        else
            oss << t.chrom << ":" << lo << "-";
    } else {
        oss << t.chrom << ":" << lo << "-" << hi;
    }
    return oss.str();
}

VcfHeaderInfo read_vcf_header_info(const std::string& path) {
    htsFile* fp = bcf_open(path.c_str(), "r");
    if (!fp) die("cannot open VCF/BCF: " + path);
    bcf_hdr_t* hdr = bcf_hdr_read(fp);
    if (!hdr) {
        bcf_close(fp);
        die("cannot read VCF header: " + path);
    }
    VcfHeaderInfo info;
    int n = bcf_hdr_nsamples(hdr);
    info.samples.reserve(static_cast<size_t>(n));
    for (int i = 0; i < n; ++i) info.samples.emplace_back(hdr->samples[i]);

    int nctg = hdr->n[BCF_DT_CTG];
    info.contigs.reserve(static_cast<size_t>(nctg));
    for (int i = 0; i < nctg; ++i) {
        const char* name = hdr->id[BCF_DT_CTG][i].key;
        if (name) info.contigs.emplace_back(name);
    }
    bcf_hdr_destroy(hdr);
    bcf_close(fp);
    return info;
}

static bool is_snp_biallelic(const bcf1_t* rec) {
    if (rec->n_allele != 2) return false;
    const char* ref = rec->d.allele[0];
    const char* alt = rec->d.allele[1];
    if (!ref || !alt) return false;
    if (std::strlen(ref) != 1 || std::strlen(alt) != 1) return false;
    if (ref[0] == '.' || alt[0] == '.') return false;
    return true;
}

static void count_pop(const int32_t* gt, int nsmpl, const std::vector<int>& idx,
                      int& alt, int& ncall, bool& multi) {
    alt = 0;
    ncall = 0;
    multi = false;
    for (int si : idx) {
        if (si < 0 || si >= nsmpl) continue;
        const int32_t a0 = gt[si * 2];
        const int32_t a1 = gt[si * 2 + 1];
        if (bcf_gt_is_missing(a0) || bcf_gt_is_missing(a1)) continue;
        const int aa = bcf_gt_allele(a0);
        const int bb = bcf_gt_allele(a1);
        if (aa < 0 || bb < 0) continue;
        if (aa > 1 || bb > 1) multi = true;
        if (aa <= 1 && bb <= 1) {
            alt += (aa == 1) + (bb == 1);
            ncall += 2;
        }
    }
}

// Dosage for LD weights: missing alleles → 0 (hardingnj unphased path).
static void dosage_pop(const int32_t* gt, int nsmpl, const std::vector<int>& idx,
                       std::vector<int8_t>& out) {
    out.resize(idx.size());
    for (size_t i = 0; i < idx.size(); ++i) {
        const int si = idx[i];
        if (si < 0 || si >= nsmpl) {
            out[i] = 0;
            continue;
        }
        const int32_t a0 = gt[si * 2];
        const int32_t a1 = gt[si * 2 + 1];
        if (bcf_gt_is_missing(a0) || bcf_gt_is_missing(a1)) {
            out[i] = 0;
            continue;
        }
        const int aa = bcf_gt_allele(a0);
        const int bb = bcf_gt_allele(a1);
        if (aa < 0 || bb < 0 || aa > 1 || bb > 1) {
            out[i] = 0;
            continue;
        }
        out[i] = static_cast<int8_t>((aa == 1) + (bb == 1));
    }
}

std::vector<SnpData> load_snps(const Options& opt, const SamplePlan& plan,
                               const RegionTarget& target) {
    if (target.chrom.empty()) die("load_snps: empty contig");

    htsFile* fp = bcf_open(opt.vcf.c_str(), "r");
    if (!fp) die("cannot open VCF/BCF: " + opt.vcf);
    set_io_threads(fp, opt.threads);

    bcf_hdr_t* hdr = bcf_hdr_read(fp);
    if (!hdr) {
        bcf_close(fp);
        die("cannot read VCF header: " + opt.vcf);
    }

    tbx_t* tbx = nullptr;
    hts_idx_t* idx = nullptr;
    hts_itr_t* itr = nullptr;
    int tid = bcf_hdr_name2id(hdr, target.chrom.c_str());
    if (tid < 0) {
        bcf_hdr_destroy(hdr);
        bcf_close(fp);
        die("chromosome not found in VCF header: " + target.chrom);
    }

    int64_t pos_lo = 1, pos_hi = 0;
    bool open_end = true;
    load_bounds(target, opt.size, pos_lo, pos_hi, open_end);
    const std::string region = hts_region_query(target, opt.size);

    if (fp->format.format == bcf) {
        idx = bcf_index_load(opt.vcf.c_str());
        if (idx) {
            itr = bcf_itr_querys(idx, hdr, region.c_str());
            if (!itr) itr = bcf_itr_querys(idx, hdr, target.chrom.c_str());
        }
    } else {
        tbx = tbx_index_load(opt.vcf.c_str());
        if (tbx) {
            itr = tbx_itr_querys(tbx, region.c_str());
            if (!itr) itr = tbx_itr_querys(tbx, target.chrom.c_str());
        }
    }

    if (itr) {
        log_info(opt, "Indexed VCF load region: " + region +
                          " (IO threads=" + std::to_string(opt.threads) + ")");
    } else {
        log_warn(opt, "no index for " + opt.vcf + "; scanning for " +
                          target.chrom + " (bcftools index -t recommended)");
    }

    bcf1_t* rec = bcf_init();
    int32_t* gt_arr = nullptr;
    int ngt_arr = 0;

    std::vector<SnpData> snps;
    snps.reserve(100000);

    int64_t n_total = 0, n_multi = 0, n_missing_pop = 0, n_fixed_p2 = 0,
            n_kept = 0;

    auto process_rec = [&](bcf1_t* r) {
        if (r->rid != tid) return;
        int64_t pos1 = static_cast<int64_t>(r->pos) + 1;
        if (pos1 < pos_lo) return;
        if (!open_end && pos1 > pos_hi) return;

        if (bcf_unpack(r, BCF_UN_STR | BCF_UN_FMT) < 0) return;
        ++n_total;

        if (!is_snp_biallelic(r)) {
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
            ++n_missing_pop;
            return;
        }

        int alt_a = 0, n_a = 0, alt_b = 0, n_b = 0;
        bool multi_a = false, multi_b = false;
        count_pop(gt_arr, nsmpl, plan.idx_a, alt_a, n_a, multi_a);
        count_pop(gt_arr, nsmpl, plan.idx_b, alt_b, n_b, multi_b);
        if (multi_a || multi_b) {
            ++n_multi;
            return;
        }
        if (n_a == 0 || n_b == 0) {
            ++n_missing_pop;
            return;
        }

        int ref_b = n_b - alt_b;
        bool fixed_p2 = false;
        if (alt_b == 0 || ref_b == 0) fixed_p2 = true;
        if (alt_b == 1 || ref_b == 1) fixed_p2 = true;
        if (fixed_p2) {
            ++n_fixed_p2;
            return;
        }

        SnpData s;
        s.pos = pos1;
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

    log_info(opt, std::to_string(n_total) + " records considered on " +
                      target.chrom + " (load " + region + ")");
    log_info(opt, std::to_string(n_multi) + " SNPs excluded as multiallelic/non-SNP");
    log_info(opt, std::to_string(n_missing_pop) +
                      " SNPs excluded as missing in all samples in a population");
    log_info(opt, std::to_string(n_fixed_p2) +
                      " SNPs excluded as invariant or singleton in population 2");
    log_info(opt, std::to_string(n_kept) + "/" + std::to_string(n_total) +
                      " SNPs included in the analysis");

    if (snps.empty()) {
        log_warn(opt, "no SNPs left after filters on " + target.chrom + " (" +
                          region + ")");
        return snps;
    }
    return snps;
}

}  // namespace xpclr
