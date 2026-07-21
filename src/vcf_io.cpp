#include "xpclr.hpp"

#include <htslib/bgzf.h>
#include <htslib/hfile.h>
#include <htslib/hts.h>
#include <htslib/tbx.h>
#include <htslib/vcf.h>

#include <cmath>
#include <cstring>
#include <sstream>
#include <string>
#include <vector>

namespace xpclr {

namespace {
constexpr size_t kSnpReserveHint = 100000;
}  // namespace

struct VcfSession {
    std::string path;
    htsFile* fp = nullptr;
    bcf_hdr_t* hdr = nullptr;
    tbx_t* tbx = nullptr;
    hts_idx_t* csi = nullptr;
    bool is_bcf = false;
    bool has_index = false;
    int64_t data_offset = -1;  // file offset of first record after header
    VcfHeaderInfo info;
};

static void set_io_threads(htsFile* fp, int n) {
    if (!fp || n <= 1) return;
    (void)hts_set_threads(fp, n);
}

static int64_t hts_tell_off(htsFile* fp) {
    if (!fp) return -1;
    if (fp->is_bgzf && fp->fp.bgzf) return bgzf_tell(fp->fp.bgzf);
    if (fp->fp.hfile) return htell(fp->fp.hfile);
    return -1;
}

static int hts_seek_off(htsFile* fp, int64_t off) {
    if (!fp || off < 0) return -1;
    if (fp->is_bgzf && fp->fp.bgzf) return bgzf_seek(fp->fp.bgzf, off, SEEK_SET);
    if (fp->fp.hfile) return hseek(fp->fp.hfile, off, SEEK_SET) < 0 ? -1 : 0;
    return -1;
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

static void load_index(VcfSession* s) {
    s->has_index = false;
    s->tbx = nullptr;
    s->csi = nullptr;
    if (s->is_bcf) {
        s->csi = bcf_index_load(s->path.c_str());
        if (s->csi) s->has_index = true;
    } else {
        s->tbx = tbx_index_load(s->path.c_str());
        if (s->tbx) s->has_index = true;
    }
}

VcfSession* vcf_session_open(const Options& opt) {
    auto* s = new VcfSession;
    s->path = opt.vcf;
    s->fp = bcf_open(opt.vcf.c_str(), "r");
    if (!s->fp) {
        delete s;
        die("cannot open VCF/BCF: " + opt.vcf);
    }
    set_io_threads(s->fp, opt.threads);
    s->hdr = bcf_hdr_read(s->fp);
    if (!s->hdr) {
        bcf_close(s->fp);
        delete s;
        die("cannot read VCF header: " + opt.vcf);
    }
    s->data_offset = hts_tell_off(s->fp);
    s->is_bcf = (s->fp->format.format == bcf);

    int n = bcf_hdr_nsamples(s->hdr);
    s->info.samples.reserve(static_cast<size_t>(n));
    for (int i = 0; i < n; ++i) s->info.samples.emplace_back(s->hdr->samples[i]);

    int nctg = s->hdr->n[BCF_DT_CTG];
    s->info.contigs.reserve(static_cast<size_t>(nctg));
    for (int i = 0; i < nctg; ++i) {
        const char* name = s->hdr->id[BCF_DT_CTG][i].key;
        if (name) s->info.contigs.emplace_back(name);
    }

    load_index(s);
    if (s->has_index) {
        log_info(opt, "VCF index loaded (shared for all regions/contigs)");
    } else {
        log_warn(opt,
                 "no CSI/TBI index for " + opt.vcf +
                     "; sequential scan is much slower. "
                     "Build an index to speed up region/contig loads: "
                     "bcftools index -t " +
                     opt.vcf + "  (or -c for CSI)");
    }
    return s;
}

void vcf_session_close(VcfSession* s) {
    if (!s) return;
    if (s->tbx) tbx_destroy(s->tbx);
    if (s->csi) hts_idx_destroy(s->csi);
    if (s->hdr) bcf_hdr_destroy(s->hdr);
    if (s->fp) bcf_close(s->fp);
    delete s;
}

const VcfHeaderInfo& vcf_session_info(const VcfSession* s) {
    if (!s) die("vcf_session_info: null session");
    return s->info;
}

bool vcf_session_has_index(const VcfSession* s) {
    return s && s->has_index;
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

static bool diploid_alleles(const int32_t* gt, int nsmpl, int si, int& aa, int& bb,
                            bool& multi) {
    if (si < 0 || si >= nsmpl) return false;
    const int32_t a0 = gt[si * 2];
    const int32_t a1 = gt[si * 2 + 1];
    if (bcf_gt_is_missing(a0) || bcf_gt_is_missing(a1)) return false;
    aa = bcf_gt_allele(a0);
    bb = bcf_gt_allele(a1);
    if (aa < 0 || bb < 0) return false;
    if (aa > 1 || bb > 1) multi = true;
    return true;
}

static void count_pop(const int32_t* gt, int nsmpl, const std::vector<int>& idx,
                      int& alt, int& ncall, bool& multi) {
    alt = 0;
    ncall = 0;
    multi = false;
    for (int si : idx) {
        int aa = 0, bb = 0;
        if (!diploid_alleles(gt, nsmpl, si, aa, bb, multi)) continue;
        if (aa <= 1 && bb <= 1) {
            alt += (aa == 1) + (bb == 1);
            ncall += 2;
        }
    }
}

// popB: counts + dosage in one pass. Dosage: missing / multi → 0 (hardingnj).
static void count_and_dosage_pop(const int32_t* gt, int nsmpl,
                                 const std::vector<int>& idx, int& alt, int& ncall,
                                 bool& multi, std::vector<int8_t>& dosage) {
    alt = 0;
    ncall = 0;
    multi = false;
    dosage.resize(idx.size());
    for (size_t i = 0; i < idx.size(); ++i) {
        const int si = idx[i];
        int aa = 0, bb = 0;
        bool m = false;
        if (!diploid_alleles(gt, nsmpl, si, aa, bb, m)) {
            dosage[i] = 0;
            continue;
        }
        if (m) multi = true;
        if (aa <= 1 && bb <= 1) {
            const int8_t d = static_cast<int8_t>((aa == 1) + (bb == 1));
            alt += d;
            ncall += 2;
            dosage[i] = d;
        } else {
            dosage[i] = 0;
        }
    }
}

std::vector<SnpData> load_snps(VcfSession* s, const Options& opt,
                               const SamplePlan& plan, const RegionTarget& target) {
    if (!s || !s->fp || !s->hdr) die("load_snps: null VCF session");
    if (target.chrom.empty()) die("load_snps: empty contig");

    htsFile* fp = s->fp;
    bcf_hdr_t* hdr = s->hdr;

    int tid = bcf_hdr_name2id(hdr, target.chrom.c_str());
    if (tid < 0) die("chromosome not found in VCF header: " + target.chrom);

    int64_t pos_lo = 1, pos_hi = 0;
    bool open_end = true;
    load_bounds(target, opt.size, pos_lo, pos_hi, open_end);
    const std::string region = hts_region_query(target, opt.size);

    hts_itr_t* itr = nullptr;
    if (s->has_index) {
        if (s->is_bcf && s->csi) {
            itr = bcf_itr_querys(s->csi, hdr, region.c_str());
            if (!itr) itr = bcf_itr_querys(s->csi, hdr, target.chrom.c_str());
        } else if (s->tbx) {
            itr = tbx_itr_querys(s->tbx, region.c_str());
            if (!itr) itr = tbx_itr_querys(s->tbx, target.chrom.c_str());
        }
    }

    if (itr) {
        log_info(opt, "Indexed load: " + region);
    } else {
        // Sequential scan of whole file, filter by tid/pos (shared session).
        if (s->data_offset >= 0 && hts_seek_off(fp, s->data_offset) != 0)
            die("cannot seek VCF to start of records for sequential scan: " + s->path);
        log_info(opt, "Sequential load: " + target.chrom + " (" + region + ")");
    }

    bcf1_t* rec = bcf_init();
    int32_t* gt_arr = nullptr;
    int ngt_arr = 0;

    std::vector<SnpData> snps;
    snps.reserve(kSnpReserveHint);
    std::vector<int8_t> dosage_b;

    int64_t n_total = 0, n_multi = 0, n_missing_pop = 0, n_fixed_p2 = 0,
            n_kept = 0, n_parse_fail = 0;

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
        if (multi_a) {
            ++n_multi;
            return;
        }
        if (n_a == 0) {
            ++n_missing_pop;
            return;
        }
        count_and_dosage_pop(gt_arr, nsmpl, plan.idx_b, alt_b, n_b, multi_b, dosage_b);
        if (multi_b) {
            ++n_multi;
            return;
        }
        if (n_b == 0) {
            ++n_missing_pop;
            return;
        }

        const int ref_b = n_b - alt_b;
        if (alt_b == 0 || ref_b == 0 || alt_b == 1 || ref_b == 1) {
            ++n_fixed_p2;
            return;
        }

        SnpData sn;
        sn.pos = pos1;
        sn.x_alt = alt_a;
        sn.n_a = n_a;
        sn.n_b = n_b;
        sn.q2 = static_cast<double>(alt_b) / static_cast<double>(n_b);
        sn.dosage_b = std::move(dosage_b);
        snps.push_back(std::move(sn));
        ++n_kept;
    };

    if (itr) {
        if (s->tbx && !s->is_bcf) {
            kstring_t str = {0, 0, nullptr};
            while (tbx_itr_next(fp, s->tbx, itr, &str) >= 0) {
                if (vcf_parse(&str, hdr, rec) < 0) {
                    ++n_parse_fail;
                    continue;
                }
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

    log_info(opt, std::to_string(n_total) + " records considered on " +
                      target.chrom + " (load " + region + ")");
    if (n_parse_fail > 0)
        log_warn(opt, std::to_string(n_parse_fail) + " VCF lines failed parse (skipped)");
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
    }
    return snps;
}

}  // namespace xpclr
