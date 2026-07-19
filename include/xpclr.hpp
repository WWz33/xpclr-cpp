#pragma once
// XP-CLR C++/htslib rewrite aligned with hardingnj/xpclr (Chen et al. 2010 logic).
// Copyright: local rewrite for IO/parallel; algorithm parity with Python reference.

#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

namespace xpclr {

struct Options {
    std::string vcf;
    std::string pop_file;
    std::string pop_a;
    std::string pop_b;
    std::string chrom;
    std::string out;
    double rrate = 1e-8;
    double ldcutoff = 0.95;
    int maxsnps = 200;
    int minsnps = 10;
    int64_t size = 20000;
    int64_t start = 1;
    int64_t stop = 0;  // 0 => last SNP pos
    int64_t step = 20000;
    int threads = 1;
    uint64_t seed = 1;
    bool phased = false;
    bool early_stop = true;  // match Python; disable after issue #115 review
    int verbose = 1;         // 0 quiet, 1 info, 2 debug
};

struct PopAssignment {
    // sample -> group from file (first-seen; duplicates handled at load)
    std::unordered_map<std::string, std::string> sample_to_group;
    // group -> ordered unique samples as listed in file
    std::unordered_map<std::string, std::vector<std::string>> group_samples;
    std::vector<std::string> group_order;
};

struct SamplePlan {
    std::vector<std::string> vcf_samples;  // all VCF sample names
    std::vector<int> idx_a;                // VCF column indices for pop A
    std::vector<int> idx_b;
    int n_input_a = 0;
    int n_input_b = 0;
    int n_matched_a = 0;
    int n_matched_b = 0;
};

// Per-SNP data kept after filters (popA/popB counts + popB dosage for LD).
struct SnpData {
    int64_t pos = 0;
    int x_alt = 0;   // popA alt allele count
    int n_a = 0;     // popA called alleles
    int n_b = 0;     // popB called alleles
    double q2 = 0.0; // popB alt frequency
    // popB dosage (n_alt per diploid sample); missing filled as 0 (Python fill=0)
    std::vector<int8_t> dosage_b;
};

struct WindowResult {
    int64_t start = 0;
    int64_t stop = 0;
    int64_t pos_start = 0;
    int64_t pos_stop = 0;
    double modelL = 0.0;
    double nullL = 0.0;
    double sel_coef = 0.0;
    int nSNPs = 0;
    int nSNPs_avail = 0;
    bool valid = false;
};

// ---- CLI / logging ----
void print_usage(const char* argv0);
Options parse_args(int argc, char** argv);
void log_info(const Options& opt, const std::string& msg);
void log_warn(const Options& opt, const std::string& msg);
void die(const std::string& msg);  // stderr + exit(1)

// ---- pop file ----
PopAssignment load_pop_file(const std::string& path, const Options& opt);
SamplePlan resolve_samples(const std::vector<std::string>& vcf_samples,
                           const PopAssignment& pop,
                           const Options& opt);

// ---- VCF IO (htslib) ----
std::vector<std::string> read_vcf_samples(const std::string& path);
// Load biallelic SNPs on chrom for planned samples; apply hardingnj filters.
std::vector<SnpData> load_snps(const Options& opt, const SamplePlan& plan);

// ---- core scan ----
std::vector<WindowResult> xpclr_scan(const std::vector<SnpData>& snps,
                                     const Options& opt);

void write_results(const std::string& path,
                   const std::string& chrom,
                   const std::vector<WindowResult>& rows);

int run_xpclr(const Options& opt);

// likelihood helpers (exposed for tests)
double determine_c(double r, double s, double ne = 20000.0, double min_rd = 1e-7,
                   int sf = 5);
double chen_likelihood(int xj, int nj, double c, double p2, double var);
double estimate_omega(const std::vector<SnpData>& snps);

}  // namespace xpclr
