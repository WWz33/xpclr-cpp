#pragma once
// XP-CLR C++/htslib rewrite aligned with hardingnj/xpclr (Chen et al. 2010 logic).

#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

namespace xpclr {

// Parsed -r/--regions (htslib/bcftools style, 1-based inclusive coords).
struct RegionTarget {
    std::string chrom;   // contig name
    int64_t beg = 0;     // 0 => from first base / first SNP
    int64_t end = 0;     // 0 => open end (to last SNP on contig)
    bool has_beg = false;
    bool has_end = false;
};

struct Options {
    std::string vcf;
    std::string pop_file;
    std::string pop_a;
    std::string pop_b;
    std::string region;  // raw -r string; empty = all contigs
    std::string out;
    double rrate = 1e-8;
    double ldcutoff = 0.95;
    int maxsnps = 200;
    int minsnps = 10;
    int64_t size = 20000;
    int64_t step = 20000;
    int threads = 1;
    uint64_t seed = 1;
    bool phased = false;
    // false = full s-grid max (default). true = unimodal early exit (python-like).
    bool unimodal_s = false;
    int verbose = 1;
};

struct PopAssignment {
    std::unordered_map<std::string, std::string> sample_to_group;
    std::unordered_map<std::string, std::vector<std::string>> group_samples;
    std::vector<std::string> group_order;
};

struct SamplePlan {
    std::vector<std::string> vcf_samples;
    std::vector<int> idx_a;
    std::vector<int> idx_b;
    int n_input_a = 0;
    int n_input_b = 0;
    int n_matched_a = 0;
    int n_matched_b = 0;
};

struct SnpData {
    int64_t pos = 0;
    int x_alt = 0;
    int n_a = 0;
    int n_b = 0;
    double q2 = 0.0;
    std::vector<int8_t> dosage_b;
};

struct WindowResult {
    std::string chrom;
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

void print_usage(const char* argv0);
Options parse_args(int argc, char** argv);
void log_info(const Options& opt, const std::string& msg);
void log_warn(const Options& opt, const std::string& msg);
void die(const std::string& msg);

// Parse -r: "", "Chr01", "Chr01:200-30000", "Chr01:200-", "Chr01:200"
RegionTarget parse_region_string(const std::string& reg);

PopAssignment load_pop_file(const std::string& path, const Options& opt);
SamplePlan resolve_samples(const std::vector<std::string>& vcf_samples,
                           const PopAssignment& pop,
                           const Options& opt);

std::vector<std::string> read_vcf_samples(const std::string& path);
std::vector<std::string> read_vcf_contigs(const std::string& path);
std::vector<SnpData> load_snps(const Options& opt, const SamplePlan& plan,
                               const RegionTarget& target);

// win_start/win_stop: window grid (stop 0 => last SNP pos on loaded set).
std::vector<WindowResult> xpclr_scan(const std::vector<SnpData>& snps,
                                     const Options& opt,
                                     const std::string& chrom,
                                     int64_t win_start, int64_t win_stop);

void write_results(const std::string& path, const std::vector<WindowResult>& rows);
int run_xpclr(const Options& opt);

double determine_c(double r, double s, double ne = 20000.0, double min_rd = 1e-7,
                   int sf = 5);
double chen_likelihood(int xj, int nj, double c, double p2, double var);
double estimate_omega(const std::vector<SnpData>& snps);

}  // namespace xpclr
