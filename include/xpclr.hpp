#pragma once
// XP-CLR C++/htslib rewrite aligned with hardingnj/xpclr (Chen et al. 2010 logic).

#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

namespace xpclr {

inline constexpr const char* kVersion = "0.3.0";

// LD weight computation mode for determine_weights.
//   0 dosage-fill  : unphased dosages; missing genotype -> 0 (hardingnj/python)
//   1 phased       : phased haplotypes; single-haplotype Pearson r (raw -p1)
//   2 em           : unphased genotypes; EM-inferred haplotype freq Pearson r (raw -p0)
//   3 pairwise     : unphased dosages; skip missing samples per SNP pair (xpclrs)
enum class LdMode {
    dosage_fill,  // -p 0
    phased,       // -p 1
    em,           // -p 2
    pairwise,     // -p 3 (default)
};

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
    std::string gmap_path;  // optional CHROM POS GDIST; empty => POS*rrate
    double rrate = 1e-8;
    double ldcutoff = 0.95;
    double ne = 20000.0;  // effective population size; enters c = 1 - exp(-ln(2*ne)*r/s)
    int maxsnps = 500;
    int minsnps = 10;
    int64_t size = 50000;
    int64_t step = 25000;
    int threads = 1;
    uint64_t seed = 1;
    // Drop highest fraction of per-SNP r=(q1-q2)^2/(q2(1-q2)) before MoM mean.
    // 0 = classic hardingnj raw mean; default 0.01 resists sweep/outlier inflation.
    double omega_trim = 0.01;
    // false = full s-grid max (default). true = unimodal early exit (python-like).
    bool unimodal_s = false;
    // LD weight computation mode (see LdMode).
    // Default pairwise: skip missing per SNP pair (no HWE assumption, robust).
    // --phased 0=dosage-fill, 1=phased, 2=EM, 3=pairwise-complete.
    LdMode ld_mode = LdMode::pairwise;
    // popB genotypes are stored as haplotypes (2*n_b columns, 0/1) iff phased.
    bool phased_input() const { return ld_mode == LdMode::phased; }
    int verbose = 1;
};

struct PopAssignment {
    std::unordered_map<std::string, std::string> sample_to_group;
    std::unordered_map<std::string, std::vector<std::string>> group_samples;
    std::vector<std::string> group_order;
};

struct SamplePlan {
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
};

// Contig/region SNP table: scalar fields + one packed popB dosage matrix
// (row-major, n_snps * n_b). Avoids per-SNP vector heap overhead.
struct SnpSet {
    std::vector<SnpData> snps;
    int n_b = 0;  // dosage columns (= matched popB samples)
    std::vector<int8_t> dosage_b;  // size snps.size()*n_b

    bool empty() const { return snps.empty(); }
    size_t size() const { return snps.size(); }
    const int8_t* dosage_row(size_t i) const {
        return dosage_b.data() + i * static_cast<size_t>(n_b);
    }
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
[[noreturn]] void die(const std::string& msg);

// Parse -r: "", "Chr01", "Chr01:200-30000", "Chr01:200-", "Chr01:200"
RegionTarget parse_region_string(const std::string& reg);

PopAssignment load_pop_file(const std::string& path, const Options& opt);
SamplePlan resolve_samples(const std::vector<std::string>& vcf_samples,
                           const PopAssignment& pop,
                           const Options& opt);

// Long-lived VCF handle: one open + header + index for all contigs/regions.
// Opaque internals live in vcf_io.cpp (htslib types).
struct VcfSession;

struct VcfHeaderInfo {
    std::vector<std::string> samples;
    std::vector<std::string> contigs;
};

// Open VCF/BCF, read samples+contigs, load CSI/TBI once if present.
// Warns once when no index (scanning is much slower).
VcfSession* vcf_session_open(const Options& opt);
void vcf_session_close(VcfSession* s);
const VcfHeaderInfo& vcf_session_info(const VcfSession* s);
bool vcf_session_has_index(const VcfSession* s);

SnpSet load_snps(VcfSession* s, const Options& opt,
                 const SamplePlan& plan, const RegionTarget& target);

// win_start/win_stop: window grid (stop 0 => last SNP pos on loaded set).
std::vector<WindowResult> xpclr_scan(const SnpSet& snps,
                                     const Options& opt,
                                     const std::string& chrom,
                                     int64_t win_start, int64_t win_stop);

void write_results(const std::string& path, const std::vector<WindowResult>& rows);
int run_xpclr(const Options& opt);

double determine_c(double r, double s, double ne, double min_rd = 1e-7,
                   int sf = 5);
double chen_likelihood(int xj, int nj, double c, double p2, double var);
// trim in [0,1): fraction of highest r_i dropped before mean (0 = no trim).
double estimate_omega(const std::vector<SnpData>& snps, double trim = 0.01);
// overload for packed set
double estimate_omega(const SnpSet& snps, double trim = 0.01);

// Internal types shared between clr_stats.cpp and scan.cpp.
namespace detail {

struct WinRow {
    // Field names mirror SnpData (x_alt/n_a/q2) so scan->math needs no translation.
    int x_alt;
    int n_a;
    double rd;     // |genetic distance to window center|
    double q2;
    double omega;
    double weight;
};

extern const double kSelCoefs[];
extern const int kNSel;

std::vector<double> determine_weights(
    const SnpSet& snps, const std::vector<int>& ix, double ldcutoff, LdMode mode);

double calculate_cl(double sc, const std::vector<WinRow>& dat, double ne);

void compute_xpclr(const std::vector<WinRow>& dat, double ne,
                   bool unimodal_s, double& modelL, double& nullL, double& sel);

}  // namespace detail

}  // namespace xpclr
