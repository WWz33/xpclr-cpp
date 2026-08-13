#include "xpclr.hpp"

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iostream>

namespace xpclr {

void die(const std::string& msg) {
    std::cerr << "[E::xpclr] " << msg << "\n";
    std::cerr << "Use -h for help.\n";
    std::exit(1);
}

void log_info(const Options& opt, const std::string& msg) {
    if (opt.verbose >= 1) std::cerr << "[I::xpclr] " << msg << "\n";
}

void log_warn(const Options& opt, const std::string& msg) {
    (void)opt;  // warnings always printed
    std::cerr << "[W::xpclr] " << msg << "\n";
}

void print_usage(const char* argv0) {
    std::cerr
        << "Usage:\n"
        << "  " << argv0
        << " -i <vcf.gz> -p <pop.txt> -a <popA> -b <popB> -o <out.tsv>\n"
        << "             [-r <region>] [options]\n"
        << "\n"
        << "Cross-population composite likelihood ratio (XP-CLR) scan.\n"
        << "C++/htslib rewrite aligned with hardingnj/xpclr (Chen et al. 2010).\n"
        << "\n"
        << "Required:\n"
        << "  -i, --input FILE     VCF/BCF (optionally bgzipped; CSI/TBI recommended)\n"
        << "  -p, --pop FILE       Population file: two columns SAMPLE  GROUP\n"
        << "  -a, --popA NAME      Name of population A (selected; target for selection)\n"
        << "  -b, --popB NAME      Name of population B (reference / non-selected)\n"
        << "  -o, --out FILE       Output TSV path\n"
        << "\n"
        << "Options:\n"
        << "  -r, --regions REG    Contig or interval (htslib style), e.g. Chr01,\n"
        << "                       Chr01:200-30000, 1:1000000-; omit = all contigs\n"
        << "  -R, --rrate FLOAT    Recombination rate per bp (default 1e-8)\n"
        << "  -L, --ld FLOAT       LD r^2 cutoff for SNP weights (default 0.95)\n"
        << "  -N, --ne FLOAT       Effective population size (default 20000)\n"
        << "  -k, --maxsnps INT    Max SNPs per window (default 500)\n"
        << "  -m, --minsnps INT    Min SNPs per window (default 10)\n"
        << "  -w, --size INT       Window size bp (default 50000)\n"
        << "  -s, --step INT       Window step bp (default 25000)\n"
        << "  -t, --threads INT    Threads (default 1)\n"
        << "  --seed INT           RNG seed for maxsnps subsample (default 1)\n"
        << "  --omega-trim FLOAT   Drop top fraction of SNP r when estimating omega\n"
        << "                       (default 0.01; 0 = hardingnj raw mean)\n"
        << "  -P, --phased INT     LD weight mode: 0=dosage fill (missing->0,\n"
        << "                       matches python; biased with missing data),\n"
        << "                       1=phased haplotype r (raw -p1; needs phased VCF),\n"
        << "                       2=EM two-locus phase inference r (raw -p0),\n"
        << "                       3=pairwise-complete dosage r (default; skip missing)\n"
        << "  --unimodal-s         Stop at first LL decline along s (hardingnj/python)\n"
        << "  -G, --gmap FILE      Genetic map CHROM POS GDIST\n"
        << "  -V, --verbose INT    0=quiet, 1=info, 2=debug (default 1)\n"
        << "  -h, --help           Show this help and exit 0\n"
        << "  -v, --version        Show version and exit 0\n"
        << "\n"
        << "Input:\n"
        << "  VCF with GT. Sample IDs matched by name to -p/--pop column 1.\n"
        << "  pop.txt: SAMPLE  GROUP (whitespace). Lines starting with # ignored.\n"
        << "  Duplicate SAMPLE with same GROUP: warn, keep once.\n"
        << "  Duplicate SAMPLE with different GROUP: error.\n"
        << "  Selected group with 0 samples in VCF: error.\n"
        << "\n"
        << "Output (TSV columns):\n"
        << "  id chrom start stop pos_start pos_stop modelL nullL sel_coef\n"
        << "  nSNPs nSNPs_avail xpclr xpclr_norm\n"
        << "\n"
        << "Examples:\n"
        << "  " << argv0
        << " -i data/smoke.vcf.gz -p data/pop_smoke.txt -a popA -b popB -o out.tsv\n"
        << "  " << argv0
        << " -i snps.vcf.gz -p pops.txt -a landrace -b wild -r Chr01 -o chr1.tsv \\\n"
        << "      -w 500000 -s 100000 -m 2 -t 8\n"
        << "  " << argv0
        << " -i snps.vcf.gz -p pops.txt -a W -b C -r Chr01:200-30000 -o sub.tsv\n";
}

static bool eq(const char* a, const char* b) { return std::strcmp(a, b) == 0; }

RegionTarget parse_region_string(const std::string& reg) {
    RegionTarget t;
    if (reg.empty()) return t;

    auto colon = reg.find(':');
    if (colon == std::string::npos) {
        t.chrom = reg;
        return t;
    }
    t.chrom = reg.substr(0, colon);
    if (t.chrom.empty()) die("invalid -r/--regions: empty contig in '" + reg + "'");

    std::string coords = reg.substr(colon + 1);
    if (coords.empty()) return t;

    auto dash = coords.find('-');
    if (dash == std::string::npos) {
        t.beg = std::stoll(coords);
        t.has_beg = true;
        if (t.beg < 1) die("invalid -r region start (<1): " + reg);
        return t;
    }

    std::string left = coords.substr(0, dash);
    std::string right = coords.substr(dash + 1);
    if (!left.empty()) {
        t.beg = std::stoll(left);
        t.has_beg = true;
        if (t.beg < 1) die("invalid -r region start (<1): " + reg);
    }
    if (!right.empty()) {
        t.end = std::stoll(right);
        t.has_end = true;
        if (t.end < 1) die("invalid -r region end (<1): " + reg);
        if (t.has_beg && t.end < t.beg)
            die("invalid -r region end < start: " + reg);
    }
    return t;
}

Options parse_args(int argc, char** argv) {
    Options opt;
    if (argc == 1) {
        print_usage(argv[0]);
        std::exit(0);
    }
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        auto need = [&](const char* name) -> std::string {
            if (i + 1 >= argc) die(std::string("missing value for ") + name);
            return argv[++i];
        };
        if (eq(argv[i], "-h") || eq(argv[i], "--help")) {
            print_usage(argv[0]);
            std::exit(0);
        }
        if (eq(argv[i], "-v") || eq(argv[i], "--version")) {
            std::cout << "xpclr-cpp " << kVersion << " (htslib; aligned hardingnj/xpclr)\n";
            std::exit(0);
        }
        if (a == "-i" || a == "--input")
            opt.vcf = need("-i");
        else if (a == "-p" || a == "--pop")
            opt.pop_file = need("-p/--pop");
        else if (a == "-a" || a == "--popA")
            opt.pop_a = need("-a");
        else if (a == "-b" || a == "--popB")
            opt.pop_b = need("-b");
        else if (a == "-r" || a == "--regions" || a == "--region")
            opt.region = need("-r");
        else if (a == "-c" || a == "--chr")
            die("removed: use -r/--regions (e.g. -r Chr01 or -r Chr01:200-30000)");
        else if (a == "--start" || a == "--stop")
            die("removed: set genomic span via -r/--regions "
                "(e.g. -r Chr01:1-5000000); --size/--step still set the window grid");
        else if (a == "-o" || a == "--out")
            opt.out = need("-o");
        else if (a == "-R" || a == "--rrate")
            opt.rrate = std::stod(need("--rrate"));
        else if (a == "-G" || a == "--gmap")
            opt.gmap_path = need("--gmap");
        else if (a == "-L" || a == "--ld")
            opt.ldcutoff = std::stod(need("--ld"));
        else if (a == "-N" || a == "--ne")
            opt.ne = std::stod(need("--ne"));
        else if (a == "-k" || a == "--maxsnps")
            opt.maxsnps = std::stoi(need("--maxsnps"));
        else if (a == "-m" || a == "--minsnps")
            opt.minsnps = std::stoi(need("--minsnps"));
        else if (a == "-w" || a == "--size")
            opt.size = std::stoll(need("--size"));
        else if (a == "-s" || a == "--step")
            opt.step = std::stoll(need("--step"));
        else if (a == "-t" || a == "--threads")
            opt.threads = std::stoi(need("--threads"));
        else if (a == "--seed")
            opt.seed = static_cast<uint64_t>(std::stoull(need("--seed")));
        else if (a == "-P" || a == "--phased") {
            int p = std::stoi(need("--phased"));
            if (p == 0) { opt.ld_mode = LdMode::dosage_fill; opt.phased_input = false; }
            else if (p == 1) { opt.ld_mode = LdMode::phased; opt.phased_input = true; }
            else if (p == 2) { opt.ld_mode = LdMode::em; opt.phased_input = false; }
            else if (p == 3) { opt.ld_mode = LdMode::pairwise; opt.phased_input = false; }
            else die("--phased must be 0 (dosage-fill), 1 (phased), 2 (EM), or 3 (pairwise)");
        }
        else if (a == "--omega-trim")
            opt.omega_trim = std::stod(need("--omega-trim"));
        else if (a == "--unimodal-s")
            opt.unimodal_s = true;
        else if (a == "-V" || a == "--verbose")
            opt.verbose = std::stoi(need("-V"));
        else if (a == "--no-early-stop")
            die("removed: full s-grid max is now default; drop this flag "
                "(or use --unimodal-s for hardingnj/python-style early exit)");
        else if (a == "--early-stop")
            die("renamed: use --unimodal-s (assume unimodal likelihood along s grid)");
        else
            die("unknown option: " + a);
    }
    if (opt.vcf.empty()) die("required: -i/--input");
    if (opt.pop_file.empty()) die("required: -p/--pop");
    if (opt.pop_a.empty()) die("required: -a/--popA");
    if (opt.pop_b.empty()) die("required: -b/--popB");
    if (opt.out.empty()) die("required: -o/--out");
    if (opt.pop_a == opt.pop_b) die("-a and -b must name different populations");
    if (opt.minsnps < 2) die("--minsnps must be >= 2");
    if (opt.maxsnps < opt.minsnps) die("--maxsnps must be >= --minsnps");
    if (opt.threads < 1) die("--threads must be >= 1");
    if (opt.size < 1 || opt.step < 1) die("--size/--step must be >= 1");
    if (!(opt.omega_trim >= 0.0 && opt.omega_trim < 1.0) || !std::isfinite(opt.omega_trim))
        die("--omega-trim must be in [0,1)");
    if (!(opt.ne > 0.0) || !std::isfinite(opt.ne))
        die("--ne must be > 0");
    if (!opt.region.empty()) (void)parse_region_string(opt.region);
    return opt;
}

}  // namespace xpclr
