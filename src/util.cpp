#include "xpclr.hpp"

#include <cstdlib>
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
    if (opt.verbose >= 0) std::cerr << "[W::xpclr] " << msg << "\n";
}

void print_usage(const char* argv0) {
    std::cerr
        << "Usage:\n"
        << "  " << argv0
        << " -i <vcf.gz> --pop <pop.txt> -a <popA> -b <popB> -c <chr> -o <out.tsv>\n"
        << "             [options]\n"
        << "\n"
        << "Cross-population composite likelihood ratio (XP-CLR) scan.\n"
        << "C++/htslib rewrite aligned with hardingnj/xpclr (Chen et al. 2010).\n"
        << "\n"
        << "Required:\n"
        << "  -i, --input FILE     VCF/BCF (optionally bgzipped; CSI/TBI recommended)\n"
        << "  --pop FILE           Population file: two columns SAMPLE  GROUP\n"
        << "  -a, --popA NAME      Name of population A (selected; target for selection)\n"
        << "  -b, --popB NAME      Name of population B (reference / non-selected)\n"
        << "  -c, --chr NAME       Contig/chromosome to scan\n"
        << "  -o, --out FILE       Output TSV path\n"
        << "\n"
        << "Options:\n"
        << "  --rrate FLOAT        Recombination rate per bp (default 1e-8)\n"
        << "  --ld FLOAT           LD r^2 cutoff for SNP weights (default 0.95)\n"
        << "  --maxsnps INT        Max SNPs per window (default 200)\n"
        << "  --minsnps INT        Min SNPs per window (default 10)\n"
        << "  --size INT           Window size bp (default 20000)\n"
        << "  --step INT           Window step bp (default 20000)\n"
        << "  --start INT          First window start (default 1)\n"
        << "  --stop INT           Last base; 0 = last SNP position (default 0)\n"
        << "  --threads INT        OpenMP threads for windows (default 1)\n"
        << "  --seed INT           RNG seed for maxsnps subsample (default 1)\n"
        << "  --phased             Use haplotype-style dosage (reserved; default off)\n"
        << "  --no-early-stop      Disable selection-grid early stop (issue #115)\n"
        << "  -V, --verbose INT    0=quiet, 1=info, 2=debug (default 1)\n"
        << "  -h, --help           Show this help and exit 0\n"
        << "  -v, --version        Show version and exit 0\n"
        << "\n"
        << "Input:\n"
        << "  VCF with GT. Sample IDs matched by name to --pop column 1.\n"
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
        << " -i demo/smoke.vcf.gz --pop demo/pop.txt -a popA -b popB -c 1 -o out.tsv\n"
        << "  " << argv0
        << " -i snps.vcf.gz --pop pops.txt -a landrace -b wild -c 5 -o chr5.tsv \\\n"
        << "      --size 500000 --step 100000 --minsnps 2 --threads 8\n";
}

static bool eq(const char* a, const char* b) { return std::string(a) == b; }

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
            std::cout << "xpclr-cpp 0.1.0 (htslib; aligned hardingnj/xpclr)\n";
            std::exit(0);
        }
        if (a == "-i" || a == "--input")
            opt.vcf = need("-i");
        else if (a == "--pop")
            opt.pop_file = need("--pop");
        else if (a == "-a" || a == "--popA")
            opt.pop_a = need("-a");
        else if (a == "-b" || a == "--popB")
            opt.pop_b = need("-b");
        else if (a == "-c" || a == "--chr")
            opt.chrom = need("-c");
        else if (a == "-o" || a == "--out")
            opt.out = need("-o");
        else if (a == "--rrate")
            opt.rrate = std::stod(need("--rrate"));
        else if (a == "--ld")
            opt.ldcutoff = std::stod(need("--ld"));
        else if (a == "--maxsnps")
            opt.maxsnps = std::stoi(need("--maxsnps"));
        else if (a == "--minsnps")
            opt.minsnps = std::stoi(need("--minsnps"));
        else if (a == "--size")
            opt.size = std::stoll(need("--size"));
        else if (a == "--step")
            opt.step = std::stoll(need("--step"));
        else if (a == "--start")
            opt.start = std::stoll(need("--start"));
        else if (a == "--stop")
            opt.stop = std::stoll(need("--stop"));
        else if (a == "--threads")
            opt.threads = std::stoi(need("--threads"));
        else if (a == "--seed")
            opt.seed = static_cast<uint64_t>(std::stoull(need("--seed")));
        else if (a == "--phased")
            opt.phased = true;
        else if (a == "--no-early-stop")
            opt.early_stop = false;
        else if (a == "-V" || a == "--verbose")
            opt.verbose = std::stoi(need("-V"));
        else
            die("unknown option: " + a);
    }
    if (opt.vcf.empty()) die("required: -i/--input");
    if (opt.pop_file.empty()) die("required: --pop");
    if (opt.pop_a.empty()) die("required: -a/--popA");
    if (opt.pop_b.empty()) die("required: -b/--popB");
    if (opt.chrom.empty()) die("required: -c/--chr");
    if (opt.out.empty()) die("required: -o/--out");
    if (opt.pop_a == opt.pop_b) die("-a and -b must name different populations");
    if (opt.minsnps < 2) die("--minsnps must be >= 2");
    if (opt.maxsnps < opt.minsnps) die("--maxsnps must be >= --minsnps");
    if (opt.threads < 1) die("--threads must be >= 1");
    if (opt.size < 1 || opt.step < 1) die("--size/--step must be >= 1");
    return opt;
}

}  // namespace xpclr
