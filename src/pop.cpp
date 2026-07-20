#include "xpclr.hpp"

#include <fstream>
#include <set>
#include <sstream>

namespace xpclr {

static std::vector<std::string> split_ws(const std::string& line) {
    std::istringstream iss(line);
    std::vector<std::string> tok;
    std::string t;
    while (iss >> t) tok.push_back(t);
    return tok;
}

PopAssignment load_pop_file(const std::string& path, const Options& opt) {
    std::ifstream in(path);
    if (!in) die("cannot open pop file: " + path);

    PopAssignment pa;
    std::unordered_map<std::string, int> input_count;  // per group raw lines
    int nline = 0;
    std::string line;
    while (std::getline(in, line)) {
        ++nline;
        if (line.empty()) continue;
        if (line[0] == '#') continue;
        // strip CR
        if (!line.empty() && line.back() == '\r') line.pop_back();
        auto tok = split_ws(line);
        if (tok.size() < 2) {
            die("pop file line " + std::to_string(nline) +
                ": need SAMPLE GROUP, got: " + line);
        }
        const std::string& sample = tok[0];
        const std::string& group = tok[1];

        auto it = pa.sample_to_group.find(sample);
        if (it != pa.sample_to_group.end()) {
            if (it->second == group) {
                log_warn(opt, "duplicate sample '" + sample + "' with same group '" +
                                  group + "' (line " + std::to_string(nline) +
                                  "); keeping once");
                input_count[group] += 1;  // still count as input occurrence
                continue;
            }
            die("sample '" + sample + "' assigned to conflicting groups: '" +
                it->second + "' vs '" + group + "' (line " + std::to_string(nline) +
                ")");
        }
        pa.sample_to_group[sample] = group;
        if (pa.group_samples.find(group) == pa.group_samples.end()) {
            pa.group_order.push_back(group);
        }
        pa.group_samples[group].push_back(sample);
        input_count[group] += 1;
    }
    if (pa.sample_to_group.empty()) die("pop file has no sample assignments: " + path);

    log_info(opt, "Reading population file: " + path);
    {
        std::ostringstream oss;
        oss << "  Groups in file:";
        for (const auto& g : pa.group_order) oss << " " << g;
        log_info(opt, oss.str());
    }
    for (const auto& g : pa.group_order) {
        log_info(opt, "  " + g + ": input N=" +
                          std::to_string(input_count[g]) +
                          ", unique samples=" +
                          std::to_string(pa.group_samples[g].size()));
    }
    return pa;
}

SamplePlan resolve_samples(const std::vector<std::string>& vcf_samples,
                           const PopAssignment& pop,
                           const Options& opt) {
    SamplePlan plan;

    std::unordered_map<std::string, int> name2idx;
    name2idx.reserve(vcf_samples.size() * 2);
    for (size_t i = 0; i < vcf_samples.size(); ++i)
        name2idx[vcf_samples[i]] = static_cast<int>(i);

    log_info(opt, "VCF samples: " + std::to_string(vcf_samples.size()));

    auto ga = pop.group_samples.find(opt.pop_a);
    auto gb = pop.group_samples.find(opt.pop_b);
    if (ga == pop.group_samples.end())
        die("population -a '" + opt.pop_a + "' not found in pop file");
    if (gb == pop.group_samples.end())
        die("population -b '" + opt.pop_b + "' not found in pop file");

    plan.n_input_a = static_cast<int>(ga->second.size());
    plan.n_input_b = static_cast<int>(gb->second.size());

    std::vector<std::string> miss_a, miss_b;
    std::set<int> used;
    for (const auto& s : ga->second) {
        auto it = name2idx.find(s);
        if (it == name2idx.end()) {
            miss_a.push_back(s);
            continue;
        }
        if (used.count(it->second)) {
            log_warn(opt, "sample '" + s + "' already selected; skip duplicate column");
            continue;
        }
        plan.idx_a.push_back(it->second);
        used.insert(it->second);
    }
    for (const auto& s : gb->second) {
        auto it = name2idx.find(s);
        if (it == name2idx.end()) {
            miss_b.push_back(s);
            continue;
        }
        // Same sample may appear in A and B (overlap allowed); warn once.
        if (used.count(it->second)) {
            log_warn(opt, "sample '" + s + "' is in both -a and -b (allowed)");
        }
        plan.idx_b.push_back(it->second);
    }

    plan.n_matched_a = static_cast<int>(plan.idx_a.size());
    plan.n_matched_b = static_cast<int>(plan.idx_b.size());

    auto log_pop = [&](const std::string& name, int n_in, int n_matched,
                       const std::vector<std::string>& miss) {
        log_info(opt, "  " + name + ": input N=" + std::to_string(n_in) +
                          ", matched in VCF n=" + std::to_string(n_matched) +
                          ", used in analysis n=" + std::to_string(n_matched));
        if (!miss.empty()) {
            log_warn(opt, "  " + name + ": " + std::to_string(miss.size()) +
                              " sample(s) not in VCF (e.g. " + miss[0] + ")");
        }
    };
    log_pop(opt.pop_a, plan.n_input_a, plan.n_matched_a, miss_a);
    log_pop(opt.pop_b, plan.n_input_b, plan.n_matched_b, miss_b);

    // other groups present
    for (const auto& g : pop.group_order) {
        if (g == opt.pop_a || g == opt.pop_b) continue;
        log_info(opt, "  (other group not selected: " + g + ", input N=" +
                          std::to_string(pop.group_samples.at(g).size()) + ")");
    }

    if (plan.n_matched_a == 0 && plan.n_matched_b == 0) {
        die("no samples matched in VCF for -a '" + opt.pop_a + "' or -b '" +
            opt.pop_b + "'");
    }
    if (plan.n_matched_a == 0) {
        die("population -a '" + opt.pop_a + "': 0 samples matched in VCF");
    }
    if (plan.n_matched_b == 0) {
        die("population -b '" + opt.pop_b + "': 0 samples matched in VCF");
    }

    log_info(opt, "Analysis: " + opt.pop_a + " (n=" +
                      std::to_string(plan.n_matched_a) + ") vs " + opt.pop_b +
                      " (n=" + std::to_string(plan.n_matched_b) + ")" +
                      (opt.region.empty() ? ", regions=all contigs"
                                          : ", region=" + opt.region));
    return plan;
}

}  // namespace xpclr
