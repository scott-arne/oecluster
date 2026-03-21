#include <CLI/CLI.hpp>
#include "oecluster/oecluster.h"
#include "MolReader.h"
#include "OutputWriter.h"

#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#ifdef OECLUSTER_HAS_SHAPE
static OECluster::ROCSScoreType ParseScoreType(const std::string& s) {
    if (s == "combo_norm") return OECluster::ROCSScoreType::ComboNorm;
    if (s == "combo")      return OECluster::ROCSScoreType::Combo;
    if (s == "shape")      return OECluster::ROCSScoreType::Shape;
    if (s == "color")      return OECluster::ROCSScoreType::Color;
    throw std::runtime_error("Unknown score type: " + s);
}

static unsigned int ParseColorFF(const std::string& s) {
    if (s == "implicit-mills-dean") return 1;
    if (s == "explicit-mills-dean") return 2;
    if (s == "implicit-no-rings")   return 3;
    if (s == "explicit-no-rings")   return 4;
    throw std::runtime_error("Unknown color force field: " + s);
}
#endif

#ifdef OECLUSTER_HAS_BIO
static unsigned int ParseAlignmentMethod(const std::string& s) {
    if (s == "identity") return 1;
    if (s == "pam250")   return 2;
    if (s == "blosum62") return 3;
    if (s == "gonnet")   return 4;
    throw std::runtime_error("Unknown alignment method: " + s);
}
#endif

static std::string JsonStr(const std::string& key, const std::string& val) {
    return "\"" + key + "\":\"" + val + "\"";
}

template<typename T>
static std::string JsonNum(const std::string& key, T val) {
    return "\"" + key + "\":" + std::to_string(val);
}

static std::string JsonBool(const std::string& key, bool val) {
    return "\"" + key + "\":" + (val ? "true" : "false");
}

template<typename MetricT>
static int run_pdist(MetricT& metric, const std::string& output_path,
                     const OEPDist::OutputMetadata& meta,
                     size_t num_threads, size_t chunk_size,
                     double cutoff, bool progress) {
    size_t n = metric.Size();
    if (n == 0) {
        std::cerr << "Error: No items loaded" << std::endl;
        return 1;
    }
    size_t n_pairs = n * (n - 1) / 2;
    OECluster::DenseStorage storage(n);
    OECluster::PDistOptions opts;
    opts.num_threads = num_threads;
    opts.chunk_size = chunk_size;
    opts.cutoff = cutoff;
    if (progress) {
        opts.progress = [](size_t done, size_t total) {
            std::cerr << "\r" << done << "/" << total << " pairs" << std::flush;
        };
    }
    OECluster::pdist(metric, storage, opts);
    if (progress) std::cerr << std::endl;

    std::vector<double> data(n_pairs);
    size_t k = 0;
    for (size_t i = 0; i < n; ++i)
        for (size_t j = i + 1; j < n; ++j)
            data[k++] = storage.Get(i, j);

    OEPDist::WritePDist(output_path, data.data(), n, meta);
    return 0;
}

template<typename MetricT>
static int run_cdist(MetricT& metric, size_t n_a,
                     const std::string& output_path,
                     const OEPDist::OutputMetadata& meta,
                     size_t num_threads, size_t chunk_size,
                     double cutoff, bool progress) {
    size_t n_b = metric.Size() - n_a;
    if (n_a == 0 || n_b == 0) {
        std::cerr << "Error: Empty input set" << std::endl;
        return 1;
    }
    std::vector<double> data(n_a * n_b);
    OECluster::CDistOptions opts;
    opts.num_threads = num_threads;
    opts.chunk_size = chunk_size;
    opts.cutoff = cutoff;
    if (progress) {
        opts.progress = [](size_t done, size_t total) {
            std::cerr << "\r" << done << "/" << total << " pairs" << std::flush;
        };
    }
    OECluster::cdist(metric, n_a, data.data(), opts);
    if (progress) std::cerr << std::endl;

    OEPDist::WriteCDist(output_path, data.data(), n_a, n_b, meta);
    return 0;
}

struct CommonOpts {
    std::string input1, input2, output;
    size_t num_threads = 0, chunk_size = 256;
    double cutoff = 0.0;
    bool progress = false, verbose = false;
};

static void AddCommonOpts(CLI::App* cmd, CommonOpts& c) {
    cmd->add_option("input", c.input1, "Input file")->required();
    cmd->add_option("input2", c.input2, "Second input file (NxM mode)");
    cmd->add_option("-o,--output", c.output, "Output file")->required();
    cmd->add_option("-t,--threads", c.num_threads, "Threads (0=auto)");
    cmd->add_option("-c,--cutoff", c.cutoff, "Distance cutoff");
    cmd->add_option("--chunk-size", c.chunk_size, "Pairs per work unit");
    cmd->add_flag("--progress", c.progress, "Show progress");
    cmd->add_flag("-v,--verbose", c.verbose, "Verbose logging");
}

int main(int argc, char** argv) {
    CLI::App app{"oepdist — pairwise and cross-distance matrices"};
    app.require_subcommand(1);

    int exit_code = 0;

#ifdef OECLUSTER_HAS_GRAPHSIM
    CommonOpts fp_co;
    std::string fp_type = "circular", fp_atom_type, fp_bond_type;
    std::string fp_similarity = "tanimoto";
    unsigned int fp_numbits = 2048, fp_min = 0, fp_max = 2;

    auto* fp_cmd = app.add_subcommand("fp", "Fingerprint distance");
    AddCommonOpts(fp_cmd, fp_co);
    fp_cmd->add_option("--fp-type", fp_type,
        "Fingerprint type: circular, tree, path, maccs, lingo");
    fp_cmd->add_option("--numbits", fp_numbits, "Number of bits");
    fp_cmd->add_option("--min-distance", fp_min, "Min radius/bond path");
    fp_cmd->add_option("--max-distance", fp_max, "Max radius/bond path");
    fp_cmd->add_option("--atom-type", fp_atom_type,
        "Pipe-delimited OEFPAtomType flags");
    fp_cmd->add_option("--bond-type", fp_bond_type,
        "Pipe-delimited OEFPBondType flags");
    fp_cmd->add_option("--similarity", fp_similarity,
        "Similarity: tanimoto, dice, cosine, manhattan, euclidean");

    fp_cmd->callback([&]() {
        OECluster::FingerprintOptions opts;
        opts.fp_type = fp_type;
        opts.numbits = fp_numbits;
        opts.min_distance = fp_min;
        opts.max_distance = fp_max;
        opts.similarity = fp_similarity;
        if (!fp_atom_type.empty())
            opts.atom_type_mask =
                OECluster::FingerprintMetric::ParseAtomTypeMask(fp_atom_type);
        if (!fp_bond_type.empty())
            opts.bond_type_mask =
                OECluster::FingerprintMetric::ParseBondTypeMask(fp_bond_type);

        std::string params = "{" + JsonStr("fp_type", fp_type) + ","
            + JsonNum("numbits", fp_numbits) + ","
            + JsonNum("min_distance", fp_min) + ","
            + JsonNum("max_distance", fp_max) + ","
            + JsonStr("similarity", fp_similarity) + "}";

        if (fp_co.input2.empty()) {
            auto ms = OEPDist::ReadMolecules(fp_co.input1, fp_co.verbose);
            if (ms.ptrs.empty()) {
                std::cerr << "Error: No molecules read from "
                          << fp_co.input1 << std::endl;
                exit_code = 1; return;
            }
            OECluster::FingerprintMetric metric(ms.ptrs, opts);
            OEPDist::OutputMetadata meta{
                "pdist", "fingerprint", params, ms.labels, ms.labels};
            exit_code = run_pdist(metric, fp_co.output, meta,
                fp_co.num_threads, fp_co.chunk_size,
                fp_co.cutoff, fp_co.progress);
        } else {
            auto sa = OEPDist::ReadMolecules(fp_co.input1, fp_co.verbose);
            auto sb = OEPDist::ReadMolecules(fp_co.input2, fp_co.verbose);
            if (sa.ptrs.empty()) {
                std::cerr << "Error: No molecules read from "
                          << fp_co.input1 << std::endl;
                exit_code = 1; return;
            }
            if (sb.ptrs.empty()) {
                std::cerr << "Error: No molecules read from "
                          << fp_co.input2 << std::endl;
                exit_code = 1; return;
            }
            std::vector<OEChem::OEMolBase*> all;
            all.insert(all.end(), sa.ptrs.begin(), sa.ptrs.end());
            all.insert(all.end(), sb.ptrs.begin(), sb.ptrs.end());
            OECluster::FingerprintMetric metric(all, opts);
            OEPDist::OutputMetadata meta{
                "cdist", "fingerprint", params, sa.labels, sb.labels};
            exit_code = run_cdist(metric, sa.ptrs.size(), fp_co.output, meta,
                fp_co.num_threads, fp_co.chunk_size,
                fp_co.cutoff, fp_co.progress);
        }
    });
#endif

#ifdef OECLUSTER_HAS_SHAPE
    CommonOpts rocs_co;
    std::string rocs_score = "combo_norm", rocs_color_ff = "implicit-mills-dean";

    auto* rocs_cmd = app.add_subcommand("rocs", "ROCS shape overlay distance");
    AddCommonOpts(rocs_cmd, rocs_co);
    rocs_cmd->add_option("--score", rocs_score,
        "Score: combo_norm, combo, shape, color");
    rocs_cmd->add_option("--color-ff", rocs_color_ff,
        "Color force field");

    rocs_cmd->callback([&]() {
        OECluster::ROCSOptions opts;
        opts.score_type = ParseScoreType(rocs_score);
        opts.color_ff_type = ParseColorFF(rocs_color_ff);

        std::string params = "{" + JsonStr("score", rocs_score) + ","
            + JsonStr("color_ff", rocs_color_ff) + "}";

        if (rocs_co.input2.empty()) {
            auto ms = OEPDist::ReadMultiConfMolecules(
                rocs_co.input1, rocs_co.verbose);
            if (ms.mols.empty()) {
                std::cerr << "Error: No molecules read from "
                          << rocs_co.input1 << std::endl;
                exit_code = 1; return;
            }
            OECluster::ROCSMetric metric(ms.mols, opts);
            OEPDist::OutputMetadata meta{
                "pdist", "rocs", params, ms.labels, ms.labels};
            exit_code = run_pdist(metric, rocs_co.output, meta,
                rocs_co.num_threads, rocs_co.chunk_size,
                rocs_co.cutoff, rocs_co.progress);
        } else {
            auto sa = OEPDist::ReadMultiConfMolecules(
                rocs_co.input1, rocs_co.verbose);
            auto sb = OEPDist::ReadMultiConfMolecules(
                rocs_co.input2, rocs_co.verbose);
            if (sa.mols.empty() || sb.mols.empty()) {
                std::cerr << "Error: No molecules read" << std::endl;
                exit_code = 1; return;
            }
            std::vector<std::shared_ptr<OEChem::OEMol>> all;
            all.insert(all.end(), sa.mols.begin(), sa.mols.end());
            all.insert(all.end(), sb.mols.begin(), sb.mols.end());
            OECluster::ROCSMetric metric(all, opts);
            OEPDist::OutputMetadata meta{
                "cdist", "rocs", params, sa.labels, sb.labels};
            exit_code = run_cdist(metric, sa.mols.size(), rocs_co.output,
                meta, rocs_co.num_threads, rocs_co.chunk_size,
                rocs_co.cutoff, rocs_co.progress);
        }
    });
#endif

#ifdef OECLUSTER_HAS_BIO
    CommonOpts sup_co;
    std::string sup_alignment = "pam250";
    int sup_gap = -10, sup_extend = -2;
    bool sup_all_atoms = false, sup_no_overlay = false;

    auto* sup_cmd = app.add_subcommand("superpose",
        "Protein superposition RMSD");
    AddCommonOpts(sup_cmd, sup_co);
    sup_cmd->add_option("--alignment", sup_alignment,
        "Alignment: pam250, blosum62, gonnet, identity");
    sup_cmd->add_option("--gap-penalty", sup_gap, "Gap open penalty");
    sup_cmd->add_option("--extend-penalty", sup_extend, "Gap extend penalty");
    sup_cmd->add_flag("--all-atoms", sup_all_atoms,
        "All backbone atoms (not just C-alpha)");
    sup_cmd->add_flag("--no-overlay", sup_no_overlay,
        "RMSD without spatial overlay");

    sup_cmd->callback([&]() {
        OECluster::SuperposeOptions opts;
        opts.alignment_method = ParseAlignmentMethod(sup_alignment);
        opts.gap_penalty = sup_gap;
        opts.extend_penalty = sup_extend;
        opts.only_calpha = !sup_all_atoms;
        opts.overlay = !sup_no_overlay;

        std::string params = "{" + JsonStr("alignment", sup_alignment) + ","
            + JsonNum("gap_penalty", sup_gap) + ","
            + JsonNum("extend_penalty", sup_extend) + ","
            + JsonBool("only_calpha", opts.only_calpha) + ","
            + JsonBool("overlay", opts.overlay) + "}";

        if (sup_co.input2.empty()) {
            try {
                auto ds = OEPDist::ReadDesignUnits(
                    sup_co.input1, sup_co.verbose);
                if (!ds.dus.empty()) {
                    OECluster::SuperposeMetric metric(ds.dus, opts);
                    OEPDist::OutputMetadata meta{
                        "pdist", "superpose", params, ds.labels, ds.labels};
                    exit_code = run_pdist(metric, sup_co.output, meta,
                        sup_co.num_threads, sup_co.chunk_size,
                        sup_co.cutoff, sup_co.progress);
                    return;
                }
            } catch (...) {}
            auto ms = OEPDist::ReadMolecules(sup_co.input1, sup_co.verbose);
            if (ms.ptrs.empty()) {
                std::cerr << "Error: No structures read from "
                          << sup_co.input1 << std::endl;
                exit_code = 1; return;
            }
            OECluster::SuperposeMetric metric(ms.ptrs, opts);
            OEPDist::OutputMetadata meta{
                "pdist", "superpose", params, ms.labels, ms.labels};
            exit_code = run_pdist(metric, sup_co.output, meta,
                sup_co.num_threads, sup_co.chunk_size,
                sup_co.cutoff, sup_co.progress);
        } else {
            auto sa = OEPDist::ReadMolecules(sup_co.input1, sup_co.verbose);
            auto sb = OEPDist::ReadMolecules(sup_co.input2, sup_co.verbose);
            if (sa.ptrs.empty() || sb.ptrs.empty()) {
                std::cerr << "Error: No structures read" << std::endl;
                exit_code = 1; return;
            }
            std::vector<OEChem::OEMolBase*> all;
            all.insert(all.end(), sa.ptrs.begin(), sa.ptrs.end());
            all.insert(all.end(), sb.ptrs.begin(), sb.ptrs.end());
            OECluster::SuperposeMetric metric(all, opts);
            OEPDist::OutputMetadata meta{
                "cdist", "superpose", params, sa.labels, sb.labels};
            exit_code = run_cdist(metric, sa.ptrs.size(), sup_co.output,
                meta, sup_co.num_threads, sup_co.chunk_size,
                sup_co.cutoff, sup_co.progress);
        }
    });

    CommonOpts sh_co;
    std::string sh_alignment = "pam250";
    int sh_gap = -10, sh_extend = -2;
    bool sh_all_atoms = false;

    auto* sh_cmd = app.add_subcommand("sitehopper",
        "Binding site RMSD");
    AddCommonOpts(sh_cmd, sh_co);
    sh_cmd->add_option("--alignment", sh_alignment,
        "Alignment: pam250, blosum62, gonnet, identity");
    sh_cmd->add_option("--gap-penalty", sh_gap, "Gap open penalty");
    sh_cmd->add_option("--extend-penalty", sh_extend, "Gap extend penalty");
    sh_cmd->add_flag("--all-atoms", sh_all_atoms,
        "All backbone atoms (not just C-alpha)");

    sh_cmd->callback([&]() {
        OECluster::SiteHopperOptions opts;
        opts.alignment_method = ParseAlignmentMethod(sh_alignment);
        opts.gap_penalty = sh_gap;
        opts.extend_penalty = sh_extend;
        opts.only_calpha = !sh_all_atoms;

        std::string params = "{" + JsonStr("alignment", sh_alignment) + ","
            + JsonNum("gap_penalty", sh_gap) + ","
            + JsonNum("extend_penalty", sh_extend) + ","
            + JsonBool("only_calpha", opts.only_calpha) + "}";

        if (sh_co.input2.empty()) {
            auto ds = OEPDist::ReadDesignUnits(
                sh_co.input1, sh_co.verbose);
            if (ds.dus.empty()) {
                std::cerr << "Error: No design units read from "
                          << sh_co.input1 << std::endl;
                exit_code = 1; return;
            }
            OECluster::SiteHopperMetric metric(ds.dus, opts);
            OEPDist::OutputMetadata meta{
                "pdist", "sitehopper", params, ds.labels, ds.labels};
            exit_code = run_pdist(metric, sh_co.output, meta,
                sh_co.num_threads, sh_co.chunk_size,
                sh_co.cutoff, sh_co.progress);
        } else {
            auto sa = OEPDist::ReadDesignUnits(
                sh_co.input1, sh_co.verbose);
            auto sb = OEPDist::ReadDesignUnits(
                sh_co.input2, sh_co.verbose);
            if (sa.dus.empty() || sb.dus.empty()) {
                std::cerr << "Error: No design units read" << std::endl;
                exit_code = 1; return;
            }
            std::vector<std::shared_ptr<OEBio::OEDesignUnit>> all;
            all.insert(all.end(), sa.dus.begin(), sa.dus.end());
            all.insert(all.end(), sb.dus.begin(), sb.dus.end());
            OECluster::SiteHopperMetric metric(all, opts);
            OEPDist::OutputMetadata meta{
                "cdist", "sitehopper", params, sa.labels, sb.labels};
            exit_code = run_cdist(metric, sa.dus.size(), sh_co.output,
                meta, sh_co.num_threads, sh_co.chunk_size,
                sh_co.cutoff, sh_co.progress);
        }
    });
#endif

    try {
        CLI11_PARSE(app, argc, argv);
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    return exit_code;
}
