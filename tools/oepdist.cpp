#include <CLI/CLI.hpp>
#include "oecluster/oecluster.h"
#include "MolReader.h"
#include "OutputWriter.h"

#include <iostream>
#include <string>
#include <vector>

#include <oesitehopper.h>

#include <atomic>
#include <csignal>
#include <cstdlib>
#include <mutex>
#include <thread>

#include <indicators/cursor_control.hpp>
#include <indicators/progress_bar.hpp>
#include <indicators/progress_spinner.hpp>

// ---------------------------------------------------------------------------
// Signal handler to restore terminal cursor on SIGINT/SIGTERM
// ---------------------------------------------------------------------------
static void signal_handler(int sig) {
    indicators::show_console_cursor(true);
    // Re-raise with default handler so the shell sees the correct exit status.
    std::signal(sig, SIG_DFL);
    std::raise(sig);
}

// ---------------------------------------------------------------------------
// Progress helpers
// ---------------------------------------------------------------------------

// Uniform prefix width so all bars/spinners align.
static constexpr size_t kPrefixWidth = 21;

static std::string PadPrefix(const std::string& label) {
    if (label.size() >= kPrefixWidth) return label;
    return label + std::string(kPrefixWidth - label.size(), ' ');
}

/// RAII spinner with background animation thread.
/// The callback updates atomic counters; a background thread renders
/// the spinner at a steady ~100ms cadence regardless of read rate.
class ReadSpinner {
public:
    explicit ReadSpinner(const std::string& label, bool enabled)
        : enabled_(enabled), label_(PadPrefix(label)) {
        if (!enabled_) return;
        indicators::show_console_cursor(false);
        spinner_ = std::make_unique<indicators::ProgressSpinner>(
            indicators::option::PrefixText{label_},
            indicators::option::PostfixText{"0 read"},
            indicators::option::ShowPercentage{false},
            indicators::option::ShowElapsedTime{true},
            indicators::option::SpinnerStates{std::vector<std::string>{
                "⠋", "⠙", "⠹", "⠸", "⠼", "⠴", "⠦", "⠧", "⠇", "⠏"}},
            indicators::option::MaxProgress{1000000},
            indicators::option::Stream{std::cerr});
        // Launch background animation thread.
        thread_ = std::thread([this] { animate(); });
    }

    /// Returns a ReadProgress callback that updates atomic counters.
    OEPDist::ReadProgress callback() {
        if (!enabled_) return nullptr;
        return [this](size_t count, double fraction) {
            count_.store(count, std::memory_order_relaxed);
            fraction_.store(fraction, std::memory_order_relaxed);
        };
    }

    void stop(const std::string& summary) {
        if (!enabled_ || stopped_) return;
        stopped_ = true;
        running_.store(false, std::memory_order_release);
        if (thread_.joinable()) thread_.join();
        spinner_->set_option(indicators::option::ShowSpinner{false});
        spinner_->set_option(indicators::option::ShowElapsedTime{false});
        spinner_->set_option(indicators::option::PostfixText{summary});
        spinner_->mark_as_completed();
        std::cerr << std::endl;
    }

    ~ReadSpinner() {
        if (enabled_ && !stopped_) {
            running_.store(false, std::memory_order_release);
            if (thread_.joinable()) thread_.join();
            indicators::show_console_cursor(true);
        }
    }

    ReadSpinner(const ReadSpinner&) = delete;
    ReadSpinner& operator=(const ReadSpinner&) = delete;

private:
    bool enabled_;
    bool stopped_ = false;
    std::string label_;
    std::unique_ptr<indicators::ProgressSpinner> spinner_;
    std::thread thread_;
    std::atomic<bool> running_{true};
    std::atomic<size_t> count_{0};
    std::atomic<double> fraction_{-1.0};

    void animate() {
        while (running_.load(std::memory_order_acquire)) {
            size_t n = count_.load(std::memory_order_relaxed);
            double frac = fraction_.load(std::memory_order_relaxed);
            std::string text = std::to_string(n) + " read";
            if (frac >= 0.0 && frac <= 1.0) {
                int pct = static_cast<int>(frac * 100.0);
                if (pct > 100) pct = 100;
                text += " (" + std::to_string(pct) + "%)";
            }
            spinner_->set_option(indicators::option::PostfixText{text});
            spinner_->tick();
            std::this_thread::sleep_for(std::chrono::milliseconds(100));
        }
    }
};

static indicators::ProgressBar MakeProgressBar(const std::string& label,
                                                size_t max_progress) {
    using namespace indicators;
    return ProgressBar{
        option::BarWidth{40},
        option::PrefixText{PadPrefix(label)},
        option::ShowPercentage{true},
        option::ShowElapsedTime{true},
        option::ShowRemainingTime{true},
        option::MaxProgress{max_progress},
        option::Stream{std::cerr},
    };
}

// ---------------------------------------------------------------------------
// Parsers
// ---------------------------------------------------------------------------

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

static OECluster::SuperposeMethod ParseSuperposeMethod(const std::string& s) {
    if (s == "global_carbon_alpha") return OECluster::SuperposeMethod::GlobalCarbonAlpha;
    if (s == "global")              return OECluster::SuperposeMethod::Global;
    if (s == "ddm")                 return OECluster::SuperposeMethod::DDM;
    if (s == "weighted")            return OECluster::SuperposeMethod::Weighted;
    if (s == "sse")                 return OECluster::SuperposeMethod::SSE;
    if (s == "sitehopper")          return OECluster::SuperposeMethod::SiteHopper;
    throw std::runtime_error("Unknown superpose method: " + s);
}

static OECluster::SuperposeScoreType ParseSuperposeScoreType(const std::string& s) {
    if (s == "auto")        return OECluster::SuperposeScoreType::Auto;
    if (s == "rmsd")        return OECluster::SuperposeScoreType::RMSD;
    if (s == "tanimoto")    return OECluster::SuperposeScoreType::Tanimoto;
    if (s == "patch_score") return OECluster::SuperposeScoreType::PatchScore;
    throw std::runtime_error("Unknown superpose score type: " + s);
}

// ---------------------------------------------------------------------------
// Patch surface generation
// ---------------------------------------------------------------------------

static void AddPatchSurfaces(OEPDist::DesignUnitSet& ds,
                              bool show_progress, bool verbose,
                              size_t num_threads) {
    size_t n = ds.dus.size();
    if (n == 0) return;
    size_t nthreads = num_threads > 0
        ? num_threads
        : std::max(1u, std::thread::hardware_concurrency());
    if (nthreads > n) nthreads = n;

    auto bar = MakeProgressBar("Generating patches", n);
    if (show_progress) bar.set_progress(0);
    std::atomic<size_t> done{0};
    std::mutex warn_mutex;

    auto worker = [&](size_t begin, size_t end) {
        for (size_t i = begin; i < end; ++i) {
            if (!OESiteHopper::OEHasPatchSurface(*ds.dus[i])) {
                if (!OESiteHopper::OEAddPatchSurface(*ds.dus[i])) {
                    std::lock_guard<std::mutex> lock(warn_mutex);
                    if (show_progress) std::cerr << "\n";
                    std::cerr << "Warning: Failed to generate patch surface for "
                              << ds.labels[i] << std::endl;
                } else if (verbose) {
                    std::lock_guard<std::mutex> lock(warn_mutex);
                    if (show_progress) std::cerr << "\n";
                    std::cerr << "Generated patch surface for "
                              << ds.labels[i] << std::endl;
                }
            }
            size_t d = done.fetch_add(1, std::memory_order_relaxed) + 1;
            if (show_progress) bar.set_progress(d);
        }
    };

    if (nthreads == 1) {
        worker(0, n);
    } else {
        std::vector<std::thread> threads;
        threads.reserve(nthreads);
        size_t chunk = n / nthreads;
        size_t remainder = n % nthreads;
        size_t start = 0;
        for (size_t t = 0; t < nthreads; ++t) {
            size_t end = start + chunk + (t < remainder ? 1 : 0);
            threads.emplace_back(worker, start, end);
            start = end;
        }
        for (auto& th : threads) th.join();
    }

    if (show_progress)
        std::cerr << std::endl;
}

// ---------------------------------------------------------------------------
// JSON helpers
// ---------------------------------------------------------------------------

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

// ---------------------------------------------------------------------------
// Distance computation
// ---------------------------------------------------------------------------

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
    auto bar = MakeProgressBar("Computing distances", n_pairs);
    if (progress) {
        bar.set_progress(0);
        opts.progress = [&bar](size_t done, size_t /*total*/) {
            bar.set_progress(done);
        };
    }
    OECluster::pdist(metric, storage, opts);
    if (progress) {
        if (!bar.is_completed())
            bar.mark_as_completed();
        std::cerr << std::endl;
    }

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
    size_t total_pairs = n_a * n_b;
    std::vector<double> data(total_pairs);
    OECluster::CDistOptions opts;
    opts.num_threads = num_threads;
    opts.chunk_size = chunk_size;
    opts.cutoff = cutoff;
    auto bar = MakeProgressBar("Computing distances", total_pairs);
    if (progress) {
        bar.set_progress(0);
        opts.progress = [&bar](size_t done, size_t /*total*/) {
            bar.set_progress(done);
        };
    }
    OECluster::cdist(metric, n_a, data.data(), opts);
    if (progress) {
        if (!bar.is_completed())
            bar.mark_as_completed();
        std::cerr << std::endl;
    }

    OEPDist::WriteCDist(output_path, data.data(), n_a, n_b, meta);
    return 0;
}

// ---------------------------------------------------------------------------
// Common CLI options
// ---------------------------------------------------------------------------

struct CommonOpts {
    std::string input1, input2, output;
    size_t num_threads = 0, chunk_size = 256;
    double cutoff = 0.0;
    bool no_progress = false, verbose = false;

    bool progress() const { return !no_progress; }
};

static void AddCommonOpts(CLI::App* cmd, CommonOpts& c) {
    cmd->add_option("input", c.input1, "Input file")->required();
    cmd->add_option("input2", c.input2, "Second input file (NxM mode)");
    cmd->add_option("-o,--output", c.output, "Output file")->required();
    cmd->add_option("-t,--threads", c.num_threads, "Threads (0=auto)");
    cmd->add_option("-c,--cutoff", c.cutoff, "Distance cutoff");
    cmd->add_option("--chunk-size", c.chunk_size, "Pairs per work unit");
    cmd->add_flag("--no-progress", c.no_progress, "Disable progress output");
    cmd->add_flag("-v,--verbose", c.verbose, "Verbose logging");
}

// ---------------------------------------------------------------------------
// Structure readers with progress
// ---------------------------------------------------------------------------

static OEPDist::MolSet ReadMolsProgress(const std::string& path,
                                          const CommonOpts& co) {
    ReadSpinner spin("Reading structures", co.progress());
    auto ms = OEPDist::ReadMolecules(path, co.verbose, spin.callback());
    spin.stop(std::to_string(ms.owned_mols.size()) + " molecules");
    return ms;
}

static OEPDist::MultiConfMolSet ReadMultiConfProgress(const std::string& path,
                                                        const CommonOpts& co) {
    ReadSpinner spin("Reading structures", co.progress());
    auto ms = OEPDist::ReadMultiConfMolecules(path, co.verbose, spin.callback());
    spin.stop(std::to_string(ms.mols.size()) + " molecules");
    return ms;
}

static OEPDist::DesignUnitSet ReadDUsProgress(const std::string& path,
                                                const CommonOpts& co) {
    ReadSpinner spin("Reading structures", co.progress());
    auto ds = OEPDist::ReadDesignUnits(path, co.verbose, spin.callback());
    spin.stop(std::to_string(ds.dus.size()) + " design units");
    return ds;
}

// ---------------------------------------------------------------------------
// main
// ---------------------------------------------------------------------------

int main(int argc, char** argv) {
    std::signal(SIGINT, signal_handler);
    std::signal(SIGTERM, signal_handler);

    CLI::App app{"oepdist — pairwise and cross-distance matrices"};
    app.require_subcommand(1);

    int exit_code = 0;

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
    bool fp_sim = false;
    fp_cmd->add_flag("--sim", fp_sim, "Return similarity instead of distance");

    fp_cmd->callback([&]() {
        OECluster::FingerprintOptions opts;
        opts.fp_type = fp_type;
        opts.numbits = fp_numbits;
        opts.min_distance = fp_min;
        opts.max_distance = fp_max;
        opts.similarity_func = fp_similarity;
        if (!fp_atom_type.empty())
            opts.atom_type_mask =
                OECluster::FingerprintMetric::ParseAtomTypeMask(fp_atom_type);
        if (!fp_bond_type.empty())
            opts.bond_type_mask =
                OECluster::FingerprintMetric::ParseBondTypeMask(fp_bond_type);
        opts.similarity = fp_sim;

        std::string params = "{" + JsonStr("fp_type", fp_type) + ","
            + JsonNum("numbits", fp_numbits) + ","
            + JsonNum("min_distance", fp_min) + ","
            + JsonNum("max_distance", fp_max) + ","
            + JsonStr("similarity_func", fp_similarity) + "}";

        if (fp_co.input2.empty()) {
            auto ms = ReadMolsProgress(fp_co.input1, fp_co);
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
                fp_co.cutoff, fp_co.progress());
        } else {
            auto sa = ReadMolsProgress(fp_co.input1, fp_co);
            auto sb = ReadMolsProgress(fp_co.input2, fp_co);
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
                fp_co.cutoff, fp_co.progress());
        }
    });

    CommonOpts rocs_co;
    std::string rocs_score = "combo_norm", rocs_color_ff = "implicit-mills-dean";

    auto* rocs_cmd = app.add_subcommand("rocs", "ROCS shape overlay distance");
    AddCommonOpts(rocs_cmd, rocs_co);
    rocs_cmd->add_option("--score", rocs_score,
        "Score: combo_norm, combo, shape, color");
    rocs_cmd->add_option("--color-ff", rocs_color_ff,
        "Color force field");
    bool rocs_sim = false;
    rocs_cmd->add_flag("--sim", rocs_sim, "Return similarity instead of distance");

    rocs_cmd->callback([&]() {
        OECluster::ROCSOptions opts;
        opts.score_type = ParseScoreType(rocs_score);
        opts.color_ff_type = ParseColorFF(rocs_color_ff);
        opts.similarity = rocs_sim;

        std::string params = "{" + JsonStr("score", rocs_score) + ","
            + JsonStr("color_ff", rocs_color_ff) + "}";

        if (rocs_co.input2.empty()) {
            auto ms = ReadMultiConfProgress(rocs_co.input1, rocs_co);
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
                rocs_co.cutoff, rocs_co.progress());
        } else {
            auto sa = ReadMultiConfProgress(rocs_co.input1, rocs_co);
            auto sb = ReadMultiConfProgress(rocs_co.input2, rocs_co);
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
                rocs_co.cutoff, rocs_co.progress());
        }
    });

    CommonOpts sup_co;
    std::string sup_method = "global_carbon_alpha";
    std::string sup_score_type = "auto";
    bool sup_sim = false;
    std::string sup_predicate, sup_ref_predicate, sup_fit_predicate;

    auto* sup_cmd = app.add_subcommand("superpose",
        "Protein superposition distance");
    AddCommonOpts(sup_cmd, sup_co);
    sup_cmd->add_option("--method", sup_method,
        "Method: global_carbon_alpha, global, ddm, weighted, sse, sitehopper");
    sup_cmd->add_option("--score-type", sup_score_type,
        "Score: auto, rmsd, tanimoto, patch_score");
    sup_cmd->add_flag("--sim", sup_sim, "Return similarity instead of distance");
    sup_cmd->add_option("--predicate", sup_predicate,
        "oeselect expression for both ref and fit");
    sup_cmd->add_option("--ref-predicate", sup_ref_predicate,
        "Override predicate for ref structures");
    sup_cmd->add_option("--fit-predicate", sup_fit_predicate,
        "Override predicate for fit structures");

    sup_cmd->callback([&]() {
        OECluster::SuperposeOptions opts;
        opts.method = ParseSuperposeMethod(sup_method);
        opts.score_type = ParseSuperposeScoreType(sup_score_type);
        opts.similarity = sup_sim;
        opts.predicate = sup_predicate;
        opts.ref_predicate = sup_ref_predicate;
        opts.fit_predicate = sup_fit_predicate;

        std::string params = "{" + JsonStr("method", sup_method) + ","
            + JsonStr("score_type", sup_score_type) + ","
            + JsonBool("similarity", sup_sim) + "}";

        if (sup_co.input2.empty()) {
            OEPDist::DesignUnitSet ds;
            try {
                ds = ReadDUsProgress(sup_co.input1, sup_co);
            } catch (...) {}
            if (!ds.dus.empty()) {
                if (opts.method == OECluster::SuperposeMethod::SiteHopper)
                    AddPatchSurfaces(ds, sup_co.progress(), sup_co.verbose,
                                     sup_co.num_threads);
                OECluster::SuperposeMetric metric(ds.dus, opts);
                OEPDist::OutputMetadata meta{
                    "pdist", metric.Name(), params, ds.labels, ds.labels};
                exit_code = run_pdist(metric, sup_co.output, meta,
                    sup_co.num_threads, sup_co.chunk_size,
                    sup_co.cutoff, sup_co.progress());
                return;
            }
            auto ms = ReadMolsProgress(sup_co.input1, sup_co);
            if (ms.ptrs.empty()) {
                std::cerr << "Error: No structures read from "
                          << sup_co.input1 << std::endl;
                exit_code = 1; return;
            }
            OECluster::SuperposeMetric metric(ms.ptrs, opts);
            OEPDist::OutputMetadata meta{
                "pdist", metric.Name(), params, ms.labels, ms.labels};
            exit_code = run_pdist(metric, sup_co.output, meta,
                sup_co.num_threads, sup_co.chunk_size,
                sup_co.cutoff, sup_co.progress());
        } else {
            // Cross-distance: try design units first
            OEPDist::DesignUnitSet sa_du, sb_du;
            try {
                sa_du = ReadDUsProgress(sup_co.input1, sup_co);
                sb_du = ReadDUsProgress(sup_co.input2, sup_co);
            } catch (...) {}
            if (!sa_du.dus.empty() && !sb_du.dus.empty()) {
                if (opts.method == OECluster::SuperposeMethod::SiteHopper) {
                    AddPatchSurfaces(sa_du, sup_co.progress(), sup_co.verbose,
                                     sup_co.num_threads);
                    AddPatchSurfaces(sb_du, sup_co.progress(), sup_co.verbose,
                                     sup_co.num_threads);
                }
                std::vector<std::shared_ptr<OEBio::OEDesignUnit>> all;
                all.insert(all.end(), sa_du.dus.begin(), sa_du.dus.end());
                all.insert(all.end(), sb_du.dus.begin(), sb_du.dus.end());
                OECluster::SuperposeMetric metric(all, opts);
                OEPDist::OutputMetadata meta{
                    "cdist", metric.Name(), params, sa_du.labels, sb_du.labels};
                exit_code = run_cdist(metric, sa_du.dus.size(), sup_co.output,
                    meta, sup_co.num_threads, sup_co.chunk_size,
                    sup_co.cutoff, sup_co.progress());
                return;
            }
            auto sa = ReadMolsProgress(sup_co.input1, sup_co);
            auto sb = ReadMolsProgress(sup_co.input2, sup_co);
            if (sa.ptrs.empty() || sb.ptrs.empty()) {
                std::cerr << "Error: No structures read" << std::endl;
                exit_code = 1; return;
            }
            std::vector<OEChem::OEMolBase*> all;
            all.insert(all.end(), sa.ptrs.begin(), sa.ptrs.end());
            all.insert(all.end(), sb.ptrs.begin(), sb.ptrs.end());
            OECluster::SuperposeMetric metric(all, opts);
            OEPDist::OutputMetadata meta{
                "cdist", metric.Name(), params, sa.labels, sb.labels};
            exit_code = run_cdist(metric, sa.ptrs.size(), sup_co.output,
                meta, sup_co.num_threads, sup_co.chunk_size,
                sup_co.cutoff, sup_co.progress());
        }
    });

    try {
        CLI11_PARSE(app, argc, argv);
    } catch (const std::exception& e) {
        indicators::show_console_cursor(true);
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    indicators::show_console_cursor(true);
    return exit_code;
}
