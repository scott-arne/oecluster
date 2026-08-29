/**
 * @file FingerprintComparison.cpp
 * @brief Implementation of fingerprint-based comparison backed by OEFP.
 */

#include "oecluster/comparisons/FingerprintComparison.h"

#include <algorithm>
#include <cctype>
#include <exception>
#include <oechem.h>
#include <oefp/oefp.h>
#include "oecluster/CDist.h"
#include "oecluster/CondensedIndex.h"
#include "oecluster/Error.h"
#include "oecluster/PDist.h"
#include "oecluster/StorageBackend.h"
#include "MetricTable.h"

namespace OECluster {

// ---------------------------------------------------------------------------
// Impl
// ---------------------------------------------------------------------------

struct FingerprintComparison::Impl {
    std::vector<OEFP::OEFP> fingerprints;
    OEFP::OEFPBatch batch;
    FingerprintOptions opts;
    OEFP::Metric metric = OEFP::Metric::Jaccard();
};

// ---------------------------------------------------------------------------
// Helper: lowercase a string in place
// ---------------------------------------------------------------------------

static std::string to_lower(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char c) { return std::tolower(c); });
    return s;
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

static OEFP::Metric make_metric(const FingerprintOptions& opts) {
    MetricParams params;
    params.p = opts.p;
    params.tversky_alpha = opts.tversky_alpha;
    params.tversky_beta = opts.tversky_beta;
    return resolve_metric(opts.metric, opts.similarity, params, MetricSurface::Fingerprint);
}

static OEFP::BatchKernelOptions make_kernel_options(size_t num_threads,
                                                    size_t chunk_size) {
    OEFP::BatchKernelOptions options;
    options.num_threads = num_threads;
    options.chunk_size = chunk_size > 0 ? chunk_size : 256;
    return options;
}

static OEFP::OEFPBatch make_batch_slice(const std::vector<OEFP::OEFP>& fingerprints,
                                        size_t begin,
                                        size_t end) {
    OEFP::OEFPBatch batch(fingerprints.front().Spec());
    for (size_t i = begin; i < end; ++i) {
        batch.Append(fingerprints[i]);
    }
    return batch;
}

static void validate_molecule(const OEChem::OEMolBase* mol, size_t index) {
    if (mol == nullptr) {
        throw ComparisonError("FingerprintComparison received null molecule pointer at index " +
                          std::to_string(index));
    }
}

/// Fold user-facing family spellings onto the three OEFP generators.
static std::string normalize_family(const std::string& fp_type) {
    const std::string lower = to_lower(fp_type);
    if (lower == "morgan") {
        return "morgan";
    }
    if (lower == "atom_pair" || lower == "atompair" || lower == "topological_atom_pair") {
        // OEFP has a single AtomPairGenerator; the "topological" spelling is the
        // 2D graph-distance model, which is the only one OEFP implements.
        return "atom_pair";
    }
    if (lower == "topological_torsions" || lower == "topological_torsion") {
        return "topological_torsions";
    }
    if (lower == "distance_atom_pair") {
        throw ComparisonError(
            "Fingerprint type 'distance_atom_pair' is not implemented by OEFP 0.3.0; "
            "use 'atom_pair' for the topological (2D) model");
    }
    if (lower == "circular" || lower == "tree" || lower == "path" || lower == "maccs" ||
        lower == "lingo") {
        throw ComparisonError("OpenEye fingerprint type '" + fp_type +
                              "' is no longer supported; use one of 'morgan', 'atom_pair', "
                              "'topological_atom_pair', 'topological_torsions'");
    }
    throw ComparisonError("Unknown OEFP fingerprint type: " + fp_type +
                          ". Supported types are 'morgan', 'atom_pair', "
                          "'topological_atom_pair', 'topological_torsions'");
}

static std::vector<OEFP::OEFP> make_binary_fingerprints(
        const std::vector<OEChem::OEMolBase*>& mols,
        const FingerprintOptions& opts,
        const std::string& family) {
    std::vector<OEFP::OEFP> fingerprints;
    fingerprints.reserve(mols.size());

    if (family == "morgan") {
        OEFP::MorganOptions generator_opts;
        generator_opts.num_bits = opts.numbits;
        generator_opts.radius = opts.radius;
        generator_opts.use_chirality = opts.use_chirality;
        OEFP::MorganGenerator generator(generator_opts);
        for (size_t i = 0; i < mols.size(); ++i) {
            validate_molecule(mols[i], i);
            fingerprints.push_back(generator.Fingerprint(*mols[i]));
        }
    } else if (family == "atom_pair") {
        OEFP::AtomPairOptions generator_opts;
        generator_opts.num_bits = opts.numbits;
        generator_opts.min_distance = opts.min_distance;
        generator_opts.max_distance = opts.max_distance;
        generator_opts.use_chirality = opts.use_chirality;
        generator_opts.use_2d = true;
        OEFP::AtomPairGenerator generator(generator_opts);
        for (size_t i = 0; i < mols.size(); ++i) {
            validate_molecule(mols[i], i);
            fingerprints.push_back(generator.Fingerprint(*mols[i]));
        }
    } else {
        // normalize_family returns exactly three values, so this is the
        // topological torsions case. Any new caller must normalize first.
        OEFP::TopologicalTorsionsOptions generator_opts;
        generator_opts.num_bits = opts.numbits;
        generator_opts.torsion_atom_count = opts.torsion_atom_count;
        generator_opts.use_chirality = opts.use_chirality;
        OEFP::TopologicalTorsionsGenerator generator(generator_opts);
        for (size_t i = 0; i < mols.size(); ++i) {
            validate_molecule(mols[i], i);
            fingerprints.push_back(generator.Fingerprint(*mols[i]));
        }
    }

    return fingerprints;
}

// ---------------------------------------------------------------------------
// Constructor
// ---------------------------------------------------------------------------

FingerprintComparison::FingerprintComparison(const std::vector<OEChem::OEMolBase*>& mols,
                                     const Options& opts) {
    auto impl = std::make_shared<Impl>();
    impl->opts = opts;
    impl->metric = make_metric(opts);

    const std::string family = normalize_family(opts.fp_type);
    try {
        impl->fingerprints = make_binary_fingerprints(mols, opts, family);
    } catch (const ComparisonError&) {
        throw;
    } catch (const std::exception& exc) {
        throw ComparisonError("Failed to compute OEFP fingerprints: " + std::string(exc.what()));
    }
    impl->batch = OEFP::OEFPBatch::FromFingerprints(impl->fingerprints);

    pimpl_ = std::move(impl);
}

// ---------------------------------------------------------------------------
// Private clone constructor
// ---------------------------------------------------------------------------

FingerprintComparison::FingerprintComparison(std::shared_ptr<const Impl> impl)
    : pimpl_(std::move(impl)) {}

// ---------------------------------------------------------------------------
// Distance
// ---------------------------------------------------------------------------

double FingerprintComparison::Compare(size_t i, size_t j) {
    const auto& fp_i = pimpl_->fingerprints[i];
    const auto& fp_j = pimpl_->fingerprints[j];
    return OEFP::Compare(fp_i, fp_j, pimpl_->metric);
}

bool FingerprintComparison::TryPDist(StorageBackend& storage,
                                     const PDistOptions& options) {
    // Asymmetric metrics have no condensed representation. OEFP refuses them
    // inside its kernel; raising here names the OECluster-level cause.
    try {
        pimpl_->metric.ValidateForPDist();
    } catch (const std::exception& exc) {
        throw ComparisonError("Metric '" + pimpl_->opts.metric +
                              "' cannot be used with pdist: " + std::string(exc.what()));
    }

    const size_t n = pimpl_->batch.Size();
    const size_t total_pairs = n * (n - 1) / 2;
    if (storage.NumSamples() != n) {
        throw ComparisonError("FingerprintComparison pdist storage size mismatch");
    }

    const OEFP::BatchKernelOptions kernel_options =
        make_kernel_options(options.num_threads, options.chunk_size);

    try {
        double* data = storage.Data();
        if (data != nullptr) {
            OEFP::PDistInto(
                pimpl_->batch, pimpl_->metric, data, storage.NumPairs(), kernel_options);
        } else {
            const std::vector<double> values =
                OEFP::PDist(pimpl_->batch, pimpl_->metric, kernel_options);
            for (size_t k = 0; k < values.size(); ++k) {
                size_t i = 0;
                size_t j = 0;
                condensed_to_pair(k, n, i, j);
                storage.Set(i, j, values[k]);
            }
        }
    } catch (const std::exception& exc) {
        throw ComparisonError("Failed to compute OEFP fingerprint pdist: " +
                              std::string(exc.what()));
    }

    if (options.progress) {
        options.progress(total_pairs, total_pairs);
    }
    return true;
}

bool FingerprintComparison::TryCDist(size_t n_a, double* output,
                                     const CDistOptions& options) {
    const size_t n_total = pimpl_->batch.Size();
    if (n_a > n_total) {
        throw ComparisonError("FingerprintComparison cdist split index is out of range");
    }
    const size_t n_b = n_total - n_a;
    const size_t total_pairs = n_a * n_b;
    if (total_pairs == 0) {
        return true;
    }

    const OEFP::BatchKernelOptions kernel_options =
        make_kernel_options(options.num_threads, options.chunk_size);

    try {
        const OEFP::OEFPBatch batch_a = make_batch_slice(pimpl_->fingerprints, 0, n_a);
        const OEFP::OEFPBatch batch_b =
            make_batch_slice(pimpl_->fingerprints, n_a, n_total);
        OEFP::CDistInto(
            batch_a, batch_b, pimpl_->metric, output, total_pairs, kernel_options);
    } catch (const std::exception& exc) {
        throw ComparisonError("Failed to compute OEFP fingerprint cdist: " +
                              std::string(exc.what()));
    }

    if (options.cutoff > 0.0) {
        for (size_t i = 0; i < total_pairs; ++i) {
            if (output[i] > options.cutoff) {
                output[i] = 0.0;
            }
        }
    }

    if (options.progress) {
        options.progress(total_pairs, total_pairs);
    }
    return true;
}

// ---------------------------------------------------------------------------
// Clone / Size / Name
// ---------------------------------------------------------------------------

std::unique_ptr<PairwiseComparison> FingerprintComparison::Clone() const {
    return std::unique_ptr<PairwiseComparison>(new FingerprintComparison(pimpl_));
}

size_t FingerprintComparison::Size() const {
    return pimpl_->fingerprints.size();
}

std::string FingerprintComparison::ComparisonName() const {
    return "fingerprint";
}

GateFacts FingerprintComparison::Facts() const {
    GateFacts facts;
    facts.zero_self =
        pimpl_->metric.HasZeroSelfDistance() ? Capability::Yes : Capability::No;
    facts.triangle =
        pimpl_->metric.SatisfiesTriangleInequality() ? Capability::Yes : Capability::No;
    facts.data_integrity = DataIntegrity::Complete;
    return facts;
}

}  // namespace OECluster
