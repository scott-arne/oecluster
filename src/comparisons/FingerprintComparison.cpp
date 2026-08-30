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
#include "KernelOptions.h"
#include "MetricTable.h"

namespace OECluster {

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

static std::string to_lower(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char c) { return std::tolower(c); });
    return s;
}

static OEFP::Metric make_metric(const FingerprintOptions& opts) {
    MetricParams params;
    params.p = opts.p;
    params.tversky_alpha = opts.tversky_alpha;
    params.tversky_beta = opts.tversky_beta;
    return resolve_metric(opts.metric, opts.similarity, params, MetricSurface::Fingerprint);
}

/**
 * Type-erased handle over one of OEFP's three batch representations.
 *
 * Immutable after construction, which is what lets ``Clone()`` keep sharing a
 * single ``shared_ptr<const Impl>`` across worker threads.
 */
struct BatchHolder {
    virtual ~BatchHolder() = default;
    virtual size_t Size() const = 0;
    virtual double ComparePair(size_t i, size_t j, const OEFP::Metric& metric) const = 0;
    virtual std::vector<double> PDist(const OEFP::Metric& metric,
                                      const OEFP::BatchKernelOptions& kernel) const = 0;
    virtual void PDistInto(const OEFP::Metric& metric, double* output, size_t length,
                           const OEFP::BatchKernelOptions& kernel) const = 0;
    virtual void CDistInto(size_t n_a, const OEFP::Metric& metric, double* output, size_t length,
                           const OEFP::BatchKernelOptions& kernel) const = 0;
};

template <typename FingerprintT, typename BatchT>
class TypedBatchHolder : public BatchHolder {
public:
    explicit TypedBatchHolder(std::vector<FingerprintT> fingerprints)
        : fingerprints_(std::move(fingerprints)),
          batch_(BatchT::FromFingerprints(fingerprints_)) {}

    size_t Size() const override { return fingerprints_.size(); }

    double ComparePair(size_t i, size_t j, const OEFP::Metric& metric) const override {
        return OEFP::Compare(fingerprints_[i], fingerprints_[j], metric);
    }

    std::vector<double> PDist(const OEFP::Metric& metric,
                              const OEFP::BatchKernelOptions& kernel) const override {
        return OEFP::PDist(batch_, metric, kernel);
    }

    void PDistInto(const OEFP::Metric& metric, double* output, size_t length,
                   const OEFP::BatchKernelOptions& kernel) const override {
        OEFP::PDistInto(batch_, metric, output, length, kernel);
    }

    void CDistInto(size_t n_a, const OEFP::Metric& metric, double* output, size_t length,
                   const OEFP::BatchKernelOptions& kernel) const override {
        const BatchT batch_a = Slice(0, n_a);
        const BatchT batch_b = Slice(n_a, fingerprints_.size());
        OEFP::CDistInto(batch_a, batch_b, metric, output, length, kernel);
    }

private:
    BatchT Slice(size_t begin, size_t end) const {
        BatchT slice(batch_.Spec());
        for (size_t i = begin; i < end; ++i) {
            slice.Append(fingerprints_[i]);
        }
        return slice;
    }

    std::vector<FingerprintT> fingerprints_;
    BatchT batch_;
};

// ---------------------------------------------------------------------------
// Impl
// ---------------------------------------------------------------------------

struct FingerprintComparison::Impl {
    std::shared_ptr<const BatchHolder> holder;
    FingerprintOptions opts;
    OEFP::Metric metric = OEFP::Metric::Jaccard();
};

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

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

static std::string normalize_storage(const std::string& storage) {
    const std::string lower = to_lower(storage);
    if (lower == "binary" || lower == "count" || lower == "sparse" || lower == "sparse_count") {
        return lower;
    }
    throw ComparisonError("Unknown fingerprint storage: " + storage +
                          ". Supported storages are 'binary', 'count', 'sparse', 'sparse_count'");
}

// OEFP's count_simulation approximates a repeat count by setting several bits
// for one feature in a *binary* fingerprint. The count generators already
// carry real counts, so it is switched off for them. On the binary path
// OEFP's per-family defaults are restated rather than unified: Task 5 pinned
// exact distances measured against those defaults, and quietly changing them
// here would move values in tests that are not about storage at all.
static OEFP::MorganOptions morgan_options(const FingerprintOptions& opts) {
    OEFP::MorganOptions generator_opts;
    generator_opts.num_bits = opts.numbits;
    generator_opts.radius = opts.radius;
    generator_opts.use_chirality = opts.use_chirality;
    generator_opts.count_simulation = false;
    return generator_opts;
}

static OEFP::AtomPairOptions atom_pair_options(const FingerprintOptions& opts, bool counted) {
    OEFP::AtomPairOptions generator_opts;
    generator_opts.num_bits = opts.numbits;
    generator_opts.min_distance = opts.min_distance;
    generator_opts.max_distance = opts.max_distance;
    generator_opts.use_chirality = opts.use_chirality;
    generator_opts.use_2d = true;
    generator_opts.count_simulation = !counted;
    return generator_opts;
}

static OEFP::TopologicalTorsionsOptions torsions_options(const FingerprintOptions& opts,
                                                         bool counted) {
    OEFP::TopologicalTorsionsOptions generator_opts;
    generator_opts.num_bits = opts.numbits;
    generator_opts.torsion_atom_count = opts.torsion_atom_count;
    generator_opts.use_chirality = opts.use_chirality;
    generator_opts.count_simulation = !counted;
    return generator_opts;
}

static std::vector<OEFP::OEFP> make_binary_fingerprints(
        const std::vector<OEChem::OEMolBase*>& mols,
        const FingerprintOptions& opts,
        const std::string& family) {
    std::vector<OEFP::OEFP> fingerprints;
    fingerprints.reserve(mols.size());

    if (family == "morgan") {
        OEFP::MorganGenerator generator(morgan_options(opts));
        for (size_t i = 0; i < mols.size(); ++i) {
            validate_molecule(mols[i], i);
            fingerprints.push_back(generator.Fingerprint(*mols[i]));
        }
    } else if (family == "atom_pair") {
        OEFP::AtomPairGenerator generator(atom_pair_options(opts, /*counted=*/false));
        for (size_t i = 0; i < mols.size(); ++i) {
            validate_molecule(mols[i], i);
            fingerprints.push_back(generator.Fingerprint(*mols[i]));
        }
    } else {
        // normalize_family returns exactly three values, so this is the
        // topological torsions case. Any new caller must normalize first.
        OEFP::TopologicalTorsionsGenerator generator(
            torsions_options(opts, /*counted=*/false));
        for (size_t i = 0; i < mols.size(); ++i) {
            validate_molecule(mols[i], i);
            fingerprints.push_back(generator.Fingerprint(*mols[i]));
        }
    }

    return fingerprints;
}

/// Build the per-molecule fingerprints for one family/storage cell and wrap
/// them in the matching batch holder.
static std::shared_ptr<const BatchHolder> make_holder(
        const std::vector<OEChem::OEMolBase*>& mols,
        const FingerprintOptions& opts,
        const std::string& family,
        const std::string& storage) {
    for (size_t i = 0; i < mols.size(); ++i) {
        validate_molecule(mols[i], i);
    }

    if (storage == "binary") {
        return std::make_shared<TypedBatchHolder<OEFP::OEFP, OEFP::OEFPBatch>>(
            make_binary_fingerprints(mols, opts, family));
    }

    if (storage == "sparse") {
        std::vector<OEFP::OEFPSparse> fingerprints;
        fingerprints.reserve(mols.size());
        for (OEChem::OEMolBase* mol : mols) {
            // Sparse storage is still a set of on-bits, so it takes the binary
            // path's count-simulation setting, not the count path's.
            if (family == "morgan") {
                fingerprints.push_back(OEFP::MakeMorganSparseFingerprint(*mol, morgan_options(opts)));
            } else if (family == "atom_pair") {
                fingerprints.push_back(OEFP::MakeAtomPairSparseFingerprint(
                    *mol, atom_pair_options(opts, /*counted=*/false)));
            } else {
                fingerprints.push_back(OEFP::MakeTopologicalTorsionsSparseFingerprint(
                    *mol, torsions_options(opts, /*counted=*/false)));
            }
        }
        return std::make_shared<TypedBatchHolder<OEFP::OEFPSparse, OEFP::OEFPSparseBatch>>(
            std::move(fingerprints));
    }

    std::vector<OEFP::OEFPCount> fingerprints;
    fingerprints.reserve(mols.size());
    const bool folded = (storage == "count");
    for (OEChem::OEMolBase* mol : mols) {
        if (family == "morgan") {
            fingerprints.push_back(
                folded ? OEFP::MakeMorganCountFingerprint(*mol, morgan_options(opts))
                       : OEFP::MakeMorganSparseCountFingerprint(*mol, morgan_options(opts)));
        } else if (family == "atom_pair") {
            const OEFP::AtomPairOptions generator_opts = atom_pair_options(opts, /*counted=*/true);
            fingerprints.push_back(
                folded ? OEFP::MakeAtomPairCountFingerprint(*mol, generator_opts)
                       : OEFP::MakeAtomPairSparseCountFingerprint(*mol, generator_opts));
        } else {
            fingerprints.push_back(OEFP::MakeTopologicalTorsionsCountFingerprint(
                *mol, torsions_options(opts, /*counted=*/true)));
        }
    }
    return std::make_shared<TypedBatchHolder<OEFP::OEFPCount, OEFP::OEFPCountBatch>>(
        std::move(fingerprints));
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
    const std::string storage = normalize_storage(opts.storage);

    // An unsupported cell cannot be fixed by changing the metric, so this guard
    // must run before the metric rule.
    if (storage == "sparse_count" && family == "topological_torsions") {
        throw ComparisonError(
            "Fingerprint type 'topological_torsions' does not support storage='sparse_count': "
            "OEFP returns OEFPCount64 for that combination and provides no batch or bulk kernel "
            "for it. Use 'binary', 'count', or 'sparse'");
    }

    // A boolean metric binarizes its inputs, discarding exactly the counts that
    // are the reason to select a counted storage.
    if ((storage == "count" || storage == "sparse_count") &&
        impl->metric.Space() == OEFP::MetricSpace::Boolean) {
        throw ComparisonError("Metric '" + opts.metric + "' is a bit-set metric and discards the "
                              "counts in storage='" + storage +
                              "'. Use a numeric metric such as 'manhattan', 'canberra', or "
                              "'bray_curtis'");
    }

    try {
        impl->holder = make_holder(mols, opts, family, storage);
    } catch (const ComparisonError&) {
        throw;
    } catch (const std::exception& exc) {
        throw ComparisonError("Failed to compute OEFP fingerprints: " + std::string(exc.what()));
    }

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
    return pimpl_->holder->ComparePair(i, j, pimpl_->metric);
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

    const size_t n = pimpl_->holder->Size();
    const size_t total_pairs = n * (n - 1) / 2;
    if (storage.NumSamples() != n) {
        throw ComparisonError("FingerprintComparison pdist storage size mismatch");
    }

    const OEFP::BatchKernelOptions kernel_options =
        make_kernel_options(options.num_threads, options.chunk_size);

    try {
        double* data = storage.Data();
        if (data != nullptr) {
            pimpl_->holder->PDistInto(pimpl_->metric, data, storage.NumPairs(), kernel_options);
        } else {
            const std::vector<double> values =
                pimpl_->holder->PDist(pimpl_->metric, kernel_options);
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
    const size_t n_total = pimpl_->holder->Size();
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
        pimpl_->holder->CDistInto(n_a, pimpl_->metric, output, total_pairs, kernel_options);
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
    return pimpl_->holder->Size();
}

std::string FingerprintComparison::ComparisonName() const {
    return "fingerprint";
}

GateFacts FingerprintComparison::Facts() const {
    GateFacts facts;
    facts.is_distance = pimpl_->metric.Type() == OEFP::MetricType::Distance
                            ? Capability::Yes
                            : Capability::No;
    facts.zero_self =
        pimpl_->metric.HasZeroSelfDistance() ? Capability::Yes : Capability::No;
    facts.triangle =
        pimpl_->metric.SatisfiesTriangleInequality() ? Capability::Yes : Capability::No;
    facts.data_integrity = DataIntegrity::Complete;
    return facts;
}

}  // namespace OECluster
