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
#include "oecluster/Error.h"
#include "oecluster/PDist.h"
#include "oecluster/StorageBackend.h"

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
    const std::string metric_lower = to_lower(opts.metric);

    if (metric_lower == "tanimoto") {
        return opts.similarity ? OEFP::Metric::Tanimoto() : OEFP::Metric::Jaccard();
    }
    if (metric_lower == "jaccard") {
        if (opts.similarity) {
            throw ComparisonError("jaccard is a distance metric; use tanimoto for similarity");
        }
        return OEFP::Metric::Jaccard();
    }
    if (metric_lower == "dice") {
        if (opts.similarity) {
            throw ComparisonError("similarity=true is incompatible with dice distance");
        }
        return OEFP::Metric::Dice();
    }
    if (metric_lower == "cosine") {
        throw ComparisonError("OEFP fingerprint metrics do not support cosine");
    }
    if (metric_lower == "manhattan") {
        if (opts.similarity) {
            throw ComparisonError("similarity=true is incompatible with manhattan distance");
        }
        return OEFP::Metric::Manhattan();
    }
    if (metric_lower == "euclidean") {
        throw ComparisonError("OEFP fingerprint metrics do not support euclidean distance");
    }

    throw ComparisonError("Unknown OEFP fingerprint metric: " + opts.metric);
}

static OEFP::BatchKernelOptions make_kernel_options(size_t num_threads,
                                                    size_t chunk_size) {
    OEFP::BatchKernelOptions options;
    options.num_threads = num_threads;
    options.chunk_size = chunk_size > 0 ? chunk_size : 256;
    return options;
}

static void condensed_to_pair(size_t k, size_t n, size_t& out_i, size_t& out_j) {
    size_t row_start = 0;
    for (size_t i = 0; i < n; ++i) {
        const size_t row_pairs = n - i - 1;
        if (k < row_start + row_pairs) {
            out_i = i;
            out_j = i + 1 + (k - row_start);
            return;
        }
        row_start += row_pairs;
    }
    throw ComparisonError("Condensed fingerprint index is out of range");
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

static std::vector<OEFP::OEFP> make_morgan_fingerprints(
        const std::vector<OEChem::OEMolBase*>& mols,
        const FingerprintOptions& opts) {
    OEFP::MorganOptions morgan_opts;
    morgan_opts.num_bits = opts.numbits;
    morgan_opts.radius = opts.max_distance;

    OEFP::MorganGenerator generator(morgan_opts);
    std::vector<OEFP::OEFP> fingerprints;
    fingerprints.reserve(mols.size());
    for (size_t i = 0; i < mols.size(); ++i) {
        validate_molecule(mols[i], i);
        fingerprints.push_back(generator.Fingerprint(*mols[i]));
    }
    return fingerprints;
}

static std::vector<OEFP::OEFP> make_atom_pair_fingerprints(
        const std::vector<OEChem::OEMolBase*>& mols,
        const FingerprintOptions& opts) {
    OEFP::AtomPairOptions atom_pair_opts;
    atom_pair_opts.num_bits = opts.numbits;
    atom_pair_opts.min_distance = opts.min_distance;
    atom_pair_opts.max_distance = opts.max_distance;

    OEFP::AtomPairGenerator generator(atom_pair_opts);
    std::vector<OEFP::OEFP> fingerprints;
    fingerprints.reserve(mols.size());
    for (size_t i = 0; i < mols.size(); ++i) {
        validate_molecule(mols[i], i);
        fingerprints.push_back(generator.Fingerprint(*mols[i]));
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

    const std::string fp_lower = to_lower(opts.fp_type);
    try {
        if (fp_lower == "morgan") {
            impl->fingerprints = make_morgan_fingerprints(mols, opts);
        } else if (fp_lower == "atom_pair" || fp_lower == "atompair") {
            impl->fingerprints = make_atom_pair_fingerprints(mols, opts);
        } else if (fp_lower == "circular" || fp_lower == "tree" ||
                   fp_lower == "path" || fp_lower == "maccs" ||
                   fp_lower == "lingo") {
            throw ComparisonError(
                "OpenEye fingerprint type '" + opts.fp_type +
                "' is no longer supported; use OEFP fingerprint type "
                "'morgan' or 'atom_pair'");
        } else {
            throw ComparisonError(
                "Unknown OEFP fingerprint type: " + opts.fp_type +
                ". Supported types are 'morgan' and 'atom_pair'");
        }
    } catch (const ComparisonError&) {
        throw;
    } catch (const std::exception& exc) {
        throw ComparisonError("Failed to compute OEFP fingerprints: " +
                          std::string(exc.what()));
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
    const size_t n = pimpl_->batch.Size();
    const size_t total_pairs = n * (n - 1) / 2;
    if (storage.NumItems() != n) {
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

}  // namespace OECluster
