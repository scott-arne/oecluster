/**
 * @file FingerprintMetric.cpp
 * @brief Implementation of fingerprint-based distance metric backed by OEFP.
 */

#include "oecluster/metrics/FingerprintMetric.h"

#include <algorithm>
#include <cctype>
#include <exception>
#include <oechem.h>
#include <oefp/oefp.h>
#include "oecluster/Error.h"

namespace OECluster {

// ---------------------------------------------------------------------------
// Impl
// ---------------------------------------------------------------------------

struct FingerprintMetric::Impl {
    std::vector<OEFP::OEFP> fingerprints;
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
    const std::string sim_lower = to_lower(opts.similarity_func);
    const OEFP::MetricMode mode = opts.similarity
        ? OEFP::MetricMode::Similarity
        : OEFP::MetricMode::Distance;

    if (sim_lower == "tanimoto") {
        return opts.similarity ? OEFP::Metric::Tanimoto() : OEFP::Metric::Jaccard();
    }
    if (sim_lower == "jaccard") {
        if (opts.similarity) {
            throw MetricError("jaccard is a distance metric; use tanimoto for similarity");
        }
        return OEFP::Metric::Jaccard();
    }
    if (sim_lower == "dice") {
        return OEFP::Metric::Dice(mode);
    }
    if (sim_lower == "cosine") {
        return OEFP::Metric::Cosine(mode);
    }
    if (sim_lower == "manhattan") {
        if (opts.similarity) {
            throw MetricError("similarity=true is incompatible with manhattan distance");
        }
        return OEFP::Metric::Manhattan();
    }
    if (sim_lower == "euclidean") {
        throw MetricError("OEFP fingerprint metrics do not support euclidean distance");
    }

    throw MetricError("Unknown similarity function: " + opts.similarity_func);
}

static void validate_molecule(const OEChem::OEMolBase* mol, size_t index) {
    if (mol == nullptr) {
        throw MetricError("FingerprintMetric received null molecule pointer at index " +
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

FingerprintMetric::FingerprintMetric(const std::vector<OEChem::OEMolBase*>& mols,
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
            throw MetricError(
                "OpenEye fingerprint type '" + opts.fp_type +
                "' is no longer supported; use OEFP fingerprint type "
                "'morgan' or 'atom_pair'");
        } else {
            throw MetricError(
                "Unknown OEFP fingerprint type: " + opts.fp_type +
                ". Supported types are 'morgan' and 'atom_pair'");
        }
    } catch (const MetricError&) {
        throw;
    } catch (const std::exception& exc) {
        throw MetricError("Failed to compute OEFP fingerprints: " +
                          std::string(exc.what()));
    }

    pimpl_ = std::move(impl);
}

// ---------------------------------------------------------------------------
// Private clone constructor
// ---------------------------------------------------------------------------

FingerprintMetric::FingerprintMetric(std::shared_ptr<const Impl> impl)
    : pimpl_(std::move(impl)) {}

// ---------------------------------------------------------------------------
// Distance
// ---------------------------------------------------------------------------

double FingerprintMetric::Distance(size_t i, size_t j) {
    const auto& fp_i = pimpl_->fingerprints[i];
    const auto& fp_j = pimpl_->fingerprints[j];
    return OEFP::Compare(fp_i, fp_j, pimpl_->metric);
}

// ---------------------------------------------------------------------------
// Clone / Size / Name
// ---------------------------------------------------------------------------

std::unique_ptr<DistanceMetric> FingerprintMetric::Clone() const {
    return std::unique_ptr<DistanceMetric>(new FingerprintMetric(pimpl_));
}

size_t FingerprintMetric::Size() const {
    return pimpl_->fingerprints.size();
}

std::string FingerprintMetric::Name() const {
    return "fingerprint";
}

}  // namespace OECluster
