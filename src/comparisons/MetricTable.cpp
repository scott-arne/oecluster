/**
 * @file MetricTable.cpp
 * @brief Implementation of user-facing metric name resolution.
 */

#include "MetricTable.h"

#include <algorithm>
#include <cctype>
#include <sstream>
#include "oecluster/Error.h"

namespace OECluster {

namespace {

std::string to_lower(const std::string& value) {
    std::string result = value;
    std::transform(result.begin(), result.end(), result.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return result;
}

/// One row of the metric table.
struct MetricEntry {
    const char* name;
    bool on_fingerprint;
    bool on_descriptor;
    bool has_similarity_form;
    OEFP::Metric (*make)(bool similarity, const MetricParams& params);
};

// Capture-less lambdas convert to plain function pointers, which keeps the
// table a constant array with no per-call allocation.
const MetricEntry METRIC_TABLE[] = {
    {"jaccard", true, false, false,
     [](bool, const MetricParams&) { return OEFP::Metric::Jaccard(); }},
    {"tanimoto", true, false, true,
     [](bool similarity, const MetricParams&) {
         return similarity ? OEFP::Metric::Tanimoto() : OEFP::Metric::Jaccard();
     }},
    {"dice", true, false, false, [](bool, const MetricParams&) { return OEFP::Metric::Dice(); }},
    {"sokal_sneath", true, false, false,
     [](bool, const MetricParams&) { return OEFP::Metric::SokalSneath(); }},
    {"matching", true, false, false,
     [](bool, const MetricParams&) { return OEFP::Metric::Matching(); }},
    {"rogers_tanimoto", true, false, false,
     [](bool, const MetricParams&) { return OEFP::Metric::RogersTanimoto(); }},
    {"russell_rao", true, false, false,
     [](bool, const MetricParams&) { return OEFP::Metric::RussellRao(); }},
    {"kulsinski", true, false, false,
     [](bool, const MetricParams&) { return OEFP::Metric::Kulsinski(); }},
    {"sokal_michener", true, false, false,
     [](bool, const MetricParams&) { return OEFP::Metric::SokalMichener(); }},
    {"euclidean", true, true, false,
     [](bool, const MetricParams&) { return OEFP::Metric::Euclidean(); }},
    {"manhattan", true, true, false,
     [](bool, const MetricParams&) { return OEFP::Metric::Manhattan(); }},
    {"chebyshev", true, true, false,
     [](bool, const MetricParams&) { return OEFP::Metric::Chebyshev(); }},
    {"hamming", true, true, false,
     [](bool, const MetricParams&) { return OEFP::Metric::Hamming(); }},
    {"canberra", true, true, false,
     [](bool, const MetricParams&) { return OEFP::Metric::Canberra(); }},
    {"bray_curtis", true, true, false,
     [](bool, const MetricParams&) { return OEFP::Metric::BrayCurtis(); }},
    {"minkowski", true, true, false,
     [](bool, const MetricParams& params) { return OEFP::Metric::Minkowski(params.p); }},
    {"tversky", true, false, true,
     [](bool, const MetricParams& params) {
         return OEFP::Metric::Tversky(params.tversky_alpha, params.tversky_beta);
     }},
    {"standardized_euclidean", false, true, false,
     [](bool, const MetricParams& params) {
         return OEFP::Metric::StandardizedEuclidean(params.variances);
     }},
    {"seuclidean", false, true, false,
     [](bool, const MetricParams& params) {
         return OEFP::Metric::StandardizedEuclidean(params.variances);
     }},
    {"mahalanobis", false, true, false,
     [](bool, const MetricParams& params) {
         return OEFP::Metric::Mahalanobis(params.inverse_covariance);
     }},
};

bool visible_on(const MetricEntry& entry, MetricSurface surface) {
    return surface == MetricSurface::Fingerprint ? entry.on_fingerprint : entry.on_descriptor;
}

/// Reject parameter values OEFP would either accept silently or fail on obscurely.
void validate_params(const std::string& name, const MetricParams& params) {
    if (name == "minkowski" && !(params.p > 0.0)) {
        std::ostringstream message;
        message << "Minkowski exponent p must be positive (got " << params.p << ")";
        throw ComparisonError(message.str());
    }
    // Mirrors OEFP's own bound (metric.cpp validate_tversky_parameter) exactly.
    // Written as a negated range rather than two comparisons so NaN — which
    // compares false against everything — is rejected here instead of leaking
    // out of OEFP as a std::invalid_argument.
    if (name == "tversky" &&
        (!(params.tversky_alpha >= 0.0 && params.tversky_alpha <= 1.0) ||
         !(params.tversky_beta >= 0.0 && params.tversky_beta <= 1.0))) {
        throw ComparisonError("Tversky alpha and beta must be in [0.0, 1.0]");
    }
}

}  // namespace

std::string supported_metric_names(MetricSurface surface) {
    std::ostringstream names;
    bool first = true;
    for (const MetricEntry& entry : METRIC_TABLE) {
        if (!visible_on(entry, surface)) {
            continue;
        }
        if (!first) {
            names << ", ";
        }
        names << entry.name;
        first = false;
    }
    return names.str();
}

OEFP::Metric resolve_metric(const std::string& name, bool similarity, const MetricParams& params,
                            MetricSurface surface) {
    const std::string key = to_lower(name);

    if (key == "haversine") {
        throw ComparisonError(
            "Metric 'haversine' is a great-circle distance over exactly two columns interpreted "
            "as latitude and longitude in radians; it has no molecular meaning and is not "
            "supported by OECluster");
    }

    const MetricEntry* found = nullptr;
    for (const MetricEntry& entry : METRIC_TABLE) {
        if (key == entry.name) {
            found = &entry;
            break;
        }
    }

    if (found == nullptr) {
        throw ComparisonError("Unknown metric '" + name + "'. Supported metrics: " +
                              supported_metric_names(surface));
    }

    if (!visible_on(*found, surface)) {
        if (surface == MetricSurface::Fingerprint) {
            throw ComparisonError("Metric '" + key +
                                  "' is a descriptor-space metric; use comparison=\"descriptor\"");
        }
        throw ComparisonError("Metric '" + key +
                              "' is a fingerprint bit-set metric; use comparison=\"fingerprint\"");
    }

    if (similarity && !found->has_similarity_form) {
        throw ComparisonError("Metric '" + key +
                              "' has no similarity form; only 'tanimoto' and 'tversky' support "
                              "similarity=True");
    }

    validate_params(key, params);
    return found->make(similarity, params);
}

}  // namespace OECluster
