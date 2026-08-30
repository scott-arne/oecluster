/**
 * @file DescriptorComparison.cpp
 * @brief Implementation of descriptor-space pairwise comparison.
 */

#include "oecluster/comparisons/DescriptorComparison.h"

#include <algorithm>
#include <atomic>
#include <cctype>
#include <cmath>
#include <exception>
#include <oechem.h>
#include <oefp/oefp.h>
#include "../DescriptorBuild.h"
#include "KernelOptions.h"
#include "MetricTable.h"
#include "oecluster/CDist.h"
#include "oecluster/CondensedIndex.h"
#include "oecluster/Error.h"
#include "oecluster/PDist.h"
#include "oecluster/StorageBackend.h"

namespace OECluster {

namespace {

std::string to_lower(const std::string& value) {
    std::string result = value;
    std::transform(result.begin(), result.end(), result.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return result;
}

/// Resolve the missing-value policy, reporting complete_case separately because
/// it uses OEFP's Propagate kernel over an input the constructor has verified.
OEFP::DescriptorMissingPolicy normalize_missing(const std::string& missing,
                                                bool& out_complete_case) {
    const std::string key = to_lower(missing);
    if (key == "complete_case") {
        out_complete_case = true;
        return OEFP::DescriptorMissingPolicy::Propagate;
    }
    out_complete_case = false;
    if (key == "propagate") {
        return OEFP::DescriptorMissingPolicy::Propagate;
    }
    if (key == "ignore") {
        return OEFP::DescriptorMissingPolicy::Ignore;
    }
    throw ComparisonError("Unknown missing-value policy: " + missing +
                          ". Supported policies are 'complete_case', 'propagate', 'ignore'");
}

bool is_fitted_metric(const std::string& metric) {
    return metric == "standardized_euclidean" || metric == "seuclidean" ||
           metric == "mahalanobis";
}

/// Resolve the numeric subset of the requested selection, recording the rest.
std::vector<size_t> numeric_selection(const OEFP::DescriptorSchema& schema,
                                      const DescriptorOptions& opts,
                                      std::vector<std::string>& dropped_columns,
                                      std::vector<std::string>& dropped_reasons) {
    const std::vector<size_t> requested =
        resolve_column_indices(schema, opts.columns, opts.groups);
    std::vector<size_t> numeric;
    for (size_t index : requested) {
        const OEFP::DescriptorDefinition& definition = schema.Definition(index);
        if (is_numeric_kind(definition.value_kind)) {
            numeric.push_back(index);
        } else {
            dropped_columns.push_back(definition.name);
            dropped_reasons.push_back("non-numeric");
        }
    }
    if (numeric.empty()) {
        throw ComparisonError("Descriptor selection resolved to zero numeric columns");
    }
    return numeric;
}

std::vector<const OEChem::OEMolBase*> checked_inputs(
        const std::vector<OEChem::OEMolBase*>& mols, const char* caller) {
    std::vector<const OEChem::OEMolBase*> inputs;
    inputs.reserve(mols.size());
    for (size_t i = 0; i < mols.size(); ++i) {
        if (mols[i] == nullptr) {
            throw ComparisonError(std::string(caller) +
                                  " received null molecule pointer at index " +
                                  std::to_string(i));
        }
        inputs.push_back(mols[i]);
    }
    return inputs;
}

/// Reject any row with an absent or non-finite value under missing='complete_case'.
void require_complete_rows(const OEFP::DescriptorNumericMatrix& matrix) {
    for (size_t row = 0; row < matrix.rows; ++row) {
        for (size_t col = 0; col < matrix.columns; ++col) {
            const size_t offset = row * matrix.columns + col;
            if (matrix.validity[offset] == 0 || !std::isfinite(matrix.values[offset])) {
                throw ComparisonError(
                    "missing='complete_case' but molecule " + std::to_string(row) +
                    " has an absent or non-finite value for descriptor '" +
                    matrix.names[col] +
                    "'; call descriptor_excluded_indices first and pass only the "
                    "retained molecules");
            }
        }
    }
}

}  // namespace

struct DescriptorComparison::Impl {
    OEFP::DescriptorNumericMatrix matrix;
    std::vector<std::string> dropped_columns;
    std::vector<std::string> dropped_reasons;
    std::vector<double> variances;
    std::vector<double> inverse_covariance;
    OEFP::DescriptorMissingPolicy missing = OEFP::DescriptorMissingPolicy::Propagate;
    bool complete_case = true;
    OEFP::Metric metric = OEFP::Metric::Euclidean();

    /// Set when a produced distance is not finite; shared across clones so the
    /// complete_case stamp cannot outlive the evidence against it.
    std::shared_ptr<std::atomic<bool>> non_finite_seen =
        std::make_shared<std::atomic<bool>>(false);
};

DescriptorComparison::DescriptorComparison(const std::vector<OEChem::OEMolBase*>& mols,
                                           const Options& opts) {
    const std::vector<const OEChem::OEMolBase*> inputs =
        checked_inputs(mols, "DescriptorComparison");

    auto impl = std::make_shared<Impl>();
    const std::string metric_name = to_lower(opts.metric);

    if (!opts.variances.empty() && !opts.inverse_covariance.empty()) {
        throw ComparisonError(
            "variances and inverse_covariance are mutually exclusive; supply the one "
            "belonging to the requested metric");
    }
    const bool is_seuclidean =
        metric_name == "standardized_euclidean" || metric_name == "seuclidean";
    if (!opts.variances.empty() && !is_seuclidean) {
        throw ComparisonError("variances applies only to metric='standardized_euclidean', not '" +
                              opts.metric + "'");
    }
    if (!opts.inverse_covariance.empty() && metric_name != "mahalanobis") {
        throw ComparisonError(
            "inverse_covariance applies only to metric='mahalanobis', not '" + opts.metric + "'");
    }

    impl->missing = normalize_missing(opts.missing, impl->complete_case);
    if (impl->missing == OEFP::DescriptorMissingPolicy::Ignore && is_fitted_metric(metric_name)) {
        throw ComparisonError(
            "missing='ignore' is not available for metric='" + opts.metric +
            "'; its whitening transform mixes columns, so a per-pair column subset is "
            "incoherent. Use missing='complete_case' or an unfitted metric such as 'euclidean'");
    }

    const std::shared_ptr<const OEFP::DescriptorCalculator> calculator =
        make_descriptor_calculator(opts.sources);
    const OEFP::DescriptorSchema& schema = calculator->Schema();
    std::vector<size_t> selection =
        numeric_selection(schema, opts, impl->dropped_columns, impl->dropped_reasons);

    OEFP::DescriptorBatch batch = calculator->CalculateBatch(inputs);
    const bool has_override = !opts.variances.empty() || !opts.inverse_covariance.empty();

    try {
        if (has_override) {
            // An override is scoped to the columns it was fitted over, so nothing
            // is dropped here: a second drop would silently desync the selection
            // from the supplied numbers. Materialize and validate before using.
            impl->matrix = batch.ToNumericMatrix(OEFP::DescriptorSelection::Indices(selection));
            if (impl->complete_case) {
                require_complete_rows(impl->matrix);
            }
            if (is_seuclidean && opts.variances.size() != selection.size()) {
                throw ComparisonError(
                    "variances has " + std::to_string(opts.variances.size()) + " entries but " +
                    std::to_string(selection.size()) +
                    " columns are selected; pass columns=stats.columns alongside "
                    "variances=stats.variance");
            }
            if (metric_name == "mahalanobis" &&
                opts.inverse_covariance.size() != selection.size() * selection.size()) {
                throw ComparisonError(
                    "inverse_covariance has " + std::to_string(opts.inverse_covariance.size()) +
                    " entries but " + std::to_string(selection.size()) +
                    " columns are selected, which needs a square " +
                    std::to_string(selection.size()) + "x" + std::to_string(selection.size()) +
                    " matrix; pass columns=stats.columns alongside "
                    "inverse_covariance=stats.inverse_covariance");
            }
            impl->variances = opts.variances;
            impl->inverse_covariance = opts.inverse_covariance;
        } else if (is_fitted_metric(metric_name)) {
            // For fitted metrics, identify zero-variance columns first so we can
            // validate only the surviving columns. A column whose gaps are
            // present-and-NaN rather than absent has a non-finite variance --
            // FractionCsp3 over a set that includes water, or RDKit's BCUT2D_*
            // over a bare ion -- and will be dropped here. Validating before the
            // drop would reject molecules that should pass, since the offending
            // column is discarded before scoring.
            const OEFP::DescriptorColumnStatistics stats =
                OEFP::ColumnStatistics(batch, OEFP::DescriptorSelection::Indices(selection));
            std::vector<size_t> surviving;
            std::vector<double> variances;
            for (size_t k = 0; k < selection.size(); ++k) {
                const double variance = stats.variance[k];
                if (!std::isfinite(variance) || variance <= 0.0) {
                    impl->dropped_columns.push_back(schema.Definition(selection[k]).name);
                    impl->dropped_reasons.push_back("zero-variance");
                    continue;
                }
                surviving.push_back(selection[k]);
                variances.push_back(variance);
            }
            if (surviving.empty()) {
                throw ComparisonError(
                    "Every selected descriptor column has zero variance over these molecules; "
                    "there is no distance to compute");
            }
            selection = surviving;
            impl->matrix =
                batch.ToNumericMatrix(OEFP::DescriptorSelection::Indices(selection));
            if (impl->complete_case) {
                require_complete_rows(impl->matrix);
            }
            if (is_seuclidean) {
                impl->variances = variances;
            } else {
                const OEFP::DescriptorInverseCovariance inverse = OEFP::InverseCovarianceMatrix(
                    batch, OEFP::DescriptorSelection::Indices(selection));
                impl->inverse_covariance = inverse.matrix;
            }
        } else {
            // Unfitted metric: no columns are dropped, so materialize and validate.
            impl->matrix = batch.ToNumericMatrix(OEFP::DescriptorSelection::Indices(selection));
            if (impl->complete_case) {
                require_complete_rows(impl->matrix);
            }
        }
    } catch (const ComparisonError&) {
        throw;
    } catch (const std::exception& exc) {
        throw ComparisonError("Failed to build descriptor comparison: " +
                              std::string(exc.what()));
    }

    MetricParams params;
    params.p = opts.p;
    params.variances = impl->variances;
    params.inverse_covariance = impl->inverse_covariance;
    impl->metric = resolve_metric(metric_name, false, params, MetricSurface::Descriptor);

    pimpl_ = std::move(impl);
}

DescriptorComparison::DescriptorComparison(std::shared_ptr<const Impl> impl)
    : pimpl_(std::move(impl)) {}

namespace {

/// Record any non-finite distance so a complete_case stamp cannot survive an
/// accumulator overflow the input check could not have predicted.
void note_non_finite(const double* values, size_t count, std::atomic<bool>& flag) {
    for (size_t i = 0; i < count; ++i) {
        if (!std::isfinite(values[i])) {
            flag.store(true);
            return;
        }
    }
}

}  // namespace

double DescriptorComparison::Compare(size_t i, size_t j) {
    const size_t columns = pimpl_->matrix.columns;
    std::vector<double> values(2 * columns);
    std::vector<std::uint8_t> validity(2 * columns);
    for (size_t col = 0; col < columns; ++col) {
        values[col] = pimpl_->matrix.values[i * columns + col];
        validity[col] = pimpl_->matrix.validity[i * columns + col];
        values[columns + col] = pimpl_->matrix.values[j * columns + col];
        validity[columns + col] = pimpl_->matrix.validity[j * columns + col];
    }

    try {
        const std::vector<double> result = OEFP::PDistNumeric(
            values.data(), validity.data(), 2, columns, pimpl_->metric, pimpl_->missing);
        return result.empty() ? 0.0 : result[0];
    } catch (const std::exception& exc) {
        throw ComparisonError("Failed to compare descriptors: " + std::string(exc.what()));
    }
}

bool DescriptorComparison::TryPDist(StorageBackend& storage, const PDistOptions& options) {
    const size_t n = pimpl_->matrix.rows;
    const size_t total_pairs = n * (n - 1) / 2;
    if (storage.NumSamples() != n) {
        throw ComparisonError("DescriptorComparison pdist storage size mismatch");
    }

    const OEFP::BatchKernelOptions kernel_options =
        make_kernel_options(options.num_threads, options.chunk_size);

    try {
        double* data = storage.Data();
        if (data != nullptr) {
            OEFP::PDistNumericInto(pimpl_->matrix.values.data(), pimpl_->matrix.validity.data(),
                                   n, pimpl_->matrix.columns, pimpl_->metric, pimpl_->missing,
                                   data, storage.NumPairs(), kernel_options);
            note_non_finite(data, storage.NumPairs(), *pimpl_->non_finite_seen);
        } else {
            const std::vector<double> values = OEFP::PDistNumeric(
                pimpl_->matrix.values.data(), pimpl_->matrix.validity.data(), n,
                pimpl_->matrix.columns, pimpl_->metric, pimpl_->missing, kernel_options);
            note_non_finite(values.data(), values.size(), *pimpl_->non_finite_seen);
            for (size_t k = 0; k < values.size(); ++k) {
                size_t i = 0;
                size_t j = 0;
                condensed_to_pair(k, n, i, j);
                storage.Set(i, j, values[k]);
            }
        }
    } catch (const std::exception& exc) {
        throw ComparisonError("Failed to compute descriptor pdist: " + std::string(exc.what()));
    }

    if (options.progress) {
        options.progress(total_pairs, total_pairs);
    }
    return true;
}

bool DescriptorComparison::TryCDist(size_t n_a, double* output, const CDistOptions& options) {
    const size_t n_total = pimpl_->matrix.rows;
    if (n_a > n_total) {
        throw ComparisonError("DescriptorComparison cdist split index is out of range");
    }
    const size_t n_b = n_total - n_a;
    const size_t total_pairs = n_a * n_b;
    if (total_pairs == 0) {
        return true;
    }

    const size_t columns = pimpl_->matrix.columns;
    const double* a_values = pimpl_->matrix.values.data();
    const std::uint8_t* a_validity = pimpl_->matrix.validity.data();

    const OEFP::BatchKernelOptions kernel_options =
        make_kernel_options(options.num_threads, options.chunk_size);

    try {
        OEFP::CDistNumericInto(a_values, a_validity, n_a, a_values + n_a * columns,
                               a_validity + n_a * columns, n_b, columns, pimpl_->metric,
                               pimpl_->missing, output, total_pairs, kernel_options);
        note_non_finite(output, total_pairs, *pimpl_->non_finite_seen);
    } catch (const std::exception& exc) {
        throw ComparisonError("Failed to compute descriptor cdist: " + std::string(exc.what()));
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

std::unique_ptr<PairwiseComparison> DescriptorComparison::Clone() const {
    return std::unique_ptr<PairwiseComparison>(new DescriptorComparison(pimpl_));
}

size_t DescriptorComparison::Size() const {
    return pimpl_->matrix.rows;
}

std::string DescriptorComparison::ComparisonName() const {
    return "descriptor";
}

GateFacts DescriptorComparison::Facts() const {
    GateFacts facts;
    // OEFP's Type() is the authoritative orientation; HasZeroSelfDistance()
    // documents that it answers about the returned number, not the intent, so
    // the two facts must be read from two different accessors.
    facts.is_distance = (pimpl_->metric.Type() == OEFP::MetricType::Distance)
                            ? Capability::Yes
                            : Capability::No;
    facts.zero_self =
        pimpl_->metric.HasZeroSelfDistance() ? Capability::Yes : Capability::No;
    facts.triangle =
        pimpl_->metric.SatisfiesTriangleInequality() ? Capability::Yes : Capability::No;
    if (pimpl_->missing == OEFP::DescriptorMissingPolicy::Ignore) {
        facts.data_integrity = DataIntegrity::SubsetScored;
    } else if (!pimpl_->complete_case || pimpl_->non_finite_seen->load()) {
        facts.data_integrity = DataIntegrity::NaNPresent;
    } else {
        facts.data_integrity = DataIntegrity::Complete;
    }
    return facts;
}

const std::vector<std::string>& DescriptorComparison::Columns() const {
    return pimpl_->matrix.names;
}

const std::vector<std::string>& DescriptorComparison::DroppedColumns() const {
    return pimpl_->dropped_columns;
}

const std::vector<std::string>& DescriptorComparison::DroppedReasons() const {
    return pimpl_->dropped_reasons;
}

const std::vector<double>& DescriptorComparison::Variances() const {
    return pimpl_->variances;
}

const std::vector<double>& DescriptorComparison::InverseCovariance() const {
    return pimpl_->inverse_covariance;
}

std::vector<size_t> descriptor_excluded_indices(const std::vector<OEChem::OEMolBase*>& mols,
                                                const DescriptorOptions& opts) {
    const std::vector<const OEChem::OEMolBase*> inputs =
        checked_inputs(mols, "descriptor_excluded_indices");

    const std::shared_ptr<const OEFP::DescriptorCalculator> calculator =
        make_descriptor_calculator(opts.sources);

    std::vector<std::string> ignored_columns;
    std::vector<std::string> ignored_reasons;
    const std::vector<size_t> selection =
        numeric_selection(calculator->Schema(), opts, ignored_columns, ignored_reasons);

    OEFP::DescriptorNumericMatrix matrix;
    try {
        matrix = calculator->CalculateBatch(inputs).ToNumericMatrix(
            OEFP::DescriptorSelection::Indices(selection));
    } catch (const std::exception& exc) {
        throw ComparisonError("Failed to compute descriptors for complete-case filtering: " +
                              std::string(exc.what()));
    }

    std::vector<size_t> excluded;
    for (size_t row = 0; row < matrix.rows; ++row) {
        for (size_t col = 0; col < matrix.columns; ++col) {
            const size_t offset = row * matrix.columns + col;
            if (matrix.validity[offset] == 0 || !std::isfinite(matrix.values[offset])) {
                excluded.push_back(row);
                break;
            }
        }
    }
    return excluded;
}

}  // namespace OECluster
