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

/// Apply the five schema-matched override rules against an already-resolved
/// selection. Other rules also govern an override -- mutual exclusion, and
/// each override's tie to its metric -- but those read nothing but the
/// options and so are settled before the schema exists.
///
/// Separated from the rules that read nothing but the options so that
/// ``validate_descriptor_options`` can show, in one glance, which half of its
/// work needs the schema resolved and which does not.
void check_override_against_selection(const OEFP::DescriptorSchema& schema,
                                      const std::vector<size_t>& selection,
                                      const DescriptorOptions& opts) {
    const std::string metric_name = to_lower(opts.metric);
    const bool is_seuclidean =
        metric_name == "standardized_euclidean" || metric_name == "seuclidean";

    // resolve_column_indices sorts into ascending schema order, but the
    // override arrives in the caller's order and is consumed positionally.
    // Equal lengths would then pair each column with the wrong number and
    // score silently wrong distances, so require the caller's own ordering to
    // already agree.
    if (!opts.columns.empty()) {
        std::vector<size_t> requested;
        requested.reserve(opts.columns.size());
        for (const std::string& name : opts.columns) {
            requested.push_back(schema.IndexOf(name));
        }
        if (!std::is_sorted(requested.begin(), requested.end())) {
            std::string expected;
            for (const size_t index : selection) {
                if (!expected.empty()) {
                    expected += ", ";
                }
                expected += schema.Definition(index).name;
            }
            throw ComparisonError(
                "columns must be in ascending schema order when variances or "
                "inverse_covariance is supplied, because those values are matched to "
                "columns by position; got the requested names in a different order. "
                "Use columns={" + expected + "}, or pass columns=stats.columns from "
                "descriptor_statistics, which is already in this order");
        }
    }

    // An override is scoped to the columns it was fitted over, so nothing is
    // dropped under it: a second drop would silently desync the selection from
    // the supplied numbers.
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
    // OEFP rejects these on the first scoring call. Checking here reports the
    // mistake where it was made, with the column named.
    for (size_t k = 0; k < opts.variances.size(); ++k) {
        if (!std::isfinite(opts.variances[k]) || opts.variances[k] <= 0.0) {
            throw ComparisonError("variances[" + std::to_string(k) + "] for column '" +
                                  schema.Definition(selection[k]).name +
                                  "' must be finite and strictly positive, got " +
                                  std::to_string(opts.variances[k]));
        }
    }
    for (size_t k = 0; k < opts.inverse_covariance.size(); ++k) {
        if (!std::isfinite(opts.inverse_covariance[k])) {
            throw ComparisonError("inverse_covariance[" + std::to_string(k) +
                                  "] must be finite, got " +
                                  std::to_string(opts.inverse_covariance[k]));
        }
    }
}

}  // namespace

void validate_descriptor_options(const DescriptorOptions& opts) {
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

    // normalize_missing validates by mapping, so calling it is what keeps the
    // policy names in one place rather than in a check and a table that drift.
    bool complete_case = false;
    const OEFP::DescriptorMissingPolicy missing = normalize_missing(opts.missing, complete_case);
    if (missing == OEFP::DescriptorMissingPolicy::Ignore && is_fitted_metric(metric_name)) {
        throw ComparisonError(
            "missing='ignore' is not available for metric='" + opts.metric +
            "'; its whitening transform mixes columns, so a per-pair column subset is "
            "incoherent. Use missing='complete_case' or an unfitted metric such as 'euclidean'");
    }

    // Only the name and the exponent are settled here; the metric itself is
    // built at the end of the constructor from the variances fitted there.
    // That is why check_metric_request exists apart from resolve_metric --
    // resolving now would succeed and hand back a metric standardized by an
    // empty vector. Pass the folded name: the constructor resolves the folded
    // spelling, and an unknown metric must be reported the same way by both.
    MetricParams params;
    params.p = opts.p;
    check_metric_request(metric_name, false, params, MetricSurface::Descriptor);

    // The remaining rules all belong to an override, and an override is matched
    // to the resolved selection by position, so they need the schema -- which
    // comes from opts.sources alone and never from a molecule. Resolving it
    // here is what lets a mistake no molecule set could rescue be refused
    // before a filter that can empty the item list. Without an override there
    // is nothing here to check, so the calculator is not built at all: that is
    // the common path and it must not start paying for a schema no rule reads.
    const bool has_override = !opts.variances.empty() || !opts.inverse_covariance.empty();
    if (!has_override) {
        return;
    }

    const std::shared_ptr<const OEFP::DescriptorCalculator> calculator =
        make_descriptor_calculator(opts.sources);
    // numeric_selection reports dropped columns through out-parameters the
    // caller is expected to keep; a verdict-only caller has nowhere to put them.
    std::vector<std::string> ignored_columns;
    std::vector<std::string> ignored_reasons;
    const std::vector<size_t> selection =
        numeric_selection(calculator->Schema(), opts, ignored_columns, ignored_reasons);
    check_override_against_selection(calculator->Schema(), selection, opts);
}

struct DescriptorComparison::Impl {
    OEFP::DescriptorNumericMatrix matrix;
    std::vector<std::string> dropped_columns;
    std::vector<std::string> dropped_reasons;
    std::vector<double> variances;
    std::vector<double> inverse_covariance;
    OEFP::DescriptorMissingPolicy missing = OEFP::DescriptorMissingPolicy::Propagate;
    bool complete_case = true;
    OEFP::Metric metric = OEFP::Metric::Euclidean();

    /// Set when a produced distance is not finite. Shared across clones through
    /// pimpl_ so the complete_case stamp cannot outlive the evidence against it,
    /// and wrapped in shared_ptr to supply mutability through the
    /// shared_ptr<const Impl> that every access path holds.
    std::shared_ptr<std::atomic<bool>> non_finite_seen =
        std::make_shared<std::atomic<bool>>(false);
};

DescriptorComparison::DescriptorComparison(const std::vector<OEChem::OEMolBase*>& mols,
                                           const Options& opts) {
    const std::vector<const OEChem::OEMolBase*> inputs =
        checked_inputs(mols, "DescriptorComparison");

    // Option rules decided with no molecules in hand, in one place.
    // validate_descriptor_options stops short of the source, column and group
    // names on a request with no override: those need a schema, and on that
    // path it builds no calculator and so resolves none. Everything else left
    // below reads the input; among those, the minimum count a fit requires,
    // the complete-case row check, the refusal when every selected column is
    // constant.
    //
    // The Python layer calls the same function before its complete-case
    // filter runs, so an unusable option value is reported as itself rather
    // than as an item list the filter emptied. When an override is supplied
    // this resolves the schema a second time, since the call above resolved
    // its own copy; that is the price of the two never disagreeing about a
    // shared rule.
    validate_descriptor_options(opts);

    auto impl = std::make_shared<Impl>();
    const std::string metric_name = to_lower(opts.metric);
    const bool has_override = !opts.variances.empty() || !opts.inverse_covariance.empty();
    const bool is_seuclidean =
        metric_name == "standardized_euclidean" || metric_name == "seuclidean";

    // The policy name was accepted above; this call is the mapping half of the
    // same function.
    impl->missing = normalize_missing(opts.missing, impl->complete_case);

    // A variance needs two observations. Without this the fit leaves every
    // variance NaN, every column is dropped, and the constructor blames the
    // descriptors for what is an input-size problem. descriptor_statistics
    // guards the same failure the same way.
    if (is_fitted_metric(metric_name) && !has_override && mols.size() < 2) {
        throw ComparisonError(
            "metric='" + metric_name + "' fits variances from the input and needs at least two "
            "molecules, got " + std::to_string(mols.size()) +
            "; pass variances= from descriptor_statistics, or use an unfitted metric such as "
            "'euclidean'");
    }

    const std::shared_ptr<const OEFP::DescriptorCalculator> calculator =
        make_descriptor_calculator(opts.sources);
    const OEFP::DescriptorSchema& schema = calculator->Schema();
    std::vector<size_t> selection =
        numeric_selection(schema, opts, impl->dropped_columns, impl->dropped_reasons);

    try {
        OEFP::DescriptorBatch batch = calculator->CalculateBatch(inputs);
        // The matrix is materialized over the selection as requested and
        // complete_case is enforced against it here, before any column is
        // dropped below. A column whose gaps are present-and-NaN rather than
        // absent has a non-finite variance -- FractionCsp3 over a set that
        // includes water, or RDKit's BCUT2D_* over a bare ion -- so the
        // zero-variance filter would remove that column and leave the offending
        // rows to pass an after-the-fact row check unchallenged. The
        // constructor would then silently score molecules that
        // descriptor_excluded_indices excludes, and the filtered and unfiltered
        // paths would disagree. Validating first is what keeps them in step.
        impl->matrix = batch.ToNumericMatrix(OEFP::DescriptorSelection::Indices(selection));
        if (impl->complete_case) {
            require_complete_rows(impl->matrix);
        }

        if (has_override) {
            // Every rule an override has to satisfy was applied by
            // validate_descriptor_options above, against the selection this
            // same schema resolves. Repeating one here would be the second copy
            // that later drifts.
            impl->variances = opts.variances;
            impl->inverse_covariance = opts.inverse_covariance;
        } else if (is_fitted_metric(metric_name)) {
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
            if (surviving.size() != selection.size()) {
                // Only rematerialize when the selection actually narrowed. The
                // rows were already validated above, so this second pass costs
                // a copy and changes no verdict.
                selection = surviving;
                impl->matrix =
                    batch.ToNumericMatrix(OEFP::DescriptorSelection::Indices(selection));
            }
            if (is_seuclidean) {
                impl->variances = variances;
            } else {
                const OEFP::DescriptorInverseCovariance inverse = OEFP::InverseCovarianceMatrix(
                    batch, OEFP::DescriptorSelection::Indices(selection));
                impl->inverse_covariance = inverse.matrix;
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
        const double value = result.empty() ? 0.0 : result[0];
        note_non_finite(&value, 1, *pimpl_->non_finite_seen);
        return value;
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
    // Observed NaN outranks the policy declaration. `ignore` skips a cell only
    // when its validity bit is clear, so a present-and-NaN value still reaches
    // the matrix under that policy; stamping SubsetScored there would downgrade
    // a hard refusal into a tier-2 one that allow_nonmetric can override.
    if (pimpl_->non_finite_seen->load()) {
        facts.data_integrity = DataIntegrity::NaNPresent;
    } else if (pimpl_->missing == OEFP::DescriptorMissingPolicy::Ignore) {
        facts.data_integrity = DataIntegrity::SubsetScored;
    } else if (!pimpl_->complete_case) {
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

    // This function is the documented way to filter an input before
    // constructing, so the options handed to it are the ones the constructor
    // will see next. Refusing them here means the caller who follows that
    // advice is told about a bad metric or a mismatched override at the call
    // where the mistake is, instead of getting a mask computed from options
    // that were never going to build.
    validate_descriptor_options(opts);

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
