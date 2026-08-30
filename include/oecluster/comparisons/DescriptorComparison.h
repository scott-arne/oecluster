/**
 * @file DescriptorComparison.h
 * @brief Pairwise comparison over numeric molecular descriptors.
 */

#ifndef OECLUSTER_COMPARISONS_DESCRIPTORCOMPARISON_H
#define OECLUSTER_COMPARISONS_DESCRIPTORCOMPARISON_H

#include <cstddef>
#include <memory>
#include <string>
#include <vector>
#include "oecluster/PairwiseComparison.h"

namespace OEChem { class OEMolBase; }

namespace OECluster {

/**
 * @brief Descriptor source, column selection, metric, and missing-value policy.
 *
 * ``variances`` and ``inverse_covariance`` are mutually exclusive overrides for
 * the two fitted metrics. Supplying one skips fitting and skips the
 * zero-variance drop, so the selection is used exactly as given.
 */
struct DescriptorOptions {
    std::vector<std::string> sources;  ///< Empty selects {"openeye"}.
    std::vector<std::string> columns;  ///< Explicit column names; empty selects all.
    std::vector<std::string> groups;   ///< Group names, unioned with ``columns``.
    std::string metric = "standardized_euclidean";
    std::vector<double> variances;          ///< standardized_euclidean override.
    std::vector<double> inverse_covariance; ///< mahalanobis override, row-major square.
    std::string missing = "complete_case";  ///< complete_case | propagate | ignore.
    double p = 2.0;                         ///< Minkowski order.
};

/**
 * @brief Descriptor-space pairwise comparison backed by OEFP numeric kernels.
 *
 * Descriptors are computed once during construction and shared across clones
 * through ``std::shared_ptr``, matching FingerprintComparison.
 *
 * Under ``missing="complete_case"`` the constructor requires every selected
 * value to be present and finite; call ``descriptor_excluded_indices`` first and
 * pass only the retained molecules. Filtering deliberately lives outside the
 * comparison so that ``Size()`` always equals the input count, which is what the
 * cdist split index is measured against.
 *
 * Rows are checked against the full requested selection *before* any
 * zero-variance column is dropped. A column whose gaps are present-and-NaN
 * rather than absent has a non-finite variance and would be dropped by the fit,
 * so validating afterwards would let the constructor silently score molecules
 * that ``descriptor_excluded_indices`` excludes. The two must not disagree.
 *
 * The integrity stamp escalates to ``NaNPresent`` whenever a non-finite
 * distance is actually produced, whatever the policy -- including under
 * ``ignore``, where the declared stamp would otherwise be ``SubsetScored``.
 */
class DescriptorComparison : public PairwiseComparison {
public:
    using Options = DescriptorOptions;

    /**
     * @brief Construct a DescriptorComparison from a set of molecules.
     *
     * :param mols: Pointers to molecules (not owned). Must not contain nulls.
     * :param opts: Descriptor options.
     * :raises ComparisonError: When a source, column, group, metric, or missing
     *     policy is invalid; when an override does not match the selection; when
     *     every selected column is constant; or when ``complete_case`` is
     *     requested and a selected value is absent or non-finite.
     */
    explicit DescriptorComparison(const std::vector<OEChem::OEMolBase*>& mols,
                                  const Options& opts = Options());

    double Compare(size_t i, size_t j) override;
    bool TryPDist(StorageBackend& storage, const PDistOptions& options) override;
    bool TryCDist(size_t n_a, double* output, const CDistOptions& options) override;
    std::unique_ptr<PairwiseComparison> Clone() const override;
    size_t Size() const override;
    std::string ComparisonName() const override;
    GateFacts Facts() const override;

    /// Surviving column names, in ascending schema order regardless of the
    /// order they were requested in. Values supplied through ``variances`` or
    /// ``inverse_covariance`` are matched to these by position.
    const std::vector<std::string>& Columns() const;
    /// Columns removed before scoring.
    const std::vector<std::string>& DroppedColumns() const;
    /// Parallel to DroppedColumns: "zero-variance" or "non-numeric".
    /// "zero-variance" covers a column variance that is zero, negative, or
    /// not finite -- a source that emits a present NaN yields the last case.
    const std::vector<std::string>& DroppedReasons() const;
    /// Fitted or supplied variances; empty unless the metric is standardized_euclidean.
    const std::vector<double>& Variances() const;
    /// Fitted or supplied inverse covariance; empty unless the metric is mahalanobis.
    const std::vector<double>& InverseCovariance() const;

private:
    struct Impl;
    std::shared_ptr<const Impl> pimpl_;

    /// Private clone constructor -- shares immutable descriptor data.
    explicit DescriptorComparison(std::shared_ptr<const Impl> impl);
};

/**
 * @brief Indices of molecules that complete-case filtering would drop.
 *
 * A molecule is excluded when any selected descriptor value is absent or not
 * finite. The selection is resolved exactly as ``DescriptorComparison`` would
 * resolve it, before any zero-variance drop, so filtering and scoring agree.
 *
 * :param mols: Pointers to molecules (not owned). Must not contain nulls.
 * :param opts: Descriptor options; only the source, column, and group fields
 *     affect the result.
 * :returns: Ascending indices into *mols*.
 * :raises ComparisonError: When a source, column, or group is unknown, or a
 *     molecule pointer is null.
 */
std::vector<size_t> descriptor_excluded_indices(const std::vector<OEChem::OEMolBase*>& mols,
                                                const DescriptorOptions& opts);

}  // namespace OECluster

#endif  // OECLUSTER_COMPARISONS_DESCRIPTORCOMPARISON_H
