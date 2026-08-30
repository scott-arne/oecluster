/**
 * @file DescriptorStatistics.h
 * @brief Fitted per-column statistics over a molecular descriptor batch.
 */

#ifndef OECLUSTER_DESCRIPTORSTATISTICS_H
#define OECLUSTER_DESCRIPTORSTATISTICS_H

#include <cstddef>
#include <string>
#include <vector>

namespace OEChem {
class OEMolBase;
}

namespace OECluster {

/**
 * @brief Which descriptors to fit statistics over.
 */
struct DescriptorStatisticsOptions {
    std::vector<std::string> sources;  ///< Source names; empty selects {"openeye"}.
    std::vector<std::string> columns;  ///< Explicit column names; empty selects all.
    std::vector<std::string> groups;   ///< Group names, unioned with ``columns``.
    bool inverse_covariance = false;   ///< Also fit the Mahalanobis inverse covariance.
};

/**
 * @brief Fitted statistics and the columns that were discarded.
 *
 * The five per-column vectors describe only the surviving columns, in the order
 * given by ``columns``. ``dropped_columns`` and ``dropped_reasons`` are parallel
 * and describe what was removed and why.
 */
struct DescriptorStatisticsResult {
    std::vector<std::string> columns;         ///< Surviving column names, schema order.
    std::vector<double> mean;                 ///< Mean of the present values.
    std::vector<double> variance;             ///< Sample variance, denominator n - 1.
    std::vector<double> minimum;              ///< Smallest present value.
    std::vector<double> maximum;              ///< Largest present value.
    std::vector<size_t> present_count;        ///< Present values per column.
    std::vector<double> inverse_covariance;   ///< Row-major k x k; empty unless requested.
    size_t inverse_covariance_rank = 0;       ///< Retained eigenvalues; 0 unless requested.
    std::vector<std::string> dropped_columns; ///< Discarded column names.
    std::vector<std::string> dropped_reasons; ///< "zero-variance" or "non-numeric".
    size_t num_rows = 0;                      ///< Molecules the fit ran over.
};

/**
 * @brief Fit descriptor statistics over a set of molecules.
 *
 * Non-numeric columns and columns whose variance is zero or undefined are
 * discarded and reported: a zero-variance column contributes nothing to a
 * distance and makes standardized Euclidean divide by zero.
 *
 * :param mols: Molecules to describe. Must not contain null pointers.
 * :param options: Source, column, group, and inverse-covariance selection.
 * :returns: The fitted statistics and the drop report.
 * :raises ComparisonError: When a source, column, or group is unknown, a
 *     molecule pointer is null, or OEFP cannot compute the requested statistics.
 */
DescriptorStatisticsResult descriptor_statistics(const std::vector<OEChem::OEMolBase*>& mols,
                                                 const DescriptorStatisticsOptions& options);

}  // namespace OECluster

#endif  // OECLUSTER_DESCRIPTORSTATISTICS_H
