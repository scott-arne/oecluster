/**
 * @file PartitionAgreement.cpp
 * @brief Agreement metrics between two labelings of the same samples.
 */

#include "oecluster/clustering/PartitionAgreement.h"

#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "ContingencyTable.h"

namespace OECluster {

namespace {

constexpr double UNDEFINED = std::numeric_limits<double>::quiet_NaN();

/// C(n, 2) in uint64_t. Exact for any n below 2^32: C(2^32, 2) is about
/// 9.2e18, inside uint64_t's 1.8e19 range.
uint64_t choose_two(uint64_t n) { return n < 2 ? 0 : n * (n - 1) / 2; }

/**
 * @brief The two sides group the samples identically, however they numbered
 *        them.
 *
 * One cell per row and one per column, each filling both its marginals, is a
 * bijection between the two sides' clusters. O(nnz), and independent of the
 * label values, so two all-singleton labelings pass whatever integers or
 * strings they used.
 */
bool partitions_identical(const detail::ContingencyTable& table) {
    if (table.cells.size() != table.marginals_a.size() ||
        table.cells.size() != table.marginals_b.size()) {
        return false;
    }
    for (const detail::ContingencyTable::Cell& cell : table.cells) {
        if (cell.count != table.marginals_a[cell.row] ||
            cell.count != table.marginals_b[cell.col]) {
            return false;
        }
    }
    return true;
}

PartitionAgreement score(const detail::ContingencyTable& table,
                         const PartitionAgreementOptions& options) {
    PartitionAgreement agreement;
    agreement.num_samples = static_cast<size_t>(table.num_samples);
    agreement.num_clusters_a = table.marginals_a.size();
    agreement.num_clusters_b = table.marginals_b.size();
    agreement.requested.adjusted_mutual_information =
        options.compute_adjusted_mutual_information;

    // Rule 1. Agreement between two labelings of one sample is not a
    // measurement, and C(1, 2) == 0 makes every pair denominator vanish.
    // Outranks rule 2, so a lone survivor reports NaN rather than 1.0.
    if (table.num_samples < 2) {
        return agreement;
    }

    // Rule 2. Outranks rule 3, which is what keeps ARI's denominator off zero:
    // the only inputs that zero it are identical partitions.
    if (partitions_identical(table)) {
        agreement.adjusted_rand_index = 1.0;
        agreement.fowlkes_mallows = 1.0;
        agreement.normalized_mutual_information = 1.0;
        agreement.homogeneity = 1.0;
        agreement.completeness = 1.0;
        agreement.v_measure = 1.0;
        if (options.compute_adjusted_mutual_information) {
            agreement.adjusted_mutual_information = 1.0;
        }
        return agreement;
    }

    // Rule 3: compute, and report NaN for any quantity whose denominator is
    // zero.
    uint64_t sum_cells = 0;
    for (const detail::ContingencyTable::Cell& cell : table.cells) {
        sum_cells += choose_two(cell.count);
    }
    uint64_t sum_a = 0;
    for (uint64_t marginal : table.marginals_a) {
        sum_a += choose_two(marginal);
    }
    uint64_t sum_b = 0;
    for (uint64_t marginal : table.marginals_b) {
        sum_b += choose_two(marginal);
    }
    const uint64_t total = choose_two(table.num_samples);

    const double expected = static_cast<double>(sum_a) *
                            static_cast<double>(sum_b) /
                            static_cast<double>(total);
    const double maximum =
        0.5 * (static_cast<double>(sum_a) + static_cast<double>(sum_b));
    agreement.adjusted_rand_index =
        (static_cast<double>(sum_cells) - expected) / (maximum - expected);

    // 0/0 when either side is all singletons: no coincident pair exists to
    // normalize against.
    agreement.fowlkes_mallows =
        (sum_a == 0 || sum_b == 0)
            ? UNDEFINED
            : static_cast<double>(sum_cells) /
                  std::sqrt(static_cast<double>(sum_a) *
                            static_cast<double>(sum_b));

    return agreement;
}

void validate_labelings(size_t size_a, size_t size_b) {
    if (size_a == 0 || size_b == 0) {
        throw std::invalid_argument(
            "partition_agreement requires non-empty labelings");
    }
    if (size_a != size_b) {
        throw std::invalid_argument(
            "b has " + std::to_string(size_b) + " entries but a has " +
            std::to_string(size_a) + " samples");
    }
}

void validate_scaffolds(size_t num_samples, size_t num_scaffolds) {
    if (num_samples == 0 || num_scaffolds == 0) {
        throw std::invalid_argument(
            "scaffold_agreement requires a non-empty clustering and a "
            "non-empty scaffold annotation");
    }
    if (num_samples != num_scaffolds) {
        throw std::invalid_argument(
            "scaffold_labels has " + std::to_string(num_scaffolds) +
            " entries but the clustering has " + std::to_string(num_samples) +
            " samples");
    }
}

}  // namespace

PartitionAgreement partition_agreement(
    const std::vector<ClusterLabel>& a, const std::vector<ClusterLabel>& b,
    const PartitionAgreementOptions& options) {
    validate_labelings(a.size(), b.size());
    return score(detail::build_contingency(a, b, options.noise_handling),
                 options);
}

PartitionAgreement partition_agreement(
    const ClusteringResult& a, const ClusteringResult& b,
    const PartitionAgreementOptions& options) {
    return partition_agreement(a.Labels(), b.Labels(), options);
}

PartitionAgreement scaffold_agreement(
    const std::vector<ClusterLabel>& labels,
    const std::vector<std::string>& scaffold_labels,
    const PartitionAgreementOptions& options) {
    validate_scaffolds(labels.size(), scaffold_labels.size());
    return score(
        detail::build_contingency(labels, scaffold_labels, options.noise_handling),
        options);
}

PartitionAgreement scaffold_agreement(
    const ClusteringResult& result,
    const std::vector<std::string>& scaffold_labels,
    const PartitionAgreementOptions& options) {
    return scaffold_agreement(result.Labels(), scaffold_labels, options);
}

}  // namespace OECluster
