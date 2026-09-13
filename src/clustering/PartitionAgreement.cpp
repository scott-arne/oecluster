/**
 * @file PartitionAgreement.cpp
 * @brief Agreement metrics between two labelings of the same samples.
 */

#include "oecluster/clustering/PartitionAgreement.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

#include "ContingencyTable.h"

namespace OECluster {

namespace {

constexpr double UNDEFINED = std::numeric_limits<double>::quiet_NaN();

/// C(n, 2) in uint64_t. Exact for any n below 2^32: C(2^32, 2) is about
/// 9.2e18, inside uint64_t's 1.8e19 range.
uint64_t choose_two(uint64_t n) { return n < 2 ? 0 : n * (n - 1) / 2; }

/// log(k!) for k in [0, n], by prefix-summing log(k). Costs 8*(n+1) bytes --
/// 800 KB at n = 100k -- and turns each inner term into nine table lookups and
/// eight adds instead of nine lgamma calls.
std::vector<double> log_factorials(uint64_t n) {
    std::vector<double> table(static_cast<size_t>(n) + 1, 0.0);
    for (uint64_t k = 2; k <= n; ++k) {
        table[static_cast<size_t>(k)] =
            table[static_cast<size_t>(k - 1)] + std::log(static_cast<double>(k));
    }
    return table;
}

/// The inner sum of Vinh et al. (2010) for one pair of marginal values. It
/// depends on (i, j) only through (u, w), which is what makes the grouping in
/// expected_mutual_information exact rather than approximate.
double expected_term(uint64_t u, uint64_t w, uint64_t n,
                     const std::vector<double>& logfact) {
    const auto lf = [&logfact](uint64_t k) {
        return logfact[static_cast<size_t>(k)];
    };
    const double total = static_cast<double>(n);
    // The hypergeometric support: a cell cannot be emptier than u + w - n, and
    // count = 0 contributes nothing because 0 * log(...) is 0.
    const uint64_t lower = (u + w > n) ? (u + w - n) : 1;
    const uint64_t upper = std::min(u, w);
    double sum = 0.0;
    for (uint64_t count = lower; count <= upper; ++count) {
        const double log_p = lf(u) + lf(w) + lf(n - u) + lf(n - w) - lf(n) -
                             lf(count) - lf(u - count) - lf(w - count) -
                             lf(n - u - w + count);
        sum += (static_cast<double>(count) / total) *
               std::log(total * static_cast<double>(count) /
                        (static_cast<double>(u) * static_cast<double>(w))) *
               std::exp(log_p);
    }
    return sum;
}

/**
 * @brief Expected mutual information under the hypergeometric model.
 *
 * The outer sums run over every (i, j) pair of marginals, including pairs
 * whose observed cell count is zero, so this cannot ride the sparse cell list.
 * Grouping equal marginal values is exact, not an approximation, and bounds the
 * work by the number of distinct cluster sizes rather than the cluster count --
 * two all-singleton sides collapse from N^2 pairs to one. This sum needs no
 * extra canonicalization: the histograms are built from the multiset of
 * cluster sizes and traversed in ascending order, both of which are unchanged
 * by a permutation of the samples.
 */
double expected_mutual_information(const detail::ContingencyTable& table,
                                   const std::vector<double>& logfact) {
    double expected = 0.0;
    detail::for_each_marginal_pair(
        detail::marginal_histogram(table.marginals_a),
        detail::marginal_histogram(table.marginals_b),
        [&](uint64_t u, uint64_t w, uint64_t multiplicity) {
            expected += static_cast<double>(multiplicity) *
                        expected_term(u, w, table.num_samples, logfact);
        });
    return expected;
}

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

    // Natural log throughout. Both sums accumulate in an order determined by
    // the values being added, never by the interned ids: ids follow first
    // appearance, so permuting the samples renumbers the clusters, and a sum
    // taken in id order would change in its last bits. Sorting by a key that
    // exactly determines each term makes equal keys mean equal terms, so the
    // emitted sequence -- and the sum -- is the same for every sample order
    // and invariant under argument swap. The MI term sort uses (min, max) of
    // the marginals rather than (row, col), making it transpose-invariant:
    // each term depends only on count and the product row*col, and (min, max,
    // count) determines both. The pair metrics above need no such treatment:
    // they are exact integer arithmetic. Cost is O(K log K + nnz log nnz),
    // the same order as building the table.
    const double n = static_cast<double>(table.num_samples);
    const auto entropy_of = [n](const std::vector<uint64_t>& marginals) {
        std::vector<uint64_t> sizes(marginals);
        std::sort(sizes.begin(), sizes.end());
        double entropy = 0.0;
        for (uint64_t size : sizes) {
            const double p = static_cast<double>(size) / n;
            if (p > 0.0) {
                entropy -= p * std::log(p);
            }
        }
        return entropy;
    };
    const double entropy_a = entropy_of(table.marginals_a);
    const double entropy_b = entropy_of(table.marginals_b);

    struct MutualInformationTerm {
        uint64_t row_size;
        uint64_t col_size;
        uint64_t count;
    };
    std::vector<MutualInformationTerm> terms;
    terms.reserve(table.cells.size());
    for (const detail::ContingencyTable::Cell& cell : table.cells) {
        terms.push_back(MutualInformationTerm{table.marginals_a[cell.row],
                                              table.marginals_b[cell.col],
                                              cell.count});
    }
    std::sort(terms.begin(), terms.end(),
              [](const MutualInformationTerm& lhs,
                 const MutualInformationTerm& rhs) {
                  const uint64_t lhs_min = std::min(lhs.row_size, lhs.col_size);
                  const uint64_t lhs_max = std::max(lhs.row_size, lhs.col_size);
                  const uint64_t rhs_min = std::min(rhs.row_size, rhs.col_size);
                  const uint64_t rhs_max = std::max(rhs.row_size, rhs.col_size);
                  return std::tie(lhs_min, lhs_max, lhs.count) <
                         std::tie(rhs_min, rhs_max, rhs.count);
              });

    double mutual_information = 0.0;
    for (const MutualInformationTerm& term : terms) {
        const double count = static_cast<double>(term.count);
        const double row = static_cast<double>(term.row_size);
        const double col = static_cast<double>(term.col_size);
        mutual_information += (count / n) * std::log((count * n) / (row * col));
    }

    // 2*MI/(H(a)+H(b)) is the definition, not the harmonic mean of the two
    // components below: the harmonic form is 0/0 both when MI is zero with two
    // positive entropies and when one entropy is zero, and the composite is
    // well defined in each case. Rules 1 and 2 have already removed the only
    // input whose denominator here is zero -- two zero-entropy sides are two
    // single-cluster partitions, which are identical -- so neither field is
    // ever NaN.
    const double shared = 2.0 * mutual_information / (entropy_a + entropy_b);
    agreement.normalized_mutual_information = shared;
    agreement.v_measure = shared;
    agreement.homogeneity =
        entropy_a > 0.0 ? mutual_information / entropy_a : UNDEFINED;
    agreement.completeness =
        entropy_b > 0.0 ? mutual_information / entropy_b : UNDEFINED;

    if (options.compute_adjusted_mutual_information) {
        const std::vector<double> logfact = log_factorials(table.num_samples);
        const double chance = expected_mutual_information(table, logfact);
        double denominator = 0.5 * (entropy_a + entropy_b) - chance;

        // The one place here where a convention replaces a NaN, matching
        // scikit-learn. Defensive rather than reachable: E[MI] never exceeds
        // min(H(a), H(b)), so the denominator is at least 0.5*|H(a) - H(b)|,
        // and over every input rule 2 does not already intercept its minimum
        // is log(2)/N -- about 3e15 samples short of epsilon.
        const double eps = std::numeric_limits<double>::epsilon();
        denominator = denominator < 0.0 ? std::min(denominator, -eps)
                                        : std::max(denominator, eps);
        agreement.adjusted_mutual_information =
            (mutual_information - chance) / denominator;
    }

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
