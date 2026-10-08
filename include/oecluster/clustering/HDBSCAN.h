/**
 * @file HDBSCAN.h
 * @brief HDBSCAN clustering over a precomputed distance matrix or a comparison.
 */

#ifndef OECLUSTER_CLUSTERING_HDBSCAN_H
#define OECLUSTER_CLUSTERING_HDBSCAN_H

#include <cstddef>
#include <string>
#include <utility>
#include <vector>

#include "oecluster/PairwiseComparison.h"
#include "oecluster/StorageBackend.h"
#include "oecluster/clustering/ClusterTypes.h"

namespace OECluster {

/**
 * @brief Cluster selection strategy for HDBSCAN condensed trees.
 */
enum class HDBSCANClusterSelectionMethod {
    EOM,
    Leaf
};

/**
 * @brief Options for HDBSCAN clustering.
 */
struct HDBSCANOptions {
    size_t min_cluster_size = 5;
    size_t min_samples = 0;
    double cluster_selection_epsilon = 0.0;
    size_t max_cluster_size = 0;
    double alpha = 1.0;
    HDBSCANClusterSelectionMethod cluster_selection_method =
        HDBSCANClusterSelectionMethod::EOM;
    bool allow_single_cluster = false;
    size_t num_threads = 0;  ///< Worker threads; 0 auto-detects hardware concurrency.
    size_t chunk_size = 64;  ///< Ceiling on rows per work unit in the matrix core-distance pass; 0 selects 64.
};

/**
 * @brief HDBSCAN result with labels, clusters, and membership probabilities.
 */
class HDBSCANResult : public ClusteringResult {
public:
    HDBSCANResult() = default;
    HDBSCANResult(std::vector<ClusterLabel> labels, Clusters members,
                  std::vector<double> probabilities)
        : ClusteringResult(std::move(labels), std::move(members)),
          probabilities_(std::move(probabilities)) {}

    /** @brief Per-item membership probabilities. */
    const std::vector<double>& Probabilities() const { return probabilities_; }

    std::string Method() const override { return "hdbscan"; }

private:
    std::vector<double> probabilities_;
};

/**
 * @brief Cluster a precomputed distance matrix with HDBSCAN.
 *
 * Every distance must be finite and non-negative; a zero is read as +0.0.
 *
 * :param storage: Complete pairwise distance storage.
 * :param options: HDBSCAN clustering options.
 * :returns: Labels, clusters, and probabilities.
 * :raises std::invalid_argument: If min_cluster_size < 2, alpha is not positive
 *     (NaN included), min_samples > NumSamples(), storage is incomplete, or a
 *     distance divided by alpha is not finite.
 * :raises std::runtime_error: If a distance is non-finite or negative, naming the pair.
 */
HDBSCANResult hdbscan_cluster(const StorageBackend& storage, const HDBSCANOptions& options);

/**
 * @brief Cluster from a comparison with HDBSCAN, holding no matrix.
 *
 * One pass compares every pair once for the core distances, holding
 * N x (min_samples - 1) doubles; a pruned Prim pass then compares at most every
 * pair again, holding O(N). With min_samples = 1 there is only the Prim pass.
 * One comparison clone per worker is held on top. The result equals the
 * storage overload's over a matrix filled through Compare(min(i, j), max(i, j)),
 * for every num_threads, provided the comparison is repeatable: Compare(i, j)
 * returns a bit-identical value on every call and every clone, whatever that
 * clone scored before. ROCS is not, and is refused.
 *
 * :param comparison: The comparison; cloned, never called itself.
 * :param options: HDBSCAN clustering options; chunk_size is not read.
 * :returns: Labels, clusters, and probabilities.
 * :raises std::invalid_argument: As the storage overload, for the options and
 *     the min_samples bound, and if a distance divided by alpha is not finite.
 * :raises ComparisonError: If the comparison is ROCS, reports similarities, a
 *     non-zero self-distance, or missing='propagate'.
 * :raises std::runtime_error: If a distance is non-finite or negative, naming the pair.
 * :raises std::length_error: If the core-distance heaps' size overflows a size_t.
 */
HDBSCANResult hdbscan_cluster(PairwiseComparison& comparison, const HDBSCANOptions& options);

namespace detail {

/**
 * @brief Each item's distance to its (min_samples - 1)-th nearest neighbor.
 *
 * Since 5.20.0 a wrapper on the pass HDBSCAN itself runs, so the distance
 * domain it enforces applies here too: a non-finite or negative distance is
 * refused rather than carried into the result, and a -0.0 is read as +0.0.
 * Earlier releases returned whatever the matrix held.
 *
 * :param storage: Complete pairwise distances; sparse storage is refused.
 * :param min_samples: Self-inclusive neighbor count, in [1, item count].
 *     A value of 1 gives every item a core distance of zero.
 * :param num_threads: Worker threads; 0 auto-detects hardware concurrency.
 * :returns: One core distance per item.
 * :raises std::invalid_argument: If min_samples is zero or above the item
 *     count, or the storage cannot provide complete distances.
 * :raises std::runtime_error: If a distance is non-finite or negative,
 *     naming the pair.
 */
std::vector<double> compute_core_distances(
    const StorageBackend& storage,
    size_t min_samples,
    size_t num_threads);

}  // namespace detail

}  // namespace OECluster

#endif  // OECLUSTER_CLUSTERING_HDBSCAN_H
