/**
 * @file HDBSCAN.h
 * @brief HDBSCAN clustering over precomputed distance matrices.
 */

#ifndef OECLUSTER_CLUSTERING_HDBSCAN_H
#define OECLUSTER_CLUSTERING_HDBSCAN_H

#include <cstddef>
#include <utility>
#include <vector>

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
    size_t num_threads = 0;
    size_t chunk_size = 4096;
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
 * :param storage: Complete pairwise distance storage.
 * :param options: HDBSCAN clustering options.
 * :returns: Labels, clusters, and probabilities.
 */
HDBSCANResult hdbscan_cluster(const StorageBackend& storage, const HDBSCANOptions& options);

namespace detail {

std::vector<double> compute_core_distances(
    const StorageBackend& storage,
    size_t min_samples,
    size_t num_threads);

}  // namespace detail

}  // namespace OECluster

#endif  // OECLUSTER_CLUSTERING_HDBSCAN_H
