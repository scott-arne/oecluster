/**
 * @file Agglomerative.h
 * @brief Hierarchical agglomerative clustering over precomputed distances.
 */

#ifndef OECLUSTER_CLUSTERING_AGGLOMERATIVE_H
#define OECLUSTER_CLUSTERING_AGGLOMERATIVE_H

#include <cstddef>
#include <vector>

#include "oecluster/StorageBackend.h"
#include "oecluster/clustering/ClusterTypes.h"

namespace OECluster {

/**
 * @brief Linkage update method for agglomerative clustering.
 */
enum class AgglomerativeLinkageMethod {
    Single,
    Complete,
    Average,
    Weighted
};

/**
 * @brief Options for hierarchical agglomerative clustering.
 */
struct AgglomerativeOptions {
    size_t n_clusters = 2;
    double distance_threshold = -1.0;
    AgglomerativeLinkageMethod linkage = AgglomerativeLinkageMethod::Average;
    bool compute_full_tree = true;
    size_t num_threads = 0;
    size_t chunk_size = 4096;
};

/**
 * @brief Agglomerative clustering result with labels and merge tree metadata.
 */
struct AgglomerativeResult : public ClusteringResult {
    std::vector<size_t> children_left;
    std::vector<size_t> children_right;
    std::vector<double> distances;
    std::vector<size_t> cluster_sizes;
};

/**
 * @brief Cluster a complete precomputed distance matrix with agglomerative clustering.
 *
 * :param storage: Complete pairwise distance storage.
 * :param options: Agglomerative clustering options.
 * :returns: Labels, clusters, and merge tree metadata.
 * :raises std::invalid_argument: If options are invalid or storage is incomplete.
 */
AgglomerativeResult agglomerative_cluster(
    const StorageBackend& storage,
    const AgglomerativeOptions& options);

}  // namespace OECluster

#endif  // OECLUSTER_CLUSTERING_AGGLOMERATIVE_H
