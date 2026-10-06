/**
 * @file Agglomerative.h
 * @brief Hierarchical agglomerative clustering over precomputed distances, or a
 *        comparison for single linkage.
 */

#ifndef OECLUSTER_CLUSTERING_AGGLOMERATIVE_H
#define OECLUSTER_CLUSTERING_AGGLOMERATIVE_H

#include <cstddef>
#include <string>
#include <utility>
#include <vector>

#include "oecluster/PairwiseComparison.h"
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
    size_t n_clusters = 2;                        ///< Number of clusters to form; ignored if distance_threshold >= 0.
    double distance_threshold = -1.0;             ///< Stop merging above this distance; negative means use n_clusters instead.
    AgglomerativeLinkageMethod linkage = AgglomerativeLinkageMethod::Average;  ///< Linkage update method.
    bool compute_full_tree = true;                ///< Build full dendrogram even when stopping early.
    size_t num_threads = 0;                       ///< Worker threads; 0 auto-detects hardware concurrency.
    size_t chunk_size = 4096;                     ///< Chunk size for parallelized linkage updates.
};

/**
 * @brief Agglomerative clustering result with labels and merge tree metadata.
 */
class AgglomerativeResult : public ClusteringResult {
public:
    AgglomerativeResult() = default;
    AgglomerativeResult(std::vector<ClusterLabel> labels, Clusters members,
                        std::vector<size_t> children_left,
                        std::vector<size_t> children_right,
                        std::vector<double> distances,
                        std::vector<size_t> cluster_sizes)
        : ClusteringResult(std::move(labels), std::move(members)),
          children_left_(std::move(children_left)),
          children_right_(std::move(children_right)),
          distances_(std::move(distances)),
          cluster_sizes_(std::move(cluster_sizes)) {}

    /** @brief Left child node index per merge. */
    const std::vector<size_t>& ChildrenLeft() const { return children_left_; }
    /** @brief Right child node index per merge. */
    const std::vector<size_t>& ChildrenRight() const { return children_right_; }
    /** @brief Merge distance per merge. */
    const std::vector<double>& Distances() const { return distances_; }
    /** @brief Merged cluster size per merge. */
    const std::vector<size_t>& ClusterSizes() const { return cluster_sizes_; }

    std::string Method() const override { return "agglomerative"; }

private:
    std::vector<size_t> children_left_;
    std::vector<size_t> children_right_;
    std::vector<double> distances_;
    std::vector<size_t> cluster_sizes_;
};

/**
 * @brief Cluster a complete precomputed distance matrix with agglomerative clustering.
 *
 * Single linkage is built from the minimum spanning tree in O(N) memory beyond
 * the matrix; its merges at a tied height come in the tree's order. Complete,
 * average and weighted linkage use the generic heap algorithm. For single
 * linkage every distance must be finite; a zero is read as +0.0.
 *
 * :param storage: Complete pairwise distance storage.
 * :param options: Agglomerative clustering options; single linkage does not
 *     read chunk_size.
 * :returns: Labels, clusters, and merge tree metadata.
 * :raises std::invalid_argument: If options are invalid or storage is incomplete.
 * :raises std::runtime_error: If single linkage reads a non-finite distance,
 *     naming the pair.
 */
AgglomerativeResult agglomerative_cluster(
    const StorageBackend& storage,
    const AgglomerativeOptions& options);

/**
 * @brief Cluster from a comparison with single linkage, holding no matrix.
 *
 * Prim's algorithm compares every pair exactly once and holds O(N), plus one
 * comparison clone per worker. The result equals the storage overload's over
 * a matrix filled through Compare(min(i, j), max(i, j)), for every
 * num_threads, provided the comparison is repeatable: Compare(i, j) returns a
 * bit-identical value on every call and every clone, whatever that clone
 * scored before. ROCS is not, and is refused.
 *
 * :param comparison: The comparison; cloned, never called itself.
 * :param options: Agglomerative clustering options; linkage must be Single.
 * :returns: Labels, clusters, and merge tree metadata.
 * :raises std::invalid_argument: If options are invalid, linkage is not Single,
 *     or n_clusters exceeds the item count without a distance_threshold.
 * :raises ComparisonError: If the comparison is ROCS, reports similarities, a
 *     non-zero self-distance, or missing='propagate'.
 * :raises std::runtime_error: If a distance is non-finite, naming the pair.
 */
AgglomerativeResult agglomerative_cluster(
    PairwiseComparison& comparison,
    const AgglomerativeOptions& options);

}  // namespace OECluster

#endif  // OECLUSTER_CLUSTERING_AGGLOMERATIVE_H
