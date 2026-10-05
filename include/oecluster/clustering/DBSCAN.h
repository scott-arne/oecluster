/**
 * @file DBSCAN.h
 * @brief DBSCAN clustering over a distance matrix or a comparison.
 */

#ifndef OECLUSTER_CLUSTERING_DBSCAN_H
#define OECLUSTER_CLUSTERING_DBSCAN_H

#include <cstddef>
#include <string>
#include <utility>
#include <vector>

#include "oecluster/PairwiseComparison.h"
#include "oecluster/StorageBackend.h"
#include "oecluster/clustering/ClusterTypes.h"

namespace OECluster {

/**
 * @brief Options for DBSCAN clustering.
 */
struct DBSCANOptions {
    double eps = 0.5;
    size_t min_samples = 5;
    size_t num_threads = 0;
    size_t chunk_size = 4096;
    /// Comparison overload only: the most memory the threshold graph may
    /// take, in bytes. 0 applies the default limit, the larger of the
    /// condensed matrix the graph replaces and 1 GiB. The storage overload
    /// refuses any other value.
    size_t max_graph_bytes = 0;
};

/**
 * @brief DBSCAN result with labels, clusters, and core sample indices.
 */
class DBSCANResult : public ClusteringResult {
public:
    DBSCANResult() = default;
    DBSCANResult(std::vector<ClusterLabel> labels, Clusters members,
                 Cluster core_sample_indices)
        : ClusteringResult(std::move(labels), std::move(members)),
          core_sample_indices_(std::move(core_sample_indices)) {}

    /** @brief Indices of core samples. */
    const Cluster& CoreSampleIndices() const { return core_sample_indices_; }

    std::string Method() const override { return "dbscan"; }

private:
    Cluster core_sample_indices_;
};

/**
 * @brief Cluster a precomputed distance matrix with DBSCAN.
 *
 * :param storage: Pairwise distance storage.
 * :param options: DBSCAN clustering options.
 * :returns: Labels, clusters, and core sample indices.
 * :raises std::invalid_argument: If options are invalid, max_graph_bytes is
 *     non-zero, or sparse storage is incomplete.
 */
DBSCANResult dbscan_cluster(const StorageBackend& storage, const DBSCANOptions& options);

/**
 * @brief Cluster a comparison with DBSCAN, holding no matrix.
 *
 * Builds the eps-neighbor graph in two passes over every pair -- one to count
 * each item's neighbors, one to record them -- so the result equals the
 * storage overload's on a matrix filled through Compare(i, j), for every
 * num_threads and chunk_size. The graph costs 16 bytes per within-eps pair
 * on a 64-bit platform, and its exact size is known before it is allocated,
 * so a graph above the limit is refused rather than attempted. A chunk_size
 * of 0 selects 4096, as on the storage overload.
 *
 * Precondition: Compare(i, j) returns a bit-identical value for a pair on
 * every call and every clone. A comparison that changes a row's neighbor
 * count between the passes is refused; one that swaps neighbors while
 * keeping every count is outside the contract.
 *
 * :param comparison: Distance comparison; cloned once per running chunk.
 * :param options: DBSCAN clustering options.
 * :returns: Labels, clusters, and core sample indices.
 * :raises std::invalid_argument: If eps is negative or min_samples is zero.
 * :raises ComparisonError: If the comparison reports similarities, a
 *     non-zero self-distance, or values that may be NaN.
 * :raises std::runtime_error: If a comparison returns NaN or infinity.
 * :raises std::length_error: If the graph would exceed its limit.
 * :raises std::logic_error: If the two passes disagree on a row's size.
 */
DBSCANResult dbscan_cluster(PairwiseComparison& comparison, const DBSCANOptions& options);

}  // namespace OECluster

#endif  // OECLUSTER_CLUSTERING_DBSCAN_H
