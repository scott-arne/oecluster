/**
 * @file DBSCAN.h
 * @brief DBSCAN clustering over precomputed distance matrices.
 */

#ifndef OECLUSTER_CLUSTERING_DBSCAN_H
#define OECLUSTER_CLUSTERING_DBSCAN_H

#include <cstddef>

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

private:
    Cluster core_sample_indices_;
};

/**
 * @brief Cluster a precomputed distance matrix with DBSCAN.
 *
 * :param storage: Pairwise distance storage.
 * :param options: DBSCAN clustering options.
 * :returns: Labels, clusters, and core sample indices.
 * :raises std::invalid_argument: If options are invalid or sparse storage is incomplete.
 */
DBSCANResult dbscan_cluster(const StorageBackend& storage, const DBSCANOptions& options);

}  // namespace OECluster

#endif  // OECLUSTER_CLUSTERING_DBSCAN_H
