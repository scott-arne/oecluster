/**
 * @file Butina.h
 * @brief Butina clustering over precomputed distance matrices.
 */

#ifndef OECLUSTER_CLUSTERING_BUTINA_H
#define OECLUSTER_CLUSTERING_BUTINA_H

#include <cstddef>

#include "oecluster/StorageBackend.h"
#include "oecluster/clustering/ClusterTypes.h"

namespace OECluster {

/**
 * @brief Options for Butina clustering.
 */
struct ButinaOptions {
    double distance_threshold = 0.0;
    bool reordering = false;
    size_t num_threads = 0;
    size_t chunk_size = 4096;
};

/**
 * @brief Cluster a precomputed distance matrix with the Butina algorithm.
 *
 * :param storage: Pairwise distance storage.
 * :param options: Butina clustering options.
 * :returns: Clusters where the first member is the highest-neighborhood representative.
 */
Clusters butina_cluster(const StorageBackend& storage, const ButinaOptions& options);

}  // namespace OECluster

#endif  // OECLUSTER_CLUSTERING_BUTINA_H
