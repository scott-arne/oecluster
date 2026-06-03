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
 * @brief Butina clustering result. Carries only labels and members; the first
 *     member of each cluster is the highest-neighborhood representative.
 */
class ButinaResult : public ClusteringResult {
public:
    using ClusteringResult::ClusteringResult;

    std::string Method() const override { return "butina"; }
};

/**
 * @brief Cluster a precomputed distance matrix with the Butina algorithm.
 *
 * :param storage: Pairwise distance storage.
 * :param options: Butina clustering options.
 * :returns: ButinaResult whose clusters are ordered with the
 *     highest-neighborhood representative first, and whose per-item labels
 *     equal each member's cluster position.
 */
ButinaResult butina_cluster(const StorageBackend& storage, const ButinaOptions& options);

}  // namespace OECluster

#endif  // OECLUSTER_CLUSTERING_BUTINA_H
