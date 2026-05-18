/**
 * @file Centroid.h
 * @brief Representative member selection for clusters.
 */

#ifndef OECLUSTER_CLUSTERING_CENTROID_H
#define OECLUSTER_CLUSTERING_CENTROID_H

#include <cstddef>

#include "oecluster/StorageBackend.h"
#include "oecluster/clustering/ClusterTypes.h"

namespace OECluster {

/**
 * @brief Algorithms for choosing a representative cluster member.
 */
enum class CentroidMethod {
    First,    ///< Return the first member stored in the cluster.
    Medoid,   ///< Return the member with the lowest total in-cluster distance.
    Minimax,  ///< Return the member with the lowest maximum in-cluster distance.
    Mean      ///< Alias for medoid over a precomputed distance matrix.
};

/**
 * @brief Select a representative member index from a cluster.
 *
 * :param cluster: Cluster member indices.
 * :param storage: Pairwise distance storage used for distance-based methods.
 * :param method: Centroid selection algorithm.
 * :returns: Member index selected as the cluster representative.
 * :raises std::invalid_argument: If the cluster is empty or sparse distances are incomplete.
 * :raises std::out_of_range: If a cluster member is outside the storage range.
 */
size_t cluster_centroid(
    const Cluster& cluster,
    const StorageBackend& storage,
    CentroidMethod method = CentroidMethod::Medoid);

}  // namespace OECluster

#endif  // OECLUSTER_CLUSTERING_CENTROID_H
