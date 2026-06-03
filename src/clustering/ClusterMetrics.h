/**
 * @file ClusterMetrics.h
 * @brief Shared per-cluster distance statistics used by representative
 *        selection and clustering-quality reporting.
 */

#ifndef OECLUSTER_SRC_CLUSTERING_CLUSTERMETRICS_H
#define OECLUSTER_SRC_CLUSTERING_CLUSTERMETRICS_H

#include <cstddef>
#include <vector>

#include "oecluster/StorageBackend.h"
#include "oecluster/clustering/ClusterTypes.h"

namespace OECluster::detail {

/** @brief Arithmetic mean of a distance list; 0.0 when empty. */
double mean_distance(const std::vector<double>& distances);

/** @brief Maximum of a distance list; 0.0 when empty. */
double max_distance(const std::vector<double>& distances);

/** @brief Median of a distance list (takes by value to sort); 0.0 when empty. */
double median_distance(std::vector<double> distances);

/** @brief Maximum pairwise distance within a cluster; 0.0 for size < 2. */
double cluster_diameter(const Cluster& cluster, const StorageBackend& storage);

}  // namespace OECluster::detail

#endif  // OECLUSTER_SRC_CLUSTERING_CLUSTERMETRICS_H
