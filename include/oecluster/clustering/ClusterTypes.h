/**
 * @file ClusterTypes.h
 * @brief Shared clustering result types.
 */

#ifndef OECLUSTER_CLUSTERING_CLUSTERTYPES_H
#define OECLUSTER_CLUSTERING_CLUSTERTYPES_H

#include <cstddef>
#include <vector>

namespace OECluster {

using Cluster = std::vector<size_t>;
using Clusters = std::vector<Cluster>;
using ClusterLabel = int;

constexpr ClusterLabel NOISE_LABEL = -1;

/**
 * @brief Shared clustering result with integer labels and grouped members.
 *
 * Labels use -1 for noise. Non-negative labels identify clusters and are
 * converted to member lists by labels_to_clusters().
 */
struct ClusteringResult {
    std::vector<ClusterLabel> labels;
    Clusters clusters;
};

/**
 * @brief Convert cluster labels into cluster member lists.
 *
 * :param labels: Per-item labels, where negative values are ignored as noise.
 * :returns: Clusters ordered by non-negative label, preserving item order.
 */
Clusters labels_to_clusters(const std::vector<ClusterLabel>& labels);

}  // namespace OECluster

#endif  // OECLUSTER_CLUSTERING_CLUSTERTYPES_H
