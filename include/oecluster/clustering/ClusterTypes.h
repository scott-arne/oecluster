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

}  // namespace OECluster

#endif  // OECLUSTER_CLUSTERING_CLUSTERTYPES_H
