/**
 * @file HDBSCANTree.h
 * @brief Internal condensed-tree utilities for HDBSCAN.
 */

#ifndef OECLUSTER_SRC_CLUSTERING_HDBSCANTREE_H
#define OECLUSTER_SRC_CLUSTERING_HDBSCANTREE_H

#include <cstddef>
#include <vector>

#include "HDBSCANLinkage.h"
#include "oecluster/clustering/ClusterTypes.h"
#include "oecluster/clustering/HDBSCAN.h"

namespace OECluster::detail {

struct CondensedNode {
    size_t parent = 0;
    size_t child = 0;
    double value = 0.0;
    size_t cluster_size = 0;
};

struct HDBSCANTreeSelection {
    std::vector<ClusterLabel> labels;
    std::vector<double> probabilities;
};

std::vector<CondensedNode> condense_tree(
    const std::vector<HDBSCANLinkageNode>& hierarchy,
    size_t min_cluster_size);

HDBSCANTreeSelection select_clusters(
    const std::vector<CondensedNode>& condensed_tree,
    HDBSCANClusterSelectionMethod cluster_selection_method,
    bool allow_single_cluster,
    double cluster_selection_epsilon,
    size_t max_cluster_size);

}  // namespace OECluster::detail

#endif  // OECLUSTER_SRC_CLUSTERING_HDBSCANTREE_H
