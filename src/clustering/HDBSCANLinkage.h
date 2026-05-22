/**
 * @file HDBSCANLinkage.h
 * @brief Internal MST and single-linkage utilities for HDBSCAN.
 */

#ifndef OECLUSTER_SRC_CLUSTERING_HDBSCANLINKAGE_H
#define OECLUSTER_SRC_CLUSTERING_HDBSCANLINKAGE_H

#include <cstddef>
#include <vector>

#include "oecluster/StorageBackend.h"

namespace OECluster::detail {

struct HDBSCANMSTEdge {
    size_t current_node = 0;
    size_t next_node = 0;
    double distance = 0.0;
};

struct HDBSCANLinkageNode {
    size_t left_node = 0;
    size_t right_node = 0;
    double value = 0.0;
    size_t cluster_size = 0;
};

std::vector<HDBSCANMSTEdge> hdbscan_mutual_reachability_mst(
    const StorageBackend& storage,
    const std::vector<double>& core_distances,
    double alpha);

std::vector<HDBSCANLinkageNode> make_hdbscan_single_linkage(
    std::vector<HDBSCANMSTEdge> mst,
    size_t n_samples);

}  // namespace OECluster::detail

#endif  // OECLUSTER_SRC_CLUSTERING_HDBSCANLINKAGE_H
