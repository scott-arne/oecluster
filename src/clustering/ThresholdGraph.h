/**
 * @file ThresholdGraph.h
 * @brief Internal threshold-neighbor graph utilities for clustering algorithms.
 */

#ifndef OECLUSTER_SRC_CLUSTERING_THRESHOLDGRAPH_H
#define OECLUSTER_SRC_CLUSTERING_THRESHOLDGRAPH_H

#include <cstddef>
#include <vector>

#include "oecluster/StorageBackend.h"

namespace OECluster {

struct ThresholdGraphOptions {
    double threshold = 0.0;
    size_t num_threads = 0;
    size_t chunk_size = 4096;
};

class ThresholdNeighborGraph {
public:
    explicit ThresholdNeighborGraph(std::vector<std::vector<size_t>> neighbors);

    size_t Size() const;
    const std::vector<size_t>& Neighbors(size_t index) const;
    std::vector<size_t>& MutableNeighbors(size_t index);

private:
    std::vector<std::vector<size_t>> neighbors_;
};

ThresholdNeighborGraph BuildThresholdNeighborGraph(
    const StorageBackend& storage,
    const ThresholdGraphOptions& options);

}  // namespace OECluster

#endif  // OECLUSTER_SRC_CLUSTERING_THRESHOLDGRAPH_H
