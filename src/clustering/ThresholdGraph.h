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

class NeighborRange {
public:
    NeighborRange(const size_t* begin, const size_t* end);

    const size_t* begin() const;
    const size_t* end() const;
    size_t size() const;
    bool empty() const;

private:
    const size_t* begin_;
    const size_t* end_;
};

class ThresholdNeighborGraph {
public:
    explicit ThresholdNeighborGraph(std::vector<std::vector<size_t>> neighbors);
    ThresholdNeighborGraph(std::vector<size_t> offsets, std::vector<size_t> indices);

    size_t Size() const;
    NeighborRange Neighbors(size_t index) const;

private:
    std::vector<size_t> offsets_;
    std::vector<size_t> indices_;
};

ThresholdNeighborGraph BuildThresholdNeighborGraph(
    const StorageBackend& storage,
    const ThresholdGraphOptions& options);

}  // namespace OECluster

#endif  // OECLUSTER_SRC_CLUSTERING_THRESHOLDGRAPH_H
