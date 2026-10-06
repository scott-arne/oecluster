/**
 * @file ThresholdGraph.h
 * @brief Internal threshold-neighbor graph utilities for clustering algorithms.
 */

#ifndef OECLUSTER_SRC_CLUSTERING_THRESHOLDGRAPH_H
#define OECLUSTER_SRC_CLUSTERING_THRESHOLDGRAPH_H

#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

#include "oecluster/PairwiseComparison.h"
#include "oecluster/StorageBackend.h"

namespace OECluster {

struct ThresholdGraphOptions {
    double threshold = 0.0;
    size_t num_threads = 0;
    size_t chunk_size = 4096;
    // Comparison builds only: the most memory the graph may take, in bytes;
    // 0 applies detail::default_threshold_graph_limit.
    size_t max_graph_bytes = 0;
    // Comparison builds only: the entry point the diagnostics name.
    const char* caller = "threshold_graph";
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

// The storage overload's graph over a matrix filled through Compare(i, j),
// built without the matrix: every pair is compared once to size each row and
// once more to fill it. Both passes must therefore see the same value for a
// pair, so Compare must be repeatable on every call and every clone; a row
// whose passes disagree in size throws std::logic_error. A non-finite
// distance throws std::runtime_error, and a graph larger than its limit
// throws std::length_error before it is allocated.
ThresholdNeighborGraph BuildThresholdNeighborGraph(
    PairwiseComparison& comparison,
    const ThresholdGraphOptions& options);

namespace detail {

// Bytes of the graph over n items with `edges` within-threshold pairs i < j:
// n + 1 offsets and n + 2 * edges indices. Throws std::length_error when the
// size does not fit a size_t.
size_t threshold_graph_bytes(size_t n, size_t edges);

// The larger of the condensed matrix the graph replaces and 1 GiB. Throws
// std::length_error when the matrix size does not fit a size_t.
size_t default_threshold_graph_limit(size_t n);

// The limit a comparison build enforces: max_graph_bytes when it is non-zero,
// the default otherwise. The builder calls this rather than choosing itself.
size_t threshold_graph_limit(size_t n, size_t max_graph_bytes);

inline std::runtime_error non_finite_distance_error(const std::string& caller,
                                                    size_t a, size_t b) {
    return std::runtime_error(
        caller + " read a non-finite distance between items " +
        std::to_string(a) + " and " + std::to_string(b));
}

// HDBSCAN's core distance treats an item's own distance as 0, which a negative
// distance would undercut; the streaming core pass has no such zero to compare.
inline std::runtime_error negative_distance_error(const std::string& caller,
                                                  size_t a, size_t b) {
    return std::runtime_error(
        caller + " read a negative distance between items " +
        std::to_string(a) + " and " + std::to_string(b));
}

}  // namespace detail

}  // namespace OECluster

#endif  // OECLUSTER_SRC_CLUSTERING_THRESHOLDGRAPH_H
