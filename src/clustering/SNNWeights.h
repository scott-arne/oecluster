/**
 * @file SNNWeights.h
 * @brief Shared-nearest-neighbor Jaccard weights over a KNNGraph, the input
 * graph of leiden.
 *
 * Kept apart from the optimizer so that a public SNN graph can later expose
 * this interface without touching the engine.
 */

#ifndef OECLUSTER_SRC_CLUSTERING_SNNWEIGHTS_H
#define OECLUSTER_SRC_CLUSTERING_SNNWEIGHTS_H

#include <climits>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

#include "oecluster/clustering/KNNGraph.h"

namespace OECluster::detail {

/// Symmetric weighted graph in CSR form, without self-loops.
struct WeightedGraph {
    size_t num_nodes = 0;
    std::vector<size_t> offsets;      // num_nodes + 1 entries
    std::vector<uint32_t> neighbors;  // ascending within each node
    std::vector<double> weights;      // aligned with neighbors
};

// ClusterLabel is int and the CSR indexes nodes as uint32_t, so INT_MAX
// items is the most that every singleton can be labeled and indexed.
inline void validate_leiden_item_count(size_t n, const std::string& caller) {
    const size_t limit = static_cast<size_t>(INT_MAX);
    if (n > limit) {
        throw std::invalid_argument(caller + " supports at most " +
                                    std::to_string(limit) + " items, got " +
                                    std::to_string(n));
    }
}

/**
 * @brief The SNN Jaccard graph over a KNNGraph's union edge set.
 *
 * :param graph: The neighbor graph.
 * :param prune: Edges whose weight is below this are dropped.
 * :param num_threads: Worker threads; 0 selects the hardware concurrency.
 * :returns: The weighted graph; it does not depend on num_threads.
 * :raises std::invalid_argument: If the graph has more than INT_MAX items.
 */
WeightedGraph snn_weights(const KNNGraph& graph, double prune,
                          size_t num_threads);

}  // namespace OECluster::detail

#endif  // OECLUSTER_SRC_CLUSTERING_SNNWEIGHTS_H
