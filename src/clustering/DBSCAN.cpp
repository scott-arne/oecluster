/**
 * @file DBSCAN.cpp
 * @brief DBSCAN clustering implementation.
 */

#include "oecluster/clustering/DBSCAN.h"

#include <stdexcept>
#include <vector>

#include "ThresholdGraph.h"

namespace OECluster {

DBSCANResult dbscan_cluster(const StorageBackend& storage, const DBSCANOptions& options) {
    if (options.eps < 0.0) {
        throw std::invalid_argument("DBSCAN eps must be non-negative");
    }
    if (options.min_samples == 0) {
        throw std::invalid_argument("DBSCAN min_samples must be at least one");
    }

    ThresholdGraphOptions graph_options;
    graph_options.threshold = options.eps;
    graph_options.num_threads = options.num_threads;
    graph_options.chunk_size = options.chunk_size;
    const ThresholdNeighborGraph graph = BuildThresholdNeighborGraph(storage, graph_options);

    std::vector<ClusterLabel> labels(graph.Size(), NOISE_LABEL);
    Cluster core_sample_indices;

    std::vector<bool> is_core(graph.Size(), false);
    for (size_t i = 0; i < graph.Size(); ++i) {
        if (graph.Neighbors(i).size() >= options.min_samples) {
            is_core[i] = true;
            core_sample_indices.push_back(i);
        }
    }

    ClusterLabel label = 0;
    std::vector<size_t> stack;
    for (size_t seed = 0; seed < graph.Size(); ++seed) {
        if (labels[seed] != NOISE_LABEL || !is_core[seed]) {
            continue;
        }

        size_t current = seed;
        while (true) {
            if (labels[current] == NOISE_LABEL) {
                labels[current] = label;
                if (is_core[current]) {
                    for (const size_t neighbor : graph.Neighbors(current)) {
                        if (labels[neighbor] == NOISE_LABEL) {
                            stack.push_back(neighbor);
                        }
                    }
                }
            }

            if (stack.empty()) {
                break;
            }
            current = stack.back();
            stack.pop_back();
        }

        ++label;
    }

    Clusters members = labels_to_clusters(labels);
    return DBSCANResult(std::move(labels), std::move(members),
                        std::move(core_sample_indices));
}

}  // namespace OECluster
