/**
 * @file DBSCAN.cpp
 * @brief DBSCAN clustering implementation.
 */

#include "oecluster/clustering/DBSCAN.h"

#include <stdexcept>
#include <vector>

#include "DiversityValidation.h"
#include "ThresholdGraph.h"

namespace OECluster {

namespace {

constexpr const char* DBSCAN_NAME = "dbscan";

void validate_dbscan_options(const DBSCANOptions& options) {
    if (options.eps < 0.0) {
        throw std::invalid_argument("DBSCAN eps must be non-negative");
    }
    if (options.min_samples == 0) {
        throw std::invalid_argument("DBSCAN min_samples must be at least one");
    }
}

// Everything after the graph is shared, so the two overloads cannot drift
// apart once the graphs agree.
DBSCANResult dbscan_from_graph(const ThresholdNeighborGraph& graph,
                               size_t min_samples) {
    std::vector<ClusterLabel> labels(graph.Size(), NOISE_LABEL);
    Cluster core_sample_indices;

    std::vector<bool> is_core(graph.Size(), false);
    for (size_t i = 0; i < graph.Size(); ++i) {
        if (graph.Neighbors(i).size() >= min_samples) {
            is_core[i] = true;
            core_sample_indices.push_back(i);
        }
    }

    // Iterative stack-based region growing avoids stack overflow on large clusters.
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

}  // namespace

DBSCANResult dbscan_cluster(const StorageBackend& storage, const DBSCANOptions& options) {
    validate_dbscan_options(options);
    // Refused rather than ignored: a budget the matrix path never reads would
    // read as protection that is not there.
    if (options.max_graph_bytes != 0) {
        throw std::invalid_argument(
            "dbscan max_graph_bytes applies only when clustering from a comparison");
    }

    ThresholdGraphOptions graph_options;
    graph_options.threshold = options.eps;
    graph_options.num_threads = options.num_threads;
    graph_options.chunk_size = options.chunk_size;
    return dbscan_from_graph(BuildThresholdNeighborGraph(storage, graph_options),
                             options.min_samples);
}

DBSCANResult dbscan_cluster(PairwiseComparison& comparison, const DBSCANOptions& options) {
    validate_dbscan_options(options);
    detail::validate_distance_facts(comparison, DBSCAN_NAME);

    ThresholdGraphOptions graph_options;
    graph_options.threshold = options.eps;
    graph_options.num_threads = options.num_threads;
    graph_options.chunk_size = options.chunk_size;
    graph_options.max_graph_bytes = options.max_graph_bytes;
    graph_options.caller = DBSCAN_NAME;
    return dbscan_from_graph(BuildThresholdNeighborGraph(comparison, graph_options),
                             options.min_samples);
}

}  // namespace OECluster
