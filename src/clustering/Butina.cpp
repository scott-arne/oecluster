/**
 * @file Butina.cpp
 * @brief Butina clustering implementation.
 */

#include "oecluster/clustering/Butina.h"

#include <stdexcept>
#include <utility>
#include <vector>

#include "DiversityValidation.h"
#include "SphereExclusionEngine.h"
#include "ThresholdGraph.h"

namespace OECluster {

namespace {

constexpr const char* BUTINA_NAME = "butina";

void validate_butina_threshold(const ButinaOptions& options) {
    if (options.distance_threshold < 0.0) {
        throw std::invalid_argument("Butina distance threshold must be non-negative");
    }
}

ButinaResult small_butina_result(size_t n) {
    if (n == 1) {
        return ButinaResult(std::vector<ClusterLabel>{0}, Clusters{Cluster{0}});
    }
    return ButinaResult();
}

ButinaResult butina_from_graph(const ThresholdNeighborGraph& graph, bool reordering) {
    detail::SphereEngineResult result =
        detail::sphere_neighbors_first(graph, reordering);
    return ButinaResult(std::move(result.labels), std::move(result.clusters));
}

}  // namespace

// Butina is sphere exclusion under the neighbor-count order. It keeps its own
// validation rather than routing through sphere_exclusion(), so it gains none
// of that function's refusals and its outputs stay exactly as they were.
ButinaResult butina_cluster(const StorageBackend& storage, const ButinaOptions& options) {
    validate_butina_threshold(options);
    // Refused rather than ignored: a budget the matrix path never reads would
    // read as protection that is not there.
    if (options.max_graph_bytes != 0) {
        throw std::invalid_argument(
            "butina max_graph_bytes applies only when clustering from a comparison");
    }

    if (storage.NumSamples() < 2) {
        return small_butina_result(storage.NumSamples());
    }

    ThresholdGraphOptions graph_options;
    graph_options.threshold = options.distance_threshold;
    graph_options.num_threads = options.num_threads;
    graph_options.chunk_size = options.chunk_size;
    return butina_from_graph(BuildThresholdNeighborGraph(storage, graph_options),
                             options.reordering);
}

ButinaResult butina_cluster(PairwiseComparison& comparison, const ButinaOptions& options) {
    validate_butina_threshold(options);
    detail::validate_distance_facts(comparison, BUTINA_NAME);

    // No small-input shortcut, unlike the storage overload: the builder's
    // trivial graph below two items still answers to an explicit budget, and
    // the engine turns it into the same result the shortcut would.
    ThresholdGraphOptions graph_options;
    graph_options.threshold = options.distance_threshold;
    graph_options.num_threads = options.num_threads;
    graph_options.chunk_size = options.chunk_size;
    graph_options.max_graph_bytes = options.max_graph_bytes;
    graph_options.caller = BUTINA_NAME;
    return butina_from_graph(BuildThresholdNeighborGraph(comparison, graph_options),
                             options.reordering);
}

}  // namespace OECluster
