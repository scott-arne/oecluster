/**
 * @file Butina.cpp
 * @brief Butina clustering implementation.
 */

#include "oecluster/clustering/Butina.h"

#include <stdexcept>
#include <utility>
#include <vector>

#include "SphereExclusionEngine.h"
#include "ThresholdGraph.h"

namespace OECluster {

// Butina is sphere exclusion under the neighbor-count order. It keeps its own
// validation rather than routing through sphere_exclusion(), so it gains none
// of that function's refusals and its outputs stay exactly as they were.
ButinaResult butina_cluster(const StorageBackend& storage, const ButinaOptions& options) {
    if (options.distance_threshold < 0.0) {
        throw std::invalid_argument("Butina distance threshold must be non-negative");
    }

    if (storage.NumSamples() < 2) {
        if (storage.NumSamples() == 1) {
            return ButinaResult(std::vector<ClusterLabel>{0}, Clusters{Cluster{0}});
        }
        return ButinaResult();
    }

    ThresholdGraphOptions graph_options;
    graph_options.threshold = options.distance_threshold;
    graph_options.num_threads = options.num_threads;
    graph_options.chunk_size = options.chunk_size;
    const ThresholdNeighborGraph graph = BuildThresholdNeighborGraph(storage, graph_options);

    detail::SphereEngineResult result =
        detail::sphere_neighbors_first(graph, options.reordering);
    return ButinaResult(std::move(result.labels), std::move(result.clusters));
}

}  // namespace OECluster
