/**
 * @file SphereExclusionEngine.h
 * @brief The sphere-exclusion engine shared by sphere_exclusion and Butina.
 */

#ifndef OECLUSTER_SRC_CLUSTERING_SPHEREEXCLUSIONENGINE_H
#define OECLUSTER_SRC_CLUSTERING_SPHEREEXCLUSIONENGINE_H

#include <cstddef>
#include <utility>
#include <vector>

#include "ThresholdGraph.h"
#include "oecluster/clustering/ClusterTypes.h"

namespace OECluster::detail {

// Clusters in center order, each center first; labels are cluster positions
// and centers[c] == clusters[c][0].
struct SphereEngineResult {
    std::vector<ClusterLabel> labels;
    Clusters clusters;
    std::vector<size_t> centers;
};

inline SphereEngineResult sphere_result_from_clusters(size_t n,
                                                      Clusters clusters) {
    SphereEngineResult result;
    result.labels.assign(n, NOISE_LABEL);
    result.centers.reserve(clusters.size());
    for (size_t c = 0; c < clusters.size(); ++c) {
        result.centers.push_back(clusters[c].front());
        for (const size_t member : clusters[c]) {
            result.labels[member] = static_cast<ClusterLabel>(c);
        }
    }
    result.clusters = std::move(clusters);
    return result;
}

// Butina's seed order over a threshold graph: descending neighbor count,
// larger index first on a tie, with optional reordering of the unclaimed
// counts after each cluster. Every item is assigned.
SphereEngineResult sphere_neighbors_first(const ThresholdNeighborGraph& graph,
                                          bool reordering);

}  // namespace OECluster::detail

#endif  // OECLUSTER_SRC_CLUSTERING_SPHEREEXCLUSIONENGINE_H
