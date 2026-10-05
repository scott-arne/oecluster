/**
 * @file SphereExclusionEngine.h
 * @brief The sphere-exclusion engine shared by sphere_exclusion and Butina.
 */

#ifndef OECLUSTER_SRC_CLUSTERING_SPHEREEXCLUSIONENGINE_H
#define OECLUSTER_SRC_CLUSTERING_SPHEREEXCLUSIONENGINE_H

#include <cmath>
#include <cstddef>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
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

// The nearest center to item, ties to the earlier center; best is a position
// in centers. A non-finite read stops the scan and is reported through
// best_distance, so the caller raises on it in a fixed order.
template <typename Read>
void nearest_center(size_t item, const std::vector<size_t>& centers,
                    Read&& read, size_t& best, double& best_distance) {
    best = 0;
    best_distance = std::numeric_limits<double>::infinity();
    for (size_t k = 0; k < centers.size(); ++k) {
        const double distance = read(centers[k], item);
        if (!std::isfinite(distance)) {
            best = k;
            best_distance = distance;
            return;
        }
        if (distance < best_distance) {
            best = k;
            best_distance = distance;
        }
    }
}

// First-claim sphere exclusion over input order (permutation null) or the
// given permutation. The source fills distances for the still-unclaimed
// items; claims are committed here, serially and in ascending item order,
// so a parallel source never touches the claim state.
template <typename Source>
SphereEngineResult sphere_ordered_first(Source& source, size_t n,
                                        const std::vector<size_t>* permutation,
                                        double threshold) {
    std::vector<bool> claimed(n, false);
    std::vector<size_t> unclaimed(n);
    std::iota(unclaimed.begin(), unclaimed.end(), size_t{0});
    std::vector<size_t> candidates;
    std::vector<size_t> remaining;
    std::vector<double> distances;
    Clusters clusters;
    for (size_t step = 0; step < n; ++step) {
        const size_t center =
            permutation == nullptr ? step : (*permutation)[step];
        if (claimed[center]) {
            continue;
        }
        claimed[center] = true;
        candidates.clear();
        for (const size_t item : unclaimed) {
            if (item != center) {
                candidates.push_back(item);
            }
        }
        source.Row(center, candidates, distances);

        Cluster cluster{center};
        remaining.clear();
        for (size_t p = 0; p < candidates.size(); ++p) {
            if (!std::isfinite(distances[p])) {
                throw non_finite_distance_error("sphere_exclusion", center,
                                                candidates[p]);
            }
            if (distances[p] <= threshold) {
                cluster.push_back(candidates[p]);
                claimed[candidates[p]] = true;
            } else {
                remaining.push_back(candidates[p]);
            }
        }
        unclaimed.swap(remaining);
        clusters.push_back(std::move(cluster));
    }
    return sphere_result_from_clusters(n, std::move(clusters));
}

// Reassigns every non-center item to its nearest center; the centers, and
// so the clusters' order, stay as the first-claim pass fixed them.
template <typename Source>
void sphere_assign_nearest(Source& source, size_t n,
                           SphereEngineResult& result) {
    const std::vector<size_t>& centers = result.centers;
    std::vector<bool> is_center(n, false);
    for (const size_t center : centers) {
        is_center[center] = true;
    }
    std::vector<size_t> others;
    for (size_t item = 0; item < n; ++item) {
        if (!is_center[item]) {
            others.push_back(item);
        }
    }
    std::vector<size_t> best;
    std::vector<double> best_distance;
    source.Nearest(others, centers, best, best_distance);

    Clusters clusters(centers.size());
    for (size_t k = 0; k < centers.size(); ++k) {
        clusters[k].push_back(centers[k]);
    }
    for (size_t p = 0; p < others.size(); ++p) {
        if (!std::isfinite(best_distance[p])) {
            throw non_finite_distance_error("sphere_exclusion",
                                            centers[best[p]], others[p]);
        }
        clusters[best[p]].push_back(others[p]);
    }
    result = sphere_result_from_clusters(n, std::move(clusters));
}

}  // namespace OECluster::detail

#endif  // OECLUSTER_SRC_CLUSTERING_SPHEREEXCLUSIONENGINE_H
