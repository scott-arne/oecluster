/**
 * @file SphereExclusion.cpp
 * @brief Sphere exclusion (leader, Butina and DISE orders) over distance
 * matrices and comparisons.
 */

#include <algorithm>
#include <cstddef>
#include <functional>
#include <utility>
#include <vector>

#include "SphereExclusionEngine.h"
#include "ThresholdGraph.h"

namespace OECluster::detail {

namespace {

using Candidate = std::pair<size_t, size_t>;  // (neighbor_count, index)

std::vector<Candidate> make_sorted_candidates(const ThresholdNeighborGraph& graph) {
    std::vector<Candidate> candidates;
    candidates.reserve(graph.Size());
    for (size_t i = 0; i < graph.Size(); ++i) {
        candidates.emplace_back(graph.Neighbors(i).size(), i);
    }
    std::sort(candidates.begin(), candidates.end(), std::greater<Candidate>());
    return candidates;
}

// After forming a cluster, recompute unseen-neighbor counts for affected
// candidates to prioritize high-connectivity nodes in remaining items. Only
// the range from the cursor on is live: the consumed prefix is never read
// again, so updating and sorting just that range gives the order the old
// erase-then-sort-everything loop gave, without its quadratic cost.
void reorder_unseen_candidates(std::vector<Candidate>& candidates,
                               size_t cursor,
                               const ThresholdNeighborGraph& graph,
                               const std::vector<bool>& seen,
                               const Cluster& new_cluster) {
    std::vector<bool> affected(graph.Size(), false);
    for (const size_t member : new_cluster) {
        for (const size_t neighbor : graph.Neighbors(member)) {
            if (!seen[neighbor]) {
                affected[neighbor] = true;
            }
        }
    }

    const auto live = candidates.begin() + static_cast<std::ptrdiff_t>(cursor);
    for (auto it = live; it != candidates.end(); ++it) {
        const size_t idx = it->second;
        if (!affected[idx]) {
            continue;
        }
        size_t unseen_neighbors = 0;
        for (const size_t neighbor : graph.Neighbors(idx)) {
            if (!seen[neighbor]) {
                ++unseen_neighbors;
            }
        }
        it->first = unseen_neighbors;
    }
    std::sort(live, candidates.end(), std::greater<Candidate>());
}

}  // namespace

SphereEngineResult sphere_neighbors_first(const ThresholdNeighborGraph& graph,
                                          bool reordering) {
    std::vector<Candidate> candidates = make_sorted_candidates(graph);
    std::vector<bool> seen(graph.Size(), false);
    Clusters clusters;
    clusters.reserve(graph.Size());

    size_t cursor = 0;
    while (cursor < candidates.size() && candidates[cursor].first > 1) {
        const size_t idx = candidates[cursor].second;
        ++cursor;
        if (seen[idx]) {
            continue;
        }

        Cluster cluster;
        cluster.reserve(graph.Neighbors(idx).size());
        cluster.push_back(idx);
        seen[idx] = true;

        for (const size_t neighbor : graph.Neighbors(idx)) {
            if (!seen[neighbor]) {
                cluster.push_back(neighbor);
                seen[neighbor] = true;
            }
        }
        clusters.push_back(cluster);

        if (reordering) {
            reorder_unseen_candidates(candidates, cursor, graph, seen, cluster);
        }
    }

    for (; cursor < candidates.size(); ++cursor) {
        const size_t idx = candidates[cursor].second;
        if (!seen[idx]) {
            clusters.push_back(Cluster{idx});
            seen[idx] = true;
        }
    }

    return sphere_result_from_clusters(graph.Size(), std::move(clusters));
}

}  // namespace OECluster::detail
