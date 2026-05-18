/**
 * @file Butina.cpp
 * @brief Butina clustering implementation.
 */

#include "oecluster/clustering/Butina.h"

#include <algorithm>
#include <functional>
#include <stdexcept>
#include <utility>
#include <vector>

#include "ThresholdGraph.h"

namespace OECluster {

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

void reorder_unseen_candidates(
    std::vector<Candidate>& candidates,
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

    for (auto& candidate : candidates) {
        const size_t idx = candidate.second;
        if (!affected[idx]) {
            continue;
        }
        size_t unseen_neighbors = 0;
        for (const size_t neighbor : graph.Neighbors(idx)) {
            if (!seen[neighbor]) {
                ++unseen_neighbors;
            }
        }
        candidate.first = unseen_neighbors;
    }
    std::sort(candidates.begin(), candidates.end(), std::greater<Candidate>());
}

}  // namespace

Clusters butina_cluster(const StorageBackend& storage, const ButinaOptions& options) {
    if (options.distance_threshold < 0.0) {
        throw std::invalid_argument("Butina distance threshold must be non-negative");
    }

    if (storage.NumItems() < 2) {
        Clusters clusters;
        if (storage.NumItems() == 1) {
            clusters.push_back(Cluster{0});
        }
        return clusters;
    }

    ThresholdGraphOptions graph_options;
    graph_options.threshold = options.distance_threshold;
    graph_options.num_threads = options.num_threads;
    graph_options.chunk_size = options.chunk_size;
    const ThresholdNeighborGraph graph = BuildThresholdNeighborGraph(storage, graph_options);

    std::vector<Candidate> candidates = make_sorted_candidates(graph);
    std::vector<bool> seen(graph.Size(), false);
    Clusters clusters;
    clusters.reserve(graph.Size());

    while (!candidates.empty() && candidates.front().first > 1) {
        const size_t idx = candidates.front().second;
        candidates.erase(candidates.begin());
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

        if (options.reordering) {
            reorder_unseen_candidates(candidates, graph, seen, cluster);
        }
    }

    while (!candidates.empty()) {
        const size_t idx = candidates.front().second;
        candidates.erase(candidates.begin());
        if (!seen[idx]) {
            clusters.push_back(Cluster{idx});
            seen[idx] = true;
        }
    }

    return clusters;
}

}  // namespace OECluster
