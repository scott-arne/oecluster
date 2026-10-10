/**
 * @file agglomerative_oracle.cpp
 * @brief The 5.20.0 heap linkage algorithm, frozen as a differential oracle.
 */

#include "agglomerative_oracle.h"

#include <algorithm>
#include <cstddef>
#include <limits>
#include <queue>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "oecluster/ThreadPool.h"

#include "../../src/clustering/AgglomerativeInternal.h"
#include "../../src/clustering/DistanceAccess.h"

namespace agglomerative_oracle {
namespace {

using OECluster::AgglomerativeLinkageMethod;
using OECluster::AgglomerativeOptions;
using OECluster::AgglomerativeResult;
using OECluster::ClusterLabel;
using OECluster::Clusters;
using OECluster::labels_to_clusters;
using OECluster::StorageBackend;
using OECluster::ThreadPool;
using OECluster::detail::labels_from_cut;
using OECluster::detail::max_node_count;

// The frozen code below is a verbatim lift from OECluster's anonymous
// namespace, where `detail::` resolved to OECluster::detail. The alias keeps
// it verbatim rather than rewriting every call site.
namespace detail = OECluster::detail;

struct MergeCandidate {
    double distance = 0.0;
    size_t left = 0;
    size_t right = 0;
};

struct MergeCandidateGreater {
    bool operator()(const MergeCandidate& lhs, const MergeCandidate& rhs) const {
        if (lhs.distance != rhs.distance) {
            return lhs.distance > rhs.distance;
        }
        if (lhs.left != rhs.left) {
            return lhs.left > rhs.left;
        }
        return lhs.right > rhs.right;
    }
};


double& cluster_distance(
    std::vector<double>& distances,
    size_t max_nodes,
    size_t left,
    size_t right) {
    return distances[detail::condensed_index(max_nodes, left, right)];
}

// Canonical ordering (left < right) ensures consistent heap key structure for the priority queue.
MergeCandidate make_candidate(double distance, size_t left, size_t right) {
    if (right < left) {
        std::swap(left, right);
    }
    return MergeCandidate{distance, left, right};
}

// Average linkage weighs by cluster sizes; Weighted uses unweighted 0.5 factor

double update_linkage_distance(
    AgglomerativeLinkageMethod linkage,
    double left_distance,
    double right_distance,
    size_t left_size,
    size_t right_size) {
    switch (linkage) {
        case AgglomerativeLinkageMethod::Single:
            return std::min(left_distance, right_distance);
        case AgglomerativeLinkageMethod::Complete:
            return std::max(left_distance, right_distance);
        case AgglomerativeLinkageMethod::Average:
            return ((static_cast<double>(left_size) * left_distance) +
                    (static_cast<double>(right_size) * right_distance)) /
                   static_cast<double>(left_size + right_size);
        case AgglomerativeLinkageMethod::Weighted:
            return 0.5 * (left_distance + right_distance);
    }

    throw std::invalid_argument("Unknown agglomerative linkage method");
}


void validate_arguments(const AgglomerativeOptions& options) {
    if (options.chunk_size == 0) {
        throw std::invalid_argument("Agglomerative chunk_size must be at least one");
    }
    if (std::isnan(options.distance_threshold)) {
        throw std::invalid_argument("Agglomerative distance_threshold must not be NaN");
    }
    if (options.distance_threshold < 0.0 && options.n_clusters == 0) {
        throw std::invalid_argument("Agglomerative n_clusters must be at least one");
    }
}

// A threshold cut ignores n_clusters, so the bound applies only without one.
void validate_cluster_bound(const AgglomerativeOptions& options, size_t n) {
    if (options.distance_threshold < 0.0 && options.n_clusters > n) {
        throw std::invalid_argument(
            "Agglomerative n_clusters must be at most the item count");
    }
}

void validate_options(
    const StorageBackend& storage,
    const AgglomerativeOptions& options) {
    detail::validate_complete_distance_storage(storage, "Agglomerative clustering");
    validate_arguments(options);
    validate_cluster_bound(options, storage.NumSamples());
}

AgglomerativeResult small_result(size_t n) {
    if (n == 0) {
        return AgglomerativeResult();
    }
    std::vector<ClusterLabel> labels{0};
    Clusters members = labels_to_clusters(labels);
    return AgglomerativeResult(std::move(labels), std::move(members),
                               {}, {}, {}, {});
}

std::vector<double> initialize_cluster_distances(
    const StorageBackend& storage,
    size_t max_nodes,
    size_t num_threads,
    size_t chunk_size) {
    const size_t n = storage.NumSamples();
    std::vector<double> distances(max_nodes * (max_nodes - 1) / 2,
                                  std::numeric_limits<double>::infinity());
    const double* data = storage.Data();

    ThreadPool pool(num_threads);
    pool.ParallelFor(0, n, chunk_size, [&](size_t begin, size_t end) {
        for (size_t i = begin; i < end; ++i) {
            for (size_t j = i + 1; j < n; ++j) {
                cluster_distance(distances, max_nodes, i, j) =
                    detail::dense_distance(data, n, i, j);
            }
        }
    });

    return distances;
}

}  // namespace

AgglomerativeResult heap_cluster(
    const StorageBackend& storage,
    const AgglomerativeOptions& options) {
    validate_options(storage, options);

    const size_t n = storage.NumSamples();
    if (n < 2) {
        return small_result(n);
    }

    const bool cut_by_threshold = options.distance_threshold >= 0.0;
    const size_t target_merges =
        (!options.compute_full_tree && !cut_by_threshold)
            ? n - options.n_clusters
            : n - 1;

    const size_t max_nodes = max_node_count(n);
    std::vector<double> distances = initialize_cluster_distances(
        storage,
        max_nodes,
        options.num_threads,
        options.chunk_size);

    std::vector<size_t> cluster_sizes(max_nodes, 0);
    std::vector<bool> active(max_nodes, false);
    for (size_t i = 0; i < n; ++i) {
        cluster_sizes[i] = 1;
        active[i] = true;
    }

    std::vector<MergeCandidate> heap_storage;
    heap_storage.reserve(storage.NumPairs() + (storage.NumPairs() / 2));
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i + 1; j < n; ++j) {
            heap_storage.push_back(
                make_candidate(cluster_distance(distances, max_nodes, i, j), i, j));
        }
    }
    std::priority_queue<
        MergeCandidate,
        std::vector<MergeCandidate>,
        MergeCandidateGreater>
        heap(MergeCandidateGreater{}, std::move(heap_storage));

    std::vector<size_t> children_left;
    std::vector<size_t> children_right;
    std::vector<double> distances_out;
    std::vector<size_t> merge_cluster_sizes;
    children_left.reserve(n - 1);
    children_right.reserve(n - 1);
    distances_out.reserve(n - 1);
    merge_cluster_sizes.reserve(n - 1);

    size_t next_node = n;
    size_t active_count = n;
    while (active_count > 1 && distances_out.size() < target_merges) {
        MergeCandidate best;
        bool found = false;
        while (!heap.empty()) {
            best = heap.top();
            heap.pop();
            if (active[best.left] && active[best.right]) {
                found = true;
                break;
            }
        }
        if (!found) {
            throw std::runtime_error("Agglomerative clustering heap was exhausted");
        }

        const size_t left = best.left;
        const size_t right = best.right;
        const size_t merged_node = next_node++;
        const size_t merged_size = cluster_sizes[left] + cluster_sizes[right];

        children_left.push_back(left);
        children_right.push_back(right);
        distances_out.push_back(best.distance);
        merge_cluster_sizes.push_back(merged_size);

        active[left] = false;
        active[right] = false;
        active[merged_node] = true;
        cluster_sizes[merged_node] = merged_size;
        --active_count;

        for (size_t node = 0; node < merged_node; ++node) {
            if (!active[node]) {
                continue;
            }

            const double updated_distance = update_linkage_distance(
                options.linkage,
                cluster_distance(distances, max_nodes, left, node),
                cluster_distance(distances, max_nodes, right, node),
                cluster_sizes[left],
                cluster_sizes[right]);
            cluster_distance(distances, max_nodes, merged_node, node) = updated_distance;
            heap.push(make_candidate(updated_distance, merged_node, node));
        }
    }

    std::vector<ClusterLabel> labels =
        labels_from_cut(children_left, children_right, distances_out, n, options);
    Clusters members = labels_to_clusters(labels);
    return AgglomerativeResult(std::move(labels), std::move(members),
                               std::move(children_left), std::move(children_right),
                               std::move(distances_out), std::move(merge_cluster_sizes));
}

}  // namespace agglomerative_oracle
