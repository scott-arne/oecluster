/**
 * @file Agglomerative.cpp
 * @brief Hierarchical agglomerative clustering implementation.
 */

#include "oecluster/clustering/Agglomerative.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <queue>
#include <stdexcept>
#include <vector>

#include "DistanceAccess.h"
#include "oecluster/ThreadPool.h"

namespace OECluster {

namespace {

class LeafUnionFind {
public:
    explicit LeafUnionFind(size_t n)
        : parent_(n),
          rank_(n, 0) {
        for (size_t i = 0; i < n; ++i) {
            parent_[i] = i;
        }
    }

    size_t Find(size_t node) {
        size_t root = node;
        while (parent_[root] != root) {
            root = parent_[root];
        }
        while (parent_[node] != root) {
            const size_t next = parent_[node];
            parent_[node] = root;
            node = next;
        }
        return root;
    }

    size_t Union(size_t left, size_t right) {
        size_t left_root = Find(left);
        size_t right_root = Find(right);
        if (left_root == right_root) {
            return left_root;
        }
        if (rank_[left_root] < rank_[right_root]) {
            std::swap(left_root, right_root);
        }
        parent_[right_root] = left_root;
        if (rank_[left_root] == rank_[right_root]) {
            ++rank_[left_root];
        }
        return left_root;
    }

private:
    std::vector<size_t> parent_;
    std::vector<size_t> rank_;
};

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

size_t max_node_count(size_t n_samples) {
    return n_samples == 0 ? 0 : (2 * n_samples - 1);
}

double& cluster_distance(
    std::vector<double>& distances,
    size_t max_nodes,
    size_t left,
    size_t right) {
    return distances[detail::condensed_index(max_nodes, left, right)];
}

double cluster_distance(
    const std::vector<double>& distances,
    size_t max_nodes,
    size_t left,
    size_t right) {
    return distances[detail::condensed_index(max_nodes, left, right)];
}

MergeCandidate make_candidate(double distance, size_t left, size_t right) {
    if (right < left) {
        std::swap(left, right);
    }
    return MergeCandidate{distance, left, right};
}

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

void validate_options(
    const StorageBackend& storage,
    const AgglomerativeOptions& options) {
    detail::validate_complete_distance_storage(storage, "Agglomerative clustering");

    if (options.chunk_size == 0) {
        throw std::invalid_argument("Agglomerative chunk_size must be at least one");
    }
    if (std::isnan(options.distance_threshold)) {
        throw std::invalid_argument("Agglomerative distance_threshold must not be NaN");
    }
    if (options.distance_threshold < 0.0) {
        if (options.n_clusters == 0) {
            throw std::invalid_argument("Agglomerative n_clusters must be at least one");
        }
        if (options.n_clusters > storage.NumItems()) {
            throw std::invalid_argument(
                "Agglomerative n_clusters must be at most the item count");
        }
    }
}

std::vector<double> initialize_cluster_distances(
    const StorageBackend& storage,
    size_t max_nodes,
    size_t num_threads,
    size_t chunk_size) {
    const size_t n = storage.NumItems();
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

std::vector<ClusterLabel> labels_from_cut(
    const std::vector<size_t>& children_left,
    const std::vector<size_t>& children_right,
    const std::vector<double>& distances,
    size_t n_samples,
    const AgglomerativeOptions& options) {
    if (n_samples == 0) {
        return {};
    }

    LeafUnionFind union_find(n_samples);
    std::vector<size_t> representatives(max_node_count(n_samples), 0);
    for (size_t i = 0; i < n_samples; ++i) {
        representatives[i] = i;
    }

    const bool cut_by_threshold = options.distance_threshold >= 0.0;
    const size_t merges_to_apply =
        cut_by_threshold ? distances.size() : n_samples - options.n_clusters;

    for (size_t merge = 0; merge < distances.size(); ++merge) {
        if (merge >= merges_to_apply) {
            break;
        }
        if (cut_by_threshold && distances[merge] > options.distance_threshold) {
            break;
        }

        const size_t left = children_left[merge];
        const size_t right = children_right[merge];
        const size_t root = union_find.Union(
            representatives[left],
            representatives[right]);
        representatives[n_samples + merge] = union_find.Find(root);
    }

    std::vector<ClusterLabel> labels(n_samples, NOISE_LABEL);
    std::vector<size_t> roots;
    roots.reserve(n_samples);
    for (size_t i = 0; i < n_samples; ++i) {
        const size_t root = union_find.Find(i);
        auto iter = std::find(roots.begin(), roots.end(), root);
        if (iter == roots.end()) {
            roots.push_back(root);
            labels[i] = static_cast<ClusterLabel>(roots.size() - 1);
        } else {
            labels[i] = static_cast<ClusterLabel>(
                static_cast<size_t>(std::distance(roots.begin(), iter)));
        }
    }

    return labels;
}

}  // namespace

AgglomerativeResult agglomerative_cluster(
    const StorageBackend& storage,
    const AgglomerativeOptions& options) {
    validate_options(storage, options);

    const size_t n = storage.NumItems();
    if (n == 0) {
        return AgglomerativeResult();
    }
    if (n == 1) {
        std::vector<ClusterLabel> labels{0};
        Clusters members = labels_to_clusters(labels);
        return AgglomerativeResult(std::move(labels), std::move(members),
                                   {}, {}, {}, {});
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

}  // namespace OECluster
