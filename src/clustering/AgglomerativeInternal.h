/**
 * @file AgglomerativeInternal.h
 * @brief Dendrogram helpers shared by the linkage kernel and its test oracle.
 *
 * `labels_from_cut` cuts a merge tree into labels. Both the production path and
 * the frozen 5.20.0 oracle must cut with the same code: two copies would let a
 * divergence between them read as a kernel defect.
 */

#ifndef OECLUSTER_CLUSTERING_AGGLOMERATIVE_INTERNAL_H
#define OECLUSTER_CLUSTERING_AGGLOMERATIVE_INTERNAL_H

#include <cstddef>
#include <vector>

#include "oecluster/clustering/Agglomerative.h"

namespace OECluster::detail {

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

/// The dendrogram's node count: n leaves plus n-1 internal nodes.
inline size_t max_node_count(size_t n_samples) {
    return n_samples == 0 ? 0 : (2 * n_samples - 1);
}

inline std::vector<ClusterLabel> labels_from_cut(
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

// Single linkage's merges are the spanning tree's edges in ascending order.
// Within a tied height they come in the tree's order, which can differ from

}  // namespace OECluster::detail

#endif  // OECLUSTER_CLUSTERING_AGGLOMERATIVE_INTERNAL_H
