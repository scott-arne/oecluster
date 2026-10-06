/**
 * @file HDBSCANLinkage.cpp
 * @brief Internal MST and single-linkage utilities for HDBSCAN.
 */

#include "HDBSCANLinkage.h"

#include <algorithm>
#include <limits>
#include <stdexcept>


namespace OECluster::detail {

namespace {

class HDBSCANLinkageUnionFind {
public:
    explicit HDBSCANLinkageUnionFind(size_t n_samples)
        : parent_(2 * n_samples - 1, kInvalid),
          size_(2 * n_samples - 1, 0),
          next_label_(n_samples) {
        for (size_t i = 0; i < n_samples; ++i) {
            size_[i] = 1;
        }
    }

    // Two-pass path compression flattens union-find trees to maintain O(α(n)) amortized time.
    size_t Find(size_t node) {
        size_t root = node;
        while (parent_[root] != kInvalid) {
            root = parent_[root];
        }
        while (parent_[node] != kInvalid && parent_[node] != root) {
            const size_t next = parent_[node];
            parent_[node] = root;
            node = next;
        }
        return root;
    }

    size_t Union(size_t left, size_t right) {
        const size_t label = next_label_++;
        parent_[left] = label;
        parent_[right] = label;
        size_[label] = size_[left] + size_[right];
        return label;
    }

    size_t Size(size_t node) const {
        return size_[node];
    }

private:
    static constexpr size_t kInvalid = std::numeric_limits<size_t>::max();

    std::vector<size_t> parent_;
    std::vector<size_t> size_;
    size_t next_label_;
};

}  // namespace

std::vector<HDBSCANLinkageNode> make_hdbscan_single_linkage(
    std::vector<HDBSCANMSTEdge> mst,
    size_t n_samples) {
    if (n_samples == 0) {
        return {};
    }
    if (mst.size() + 1 != n_samples) {
        throw std::invalid_argument("MST edge count must be n_samples - 1");
    }

    std::sort(mst.begin(), mst.end(), [](const HDBSCANMSTEdge& lhs, const HDBSCANMSTEdge& rhs) {
        if (lhs.distance != rhs.distance) {
            return lhs.distance < rhs.distance;
        }
        if (lhs.current_node != rhs.current_node) {
            return lhs.current_node < rhs.current_node;
        }
        return lhs.next_node < rhs.next_node;
    });

    HDBSCANLinkageUnionFind union_find(n_samples);
    std::vector<HDBSCANLinkageNode> linkage;
    linkage.reserve(mst.size());

    for (const HDBSCANMSTEdge& edge : mst) {
        const size_t left = union_find.Find(edge.current_node);
        const size_t right = union_find.Find(edge.next_node);
        const size_t cluster_size = union_find.Size(left) + union_find.Size(right);
        linkage.push_back(HDBSCANLinkageNode{left, right, edge.distance, cluster_size});
        union_find.Union(left, right);
    }

    return linkage;
}

}  // namespace OECluster::detail
