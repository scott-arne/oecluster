/**
 * @file HDBSCANLinkage.cpp
 * @brief Internal MST and single-linkage utilities for HDBSCAN.
 */

#include "HDBSCANLinkage.h"

#include <algorithm>
#include <limits>
#include <stdexcept>

#include "DistanceAccess.h"

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

double mutual_reachability(
    const double* data,
    size_t n,
    const std::vector<double>& core_distances,
    size_t i,
    size_t j,
    double alpha) {
    return std::max({core_distances[i], core_distances[j],
                     dense_distance(data, n, i, j) / alpha});
}

}  // namespace

std::vector<HDBSCANMSTEdge> hdbscan_mutual_reachability_mst(
    const StorageBackend& storage,
    const std::vector<double>& core_distances,
    double alpha) {
    validate_complete_distance_storage(storage, "HDBSCAN");
    if (alpha <= 0.0) {
        throw std::invalid_argument("HDBSCAN alpha must be positive");
    }

    const size_t n = storage.NumItems();
    if (core_distances.size() != n) {
        throw std::invalid_argument("Core distance count must match storage item count");
    }
    if (n < 2) {
        return {};
    }

    std::vector<HDBSCANMSTEdge> mst;
    mst.reserve(n - 1);
    std::vector<bool> in_tree(n, false);
    std::vector<double> min_reachability(n, std::numeric_limits<double>::infinity());
    std::vector<size_t> current_sources(n, 0);
    const double* data = storage.Data();
    size_t current_node = 0;

    for (size_t step = 0; step < n - 1; ++step) {
        in_tree[current_node] = true;

        double new_reachability = std::numeric_limits<double>::infinity();
        size_t source_node = 0;
        size_t new_node = 0;

        for (size_t candidate = 0; candidate < n; ++candidate) {
            if (in_tree[candidate]) {
                continue;
            }

            const double next_min_reach = min_reachability[candidate];
            const size_t next_source = current_sources[candidate];
            const double distance = mutual_reachability(
                data, n, core_distances, current_node, candidate, alpha);

            if (distance < next_min_reach) {
                min_reachability[candidate] = distance;
                current_sources[candidate] = current_node;
                if (distance < new_reachability) {
                    new_reachability = distance;
                    source_node = current_node;
                    new_node = candidate;
                }
            } else if (next_min_reach < new_reachability) {
                new_reachability = next_min_reach;
                source_node = next_source;
                new_node = candidate;
            }
        }

        mst.push_back(HDBSCANMSTEdge{source_node, new_node, new_reachability});
        current_node = new_node;
    }

    return mst;
}

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
