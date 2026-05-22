/**
 * @file HDBSCANTree.cpp
 * @brief Internal condensed-tree utilities for HDBSCAN.
 */

#include "HDBSCANTree.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace OECluster::detail {

namespace {

constexpr double INFTY = std::numeric_limits<double>::infinity();

std::vector<size_t> bfs_from_hierarchy(
    const std::vector<HDBSCANLinkageNode>& hierarchy,
    size_t root) {
    const size_t n_samples = hierarchy.size() + 1;
    std::vector<size_t> result;
    std::vector<size_t> process_queue{root};

    while (!process_queue.empty()) {
        result.insert(result.end(), process_queue.begin(), process_queue.end());

        std::vector<size_t> next_queue;
        for (const size_t node : process_queue) {
            if (node >= n_samples) {
                const HDBSCANLinkageNode& children = hierarchy.at(node - n_samples);
                next_queue.push_back(children.left_node);
                next_queue.push_back(children.right_node);
            }
        }
        process_queue = std::move(next_queue);
    }

    return result;
}

size_t child_count(
    const std::vector<HDBSCANLinkageNode>& hierarchy,
    size_t n_samples,
    size_t node) {
    return node >= n_samples ? hierarchy.at(node - n_samples).cluster_size : 1;
}

std::unordered_map<size_t, double> compute_stability(
    const std::vector<CondensedNode>& condensed_tree) {
    std::unordered_map<size_t, double> stability;
    if (condensed_tree.empty()) {
        return stability;
    }

    size_t smallest_parent = condensed_tree.front().parent;
    size_t largest_parent = condensed_tree.front().parent;
    size_t largest_child = condensed_tree.front().child;
    for (const CondensedNode& node : condensed_tree) {
        smallest_parent = std::min(smallest_parent, node.parent);
        largest_parent = std::max(largest_parent, node.parent);
        largest_child = std::max(largest_child, node.child);
    }

    std::vector<double> births(std::max(largest_child, smallest_parent) + 1, 0.0);
    for (const CondensedNode& node : condensed_tree) {
        births[node.child] = node.value;
    }
    births[smallest_parent] = 0.0;

    for (size_t parent = smallest_parent; parent <= largest_parent; ++parent) {
        stability[parent] = 0.0;
    }
    for (const CondensedNode& node : condensed_tree) {
        stability[node.parent] +=
            (node.value - births[node.parent]) * static_cast<double>(node.cluster_size);
    }

    return stability;
}

std::vector<CondensedNode> cluster_tree_nodes(
    const std::vector<CondensedNode>& condensed_tree) {
    std::vector<CondensedNode> result;
    for (const CondensedNode& node : condensed_tree) {
        if (node.cluster_size > 1) {
            result.push_back(node);
        }
    }
    return result;
}

std::vector<size_t> bfs_from_cluster_tree(
    const std::vector<CondensedNode>& cluster_tree,
    size_t root) {
    std::vector<size_t> result;
    std::vector<size_t> process_queue{root};

    while (!process_queue.empty()) {
        result.insert(result.end(), process_queue.begin(), process_queue.end());
        std::vector<size_t> next_queue;
        for (const size_t parent : process_queue) {
            for (const CondensedNode& node : cluster_tree) {
                if (node.parent == parent) {
                    next_queue.push_back(node.child);
                }
            }
        }
        process_queue = std::move(next_queue);
    }

    return result;
}

std::vector<double> max_lambdas(const std::vector<CondensedNode>& condensed_tree) {
    size_t largest_parent = 0;
    for (const CondensedNode& node : condensed_tree) {
        largest_parent = std::max(largest_parent, node.parent);
    }

    std::vector<double> deaths(largest_parent + 1, 0.0);
    for (const CondensedNode& node : condensed_tree) {
        deaths[node.parent] = std::max(deaths[node.parent], node.value);
    }
    return deaths;
}

class TreeUnionFind {
public:
    explicit TreeUnionFind(size_t size)
        : parent_(size), rank_(size, 0) {
        std::iota(parent_.begin(), parent_.end(), 0);
    }

    size_t Find(size_t node) {
        if (parent_[node] != node) {
            parent_[node] = Find(parent_[node]);
        }
        return parent_[node];
    }

    void Union(size_t left, size_t right) {
        size_t left_root = Find(left);
        size_t right_root = Find(right);
        if (left_root == right_root) {
            return;
        }
        if (rank_[left_root] < rank_[right_root]) {
            parent_[left_root] = right_root;
        } else if (rank_[left_root] > rank_[right_root]) {
            parent_[right_root] = left_root;
        } else {
            parent_[right_root] = left_root;
            ++rank_[left_root];
        }
    }

private:
    std::vector<size_t> parent_;
    std::vector<size_t> rank_;
};

std::vector<size_t> cluster_tree_leaves(
    const std::vector<CondensedNode>& cluster_tree) {
    if (cluster_tree.empty()) {
        return {};
    }

    size_t root = cluster_tree.front().parent;
    for (const CondensedNode& node : cluster_tree) {
        root = std::min(root, node.parent);
    }

    std::vector<size_t> leaves;
    std::vector<size_t> stack{root};
    while (!stack.empty()) {
        const size_t current = stack.back();
        stack.pop_back();

        bool has_children = false;
        for (const CondensedNode& node : cluster_tree) {
            if (node.parent == current) {
                has_children = true;
                stack.push_back(node.child);
            }
        }
        if (!has_children) {
            leaves.push_back(current);
        }
    }
    std::sort(leaves.begin(), leaves.end());
    return leaves;
}

double parent_distance_for_child(
    const std::vector<CondensedNode>& cluster_tree,
    size_t child) {
    for (const CondensedNode& node : cluster_tree) {
        if (node.child == child) {
            return node.value == 0.0 ? INFTY : 1.0 / node.value;
        }
    }
    return INFTY;
}

size_t parent_for_child(
    const std::vector<CondensedNode>& cluster_tree,
    size_t child) {
    for (const CondensedNode& node : cluster_tree) {
        if (node.child == child) {
            return node.parent;
        }
    }
    return child;
}

size_t traverse_upwards(
    const std::vector<CondensedNode>& cluster_tree,
    double cluster_selection_epsilon,
    size_t leaf,
    bool allow_single_cluster) {
    size_t root = cluster_tree.front().parent;
    for (const CondensedNode& node : cluster_tree) {
        root = std::min(root, node.parent);
    }

    const size_t parent = parent_for_child(cluster_tree, leaf);
    if (parent == root) {
        return allow_single_cluster ? parent : leaf;
    }

    const double parent_eps = parent_distance_for_child(cluster_tree, parent);
    if (parent_eps > cluster_selection_epsilon) {
        return parent;
    }
    return traverse_upwards(
        cluster_tree, cluster_selection_epsilon, parent, allow_single_cluster);
}

std::unordered_set<size_t> epsilon_search(
    const std::vector<size_t>& leaves,
    const std::vector<CondensedNode>& cluster_tree,
    double cluster_selection_epsilon,
    bool allow_single_cluster) {
    std::unordered_set<size_t> selected_clusters;
    std::unordered_set<size_t> processed;

    for (const size_t leaf : leaves) {
        const double eps = parent_distance_for_child(cluster_tree, leaf);
        if (eps < cluster_selection_epsilon) {
            if (processed.find(leaf) == processed.end()) {
                const size_t epsilon_child = traverse_upwards(
                    cluster_tree, cluster_selection_epsilon, leaf, allow_single_cluster);
                selected_clusters.insert(epsilon_child);
                for (const size_t sub_node : bfs_from_cluster_tree(cluster_tree, epsilon_child)) {
                    if (sub_node != epsilon_child) {
                        processed.insert(sub_node);
                    }
                }
            }
        } else {
            selected_clusters.insert(leaf);
        }
    }

    return selected_clusters;
}

size_t infer_sample_count(const std::vector<CondensedNode>& condensed_tree) {
    size_t root = condensed_tree.front().parent;
    size_t max_point = 0;
    bool has_point = false;
    for (const CondensedNode& node : condensed_tree) {
        root = std::min(root, node.parent);
    }
    for (const CondensedNode& node : condensed_tree) {
        if (node.cluster_size == 1 && node.child < root) {
            max_point = std::max(max_point, node.child);
            has_point = true;
        }
    }
    return has_point ? max_point + 1 : root;
}

std::vector<ClusterLabel> do_labelling(
    const std::vector<CondensedNode>& condensed_tree,
    const std::unordered_set<size_t>& clusters,
    const std::unordered_map<size_t, ClusterLabel>& cluster_label_map,
    bool allow_single_cluster,
    double cluster_selection_epsilon) {
    size_t root_cluster = condensed_tree.front().parent;
    size_t max_node = 0;
    for (const CondensedNode& node : condensed_tree) {
        root_cluster = std::min(root_cluster, node.parent);
        max_node = std::max({max_node, node.parent, node.child});
    }

    std::vector<ClusterLabel> result(root_cluster, NOISE_LABEL);
    TreeUnionFind union_find(max_node + 1);
    for (const CondensedNode& node : condensed_tree) {
        if (clusters.find(node.child) == clusters.end()) {
            union_find.Union(node.parent, node.child);
        }
    }

    double root_threshold = 0.0;
    for (const CondensedNode& node : condensed_tree) {
        if (node.parent == root_cluster) {
            root_threshold = std::max(root_threshold, node.value);
        }
    }

    for (size_t point = 0; point < root_cluster; ++point) {
        const size_t cluster = union_find.Find(point);
        ClusterLabel label = NOISE_LABEL;
        if (cluster != root_cluster) {
            label = cluster_label_map.at(cluster);
        } else if (clusters.size() == 1 && allow_single_cluster) {
            double parent_lambda = 0.0;
            for (const CondensedNode& node : condensed_tree) {
                if (node.child == point) {
                    parent_lambda = node.value;
                    break;
                }
            }
            const double threshold = cluster_selection_epsilon != 0.0
                                         ? 1.0 / cluster_selection_epsilon
                                         : root_threshold;
            if (parent_lambda >= threshold) {
                label = cluster_label_map.at(cluster);
            }
        }
        result[point] = label;
    }

    return result;
}

std::vector<double> get_probabilities(
    const std::vector<CondensedNode>& condensed_tree,
    const std::unordered_map<ClusterLabel, size_t>& reverse_cluster_map,
    const std::vector<ClusterLabel>& labels) {
    std::vector<double> result(labels.size(), 0.0);
    const std::vector<double> deaths = max_lambdas(condensed_tree);

    size_t root_cluster = condensed_tree.front().parent;
    for (const CondensedNode& node : condensed_tree) {
        root_cluster = std::min(root_cluster, node.parent);
    }

    for (const CondensedNode& node : condensed_tree) {
        const size_t point = node.child;
        if (point >= root_cluster) {
            continue;
        }
        const ClusterLabel label = labels[point];
        if (label == NOISE_LABEL) {
            continue;
        }
        const size_t cluster = reverse_cluster_map.at(label);
        const double max_lambda = deaths[cluster];
        if (max_lambda == 0.0 || std::isinf(node.value)) {
            result[point] = 1.0;
        } else {
            result[point] = std::min(node.value, max_lambda) / max_lambda;
        }
    }

    return result;
}

}  // namespace

std::vector<CondensedNode> condense_tree(
    const std::vector<HDBSCANLinkageNode>& hierarchy,
    size_t min_cluster_size) {
    if (hierarchy.empty()) {
        return {};
    }

    const size_t root = 2 * hierarchy.size();
    const size_t n_samples = hierarchy.size() + 1;
    size_t next_label = n_samples + 1;
    const std::vector<size_t> node_list = bfs_from_hierarchy(hierarchy, root);

    std::vector<size_t> relabel(root + 1, 0);
    relabel[root] = n_samples;
    std::vector<bool> ignore(root + 1, false);
    std::vector<CondensedNode> condensed;
    condensed.reserve(n_samples);

    for (const size_t node : node_list) {
        if (ignore[node] || node < n_samples) {
            continue;
        }

        const HDBSCANLinkageNode& children = hierarchy.at(node - n_samples);
        const size_t left = children.left_node;
        const size_t right = children.right_node;
        const double lambda_value = children.value > 0.0 ? 1.0 / children.value : INFTY;
        const size_t left_count = child_count(hierarchy, n_samples, left);
        const size_t right_count = child_count(hierarchy, n_samples, right);

        if (left_count >= min_cluster_size && right_count >= min_cluster_size) {
            relabel[left] = next_label++;
            condensed.push_back(
                CondensedNode{relabel[node], relabel[left], lambda_value, left_count});

            relabel[right] = next_label++;
            condensed.push_back(
                CondensedNode{relabel[node], relabel[right], lambda_value, right_count});
        } else if (left_count < min_cluster_size && right_count < min_cluster_size) {
            for (const size_t sub_node : bfs_from_hierarchy(hierarchy, left)) {
                if (sub_node < n_samples) {
                    condensed.push_back(CondensedNode{relabel[node], sub_node, lambda_value, 1});
                }
                ignore[sub_node] = true;
            }
            for (const size_t sub_node : bfs_from_hierarchy(hierarchy, right)) {
                if (sub_node < n_samples) {
                    condensed.push_back(CondensedNode{relabel[node], sub_node, lambda_value, 1});
                }
                ignore[sub_node] = true;
            }
        } else if (left_count < min_cluster_size) {
            relabel[right] = relabel[node];
            for (const size_t sub_node : bfs_from_hierarchy(hierarchy, left)) {
                if (sub_node < n_samples) {
                    condensed.push_back(CondensedNode{relabel[node], sub_node, lambda_value, 1});
                }
                ignore[sub_node] = true;
            }
        } else {
            relabel[left] = relabel[node];
            for (const size_t sub_node : bfs_from_hierarchy(hierarchy, right)) {
                if (sub_node < n_samples) {
                    condensed.push_back(CondensedNode{relabel[node], sub_node, lambda_value, 1});
                }
                ignore[sub_node] = true;
            }
        }
    }

    return condensed;
}

HDBSCANTreeSelection select_clusters(
    const std::vector<CondensedNode>& condensed_tree,
    HDBSCANClusterSelectionMethod cluster_selection_method,
    bool allow_single_cluster,
    double cluster_selection_epsilon,
    size_t max_cluster_size) {
    HDBSCANTreeSelection selection;
    if (condensed_tree.empty()) {
        return selection;
    }

    std::unordered_map<size_t, double> stability = compute_stability(condensed_tree);
    std::vector<size_t> node_list;
    node_list.reserve(stability.size());
    for (const auto& [node, _] : stability) {
        node_list.push_back(node);
    }
    std::sort(node_list.begin(), node_list.end(), std::greater<size_t>());
    if (!allow_single_cluster && !node_list.empty()) {
        node_list.pop_back();
    }

    const std::vector<CondensedNode> cluster_tree = cluster_tree_nodes(condensed_tree);
    std::unordered_map<size_t, bool> is_cluster;
    for (const size_t node : node_list) {
        is_cluster[node] = true;
    }

    const size_t n_samples = infer_sample_count(condensed_tree);
    if (max_cluster_size == 0) {
        max_cluster_size = n_samples + 1;
    }

    std::unordered_map<size_t, size_t> cluster_sizes;
    for (const CondensedNode& node : cluster_tree) {
        cluster_sizes[node.child] = node.cluster_size;
    }
    if (allow_single_cluster && !node_list.empty()) {
        const size_t root = node_list.back();
        size_t root_size = 0;
        for (const CondensedNode& node : cluster_tree) {
            if (node.parent == root) {
                root_size += node.cluster_size;
            }
        }
        cluster_sizes[root] = root_size;
    }

    if (cluster_selection_method == HDBSCANClusterSelectionMethod::EOM) {
        for (const size_t node : node_list) {
            double subtree_stability = 0.0;
            for (const CondensedNode& child : cluster_tree) {
                if (child.parent == node) {
                    subtree_stability += stability[child.child];
                }
            }

            const size_t cluster_size = cluster_sizes.count(node) ? cluster_sizes[node] : 0;
            if (subtree_stability > stability[node] || cluster_size > max_cluster_size) {
                is_cluster[node] = false;
                stability[node] = subtree_stability;
            } else {
                for (const size_t sub_node : bfs_from_cluster_tree(cluster_tree, node)) {
                    if (sub_node != node) {
                        is_cluster[sub_node] = false;
                    }
                }
            }
        }

        if (cluster_selection_epsilon != 0.0 && !cluster_tree.empty()) {
            std::vector<size_t> eom_clusters;
            for (const auto& [cluster, selected] : is_cluster) {
                if (selected) {
                    eom_clusters.push_back(cluster);
                }
            }
            std::unordered_set<size_t> selected_clusters;
            size_t root = cluster_tree.front().parent;
            for (const CondensedNode& node : cluster_tree) {
                root = std::min(root, node.parent);
            }
            if (eom_clusters.size() == 1 && eom_clusters.front() == root) {
                if (allow_single_cluster) {
                    selected_clusters.insert(root);
                }
            } else {
                selected_clusters = epsilon_search(
                    eom_clusters,
                    cluster_tree,
                    cluster_selection_epsilon,
                    allow_single_cluster);
            }
            for (auto& [cluster, selected] : is_cluster) {
                selected = selected_clusters.find(cluster) != selected_clusters.end();
            }
        }
    } else {
        std::vector<size_t> leaves = cluster_tree_leaves(cluster_tree);
        if (leaves.empty()) {
            for (auto& [_, selected] : is_cluster) {
                selected = false;
            }
            size_t root = condensed_tree.front().parent;
            for (const CondensedNode& node : condensed_tree) {
                root = std::min(root, node.parent);
            }
            is_cluster[root] = true;
        }

        std::unordered_set<size_t> selected_clusters;
        if (cluster_selection_epsilon != 0.0) {
            selected_clusters = epsilon_search(
                leaves,
                cluster_tree,
                cluster_selection_epsilon,
                allow_single_cluster);
        } else {
            selected_clusters.insert(leaves.begin(), leaves.end());
        }
        for (auto& [cluster, selected] : is_cluster) {
            selected = selected_clusters.find(cluster) != selected_clusters.end();
        }
    }

    std::vector<size_t> clusters;
    for (const auto& [cluster, selected] : is_cluster) {
        if (selected) {
            clusters.push_back(cluster);
        }
    }
    std::sort(clusters.begin(), clusters.end());

    std::unordered_map<size_t, ClusterLabel> cluster_label_map;
    std::unordered_map<ClusterLabel, size_t> reverse_cluster_map;
    ClusterLabel label = 0;
    for (const size_t cluster : clusters) {
        cluster_label_map[cluster] = label;
        reverse_cluster_map[label] = cluster;
        ++label;
    }

    std::unordered_set<size_t> cluster_set(clusters.begin(), clusters.end());
    selection.labels = do_labelling(
        condensed_tree,
        cluster_set,
        cluster_label_map,
        allow_single_cluster,
        cluster_selection_epsilon);
    selection.probabilities =
        get_probabilities(condensed_tree, reverse_cluster_map, selection.labels);
    return selection;
}

}  // namespace OECluster::detail
