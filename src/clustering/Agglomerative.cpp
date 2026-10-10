/**
 * @file Agglomerative.cpp
 * @brief Hierarchical agglomerative clustering implementation.
 */

#include "oecluster/clustering/Agglomerative.h"

#include <algorithm>
#include <cmath>
#include <optional>
#include <stdexcept>
#include <vector>

#include "AgglomerativeInternal.h"
#include "AgglomerativeRowCache.h"
#include "DistanceAccess.h"
#include "DiversityValidation.h"
#include "HDBSCANLinkage.h"
#include "PrimMST.h"

namespace OECluster {

namespace {

constexpr const char* AGGLOMERATIVE_NAME = "agglomerative";

// The checks that need no item count, in the order 5.19.0 applied them.
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

// Single linkage's merges are the spanning tree's edges in ascending order.
// Within a tied height they come in the tree's order, which can differ from
// 5.20.0's heap; the heights, and so every distance_threshold cut, cannot.
AgglomerativeResult single_linkage_result(std::vector<detail::HDBSCANMSTEdge> mst,
                                          size_t n,
                                          const AgglomerativeOptions& options) {
    const std::vector<detail::HDBSCANLinkageNode> linkage =
        detail::make_hdbscan_single_linkage(std::move(mst), n);
    const bool cut_by_threshold = options.distance_threshold >= 0.0;
    const size_t merges =
        (!options.compute_full_tree && !cut_by_threshold) ? n - options.n_clusters
                                                          : n - 1;

    std::vector<size_t> children_left;
    std::vector<size_t> children_right;
    std::vector<double> distances;
    std::vector<size_t> cluster_sizes;
    children_left.reserve(merges);
    children_right.reserve(merges);
    distances.reserve(merges);
    cluster_sizes.reserve(merges);
    for (size_t merge = 0; merge < merges; ++merge) {
        const detail::HDBSCANLinkageNode& node = linkage[merge];
        children_left.push_back(std::min(node.left_node, node.right_node));
        children_right.push_back(std::max(node.left_node, node.right_node));
        distances.push_back(node.value);
        cluster_sizes.push_back(node.cluster_size);
    }

    std::vector<ClusterLabel> labels =
        detail::labels_from_cut(children_left, children_right, distances, n, options);
    Clusters members = labels_to_clusters(labels);
    return AgglomerativeResult(std::move(labels), std::move(members),
                               std::move(children_left), std::move(children_right),
                               std::move(distances), std::move(cluster_sizes));
}

detail::PrimOptions prim_options(const AgglomerativeOptions& options) {
    detail::PrimOptions prim;
    prim.num_threads = options.num_threads;
    prim.caller = AGGLOMERATIVE_NAME;
    return prim;
}

}  // namespace

AgglomerativeResult agglomerative_cluster(
    const StorageBackend& storage,
    const AgglomerativeOptions& options) {
    if (options.linkage != AgglomerativeLinkageMethod::Single) {
        return detail::row_cache_result(storage, options);
    }
    validate_options(storage, options);
    const size_t n = storage.NumSamples();
    if (n < 2) {
        return small_result(n);
    }
    return single_linkage_result(
        detail::prim_mst(storage, detail::PrimWeights(), prim_options(options)), n,
        options);
}

AgglomerativeResult agglomerative_cluster(
    PairwiseComparison& comparison,
    const AgglomerativeOptions& options) {
    validate_arguments(options);
    if (options.linkage != AgglomerativeLinkageMethod::Single) {
        throw std::invalid_argument(
            "Agglomerative clustering from a comparison supports only single "
            "linkage; complete, average and weighted linkage need a matrix from "
            "pdist()");
    }
    detail::refuse_unrepeatable(comparison, AGGLOMERATIVE_NAME, "cluster",
                                "the tree would depend on the order the pairs "
                                "are scored");
    detail::validate_distance_facts(comparison, AGGLOMERATIVE_NAME);
    const size_t n = comparison.Size();
    validate_cluster_bound(options, n);
    if (n < 2) {
        return small_result(n);
    }
    return single_linkage_result(
        detail::prim_mst(comparison, detail::PrimWeights(), prim_options(options)), n,
        options);
}

namespace detail {

AgglomerativeResult row_cache_result(
    const StorageBackend& storage,
    const AgglomerativeOptions& options,
    std::optional<size_t> serial_cutoff,
    RowCacheStats* stats) {
    validate_options(storage, options);
    const size_t n = storage.NumSamples();
    if (n < 2) {
        return small_result(n);
    }

    const bool cut_by_threshold = options.distance_threshold >= 0.0;
    RowCacheOptions kernel;
    kernel.linkage = options.linkage;
    kernel.target_merges = (!options.compute_full_tree && !cut_by_threshold)
                               ? n - options.n_clusters
                               : n - 1;
    kernel.num_threads = options.num_threads;
    kernel.chunk_size = options.chunk_size;
    kernel.serial_cutoff = serial_cutoff;
    kernel.stats = stats;

    LinkageTree tree = agglomerative_row_cache(storage, kernel);
    std::vector<ClusterLabel> labels =
        labels_from_cut(tree.children_left, tree.children_right, tree.distances, n,
                        options);
    Clusters members = labels_to_clusters(labels);
    return AgglomerativeResult(std::move(labels), std::move(members),
                               std::move(tree.children_left),
                               std::move(tree.children_right),
                               std::move(tree.distances),
                               std::move(tree.cluster_sizes));
}

}  // namespace detail

}  // namespace OECluster
