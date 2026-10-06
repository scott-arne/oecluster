/**
 * @file HDBSCAN.cpp
 * @brief HDBSCAN clustering implementation and dense infrastructure.
 */

#include "oecluster/clustering/HDBSCAN.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>
#include <vector>

#include "DistanceAccess.h"
#include "DiversityValidation.h"
#include "HDBSCANLinkage.h"
#include "HDBSCANTree.h"
#include "PrimMST.h"
#include "StreamingCoreDistances.h"

namespace OECluster {

namespace {

constexpr const char* HDBSCAN_NAME = "hdbscan";

// The checks that need no item count, in the order 5.19.0 applied them.
size_t validated_min_samples(const HDBSCANOptions& options) {
    if (options.min_cluster_size < 2) {
        throw std::invalid_argument("HDBSCAN min_cluster_size must be at least two");
    }
    // Written as !(alpha > 0) so that NaN, which alpha <= 0 let through, is refused.
    if (!(options.alpha > 0.0)) {
        throw std::invalid_argument("HDBSCAN alpha must be positive");
    }
    return options.min_samples == 0 ? options.min_cluster_size : options.min_samples;
}

void validate_min_samples_bound(size_t min_samples, size_t n) {
    if (min_samples > n) {
        throw std::invalid_argument("HDBSCAN min_samples must be at most the item count");
    }
}

// With a core pass every pair was read, so the largest distance decides whether
// any quotient overflows; with min_samples = 1 the unpruned Prim pass checks
// each quotient as it reads it.
detail::PrimWeights mutual_reachability_weights(detail::CoreDistances core,
                                                size_t min_samples, double alpha) {
    if (min_samples > 1 && !std::isfinite(core.max_distance / alpha)) {
        throw detail::alpha_overflow_error(HDBSCAN_NAME, alpha);
    }
    detail::PrimWeights weights;
    weights.core = std::move(core.values);
    weights.alpha = alpha;
    weights.prune = min_samples > 1;
    return weights;
}

detail::PrimOptions prim_options(const HDBSCANOptions& options) {
    detail::PrimOptions prim;
    prim.num_threads = options.num_threads;
    prim.caller = HDBSCAN_NAME;
    return prim;
}

HDBSCANResult hdbscan_from_tree(std::vector<detail::HDBSCANMSTEdge> mst, size_t n,
                                const HDBSCANOptions& options) {
    const std::vector<detail::HDBSCANLinkageNode> linkage =
        detail::make_hdbscan_single_linkage(std::move(mst), n);
    const std::vector<detail::CondensedNode> condensed_tree =
        detail::condense_tree(linkage, options.min_cluster_size);
    const detail::HDBSCANTreeSelection selection =
        detail::select_clusters(
            condensed_tree,
            options.cluster_selection_method,
            options.allow_single_cluster,
            options.cluster_selection_epsilon,
            options.max_cluster_size);

    std::vector<ClusterLabel> labels = selection.labels;
    std::vector<double> probabilities = selection.probabilities;
    if (labels.empty()) {
        labels.assign(n, NOISE_LABEL);
    }
    if (probabilities.empty()) {
        probabilities.assign(n, 0.0);
    }
    Clusters members = labels_to_clusters(labels);
    return HDBSCANResult(std::move(labels), std::move(members),
                         std::move(probabilities));
}

}  // namespace

namespace detail {

std::vector<double> compute_core_distances(
    const StorageBackend& storage,
    size_t min_samples,
    size_t num_threads) {
    return matrix_core_distances(storage, min_samples, num_threads, 0, HDBSCAN_NAME)
        .values;
}

}  // namespace detail

HDBSCANResult hdbscan_cluster(const StorageBackend& storage, const HDBSCANOptions& options) {
    const size_t min_samples = validated_min_samples(options);
    validate_min_samples_bound(min_samples, storage.NumSamples());
    detail::validate_complete_distance_storage(storage, "HDBSCAN");

    detail::PrimWeights weights = mutual_reachability_weights(
        detail::matrix_core_distances(storage, min_samples, options.num_threads,
                                      options.chunk_size, HDBSCAN_NAME),
        min_samples, options.alpha);
    return hdbscan_from_tree(detail::prim_mst(storage, weights, prim_options(options)),
                             storage.NumSamples(), options);
}

HDBSCANResult hdbscan_cluster(PairwiseComparison& comparison, const HDBSCANOptions& options) {
    const size_t min_samples = validated_min_samples(options);
    detail::refuse_unrepeatable(comparison, HDBSCAN_NAME, "cluster",
                                "its core-distance and spanning-tree passes can "
                                "disagree");
    detail::validate_distance_facts(comparison, HDBSCAN_NAME);
    const size_t n = comparison.Size();
    validate_min_samples_bound(min_samples, n);

    detail::PrimWeights weights = mutual_reachability_weights(
        detail::streaming_core_distances(comparison, min_samples, options.num_threads,
                                         HDBSCAN_NAME),
        min_samples, options.alpha);
    return hdbscan_from_tree(
        detail::prim_mst(comparison, weights, prim_options(options)), n, options);
}

}  // namespace OECluster
