/**
 * @file HDBSCAN.cpp
 * @brief HDBSCAN clustering implementation and dense infrastructure.
 */

#include "oecluster/clustering/HDBSCAN.h"

#include <algorithm>
#include <stdexcept>
#include <vector>

#include "DistanceAccess.h"
#include "HDBSCANLinkage.h"
#include "HDBSCANTree.h"
#include "oecluster/ThreadPool.h"

namespace OECluster {

namespace detail {

std::vector<double> compute_core_distances(
    const StorageBackend& storage,
    size_t min_samples,
    size_t num_threads) {
    validate_complete_distance_storage(storage, "HDBSCAN");
    const size_t n = storage.NumSamples();
    if (min_samples == 0) {
        throw std::invalid_argument("HDBSCAN min_samples must be at least one");
    }
    if (min_samples > n) {
        throw std::invalid_argument("HDBSCAN min_samples must be at most the item count");
    }

    std::vector<double> core_distances(n, 0.0);
    if (min_samples == 1 || n == 0) {
        return core_distances;
    }

    const double* data = storage.Data();
    ThreadPool pool(num_threads);
    pool.ParallelFor(0, n, 64, [&](size_t begin, size_t end) {
        std::vector<double> row(n);
        for (size_t i = begin; i < end; ++i) {
            row[0] = 0.0;
            size_t write_index = 1;
            for (size_t j = 0; j < n; ++j) {
                if (i != j) {
                    row[write_index++] = dense_distance(data, n, i, j);
                }
            }
            const auto nth = row.begin() + static_cast<std::ptrdiff_t>(min_samples - 1);
            std::nth_element(row.begin(), nth, row.end());
            core_distances[i] = *nth;
        }
    });

    return core_distances;
}

}  // namespace detail

HDBSCANResult hdbscan_cluster(const StorageBackend& storage, const HDBSCANOptions& options) {
    if (options.min_cluster_size < 2) {
        throw std::invalid_argument("HDBSCAN min_cluster_size must be at least two");
    }
    if (options.alpha <= 0.0) {
        throw std::invalid_argument("HDBSCAN alpha must be positive");
    }

    const size_t min_samples =
        options.min_samples == 0 ? options.min_cluster_size : options.min_samples;
    if (min_samples > storage.NumSamples()) {
        throw std::invalid_argument("HDBSCAN min_samples must be at most the item count");
    }

    detail::validate_complete_distance_storage(storage, "HDBSCAN");

    const std::vector<double> core_distances =
        detail::compute_core_distances(storage, min_samples, options.num_threads);
    const std::vector<detail::HDBSCANMSTEdge> mst =
        detail::hdbscan_mutual_reachability_mst(storage, core_distances, options.alpha);
    const std::vector<detail::HDBSCANLinkageNode> linkage =
        detail::make_hdbscan_single_linkage(mst, storage.NumSamples());
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
        labels.assign(storage.NumSamples(), NOISE_LABEL);
    }
    if (probabilities.empty()) {
        probabilities.assign(storage.NumSamples(), 0.0);
    }
    Clusters members = labels_to_clusters(labels);
    return HDBSCANResult(std::move(labels), std::move(members),
                         std::move(probabilities));
}

}  // namespace OECluster
