/**
 * @file SphereExclusion.cpp
 * @brief Sphere exclusion (leader, Butina and DISE orders) over distance
 * matrices and comparisons.
 */

#include "oecluster/clustering/SphereExclusion.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "ChunkedComparisons.h"
#include "DistanceAccess.h"
#include "DiversityValidation.h"
#include "SphereExclusionEngine.h"
#include "ThresholdGraph.h"

namespace OECluster::detail {

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

// After forming a cluster, recompute unseen-neighbor counts for affected
// candidates to prioritize high-connectivity nodes in remaining items. Only
// the range from the cursor on is live: the consumed prefix is never read
// again, so updating and sorting just that range gives the order the old
// erase-then-sort-everything loop gave, without its quadratic cost.
void reorder_unseen_candidates(std::vector<Candidate>& candidates,
                               size_t cursor,
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

    const auto live = candidates.begin() + static_cast<std::ptrdiff_t>(cursor);
    for (auto it = live; it != candidates.end(); ++it) {
        const size_t idx = it->second;
        if (!affected[idx]) {
            continue;
        }
        size_t unseen_neighbors = 0;
        for (const size_t neighbor : graph.Neighbors(idx)) {
            if (!seen[neighbor]) {
                ++unseen_neighbors;
            }
        }
        it->first = unseen_neighbors;
    }
    std::sort(live, candidates.end(), std::greater<Candidate>());
}

}  // namespace

SphereEngineResult sphere_neighbors_first(const ThresholdNeighborGraph& graph,
                                          bool reordering) {
    std::vector<Candidate> candidates = make_sorted_candidates(graph);
    std::vector<bool> seen(graph.Size(), false);
    Clusters clusters;
    clusters.reserve(graph.Size());

    size_t cursor = 0;
    while (cursor < candidates.size() && candidates[cursor].first > 1) {
        const size_t idx = candidates[cursor].second;
        ++cursor;
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

        if (reordering) {
            reorder_unseen_candidates(candidates, cursor, graph, seen, cluster);
        }
    }

    for (; cursor < candidates.size(); ++cursor) {
        const size_t idx = candidates[cursor].second;
        if (!seen[idx]) {
            clusters.push_back(Cluster{idx});
            seen[idx] = true;
        }
    }

    return sphere_result_from_clusters(graph.Size(), std::move(clusters));
}

}  // namespace OECluster::detail

namespace OECluster {

namespace {

constexpr const char* SPHERE_NAME = "sphere_exclusion";

// Reads a complete condensed matrix serially.
class MatrixSphereSource {
public:
    MatrixSphereSource(const double* data, size_t n) : data_(data), n_(n) {}

    void Row(size_t center, const std::vector<size_t>& targets,
             std::vector<double>& distances) const {
        distances.resize(targets.size());
        for (size_t p = 0; p < targets.size(); ++p) {
            distances[p] = detail::dense_distance(data_, n_, center, targets[p]);
        }
    }

    void Nearest(const std::vector<size_t>& items,
                 const std::vector<size_t>& centers, std::vector<size_t>& best,
                 std::vector<double>& best_distance) const {
        best.resize(items.size());
        best_distance.resize(items.size());
        const auto read = [this](size_t a, size_t b) {
            return detail::dense_distance(data_, n_, a, b);
        };
        for (size_t p = 0; p < items.size(); ++p) {
            detail::nearest_center(items[p], centers, read, best[p],
                                   best_distance[p]);
        }
    }

private:
    const double* data_;
    size_t n_;
};

// Runs a comparison in parallel chunks. Each chunk writes only its own slots
// of the output vectors, which are sized before Run() starts; the engine
// reads them, and commits claims, on the caller thread after Run() returns.
class ComparisonSphereSource {
public:
    explicit ComparisonSphereSource(detail::ChunkedComparisons& work)
        : work_(work) {}

    void Row(size_t center, const std::vector<size_t>& targets,
             std::vector<double>& distances) {
        distances.resize(targets.size());
        work_.Run(targets.size(), [&](PairwiseComparison& local, size_t begin,
                                      size_t end) {
            for (size_t p = begin; p < end; ++p) {
                const size_t item = targets[p];
                distances[p] = local.Compare(std::min(center, item),
                                             std::max(center, item));
            }
        });
    }

    void Nearest(const std::vector<size_t>& items,
                 const std::vector<size_t>& centers, std::vector<size_t>& best,
                 std::vector<double>& best_distance) {
        best.resize(items.size());
        best_distance.resize(items.size());
        work_.Run(items.size(), [&](PairwiseComparison& local, size_t begin,
                                    size_t end) {
            const auto read = [&local](size_t a, size_t b) {
                return local.Compare(std::min(a, b), std::max(a, b));
            };
            for (size_t p = begin; p < end; ++p) {
                detail::nearest_center(items[p], centers, read, best[p],
                                       best_distance[p]);
            }
        });
    }

private:
    detail::ChunkedComparisons& work_;
};

void validate_sphere_options(const SphereExclusionOptions& options) {
    if (!std::isfinite(options.distance_threshold)) {
        throw std::invalid_argument(
            "sphere_exclusion distance_threshold must be finite");
    }
    if (options.distance_threshold < 0.0) {
        throw std::invalid_argument(
            "sphere_exclusion distance_threshold must be non-negative");
    }
    switch (options.order) {
        case SphereOrder::Input:
        case SphereOrder::Neighbors:
        case SphereOrder::Permutation:
            break;
        default:
            throw std::invalid_argument("Unknown sphere_exclusion order");
    }
    switch (options.assignment) {
        case SphereAssignment::First:
        case SphereAssignment::Nearest:
            break;
        default:
            throw std::invalid_argument("Unknown sphere_exclusion assignment");
    }
    if (options.reordering && options.order != SphereOrder::Neighbors) {
        throw std::invalid_argument(
            "sphere_exclusion reordering requires the Neighbors order");
    }
    if (!options.permutation.empty() &&
        options.order != SphereOrder::Permutation) {
        throw std::invalid_argument(
            "sphere_exclusion permutation requires the Permutation order");
    }
    detail::validate_chunk_size(options.chunk_size, SPHERE_NAME);
}

void validate_sphere_permutation(const SphereExclusionOptions& options,
                                 size_t n) {
    if (options.order != SphereOrder::Permutation) {
        return;
    }
    if (options.permutation.size() != n) {
        throw std::invalid_argument(
            "sphere_exclusion permutation has " +
            std::to_string(options.permutation.size()) + " entries for " +
            std::to_string(n) + " items");
    }
    std::vector<bool> used(n, false);
    for (const size_t index : options.permutation) {
        if (index >= n) {
            throw std::invalid_argument("sphere_exclusion permutation entry " +
                                        std::to_string(index) +
                                        " is outside the item range");
        }
        if (used[index]) {
            throw std::invalid_argument(
                "sphere_exclusion permutation repeats item " +
                std::to_string(index));
        }
        used[index] = true;
    }
}

// The clustering convention (butina, dbscan): no items give an empty result,
// and one item is its own center.
SphereExclusionResult small_sphere_result(size_t n) {
    if (n == 1) {
        return SphereExclusionResult(std::vector<ClusterLabel>{0},
                                     Clusters{Cluster{0}},
                                     std::vector<size_t>{0});
    }
    return SphereExclusionResult();
}

SphereExclusionResult to_public_result(detail::SphereEngineResult result) {
    return SphereExclusionResult(std::move(result.labels),
                                 std::move(result.clusters),
                                 std::move(result.centers));
}

const std::vector<size_t>* permutation_or_null(
    const SphereExclusionOptions& options) {
    return options.order == SphereOrder::Permutation ? &options.permutation
                                                     : nullptr;
}

}  // namespace

SphereExclusionResult sphere_exclusion(const StorageBackend& storage,
                                       const SphereExclusionOptions& options) {
    validate_sphere_options(options);
    detail::validate_complete_distance_storage(storage, SPHERE_NAME);
    const size_t n = storage.NumSamples();
    validate_sphere_permutation(options, n);
    if (n < 2) {
        return small_sphere_result(n);
    }

    const double* data = storage.Data();
    const MatrixSphereSource source(data, n);
    detail::SphereEngineResult result;
    if (options.order == SphereOrder::Neighbors) {
        // The graph would read a NaN as "not a neighbor" and carry on, while
        // the ordered sources raise on the same matrix; scanning first makes
        // every order refuse it.
        size_t k = 0;
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = i + 1; j < n; ++j, ++k) {
                if (!std::isfinite(data[k])) {
                    throw detail::sphere_non_finite_error(i, j);
                }
            }
        }
        ThresholdGraphOptions graph_options;
        graph_options.threshold = options.distance_threshold;
        // Capped at the item count, as ChunkedComparisons and k-medoids cap
        // theirs: the graph builder parallelizes over pairs, so a huge request
        // could otherwise launch pair-count-scale workers.
        graph_options.num_threads =
            detail::capped_threads(options.num_threads, n);
        graph_options.chunk_size = options.chunk_size;
        result = detail::sphere_neighbors_first(
            BuildThresholdNeighborGraph(storage, graph_options),
            options.reordering);
    } else {
        result = detail::sphere_ordered_first(source, n,
                                              permutation_or_null(options),
                                              options.distance_threshold);
    }
    if (options.assignment == SphereAssignment::Nearest) {
        detail::sphere_assign_nearest(source, n, result);
    }
    return to_public_result(std::move(result));
}

SphereExclusionResult sphere_exclusion(PairwiseComparison& comparison,
                                       const SphereExclusionOptions& options) {
    validate_sphere_options(options);
    detail::validate_comparison_facts(comparison, SPHERE_NAME);
    if (options.order == SphereOrder::Neighbors) {
        throw std::invalid_argument(
            "sphere_exclusion with the Neighbors order needs every pairwise "
            "distance; pass a precomputed distance matrix");
    }
    const size_t n = comparison.Size();
    validate_sphere_permutation(options, n);
    if (n < 2) {
        return small_sphere_result(n);
    }

    detail::ChunkedComparisons work(comparison, n, options.num_threads,
                                    options.chunk_size);
    ComparisonSphereSource source(work);
    detail::SphereEngineResult result = detail::sphere_ordered_first(
        source, n, permutation_or_null(options), options.distance_threshold);
    if (options.assignment == SphereAssignment::Nearest) {
        detail::sphere_assign_nearest(source, n, result);
    }
    return to_public_result(std::move(result));
}

}  // namespace OECluster
