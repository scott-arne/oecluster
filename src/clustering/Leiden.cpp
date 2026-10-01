/**
 * @file Leiden.cpp
 * @brief Leiden community detection over a KNNGraph and over raw input.
 */

#include "oecluster/clustering/Leiden.h"

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "DiversityValidation.h"
#include "KNNGraphBuild.h"
#include "LeidenEngine.h"
#include "SNNWeights.h"

namespace OECluster {

namespace {

constexpr const char* LEIDEN_NAME = "leiden";

// A forged enum value falls through the switch instead of reaching the
// engine's two-way branches as an unchecked third case.
bool known_objective(LeidenObjective objective) {
    switch (objective) {
        case LeidenObjective::Modularity:
        case LeidenObjective::CPM:
            return true;
    }
    return false;
}

void validate_options(const LeidenOptions& options) {
    const std::string name = LEIDEN_NAME;
    if (!known_objective(options.objective)) {
        throw std::invalid_argument(name +
                                    " objective is not a known LeidenObjective");
    }
    if (!std::isfinite(options.resolution) || options.resolution < 0.0) {
        throw std::invalid_argument(name +
                                    " resolution must be finite and non-negative");
    }
    if (!std::isfinite(options.prune) || options.prune < 0.0 ||
        options.prune >= 1.0) {
        throw std::invalid_argument(name + " prune must be finite and in [0, 1)");
    }
    if (!std::isfinite(options.theta) || options.theta <= 0.0) {
        throw std::invalid_argument(name + " theta must be finite and positive");
    }
    if (options.n_iterations < -1) {
        throw std::invalid_argument(name +
                                    " n_iterations must be -1 or non-negative, got " +
                                    std::to_string(options.n_iterations));
    }
}

// Raw-input validation: chunk_size, the options, zero items, the item
// count, then k, so no refusal costs a comparison. Returns false for zero
// items.
bool validate_raw(size_t n, const LeidenOptions& options) {
    detail::validate_chunk_size(options.chunk_size, LEIDEN_NAME);
    validate_options(options);
    if (n == 0) {
        return false;
    }
    detail::validate_leiden_item_count(n, LEIDEN_NAME);
    detail::validate_knn_k(n, options.k, LEIDEN_NAME);
    return true;
}

KNNGraphOptions graph_options(const LeidenOptions& options) {
    KNNGraphOptions graph;
    graph.k = options.k;
    graph.num_threads = options.num_threads;
    graph.chunk_size = options.chunk_size;
    return graph;
}

LeidenResult empty_result(const LeidenOptions& options, size_t k) {
    return LeidenResult({}, {}, 0.0, 0, options.objective, options.resolution, k);
}

// The weights arrive by value and move into the engine, so neither a second
// CSR nor, on the raw overloads, the kNN graph is alive during optimization.
LeidenResult optimize(detail::WeightedGraph weighted, size_t k,
                      const LeidenOptions& options) {
    detail::LeidenParams params;
    params.objective = options.objective;
    params.resolution = options.resolution;
    params.theta = options.theta;
    detail::LeidenRun run = detail::run_leiden(
        std::move(weighted), params, options.n_iterations, options.seed);
    Clusters clusters = labels_to_clusters(run.labels);
    return LeidenResult(std::move(run.labels), std::move(clusters), run.quality,
                        run.iterations, options.objective, options.resolution, k);
}

}  // namespace

LeidenResult leiden(const KNNGraph& graph, const LeidenOptions& options) {
    validate_options(options);
    const size_t n = graph.NumItems();
    if (n == 0) {
        return empty_result(options, graph.K());
    }
    detail::validate_leiden_item_count(n, LEIDEN_NAME);
    return optimize(detail::snn_weights(graph, options.prune, options.num_threads),
                    graph.K(), options);
}

LeidenResult leiden(const StorageBackend& storage, const LeidenOptions& options) {
    if (!validate_raw(storage.NumSamples(), options)) {
        return empty_result(options, options.k);
    }
    // The temporary graph is destroyed once its weights are built.
    detail::WeightedGraph weighted = detail::snn_weights(
        detail::build_knn_graph(storage, graph_options(options), LEIDEN_NAME),
        options.prune, options.num_threads);
    return optimize(std::move(weighted), options.k, options);
}

LeidenResult leiden(PairwiseComparison& comparison, const LeidenOptions& options) {
    if (!validate_raw(comparison.Size(), options)) {
        return empty_result(options, options.k);
    }
    detail::WeightedGraph weighted = detail::snn_weights(
        detail::build_knn_graph(comparison, graph_options(options), LEIDEN_NAME),
        options.prune, options.num_threads);
    return optimize(std::move(weighted), options.k, options);
}

}  // namespace OECluster
