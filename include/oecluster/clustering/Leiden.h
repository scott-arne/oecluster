/**
 * @file Leiden.h
 * @brief Leiden community detection over a shared-nearest-neighbor graph.
 */

#ifndef OECLUSTER_CLUSTERING_LEIDEN_H
#define OECLUSTER_CLUSTERING_LEIDEN_H

#include <cstddef>
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

#include "oecluster/PairwiseComparison.h"
#include "oecluster/StorageBackend.h"
#include "oecluster/clustering/ClusterTypes.h"
#include "oecluster/clustering/KNNGraph.h"

namespace OECluster {

/** @brief Quality function optimized by leiden. */
enum class LeidenObjective {
    Modularity,  ///< Reichardt-Bornholdt modularity with a resolution.
    CPM,         ///< Constant Potts Model.
};

/** @brief Options for leiden. */
struct LeidenOptions {
    /// Neighbors per item for the graph (raw overloads only); 1 <= k <= n - 1.
    size_t k = 0;
    LeidenObjective objective = LeidenObjective::Modularity;
    /// Resolution; finite and >= 0.
    double resolution = 1.0;
    /// SNN Jaccard weights below this are dropped; 0 <= prune < 1.
    double prune = 1.0 / 15.0;
    /// Refinement randomness; finite and > 0.
    double theta = 0.01;
    /// -1 iterates until a pass changes nothing; otherwise >= 0 passes.
    int64_t n_iterations = -1;
    /// Seed for the random stream.
    uint64_t seed = 0;
    /// Worker threads for the graph and the SNN weights; 0 selects the
    /// hardware concurrency.
    size_t num_threads = 0;
    /// Pairwise distances per work unit for the graph; at least one.
    size_t chunk_size = 4096;
};

/**
 * @brief Leiden result: labels, members, the partition's quality, and the
 * options that shaped it.
 *
 * Clusters are ordered by their smallest member and list members in
 * ascending index. Every item is assigned; an item with no retained edge is
 * a singleton.
 */
class LeidenResult : public ClusteringResult {
public:
    LeidenResult() = default;
    LeidenResult(std::vector<ClusterLabel> labels, Clusters clusters,
                 double quality, size_t iterations, LeidenObjective objective,
                 double resolution, size_t k)
        : ClusteringResult(std::move(labels), std::move(clusters)),
          quality_(quality),
          iterations_(iterations),
          objective_(objective),
          resolution_(resolution),
          k_(k) {}

    /// Objective value of the returned partition on the SNN graph.
    double Quality() const { return quality_; }
    /// Full passes run, including a final pass that changed nothing.
    size_t Iterations() const { return iterations_; }
    LeidenObjective Objective() const { return objective_; }
    double Resolution() const { return resolution_; }
    /// Neighbors per item in the graph; 0 on a default-constructed result.
    size_t K() const { return k_; }
    std::string Method() const override { return "leiden"; }

private:
    double quality_ = 0.0;
    size_t iterations_ = 0;
    LeidenObjective objective_ = LeidenObjective::Modularity;
    double resolution_ = 0.0;
    size_t k_ = 0;
};

/**
 * @brief Leiden community detection over a prebuilt k-nearest-neighbor graph.
 *
 * The graph is reweighted by shared-nearest-neighbor Jaccard: with N+(i) row
 * i plus i itself, an edge {i, j} for j in row i or i in row j weighs
 * s / (2(k + 1) - s), s = |N+(i) intersect N+(j)|, and edges below
 * options.prune are dropped. The partition then maximizes options.objective
 * by the Leiden algorithm (Traag, Waltman and van Eck, 2019), so every
 * cluster induces a connected subgraph. The same options on the same build
 * give the same result, whatever num_threads is.
 *
 * :param graph: The neighbor graph; options.k and options.chunk_size are
 *     ignored.
 * :param options: Objective, resolution, prune, theta, n_iterations, seed
 *     and num_threads (for the weights).
 * :returns: The clustering; empty for a zero-item graph.
 * :raises std::invalid_argument: On an unknown objective, a resolution that
 *     is negative or not finite, a prune outside [0, 1), a theta that is not
 *     a positive finite number, n_iterations < -1, or more than INT_MAX items.
 */
LeidenResult leiden(const KNNGraph& graph, const LeidenOptions& options);

/**
 * @brief Leiden community detection over a precomputed distance matrix.
 *
 * Validates every option before building the graph, then builds it as
 * knn_graph() does and clusters it as the graph overload does.
 *
 * :param storage: Dense, memory-mapped or finalized sparse storage.
 * :param options: Graph, weighting and optimization options.
 * :returns: The clustering; empty for zero items.
 * :raises std::invalid_argument: On the graph overload's refusals or
 *     knn_graph()'s.
 * :raises std::runtime_error: If a distance read is NaN or infinite.
 */
LeidenResult leiden(const StorageBackend& storage, const LeidenOptions& options);

/**
 * @brief Leiden community detection over a comparison, evaluated lazily.
 *
 * :param comparison: Distance comparison; cloned once per running unit.
 * :param options: Graph, weighting and optimization options.
 * :returns: The clustering; empty for zero items.
 * :raises std::invalid_argument: On the graph overload's refusals or
 *     knn_graph()'s; all are checked before any comparison runs.
 * :raises ComparisonError: If the comparison's facts rule out ranking.
 * :raises std::runtime_error: If a comparison returns NaN or infinity.
 */
LeidenResult leiden(PairwiseComparison& comparison, const LeidenOptions& options);

}  // namespace OECluster

#endif  // OECLUSTER_CLUSTERING_LEIDEN_H
