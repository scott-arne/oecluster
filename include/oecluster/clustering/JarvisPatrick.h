/**
 * @file JarvisPatrick.h
 * @brief Jarvis-Patrick shared-nearest-neighbor clustering.
 */

#ifndef OECLUSTER_CLUSTERING_JARVISPATRICK_H
#define OECLUSTER_CLUSTERING_JARVISPATRICK_H

#include <cstddef>
#include <string>
#include <utility>
#include <vector>

#include "oecluster/PairwiseComparison.h"
#include "oecluster/StorageBackend.h"
#include "oecluster/clustering/ClusterTypes.h"
#include "oecluster/clustering/KNNGraph.h"

namespace OECluster {

/** @brief Options for jarvis_patrick over raw input. */
struct JarvisPatrickOptions {
    /// Neighbors per item for the graph; 1 <= k <= n - 1 for n >= 1 items.
    size_t k = 0;
    /// Shared neighbors required to link a mutual pair; kmin < k.
    size_t kmin = 0;
    /// Worker threads for the graph; 0 selects the hardware concurrency.
    size_t num_threads = 0;
    /// Pairwise distances per work unit for the graph; at least one.
    size_t chunk_size = 4096;
};

/**
 * @brief Jarvis-Patrick result: labels, members, and the k and kmin used.
 *
 * Clusters are ordered by their smallest member and list members in
 * ascending index. Every item is assigned; an unlinked item is a singleton.
 */
class JarvisPatrickResult : public ClusteringResult {
public:
    JarvisPatrickResult() = default;
    JarvisPatrickResult(std::vector<ClusterLabel> labels, Clusters clusters,
                        size_t k, size_t kmin)
        : ClusteringResult(std::move(labels), std::move(clusters)),
          k_(k),
          kmin_(kmin) {}

    /// Neighbors per item in the graph; 0 on a default-constructed result.
    size_t K() const { return k_; }
    /// Shared neighbors required to link; 0 on a default-constructed result.
    size_t KMin() const { return kmin_; }
    std::string Method() const override { return "jarvis_patrick"; }

private:
    size_t k_ = 0;
    size_t kmin_ = 0;
};

/**
 * @brief Jarvis-Patrick clustering over a prebuilt k-nearest-neighbor graph.
 *
 * Items i and j are linked iff each is in the other's row and the two rows
 * share at least kmin items. Clusters are the connected components of the
 * links. Rows exclude their own item, so a self-inclusive formulation's k is
 * this k + 1 and its kmin for a mutual pair is this kmin + 2.
 *
 * :param graph: The neighbor graph.
 * :param kmin: Shared neighbors required to link a mutual pair.
 * :returns: The clustering; empty for a zero-item graph whatever kmin is.
 * :raises std::invalid_argument: If kmin >= graph.K() on a non-empty graph.
 */
JarvisPatrickResult jarvis_patrick(const KNNGraph& graph, size_t kmin);

/**
 * @brief Jarvis-Patrick clustering over a precomputed distance matrix.
 *
 * Checks kmin < k before building the graph, then builds it as knn_graph()
 * does and clusters it as the graph overload does.
 *
 * :param storage: Dense, memory-mapped or finalized sparse storage.
 * :param options: k, kmin and threading options.
 * :returns: The clustering; empty for zero items.
 * :raises std::invalid_argument: On knn_graph()'s refusals, or kmin >= k.
 * :raises std::runtime_error: If a distance read is NaN or infinite.
 */
JarvisPatrickResult jarvis_patrick(const StorageBackend& storage,
                                   const JarvisPatrickOptions& options);

/**
 * @brief Jarvis-Patrick clustering over a comparison, evaluated lazily.
 *
 * :param comparison: Distance comparison; cloned once per running unit.
 * :param options: k, kmin and threading options.
 * :returns: The clustering; empty for zero items.
 * :raises std::invalid_argument: On knn_graph()'s refusals, or kmin >= k;
 *     both are checked before any comparison runs.
 * :raises ComparisonError: If the comparison's facts rule out ranking.
 * :raises std::runtime_error: If a comparison returns NaN or infinity.
 */
JarvisPatrickResult jarvis_patrick(PairwiseComparison& comparison,
                                   const JarvisPatrickOptions& options);

}  // namespace OECluster

#endif  // OECLUSTER_CLUSTERING_JARVISPATRICK_H
