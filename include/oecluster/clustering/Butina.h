/**
 * @file Butina.h
 * @brief Butina clustering over a distance matrix or a comparison.
 */

#ifndef OECLUSTER_CLUSTERING_BUTINA_H
#define OECLUSTER_CLUSTERING_BUTINA_H

#include <cstddef>
#include <string>

#include "oecluster/PairwiseComparison.h"
#include "oecluster/StorageBackend.h"
#include "oecluster/clustering/ClusterTypes.h"

namespace OECluster {

/**
 * @brief Options for Butina clustering.
 */
struct ButinaOptions {
    double distance_threshold = 0.0;
    bool reordering = false;
    size_t num_threads = 0;
    size_t chunk_size = 4096;
    /// Comparison overload only: the most memory the threshold graph may
    /// take, in bytes. 0 applies the default limit, the larger of the
    /// condensed matrix the graph replaces and 1 GiB. The storage overload
    /// refuses any other value.
    size_t max_graph_bytes = 0;
};

/**
 * @brief Butina clustering result. Carries only labels and members; the first
 *     member of each cluster is the highest-neighborhood representative.
 */
class ButinaResult : public ClusteringResult {
public:
    using ClusteringResult::ClusteringResult;

    std::string Method() const override { return "butina"; }
};

/**
 * @brief Cluster a precomputed distance matrix with the Butina algorithm.
 *
 * :param storage: Pairwise distance storage.
 * :param options: Butina clustering options.
 * :returns: ButinaResult whose clusters are ordered with the
 *     highest-neighborhood representative first, and whose per-item labels
 *     equal each member's cluster position.
 * :raises std::invalid_argument: On a negative threshold or a non-zero
 *     max_graph_bytes.
 */
ButinaResult butina_cluster(const StorageBackend& storage, const ButinaOptions& options);

/**
 * @brief Cluster a comparison with the Butina algorithm, holding no matrix.
 *
 * Builds the threshold graph in two passes over every pair -- one to count
 * each item's neighbors, one to record them -- so the result equals the
 * storage overload's on a matrix filled through Compare(i, j), for every
 * num_threads and chunk_size. The graph costs 16 bytes per within-threshold
 * pair plus 16 per item on a 64-bit platform, and its exact size is known
 * before it is allocated, so a graph above the limit is refused rather than
 * attempted.
 * A chunk_size of 0 selects 4096, as on the storage overload.
 *
 * Precondition: Compare(i, j) returns a bit-identical value for a pair on
 * every call and every clone. A comparison that changes a row's neighbor
 * count between the passes is refused; one that swaps neighbors while
 * keeping every count is outside the contract.
 *
 * :param comparison: Distance comparison; cloned once per running chunk.
 * :param options: Butina clustering options.
 * :returns: As the storage overload.
 * :raises std::invalid_argument: On a negative threshold.
 * :raises ComparisonError: If the comparison reports similarities, a
 *     non-zero self-distance, or values that may be NaN; or if it is a ROCS
 *     comparison, whose scores depend on what its overlay scored before.
 * :raises std::runtime_error: If a comparison returns NaN or infinity.
 * :raises std::length_error: If the graph would exceed its limit.
 * :raises std::logic_error: If the two passes disagree on a row's size.
 */
ButinaResult butina_cluster(PairwiseComparison& comparison, const ButinaOptions& options);

}  // namespace OECluster

#endif  // OECLUSTER_CLUSTERING_BUTINA_H
