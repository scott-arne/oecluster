/**
 * @file SphereExclusion.h
 * @brief Sphere-exclusion clustering: leader, Butina and DISE orders.
 */

#ifndef OECLUSTER_CLUSTERING_SPHEREEXCLUSION_H
#define OECLUSTER_CLUSTERING_SPHEREEXCLUSION_H

#include <cstddef>
#include <string>
#include <utility>
#include <vector>

#include "oecluster/PairwiseComparison.h"
#include "oecluster/StorageBackend.h"
#include "oecluster/clustering/ClusterTypes.h"

namespace OECluster {

/** @brief Order in which sphere exclusion takes centers. */
enum class SphereOrder {
    Input,        ///< Ascending item index (leader clustering).
    Neighbors,    ///< Descending neighbor count (Butina).
    Permutation,  ///< SphereExclusionOptions::permutation (DISE).
};

/** @brief How non-center items are assigned to centers. */
enum class SphereAssignment {
    First,    ///< The first sphere that reaches the item.
    Nearest,  ///< The nearest center; ties to the earlier center.
};

/** @brief Options for sphere_exclusion. */
struct SphereExclusionOptions {
    /// Items at or within this distance of a center join its sphere.
    double distance_threshold = 0.0;
    /// Order in which centers are taken.
    SphereOrder order = SphereOrder::Input;
    /// A complete permutation of the item indices; required with, and only
    /// with, SphereOrder::Permutation.
    std::vector<size_t> permutation;
    /// Recompute unclaimed neighbor counts after each cluster, as Butina's
    /// reordering does; SphereOrder::Neighbors only.
    bool reordering = false;
    /// How non-center items are assigned once the centers are fixed.
    SphereAssignment assignment = SphereAssignment::First;
    /// Worker threads: the threshold graph under the neighbor order, and
    /// the comparisons on a comparison; 0 selects the hardware concurrency.
    size_t num_threads = 0;
    /// Pairs or items per work unit; at least one.
    size_t chunk_size = 4096;
    /// Comparison overload under SphereOrder::Neighbors only: the most
    /// memory the threshold graph may take, in bytes. 0 applies the default
    /// limit, the larger of the condensed matrix the graph replaces and
    /// 1 GiB. Any other value is refused everywhere else.
    size_t max_graph_bytes = 0;
};

/**
 * @brief Sphere-exclusion result: labels, members and one center per cluster.
 *
 * Clusters are in center order and list their center first, then the other
 * members in ascending index. Every item is assigned, and an item's label is
 * its cluster's position.
 */
class SphereExclusionResult : public ClusteringResult {
public:
    SphereExclusionResult() = default;
    SphereExclusionResult(std::vector<ClusterLabel> labels, Clusters clusters,
                          std::vector<size_t> centers)
        : ClusteringResult(std::move(labels), std::move(clusters)),
          centers_(std::move(centers)) {}

    /// One center per cluster, in cluster order; Centers()[i] is
    /// Members()[i][0].
    const std::vector<size_t>& Centers() const { return centers_; }
    std::string Method() const override { return "sphere_exclusion"; }

private:
    std::vector<size_t> centers_;
};

/**
 * @brief Sphere exclusion over a precomputed distance matrix.
 *
 * Takes centers in the configured order. Each center claims every unclaimed
 * item j with d(center, j) <= distance_threshold. Under
 * SphereOrder::Neighbors with SphereAssignment::First the result equals
 * butina_cluster() with the same threshold and reordering. Nearest
 * assignment then moves each non-center item to its nearest center, with
 * ties going to the earlier center; the centers do not change.
 *
 * :param storage: Complete dense or memory-mapped distance storage.
 * :param options: Threshold, order, assignment and threading options.
 * :returns: The clustering; empty for zero items.
 * :raises std::invalid_argument: On a non-finite or negative threshold, an
 *     unknown order or assignment, reordering without the neighbor order, a
 *     permutation that is not a complete permutation of the items (or is
 *     given with another order), a zero chunk_size, a non-zero
 *     max_graph_bytes, or sparse or data-less storage.
 * :raises std::runtime_error: If a distance read is NaN or infinite. Under
 *     the neighbor order every distance is read up front.
 */
SphereExclusionResult sphere_exclusion(const StorageBackend& storage,
                                       const SphereExclusionOptions& options);

/**
 * @brief Sphere exclusion over a comparison, evaluated lazily.
 *
 * Under SphereOrder::Input and SphereOrder::Permutation each center
 * compares against the still-unclaimed items in parallel chunks, as
 * Compare(min, max). Under SphereOrder::Neighbors the threshold graph is
 * built in two passes over every pair, one to count each item's neighbors
 * and one to record them, at 16 bytes per within-threshold pair plus 16 per
 * item on a 64-bit platform; its exact size is known before it is allocated,
 * so a graph above max_graph_bytes (or the default limit) is refused rather
 * than attempted.
 * Nearest assignment compares every non-center item with every center. The
 * result equals the matrix overload's on the same distances, for every
 * num_threads and chunk_size.
 *
 * Precondition under SphereOrder::Neighbors: Compare(i, j) returns a
 * bit-identical value for a pair on every call and every clone. A comparison
 * that changes a row's neighbor count between the passes is refused; one
 * that swaps neighbors while keeping every count is outside the contract.
 *
 * :param comparison: Distance comparison; cloned once per running chunk.
 * :param options: Threshold, order, assignment and threading options.
 * :returns: The clustering; empty for zero items.
 * :raises std::invalid_argument: On the matrix overload's option refusals
 *     other than max_graph_bytes, or a non-zero max_graph_bytes with an
 *     order other than SphereOrder::Neighbors.
 * :raises ComparisonError: If the comparison's facts rule out ranking its
 *     distances; or, under SphereOrder::Neighbors, if it is a ROCS
 *     comparison, whose scores depend on what its overlay scored before.
 * :raises std::runtime_error: If a comparison returns NaN or infinity.
 * :raises std::length_error: If the threshold graph would exceed its limit.
 * :raises std::logic_error: If the graph's two passes disagree on a row's
 *     size.
 */
SphereExclusionResult sphere_exclusion(PairwiseComparison& comparison,
                                       const SphereExclusionOptions& options);

}  // namespace OECluster

#endif  // OECLUSTER_CLUSTERING_SPHEREEXCLUSION_H
