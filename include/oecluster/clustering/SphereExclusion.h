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

#include "oecluster/StorageBackend.h"
#include "oecluster/clustering/ClusterTypes.h"

namespace OECluster {

/** @brief Order in which sphere exclusion takes centers. */
enum class SphereOrder {
    Input,        ///< Ascending item index (leader clustering).
    Neighbors,    ///< Descending neighbor count (Butina); matrix only.
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
 *     given with another order), a zero chunk_size, or sparse or data-less
 *     storage.
 * :raises std::runtime_error: If a distance read is NaN or infinite. Under
 *     the neighbor order every distance is read up front.
 */
SphereExclusionResult sphere_exclusion(const StorageBackend& storage,
                                       const SphereExclusionOptions& options);

}  // namespace OECluster

#endif  // OECLUSTER_CLUSTERING_SPHEREEXCLUSION_H
