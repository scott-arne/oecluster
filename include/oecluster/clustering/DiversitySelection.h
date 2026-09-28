/**
 * @file DiversitySelection.h
 * @brief Farthest-first subset selection and the #Circles coverage measure.
 *
 * Both entry points run on one deterministic kernel and take either a
 * precomputed distance matrix or a comparison evaluated lazily, so a library
 * too large for a full matrix can still be subset. Neither is a clustering:
 * each result names a subset, not a partition.
 */

#ifndef OECLUSTER_CLUSTERING_DIVERSITYSELECTION_H
#define OECLUSTER_CLUSTERING_DIVERSITYSELECTION_H

#include <cstddef>
#include <limits>
#include <vector>

#include "oecluster/PairwiseComparison.h"
#include "oecluster/StorageBackend.h"

namespace OECluster {

/** @brief Where maxmin_select starts when no initial selection is given. */
enum class MaxMinSeed {
    Index,     ///< MaxMinOptions::seed.
    Medoid,    ///< The global medoid; the distance-matrix overload only.
    Farthest,  ///< The item farthest from item 0, ties to the smaller index.
};

/** @brief Why maxmin_select stopped. */
enum class MaxMinStop {
    Count,      ///< The selection reached MaxMinOptions::count.
    Threshold,  ///< The best remaining candidate was at or within the threshold.
    Exhausted,  ///< No unselected item remained.
};

/** @brief How circles builds its packing. */
enum class CirclesMethod {
    MaxMin,      ///< Farthest-first from item 0.
    Sequential,  ///< The reference greedy pass over input order.
};

/** @brief Options for maxmin_select. */
struct MaxMinOptions {
    /// Total selection size, initial entries included; 0 means no count limit.
    size_t count = 0;
    /// Stop before a candidate at or within this distance of the selection.
    /// NaN means unset; at least one of count and threshold is required.
    double threshold = std::numeric_limits<double>::quiet_NaN();
    /// Where to start when initial is empty.
    MaxMinSeed seed_mode = MaxMinSeed::Index;
    /// The starting item; read only when seed_mode is Index.
    size_t seed = 0;
    /// An existing selection to extend. Non-empty only with the default seed.
    std::vector<size_t> initial;
    /// Worker threads for the comparison overload; 0 auto-detects hardware
    /// concurrency. An explicit value is capped at the item count.
    size_t num_threads = 0;
    /// Items per work unit, with the same scope as num_threads; at least one.
    size_t chunk_size = 256;
};

/** @brief The result of maxmin_select. */
struct MaxMinSelection {
    /// Selected items in selection order, initial entries first.
    std::vector<size_t> indices;
    /// Each pick's distance to the earlier selection when it was picked; NaN
    /// for the seed and every initial entry. Never increases after them.
    std::vector<double> pick_distances;
    /// Why the selection stopped.
    MaxMinStop stop = MaxMinStop::Exhausted;
};

/** @brief Options for circles. */
struct CirclesOptions {
    /// Packing construction.
    CirclesMethod method = CirclesMethod::MaxMin;
    /// Worker threads for the comparison overload; 0 auto-detects hardware
    /// concurrency. An explicit value is capped at the item count.
    size_t num_threads = 0;
    /// Items per work unit, with the same scope as num_threads; at least one.
    size_t chunk_size = 256;
};

/**
 * @brief The result of circles: a packing whose members are pairwise more
 * than the threshold apart.
 *
 * Any valid packing is a lower bound on the packing number, so count is a
 * lower bound under either method, and the two methods can disagree.
 */
struct CirclesResult {
    /// Number of members; the #Circles value.
    size_t count = 0;
    /// Members in pick order (MaxMin) or input order (Sequential).
    std::vector<size_t> members;
    /// The threshold the packing was built at.
    double threshold = std::numeric_limits<double>::quiet_NaN();
    /// The method that built it.
    CirclesMethod method = CirclesMethod::MaxMin;
};

/**
 * @brief Farthest-first (MaxMin) selection over a precomputed distance matrix.
 *
 * :param storage: Complete dense or memory-mapped distances.
 * :param options: Stop conditions, seed and threading.
 * :returns: The selection, its pick distances and the stop reason.
 * :raises std::invalid_argument: On an out-of-range enum, incomplete storage,
 *     no items, a zero chunk_size, an invalid count, threshold, seed or
 *     initial set, a non-finite distance read, or, with the Medoid seed, any
 *     non-finite entry or a row sum that overflows.
 */
MaxMinSelection maxmin_select(const StorageBackend& storage,
                              const MaxMinOptions& options);

/**
 * @brief The #Circles packing of a precomputed distance matrix.
 *
 * :param storage: Complete dense or memory-mapped distances.
 * :param threshold: Members are pairwise strictly farther apart than this.
 * :param options: Method and threading.
 * :returns: The packing.
 * :raises std::invalid_argument: On an out-of-range method, incomplete
 *     storage, no items, a zero chunk_size, a NaN, infinite or negative
 *     threshold, or a non-finite distance read.
 */
CirclesResult circles(const StorageBackend& storage, double threshold,
                      const CirclesOptions& options);

}  // namespace OECluster

#endif  // OECLUSTER_CLUSTERING_DIVERSITYSELECTION_H
