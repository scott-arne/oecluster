/**
 * @file StreamingCoreDistances.h
 * @brief HDBSCAN core distances from a matrix or from one pass over a comparison.
 */

#ifndef OECLUSTER_SRC_CLUSTERING_STREAMINGCOREDISTANCES_H
#define OECLUSTER_SRC_CLUSTERING_STREAMINGCOREDISTANCES_H

#include <cstddef>
#include <string>
#include <utility>
#include <vector>

#include "oecluster/PairwiseComparison.h"
#include "oecluster/StorageBackend.h"

namespace OECluster::detail {

/** @brief Core distances and the largest distance read while computing them. */
struct CoreDistances {
    std::vector<double> values;
    // 0 when no pair was read (min_samples = 1).
    double max_distance = 0.0;
};

/**
 * @brief The blocks and rounds of the one-pass core schedule.
 *
 * Block b holds items [bounds[b], bounds[b + 1]). rounds[0] pairs every block
 * with itself; every later round pairs blocks by the circle method, so the
 * tiles of one round touch disjoint blocks and every unordered pair of items
 * falls in exactly one tile of one round.
 */
struct CoreBlockSchedule {
    std::vector<size_t> bounds;
    std::vector<std::vector<std::pair<size_t, size_t>>> rounds;
};

/**
 * @brief The schedule for n items and the given participant count.
 *
 * :param n: Item count.
 * :param participants: Workers, at least 1 when n > 0.
 * :returns: min(n, 4 x participants) blocks; an empty schedule when n is 0.
 */
CoreBlockSchedule core_block_schedule(size_t n, size_t participants);

/**
 * @brief Core distances from one pass over a comparison, each pair compared once.
 *
 * Item i's core distance is its (min_samples - 1)-th smallest distance to
 * another item, the value matrix_core_distances takes over a matrix filled
 * through the same Compare(min(i, j), max(i, j)). Every distance read must be
 * finite and non-negative; a zero is read as +0.0.
 *
 * :param comparison: The comparison; cloned once per running tile.
 * :param min_samples: Self-inclusive neighbor count, 1 to comparison.Size().
 * :param num_threads: Worker threads; 0 selects the hardware count.
 * :param caller: Entry point name for error messages.
 * :raises std::invalid_argument: If min_samples is 0 or above the item count.
 * :raises std::length_error: If the per-item heaps' size overflows a size_t.
 * :raises std::runtime_error: For a non-finite or negative distance, naming the pair.
 */
CoreDistances streaming_core_distances(PairwiseComparison& comparison,
                                       size_t min_samples, size_t num_threads,
                                       const std::string& caller);

/**
 * @brief Core distances from a dense or memory-mapped matrix, row by row.
 *
 * The rows include the item's own distance as 0, as in 5.19.0. Every distance
 * read must be finite and non-negative; a zero is read as +0.0. When
 * min_samples is 1 nothing is read.
 *
 * :param chunk_size: Ceiling on rows per work unit; 0 selects 64.
 * :raises std::invalid_argument: If the storage is sparse or not contiguous, or
 *     min_samples is 0 or above the item count.
 * :raises std::runtime_error: For a non-finite or negative distance, naming the pair.
 */
CoreDistances matrix_core_distances(const StorageBackend& storage,
                                    size_t min_samples, size_t num_threads,
                                    size_t chunk_size, const std::string& caller);

}  // namespace OECluster::detail

#endif  // OECLUSTER_SRC_CLUSTERING_STREAMINGCOREDISTANCES_H
