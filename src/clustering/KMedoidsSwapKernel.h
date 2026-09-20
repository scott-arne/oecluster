/**
 * @file KMedoidsSwapKernel.h
 * @brief FastPAM1 swap scan, verification pass and the assignment cache they share.
 */

#ifndef OECLUSTER_SRC_CLUSTERING_KMEDOIDSSWAPKERNEL_H
#define OECLUSTER_SRC_CLUSTERING_KMEDOIDSSWAPKERNEL_H

#include <algorithm>
#include <cstddef>
#include <limits>
#include <vector>

#include "DistanceAccess.h"
#include "oecluster/ThreadPool.h"

namespace OECluster::detail {

/**
 * @brief Cached nearest and second-nearest medoid facts for one item.
 */
struct Assignment {
    size_t nearest_slot = 0;          ///< Slot index, not an item index.
    double nearest_distance = 0.0;
    double second_nearest_distance =  ///< Where the item lands if its own medoid leaves.
        std::numeric_limits<double>::infinity();
};

/**
 * @brief Rewrite ``slot_of`` so every medoid item maps to its slot.
 *
 * Non-medoid items are marked with ``slot_of.size()``, which is the item count,
 * so "is h a medoid" is a single comparison in the scans below.
 *
 * :param slot_of: Item-indexed slot map, sized to the item count, overwritten.
 * :param medoids: Current medoid item indices, one per slot.
 */
inline void refresh_slot_map(std::vector<size_t>& slot_of,
                             const std::vector<size_t>& medoids) {
    std::fill(slot_of.begin(), slot_of.end(), slot_of.size());
    for (size_t slot = 0; slot < medoids.size(); ++slot) {
        slot_of[medoids[slot]] = slot;
    }
}

/**
 * @brief Clamp a chunk size to the item count.
 *
 * ThreadPool derives its chunk count as ``(range + chunk_size - 1) /
 * chunk_size``, which wraps to zero for a chunk size near SIZE_MAX and then
 * runs no chunk at all: every scan below would return its default-initialized
 * output and the call would report a zero-cost single cluster as converged.
 * Clamping to the item count is observationally free -- a chunk at least as
 * wide as the range is one chunk either way -- and keeps both this header's own
 * ceiling arithmetic and the pool's below the overflow. Validation guarantees
 * ``chunk_size >= 1``, and every call site has ``n >= n_clusters >= 1``, so the
 * result is never zero.
 *
 * :param n: Number of items.
 * :param chunk_size: Caller-requested chunk size; at least one.
 * :returns: The chunk size to hand ThreadPool, never zero and never above n.
 */
inline size_t effective_chunk_size(size_t n, size_t chunk_size) {
    return std::min(chunk_size, n);
}

/**
 * @brief Build the nearest and second-nearest cache for every item.
 *
 * :param data: Condensed distance array from ``StorageBackend::Data()``.
 * :param n: Number of items.
 * :param medoids: Current medoid item indices, one per slot; must be distinct.
 * :param num_threads: Worker threads; 0 auto-detects hardware concurrency.
 * :param chunk_size: Chunk size for the parallel scan; at least one.
 * :returns: One Assignment per item, indexed by item.
 */
inline std::vector<Assignment> build_assignments(
    const double* data, size_t n, const std::vector<size_t>& medoids,
    size_t num_threads, size_t chunk_size) {
    const size_t k = medoids.size();
    const double infinity = std::numeric_limits<double>::infinity();

    std::vector<size_t> slot_of(n, n);
    refresh_slot_map(slot_of, medoids);

    std::vector<Assignment> assignments(n);
    ThreadPool pool(num_threads);
    pool.ParallelFor(0, n, effective_chunk_size(n, chunk_size),
                     [&](size_t begin, size_t end) {
        for (size_t j = begin; j < end; ++j) {
            Assignment entry;

            if (slot_of[j] != n) {
                // The self-assignment rule: a medoid always belongs to its own
                // slot. Ranking it against the others would hand two duplicate
                // medoids at distance 0 to the same slot and leave the other
                // cluster empty, breaking the exactly-k guarantee.
                entry.nearest_slot = slot_of[j];
                entry.nearest_distance = 0.0;
                entry.second_nearest_distance = infinity;
                for (size_t slot = 0; slot < k; ++slot) {
                    if (slot == entry.nearest_slot) {
                        continue;
                    }
                    const double distance =
                        dense_distance(data, n, j, medoids[slot]);
                    if (distance < entry.second_nearest_distance) {
                        entry.second_nearest_distance = distance;
                    }
                }
            } else {
                entry.nearest_slot = 0;
                entry.nearest_distance = infinity;
                entry.second_nearest_distance = infinity;
                for (size_t slot = 0; slot < k; ++slot) {
                    const double distance =
                        dense_distance(data, n, j, medoids[slot]);
                    // Slots are visited in slot order, but the tie rule ranks
                    // on the medoid's item index, so an equal distance only
                    // displaces the incumbent when its item index is smaller.
                    const bool wins =
                        distance < entry.nearest_distance ||
                        (distance == entry.nearest_distance &&
                         medoids[slot] < medoids[entry.nearest_slot]);
                    if (wins) {
                        entry.second_nearest_distance = entry.nearest_distance;
                        entry.nearest_distance = distance;
                        entry.nearest_slot = slot;
                    } else if (distance < entry.second_nearest_distance) {
                        entry.second_nearest_distance = distance;
                    }
                }
            }

            assignments[j] = entry;
        }
    });

    return assignments;
}

/**
 * @brief Sum every item's distance to its assigned medoid.
 *
 * :param assignments: Cache from ``build_assignments``.
 * :returns: The k-medoids objective for the configuration the cache describes.
 */
inline double total_cost(const std::vector<Assignment>& assignments) {
    // Ascending item order, so the sum is bit-reproducible whatever the
    // chunking that filled the cache.
    double cost = 0.0;
    for (const Assignment& entry : assignments) {
        cost += entry.nearest_distance;
    }
    return cost;
}

/**
 * @brief One evaluated (leaving slot, entering item) swap.
 *
 * ``score`` is the predicted delta in the fast loop and the recomputed total
 * in the verification pass. Both are minimized under the same key, so one
 * comparator serves both.
 */
struct SwapCandidate {
    double score = 0.0;
    size_t entering_item = 0;
    size_t leaving_item = 0;
    size_t leaving_slot = 0;
    bool valid = false;
};

/**
 * @brief Rank two candidates by (score, entering item, leaving medoid item).
 *
 * The key is a strict total order over distinct candidate pairs, so reducing
 * chunk-local winners in any order gives the same answer.
 *
 * :param candidate: Challenger; an invalid candidate never wins.
 * :param incumbent: Current best; an invalid incumbent always loses.
 * :returns: True when the challenger should displace the incumbent.
 */
inline bool is_better(const SwapCandidate& candidate,
                      const SwapCandidate& incumbent) {
    if (!candidate.valid) {
        return false;
    }
    if (!incumbent.valid) {
        return true;
    }
    if (candidate.score != incumbent.score) {
        return candidate.score < incumbent.score;
    }
    if (candidate.entering_item != incumbent.entering_item) {
        return candidate.entering_item < incumbent.entering_item;
    }
    return candidate.leaving_item < incumbent.leaving_item;
}

/**
 * @brief Best swap over every candidate, scored by the FastPAM1 predicted delta.
 *
 * The returned ``score`` is a **delta** -- the predicted change in the total,
 * negative when the swap helps -- and not a total. ``verification_pass`` fills
 * the same field with a recomputed total, so the two producers of a
 * ``SwapCandidate`` mean different things by it and only the comparator is
 * shared; reading one for the other is the easiest mistake to make here.
 *
 * The delta splits into a ``shared`` part that every leaving slot sees and a
 * per-slot ``correction``, which is why one O(n) pass per entering item scores
 * all k slots at once instead of k passes. The cached d1 and d2 are the entire
 * speedup: no step appeals to the triangle inequality, because this library's
 * dissimilarities are free to violate it, and every candidate that is not
 * already a medoid is evaluated.
 *
 * :param data: Condensed distance array from ``StorageBackend::Data()``.
 * :param n: Number of items.
 * :param medoids: Current medoid item indices, one per slot.
 * :param assignments: Cache from ``build_assignments`` for those medoids.
 * :param slot_of: Slot map from ``refresh_slot_map`` for those medoids.
 * :param num_threads: Worker threads; 0 auto-detects hardware concurrency.
 * :param chunk_size: Chunk size for the parallel scan; at least one.
 * :returns: The best candidate, or an invalid candidate when none exists.
 */
inline SwapCandidate best_predicted_swap(
    const double* data, size_t n, const std::vector<size_t>& medoids,
    const std::vector<Assignment>& assignments,
    const std::vector<size_t>& slot_of, size_t num_threads,
    size_t chunk_size) {
    const size_t k = medoids.size();
    const size_t chunk = effective_chunk_size(n, chunk_size);
    const size_t total_chunks = (n + chunk - 1) / chunk;
    std::vector<SwapCandidate> chunk_best(total_chunks);

    ThreadPool pool(num_threads);
    pool.ParallelFor(0, n, chunk, [&](size_t begin, size_t end) {
        std::vector<double> correction(k, 0.0);
        SwapCandidate local;

        for (size_t h = begin; h < end; ++h) {
            if (slot_of[h] != n) {
                continue;
            }

            std::fill(correction.begin(), correction.end(), 0.0);
            double shared = 0.0;
            // A candidate's entire scan runs inside one chunk in ascending item
            // order, so the summation order is fixed by the indexing rather
            // than by the schedule and every delta is bit-identical whatever
            // num_threads and chunk_size are.
            for (size_t j = 0; j < n; ++j) {
                const double d_hj = dense_distance(data, n, h, j);
                const Assignment& entry = assignments[j];
                const double surviving =
                    std::min(d_hj - entry.nearest_distance, 0.0);
                shared += surviving;
                const double leaving =
                    std::min(entry.second_nearest_distance, d_hj) -
                    entry.nearest_distance;
                correction[entry.nearest_slot] += leaving - surviving;
            }

            for (size_t slot = 0; slot < k; ++slot) {
                SwapCandidate candidate;
                candidate.score = shared + correction[slot];
                candidate.entering_item = h;
                candidate.leaving_item = medoids[slot];
                candidate.leaving_slot = slot;
                candidate.valid = true;
                if (is_better(candidate, local)) {
                    local = candidate;
                }
            }
        }

        chunk_best[begin / chunk] = local;
    });

    SwapCandidate winner;
    for (const SwapCandidate& candidate : chunk_best) {
        if (is_better(candidate, winner)) {
            winner = candidate;
        }
    }
    return winner;
}

/**
 * @brief Total cost after replacing one slot's medoid, from the cache alone.
 *
 * Summed in ascending item order over the same summands ``total_cost`` uses, so
 * the comparison against the reported cost is exact. The cache supplies each
 * item's nearest surviving medoid distance: for an item whose medoid survives
 * that is d1, and for an item whose medoid leaves it is d2 by definition.
 *
 * :param data: Condensed distance array from ``StorageBackend::Data()``.
 * :param n: Number of items.
 * :param assignments: Cache from ``build_assignments`` for the current medoids.
 * :param leaving_slot: Slot whose medoid is replaced.
 * :param entering_item: Item index taking that slot.
 * :returns: The objective of the resulting configuration.
 */
inline double recomputed_total(const double* data, size_t n,
                               const std::vector<Assignment>& assignments,
                               size_t leaving_slot, size_t entering_item) {
    double total = 0.0;
    for (size_t j = 0; j < n; ++j) {
        const Assignment& entry = assignments[j];
        const double surviving = entry.nearest_slot == leaving_slot
                                     ? entry.second_nearest_distance
                                     : entry.nearest_distance;
        total += std::min(surviving, dense_distance(data, n, j, entering_item));
    }
    return total;
}

/**
 * @brief Textbook PAM's own scan, on recomputed totals rather than predicted deltas.
 *
 * The fast loop is not allowed to declare convergence: a swap that genuinely
 * lowers the recomputed cost can be predicted at exactly zero, and exiting on
 * that prediction would assert a local optimum nothing had checked.
 *
 * Note that the returned ``score`` is a **total**, not a delta, unlike
 * ``best_predicted_swap``'s.
 *
 * :param data: Condensed distance array from ``StorageBackend::Data()``.
 * :param n: Number of items.
 * :param medoids: Current medoid item indices, one per slot.
 * :param assignments: Cache from ``build_assignments`` for those medoids.
 * :param slot_of: Slot map from ``refresh_slot_map`` for those medoids.
 * :param current_cost: Objective to beat; only strict improvements qualify.
 * :param num_threads: Worker threads; 0 auto-detects hardware concurrency.
 * :param chunk_size: Chunk size for the parallel scan; at least one.
 * :returns: The improving candidate with the smallest total, or an invalid
 *     candidate when the configuration is a verified local optimum.
 */
inline SwapCandidate verification_pass(
    const double* data, size_t n, const std::vector<size_t>& medoids,
    const std::vector<Assignment>& assignments,
    const std::vector<size_t>& slot_of, double current_cost,
    size_t num_threads, size_t chunk_size) {
    const size_t k = medoids.size();
    const size_t chunk = effective_chunk_size(n, chunk_size);
    const size_t total_chunks = (n + chunk - 1) / chunk;
    std::vector<SwapCandidate> chunk_best(total_chunks);

    ThreadPool pool(num_threads);
    pool.ParallelFor(0, n, chunk, [&](size_t begin, size_t end) {
        SwapCandidate local;

        for (size_t h = begin; h < end; ++h) {
            if (slot_of[h] != n) {
                continue;
            }
            for (size_t slot = 0; slot < k; ++slot) {
                const double total =
                    recomputed_total(data, n, assignments, slot, h);
                if (!(total < current_cost)) {
                    continue;
                }
                SwapCandidate candidate;
                candidate.score = total;
                candidate.entering_item = h;
                candidate.leaving_item = medoids[slot];
                candidate.leaving_slot = slot;
                candidate.valid = true;
                if (is_better(candidate, local)) {
                    local = candidate;
                }
            }
        }

        chunk_best[begin / chunk] = local;
    });

    SwapCandidate winner;
    for (const SwapCandidate& candidate : chunk_best) {
        if (is_better(candidate, winner)) {
            winner = candidate;
        }
    }
    return winner;
}

}  // namespace OECluster::detail

#endif  // OECLUSTER_SRC_CLUSTERING_KMEDOIDSSWAPKERNEL_H
