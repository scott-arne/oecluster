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
// Owns the pool and the partition these kernels used to build for themselves.
#include "KMedoidsDistanceSource.h"

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
 * @brief Distances from every item to every medoid, indexed [slot * n + item].
 *
 * The streaming path pays a comparison for each of these, and an accepted swap
 * changes one medoid, so holding them turns a rebuild from ``n * k``
 * comparisons into ``n``. The matrix path builds one too rather than
 * branching: reading it is no dearer than reading the matrix, and one code
 * path is worth more than the bytes.
 *
 * ``local(medoids[slot], item)`` canonicalises internally, so this reads the
 * same value as ``local(item, medoids[slot])``; one spelling is kept here and
 * in refresh_medoid_column so a reader does not have to re-derive that.
 *
 * :param distances: Distance source; see MatrixDistanceRunner.
 * :param medoids: Current medoid item indices, one per slot.
 * :returns: A k * n table of distances to each medoid.
 */
template <typename Runner>
inline std::vector<double> build_medoid_columns(
    Runner& distances, const std::vector<size_t>& medoids) {
    const size_t n = distances.Size();
    const size_t k = medoids.size();
    std::vector<double> columns(k * n);
    distances.Run(n, [&](const auto& local, size_t begin, size_t end) {
        for (size_t item = begin; item < end; ++item) {
            for (size_t slot = 0; slot < k; ++slot) {
                columns[slot * n + item] = local(medoids[slot], item);
            }
        }
    });
    return columns;
}

/**
 * @brief Refill one slot's column after its medoid changed.
 *
 * :param distances: Distance source; see MatrixDistanceRunner.
 * :param columns: Table from build_medoid_columns, one slot overwritten.
 * :param slot: Slot whose medoid changed.
 * :param medoid: Item index now occupying that slot.
 */
template <typename Runner>
inline void refresh_medoid_column(Runner& distances,
                                  std::vector<double>& columns, size_t slot,
                                  size_t medoid) {
    const size_t n = distances.Size();
    distances.Run(n, [&](const auto& local, size_t begin, size_t end) {
        for (size_t item = begin; item < end; ++item) {
            columns[slot * n + item] = local(medoid, item);
        }
    });
}

/**
 * @brief Sort the medoids ascending, carrying their cached columns with them.
 *
 * Permuted rather than refilled: sorting renames slots without changing a
 * single distance, so a streamed caller would be paying ``n * k`` comparisons
 * for values it already holds.
 *
 * :param medoids: Medoid item indices, sorted in place; must be distinct.
 * :param columns: Table from build_medoid_columns, permuted to match.
 */
inline void sort_medoids_with_columns(std::vector<size_t>& medoids,
                                      std::vector<double>& columns) {
    const size_t k = medoids.size();
    const size_t n = k == 0 ? 0 : columns.size() / k;

    std::vector<size_t> order(k);
    for (size_t slot = 0; slot < k; ++slot) {
        order[slot] = slot;
    }
    std::sort(order.begin(), order.end(), [&](size_t a, size_t b) {
        return medoids[a] < medoids[b];
    });

    std::vector<size_t> sorted(k);
    std::vector<double> permuted(columns.size());
    for (size_t slot = 0; slot < k; ++slot) {
        const size_t source = order[slot];
        sorted[slot] = medoids[source];
        std::copy(columns.begin() + static_cast<std::ptrdiff_t>(source * n),
                  columns.begin() + static_cast<std::ptrdiff_t>((source + 1) * n),
                  permuted.begin() + static_cast<std::ptrdiff_t>(slot * n));
    }
    medoids = std::move(sorted);
    columns = std::move(permuted);
}

/**
 * @brief Build the nearest and second-nearest cache for every item.
 *
 * :param distances: Distance source; see MatrixDistanceRunner.
 * :param medoids: Current medoid item indices, one per slot; must be distinct.
 * :param columns: Table from build_medoid_columns for exactly those medoids.
 * :returns: One Assignment per item, indexed by item.
 */
template <typename Runner>
inline std::vector<Assignment> build_assignments(
    Runner& distances, const std::vector<size_t>& medoids,
    const std::vector<double>& columns) {
    const size_t n = distances.Size();
    const size_t k = medoids.size();
    const double infinity = std::numeric_limits<double>::infinity();

    std::vector<size_t> slot_of(n, n);
    refresh_slot_map(slot_of, medoids);

    std::vector<Assignment> assignments(n);
    distances.Run(n, [&](const auto&, size_t begin, size_t end) {
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
                    const double distance = columns[slot * n + j];
                    if (distance < entry.second_nearest_distance) {
                        entry.second_nearest_distance = distance;
                    }
                }
            } else {
                entry.nearest_slot = 0;
                entry.nearest_distance = infinity;
                entry.second_nearest_distance = infinity;
                for (size_t slot = 0; slot < k; ++slot) {
                    const double distance = columns[slot * n + j];
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
 * @brief Build the assignment cache for callers that hold no column cache.
 *
 * Fills the columns and discards them, which is what the swap loop did before
 * it started carrying one. A caller inside the loop holds the table instead:
 * refilling it per call would cost ``n * k`` comparisons an accepted swap only
 * invalidates ``n`` of.
 *
 * :param distances: Distance source; see MatrixDistanceRunner.
 * :param medoids: Current medoid item indices, one per slot; must be distinct.
 * :returns: One Assignment per item, indexed by item.
 */
template <typename Runner>
inline std::vector<Assignment> build_assignments(
    Runner& distances, const std::vector<size_t>& medoids) {
    return build_assignments(distances, medoids,
                             build_medoid_columns(distances, medoids));
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
 * :param distances: Distance source; see MatrixDistanceRunner.
 * :param medoids: Current medoid item indices, one per slot.
 * :param assignments: Cache from ``build_assignments`` for those medoids.
 * :param slot_of: Slot map from ``refresh_slot_map`` for those medoids.
 * :returns: The best candidate, or an invalid candidate when none exists.
 */
template <typename Runner>
inline SwapCandidate best_predicted_swap(
    Runner& distances, const std::vector<size_t>& medoids,
    const std::vector<Assignment>& assignments,
    const std::vector<size_t>& slot_of) {
    const size_t n = distances.Size();
    const size_t k = medoids.size();
    const size_t chunk = distances.ChunkSize(n);
    const size_t total_chunks = (n + chunk - 1) / chunk;
    std::vector<SwapCandidate> chunk_best(total_chunks);

    distances.Run(n, [&](const auto& distance_of, size_t begin, size_t end) {
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
                const double d_hj = distance_of(h, j);
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

        // Index arithmetic assumes ParallelFor hands out uniform chunk-start
        // offsets; guided or split-tail chunking would write past the end.
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
 * Takes an accessor rather than a runner: it is the inner computation of one
 * verification candidate, so it runs on whatever worker already holds the
 * chunk, and giving it a runner would invite a nested pass per candidate.
 *
 * :param distance: Distance accessor; see MatrixDistances.
 * :param n: Number of items.
 * :param assignments: Cache from ``build_assignments`` for the current medoids.
 * :param leaving_slot: Slot whose medoid is replaced.
 * :param entering_item: Item index taking that slot.
 * :returns: The objective of the resulting configuration.
 */
template <typename Accessor>
inline double recomputed_total(const Accessor& distance, size_t n,
                               const std::vector<Assignment>& assignments,
                               size_t leaving_slot, size_t entering_item) {
    double total = 0.0;
    for (size_t j = 0; j < n; ++j) {
        const Assignment& entry = assignments[j];
        const double surviving = entry.nearest_slot == leaving_slot
                                     ? entry.second_nearest_distance
                                     : entry.nearest_distance;
        total += std::min(surviving, distance(j, entering_item));
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
 * :param distances: Distance source; see MatrixDistanceRunner.
 * :param medoids: Current medoid item indices, one per slot.
 * :param assignments: Cache from ``build_assignments`` for those medoids.
 * :param slot_of: Slot map from ``refresh_slot_map`` for those medoids.
 * :param current_cost: Objective to beat; only strict improvements qualify.
 * :returns: The improving candidate with the smallest total, or an invalid
 *     candidate when the configuration is a verified local optimum.
 */
template <typename Runner>
inline SwapCandidate verification_pass(
    Runner& distances, const std::vector<size_t>& medoids,
    const std::vector<Assignment>& assignments,
    const std::vector<size_t>& slot_of, double current_cost) {
    const size_t n = distances.Size();
    const size_t k = medoids.size();
    const size_t chunk = distances.ChunkSize(n);
    const size_t total_chunks = (n + chunk - 1) / chunk;
    std::vector<SwapCandidate> chunk_best(total_chunks);

    distances.Run(n, [&](const auto& distance_of, size_t begin, size_t end) {
        SwapCandidate local;

        for (size_t h = begin; h < end; ++h) {
            if (slot_of[h] != n) {
                continue;
            }
            for (size_t slot = 0; slot < k; ++slot) {
                const double total =
                    recomputed_total(distance_of, n, assignments, slot, h);
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

        // Index arithmetic assumes ParallelFor hands out uniform chunk-start
        // offsets; guided or split-tail chunking would write past the end.
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
