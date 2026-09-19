/**
 * @file MaxMinKernel.h
 * @brief Deterministic farthest-first (MaxMin) selection over dense distances.
 */

#ifndef OECLUSTER_SRC_CLUSTERING_MAXMINKERNEL_H
#define OECLUSTER_SRC_CLUSTERING_MAXMINKERNEL_H

#include <cstddef>
#include <stdexcept>
#include <vector>

#include "DistanceAccess.h"

namespace OECluster::detail {

/**
 * @brief Deterministic farthest-first (MaxMin) selection over dense distances.
 *
 * Ties resolve to the smaller item index at every step, so the returned
 * sequence depends only on the distances and the seed, never on thread count
 * or iteration order.
 *
 * Already-selected items are excluded from every subsequent pick. That is not
 * an optimization: on a matrix where every distance is equal, an unmasked scan
 * would tie every candidate at zero and the smaller-index rule would return
 * the seed ``count`` times.
 *
 * :param data: Condensed distance array from ``StorageBackend::Data()``.
 * :param n: Number of items.
 * :param count: Number of items to select; must be in [1, n].
 * :param seed: Index of the first selection; must be less than n.
 * :returns: ``count`` item indices in selection order.
 * :raises std::invalid_argument: If count is outside [1, n] or seed is not a
 *     valid item index.
 */
inline std::vector<size_t> maxmin_select_from(const double* data, size_t n,
                                              size_t count, size_t seed) {
    if (count == 0 || count > n) {
        throw std::invalid_argument(
            "MaxMin selection count must be in [1, item count]");
    }
    if (seed >= n) {
        throw std::invalid_argument(
            "MaxMin selection seed is outside the item range");
    }

    std::vector<size_t> selection;
    selection.reserve(count);
    selection.push_back(seed);

    std::vector<bool> selected(n, false);
    selected[seed] = true;

    std::vector<double> nearest(n);
    for (size_t j = 0; j < n; ++j) {
        nearest[j] = dense_distance(data, n, j, seed);
    }

    while (selection.size() < count) {
        // Ascending scan with a strict improvement is the tie rule: the
        // smallest index wins any tie because a later equal never displaces it.
        size_t best = n;
        for (size_t j = 0; j < n; ++j) {
            if (selected[j]) {
                continue;
            }
            if (best == n || nearest[j] > nearest[best]) {
                best = j;
            }
        }

        selected[best] = true;
        selection.push_back(best);
        for (size_t j = 0; j < n; ++j) {
            const double distance = dense_distance(data, n, j, best);
            if (distance < nearest[j]) {
                nearest[j] = distance;
            }
        }
    }

    return selection;
}

}  // namespace OECluster::detail

#endif  // OECLUSTER_SRC_CLUSTERING_MAXMINKERNEL_H
