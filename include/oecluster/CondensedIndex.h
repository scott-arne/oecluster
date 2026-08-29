/**
 * @file CondensedIndex.h
 * @brief Index arithmetic for the condensed upper-triangle pairwise layout.
 */

#ifndef OECLUSTER_CONDENSEDINDEX_H
#define OECLUSTER_CONDENSEDINDEX_H

#include <cstddef>

#include "oecluster/Error.h"

namespace OECluster {

/**
 * @brief Map a row-major upper-triangle pair to its condensed offset.
 *
 * The two indices may be given in either order.
 *
 * :param i: First item index.
 * :param j: Second item index, distinct from ``i``.
 * :param n: Number of items.
 * :returns: The offset of the pair within the condensed vector.
 */
inline size_t pair_to_condensed(size_t i, size_t j, size_t n) {
    const size_t lo = i < j ? i : j;
    const size_t hi = i < j ? j : i;
    return n * lo + hi - ((lo + 2) * (lo + 1)) / 2;
}

/**
 * @brief Recover the pair of item indices addressed by a condensed offset.
 *
 * :param index: Offset within the condensed vector.
 * :param n: Number of items.
 * :param i: Receives the lower item index.
 * :param j: Receives the higher item index.
 * :raises ComparisonError: When ``index`` is not a valid condensed offset for ``n``.
 */
inline void condensed_to_pair(size_t index, size_t n, size_t& i, size_t& j) {
    size_t remaining = index;
    for (size_t row = 0; row + 1 < n; ++row) {
        const size_t row_length = n - row - 1;
        if (remaining < row_length) {
            i = row;
            j = row + 1 + remaining;
            return;
        }
        remaining -= row_length;
    }
    throw ComparisonError("Condensed pair index is out of range");
}

}  // namespace OECluster

#endif  // OECLUSTER_CONDENSEDINDEX_H
