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
 * :raises ComparisonError: When the indices are identical or either is not less than ``n``.
 */
inline size_t pair_to_condensed(size_t i, size_t j, size_t n) {
    if (i == j) {
        throw ComparisonError("Condensed pair indices must be distinct");
    }
    const size_t lo = i < j ? i : j;
    const size_t hi = i < j ? j : i;
    if (hi >= n) {
        throw ComparisonError("Condensed pair index is out of range");
    }
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
    if (n < 2 || index >= n * (n - 1) / 2) {
        throw ComparisonError("Condensed pair index is out of range");
    }

    // The number of pairs above row r rises monotonically with r, so the row
    // owning an offset is the largest r whose row_start does not exceed it.
    // Integer-only on purpose: a floating-point quadratic solve is O(1) but
    // its rounding is a correctness risk at large n that no test can cheaply
    // pin. r * (2n - r - 1) is always even, so the division is exact.
    auto row_start = [n](size_t row) { return row * (2 * n - row - 1) / 2; };

    size_t low = 0;
    size_t high = n - 2;
    while (low < high) {
        const size_t mid = low + (high - low + 1) / 2;
        if (row_start(mid) <= index) {
            low = mid;
        } else {
            high = mid - 1;
        }
    }

    i = low;
    j = low + 1 + (index - row_start(low));
}

}  // namespace OECluster

#endif  // OECLUSTER_CONDENSEDINDEX_H
