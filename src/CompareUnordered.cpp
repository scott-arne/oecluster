/**
 * @file CompareUnordered.cpp
 * @brief The one call through the protected unordered-pair seam.
 */

#include "CompareUnordered.h"

#include <algorithm>

namespace OECluster::detail {

// Declared in CompareUnordered.h, which is where its contract is written.
// PairwiseComparison names this function, and only this function, a friend.
//
// Ascending order, and the seam's contract depends on it: the base
// CompareUnordered forwards to Compare without reordering, so a directional
// comparison that does not override the seam would answer a caller's (j, i)
// with the other triangle's value. Canonicalizing here is what makes the base
// default a function of the unordered pair, and it is also the cell pdist
// stores, since condensed_to_pair yields i < j.
double compare_unordered(PairwiseComparison& comparison, size_t i, size_t j) {
    return comparison.CompareUnordered(std::min(i, j), std::max(i, j));
}

}  // namespace OECluster::detail
