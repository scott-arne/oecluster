/**
 * @file CompareUnordered.h
 * @brief Namespace-scope access to the unordered-pair comparison seam.
 *
 * Private to the library: not installed, not exposed through SWIG.
 *
 * ``PairwiseComparison::CompareUnordered`` is protected so that it cannot
 * become a second public way to compare two items; making it public would
 * also hand it to SWIG, which wraps public members by default. Every caller
 * is a free function, so reaching a protected member means friendship, and
 * without this accessor there would be one friend declaration per caller.
 * There is one, and it names this function.
 *
 * DescriptorBuild.h, StorageOverflow.h and ROCSScoreCombine.h are the
 * precedent for declaring a private helper in a header of its own rather
 * than in the anonymous namespace of the one file that uses it.
 */

#ifndef OECLUSTER_COMPAREUNORDERED_H
#define OECLUSTER_COMPAREUNORDERED_H

#include <cstddef>

#include "oecluster/PairwiseComparison.h"

namespace OECluster::detail {

/**
 * @brief The comparison's value for the unordered pair {i, j}.
 *
 * The one caller-side name for ``PairwiseComparison::CompareUnordered``, whose
 * contract is written on that member. Use it wherever a caller reads d(i, j)
 * as a symmetric distance -- a silhouette, a cluster diameter, a nearest
 * neighbour, a kernel entry -- so that the lazy route and a matrix from
 * ``pdist()`` answer with the same number. ``Compare`` stays the right call
 * where the direction is the question being asked, which is ``cdist()`` and a
 * caller going straight to the comparison.
 *
 * The two indices are passed on in ascending order, so a comparison that does
 * not override the seam reads exactly the ``Compare(min(i, j), max(i, j))``
 * it did before. ``ROCSComparison`` does override it, at the cost of a second
 * overlay search per pair.
 *
 * :param comparison: Comparison to read; typically a per-worker clone.
 * :param i: Index of one item of the pair; either order is accepted.
 * :param j: Index of the other.
 * :returns: The value the comparison reports for the unordered pair.
 */
double compare_unordered(PairwiseComparison& comparison, size_t i, size_t j);

}  // namespace OECluster::detail

#endif  // OECLUSTER_COMPAREUNORDERED_H
