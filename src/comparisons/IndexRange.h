/**
 * @file IndexRange.h
 * @brief Index-range refusal shared by the concrete comparison classes.
 *
 * Private to the comparison implementations: not installed, not exposed through SWIG.
 */

#ifndef OECLUSTER_COMPARISONS_INDEXRANGE_H
#define OECLUSTER_COMPARISONS_INDEXRANGE_H

#include <cstddef>
#include <string>

#include "oecluster/Error.h"

namespace OECluster::detail {

/**
 * @brief Build and throw the out-of-range diagnostic for a bad Compare index.
 *
 * Deliberately out of line, for the reason ``report_index_out_of_range`` in
 * StorageBackend.cpp records: the message temporaries would otherwise force
 * exception cleanup paths into every Compare.
 *
 * :param comparison: Comparison class name, used in the message.
 * :param index: The offending index.
 * :param n: Number of items the comparison holds.
 * :raises ComparisonError: Always.
 */
[[noreturn]] inline void report_compare_index_out_of_range(const char* comparison, size_t index,
                                                           size_t n) {
    throw ComparisonError(std::string(comparison) + " index " + std::to_string(index) +
                          " is outside the comparison range of " + std::to_string(n) + " items");
}

/**
 * @brief Refuse a Compare whose indices fall outside the held item range.
 *
 * ``Compare`` is the only one of these classes' index-taking methods SWIG
 * exports, so an out-of-range index is reachable from Python. Every
 * implementation indexes its container without checking first, and the
 * unchecked read is worse than a crash near the end of the container: a
 * two-molecule fingerprint comparison returned 0.0 -- which reads as "these two
 * are identical" -- for ``Compare(1000000, 1000001)``, and only took the
 * process down with SIGSEGV once the index was far enough out to leave the
 * mapping. The bulk paths need no such guard: ``pdist`` and ``cdist`` derive
 * every index they pass from ``Size()``.
 *
 * The message follows ``check_index_range`` in StorageBackend.cpp, which
 * refuses the same class of mistake one layer down.
 *
 * :param comparison: Comparison class name, used in the message.
 * :param i: Index of first item.
 * :param j: Index of second item.
 * :param n: Number of items the comparison holds.
 * :raises ComparisonError: If either index is at or beyond n.
 */
inline void check_compare_index_range(const char* comparison, size_t i, size_t j, size_t n) {
    if (i >= n || j >= n) {
        report_compare_index_out_of_range(comparison, i >= n ? i : j, n);
    }
}

}  // namespace OECluster::detail

#endif  // OECLUSTER_COMPARISONS_INDEXRANGE_H
