/**
 * @file KMedoidsBudget.h
 * @brief How much comparison work a streamed k-medoids run costs.
 */

#ifndef OECLUSTER_SRC_CLUSTERING_KMEDOIDSBUDGET_H
#define OECLUSTER_SRC_CLUSTERING_KMEDOIDSBUDGET_H

#include <cstddef>

#include "oecluster/clustering/KMedoids.h"

namespace OECluster::detail {

/// Default ceiling on comparisons for a streamed run: about three hours at a
/// microsecond each. A caller whose comparison is far cheaper or far dearer
/// than that sets max_comparisons rather than living with this number.
constexpr size_t KMEDOIDS_DEFAULT_MAX_COMPARISONS = 10'000'000'000;

/**
 * @brief Comparisons in one pass over every pair.
 *
 * Every ordered pair except the diagonal: a candidate's scan reads its whole
 * row, and nothing is retained between candidates, so each unordered pair is
 * read twice per sweep. The diagonal is free.
 *
 * :param n: Item count.
 * :returns: n * (n - 1), saturating rather than wrapping.
 */
inline size_t sweep_comparisons(size_t n) {
    if (n < 2) {
        return 0;
    }
    // Saturate for the same reason unavoidable_comparisons() does, and it has
    // to happen here too: the product passes 2^64 just above n = 2^32, and a
    // wrap there would answer 2^32 -- cheap enough for any budget to accept.
    const size_t limit = static_cast<size_t>(-1);
    return n > limit / (n - 1) ? limit : n * (n - 1);
}

/**
 * @brief Comparisons a run must make whatever the swap loop does.
 *
 * The initializer, the one-off medoid column cache fill, and the final
 * verification sweep. Nothing here depends on the iteration count, which is
 * what makes refusing on it honest.
 *
 * :param n: Item count.
 * :param k: Medoid count.
 * :param init: Initialization strategy.
 * :returns: The comparison count, saturating rather than wrapping.
 */
inline size_t unavoidable_comparisons(size_t n, size_t k, KMedoidsInit init) {
    const size_t sweep = sweep_comparisons(n);
    size_t init_sweeps = 0;
    switch (init) {
        case KMedoidsInit::Build:
            init_sweeps = k > 0 ? k - 1 : 0;
            break;
        case KMedoidsInit::FarthestFirst:
            init_sweeps = 1;
            break;
        case KMedoidsInit::Explicit:
            init_sweeps = 0;
            break;
    }
    // Saturate: these products reach 1e10 at n = 100000 and the caller
    // compares against a budget, so a wrap would read as cheap work.
    const size_t limit = static_cast<size_t>(-1);
    if (sweep != 0 && init_sweeps + 1 > limit / sweep) {
        return limit;
    }
    const size_t sweeps_total = (init_sweeps + 1) * sweep;
    if (n != 0 && k > (limit - sweeps_total) / n) {
        return limit;
    }
    return sweeps_total + n * k;
}

}  // namespace OECluster::detail

#endif  // OECLUSTER_SRC_CLUSTERING_KMEDOIDSBUDGET_H
