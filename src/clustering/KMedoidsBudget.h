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
 * @brief Comparisons in one verification pass.
 *
 * The pass scores every (candidate, slot) pair and rescans every item for
 * each, so it costs ``k`` sweeps where the predicted-delta scan costs one.
 * Pricing it as a single sweep is what let a budgeted run overrun its limit
 * by a factor of ``k``.
 *
 * :param n: Item count.
 * :param k: Medoid count.
 * :returns: k sweeps, saturating rather than wrapping.
 */
inline size_t verification_comparisons(size_t n, size_t k) {
    const size_t sweep = sweep_comparisons(n);
    const size_t limit = static_cast<size_t>(-1);
    if (sweep != 0 && k > limit / sweep) {
        return limit;
    }
    return k * sweep;
}

/**
 * @brief Comparisons a run must make whatever the swap loop does.
 *
 * The initializer, the one-off medoid column cache fill, and one verification
 * pass: ``(init sweeps + k) * n * (n - 1) + n * k``. Nothing here depends on
 * the iteration count, which is what makes refusing on it honest.
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
            // k, not k - 1: build_initialize() opens with a global-medoid
            // pass to pick the first seed and only then loops k - 1 times.
            init_sweeps = k;
            break;
        case KMedoidsInit::FarthestFirst:
            init_sweeps = 1;
            break;
        case KMedoidsInit::Explicit:
            init_sweeps = 0;
            break;
    }
    // Saturate: these products reach 1e10 at n = 100000 and the caller
    // compares against a budget, so a wrap would read as cheap work. The
    // sweep count is summed before it is multiplied, so it is guarded too --
    // Build's total is 2k, which wraps from k = 2^63 up, and this fires at
    // exactly that boundary.
    const size_t limit = static_cast<size_t>(-1);
    if (init_sweeps > limit - k) {
        return limit;
    }
    const size_t sweeps = init_sweeps + k;
    if (sweep != 0 && sweeps > limit / sweep) {
        return limit;
    }
    const size_t sweeps_total = sweeps * sweep;
    if (n != 0 && k > (limit - sweeps_total) / n) {
        return limit;
    }
    return sweeps_total + n * k;
}

}  // namespace OECluster::detail

#endif  // OECLUSTER_SRC_CLUSTERING_KMEDOIDSBUDGET_H
