/**
 * @file InternalIndices.h
 * @brief Numeric helpers for the internal cluster-validity indices.
 *
 * Header-only and inline, matching DistanceAccess.h. These live outside
 * ClusterReport.cpp's anonymous namespace so the C++ tests can drive them
 * directly: the arithmetic they decide -- a wrapping couple counter, a
 * cancelling variance, a sign taken after a width was dropped -- fails by
 * returning a plausible number rather than by crashing.
 */

#ifndef OECLUSTER_SRC_CLUSTERING_INTERNALINDICES_H
#define OECLUSTER_SRC_CLUSTERING_INTERNALINDICES_H

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "oecluster/StorageBackend.h"

namespace OECluster::detail {

/**
 * @brief Running count, mean and sum of squared deviations (Welford).
 *
 * The point-biserial index needs a standard deviation over every pairwise
 * distance. Distances live in [0, 1] and are frequently near-constant on
 * fingerprint data, which is exactly where sqrt(sum_sq/P - mean^2) cancels: at
 * P on the order of 1e9 the two terms agree to more significant digits than a
 * double carries, and the difference can come out negative. Welford costs one
 * extra multiply per value and is unconditionally stable.
 */
struct DistanceMoments {
    size_t count = 0;
    double mean = 0.0;
    double m2 = 0.0;

    void Add(const double value) {
        ++count;
        const double delta = value - mean;
        mean += delta / static_cast<double>(count);
        m2 += delta * (value - mean);
    }
};

/// Combines two independently accumulated streams (Chan's parallel variance).
inline DistanceMoments merge_moments(const DistanceMoments& a, const DistanceMoments& b) {
    if (a.count == 0) {
        return b;
    }
    if (b.count == 0) {
        return a;
    }
    DistanceMoments merged;
    merged.count = a.count + b.count;
    const double count_a = static_cast<double>(a.count);
    const double count_b = static_cast<double>(b.count);
    const double total = static_cast<double>(merged.count);
    const double delta = b.mean - a.mean;
    merged.mean = a.mean + delta * count_b / total;
    merged.m2 = a.m2 + b.m2 + delta * delta * count_a * count_b / total;
    return merged;
}

/// Population standard deviation; 0.0 for an empty or single-value stream.
inline double population_stddev(const DistanceMoments& moments) {
    if (moments.count == 0) {
        return 0.0;
    }
    return std::sqrt(moments.m2 / static_cast<double>(moments.count));
}

/**
 * @brief Adds a Baker-Hubert couple count, refusing rather than wrapping.
 *
 * A pre-walk guard on P_w * P_b is the obvious alternative and is wrong: that
 * product is an upper bound ties never reach, so a large all-tied input whose
 * Gamma is a defined NaN would be refused. This fires only when the counter
 * that actually exists would wrap.
 *
 * length_error rather than a bespoke type: a caller at this scale is already in
 * memory trouble, and the SWIG layer maps length_error to Python MemoryError,
 * which tells them the truthful thing -- the dataset is too large for this
 * stage, not that their clustering is bad.
 *
 * :raises std::length_error: if the addition would exceed 2^64 - 1.
 */
inline void add_couples(unsigned long long& counter, const unsigned long long increment) {
    if (counter > std::numeric_limits<unsigned long long>::max() - increment) {
        throw std::length_error(
            "cluster_report: Baker-Hubert Gamma couple count would exceed the "
            "64-bit counter (" +
            std::to_string(counter) + " + " + std::to_string(increment) + ")");
    }
    counter += increment;
}

/// The two pair-rank indices, both NaN when undefined.
struct PairRankIndices {
    double c_index = std::numeric_limits<double>::quiet_NaN();
    double baker_hubert_gamma = std::numeric_limits<double>::quiet_NaN();
};

/**
 * @brief C-index and Baker-Hubert Gamma from within- and between-pair distances.
 *
 * Both arrays are taken by value and sorted in place. Keeping them separate
 * rather than tagging each distance drops the per-element tag: the allocation
 * is P doubles rather than P 16-byte structs, and both indices come off the
 * same two sorted arrays.
 *
 * Both arrays must contain only finite values -- read them through
 * checked_distance. A NaN element does not merely give a meaningless index: the
 * sorts below have no strict weak ordering over such a range, which is already
 * undefined behaviour, and even past them the run scans advance on `== value`,
 * which is false for a NaN against itself, so the walk would not terminate.
 *
 * :param within: every pairwise distance inside a cluster.
 * :param between: every pairwise distance across two distinct clusters.
 * :returns: c_index NaN when S_max == S_min, there are no within-pairs, or
 *     there are no between-pairs; baker_hubert_gamma NaN when no couple is
 *     concordant or discordant.
 * :raises std::length_error: if a couple counter would wrap (see add_couples).
 */
inline PairRankIndices pair_rank_indices(
    std::vector<double> within,
    std::vector<double> between) {
    std::sort(within.begin(), within.end());
    std::sort(between.begin(), between.end());
    const size_t within_count = within.size();
    const size_t between_count = between.size();

    PairRankIndices result;

    // With no within-pairs or no between-pairs the index is undefined by
    // definition. Naming both here tells the reader which input to fix rather
    // than making them derive it from a zero denominator two screens down.
    if (within_count > 0 && between_count > 0) {
        // The two merge walks pick the within_count smallest and the
        // within_count largest of the pooled distances, but only their index
        // splits are needed: what the index is defined from is S_w - S_min and
        // S_max - S_w, and each of those is a sum of differences between
        // elements the walk left behind and elements it took in their place.
        //
        // Accumulating those differences rather than three separate sums is
        // what makes this correct on real data. Distances are frequently
        // near-constant, so S_w, S_min and S_max agree to nearly every digit a
        // double carries; differencing them afterwards cancels away the very
        // quantity the index measures, and when the cancellation is total the
        // function used to report a defined index as undefined.
        size_t low_within = 0;
        size_t low_between = 0;
        for (size_t taken = 0; taken < within_count; ++taken) {
            if (low_within < within_count &&
                (low_between >= between_count ||
                 within[low_within] <= between[low_between])) {
                ++low_within;
            } else {
                ++low_between;
            }
        }
        // Every element the walk took is <= every element it left, so pairing
        // the within-elements it left against the between-elements it took in
        // ascending order gives one non-negative term per swap.
        double lower_gap = 0.0;
        for (size_t k = 0; k < low_between; ++k) {
            lower_gap += within[low_within + k] - between[k];
        }

        size_t high_within = within_count;
        size_t high_between = between_count;
        for (size_t taken = 0; taken < within_count; ++taken) {
            if (high_within > 0 &&
                (high_between == 0 ||
                 within[high_within - 1] >= between[high_between - 1])) {
                --high_within;
            } else {
                --high_between;
            }
        }
        const size_t taken_between = between_count - high_between;
        double upper_gap = 0.0;
        for (size_t k = 0; k < taken_between; ++k) {
            upper_gap +=
                between[between_count - 1 - k] - within[taken_between - 1 - k];
        }

        // Both gaps are sums of non-negative terms, so the sum is zero only
        // when every term is, which is exactly S_max == S_min -- no longer a
        // floating-point approximation of that question.
        const double denominator = lower_gap + upper_gap;
        if (denominator != 0.0) {
            result.c_index = lower_gap / denominator;
        }
    }

    unsigned long long concordant = 0;  // d_within < d_between
    unsigned long long discordant = 0;  // d_between < d_within
    size_t next_within = 0;
    size_t next_between = 0;
    unsigned long long within_seen = 0;
    unsigned long long between_seen = 0;
    while (next_within < within_count || next_between < between_count) {
        const double value =
            (next_within < within_count &&
             (next_between >= between_count ||
              within[next_within] <= between[next_between]))
                ? within[next_within]
                : between[next_between];

        // The whole run of `value` is consumed from both arrays at once, so
        // within/between couples inside the run are ties and score neither
        // concordant nor discordant.
        size_t run_within = 0;
        while (next_within + run_within < within_count &&
               within[next_within + run_within] == value) {
            ++run_within;
        }
        size_t run_between = 0;
        while (next_between + run_between < between_count &&
               between[next_between + run_between] == value) {
            ++run_between;
        }

        // Accumulated one element at a time rather than as a product, so the
        // guard sees every value the counter actually takes -- a multiply could
        // wrap before add_couples ever looked at it.
        for (size_t i = 0; i < run_between; ++i) {
            add_couples(concordant, within_seen);
        }
        for (size_t i = 0; i < run_within; ++i) {
            add_couples(discordant, between_seen);
        }

        within_seen += static_cast<unsigned long long>(run_within);
        between_seen += static_cast<unsigned long long>(run_between);
        next_within += run_within;
        next_between += run_between;
    }

    const double denominator =
        static_cast<double>(concordant) + static_cast<double>(discordant);
    if (denominator > 0.0) {
        // The sign is taken before the width is dropped. Gamma is negative
        // whenever discordant couples dominate, which is an ordinary outcome for
        // a bad clustering and the case the index exists to report;
        // concordant - discordant on unsigned long long wraps to a value near
        // 2^64 and the ratio comes out near +1.
        const double numerator = concordant >= discordant
            ? static_cast<double>(concordant - discordant)
            : -static_cast<double>(discordant - concordant);
        result.baker_hubert_gamma = numerator / denominator;
    }

    return result;
}

/**
 * @brief Refuses a non-finite distance with cluster_report's message.
 *
 * Shared by the storage read below and the comparison source, so both paths
 * name a bad pair identically.
 *
 * :raises std::invalid_argument: if the distance is NaN or infinite.
 */
inline double finite_report_distance(const double distance, const size_t i,
                                     const size_t j) {
    if (!std::isfinite(distance)) {
        throw std::invalid_argument(
            "cluster_report: distance between samples " + std::to_string(i) +
            " and " + std::to_string(j) + " is not finite");
    }
    return distance;
}

/**
 * @brief Reads one distance, refusing a non-finite value.
 *
 * A NaN distance is not merely an undefined metric: std::sort over a range
 * containing one has no strict weak ordering, so both the pair-rank sorts and
 * the pre-existing median_distance(intra_pairs) call are undefined behaviour on
 * such input. Refusing replaces undefined behaviour with a named error rather
 * than removing a defined result.
 *
 * :raises std::invalid_argument: if the distance is NaN or infinite.
 */
inline double checked_distance(
    const StorageBackend& storage,
    const size_t i,
    const size_t j) {
    return finite_report_distance(storage.Get(i, j), i, j);
}

}  // namespace OECluster::detail

#endif  // OECLUSTER_SRC_CLUSTERING_INTERNALINDICES_H
