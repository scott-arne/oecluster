/**
 * @file ActivityMetrics.h
 * @brief Numeric helpers shared by the SAR-coherence metrics.
 *
 * Header-only and inline, matching ContingencyTable.h. These live outside
 * SARCoherence.cpp's anonymous namespace so the C++ tests can drive the
 * arithmetic directly: a sums-of-squares decomposition, an effect size, and a
 * standard deviation all fail by returning a plausible number rather than by
 * crashing, and the overflow guards can only be exercised at inputs no public
 * entry point should have to be called with twice.
 */

#ifndef OECLUSTER_SRC_CLUSTERING_ACTIVITYMETRICS_H
#define OECLUSTER_SRC_CLUSTERING_ACTIVITYMETRICS_H

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace OECluster::detail {

/// Sets drop[i] for every sample whose activity is missing (NaN). Every other
/// bit is left as found, so this composes with detail::mark_excluded in either
/// order. Throws std::invalid_argument on an infinite value, naming caller.
/// drop must already be sized to activity.size().
inline void mark_missing_activity(const std::vector<double>& activity,
                                  const std::string& caller,
                                  std::vector<bool>& drop) {
    for (std::size_t i = 0; i < activity.size(); ++i) {
        const double value = activity[i];
        if (std::isnan(value)) {
            drop[i] = true;
        } else if (std::isinf(value)) {
            throw std::invalid_argument(
                caller + ": activity[" + std::to_string(i) +
                "] is infinite; use NaN for a missing measurement");
        }
    }
}

/// Original indices of the samples drop does not mark, in input order.
inline std::vector<std::size_t> gather_indices(const std::vector<bool>& drop) {
    std::vector<std::size_t> indices;
    for (std::size_t i = 0; i < drop.size(); ++i) {
        if (!drop[i]) {
            indices.push_back(i);
        }
    }
    return indices;
}

/// Scored samples: original indices and their finite activity values, in
/// input order.
struct ScoredActivity {
    std::vector<std::size_t> indices;
    std::vector<double> values;
};

/// Gathers the samples drop does not mark. drop must be final: every exclusion
/// reason has already been ORed into it.
inline ScoredActivity gather_scored(const std::vector<double>& activity,
                                    const std::vector<bool>& drop) {
    ScoredActivity scored;
    for (std::size_t i = 0; i < activity.size(); ++i) {
        if (!drop[i]) {
            scored.indices.push_back(i);
            scored.values.push_back(activity[i]);
        }
    }
    return scored;
}

/// Population standard deviation, two-pass. NaN for fewer than 2 values.
/// Returns +inf when the mean or the deviation sum overflows, which only
/// happens outside the supported numeric domain; the caller turns that into a
/// throw rather than reporting it.
inline double population_stddev(const std::vector<double>& values) {
    if (values.size() < 2) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    double sum = 0.0;
    for (const double value : values) {
        sum += value;
    }
    const double mean = sum / static_cast<double>(values.size());
    if (std::isinf(mean)) {
        // Report the infinity rather than the exact zero the deviations would
        // suggest: the caller turns this into a refusal, and it can only do
        // that if the unformable mean is visible in the return value.
        return std::numeric_limits<double>::infinity();
    }

    // Refine the mean with a correction pass. A naive sum's error scales with
    // the magnitude of the values, while the deviations it is then subtracted
    // from scale with their spread. When the spread is near the values' rounding
    // floor the naive mean can land a full ulp away from the true mean, driving
    // the spread to zero or reversing the effect sizes. The residuals (value -
    // mean) are on the order of the spread rather than the magnitude, so their
    // sum is accurate where the raw sum is not.
    double correction = 0.0;
    for (const double value : values) {
        correction += value - mean;
    }
    const double refined_mean = mean + correction / static_cast<double>(values.size());

    // Squaring before dividing is what loses the small end. A deviation of
    // 1e-163 squares to zero, so a naive accumulation reports a spread of
    // exactly 0.0 for data that has one -- and 0.0 is not a small answer here,
    // it is the answer that collapses the RMODI band and makes every sample
    // in-band. Dividing each deviation by the largest one first puts the
    // biggest term at exactly 1.0, so the sum can neither underflow nor
    // overflow, and multiplying the scale back in after the square root keeps
    // it out of the squaring entirely.
    double scale = 0.0;
    for (const double value : values) {
        const double magnitude = std::fabs(value - refined_mean);
        if (magnitude > scale) {
            scale = magnitude;
        }
    }
    if (std::isinf(scale)) {
        return std::numeric_limits<double>::infinity();
    }
    if (scale == 0.0) {
        return 0.0;
    }

    double deviation = 0.0;
    for (const double value : values) {
        const double difference = (value - refined_mean) / scale;
        deviation += difference * difference;
    }
    // The scaling rescues the small end, and it would rescue the large end too
    // -- {-1e300, 1e300} has a perfectly representable spread of 1e300 even
    // though its sum of squared deviations does not fit. That is a wider domain
    // than §3.3 declares, and widening it is not this task's decision to make:
    // sums_of_squares refuses the same input, and the two must not disagree
    // about what is in range. deviation is at least 1.0 by construction, so this
    // product is infinite exactly when the unscaled sum would have been.
    if (std::isinf(scale * scale * deviation)) {
        return std::numeric_limits<double>::infinity();
    }
    return scale * std::sqrt(deviation / static_cast<double>(values.size()));
}

/// A one-way sum-of-squares decomposition. total and within are computed;
/// between is derived as clamp(total - within, 0, total), so
/// 0 <= between <= total holds exactly while between + within reproduces total
/// only to within rounding. Deriving in that direction is a correctness
/// decision: computing between instead would let rounding report a perfect
/// effect size on unseparated data, whereas this way rounding reports none.
///
/// The three sums are held **in units of scale squared**, where scale is the
/// largest absolute deviation from the grand mean. Multiply by scale * scale
/// -- or call ss_total, ss_within, ss_between -- for the caller's own units.
/// Holding them scaled is what makes the ratios trustworthy across the whole
/// supported domain: an unscaled squared deviation flushes to zero below about
/// 1e-162, and it does so unevenly, so SS_within can reach zero while SS_total
/// has not. That reports eta squared as 1.0 -- perfect separation -- for data
/// that is not separated at all. Scaled, the largest term is exactly 1.0 by
/// construction, the sums lie in [1, num_scored], and scale cancels out of
/// every ratio, so no ratio ever sees the bottom of the double range.
struct SumsOfSquares {
    double total = 0.0;
    double between = 0.0;
    double within = 0.0;
    double scale = 0.0;
    std::size_t num_scored = 0;
    std::size_t num_groups = 0;
};

/// The decomposition's sums in the caller's own units. These can overflow or
/// underflow where the scaled sums cannot, which is exactly why eta squared
/// and omega squared are not computed through them.
inline double ss_total(const SumsOfSquares& ss) {
    return ss.scale * ss.scale * ss.total;
}
inline double ss_within(const SumsOfSquares& ss) {
    return ss.scale * ss.scale * ss.within;
}
inline double ss_between(const SumsOfSquares& ss) {
    return ss.scale * ss.scale * ss.between;
}

/// Decomposes scored values by group id. Group ids must be dense in
/// [0, num_groups) and parallel to values. Throws std::invalid_argument if the
/// lengths disagree, or if the grand mean, SS_total, or any group mean
/// overflows to an infinity, naming whichever it was.
inline SumsOfSquares sums_of_squares(const std::vector<std::uint32_t>& group_ids,
                                     const std::vector<double>& values,
                                     std::uint32_t num_groups,
                                     const std::string& caller) {
    if (group_ids.size() != values.size()) {
        throw std::invalid_argument(
            caller + ": group ids and values have different lengths");
    }

    SumsOfSquares ss;
    ss.num_scored = values.size();
    ss.num_groups = num_groups;
    if (values.empty()) {
        return ss;
    }

    double grand_sum = 0.0;
    for (const double value : values) {
        grand_sum += value;
    }
    const double grand_mean = grand_sum / static_cast<double>(values.size());
    if (std::isinf(grand_mean)) {
        throw std::invalid_argument(
            caller +
            ": the activity mean overflows to infinity; the supported range is "
            "values whose sum is finite in double precision");
    }

    // Refine the grand mean with a correction pass. A naive sum's error scales
    // with the magnitude of the values, while the deviations it is then
    // subtracted from scale with their spread. When the spread is near the
    // values' rounding floor the naive mean can land a full ulp away from the
    // true mean, driving SS_between to zero or reversing the effect sizes.
    double grand_correction = 0.0;
    for (const double value : values) {
        grand_correction += value - grand_mean;
    }
    const double refined_grand_mean = grand_mean + grand_correction / static_cast<double>(values.size());

    for (const double value : values) {
        const double magnitude = std::fabs(value - refined_grand_mean);
        if (magnitude > ss.scale) {
            ss.scale = magnitude;
        }
    }
    if (ss.scale == 0.0) {
        // Every value equals the grand mean. The sums are all exactly zero and
        // there is no scale to divide by; returning here keeps the loop below
        // from forming 0/0.
        return ss;
    }
    if (std::isinf(ss.scale)) {
        throw std::invalid_argument(
            caller +
            ": SS_total overflows to infinity; the supported range is values "
            "whose squared deviations sum finitely in double precision");
    }

    for (const double value : values) {
        const double difference = (value - refined_grand_mean) / ss.scale;
        ss.total += difference * difference;
    }
    // The guard is on the sum in the caller's units, which is the quantity
    // §3.3's domain is stated over. Squaring the scale first is deliberate:
    // the scaled total is at least 1.0, so the true total is at least the
    // scale squared, and an overflow in that product is therefore an overflow
    // in the truth rather than an artifact of the multiplication order.
    if (std::isinf(ss_total(ss))) {
        throw std::invalid_argument(
            caller +
            ": SS_total overflows to infinity; the supported range is values "
            "whose squared deviations sum finitely in double precision");
    }

    std::vector<double> group_sums(num_groups, 0.0);
    std::vector<std::size_t> group_counts(num_groups, 0);
    for (std::size_t i = 0; i < values.size(); ++i) {
        group_sums[group_ids[i]] += values[i];
        ++group_counts[group_ids[i]];
    }
    std::vector<double> group_means(num_groups, 0.0);
    for (std::uint32_t group = 0; group < num_groups; ++group) {
        if (group_counts[group] == 0) {
            continue;
        }
        group_means[group] =
            group_sums[group] / static_cast<double>(group_counts[group]);
        if (std::isinf(group_means[group])) {
            throw std::invalid_argument(
                caller + ": the mean of cluster " + std::to_string(group) +
                " overflows to infinity; the supported range is values whose "
                "sum is finite in double precision");
        }
    }

    // Refine the group means with a correction pass, for the same reason as
    // the grand mean.
    std::vector<double> group_corrections(num_groups, 0.0);
    for (std::size_t i = 0; i < values.size(); ++i) {
        group_corrections[group_ids[i]] += values[i] - group_means[group_ids[i]];
    }
    for (std::uint32_t group = 0; group < num_groups; ++group) {
        if (group_counts[group] == 0) {
            continue;
        }
        group_means[group] += group_corrections[group] / static_cast<double>(group_counts[group]);
    }

    // Divided by the same scale as SS_total, so the two remain comparable and
    // their difference is meaningful.
    for (std::size_t i = 0; i < values.size(); ++i) {
        const double difference =
            (values[i] - group_means[group_ids[i]]) / ss.scale;
        ss.within += difference * difference;
    }

    // Derive between rather than computing it directly, and clamp to [0, total].
    // Mathematically within <= total because each value is at least as close to
    // its group mean as to the grand mean, but the two sums accumulate in
    // separate loops and round independently. On tightly spaced input where the
    // group means straddle the grand mean, within can finish one ulp above total,
    // leaving total - within strictly negative. The clamp keeps between and the
    // effect sizes from going negative there. See
    // BetweenClampsToZeroWhenWithinExceedsTotal for a reproducing fixture.
    const double between = ss.total - ss.within;
    ss.between = between < 0.0 ? 0.0 : (between > ss.total ? ss.total : between);
    return ss;
}

/// SS_between / SS_total, read off the scaled sums so the ratio is the same
/// whatever units the caller measured in. NaN when fewer than two values were
/// scored or when SS_total is zero, where there is no variance to explain --
/// and scaled, SS_total is zero exactly when every value equals the grand
/// mean, never merely because the deviations were small.
inline double eta_squared(const SumsOfSquares& ss) {
    if (ss.num_scored < 2 || ss.total == 0.0) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return ss.between / ss.total;
}

/// Scale-normalised: (eta2 - (K-1)*m) / (1 + m) with
/// m = (within / total) / df_within, evaluated in that order. This is the
/// textbook expression rearranged so that the first operation cancels the
/// input's scale; the two other orderings lose the result to overflow and to
/// underflow respectively. Do not rearrange. The sums it reads are themselves
/// scaled, which is a second and independent line of defence: the rearranged
/// order protects the division, the scaled sums protect the squaring that
/// produced them.
inline double omega_squared(const SumsOfSquares& ss) {
    if (ss.num_scored < 2 || ss.total == 0.0 ||
        ss.num_scored == ss.num_groups) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    const double eta2 = ss.between / ss.total;
    const double df_within =
        static_cast<double>(ss.num_scored) - static_cast<double>(ss.num_groups);
    const double m = (ss.within / ss.total) / df_within;
    const double groups_minus_one = static_cast<double>(ss.num_groups) - 1.0;
    return (eta2 - groups_minus_one * m) / (1.0 + m);
}

}  // namespace OECluster::detail

#endif  // OECLUSTER_SRC_CLUSTERING_ACTIVITYMETRICS_H
