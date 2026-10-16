/**
 * @file ROCSScoreCombine.h
 * @brief Combining the two overlay orientations ROCS scores for one pair.
 *
 * Private to the ROCS comparison: not installed, not exposed through SWIG.
 *
 * Declaring the helper here gives it external linkage, which it would not have
 * in the anonymous namespace ROCSComparison.cpp already opens. That is the
 * trade being made rather than an oversight: the rule below decides what a
 * non-finite orientation contributes, and the inputs that exercise it are ones
 * no overlay fixture can be relied on to produce -- a molecule pair that scores
 * NaN in exactly one direction and a finite value in the other may not exist.
 * Reaching the rule from the test suite is what makes those cases testable at
 * all. DescriptorBuild.h and StorageOverflow.h declare their private helpers
 * the same way, and one shared pattern is worth more than one fewer symbol.
 */

#ifndef OECLUSTER_COMPARISONS_ROCSSCORECOMBINE_H
#define OECLUSTER_COMPARISONS_ROCSSCORECOMBINE_H

namespace OECluster {

/**
 * @brief The better of two raw overlay scores for the same unordered pair.
 *
 * ``pdist`` stores one triangle, so ROCS scores a pair in both orientations and
 * keeps the better overlay. Every ROCS raw score is a Tanimoto that rises with
 * overlay quality, whatever the configured score type or direction flag, so
 * "better" is "larger" here and the caller converts once afterwards. Taking the
 * maximum of the converted values instead would need a branch on the direction
 * flag and would pick the worse overlay on the distance path.
 *
 * A non-finite operand collapses to one canonical quiet NaN rather than being
 * returned as it stands. ``std::max`` compares false against NaN and so returns
 * whichever operand it was handed first, which would leave a ``{NaN, -inf}``
 * pair scoring differently depending on the argument order -- the same
 * order-dependence this function exists to remove, one level down. The test is
 * ``!std::isfinite``, not ``std::isnan``: an infinity is not a score either,
 * and a NaN-only test cannot tell the two rules apart.
 *
 * :param a: Raw score for one orientation.
 * :param b: Raw score for the other orientation.
 * :returns: The larger of the two when both are finite; a quiet NaN otherwise.
 */
double better_raw_score(double a, double b);

}  // namespace OECluster

#endif  // OECLUSTER_COMPARISONS_ROCSSCORECOMBINE_H
