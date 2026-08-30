/**
 * @file GateFacts.h
 * @brief Machine-checkable facts a comparison reports about its own metric behavior.
 */

#ifndef OECLUSTER_GATEFACTS_H
#define OECLUSTER_GATEFACTS_H

namespace OECluster {

/**
 * @brief Tri-state answer to a metric capability question.
 *
 * ``Unknown`` means the comparison cannot prove the property either way; the
 * gate treats it as permissive so that new comparisons are not silently
 * refused before they declare themselves.
 */
enum class Capability { Unknown, No, Yes };

/**
 * @brief How completely the scored values cover the requested pairs.
 *
 * This is a conservative declaration of what the configured policy admits,
 * not a measurement of the data. ``Complete`` means every requested pair
 * was scored on the full data. ``NaNPresent`` is reported when the policy
 * permits non-finite distances, even if this particular input produced none;
 * it also escalates to ``NaNPresent`` whenever a non-finite distance is
 * actually produced, whatever the policy -- including under ``ignore``, where
 * the declared stamp would otherwise be the weaker ``SubsetScored``.
 * ``SubsetScored`` means values were computed from a per-pair subset of the
 * available dimensions, which makes them mutually incomparable.
 */
enum class DataIntegrity { Complete, NaNPresent, SubsetScored };

/**
 * @brief The four facts the capability gate reads off a comparison.
 *
 * Defaults are deliberately permissive: a comparison that does not override
 * ``PairwiseComparison::Facts`` reports ``Unknown`` capabilities and
 * ``Complete`` integrity, so the gate lets its results through.
 *
 * ``is_distance`` and ``zero_self`` answer different questions and must not be
 * conflated. ``is_distance`` is about orientation -- whether a small value
 * means "close" -- and it is what the gate refuses a similarity matrix on.
 * ``zero_self`` is about the number the comparison actually returns on the
 * diagonal. A comparison can be a similarity whose self-value happens to be
 * zero, and reporting that honestly is the point of keeping the two apart.
 */
struct GateFacts {
    Capability is_distance = Capability::Unknown;  ///< Do small values mean "close"?
    Capability zero_self = Capability::Unknown;    ///< Is d(x, x) == 0?
    Capability triangle = Capability::Unknown;     ///< Does d obey the triangle inequality?
    DataIntegrity data_integrity = DataIntegrity::Complete;  ///< Coverage of the scored values.
};

}  // namespace OECluster

#endif  // OECLUSTER_GATEFACTS_H
