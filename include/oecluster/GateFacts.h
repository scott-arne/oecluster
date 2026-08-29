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
 * ``Complete`` means every requested pair was scored on the full data.
 * ``NaNPresent`` means at least one value is not a number. ``SubsetScored``
 * means values were computed from a per-pair subset of the available
 * dimensions, which makes them mutually incomparable.
 */
enum class DataIntegrity { Complete, NaNPresent, SubsetScored };

/**
 * @brief The three facts the capability gate reads off a comparison.
 *
 * Defaults are deliberately permissive: a comparison that does not override
 * ``PairwiseComparison::Facts`` reports ``Unknown`` capabilities and
 * ``Complete`` integrity, so the gate lets its results through.
 */
struct GateFacts {
    Capability zero_self = Capability::Unknown;  ///< Is d(x, x) == 0?
    Capability triangle = Capability::Unknown;   ///< Does d obey the triangle inequality?
    DataIntegrity data_integrity = DataIntegrity::Complete;  ///< Coverage of the scored values.
};

}  // namespace OECluster

#endif  // OECLUSTER_GATEFACTS_H
