/**
 * @file Consensus.h
 * @brief Co-association consensus over an ensemble of partitions.
 */
#ifndef OECLUSTER_CONSENSUS_H
#define OECLUSTER_CONSENSUS_H

#include <cstddef>
#include <vector>

#include "oecluster/StorageBackend.h"

namespace OECluster {

/**
 * @brief Shared tuning for the consensus kernels.
 */
struct ConsensusOptions {
    size_t num_threads = 0;   ///< Worker threads; 0 auto-detects hardware concurrency.
    size_t chunk_size = 4096; ///< Rows or clusters per work unit.
};

/**
 * @brief What `coassociation_distances` observed while it ran.
 */
struct ConsensusMatrixSummary {
    size_t num_partitions = 0;   ///< Members accumulated.
    size_t unobserved_pairs = 0; ///< Pairs no member observed together.
};

/**
 * @brief Per-item and per-cluster mean co-association.
 */
struct ConsensusStrength {
    std::vector<double> item_consensus;    ///< One per item; NaN for a singleton or noise.
    std::vector<double> cluster_consensus; ///< One per distinct non-negative label, ascending.
};

/**
 * @brief Fill `destination` with the co-association distance of every pair.
 *
 * The distance is ``1 - co(i, j) / obs(i, j)``, where ``co`` counts the
 * members that placed both items in one cluster and ``obs`` the members that
 * observed both. A pair no member observed together has no evidence either
 * way and is given distance 1.0, counted in the returned summary.
 *
 * Members arrive concatenated: ``offsets`` holds ``R + 1`` entries and member
 * ``r`` owns ``positions[offsets[r] .. offsets[r + 1])`` with its labels at
 * the same slots. A negative label is noise and joins no cluster.
 *
 * The destination is caller-allocated with ``NumSamples() == num_items``, the
 * pdist arrangement, and is zeroed before accumulation: `MMapStorage` reuses
 * an existing file of the right size without clearing it, so a second run
 * over one path would otherwise accumulate onto the first run's values.
 *
 * :param num_items: Number of items N.
 * :param offsets: Member boundaries into `positions` and `labels`.
 * :param positions: Concatenated item positions, distinct within a member.
 * :param labels: Concatenated labels, one per position.
 * :param destination: Storage for the distances, sized to `num_items`.
 * :param options: Threads and chunk size.
 * :returns: The member count and the unobserved-pair count.
 * :raises std::invalid_argument: For fewer than two items, a destination of
 *     another size or without contiguous data, a malformed offset vector, an
 *     empty member, an out-of-range or repeated position within one member,
 *     or a zero chunk size.
 */
ConsensusMatrixSummary coassociation_distances(
    size_t num_items,
    const std::vector<size_t>& offsets,
    const std::vector<size_t>& positions,
    const std::vector<int>& labels,
    StorageBackend& destination,
    const ConsensusOptions& options = ConsensusOptions());

/**
 * @brief Union every pair whose co-association is at least `threshold`.
 *
 * The comparison runs in the distance domain against
 * ``cutoff = 1 - threshold``: both sides lose the same rounding that way,
 * so a pair whose support equals the threshold merges. Components are
 * numbered by their smallest member position, ascending, so the labels
 * depend only on the matrix and the threshold.
 *
 * :param matrix: A consensus matrix with contiguous data.
 * :param threshold: Co-association fraction in [0, 1].
 * :param options: Unused by this sequential pass; accepted for symmetry.
 * :returns: One label per item, `0..K-1`.
 * :raises std::invalid_argument: For a matrix without contiguous data or with
 *     fewer than two items, or a threshold outside [0, 1] or not finite.
 */
std::vector<int> consensus_components(
    const StorageBackend& matrix,
    double threshold,
    const ConsensusOptions& options = ConsensusOptions());

/**
 * @brief Mean co-association within each cluster, per item and per cluster.
 *
 * :param matrix: A consensus matrix with contiguous data.
 * :param labels: One label per item; negative labels are noise.
 * :param options: Threads and chunk size.
 * :returns: `item_consensus` of length N and one `cluster_consensus` per
 *     distinct non-negative label in ascending order; NaN where a cluster
 *     has one member and for an item that is noise.
 * :raises std::invalid_argument: For a matrix without contiguous data, a
 *     label count that disagrees with the matrix, or a zero chunk size.
 */
ConsensusStrength consensus_strength(
    const StorageBackend& matrix,
    const std::vector<int>& labels,
    const ConsensusOptions& options = ConsensusOptions());

}  // namespace OECluster

#endif  // OECLUSTER_CONSENSUS_H
