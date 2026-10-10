/**
 * @file AgglomerativeRowCache.h
 * @brief Complete, average and weighted linkage over an active-slot condensed
 *        matrix with cached row minima.
 */

#ifndef OECLUSTER_SRC_CLUSTERING_AGGLOMERATIVEROWCACHE_H
#define OECLUSTER_SRC_CLUSTERING_AGGLOMERATIVEROWCACHE_H

#include <cstddef>
#include <cstdint>
#include <optional>
#include <string>
#include <vector>

#include "oecluster/StorageBackend.h"
#include "oecluster/clustering/Agglomerative.h"

namespace OECluster::detail {

/**
 * @brief Below this many rows to update a merge step runs on the caller's thread.
 *
 * A step touches one row of the condensed workspace per active slot, so a team
 * only pays for its barrier on a long row list. Set from the plan's measurement.
 */
constexpr size_t ROW_CACHE_SERIAL_CUTOFF = 512;

/**
 * @brief The default participant cap for a merge step.
 *
 * Every step streams one scattered double per active slot and then, for the
 * rows it invalidated, a whole row again, so the pass is bound by memory
 * traffic; and every step ends at a barrier, so each extra participant is one
 * more thread a busy machine can deschedule mid-step.
 */
constexpr size_t ROW_CACHE_DEFAULT_PARTICIPANTS = 8;

/**
 * @brief Participants for a merge step over n items.
 *
 * :param num_threads: Requested count; 0 selects min(hardware, 8).
 * :param n: Item count.
 * :param hardware: What std::thread::hardware_concurrency() reported.
 * :returns: resolve_participants() over that request.
 */
size_t row_cache_participants(size_t num_threads, size_t n, size_t hardware);

/**
 * @brief Whether any merge step can reach the team path.
 *
 * The largest step updates `n - 2` rows, so when that is below the cutoff every
 * step runs serially and a team would start workers that are joined without
 * having updated a row. Starting them costs more than the whole pass at small
 * item counts, which is why this is a question rather than an assumption.
 *
 * :param participants: Resolved participant count, 1 or more.
 * :param n: Item count, 2 or more.
 * :param cutoff: Rows below which a step runs on the calling thread.
 * :returns: True when a team should be constructed.
 */
inline bool row_cache_team_is_useful(size_t participants, size_t n, size_t cutoff) {
    return participants > 1 && n >= 2 && n - 2 >= cutoff;
}

/** @brief The merge tree a linkage pass produces, before any cut is applied. */
struct LinkageTree {
    std::vector<size_t> children_left;   ///< Lower child node id per merge.
    std::vector<size_t> children_right;  ///< Higher child node id per merge.
    std::vector<double> distances;       ///< Merge height per merge.
    std::vector<size_t> cluster_sizes;   ///< Merged cluster size per merge.
};

/**
 * @brief Counters the degenerate-case benchmark reads; unused in production.
 *
 * The rescan is the pass's only super-quadratic term, so its count is the
 * single number that says whether a tie-heavy input has reached the cubic
 * worst case.
 */
struct RowCacheStats {
    uint64_t rescans = 0;          ///< Rows whose cached partner was merged away.
    uint64_t rescanned_slots = 0;  ///< Distances a rescan read, summed.
    uint64_t updates = 0;          ///< Rows a merge step touched, summed.
};

/** @brief What the kernel needs from AgglomerativeOptions, plus test hooks. */
struct RowCacheOptions {
    AgglomerativeLinkageMethod linkage = AgglomerativeLinkageMethod::Average;
    // Merges to record; the pass stops early once it has that many.
    size_t target_merges = 0;
    size_t num_threads = 0;
    // Rows per chunk of the initial copy, clamped to the row count.
    size_t chunk_size = 4096;
    // Steps with fewer rows to update run on the calling thread alone. Unset
    // selects ROW_CACHE_SERIAL_CUTOFF.
    std::optional<size_t> serial_cutoff;
    std::string caller = "agglomerative";
    RowCacheStats* stats = nullptr;
};

/**
 * @brief The merge tree for complete, average or weighted linkage.
 *
 * Holds one condensed matrix over the N slots a live cluster can occupy, plus
 * six O(N) arrays, and never the 2N-1 node table or the all-pairs heap the
 * 5.20.0 path allocated. Each merge picks the active pair that is smallest
 * under the ascending key (distance, lower node id, higher node id), the key
 * 5.20.0's heap popped under, so the tree does not depend on the thread count
 * or the chunk size and is bit-identical to 5.20.0's for a finite non-negative
 * input.
 *
 * :param storage: Complete pairwise distance storage over N >= 2 items.
 * :param options: Linkage, stopping point and parallelism.
 * :returns: The merges in the order they were made; empty for N < 2.
 * :raises std::runtime_error: If an input distance is not finite, naming the
 *     pair of items. Derived cluster distances are not checked: a finite input
 *     can overflow average or weighted linkage to +inf, which 5.20.0 ordered
 *     normally and this path reproduces.
 */
LinkageTree agglomerative_row_cache(const StorageBackend& storage,
                                    const RowCacheOptions& options);

/**
 * @brief agglomerative_cluster()'s complete/average/weighted path, in full.
 *
 * Defined in Agglomerative.cpp beside the single-linkage assembly. Production
 * calls it with both overrides defaulted; the differential test supplies them
 * so that it exercises the one production path rather than a copy of it.
 *
 * :param storage: Complete pairwise distance storage.
 * :param options: Public agglomerative options; validated here.
 * :param serial_cutoff: Overrides ROW_CACHE_SERIAL_CUTOFF when set.
 * :param stats: Receives rescan counters when non-null.
 * :returns: Labels, members and the merge tree.
 */
AgglomerativeResult row_cache_result(
    const StorageBackend& storage,
    const AgglomerativeOptions& options,
    std::optional<size_t> serial_cutoff = std::nullopt,
    RowCacheStats* stats = nullptr);

}  // namespace OECluster::detail

#endif  // OECLUSTER_SRC_CLUSTERING_AGGLOMERATIVEROWCACHE_H
