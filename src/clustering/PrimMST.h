/**
 * @file PrimMST.h
 * @brief Dense Prim's algorithm over a matrix or a comparison, shared by HDBSCAN
 *        and single-linkage agglomerative clustering.
 */

#ifndef OECLUSTER_SRC_CLUSTERING_PRIMMST_H
#define OECLUSTER_SRC_CLUSTERING_PRIMMST_H

#include <cstddef>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

#include "HDBSCANLinkage.h"
#include "oecluster/PairwiseComparison.h"
#include "oecluster/StorageBackend.h"

namespace OECluster::detail {

/**
 * @brief Below this many remaining candidates a matrix step runs serially.
 *
 * A matrix read costs nanoseconds, so a team step only pays for its
 * synchronization on a long candidate list. Set from the plan's measurement.
 */
constexpr size_t PRIM_SERIAL_CUTOFF_MATRIX = 8192;

/**
 * @brief The default participant cap for the Prim pass.
 *
 * Every step re-reads every remaining candidate, so the pass is bound by memory
 * traffic and measured no faster beyond 8 threads; and every step ends at a
 * barrier, so each extra participant is one more thread a busy machine can
 * deschedule mid-step. An explicit num_threads is honored as given.
 */
constexpr size_t PRIM_DEFAULT_PARTICIPANTS = 8;

/**
 * @brief Participants for the Prim pass over n items.
 *
 * :param num_threads: Requested count; 0 selects min(hardware, 8).
 * :param hardware: What std::thread::hardware_concurrency() reported.
 * :returns: resolve_participants() over that request.
 */
size_t prim_participants(size_t num_threads, size_t n, size_t hardware);

/** @brief How an edge's weight is derived from a distance. */
struct PrimWeights {
    // Empty for single linkage, whose weight is the distance itself. Otherwise
    // HDBSCAN's core distances, and the weight is max(core_i, core_j, d / alpha).
    std::vector<double> core;
    double alpha = 1.0;
    // Skip a pair whose weight cannot lower its candidate's reach. Valid only
    // when a core pass has already read, and so checked, every pair.
    bool prune = false;
};

struct PrimOptions {
    size_t num_threads = 0;
    // 0 selects min(hardware, PRIM_DEFAULT_PARTICIPANTS) participants.
    // Steps with fewer remaining candidates than serial_cutoff run on the
    // calling thread alone. Unset selects the provider's default:
    // PRIM_SERIAL_CUTOFF_MATRIX for a matrix, twice the participant count for
    // a comparison.
    std::optional<size_t> serial_cutoff;
    std::string caller = "prim";
};

/**
 * @brief The error for a distance whose quotient by alpha is not finite.
 */
std::invalid_argument alpha_overflow_error(const std::string& caller, double alpha);

/**
 * @brief Minimum spanning tree from node 0, edges in the order nodes join.
 *
 * Each step picks the remaining candidate with the smallest (reach, index), the
 * choice the sequential scan of 5.19.0 made, so the edges do not depend on the
 * thread count or the cutoff. Every distance read must be finite, and also
 * non-negative when weights.core is set; a zero is read as +0.0.
 *
 * :raises std::runtime_error: For a distance outside that domain, naming the pair.
 * :raises std::invalid_argument: If a distance divided by alpha is not finite.
 */
std::vector<HDBSCANMSTEdge> prim_mst(const StorageBackend& storage,
                                     const PrimWeights& weights,
                                     const PrimOptions& options);

/**
 * @brief The same tree, reading Compare(min(i, j), max(i, j)) on one clone per
 *        participant, built serially before any thread starts.
 */
std::vector<HDBSCANMSTEdge> prim_mst(PairwiseComparison& comparison,
                                     const PrimWeights& weights,
                                     const PrimOptions& options);

}  // namespace OECluster::detail

#endif  // OECLUSTER_SRC_CLUSTERING_PRIMMST_H
