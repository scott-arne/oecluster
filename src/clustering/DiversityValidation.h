/**
 * @file DiversityValidation.h
 * @brief Validators shared by the set-level diversity entry points.
 *
 * Shared by the selection entry points (DiversitySelection.cpp) and the set
 * diversity scores (SetDiversity.cpp), so that one comparison is refused, and
 * one thread request capped, the same way everywhere.
 */

#ifndef OECLUSTER_SRC_CLUSTERING_DIVERSITYVALIDATION_H
#define OECLUSTER_SRC_CLUSTERING_DIVERSITYVALIDATION_H

#include <algorithm>
#include <cstddef>
#include <stdexcept>
#include <string>

#include "oecluster/Error.h"
#include "oecluster/GateFacts.h"
#include "oecluster/PairwiseComparison.h"

namespace OECluster::detail {

// The pre-scoring half of the Python gate's require_comparable: every refusal
// here is decided by the comparison's declared facts, before any pair is
// scored, so a count-limited call cannot succeed merely because it stopped
// short of a pair the comparison already said was bad. Unknown is accepted.
inline void validate_comparison_facts(const PairwiseComparison& comparison,
                                      const std::string& caller) {
    const GateFacts facts = comparison.Facts();
    if (facts.is_distance == Capability::No) {
        throw ComparisonError(
            caller + " requires distances, but the comparison reports "
            "similarities");
    }
    if (facts.zero_self == Capability::No) {
        throw ComparisonError(
            caller + " requires a zero self-distance, but the comparison "
            "reports that d(x, x) is not zero");
    }
    if (facts.data_integrity == DataIntegrity::NaNPresent) {
        throw ComparisonError(
            caller + " cannot rank distances the comparison declares may be "
            "non-finite (missing='propagate')");
    }
    if (facts.data_integrity == DataIntegrity::SubsetScored) {
        throw ComparisonError(
            caller + " cannot rank distances scored on per-pair feature "
            "subsets (missing='ignore'); they are not mutually comparable");
    }
}

inline void validate_item_count(size_t n, const std::string& caller) {
    if (n == 0) {
        throw std::invalid_argument(caller + " requires at least one item");
    }
}

// ParallelFor silently does no work for a zero chunk, which would read as a
// result with nothing computed rather than as an error.
inline void validate_chunk_size(size_t chunk_size, const std::string& caller) {
    if (chunk_size == 0) {
        throw std::invalid_argument(caller + " chunk_size must be at least one");
    }
}

// The selection entry points' order: the item count before the chunk size.
inline void validate_size_and_chunk(size_t n, size_t chunk_size,
                                    const std::string& caller) {
    validate_item_count(n, caller);
    validate_chunk_size(chunk_size, caller);
}

// Follows the k-medoids precedent (KMedoids.cpp): a huge explicit request can
// terminate the process, so it is capped at the item count. Zero passes
// through and keeps meaning hardware concurrency; n >= 1 here, so the minimum
// never manufactures a zero from a nonzero request.
inline size_t capped_threads(size_t num_threads, size_t n) {
    return std::min(num_threads, n);
}

}  // namespace OECluster::detail

#endif  // OECLUSTER_SRC_CLUSTERING_DIVERSITYVALIDATION_H
