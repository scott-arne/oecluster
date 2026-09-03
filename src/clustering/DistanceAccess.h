/**
 * @file DistanceAccess.h
 * @brief Internal helpers for dense precomputed distance access.
 */

#ifndef OECLUSTER_SRC_CLUSTERING_DISTANCEACCESS_H
#define OECLUSTER_SRC_CLUSTERING_DISTANCEACCESS_H

#include <cstddef>
#include <stdexcept>
#include <string>
#include <unordered_set>

#include "oecluster/StorageBackend.h"
#include "oecluster/clustering/ClusterTypes.h"

namespace OECluster::detail {

inline size_t condensed_index(size_t n, size_t i, size_t j) {
    if (i == j) {
        throw std::invalid_argument("condensed_index requires distinct indices");
    }
    if (i > j) {
        const size_t tmp = i;
        i = j;
        j = tmp;
    }
    return n * i - i * (i + 1) / 2 + j - i - 1;
}

inline double dense_distance(const double* data, size_t n, size_t i, size_t j) {
    if (i == j) {
        return 0.0;
    }
    if (i > j) {
        const size_t tmp = i;
        i = j;
        j = tmp;
    }
    return data[condensed_index(n, i, j)];
}

/**
 * @brief Refuse a cluster whose members cannot name items in the storage.
 *
 * Shared so that every path which reads storage for a caller-supplied cluster
 * runs it *before* its first ``storage.Get``. A bad member index otherwise
 * trips the backend's own range check first, which reports a storage class
 * rather than the cluster the caller has to fix. A path that returns without
 * reading storage skips it -- ``select_representatives`` with ``k == 0``
 * answers an unusable cluster with an empty result rather than a diagnosis.
 *
 * The messages name the cluster, not the operation that happens to be running:
 * ``rank_representatives`` and ``cluster_report`` are both callers, and the
 * report's caller never asked for a representative.
 *
 * :param cluster: Member indices to validate.
 * :param num_samples: Number of samples the storage holds.
 * :raises std::invalid_argument: If the cluster is empty or repeats a member.
 * :raises std::out_of_range: If a member is at or beyond num_samples.
 */
inline void validate_cluster_members(const Cluster& cluster, const size_t num_samples) {
    if (cluster.empty()) {
        throw std::invalid_argument("Cluster must contain at least one member");
    }

    std::unordered_set<size_t> seen;
    seen.reserve(cluster.size());
    for (const size_t member : cluster) {
        if (member >= num_samples) {
            throw std::out_of_range("Cluster member index is outside the storage range");
        }
        if (!seen.insert(member).second) {
            throw std::invalid_argument("Cluster members must be unique");
        }
    }
}

inline void validate_complete_distance_storage(
    const StorageBackend& storage,
    const std::string& algorithm_name) {
    if (dynamic_cast<const SparseStorage*>(&storage) != nullptr) {
        throw std::invalid_argument(
            algorithm_name + " requires complete pairwise distances; "
            "SparseStorage is not supported");
    }
    if (storage.NumPairs() > 0 && storage.Data() == nullptr) {
        throw std::invalid_argument(
            algorithm_name + " requires contiguous dense or memory-mapped storage");
    }
}

}  // namespace OECluster::detail

#endif  // OECLUSTER_SRC_CLUSTERING_DISTANCEACCESS_H
