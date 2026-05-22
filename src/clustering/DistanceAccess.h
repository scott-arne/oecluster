/**
 * @file DistanceAccess.h
 * @brief Internal helpers for dense precomputed distance access.
 */

#ifndef OECLUSTER_SRC_CLUSTERING_DISTANCEACCESS_H
#define OECLUSTER_SRC_CLUSTERING_DISTANCEACCESS_H

#include <cstddef>
#include <stdexcept>
#include <string>

#include "oecluster/StorageBackend.h"

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
