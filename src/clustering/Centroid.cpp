/**
 * @file Centroid.cpp
 * @brief Representative member selection for clusters.
 */

#include "oecluster/clustering/Centroid.h"

#include <algorithm>
#include <stdexcept>
#include <string>
#include <unordered_set>

namespace OECluster {

namespace {

void validate_cluster_members(const Cluster& cluster, size_t num_items) {
    if (cluster.empty()) {
        throw std::invalid_argument("Cluster centroid requires at least one member");
    }

    std::unordered_set<size_t> seen;
    seen.reserve(cluster.size());
    for (const size_t member : cluster) {
        if (member >= num_items) {
            throw std::out_of_range("Cluster member index is outside the storage range");
        }
        if (!seen.insert(member).second) {
            throw std::invalid_argument("Cluster members must be unique");
        }
    }
}

void validate_distance_storage(const StorageBackend& storage, CentroidMethod method) {
    if (method == CentroidMethod::First) {
        return;
    }
    if (dynamic_cast<const SparseStorage*>(&storage) != nullptr) {
        throw std::invalid_argument(
            "SparseStorage cannot provide complete distances for centroid methods");
    }
}

size_t medoid_centroid(const Cluster& cluster, const StorageBackend& storage) {
    size_t best_member = cluster.front();
    double best_score = 0.0;
    bool has_best = false;

    for (const size_t candidate : cluster) {
        double score = 0.0;
        for (const size_t member : cluster) {
            if (member != candidate) {
                score += storage.Get(candidate, member);
            }
        }
        if (!has_best || score < best_score) {
            best_member = candidate;
            best_score = score;
            has_best = true;
        }
    }

    return best_member;
}

size_t minimax_centroid(const Cluster& cluster, const StorageBackend& storage) {
    size_t best_member = cluster.front();
    double best_score = 0.0;
    bool has_best = false;

    for (const size_t candidate : cluster) {
        double score = 0.0;
        for (const size_t member : cluster) {
            if (member != candidate) {
                score = std::max(score, storage.Get(candidate, member));
            }
        }
        if (!has_best || score < best_score) {
            best_member = candidate;
            best_score = score;
            has_best = true;
        }
    }

    return best_member;
}

}  // namespace

size_t cluster_centroid(
    const Cluster& cluster,
    const StorageBackend& storage,
    CentroidMethod method) {
    validate_cluster_members(cluster, storage.NumItems());
    if (cluster.size() == 1) {
        return cluster.front();
    }
    validate_distance_storage(storage, method);

    switch (method) {
        case CentroidMethod::First:
            return cluster.front();
        case CentroidMethod::Medoid:
        case CentroidMethod::Mean:
            return medoid_centroid(cluster, storage);
        case CentroidMethod::Minimax:
            return minimax_centroid(cluster, storage);
    }

    throw std::invalid_argument("Unknown centroid method");
}

}  // namespace OECluster
