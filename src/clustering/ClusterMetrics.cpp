/**
 * @file ClusterMetrics.cpp
 * @brief Shared per-cluster distance statistics implementation.
 */

#include "ClusterMetrics.h"

#include <algorithm>

namespace OECluster::detail {

double mean_distance(const std::vector<double>& distances) {
    if (distances.empty()) {
        return 0.0;
    }
    double total = 0.0;
    for (const double distance : distances) {
        total += distance;
    }
    return total / static_cast<double>(distances.size());
}

double max_distance(const std::vector<double>& distances) {
    if (distances.empty()) {
        return 0.0;
    }
    return *std::max_element(distances.begin(), distances.end());
}

double median_distance(std::vector<double> distances) {
    if (distances.empty()) {
        return 0.0;
    }
    std::sort(distances.begin(), distances.end());
    const size_t middle = distances.size() / 2;
    if (distances.size() % 2 == 1) {
        return distances[middle];
    }
    return (distances[middle - 1] + distances[middle]) / 2.0;
}

double cluster_diameter(const Cluster& cluster, const StorageBackend& storage) {
    double diameter = 0.0;
    for (size_t i = 0; i < cluster.size(); ++i) {
        for (size_t j = i + 1; j < cluster.size(); ++j) {
            diameter = std::max(diameter, storage.Get(cluster[i], cluster[j]));
        }
    }
    return diameter;
}

}  // namespace OECluster::detail
