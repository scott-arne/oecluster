/**
 * @file ClusterTypes.cpp
 * @brief Shared clustering result helpers.
 */

#include "oecluster/clustering/ClusterTypes.h"

namespace OECluster {

Clusters labels_to_clusters(const std::vector<ClusterLabel>& labels) {
    ClusterLabel max_label = NOISE_LABEL;
    for (const ClusterLabel label : labels) {
        if (label > max_label) {
            max_label = label;
        }
    }
    if (max_label < 0) {
        return Clusters{};
    }

    const size_t cluster_count = static_cast<size_t>(max_label) + 1;
    Clusters clusters(cluster_count);
    for (size_t index = 0; index < labels.size(); ++index) {
        const ClusterLabel label = labels[index];
        if (label >= 0) {
            clusters[static_cast<size_t>(label)].push_back(index);
        }
    }

    return clusters;
}

}  // namespace OECluster
