/**
 * @file ClusterTypes.h
 * @brief Shared clustering result types.
 */

#ifndef OECLUSTER_CLUSTERING_CLUSTERTYPES_H
#define OECLUSTER_CLUSTERING_CLUSTERTYPES_H

#include <cstddef>
#include <utility>
#include <vector>

namespace OECluster {

using Cluster = std::vector<size_t>;
using Clusters = std::vector<Cluster>;
using ClusterLabel = int;

constexpr ClusterLabel NOISE_LABEL = -1;

/**
 * @brief Shared clustering result with integer labels and grouped members.
 *
 * Labels use -1 for noise. Non-negative labels identify clusters. Results are
 * immutable: construct with labels and members, then read via const getters.
 */
class ClusteringResult {
public:
    ClusteringResult() = default;
    ClusteringResult(std::vector<ClusterLabel> labels, Clusters members)
        : labels_(std::move(labels)), members_(std::move(members)) {}
    virtual ~ClusteringResult() = default;

    /** @brief Per-item labels; -1 marks noise. */
    const std::vector<ClusterLabel>& Labels() const { return labels_; }
    /** @brief Cluster member-index lists. */
    const Clusters& Members() const { return members_; }
    /** @brief Number of clusters. */
    size_t NumClusters() const { return members_.size(); }
    /** @brief Number of items (length of the labels vector). */
    size_t NumSamples() const { return labels_.size(); }

protected:
    std::vector<ClusterLabel> labels_;
    Clusters members_;
};

/**
 * @brief Convert cluster labels into cluster member lists.
 *
 * :param labels: Per-item labels, where negative values are ignored as noise.
 * :returns: Clusters ordered by non-negative label, preserving item order.
 */
Clusters labels_to_clusters(const std::vector<ClusterLabel>& labels);

}  // namespace OECluster

#endif  // OECLUSTER_CLUSTERING_CLUSTERTYPES_H
