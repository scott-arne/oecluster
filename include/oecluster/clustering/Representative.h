/**
 * @file Representative.h
 * @brief Representative member selection and ranking for clusters.
 */

#ifndef OECLUSTER_CLUSTERING_REPRESENTATIVE_H
#define OECLUSTER_CLUSTERING_REPRESENTATIVE_H

#include <cstddef>
#include <string>
#include <vector>

#include "oecluster/StorageBackend.h"
#include "oecluster/clustering/ClusterTypes.h"

namespace OECluster {

/**
 * @brief Algorithms for scoring cluster-member representatives.
 */
enum class RepresentativeMethod {
    Medoid,               ///< Minimize mean in-cluster distance.
    Minimax,              ///< Minimize maximum in-cluster distance.
    HighestNeighborhood,  ///< Maximize threshold-neighbor fraction.
    WeightedMedoid        ///< Combine centrality with penalty and priority metadata.
};

/**
 * @brief Strategy for selecting k representatives from a ranked cluster.
 */
enum class RepresentativeSelection {
    Score,    ///< Return the top k representatives by score.
    Diversity ///< Greedily spread representatives across the cluster.
};

/**
 * @brief Weights for weighted-medoid scoring.
 */
struct RepresentativeWeights {
    double alpha = 1.0;
    double beta = 1.0;
    double gamma = 1.0;
};

/**
 * @brief Quality metrics for one representative member.
 */
struct RepresentativeMetrics {
    double mean_distance_to_cluster = 0.0;
    double max_distance_to_cluster = 0.0;
    double median_distance_to_cluster = 0.0;
    double neighbor_fraction_at_threshold = 0.0;
    double nearest_external_distance = 0.0;
    double cluster_radius = 0.0;
    double cluster_diameter = 0.0;
    double silhouette_like_score = 0.0;
    double scaffold_purity = 0.0;
    size_t representative_rank = 0;
};

/**
 * @brief A scored representative member and its quality metrics.
 */
struct ClusterRepresentative {
    size_t member = 0;
    double score = 0.0;
    RepresentativeMetrics metrics;
};

/**
 * @brief Options for ranking and selecting representative members.
 */
struct RepresentativeOptions {
    RepresentativeMethod method = RepresentativeMethod::Medoid;
    RepresentativeSelection selection = RepresentativeSelection::Score;
    double neighbor_threshold = -1.0;
    RepresentativeWeights weights;
    std::vector<double> liability_penalties;
    std::vector<double> priority_scores;
    std::vector<std::string> scaffold_labels;
};

/**
 * @brief Rank every member of a cluster as a possible representative.
 *
 * :param cluster: Cluster member indices.
 * :param storage: Pairwise distance storage used for scoring and metrics.
 * :param options: Representative scoring options.
 * :returns: Representatives ordered from best to worst under the selected method.
 * :raises std::invalid_argument: If inputs are invalid or sparse distances are incomplete.
 * :raises std::out_of_range: If a cluster member is outside the storage range.
 */
std::vector<ClusterRepresentative> rank_representatives(
    const Cluster& cluster,
    const StorageBackend& storage,
    const RepresentativeOptions& options = RepresentativeOptions());

/**
 * @brief Select up to k representatives from a cluster.
 *
 * :param cluster: Cluster member indices.
 * :param storage: Pairwise distance storage used for scoring and metrics.
 * :param k: Maximum number of representatives to return.
 * :param options: Representative scoring and selection options.
 * :returns: Selected representatives in selection order.
 */
std::vector<ClusterRepresentative> select_representatives(
    const Cluster& cluster,
    const StorageBackend& storage,
    size_t k,
    const RepresentativeOptions& options = RepresentativeOptions());

/**
 * @brief Select the best representative member index from a cluster.
 *
 * :param cluster: Cluster member indices.
 * :param storage: Pairwise distance storage used for scoring.
 * :param options: Representative scoring options.
 * :returns: Member index selected as the cluster representative.
 */
size_t cluster_representative(
    const Cluster& cluster,
    const StorageBackend& storage,
    const RepresentativeOptions& options = RepresentativeOptions());

/**
 * @brief Select the best representative member index with a scoring method.
 *
 * :param cluster: Cluster member indices.
 * :param storage: Pairwise distance storage used for scoring.
 * :param method: Representative scoring method.
 * :returns: Member index selected as the cluster representative.
 */
size_t cluster_representative(
    const Cluster& cluster,
    const StorageBackend& storage,
    RepresentativeMethod method);

}  // namespace OECluster

#endif  // OECLUSTER_CLUSTERING_REPRESENTATIVE_H
