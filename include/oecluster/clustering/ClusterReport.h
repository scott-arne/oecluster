/**
 * @file ClusterReport.h
 * @brief Method-agnostic clustering-quality report and comparison.
 */

#ifndef OECLUSTER_CLUSTERING_CLUSTERREPORT_H
#define OECLUSTER_CLUSTERING_CLUSTERREPORT_H

#include <cstddef>
#include <vector>

#include "oecluster/StorageBackend.h"
#include "oecluster/clustering/ClusterTypes.h"
#include "oecluster/clustering/Representative.h"

namespace OECluster {

/**
 * @brief Threshold preset for clustering-quality reports.
 *
 * Distances are Tanimoto/Jaccard (distance = 1 - similarity).
 */
enum class ClusterThreshold {
    Default,    ///< coverage {0.25,0.35,0.45}, boundary 0.30.
    Tight,      ///< coverage {0.20,0.30,0.40}, boundary 0.25.
    Diversity   ///< coverage {0.40,0.50,0.60}, boundary 0.40.
};

/**
 * @brief Options controlling clustering-quality report computation.
 */
struct ClusterReportOptions {
    std::vector<double> coverage_thresholds;
    double boundary_threshold = 0.30;
    RepresentativeMethod representative_method = RepresentativeMethod::Medoid;
    bool treat_noise_as_singletons = true;
    size_t num_threads = 0;

    /** @brief Seed coverage_thresholds and boundary_threshold from a preset. */
    explicit ClusterReportOptions(ClusterThreshold preset = ClusterThreshold::Default);
};

/**
 * @brief Clustering-quality scorecard. NaN marks an undefined metric.
 */
struct ClusterReport {
    // Basic profile.
    size_t num_samples = 0;
    size_t num_clusters = 0;
    size_t num_noise = 0;
    size_t num_singletons = 0;
    double noise_fraction = 0.0;
    double singleton_fraction = 0.0;
    double largest_cluster_fraction = 0.0;
    double cluster_size_median = 0.0;
    double cluster_size_p90 = 0.0;
    double size_gini = 0.0;
    double size_entropy = 0.0;
    // Compactness / separation.
    double mean_intra_distance = 0.0;
    double median_intra_distance = 0.0;
    double median_radius = 0.0;
    double p95_diameter = 0.0;
    double silhouette = 0.0;
    double dunn_index = 0.0;
    size_t boundary_violations = 0;
    // Representative / coverage.
    double median_medoid_member_distance = 0.0;
    double representative_redundancy = 0.0;
    std::vector<double> coverage_thresholds;
    std::vector<double> coverage_at;
};

/**
 * @brief Two reports aligned for side-by-side comparison.
 */
struct ClusterReportComparison {
    ClusterReport a;
    ClusterReport b;
};

/**
 * @brief Compute a clustering-quality report.
 *
 * :param result: Any algorithm's clustering result (base ClusteringResult).
 * :param storage: Complete pairwise distance storage.
 * :param options: Report options (thresholds, representative method, flags).
 * :returns: A ClusterReport scorecard.
 * :raises std::invalid_argument: If storage cannot provide complete distances.
 */
ClusterReport cluster_report(
    const ClusteringResult& result,
    const StorageBackend& storage,
    const ClusterReportOptions& options = ClusterReportOptions());

/**
 * @brief Pair two reports for side-by-side reading. No agreement math.
 */
ClusterReportComparison compare_reports(const ClusterReport& a, const ClusterReport& b);

}  // namespace OECluster

#endif  // OECLUSTER_CLUSTERING_CLUSTERREPORT_H
