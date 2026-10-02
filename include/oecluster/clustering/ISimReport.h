/**
 * @file ISimReport.h
 * @brief iSIM set similarity and the approximate fingerprint-native cluster report.
 *
 * iSIM values are union-weighted: summed pairwise intersections over summed
 * pairwise unions. They are exact for that quantity and are not the mean of
 * per-pair Tanimoto values.
 */

#ifndef OECLUSTER_CLUSTERING_ISIMREPORT_H
#define OECLUSTER_CLUSTERING_ISIMREPORT_H

#include <cstddef>
#include <limits>
#include <string>
#include <vector>

#include "oefp/batch.h"
#include "oecluster/clustering/ClusterReport.h"
#include "oecluster/clustering/ClusterTypes.h"

namespace OECluster {

/// Options for isim(). metric accepts only "tanimoto"; it is reserved so
/// metrics whose iSIM is the exact mean pair similarity can be added later.
struct ISimOptions {
    std::string metric = "tanimoto";
};

/**
 * @brief iSIM Tanimoto similarity of a set of binary fingerprints.
 *
 * NaN for fewer than two fingerprints; 1.0 when every fingerprint is all-zero.
 *
 * @throws std::invalid_argument for a metric other than "tanimoto", a
 *         non-empty batch of zero-width fingerprints, fingerprints of
 *         2^31 or more bits, or 2^32 or more fingerprints.
 */
double isim(const OEFP::OEFPBatch& fingerprints, const ISimOptions& options = ISimOptions());

/**
 * @brief Options for isim_report().
 *
 * coverage_thresholds is seeded from the preset with the same table
 * cluster_report uses. A NaN entry is refused; a negative one is accepted and
 * covers no sample.
 */
struct ISimReportOptions {
    std::string metric = "tanimoto";
    std::vector<double> coverage_thresholds;
    bool treat_noise_as_singletons = true;
    /// Runs the O(N K) centroid stage: silhouette, nearest cluster, medoid
    /// Davies-Bouldin and Dunn, and coverage.
    bool compute_centroid_indices = false;
    bool compute_per_cluster_records = false;
    /// 0 means hardware concurrency; capped at the cluster count.
    size_t num_threads = 0;

    explicit ISimReportOptions(ClusterThreshold preset = ClusterThreshold::Default);
};

/// What the caller asked for, set even when the answer is undefined.
struct ISimReportRequested {
    bool centroid_indices = false;
    bool per_cluster_records = false;
};

/**
 * @brief One cluster's row, in Members() order. label and nearest_cluster are
 *        ordinals into Members().
 */
struct ISimClusterRecord {
    ClusterLabel label = 0;
    size_t size = 0;
    size_t medoid = 0;  ///< Sample index of the iSIM medoid.
    /// 1 - iSIM of the cluster; NaN for a singleton.
    double isim_intra_distance = std::numeric_limits<double>::quiet_NaN();
    /// 1 - iSIM of the cluster against every clustered point outside it; NaN when K < 2.
    double isim_separation = std::numeric_limits<double>::quiet_NaN();
    /// Largest medoid-to-member distance; 0.0 for a singleton.
    double radius = std::numeric_limits<double>::quiet_NaN();
    /// Mean medoid-to-other-member distance (n_k - 1 denominator); 0.0 for a singleton.
    double mean_medoid_distance = std::numeric_limits<double>::quiet_NaN();
    /// Centroid stage only; NaN when K < 2, 0 for a singleton cluster.
    double isim_silhouette = std::numeric_limits<double>::quiet_NaN();
    /// Centroid stage only: highest iSIM cluster-to-cluster similarity, ties to
    /// the lowest ordinal. Not cluster_report's single-linkage nearest.
    ClusterLabel nearest_cluster = NO_NEAREST_CLUSTER;
    double nearest_cluster_similarity = std::numeric_limits<double>::quiet_NaN();
};

/**
 * @brief Approximate clustering-quality report from binary fingerprints.
 *
 * Profile fields are exact; isim_* fields are iSIM ratios; median_radius,
 * median_medoid_member_distance, calinski_harabasz_medoid, the medoid
 * Davies-Bouldin and Dunn and coverage are exact relative to the iSIM medoids.
 * Undefined values are NaN. Noise is excluded from every quality field except
 * coverage_at and noise_coverage_at, which measure every sample against the
 * cluster medoids.
 */
struct ISimReport {
    size_t num_samples = 0;
    size_t num_clusters = 0;
    size_t num_noise = 0;
    size_t num_singletons = 0;
    double noise_fraction = std::numeric_limits<double>::quiet_NaN();
    double singleton_fraction = std::numeric_limits<double>::quiet_NaN();
    double largest_cluster_fraction = std::numeric_limits<double>::quiet_NaN();
    double cluster_size_median = std::numeric_limits<double>::quiet_NaN();
    double cluster_size_p90 = std::numeric_limits<double>::quiet_NaN();
    double size_gini = std::numeric_limits<double>::quiet_NaN();
    double size_entropy = std::numeric_limits<double>::quiet_NaN();

    double isim_intra_distance = std::numeric_limits<double>::quiet_NaN();
    double isim_inter_distance = std::numeric_limits<double>::quiet_NaN();

    double median_radius = std::numeric_limits<double>::quiet_NaN();
    double median_medoid_member_distance = std::numeric_limits<double>::quiet_NaN();
    double calinski_harabasz_medoid = std::numeric_limits<double>::quiet_NaN();

    double isim_silhouette = std::numeric_limits<double>::quiet_NaN();
    double davies_bouldin_medoid = std::numeric_limits<double>::quiet_NaN();
    double dunn_medoid_separation_medoid_spread = std::numeric_limits<double>::quiet_NaN();
    /// Echoed whether or not the centroid stage runs.
    std::vector<double> coverage_thresholds;
    /// Centroid stage only: fraction of all samples within each threshold of
    /// their nearest medoid.
    std::vector<double> coverage_at;
    /// Centroid stage only: the same over noise samples; NaN entries when there is no noise.
    std::vector<double> noise_coverage_at;

    std::vector<ISimClusterRecord> records;
    ISimReportRequested requested;
};

/**
 * @brief Approximate clustering-quality report over binary fingerprints.
 *
 * Core cost O(N words + sum of popcounts); the centroid stage adds
 * O(sum of popcounts * K + K^2 words + N K words) time and K * bits * 4 bytes.
 *
 * @throws std::invalid_argument for a metric other than "tanimoto", a
 *         non-empty zero-width batch, 2^31 or more bits, 2^32 or more fingerprints, a NaN
 *         coverage threshold, or a malformed partition.
 * @throws std::out_of_range when the result has clusters and more labels than
 *         the batch has fingerprints, or a member indexes past the labels.
 */
ISimReport isim_report(const ClusteringResult& result, const OEFP::OEFPBatch& fingerprints,
                       const ISimReportOptions& options = ISimReportOptions());

}  // namespace OECluster

#endif  // OECLUSTER_CLUSTERING_ISIMREPORT_H
