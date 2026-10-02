/**
 * @file ReportCommon.h
 * @brief Profile, partition validation and helpers shared by the cluster reports.
 *
 * cluster_report and isim_report must agree on the exact profile, on which
 * partitions they refuse and on how they refuse them. One definition keeps the
 * two from drifting; the caller and noun parameters keep each message naming
 * the entry point the caller actually called.
 */

#ifndef OECLUSTER_CLUSTERING_REPORTCOMMON_H
#define OECLUSTER_CLUSTERING_REPORTCOMMON_H

#include <cstddef>
#include <limits>
#include <string>
#include <vector>

#include "oecluster/clustering/ClusterReport.h"
#include "oecluster/clustering/ClusterTypes.h"

namespace OECluster::detail {

/// The coverage thresholds a preset seeds.
std::vector<double> preset_coverage_thresholds(ClusterThreshold preset);

/// Fractional-rank percentile (NumPy/Pandas default interpolation); NaN on empty.
double percentile(std::vector<double> values, double q);

/// Rousseeuw's s(i): 0 for a member of a singleton cluster and 0 when
/// max(a, b) == 0.
double silhouette_term(double a_term, double b_term, size_t cluster_size);

/// The exact size profile both reports publish.
struct ReportProfile {
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
};

ReportProfile report_profile(const ClusteringResult& result, bool treat_noise_as_singletons);

/// Copies a profile into any report carrying the eleven profile fields.
template <class Report>
void assign_profile(Report& report, const ReportProfile& profile) {
    report.num_samples = profile.num_samples;
    report.num_clusters = profile.num_clusters;
    report.num_noise = profile.num_noise;
    report.num_singletons = profile.num_singletons;
    report.noise_fraction = profile.noise_fraction;
    report.singleton_fraction = profile.singleton_fraction;
    report.largest_cluster_fraction = profile.largest_cluster_fraction;
    report.cluster_size_median = profile.cluster_size_median;
    report.cluster_size_p90 = profile.cluster_size_p90;
    report.size_gini = profile.size_gini;
    report.size_entropy = profile.size_entropy;
}

/// Owner entry of a sample no cluster lists (noise).
constexpr size_t NO_OWNER = std::numeric_limits<size_t>::max();

/**
 * @brief Refuses a result whose labels and members are not one partition of
 *        num_items samples, before any distance or fingerprint is read.
 *
 * @returns owner[i], the Members() ordinal holding sample i, or NO_OWNER.
 * @throws std::out_of_range when members is non-empty and the label count
 *         exceeds num_items, or a member is at or beyond the label count.
 * @throws std::invalid_argument for a malformed cluster or a partition that
 *         disagrees with the labels.
 */
std::vector<size_t> validate_report_partition(const ClusteringResult& result,
                                              size_t num_items,
                                              const std::string& caller,
                                              const std::string& noun);

}  // namespace OECluster::detail

#endif  // OECLUSTER_CLUSTERING_REPORTCOMMON_H
