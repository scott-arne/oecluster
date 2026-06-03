/**
 * @file ClusterReport.cpp
 * @brief Clustering-quality report implementation.
 */

#include "oecluster/clustering/ClusterReport.h"

#include <algorithm>
#include <cmath>
#include <limits>

#include "ClusterMetrics.h"
#include "DistanceAccess.h"

namespace OECluster {

ClusterReportOptions::ClusterReportOptions(ClusterThreshold preset) {
    switch (preset) {
        case ClusterThreshold::Tight:
            coverage_thresholds = {0.20, 0.30, 0.40};
            boundary_threshold = 0.25;
            break;
        case ClusterThreshold::Diversity:
            coverage_thresholds = {0.40, 0.50, 0.60};
            boundary_threshold = 0.40;
            break;
        case ClusterThreshold::Default:
        default:
            coverage_thresholds = {0.25, 0.35, 0.45};
            boundary_threshold = 0.30;
            break;
    }
}

namespace {

double percentile(std::vector<double> values, const double q) {
    if (values.empty()) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    std::sort(values.begin(), values.end());
    if (values.size() == 1) {
        return values.front();
    }
    const double rank = q * static_cast<double>(values.size() - 1);
    const size_t lo = static_cast<size_t>(std::floor(rank));
    const size_t hi = static_cast<size_t>(std::ceil(rank));
    const double frac = rank - static_cast<double>(lo);
    return values[lo] + frac * (values[hi] - values[lo]);
}

double size_gini(const std::vector<double>& sizes) {
    const size_t n = sizes.size();
    if (n <= 1) {
        return 0.0;
    }
    std::vector<double> sorted = sizes;
    std::sort(sorted.begin(), sorted.end());
    double total = 0.0;
    double weighted = 0.0;
    for (size_t i = 0; i < n; ++i) {
        total += sorted[i];
        weighted += static_cast<double>(i + 1) * sorted[i];
    }
    if (total == 0.0) {
        return 0.0;
    }
    return (2.0 * weighted) / (static_cast<double>(n) * total) -
           (static_cast<double>(n) + 1.0) / static_cast<double>(n);
}

double size_entropy(const std::vector<double>& sizes) {
    if (sizes.size() <= 1) {
        return 0.0;
    }
    double total = 0.0;
    for (const double s : sizes) {
        total += s;
    }
    if (total == 0.0) {
        return 0.0;
    }
    double entropy = 0.0;
    for (const double s : sizes) {
        if (s > 0.0) {
            const double p = s / total;
            entropy -= p * std::log2(p);
        }
    }
    return entropy;
}

}  // namespace

ClusterReport cluster_report(
    const ClusteringResult& result,
    const StorageBackend& storage,
    const ClusterReportOptions& options) {
    detail::validate_complete_distance_storage(storage, "cluster_report");

    ClusterReport report;
    report.coverage_thresholds = options.coverage_thresholds;

    const std::vector<ClusterLabel>& labels = result.Labels();
    const Clusters& members = result.Members();

    report.num_samples = labels.size();
    report.num_clusters = members.size();

    size_t num_noise = 0;
    for (const ClusterLabel label : labels) {
        if (label < 0) {
            ++num_noise;
        }
    }
    report.num_noise = num_noise;

    std::vector<double> sizes;
    sizes.reserve(members.size());
    size_t largest = 0;
    size_t num_singletons = 0;
    for (const Cluster& cluster : members) {
        sizes.push_back(static_cast<double>(cluster.size()));
        largest = std::max(largest, cluster.size());
        if (cluster.size() == 1) {
            ++num_singletons;
        }
    }
    report.num_singletons = num_singletons;

    const double n = static_cast<double>(report.num_samples);
    report.noise_fraction = n > 0.0 ? static_cast<double>(num_noise) / n : 0.0;
    report.largest_cluster_fraction = n > 0.0 ? static_cast<double>(largest) / n : 0.0;

    if (options.treat_noise_as_singletons) {
        const double denom = static_cast<double>(report.num_clusters + num_noise);
        report.singleton_fraction =
            denom > 0.0 ? static_cast<double>(num_singletons + num_noise) / denom : 0.0;
    } else {
        const double denom = static_cast<double>(report.num_clusters);
        report.singleton_fraction =
            denom > 0.0 ? static_cast<double>(num_singletons) / denom : 0.0;
    }

    if (sizes.empty()) {
        report.cluster_size_median = std::numeric_limits<double>::quiet_NaN();
        report.cluster_size_p90 = std::numeric_limits<double>::quiet_NaN();
    } else {
        report.cluster_size_median = percentile(sizes, 0.5);
        report.cluster_size_p90 = percentile(sizes, 0.9);
    }
    report.size_gini = size_gini(sizes);
    report.size_entropy = size_entropy(sizes);

    return report;
}

ClusterReportComparison compare_reports(const ClusterReport& a, const ClusterReport& b) {
    ClusterReportComparison comparison;
    comparison.a = a;
    comparison.b = b;
    return comparison;
}

}  // namespace OECluster
