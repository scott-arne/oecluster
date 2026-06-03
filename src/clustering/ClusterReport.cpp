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

ClusterReport cluster_report(
    const ClusteringResult& result,
    const StorageBackend& storage,
    const ClusterReportOptions& options) {
    detail::validate_complete_distance_storage(storage, "cluster_report");

    ClusterReport report;
    report.coverage_thresholds = options.coverage_thresholds;
    return report;
}

ClusterReportComparison compare_reports(const ClusterReport& a, const ClusterReport& b) {
    ClusterReportComparison comparison;
    comparison.a = a;
    comparison.b = b;
    return comparison;
}

}  // namespace OECluster
