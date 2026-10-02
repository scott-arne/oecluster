/**
 * @file ClusterReportTuned.h
 * @brief cluster_report with its memory budgets exposed for tests.
 */

#ifndef OECLUSTER_SRC_CLUSTERING_CLUSTERREPORTTUNED_H
#define OECLUSTER_SRC_CLUSTERING_CLUSTERREPORTTUNED_H

#include <cstddef>
#include <functional>

#include "ExactMedian.h"
#include "ReportDistanceSource.h"
#include "oecluster/StorageBackend.h"
#include "oecluster/clustering/ClusterReport.h"

namespace OECluster::detail {

/**
 * @brief Budgets the report engine normally fixes at 2^20.
 *
 * Both budgets only bite above a million pairs, far beyond what a unit test
 * can afford, so tests shrink them to drive the radix median route and the
 * block cap on small inputs. The public entry points pass the defaults.
 */
struct ReportTuning {
    size_t median_direct_budget = MEDIAN_DIRECT_BUDGET;
    size_t fill_block_distances = FILL_BLOCK_DISTANCES;
    std::function<void(size_t, size_t)> on_block;
};

ClusterReport cluster_report_tuned(
    const ClusteringResult& result,
    const StorageBackend& storage,
    const ClusterReportOptions& options,
    const ReportTuning& tuning);

}  // namespace OECluster::detail

#endif  // OECLUSTER_SRC_CLUSTERING_CLUSTERREPORTTUNED_H
