#include <gtest/gtest.h>

#include <cmath>
#include <vector>

#include "oecluster/StorageBackend.h"
#include "oecluster/clustering/ClusterReport.h"
#include "oecluster/clustering/ClusterTypes.h"

using namespace OECluster;

namespace {

// Two clusters {0,1} and {2,3}; intra 0.2, all cross pairs 0.8.
DenseStorage MakeTwoClusterStorage() {
    DenseStorage storage(4);
    storage.Set(0, 1, 0.2);
    storage.Set(2, 3, 0.2);
    storage.Set(0, 2, 0.8);
    storage.Set(0, 3, 0.8);
    storage.Set(1, 2, 0.8);
    storage.Set(1, 3, 0.8);
    return storage;
}

// A ClusteringResult with explicit labels/members for testing.
ClusteringResult MakeResult(std::vector<ClusterLabel> labels) {
    Clusters members = labels_to_clusters(labels);
    return ClusteringResult(std::move(labels), std::move(members));
}

}  // namespace

TEST(ClusterReportOptionsTest, PresetsSeedDocumentedThresholds) {
    const ClusterReportOptions def(ClusterThreshold::Default);
    EXPECT_EQ(def.coverage_thresholds, std::vector<double>({0.25, 0.35, 0.45}));
    EXPECT_DOUBLE_EQ(def.boundary_threshold, 0.30);

    const ClusterReportOptions tight(ClusterThreshold::Tight);
    EXPECT_EQ(tight.coverage_thresholds, std::vector<double>({0.20, 0.30, 0.40}));
    EXPECT_DOUBLE_EQ(tight.boundary_threshold, 0.25);

    const ClusterReportOptions diversity(ClusterThreshold::Diversity);
    EXPECT_EQ(diversity.coverage_thresholds, std::vector<double>({0.40, 0.50, 0.60}));
    EXPECT_DOUBLE_EQ(diversity.boundary_threshold, 0.40);

    // Default-constructed equals the Default preset.
    const ClusterReportOptions implicit;
    EXPECT_EQ(implicit.coverage_thresholds, def.coverage_thresholds);
    EXPECT_DOUBLE_EQ(implicit.boundary_threshold, def.boundary_threshold);
}

TEST(ClusterReportTest, SparseStorageThrows) {
    SparseStorage storage(3, 0.5);
    storage.Set(0, 1, 0.1);
    storage.Finalize();
    const ClusteringResult result = MakeResult({0, 0, -1});
    EXPECT_THROW(cluster_report(result, storage, ClusterReportOptions()),
                 std::invalid_argument);
}
