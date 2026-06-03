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

TEST(ClusterReportTest, BasicProfileTwoEqualClusters) {
    const DenseStorage storage = MakeTwoClusterStorage();
    const ClusteringResult result = MakeResult({0, 0, 1, 1});

    const ClusterReport r = cluster_report(result, storage, ClusterReportOptions());

    EXPECT_EQ(r.num_samples, 4u);
    EXPECT_EQ(r.num_clusters, 2u);
    EXPECT_EQ(r.num_noise, 0u);
    EXPECT_EQ(r.num_singletons, 0u);
    EXPECT_DOUBLE_EQ(r.noise_fraction, 0.0);
    EXPECT_DOUBLE_EQ(r.singleton_fraction, 0.0);
    EXPECT_DOUBLE_EQ(r.largest_cluster_fraction, 0.5);
    EXPECT_DOUBLE_EQ(r.cluster_size_median, 2.0);
    EXPECT_DOUBLE_EQ(r.cluster_size_p90, 2.0);
    EXPECT_DOUBLE_EQ(r.size_gini, 0.0);
    EXPECT_DOUBLE_EQ(r.size_entropy, 1.0);
}

TEST(ClusterReportTest, BasicProfileSkewedSizesGiniEntropy) {
    // Cluster 0 = {0,1,2}, cluster 1 = {3}. sizes [3,1].
    DenseStorage storage(4);
    for (size_t i = 0; i < 4; ++i) {
        for (size_t j = i + 1; j < 4; ++j) {
            storage.Set(i, j, 0.5);
        }
    }
    const ClusteringResult result = MakeResult({0, 0, 0, 1});

    const ClusterReport r = cluster_report(result, storage, ClusterReportOptions());

    EXPECT_EQ(r.num_clusters, 2u);
    EXPECT_EQ(r.num_singletons, 1u);
    EXPECT_DOUBLE_EQ(r.largest_cluster_fraction, 0.75);
    // sizes sorted [1,3]: gini = (2*(1*1 + 2*3))/(2*4) - 3/2 = 1.75 - 1.5 = 0.25
    EXPECT_NEAR(r.size_gini, 0.25, 1e-12);
    // p=[0.75,0.25]: H = -(0.75*log2 0.75 + 0.25*log2 0.25) bits
    const double expected_h = -(0.75 * std::log2(0.75) + 0.25 * std::log2(0.25));
    EXPECT_NEAR(r.size_entropy, expected_h, 1e-12);
}

TEST(ClusterReportTest, NoiseTreatedAsSingletonsToggle) {
    // One real cluster {0,1}; points 2,3 are noise.
    const DenseStorage storage = MakeTwoClusterStorage();
    const ClusteringResult result = MakeResult({0, 0, -1, -1});

    ClusterReportOptions fold;            // default treat_noise_as_singletons = true
    const ClusterReport rf = cluster_report(result, storage, fold);
    EXPECT_EQ(rf.num_clusters, 1u);
    EXPECT_EQ(rf.num_noise, 2u);
    EXPECT_EQ(rf.num_singletons, 0u);
    EXPECT_DOUBLE_EQ(rf.noise_fraction, 0.5);
    // folded: (0 + 2) / (1 + 2) = 2/3
    EXPECT_NEAR(rf.singleton_fraction, 2.0 / 3.0, 1e-12);

    ClusterReportOptions split(ClusterThreshold::Default);
    split.treat_noise_as_singletons = false;
    const ClusterReport rs = cluster_report(result, storage, split);
    // not folded: 0 / 1 = 0
    EXPECT_DOUBLE_EQ(rs.singleton_fraction, 0.0);
    EXPECT_EQ(rs.num_noise, 2u);  // still reported
}
