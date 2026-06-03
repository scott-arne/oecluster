#include <gtest/gtest.h>

#include <cmath>
#include <vector>

#include "oecluster/StorageBackend.h"
#include "oecluster/clustering/Agglomerative.h"
#include "oecluster/clustering/Butina.h"
#include "oecluster/clustering/ClusterReport.h"
#include "oecluster/clustering/ClusterTypes.h"
#include "oecluster/clustering/DBSCAN.h"
#include "oecluster/clustering/HDBSCAN.h"

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

TEST(ClusterReportTest, CompactnessTwoClusters) {
    const DenseStorage storage = MakeTwoClusterStorage();
    const ClusteringResult result = MakeResult({0, 0, 1, 1});

    const ClusterReport r = cluster_report(result, storage, ClusterReportOptions());

    EXPECT_DOUBLE_EQ(r.mean_intra_distance, 0.2);
    EXPECT_DOUBLE_EQ(r.median_intra_distance, 0.2);
    EXPECT_DOUBLE_EQ(r.median_radius, 0.2);
    EXPECT_DOUBLE_EQ(r.p95_diameter, 0.2);
    EXPECT_NEAR(r.silhouette, 0.75, 1e-12);
    EXPECT_NEAR(r.dunn_index, 4.0, 1e-12);
}

TEST(ClusterReportTest, BoundaryViolationsThreshold) {
    const DenseStorage storage = MakeTwoClusterStorage();
    const ClusteringResult result = MakeResult({0, 0, 1, 1});

    ClusterReportOptions strict;
    strict.boundary_threshold = 0.3;  // cross pairs are 0.8 > 0.3
    EXPECT_EQ(cluster_report(result, storage, strict).boundary_violations, 0u);

    ClusterReportOptions loose;
    loose.boundary_threshold = 0.9;  // all 4 cross pairs <= 0.9
    EXPECT_EQ(cluster_report(result, storage, loose).boundary_violations, 4u);
}

TEST(ClusterReportTest, SeparationMetricsNaNForSingleCluster) {
    const DenseStorage storage = MakeTwoClusterStorage();
    const ClusteringResult result = MakeResult({0, 0, 0, 0});

    const ClusterReport r = cluster_report(result, storage, ClusterReportOptions());

    EXPECT_TRUE(std::isnan(r.silhouette));
    EXPECT_TRUE(std::isnan(r.dunn_index));
    // compactness still defined
    EXPECT_FALSE(std::isnan(r.median_radius));
}

TEST(ClusterReportTest, RepresentativeAndCoverage) {
    const DenseStorage storage = MakeTwoClusterStorage();
    const ClusteringResult result = MakeResult({0, 0, 1, 1});

    ClusterReportOptions options;  // coverage {0.25,0.35,0.45}
    const ClusterReport r = cluster_report(result, storage, options);

    EXPECT_DOUBLE_EQ(r.median_medoid_member_distance, 0.2);
    EXPECT_DOUBLE_EQ(r.representative_redundancy, 0.8);
    ASSERT_EQ(r.coverage_thresholds.size(), 3u);
    ASSERT_EQ(r.coverage_at.size(), 3u);
    // Every point is within 0.2 of its cluster medoid, so coverage is 1.0 at
    // all thresholds >= 0.2.
    for (const double c : r.coverage_at) {
        EXPECT_DOUBLE_EQ(c, 1.0);
    }
}

TEST(ClusterReportTest, CoverageCountsNoiseInDenominator) {
    // Cluster {0,1}; points 2,3 noise and far (0.8) from medoid.
    const DenseStorage storage = MakeTwoClusterStorage();
    const ClusteringResult result = MakeResult({0, 0, -1, -1});

    ClusterReportOptions options;
    options.coverage_thresholds = {0.3};
    const ClusterReport r = cluster_report(result, storage, options);

    ASSERT_EQ(r.coverage_at.size(), 1u);
    // Only points 0,1 are within 0.3 of the single medoid; 2,3 (0.8) are not.
    EXPECT_DOUBLE_EQ(r.coverage_at[0], 0.5);
}

TEST(ClusterReportTest, CompareReportsHoldsBothScorecards) {
    const DenseStorage storage = MakeTwoClusterStorage();
    const ClusterReport a = cluster_report(MakeResult({0, 0, 1, 1}), storage, ClusterReportOptions());
    const ClusterReport b = cluster_report(MakeResult({0, 0, 0, 0}), storage, ClusterReportOptions());

    const ClusterReportComparison cmp = compare_reports(a, b);
    EXPECT_EQ(cmp.a.num_clusters, 2u);
    EXPECT_EQ(cmp.b.num_clusters, 1u);
}

TEST(ClusterReportTest, ResultMethodNames) {
    EXPECT_EQ(ClusteringResult().Method(), "");
    EXPECT_EQ(ButinaResult().Method(), "butina");
    EXPECT_EQ(DBSCANResult().Method(), "dbscan");
    EXPECT_EQ(HDBSCANResult().Method(), "hdbscan");
    EXPECT_EQ(AgglomerativeResult().Method(), "agglomerative");
}
