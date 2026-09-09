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

// Six points, clusters {0,1,2} and {3,4,5}. Every cross pair is 0.8, so each
// new index in section 5.3 has a closed form worked out by hand; see
// ClusterReportTest.HandComputedInternalIndices.
//
// Note that Medoid and Minimax pick the SAME representative here -- sample 0
// has both the lowest total (0.4 against 0.6) and the lowest maximum (0.2
// against 0.4). That is why the Minimax baseline below uses the other fixture.
DenseStorage MakeSixPointStorage() {
    DenseStorage storage(6);
    storage.Set(0, 1, 0.2);
    storage.Set(0, 2, 0.2);
    storage.Set(1, 2, 0.4);
    storage.Set(3, 4, 0.2);
    storage.Set(3, 5, 0.2);
    storage.Set(4, 5, 0.4);
    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 3; j < 6; ++j) {
            storage.Set(i, j, 0.8);
        }
    }
    return storage;
}

// Eight points, clusters {0,1,2,3} and {4,5,6,7}, built so the two
// representative methods disagree: Medoid picks 0 and 4, Minimax picks 1 and 5.
// Task 6 reuses it to show the medoid-named indices ignore the method.
DenseStorage MakeDivergentRepresentativeStorage() {
    DenseStorage storage(8);
    // Cluster {0,1,2,3}: totals are 0.9, 1.1, 1.1, 1.7, so the medoid is 0
    // outright. Maxima are 0.7, 0.5, 0.5, 0.7, so minimax is 1 -- it ties with
    // 2 and wins on the earliest-member rule.
    //
    // The winners must differ in TOTAL, not only in maximum. Section 3.7 makes
    // median_medoid_member_distance the one legacy scalar keyed to the
    // configured representative, and it is the mean of that winner's total. If
    // 0 and 1 tied at 1.1 the field would read the same under both methods, and
    // no baseline could catch a later task computing it from the true-medoid
    // vector instead -- exactly the identity mix-up section 3.7 warns about.
    // At 0.7 the means separate: 0.3 for the medoid, 1.1/3 for the minimax.
    storage.Set(0, 1, 0.1);
    storage.Set(0, 2, 0.1);
    storage.Set(0, 3, 0.7);
    storage.Set(1, 2, 0.5);
    storage.Set(1, 3, 0.5);
    storage.Set(2, 3, 0.5);
    // Cluster {4,5,6,7}: the same shape, shifted. Medoid 4, minimax 5. The
    // symmetry is what keeps median_radius and median_medoid_member_distance
    // single-valued rather than an average of two unequal cluster values.
    storage.Set(4, 5, 0.1);
    storage.Set(4, 6, 0.1);
    storage.Set(4, 7, 0.7);
    storage.Set(5, 6, 0.5);
    storage.Set(5, 7, 0.5);
    storage.Set(6, 7, 0.5);
    for (size_t i = 0; i < 4; ++i) {
        for (size_t j = 4; j < 8; ++j) {
            storage.Set(i, j, 0.95);
        }
    }
    // The one asymmetric cross distance, and it is load-bearing. With every
    // cross distance at 0.95, representative_redundancy is the median of
    // {0.95, 0.95} whichever pair of representatives is chosen -- so the
    // baseline would record the same number for both methods and Task 6's
    // EXPECT_NE could never fire. 1 and 5 are the minimax representatives, so
    // only the minimax run sees it.
    storage.Set(1, 5, 0.85);
    return storage;
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

// The bounds check on Get would otherwise answer first, naming the storage
// class rather than the cluster member the caller has to fix. Both messages
// contain "outside the storage range", so the assertion is on the part that
// distinguishes them.
TEST(ClusterReportTest, OutOfRangeClusterMemberNamesTheMemberNotTheBackend) {
    const DenseStorage storage = MakeTwoClusterStorage();
    const ClusteringResult result(
        std::vector<ClusterLabel>{0, 0, 0, 0}, Clusters{{0, 1, 99}});

    try {
        cluster_report(result, storage, ClusterReportOptions());
        FAIL() << "expected an out-of-range refusal";
    } catch (const std::out_of_range& e) {
        EXPECT_STREQ(e.what(), "Cluster member index is outside the storage range");
    }
}

// validate_cluster_members is cluster_report's first *cluster* diagnostic: the
// completeness check exercised by SparseStorageThrows above answers earlier. So
// its wording has to describe the cluster. It used to say "Cluster
// representative requires at least one member", naming an operation -- the same
// misdirection as naming the storage class above.
TEST(ClusterReportTest, EmptyClusterNamesTheClusterNotARepresentative) {
    const DenseStorage storage = MakeTwoClusterStorage();
    const ClusteringResult result(
        std::vector<ClusterLabel>{0, 0, 0, 0}, Clusters{{0, 1}, {}});

    try {
        cluster_report(result, storage, ClusterReportOptions());
        FAIL() << "expected an empty-cluster refusal";
    } catch (const std::invalid_argument& e) {
        EXPECT_STREQ(e.what(), "Cluster must contain at least one member");
    }
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

// Spec section 7.1 item 5. These literals were captured from the 5.0.0 build
// before the passes of section 5.1 were fused. Collapsing three cross-cluster
// walks into one must not move a single existing number, and a test that
// recomputed the expectation with the new code could not tell if it did.
//
// EXPECT_EQ on doubles is intentional: EXPECT_DOUBLE_EQ tolerates four ULPs,
// which is exactly the drift a reassociating fusion would introduce.
TEST(ClusterReportTest, FusionEquivalenceMedoidBaseline) {
    const DenseStorage storage = MakeSixPointStorage();
    const ClusteringResult result = MakeResult({0, 0, 0, 1, 1, 1});
    ClusterReportOptions options;
    options.representative_method = RepresentativeMethod::Medoid;
    const ClusterReport r = cluster_report(result, storage, options);

    EXPECT_EQ(r.num_samples, 6u);
    EXPECT_EQ(r.num_clusters, 2u);
    EXPECT_EQ(r.num_noise, 0u);
    EXPECT_EQ(r.num_singletons, 0u);
    EXPECT_EQ(r.boundary_violations, 0u);
    EXPECT_EQ(r.noise_fraction, 0);
    EXPECT_EQ(r.largest_cluster_fraction, 0.5);
    EXPECT_EQ(r.singleton_fraction, 0);
    EXPECT_EQ(r.cluster_size_median, 3);
    EXPECT_EQ(r.cluster_size_p90, 3);
    EXPECT_EQ(r.size_gini, 0);
    EXPECT_EQ(r.size_entropy, 1);
    EXPECT_EQ(r.mean_intra_distance, 0.26666666666666666);
    EXPECT_EQ(r.median_intra_distance, 0.20000000000000001);
    EXPECT_EQ(r.median_radius, 0.20000000000000001);
    EXPECT_EQ(r.p95_diameter, 0.40000000000000002);
    EXPECT_EQ(r.silhouette, 0.66666666666666663);
    EXPECT_EQ(r.dunn_index, 2);
    EXPECT_EQ(r.median_medoid_member_distance, 0.20000000000000001);
    EXPECT_EQ(r.representative_redundancy, 0.80000000000000004);
    ASSERT_EQ(r.coverage_at.size(), 3u);
    EXPECT_EQ(r.coverage_at[0], 1);
    EXPECT_EQ(r.coverage_at[1], 1);
    EXPECT_EQ(r.coverage_at[2], 1);
}

// The divergent fixture, not MakeSixPointStorage. There Minimax and Medoid
// select the same representatives, so a second baseline over it would pin the
// same twenty numbers twice and leave the minimax path of the fused walk
// untested. Here Minimax picks 1 and 5 where Medoid picks 0 and 4, so
// median_radius, median_medoid_member_distance, representative_redundancy and
// coverage_at all take values the Medoid baseline never sees.
//
// EXPECT_EQ on doubles is intentional: EXPECT_DOUBLE_EQ tolerates four ULPs,
// which is exactly the drift a reassociating fusion would introduce.
TEST(ClusterReportTest, FusionEquivalenceMinimaxBaseline) {
    const DenseStorage storage = MakeDivergentRepresentativeStorage();
    const ClusteringResult result = MakeResult({0, 0, 0, 0, 1, 1, 1, 1});
    ClusterReportOptions options;
    options.representative_method = RepresentativeMethod::Minimax;
    const ClusterReport r = cluster_report(result, storage, options);

    EXPECT_EQ(r.num_samples, 8u);
    EXPECT_EQ(r.num_clusters, 2u);
    EXPECT_EQ(r.num_noise, 0u);
    EXPECT_EQ(r.num_singletons, 0u);
    EXPECT_EQ(r.boundary_violations, 0u);
    EXPECT_EQ(r.noise_fraction, 0);
    EXPECT_EQ(r.largest_cluster_fraction, 0.5);
    EXPECT_EQ(r.singleton_fraction, 0);
    EXPECT_EQ(r.cluster_size_median, 4);
    EXPECT_EQ(r.cluster_size_p90, 4);
    EXPECT_EQ(r.size_gini, 0);
    EXPECT_EQ(r.size_entropy, 1);
    EXPECT_EQ(r.mean_intra_distance, 0.39999999999999997);
    EXPECT_EQ(r.median_intra_distance, 0.5);
    EXPECT_EQ(r.median_radius, 0.5);
    EXPECT_EQ(r.p95_diameter, 0.69999999999999996);
    EXPECT_EQ(r.silhouette, 0.57633949739212897);
    EXPECT_EQ(r.dunn_index, 1.2142857142857144);
    EXPECT_EQ(r.median_medoid_member_distance, 0.3666666666666667);
    EXPECT_EQ(r.representative_redundancy, 0.84999999999999998);
    ASSERT_EQ(r.coverage_at.size(), 3u);
    EXPECT_EQ(r.coverage_at[0], 0.5);
    EXPECT_EQ(r.coverage_at[1], 0.5);
    EXPECT_EQ(r.coverage_at[2], 0.5);
}

// The new surface defaults to "nobody asked for anything", so an existing
// caller who never mentions the flags gets an empty table and a requested
// struct that says so. Section 4.3: the struct records the request, not the
// outcome.
TEST(ClusterReportTest, NewSurfaceDefaultsToUnrequested) {
    const ClusterReportOptions options;
    EXPECT_FALSE(options.compute_pair_rank_indices);
    EXPECT_FALSE(options.compute_per_cluster_records);

    const DenseStorage storage = MakeSixPointStorage();
    const ClusteringResult result = MakeResult({0, 0, 0, 1, 1, 1});
    const ClusterReport r = cluster_report(result, storage, options);

    EXPECT_TRUE(r.records.empty());
    EXPECT_FALSE(r.requested.pair_rank_indices);
    EXPECT_FALSE(r.requested.per_cluster_records);
    EXPECT_EQ(NO_NEAREST_CLUSTER, -1);

    // A default-constructed record carries the section 4.4 initialisers.
    const ClusterRecord record;
    EXPECT_EQ(record.label, 0);
    EXPECT_EQ(record.size, 0u);
    EXPECT_EQ(record.nearest_cluster, NO_NEAREST_CLUSTER);
}
