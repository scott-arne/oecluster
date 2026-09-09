#include <gtest/gtest.h>

#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

#include "oecluster/StorageBackend.h"
#include "oecluster/clustering/Agglomerative.h"
#include "oecluster/clustering/Butina.h"
#include "oecluster/clustering/ClusterReport.h"
#include "oecluster/clustering/ClusterTypes.h"
#include "oecluster/clustering/DBSCAN.h"
#include "oecluster/clustering/HDBSCAN.h"
#include "../../src/clustering/InternalIndices.h"

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
    EXPECT_EQ(record.radius, 0.0);    // 0.0 IS the singleton value here
    EXPECT_EQ(record.diameter, 0.0);  // likewise
    // The four undefined-valued fields default to NaN. The SWIG surface
    // exports ClusterRecord's default constructor, so these defaults are
    // reachable from Python and must not read as measurements.
    EXPECT_TRUE(std::isnan(record.mean_intra_distance));
    EXPECT_TRUE(std::isnan(record.median_intra_distance));
    EXPECT_TRUE(std::isnan(record.nearest_cluster_distance));
    EXPECT_TRUE(std::isnan(record.silhouette));

    // A default-constructed report likewise: the seven new scalars are NaN
    // before anything populates them.
    const ClusterReport blank;
    EXPECT_TRUE(std::isnan(blank.calinski_harabasz_medoid));
    EXPECT_TRUE(std::isnan(blank.davies_bouldin_medoid));
    EXPECT_TRUE(std::isnan(blank.dunn_mean_separation_mean_diameter));
    EXPECT_TRUE(std::isnan(blank.dunn_medoid_separation_medoid_spread));
    EXPECT_TRUE(std::isnan(blank.point_biserial));
    EXPECT_TRUE(std::isnan(blank.c_index));
    EXPECT_TRUE(std::isnan(blank.baker_hubert_gamma));
    // The pre-existing fields keep their 0.0 defaults -- pinned so that a
    // later "consistency" cleanup of the whole struct fails loudly here.
    EXPECT_EQ(blank.silhouette, 0.0);
    EXPECT_EQ(blank.dunn_index, 0.0);
}

// Section 5.2 precondition 1, new invalid_argument refusals. Each message must
// name the offending index: a caller with a hand-built result needs to know
// which sample to fix, not that "something is wrong".
TEST(ClusterReportTest, SampleInTwoClustersIsRefused) {
    const DenseStorage storage = MakeSixPointStorage();
    const ClusteringResult result(
        std::vector<ClusterLabel>{0, 0, 0, 1, 1, 1}, Clusters{{0, 1, 2}, {2, 3, 4}});
    try {
        cluster_report(result, storage, ClusterReportOptions());
        FAIL() << "expected std::invalid_argument";
    } catch (const std::invalid_argument& e) {
        const std::string message(e.what());
        EXPECT_NE(message.find("sample 2"), std::string::npos) << message;
        EXPECT_NE(message.find("appears in clusters"), std::string::npos) << message;
    }
}

TEST(ClusterReportTest, LabelDisagreeingWithClusterOrdinalIsRefused) {
    const DenseStorage storage = MakeSixPointStorage();
    const ClusteringResult result(
        std::vector<ClusterLabel>{0, 0, 1, 1, 1, 1}, Clusters{{0, 1, 2}, {3, 4, 5}});
    try {
        cluster_report(result, storage, ClusterReportOptions());
        FAIL() << "expected std::invalid_argument";
    } catch (const std::invalid_argument& e) {
        const std::string message(e.what());
        EXPECT_NE(message.find("sample 2"), std::string::npos) << message;
        EXPECT_NE(message.find("has label 1 but appears in cluster 0"),
                  std::string::npos)
            << message;
    }
}

TEST(ClusterReportTest, ClusteredSampleOmittedFromEveryMemberListIsRefused) {
    const DenseStorage storage = MakeTwoClusterStorage();
    const ClusteringResult result(
        std::vector<ClusterLabel>{0, 0, -1, -1}, Clusters{{0}});
    try {
        cluster_report(result, storage, ClusterReportOptions());
        FAIL() << "expected std::invalid_argument";
    } catch (const std::invalid_argument& e) {
        const std::string message(e.what());
        EXPECT_NE(message.find("sample 1"), std::string::npos) << message;
        EXPECT_NE(message.find("appears in no cluster"), std::string::npos) << message;
    }
}

// The guard at ClusterReport.cpp:193 used to skip the pre-pass entirely when
// members was empty, so this shape took the empty-clustering branch and
// returned an all-NaN report instead of refusing.
TEST(ClusterReportTest, NonNoiseLabelNamingNoClusterIsRefused) {
    const DenseStorage storage = MakeTwoClusterStorage();
    const ClusteringResult result(
        std::vector<ClusterLabel>{0, -1, -1, -1}, Clusters{});
    try {
        cluster_report(result, storage, ClusterReportOptions());
        FAIL() << "expected std::invalid_argument";
    } catch (const std::invalid_argument& e) {
        const std::string message(e.what());
        EXPECT_NE(message.find("sample 0"), std::string::npos) << message;
        EXPECT_NE(message.find("no cluster"), std::string::npos) << message;
    }
}

TEST(ClusterReportTest, NoiseLabelledSampleInsideAClusterIsRefused) {
    const DenseStorage storage = MakeTwoClusterStorage();
    const ClusteringResult result(
        std::vector<ClusterLabel>{0, -1, 1, 1}, Clusters{{0, 1}, {2, 3}});
    try {
        cluster_report(result, storage, ClusterReportOptions());
        FAIL() << "expected std::invalid_argument";
    } catch (const std::invalid_argument& e) {
        const std::string message(e.what());
        EXPECT_NE(message.find("sample 1"), std::string::npos) << message;
        EXPECT_NE(message.find("labelled noise"), std::string::npos) << message;
    }
}

// A member index can sit inside storage.NumSamples() and past the end of a
// shorter label vector. Without this check the label pass reads out of bounds,
// so the test is the difference between a diagnosis and undefined behaviour.
TEST(ClusterReportTest, MemberBeyondLabelCountIsOutOfRange) {
    const DenseStorage storage = MakeSixPointStorage();
    const ClusteringResult result(
        std::vector<ClusterLabel>{0, 0}, Clusters{{0, 1, 4}});
    try {
        cluster_report(result, storage, ClusterReportOptions());
        FAIL() << "expected std::out_of_range";
    } catch (const std::out_of_range& e) {
        const std::string message(e.what());
        // 4 is inside storage (six samples) and past the two-entry label
        // vector, so naming the member is the whole point of the message.
        EXPECT_NE(message.find("member 4"), std::string::npos) << message;
        EXPECT_NE(message.find("label count 2"), std::string::npos) << message;
    }

    // The boundary itself. member == labels.size() is the first index that
    // reads past the end, so `>=` rather than `>` is the whole guard; without
    // this arm a weakened comparison stays green while owner[2] runs off a
    // two-element vector.
    const ClusteringResult at_boundary(
        std::vector<ClusterLabel>{0, 0}, Clusters{{0, 1, 2}});
    try {
        cluster_report(at_boundary, storage, ClusterReportOptions());
        FAIL() << "expected std::out_of_range";
    } catch (const std::out_of_range& e) {
        const std::string message(e.what());
        EXPECT_NE(message.find("member 2"), std::string::npos) << message;
        EXPECT_NE(message.find("label count 2"), std::string::npos) << message;
    }
}

// Section 5.2. The header promises out_of_range when a result that has clusters
// labels more samples than storage holds. Before the bijection pre-pass the
// backend produced that type incidentally, from deep inside the coverage loop;
// the pre-pass now reaches the surplus sample first, so without an explicit
// check the caller is told to fix a cluster list when the real error is a
// mismatched storage. INVARIANT 1.
TEST(ClusterReportTest, LabelsLongerThanStorageIsOutOfRange) {
    const DenseStorage storage = MakeTwoClusterStorage();  // four samples

    const ClusteringResult surplus_clustered(
        std::vector<ClusterLabel>{0, 0, 1, 1, 0, 1}, Clusters{{0, 1}, {2, 3}});
    try {
        cluster_report(surplus_clustered, storage, ClusterReportOptions());
        FAIL() << "expected std::out_of_range";
    } catch (const std::out_of_range& e) {
        const std::string message(e.what());
        EXPECT_NE(message.find("label count 6"), std::string::npos) << message;
        EXPECT_NE(message.find("storage sample count 4"), std::string::npos)
            << message;
    }

    // Surplus labelled noise: the one subset for which the documented type
    // survived on its own. Asserting the type alone would also pass on the old
    // downstream DenseStorage refusal, so assert the message that only the new
    // check can produce.
    const ClusteringResult surplus_noise(
        std::vector<ClusterLabel>{0, 0, 1, 1, -1, -1}, Clusters{{0, 1}, {2, 3}});
    try {
        cluster_report(surplus_noise, storage, ClusterReportOptions());
        FAIL() << "expected std::out_of_range";
    } catch (const std::out_of_range& e) {
        const std::string message(e.what());
        EXPECT_NE(message.find("label count 6"), std::string::npos) << message;
        EXPECT_NE(message.find("storage sample count 4"), std::string::npos)
            << message;
    }

    // With no clusters the header's clause does not apply, and a long all-noise
    // label vector stays accepted. INVARIANT 3: this half is the over-refusal
    // guard on the new check.
    const ClusteringResult no_clusters(
        std::vector<ClusterLabel>{-1, -1, -1, -1, -1, -1}, Clusters{});
    EXPECT_NO_THROW(cluster_report(no_clusters, storage, ClusterReportOptions()));

    // Both invalid at once. INVARIANT 1: the caller paired a result with the
    // wrong storage, and the out-of-range member is a symptom of that, so the
    // mismatch has to be named first. This arm is what makes the ordering of
    // the guard against validate_cluster_members a tested property rather than
    // a comment.
    const ClusteringResult mismatched_and_bad_member(
        std::vector<ClusterLabel>{0, 0, 1, 1, 0, 1},
        Clusters{{0, 1}, {2, 99}});
    try {
        cluster_report(mismatched_and_bad_member, storage, ClusterReportOptions());
        FAIL() << "expected std::out_of_range";
    } catch (const std::out_of_range& e) {
        const std::string message(e.what());
        EXPECT_NE(message.find("label count 6"), std::string::npos) << message;
        EXPECT_NE(message.find("storage sample count 4"), std::string::npos)
            << message;
    }
}

// The extended pre-pass must not reclassify the three refusals
// validate_cluster_members already owns.
TEST(ClusterReportTest, ExistingRefusalsKeepTheirTypes) {
    const DenseStorage storage = MakeTwoClusterStorage();

    const ClusteringResult empty_cluster(
        std::vector<ClusterLabel>{-1, -1, -1, -1}, Clusters{{}});
    try {
        cluster_report(empty_cluster, storage, ClusterReportOptions());
        FAIL() << "expected std::invalid_argument";
    } catch (const std::invalid_argument& e) {
        EXPECT_NE(std::string(e.what()).find("Cluster must contain at least one member"),
                  std::string::npos)
            << e.what();
    }

    const ClusteringResult duplicate(
        std::vector<ClusterLabel>{0, -1, -1, -1}, Clusters{{0, 0}});
    try {
        cluster_report(duplicate, storage, ClusterReportOptions());
        FAIL() << "expected std::invalid_argument";
    } catch (const std::invalid_argument& e) {
        EXPECT_NE(std::string(e.what()).find("Cluster members must be unique"),
                  std::string::npos)
            << e.what();
    }

    const ClusteringResult beyond_storage(
        std::vector<ClusterLabel>{0, 0, 0, 0}, Clusters{{0, 1, 99}});
    try {
        cluster_report(beyond_storage, storage, ClusterReportOptions());
        FAIL() << "expected std::out_of_range";
    } catch (const std::out_of_range& e) {
        EXPECT_NE(
            std::string(e.what()).find("Cluster member index is outside the storage range"),
            std::string::npos)
            << e.what();
    }
}

// When a result is invalid in two different LAYERS at once -- a malformed
// cluster and a partition that double-counts a sample -- which refusal fires is
// a property of the errors, not of which cluster happens to hold them.
// INVARIANT 1: name the reason the caller must fix first, and name the same one
// whichever order the clusters arrive in. Within a single layer the selection is
// deliberately not canonicalised; see
// StructuralErrorSelectionIsNotCanonicalised below.
TEST(ClusterReportTest, CompoundInvalidResultsRefuseInAStableOrder) {
    const DenseStorage storage = MakeTwoClusterStorage();  // four samples

    // A malformed cluster outranks a partition that double-counts a sample, and
    // it does so whichever cluster holds which. Each pair below is the same two
    // errors with the cluster order swapped; both halves must answer the same
    // way. This is the boundary the staged passes exist to guarantee -- see the
    // comment in ClusterReport.cpp for what is deliberately NOT guaranteed
    // inside a single cluster.
    const std::vector<ClusterLabel> one_clustered{0, -1, -1, -1};
    const struct {
        const char* what;
        Clusters clusters;
    } structural_beats_ownership[] = {
        {"Cluster must contain at least one member", Clusters{{0}, {0}, {}}},
        {"Cluster must contain at least one member", Clusters{{}, {0}, {0}}},
        {"Cluster members must be unique", Clusters{{0}, {0}, {2, 2}}},
        {"Cluster members must be unique", Clusters{{2, 2}, {0}, {0}}},
    };
    for (const auto& arm : structural_beats_ownership) {
        const ClusteringResult result(one_clustered, arm.clusters);
        try {
            cluster_report(result, storage, ClusterReportOptions());
            FAIL() << "expected std::invalid_argument for " << arm.what;
        } catch (const std::invalid_argument& e) {
            EXPECT_NE(std::string(e.what()).find(arm.what), std::string::npos)
                << e.what();
        }
    }

    // The storage-range half of the same boundary. Its type is out_of_range, so
    // it cannot share the loop above.
    for (const Clusters& clusters :
         {Clusters{{0}, {0}, {99}}, Clusters{{99}, {0}, {0}}}) {
        const ClusteringResult result(one_clustered, clusters);
        try {
            cluster_report(result, storage, ClusterReportOptions());
            FAIL() << "expected std::out_of_range";
        } catch (const std::out_of_range& e) {
            EXPECT_NE(std::string(e.what())
                          .find("Cluster member index is outside the storage range"),
                      std::string::npos)
                << e.what();
        }
    }

    // A member past the label vector beats a duplicate regardless of which
    // cluster holds which. Before the passes were staged these two shapes
    // disagreed with each other.
    const ClusteringResult bad_index_first(
        std::vector<ClusterLabel>{0, 0}, Clusters{{0, 3}, {0}});
    EXPECT_THROW(cluster_report(bad_index_first, storage, ClusterReportOptions()),
                 std::out_of_range);

    const ClusteringResult duplicate_first(
        std::vector<ClusterLabel>{0, 0}, Clusters{{0}, {0, 3}});
    EXPECT_THROW(cluster_report(duplicate_first, storage, ClusterReportOptions()),
                 std::out_of_range);
}

// The other half of the round-3 decision, kept visible. The structural pass
// short-circuits twice over: it stops at the first malformed cluster, and the
// shared validator stops at the first bad member inside that cluster. So both
// the order of the clusters and the order of the members within one change
// which refusal a caller sees, and the intra-cluster pair changes the exception
// type too. Declined rather than fixed: ordering these from cluster_report
// would mean duplicating checks that belong in DistanceAccess.h, and no message
// misleads the caller, who has a malformed cluster list either way. This test
// exists so the asymmetry is a recorded decision rather than a surprise.
//
// Both pairs flip under canonicalisation, but not through the same arm. The
// across-cluster pair flips if pass 1 stops halting at the first malformed
// cluster. The intra-cluster pair flips whichever way the member scan is split
// into complete passes -- uniqueness first reddens {99, 0, 0}, range first
// reddens {0, 0, 99} -- so neither arm alone covers it. Note that merely
// swapping the two checks inside the validator's existing single loop is not a
// canonicalisation and changes no answer: the first member to fault has by
// definition no earlier occurrence, so no member reaches either check both out
// of range and already seen.
TEST(ClusterReportTest, StructuralErrorSelectionIsNotCanonicalised) {
    const DenseStorage storage = MakeTwoClusterStorage();  // four samples
    const std::vector<ClusterLabel> one_clustered{0, -1, -1, -1};

    // Same two errors, one per cluster, clusters swapped.
    const struct {
        const char* what;
        Clusters clusters;
    } across_clusters[] = {
        {"Cluster members must be unique", Clusters{{0, 0}, {}}},
        {"Cluster must contain at least one member", Clusters{{}, {0, 0}}},
    };
    for (const auto& arm : across_clusters) {
        const ClusteringResult result(one_clustered, arm.clusters);
        try {
            cluster_report(result, storage, ClusterReportOptions());
            FAIL() << "expected std::invalid_argument for " << arm.what;
        } catch (const std::invalid_argument& e) {
            EXPECT_NE(std::string(e.what()).find(arm.what), std::string::npos)
                << e.what();
        }
    }

    // Same two errors inside ONE cluster, members swapped. This is the pair the
    // across-cluster arms above cannot see: it turns on the validator's own
    // loop order, and it crosses exception types.
    const ClusteringResult duplicate_before_out_of_range(
        one_clustered, Clusters{{0, 0, 99}});
    try {
        cluster_report(duplicate_before_out_of_range, storage, ClusterReportOptions());
        FAIL() << "expected std::invalid_argument";
    } catch (const std::invalid_argument& e) {
        EXPECT_NE(std::string(e.what()).find("Cluster members must be unique"),
                  std::string::npos)
            << e.what();
    }

    const ClusteringResult out_of_range_before_duplicate(
        one_clustered, Clusters{{99, 0, 0}});
    try {
        cluster_report(out_of_range_before_duplicate, storage, ClusterReportOptions());
        FAIL() << "expected std::out_of_range";
    } catch (const std::out_of_range& e) {
        EXPECT_NE(std::string(e.what())
                      .find("Cluster member index is outside the storage range"),
                  std::string::npos)
            << e.what();
    }
}

// Section 5.5 defines both of these as answerable, so the bijection check must
// not turn them into refusals. INVARIANT 3: over-refusal is as severe as the
// wrong number it was meant to prevent.
TEST(ClusterReportTest, DegenerateShapesAreStillAccepted) {
    const DenseStorage storage = MakeTwoClusterStorage();
    const ClusteringResult all_noise(
        std::vector<ClusterLabel>{-1, -1, -1, -1}, Clusters{});
    EXPECT_NO_THROW(cluster_report(all_noise, storage, ClusterReportOptions()));

    const DenseStorage empty_storage(0);
    const ClusteringResult zero_samples(std::vector<ClusterLabel>{}, Clusters{});
    EXPECT_NO_THROW(
        cluster_report(zero_samples, empty_storage, ClusterReportOptions()));
}

// Section 7.1 item 9. A new refusal that rejects the library's own output is a
// worse defect than the one it prevents. Options are set explicitly rather than
// defaulted: on six points the default min_samples of 5 and Butina's default
// threshold of 0.0 both degenerate, and a shape with no clusters would not
// exercise the check.
TEST(ClusterReportTest, PreconditionsAcceptEveryShippedAlgorithm) {
    const DenseStorage storage = MakeSixPointStorage();
    const ClusterReportOptions options;

    ButinaOptions butina_options;
    butina_options.distance_threshold = 0.5;
    const auto butina = butina_cluster(storage, butina_options);
    ASSERT_GE(butina.NumClusters(), 2u);
    EXPECT_NO_THROW(cluster_report(butina, storage, options));

    AgglomerativeOptions agglomerative_options;
    agglomerative_options.n_clusters = 2;
    const auto agglomerative = agglomerative_cluster(storage, agglomerative_options);
    ASSERT_GE(agglomerative.NumClusters(), 2u);
    EXPECT_NO_THROW(cluster_report(agglomerative, storage, options));

    DBSCANOptions dbscan_options;
    dbscan_options.eps = 0.5;
    dbscan_options.min_samples = 2;
    const auto dbscan = dbscan_cluster(storage, dbscan_options);
    ASSERT_GE(dbscan.NumClusters(), 2u);
    EXPECT_NO_THROW(cluster_report(dbscan, storage, options));

    HDBSCANOptions hdbscan_options;
    hdbscan_options.min_cluster_size = 2;
    const auto hdbscan = hdbscan_cluster(storage, hdbscan_options);
    ASSERT_GE(hdbscan.NumClusters(), 2u);
    EXPECT_NO_THROW(cluster_report(hdbscan, storage, options));
}

// Section 7.1 item 2. The tie-aware merge is the riskiest code in A1, so it is
// checked against the definition it implements rather than against itself.
namespace {

detail::PairRankIndices BruteForcePairRank(
    std::vector<double> within,
    std::vector<double> between) {
    unsigned long long s_plus = 0;
    unsigned long long s_minus = 0;
    for (const double w : within) {
        for (const double b : between) {
            if (w < b) {
                ++s_plus;
            } else if (b < w) {
                ++s_minus;
            }
        }
    }
    detail::PairRankIndices out{
        std::numeric_limits<double>::quiet_NaN(),
        std::numeric_limits<double>::quiet_NaN()};
    const double denominator =
        static_cast<double>(s_plus) + static_cast<double>(s_minus);
    if (denominator > 0.0) {
        out.baker_hubert_gamma =
            (static_cast<double>(s_plus) - static_cast<double>(s_minus)) / denominator;
    }

    // C-index by direct summation of the w smallest and w largest of all P.
    std::vector<double> all = within;
    all.insert(all.end(), between.begin(), between.end());
    std::sort(all.begin(), all.end());
    const size_t w_count = within.size();
    if (w_count > 0) {
        double s_w = 0.0;
        for (const double d : within) {
            s_w += d;
        }
        double s_min = 0.0;
        double s_max = 0.0;
        for (size_t i = 0; i < w_count; ++i) {
            s_min += all[i];
            s_max += all[all.size() - 1 - i];
        }
        if (s_max != s_min) {
            out.c_index = (s_w - s_min) / (s_max - s_min);
        }
    }
    return out;
}

}  // namespace

// Section 7.1 items 2 and 4. The brute-force reference is the whole point: it
// recomputes gamma by the O(P_w * P_b) definition and c_index by explicitly
// sorting all pooled pairs, so it shares no code path with the run-at-a-time
// walk it checks. Item 4 asks for exactly this on c_index.
TEST(InternalIndicesTest, PairRankMatchesBruteForce) {
    const std::vector<double> within{0.2, 0.35, 0.5, 0.1};
    const std::vector<double> between{0.6, 0.15, 0.9, 0.45, 0.7};
    const detail::PairRankIndices got = detail::pair_rank_indices(within, between);
    const detail::PairRankIndices want = BruteForcePairRank(within, between);
    // Gamma is a ratio of two exactly-represented integer counts, so both
    // routes land on the identical double and exact equality is honest.
    EXPECT_DOUBLE_EQ(got.baker_hubert_gamma, want.baker_hubert_gamma);
    // C-index is not, and must not be compared that way. The reference sums
    // `within` in its original order while pair_rank_indices sums it after
    // sorting; on this fixture the two sums differ by one ULP, which the
    // division amplifies to five -- past EXPECT_DOUBLE_EQ's four-ULP budget.
    // Measured, not estimated: 0.18421052631578951897 against
    // 0.18421052631578938019 under clang -O2. The tolerance is absolute and
    // still four orders of magnitude tighter than any real disagreement the
    // merge-walk could produce, since taking a wrong element shifts S_min or
    // S_max by at least 0.05.
    EXPECT_NEAR(got.c_index, want.c_index, 1e-12);
}

// Section 7.1 item 3. A couple whose two distances are equal contributes to
// neither counter; the run-at-a-time walk is what makes that true.
TEST(InternalIndicesTest, PairRankExcludesTiedCouples) {
    const std::vector<double> within{0.5};
    const std::vector<double> between{0.5, 0.9};
    const detail::PairRankIndices got = detail::pair_rank_indices(within, between);
    // One concordant couple (0.5 < 0.9), one tie, no discordant couples.
    EXPECT_DOUBLE_EQ(got.baker_hubert_gamma, 1.0);
    EXPECT_DOUBLE_EQ(got.baker_hubert_gamma,
                     BruteForcePairRank(within, between).baker_hubert_gamma);
}

// Section 7.1 item 16. An unsigned subtraction wraps here and returns a value
// near +1, scoring the worst clustering as the best one.
TEST(InternalIndicesTest, PairRankGammaGoesNegative) {
    const std::vector<double> within{0.9, 0.9};
    const std::vector<double> between{0.1, 0.1, 0.1, 0.1};
    const detail::PairRankIndices got = detail::pair_rank_indices(within, between);
    EXPECT_LT(got.baker_hubert_gamma, 0.0);
    EXPECT_DOUBLE_EQ(got.baker_hubert_gamma, -1.0);
    EXPECT_DOUBLE_EQ(got.baker_hubert_gamma,
                     BruteForcePairRank(within, between).baker_hubert_gamma);
}

TEST(InternalIndicesTest, PairRankAllTiedIsNaNNotARefusal) {
    const std::vector<double> within{0.5, 0.5};
    const std::vector<double> between{0.5, 0.5, 0.5, 0.5};
    detail::PairRankIndices got{0.0, 0.0};
    EXPECT_NO_THROW(got = detail::pair_rank_indices(within, between));
    EXPECT_TRUE(std::isnan(got.baker_hubert_gamma));
}

// Section 7.1 item 17. The refusal's real-scale trigger needs about 69 GB and
// is untestable; the arithmetic that decides it is not, which is why the guard
// is a named function rather than an inline condition.
TEST(InternalIndicesTest, AddCouplesRefusesRatherThanWrapping) {
    constexpr unsigned long long MAXIMUM =
        std::numeric_limits<unsigned long long>::max();

    unsigned long long counter = MAXIMUM - 5;
    EXPECT_NO_THROW(detail::add_couples(counter, 5));
    EXPECT_EQ(counter, MAXIMUM);

    EXPECT_THROW(detail::add_couples(counter, 1), std::length_error);
    EXPECT_EQ(counter, MAXIMUM);

    EXPECT_NO_THROW(detail::add_couples(counter, 0));
    EXPECT_EQ(counter, MAXIMUM);
}

// Section 5.3. The naive sum-of-squares form cancels on near-constant
// fingerprint distances; Welford does not.
TEST(InternalIndicesTest, WelfordMergeMatchesSinglePass) {
    const std::vector<double> left{0.2, 0.2, 0.4};
    const std::vector<double> right{0.8, 0.8, 0.8, 0.8};

    detail::DistanceMoments a;
    for (const double d : left) {
        a.Add(d);
    }
    detail::DistanceMoments b;
    for (const double d : right) {
        b.Add(d);
    }
    const detail::DistanceMoments merged = detail::merge_moments(a, b);

    detail::DistanceMoments single;
    for (const double d : left) {
        single.Add(d);
    }
    for (const double d : right) {
        single.Add(d);
    }

    EXPECT_EQ(merged.count, single.count);
    EXPECT_NEAR(merged.mean, single.mean, 1e-15);
    EXPECT_NEAR(detail::population_stddev(merged),
                detail::population_stddev(single), 1e-15);
}

TEST(InternalIndicesTest, PopulationStddevIsZeroBelowTwoValues) {
    const detail::DistanceMoments empty;
    EXPECT_DOUBLE_EQ(detail::population_stddev(empty), 0.0);

    detail::DistanceMoments one;
    one.Add(0.42);
    EXPECT_DOUBLE_EQ(detail::population_stddev(one), 0.0);
}
