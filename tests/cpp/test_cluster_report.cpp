#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <limits>
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
#include "../../src/clustering/ClusterMetrics.h"
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

// Nine samples in three clusters of unequal size -- A = {0,1}, B = {2,3,4},
// C = {5,6,7,8}. Most partitions in this file are K <= 2, which leaves a report
// that aggregates only the first cluster indistinguishable from one that
// aggregates all of them, and a cross pass that stops after pair (0,1)
// indistinguishable from one that walks all three pairs. This is the fixture
// built to separate those: three clusters, three sizes and three separations,
// no two of them shared. Other fixtures do reach K >= 3 -- the singleton ones
// below, and MakeSixPointStorage under {0,0,1,1,2,-1} or an all-singleton
// partition -- but every one of them repeats a cluster size, so a size-weighted
// aggregate and an unweighted one can still coincide there. Reaching K >= 3 is
// not the same as separating what K >= 3 makes separable.
//
// Three properties are load-bearing and none of them is incidental:
//
//   * The three intra multisets differ, so dropping any one cluster moves the
//     intra median as well as the mean.
//   * Each cluster pair has its own separation -- A-B 0.80, A-C 0.90,
//     B-C 0.60 -- so the Dunn numerator has a unique minimum and the third
//     pair is not a repeat of the first two. Within a pair the separation is
//     flat, which makes every cross mean equal to the separation itself and
//     the silhouette b terms exact by inspection.
//   * The clusters have distinct radii (0.32, 0.24, 0.08) and distinct
//     medoid-member means (0.32, 0.18, 0.06), so neither median sits on the
//     first cluster's value.
//
// Medoid representatives, by lowest mean distance to the rest of the cluster:
//   A: both members are at 0.32, and the earliest-member tie rule picks 0.
//   B: 0.18, 0.24, 0.30      -> 2
//   C: 0.06, 0.18, 0.30, 0.42 -> 5
DenseStorage MakeThreeClusterStorage() {
    DenseStorage storage(9);
    // Cluster A = {0,1}: a single intra pair.
    storage.Set(0, 1, 0.32);
    // Cluster B = {2,3,4}: member means 0.18, 0.24, 0.30.
    storage.Set(2, 3, 0.12);
    storage.Set(2, 4, 0.24);
    storage.Set(3, 4, 0.36);
    // Cluster C = {5,6,7,8}: member means 0.06, 0.18, 0.30, 0.42. Sample 5 sits
    // very close to all three others while 7 and 8 are far apart, so C is the
    // only cluster whose radius and diameter diverge sharply.
    storage.Set(5, 6, 0.04);
    storage.Set(5, 7, 0.06);
    storage.Set(5, 8, 0.08);
    storage.Set(6, 7, 0.08);
    storage.Set(6, 8, 0.42);
    storage.Set(7, 8, 0.76);
    for (size_t i = 0; i < 2; ++i) {
        for (size_t j = 2; j < 5; ++j) {
            storage.Set(i, j, 0.80);
        }
        for (size_t j = 5; j < 9; ++j) {
            storage.Set(i, j, 0.90);
        }
    }
    for (size_t i = 2; i < 5; ++i) {
        for (size_t j = 5; j < 9; ++j) {
            storage.Set(i, j, 0.60);
        }
    }
    return storage;
}

// Four points, clusters {0,1} and {2,3}, with the cross block deliberately
// lopsided: sample 1 sits close to both members of the other cluster while
// sample 0 sits far from them. Every other multi-cluster fixture in this file
// has a flat cross block, and that flatness is what hides the two distinctions
// below.
//
//   * Both clusters tie internally, so the earliest-member rule elects medoids
//     0 and 2, while the smallest total over the clustered set belongs to
//     sample 1. The global medoid M is therefore a point that no cluster
//     elected -- the case a search restricted to the elected medoids misses.
//   * The smallest medoid-to-medoid distance is d(0,2) = 7/8, but the smallest
//     cross-cluster MEMBER distance is d(1,2) = 1/8. The two Dunn numerators
//     are a factor of seven apart rather than coincidentally equal.
//
// Every distance is a dyadic fraction, so each quantity derived from them is
// exact in binary floating point and can be asserted with EXPECT_DOUBLE_EQ.
DenseStorage MakeOffMedoidGlobalStorage() {
    DenseStorage storage(4);
    storage.Set(0, 1, 0.5);
    storage.Set(2, 3, 0.5);
    storage.Set(0, 2, 0.875);
    storage.Set(0, 3, 0.875);
    storage.Set(1, 2, 0.125);
    storage.Set(1, 3, 0.125);
    return storage;
}

// Three singleton clusters, one per sample. A one-member cluster is its own
// representative under every method, so this fixture removes representative
// selection from the question entirely and leaves only the reduction that
// turns K representatives into representative_redundancy.
//
// All three pairwise distances differ, and that is load-bearing. Almost every
// other partition this file feeds to the redundancy reduction has exactly two
// clusters, where the row-minima vector is [d, d] and every candidate
// reduction returns d. Three distinct distances are what separate the loop
// bounds and the median-versus-mean reduction.
//
// They do not separate everything. At K = 3 the smallest distance is a row
// minimum for two of the three rows, so the median of the row minima is always
// their minimum and always equals the median of their prefix minima -- neither
// a min-instead-of-median reduction nor a running minimum hoisted above the
// outer loop can be seen here, under any label permutation. That is what
// MakeFourSingletonStorage below is for.
//
// Which number each candidate lands on depends on the label order as well as
// on the distances, because the reduction walks the representatives in cluster
// ordinal order. The worked arithmetic therefore lives with the test that fixes
// that order -- see
// NonMedoidRepresentativeRedundancyIsTheMedianOfThreeRowMinima, which permutes
// the labels deliberately.
DenseStorage MakeThreeSingletonStorage() {
    DenseStorage storage(3);
    storage.Set(0, 1, 0.2);
    storage.Set(0, 2, 0.9);
    storage.Set(1, 2, 0.8);
    return storage;
}

// Four singleton clusters, and the companion the fixture above cannot be
// rewritten into. At K = 3 the smallest of the three distances is a row minimum
// for two of the three rows, so the row-minima vector is always [m, m, d] with
// m the global minimum: its median IS its minimum, and its prefix minima have
// that same median. No label permutation changes either fact, so K = 3 cannot
// tell the shipped median of the row minima from their minimum, nor from a
// running minimum carried across rows. K = 4 is the smallest partition that
// can.
//
// A one-member cluster is its own representative under every method, so this
// one storage drives both the Medoid and the non-Medoid reduction -- two
// separate loops over what are here the same points.
//
// With MakeResult({0, 1, 2, 3}) the clusters come out in storage-index order,
// so the representatives are [0, 1, 2, 3] and the row minima are
//
//   row 0: min(0.9, 0.7, 0.2) = 0.2
//   row 1: min(0.9, 0.6, 0.8) = 0.6
//   row 2: min(0.7, 0.6, 0.5) = 0.5
//   row 3: min(0.2, 0.8, 0.5) = 0.2
//
// giving [0.2, 0.6, 0.5, 0.2], sorted [0.2, 0.2, 0.5, 0.6], median
// (0.2 + 0.5) / 2 = 0.35. That quotient is exact in IEEE double, so the
// assertions need no tolerance.
//
// Which index holds the global minimum is load-bearing, the way the label order
// is load-bearing for MakeThreeSingletonStorage. The smallest distance in the
// matrix is d(0, 3) = 0.2 and ROW 0 reaches it, so a running minimum hoisted
// above the outer loop collapses the vector to its prefix minima,
// [0.2, 0.2, 0.2, 0.2]. Any assignment whose row minima come out non-increasing
// is its own prefix-min vector and would leave that hoist green -- reordering
// the matrix disarms the test without reddening anything. Re-derive every
// candidate, the prefix minima included, before changing a distance; the worked
// table lives with the tests that fix the order, see
// NonMedoidRepresentativeRedundancyIsTheMedianOfFourRowMinima.
DenseStorage MakeFourSingletonStorage() {
    DenseStorage storage(4);
    storage.Set(0, 1, 0.9);
    storage.Set(0, 2, 0.7);
    storage.Set(0, 3, 0.2);
    storage.Set(1, 2, 0.6);
    storage.Set(1, 3, 0.8);
    storage.Set(2, 3, 0.5);
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

// The b term's divisor is the size of the OTHER cluster, and the fused cross
// pass has to supply it by hand -- `size_b` for a point in cluster a and
// `size_a` for a point in cluster b. Every other fixture that asserts
// silhouette has equal-sized clusters (2+2, 3+3, 4+4), which makes exchanging
// the two divisors arithmetically invisible rather than merely undetected.
// Here the clusters are 3 and 2, and the exchange moves the result to
// 0.63333333333333341.
//
// Points 0, 3 and 4 score 0.75 and points 1 and 2 score 0.625, so the mean
// over the five clustered points is 3.5 / 5. Sample 5 is noise and is not
// scored. EXPECT_NEAR rather than an exact literal because the value is
// hand-derived and each per-point term is a division: the accumulated rounding
// is around 1e-16, while the defect this pins is 0.067 away.
TEST(ClusterReportTest, SilhouetteBTermDividesByTheOtherClusterSize) {
    const DenseStorage storage = MakeSixPointStorage();
    const ClusteringResult result = MakeResult({0, 0, 0, 1, 1, -1});

    const ClusterReport r = cluster_report(result, storage, ClusterReportOptions());

    EXPECT_EQ(r.num_clusters, 2u);
    EXPECT_EQ(r.num_noise, 1u);
    EXPECT_NEAR(r.silhouette, 0.7, 1e-12);
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

    // The two cases above bracket 0.8 without ever landing on it, so neither
    // can tell `distance <= threshold` from `distance <`. Every threshold
    // elsewhere in this file brackets its distances the same way. A threshold
    // set to the cross distance exactly is the only input that separates the
    // two, and it is exact: storage.Set writes the literal 0.8 and this reads
    // the same literal, so both sides are bit-identical doubles regardless of
    // 0.8 being inexact in binary.
    ClusterReportOptions exact;
    exact.boundary_threshold = 0.8;
    exact.compute_per_cluster_records = true;
    const ClusterReport on_boundary = cluster_report(result, storage, exact);
    EXPECT_EQ(on_boundary.boundary_violations, 4u);

    // The per-record counts too: the single cluster pair contributes all four
    // to both records, so a strict comparison zeroes all three numbers.
    ASSERT_EQ(on_boundary.records.size(), 2u);
    EXPECT_EQ(on_boundary.records[0].boundary_violations, 4u);
    EXPECT_EQ(on_boundary.records[1].boundary_violations, 4u);
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

// Every distance-derived legacy scalar over three unequal clusters, hand
// derived from MakeThreeClusterStorage. This is the only assertion in the file
// that can tell whole-partition aggregation apart from first-cluster-only
// aggregation, or a complete cross pass apart from one that stops after the
// first cluster pair.
TEST(ClusterReportTest, HandComputedThreeClusterReport) {
    const DenseStorage storage = MakeThreeClusterStorage();
    const ClusteringResult result = MakeResult({0, 0, 1, 1, 1, 2, 2, 2, 2});

    ClusterReportOptions options;
    // 0.80 and 0.60 fall inside, 0.90 does not, so the count spans two of the
    // three cluster pairs rather than resting on a single one.
    options.boundary_threshold = 0.85;
    const ClusterReport r = cluster_report(result, storage, options);

    // The ten intra pairs are 0.32 | 0.12 0.24 0.36 | 0.04 0.06 0.08 0.08 0.42
    // 0.76, summing to 2.48.
    EXPECT_DOUBLE_EQ(r.mean_intra_distance, 0.248);
    // Sorted they are 0.04 0.06 0.08 0.08 0.12 0.24 0.32 0.36 0.42 0.76, and
    // the even-count median averages the fifth and sixth: (0.12 + 0.24) / 2.
    EXPECT_DOUBLE_EQ(r.median_intra_distance, 0.18);
    // Radii, each the medoid's farthest member: 0.32, max(0.12, 0.24) = 0.24,
    // max(0.04, 0.06, 0.08) = 0.08. The median of three is the middle one.
    EXPECT_DOUBLE_EQ(r.median_radius, 0.24);
    // Diameters 0.32, 0.36, 0.76. Fractional-rank p95 over three values lands
    // at rank 0.95 * 2 = 1.9, i.e. 0.36 + 0.9 * (0.76 - 0.36).
    EXPECT_DOUBLE_EQ(r.p95_diameter, 0.72);
    // Each cross group is flat, so b is 0.80 for A and 0.60 for B and C, and
    // every a term is the member mean listed on the fixture. Per point:
    //   A: (0.80 - 0.32) / 0.80 = 0.60, twice                       -> 1.20
    //   B: 0.70, 0.60, 0.50 from a = 0.18, 0.24, 0.30               -> 1.80
    //   C: 0.90, 0.70, 0.50, 0.30 from a = 0.06, 0.18, 0.30, 0.42   -> 2.40
    // 5.40 over nine points.
    EXPECT_DOUBLE_EQ(r.silhouette, 0.6);
    // Closest cross pair is B-C at 0.60; largest diameter is C's 0.76.
    EXPECT_DOUBLE_EQ(r.dunn_index, 0.6 / 0.76);
    // A-B contributes 2 * 3 = 6 pairs at 0.80 and B-C contributes 3 * 4 = 12 at
    // 0.60; A-C's 0.90 is outside the threshold.
    EXPECT_EQ(r.boundary_violations, 18u);
    // Medoid-to-member means: 0.32, (0.12 + 0.24) / 2, (0.04 + 0.06 + 0.08) / 3.
    EXPECT_DOUBLE_EQ(r.median_medoid_member_distance, 0.18);
    // Representatives 0, 2 and 5, whose nearest-other distances are 0.80, 0.60
    // and 0.60. Dropping cluster C would leave only 0 and 2 and read 0.80.
    EXPECT_DOUBLE_EQ(r.representative_redundancy, 0.6);
    // Sizes 2, 3 and 4. The fractional rank for the median is 0.5 * 2 = 1,
    // which lands on the middle size outright; for p90 it is 0.9 * 2 = 1.8,
    // i.e. 3 + 0.8 * (4 - 3). This is the only fixture in the file that
    // asserts either field on unequal sizes -- everywhere else the clusters
    // are the same size and the two percentiles, the minimum and the maximum
    // all coincide.
    EXPECT_DOUBLE_EQ(r.cluster_size_median, 3.0);
    EXPECT_DOUBLE_EQ(r.cluster_size_p90, 3.8);
}

// The same fixture at a threshold wide enough to admit every cross pair, which
// is the only assertion here that can see the A-C iteration at all.
TEST(ClusterReportTest, BoundaryCountIncludesTheWidestClusterPair) {
    const DenseStorage storage = MakeThreeClusterStorage();
    const ClusteringResult result = MakeResult({0, 0, 1, 1, 1, 2, 2, 2, 2});

    ClusterReportOptions options;
    options.boundary_threshold = 0.95;
    const ClusterReport r = cluster_report(result, storage, options);

    // All 26 cross-cluster pairs: A-B 6 at 0.80, A-C 8 at 0.90, B-C 12 at 0.60.
    // HandComputedThreeClusterReport's 0.85 threshold excludes A-C's eight pairs
    // entirely, so that assertion alone cannot see the A-C iteration being
    // skipped. Holding both thresholds pins the count in both directions: drop
    // the A-C iteration and this 26 breaks; count every cross pair regardless of
    // the threshold and the other test's 18 breaks.
    //
    // 0.95 rather than 0.90 keeps the A-C pairs strictly inside the threshold
    // instead of exactly on it. That is a free choice here -- either value
    // admits all eight pairs, so neither is what pins the A-C iteration -- and
    // it is not a stronger assertion. A distance and a threshold written as
    // the same decimal literal are bit-identical doubles whatever that literal
    // rounds to, so an exact-boundary comparison is reliable rather than
    // fragile. BoundaryViolationsThreshold pins that case directly, and it is
    // the only shape that separates an inclusive predicate from a strict one.
    EXPECT_EQ(r.boundary_violations, 26u);
}

// A coverage curve that actually rises. Every other multi-threshold coverage
// assertion in this file is flat, which cannot distinguish indexing
// coverage_thresholds[t] from reading the same entry on every iteration.
TEST(ClusterReportTest, CoverageCurveRisesWithEachThreshold) {
    const DenseStorage storage = MakeThreeClusterStorage();
    const ClusteringResult result = MakeResult({0, 0, 1, 1, 1, 2, 2, 2, 2});

    ClusterReportOptions options;
    options.coverage_thresholds = {0.05, 0.10, 0.15};
    const ClusterReport r = cluster_report(result, storage, options);

    // Distance from each sample to its nearest representative (0, 2 or 5), in
    // sample order: 0, 0.32, 0, 0.12, 0.24, 0, 0.04, 0.06, 0.08. Sorted that is
    // 0, 0, 0, 0.04, 0.06, 0.08, 0.12, 0.24, 0.32, so each threshold admits a
    // different number of samples out of nine.
    ASSERT_EQ(r.coverage_at.size(), 3u);
    EXPECT_DOUBLE_EQ(r.coverage_at[0], 4.0 / 9.0);
    EXPECT_DOUBLE_EQ(r.coverage_at[1], 6.0 / 9.0);
    EXPECT_DOUBLE_EQ(r.coverage_at[2], 7.0 / 9.0);
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

// Section 7.1 item 1. Every value below is worked out by hand on
// MakeSixPointStorage; see the plan's Task 6 for the derivations.
// Medoids are samples 0 and 3; the global medoid M is sample 0.
TEST(ClusterReportTest, HandComputedInternalIndices) {
    const DenseStorage storage = MakeSixPointStorage();
    const ClusteringResult result = MakeResult({0, 0, 0, 1, 1, 1});
    const ClusterReport r = cluster_report(result, storage, ClusterReportOptions());

    // CH = [3*0 + 3*0.8^2] / 1  /  [(0.2^2 + 0.2^2)*2 / (6 - 2)] = 1.92 / 0.04.
    EXPECT_NEAR(r.calinski_harabasz_medoid, 48.0, 1e-9);
    // S_k = (0 + 0.2 + 0.2)/3 for both; DB = (S_a + S_b)/0.8.
    EXPECT_NEAR(r.davies_bouldin_medoid, 1.0 / 3.0, 1e-12);
    // 0.8 / ((0.2 + 0.2 + 0.4)/3).
    EXPECT_NEAR(r.dunn_mean_separation_mean_diameter, 3.0, 1e-12);
    // 0.8 / (2 * 0.4/3).
    EXPECT_NEAR(r.dunn_medoid_separation_medoid_spread, 3.0, 1e-12);
    // Six within-pairs, nine between-pairs; r_pb^2 = 96/101 exactly.
    EXPECT_NEAR(r.point_biserial, std::sqrt(96.0 / 101.0), 1e-12);
}

// The companion to the test above, and the one that actually separates the
// five formulas. MakeSixPointStorage is too symmetric to do that: every
// cluster there has two members and every cross distance is 0.8, so a
// between-scatter that overwrites instead of accumulating, a Davies-Bouldin
// row that sums instead of maximising, a mean-within divided by member count
// instead of pair count, a medoid-Dunn wired to the wrong denominator, and a
// point-biserial pair count taken from the sample count all produce the
// correct answer there. Each of those five reads a different wrong number off
// this fixture.
//
// Three properties are load-bearing and none is incidental:
//
//   * Cluster B has THREE members and a unique medoid (intra sums 1.4375,
//     1.0, 1.0625, so sample 3 wins outright rather than by a tiebreak).
//     Three members is what puts a cluster in the fixture whose pair count
//     C(3,2) equals its member count, so the pair-versus-member denominator
//     error has to be separated by A and C, where the two differ.
//   * d(0,3) = 0.75 is the single asymmetry in the A-B block. It pulls the
//     minimum mean separation (41/48) off the minimum medoid separation
//     (3/4), so swapping the two Dunn numerators would fail too.
//   * d(2,4) = 0.75 rather than 0.875. At 0.875 cluster B's mean-within rises
//     to exactly A's 5/8, the wrong denominator's maximum coincides with the
//     right one, and the mutation survives. 0.75 keeps B at 7/12, strictly
//     below A, while leaving B's medoid unique. Anyone retuning this fixture
//     has to re-check both of those at once.
//
// Every distance is a dyadic multiple of 1/16, so all 21 pairs, both scatters
// and every quotient below are exact in binary.
TEST(ClusterReportTest, HandComputedInternalIndicesOnUnequalClusterSizes) {
    DenseStorage storage(7);
    // A = {0,1}, B = {2,3,4}, C = {5,6}.
    storage.Set(0, 1, 0.625);
    storage.Set(2, 3, 0.6875);
    storage.Set(2, 4, 0.75);
    storage.Set(3, 4, 0.3125);
    storage.Set(5, 6, 0.25);
    for (size_t i = 0; i < 2; ++i) {
        for (size_t j = 2; j < 5; ++j) {
            storage.Set(i, j, 0.875);
        }
        for (size_t j = 5; j < 7; ++j) {
            storage.Set(i, j, 0.9375);
        }
    }
    storage.Set(0, 3, 0.75);  // The one asymmetry; see the note above.
    for (size_t i = 2; i < 5; ++i) {
        for (size_t j = 5; j < 7; ++j) {
            storage.Set(i, j, 0.9375);
        }
    }

    const ClusteringResult result = MakeResult({0, 0, 1, 1, 1, 2, 2});
    const ClusterReport r = cluster_report(result, storage, ClusterReportOptions());

    // True medoids 0, 3 and 5. Totals to all clustered points are 5, 41/8,
    // 81/16, 9/2, 75/16, 79/16, 79/16, so M = 3 is the unique minimum and no
    // tiebreak is involved.
    //
    // Between-scatter = 2*d(0,3)^2 + 3*0 + 2*d(5,3)^2 = 9/8 + 225/128
    //                 = 369/128.
    // Within-scatter  = 25/64 + 73/128 + 1/16 = 131/128.
    // CH = (369/128 / 2) / (131/128 / 4) = 738/131.
    EXPECT_NEAR(r.calinski_harabasz_medoid, 738.0 / 131.0, 1e-12);

    // Scatters S = 5/16, 1/3, 1/8; separations d(0,3) = 3/4 and
    // d(0,5) = d(3,5) = 15/16. Row maxima are 31/36, 31/36 and 22/45, so
    // DB = (31/36 + 31/36 + 22/45) / 3 = (199/90) / 3 = 199/270. Summing each
    // row instead of maximising it gives 109/90.
    EXPECT_NEAR(r.davies_bouldin_medoid, 199.0 / 270.0, 1e-12);

    // Mean separations are 41/48 (A-B), 15/16 and 15/16, so the minimum is
    // 41/48. Mean within-pair distances are 5/8, 7/12 and 1/4 over PAIR
    // counts 1, 3 and 1, so the maximum is A's 5/8 and the index is
    // (41/48)/(5/8) = 41/30. Dividing by member count instead would make B's
    // 7/12 the maximum and read 41/28.
    EXPECT_NEAR(r.dunn_mean_separation_mean_diameter, 41.0 / 30.0, 1e-12);

    // The medoid variant takes its numerator from the medoid separations, so
    // 3/4 rather than 41/48, and its denominator from twice the scatters --
    // 5/8, 2/3, 1/4 -- so 2/3 rather than the 5/8 the mean variant uses.
    // (3/4)/(2/3) = 9/8. Wiring it to max_mean_within instead reads 6/5.
    EXPECT_NEAR(r.dunn_medoid_separation_medoid_spread, 9.0 / 8.0, 1e-12);

    // Five within-pairs against sixteen between-pairs, twenty-one in all.
    // Mean difference 29/32 - 21/40 = 61/160; population variance 269/7056.
    // The cluster sizes are unequal, so P_w = 5 is no longer the clustered
    // sample count the way it was on the six-point fixture: reading the count
    // from there gives sqrt(112/269) in place of sqrt(80/269) and the index
    // reads 0.984 instead. EXPECT_NEAR rather than EXPECT_DOUBLE_EQ because
    // the implementation reaches this through Welford plus merge_moments and
    // lands two ULP off the closed form.
    EXPECT_NEAR(
        r.point_biserial,
        (61.0 / 160.0) / std::sqrt(269.0 / 7056.0) * std::sqrt(80.0) / 21.0,
        1e-12);

    // The silhouette b term is looked up per sample, and this is the only
    // fixture in the file where two members of one cluster disagree about it.
    // The d(0,3) = 0.75 asymmetry pulls sample 0's mean distance to B down to
    // 2.5/3 = 5/6 while sample 1 stays at 2.625/3 = 7/8, and both sit below
    // the 15/16 each sees to C, so A's two b terms are 5/6 and 7/8. Both a
    // terms are the lone intra distance, 5/8. The per-point silhouettes are
    // (5/6 - 5/8)/(5/6) = 1/4 and (7/8 - 5/8)/(7/8) = 2/7, so the record
    // averages them to 15/56.
    //
    // Sample 0 is also cluster A's ordinal, and it holds the smaller of the
    // two b terms. A lookup that reaches for the ordinal's slot rather than
    // the point's -- either where the record loop reads the b term or where
    // the cross pass carries its running minimum from one cluster pair to the
    // next -- therefore gives sample 1 the 5/6 as well, and the record reads
    // 1/4. Every other multi-cluster fixture here has a flat cross block
    // within each pair, or puts the larger b term on the first member where
    // the running minimum discards it; either shape hides the substitution.
    //
    // EXPECT_NEAR rather than EXPECT_DOUBLE_EQ because the division by
    // cluster B's three members is the one quotient on this fixture that
    // leaves the dyadic rationals.
    ClusterReportOptions with_records;
    with_records.compute_per_cluster_records = true;
    const ClusterReport detailed = cluster_report(result, storage, with_records);
    ASSERT_EQ(detailed.records.size(), 3u);
    EXPECT_NEAR(detailed.records[0].silhouette, 15.0 / 56.0, 1e-12);
}

// Section 7.1 item 7. This is the one arithmetic error in the design that would
// otherwise produce entirely plausible numbers: S_k divides by n_k, while
// medoid_member_means divides by n_k - 1, and the two are a factor of two apart
// at n_k == 2.
TEST(ClusterReportTest, MedoidScatterUsesTheClusterSizeDenominator) {
    DenseStorage storage(4);
    storage.Set(0, 1, 0.4);
    storage.Set(2, 3, 0.4);
    storage.Set(0, 2, 0.8);
    storage.Set(0, 3, 0.8);
    storage.Set(1, 2, 0.8);
    storage.Set(1, 3, 0.8);
    const ClusteringResult result = MakeResult({0, 0, 1, 1});
    const ClusterReport r = cluster_report(result, storage, ClusterReportOptions());

    // S_k = (0 + 0.4)/2 = 0.2 for both clusters; DB = (0.2 + 0.2)/0.8 = 0.5.
    // With the n_k - 1 denominator it would be (0.4 + 0.4)/0.8 = 1.0.
    EXPECT_NEAR(r.davies_bouldin_medoid, 0.5, 1e-12);
    // The existing field keeps its own denominator and its 5.0.0 value.
    EXPECT_DOUBLE_EQ(r.median_medoid_member_distance, 0.4);
}

// Section 7.1 item 6. A fixture where the two representatives coincide would
// pass whatever the implementation did, so MakeDivergentRepresentativeStorage
// (Task 1) is built to make them differ, and the test asserts that they differ
// before asserting anything else.
TEST(ClusterReportTest, MedoidNamedIndicesIgnoreRepresentativeMethod) {
    const DenseStorage storage = MakeDivergentRepresentativeStorage();
    const ClusteringResult result = MakeResult({0, 0, 0, 0, 1, 1, 1, 1});

    ClusterReportOptions medoid_options;
    medoid_options.representative_method = RepresentativeMethod::Medoid;
    ClusterReportOptions minimax_options;
    minimax_options.representative_method = RepresentativeMethod::Minimax;

    const ClusterReport by_medoid = cluster_report(result, storage, medoid_options);
    const ClusterReport by_minimax = cluster_report(result, storage, minimax_options);

    // The guard that the fixture actually diverged. Radius is measured from the
    // configured representative, so 0.7 against 0.5 is the two runs naming
    // different points -- checked here rather than through records[], so that
    // this test needs nothing from Task 8.
    ASSERT_DOUBLE_EQ(by_medoid.median_radius, 0.7);
    ASSERT_DOUBLE_EQ(by_minimax.median_radius, 0.5);

    EXPECT_DOUBLE_EQ(by_medoid.calinski_harabasz_medoid,
                     by_minimax.calinski_harabasz_medoid);
    EXPECT_DOUBLE_EQ(by_medoid.davies_bouldin_medoid,
                     by_minimax.davies_bouldin_medoid);
    EXPECT_DOUBLE_EQ(by_medoid.dunn_medoid_separation_medoid_spread,
                     by_minimax.dunn_medoid_separation_medoid_spread);

    // Both configured-representative fields move, and neither medoid-named
    // field above did.
    EXPECT_DOUBLE_EQ(by_medoid.representative_redundancy, 0.95);
    EXPECT_DOUBLE_EQ(by_minimax.representative_redundancy, 0.85);
}

// The non-Medoid representative_redundancy branch reduced over three
// representatives rather than two. Every other test of that branch runs on a
// two-cluster fixture, and with two representatives the row-minima vector is
// [d, d]: the reduction returns d whether it starts at row 0 or stops before
// the last one, whether it keeps the least partner distance or the last one,
// and whether it finishes with a median or a mean. K = 3 is the smallest
// partition that separates them.
//
// The labels are permuted, and the permutation is as load-bearing as the
// distances. At K = 3 the smallest of the three distances is a row minimum for
// two of the three rows, so two minima always tie at the global minimum and the
// median is always that minimum. If the lone large row minimum sits at either
// end of the vector, an off-by-one that drops the row at the OTHER end leaves a
// two-entry vector whose even-count median is still that same minimum, and the
// assertion below cannot see it. {2, 0, 1} puts the large minimum in the middle,
// where dropping either end moves the answer. Restoring the labels to {0, 1, 2}
// would keep this test green while disarming half of what it claims to pin.
//
// The permutation gives members [{1}, {2}, {0}] and therefore representatives
// [1, 2, 0], whose row minima are
//
//   [min(0.8, 0.2), min(0.8, 0.9), min(0.2, 0.9)] = [0.2, 0.8, 0.2]
//
// and each candidate reduction over that vector names its own number: the
// median is 0.2, the mean is 0.4, dropping either the first or the last row
// gives 0.5, dropping the first or the last partner column gives 0.9 and 0.8,
// and folding to the largest partner distance or to the last one gives 0.9.
TEST(ClusterReportTest, NonMedoidRepresentativeRedundancyIsTheMedianOfThreeRowMinima) {
    const DenseStorage storage = MakeThreeSingletonStorage();
    const ClusteringResult result = MakeResult({2, 0, 1});
    ClusterReportOptions options;
    options.representative_method = RepresentativeMethod::Minimax;
    const ClusterReport r = cluster_report(result, storage, options);

    // The guard that the reduction really saw three rows. If the partition
    // collapsed to two clusters the assertion below would hold for the wrong
    // reason, which is exactly the degeneracy this test exists to remove.
    ASSERT_EQ(r.num_clusters, 3u);
    ASSERT_EQ(r.num_singletons, 3u);

    // Median of [0.2, 0.8, 0.2]. The mean of the same three is 0.4, and an
    // off-by-one at either end of the outer loop gives 0.5.
    EXPECT_DOUBLE_EQ(r.representative_redundancy, 0.2);
}

// The K = 4 companion to the test above, and the only place in this file where
// the median of the row minima is separated from their minimum. Two one-token
// reductions of the redundancy loops ship green against every other assertion
// of representative_redundancy here, because every one of them is at K = 2 or
// K = 3 where the median of the row minima IS their minimum: reducing with
// std::min_element instead of detail::median_distance, and hoisting the per-row
// `smallest` above the outer loop so the pushed vector becomes the prefix
// minima of the row minima. Both land on 0.2 against the shipped 0.35.
//
// Over the fixture's [0.2, 0.6, 0.5, 0.2] every one-line reduction names its
// own number:
//
//   median of the row minima (shipped)     0.35
//   minimum instead of median              0.2
//   prefix minima (hoisted running min)    0.2
//   mean instead of median                 0.375
//   maximum instead of median              0.6
//   outer loop drops the first row         0.5
//   outer loop drops the last row          0.5
//   inner loop drops the first column      0.5
//   inner loop drops the last column       0.6
//   i != j narrowed to i < j               0.55
//   i != j dropped, self at distance 0     0.0
//   the min fold written as max            0.85
//   last partner kept instead of the least 0.5
//
// The mean at 0.375 is the closest of these to the shipped 0.35, and
// EXPECT_DOUBLE_EQ still separates the two.
//
// Both paths need their own test. The two reductions live in separate loops --
// the non-Medoid one over the configured representatives, the Medoid one folded
// into the Davies-Bouldin walk over the true medoids -- so a mutation of either
// is invisible to the other's test. Singletons make the two loops read the same
// points, which is what lets one fixture pin both.
TEST(ClusterReportTest, NonMedoidRepresentativeRedundancyIsTheMedianOfFourRowMinima) {
    const DenseStorage storage = MakeFourSingletonStorage();
    const ClusteringResult result = MakeResult({0, 1, 2, 3});
    ClusterReportOptions options;
    options.representative_method = RepresentativeMethod::Minimax;
    const ClusterReport r = cluster_report(result, storage, options);

    // The guard that the reduction really saw four rows.
    ASSERT_EQ(r.num_clusters, 4u);
    ASSERT_EQ(r.num_singletons, 4u);

    EXPECT_DOUBLE_EQ(r.representative_redundancy, 0.35);
}

TEST(ClusterReportTest, MedoidRepresentativeRedundancyIsTheMedianOfFourRowMinima) {
    const DenseStorage storage = MakeFourSingletonStorage();
    const ClusteringResult result = MakeResult({0, 1, 2, 3});
    ClusterReportOptions options;
    options.representative_method = RepresentativeMethod::Medoid;
    const ClusterReport r = cluster_report(result, storage, options);

    ASSERT_EQ(r.num_clusters, 4u);
    ASSERT_EQ(r.num_singletons, 4u);

    // The guard that the answer came from the medoid walk and not from the
    // non-Medoid loop, which this method switches off. Davies-Bouldin is
    // computed in the same walk, and every singleton scatter is zero, so each
    // ratio is 0 / separation with a non-zero separation.
    ASSERT_DOUBLE_EQ(r.davies_bouldin_medoid, 0.0);

    EXPECT_DOUBLE_EQ(r.representative_redundancy, 0.35);
}

// Section 7.1 item 15. Ties in M resolve to the lowest sample index, which is
// not the same rule as the "earliest member" tie-break for m_k: Butina emits
// members in representative-first order, so the two can differ. This fixture
// makes them differ, and pins a single Calinski-Harabasz value that only the
// correct pair of rules produces.
//
// Every distance is a dyadic fraction. That is required, not stylistic: the
// two tied totals are accumulated in different orders (0 sums intra-then-cross
// over {4,3}; 3 sums one intra then three cross), and with values like 0.2 the
// two sums land two ULPs apart, so the tie the test depends on would not be a
// tie at all and the strict < would pick a winner for the wrong reason.
TEST(ClusterReportTest, GlobalMedoidTieResolvesToTheLowestSampleIndex) {
    DenseStorage storage(5);
    storage.Set(0, 1, 0.25);
    storage.Set(0, 2, 0.25);
    storage.Set(1, 2, 0.5);
    storage.Set(3, 4, 0.25);
    for (size_t i = 0; i < 3; ++i) {
        storage.Set(i, 3, 0.5);
        storage.Set(i, 4, 0.75);
    }

    // Members {4,3} for the second cluster, the order Butina would emit if 4
    // were its representative. The earliest-member tie-break therefore makes
    // m_b = 4, not 3.
    const ClusteringResult result(
        std::vector<ClusterLabel>{0, 0, 0, 1, 1}, Clusters{{0, 1, 2}, {4, 3}});
    const ClusterReport r = cluster_report(result, storage, ClusterReportOptions());

    // Totals to all clustered points: 0 -> 1.75, 1 -> 2.0, 2 -> 2.0,
    // 3 -> 1.75, 4 -> 2.5. Samples 0 and 3 tie at the minimum, so M = 0.
    // Between-scatter = 3*d(0,0)^2 + 2*d(4,0)^2 = 2 * 0.5625 = 1.125.
    // Within-scatter  = (0 + 0.0625 + 0.0625) + (0 + 0.0625) = 0.1875.
    // CH = (1.125 / 1) / (0.1875 / 3) = 18.
    //
    // Each wrong rule lands somewhere else and is caught: M = 3 gives 14,
    // and m_b = 3 (lowest index rather than earliest member) gives 8.
    EXPECT_DOUBLE_EQ(r.calinski_harabasz_medoid, 18.0);
}

// The companion to the tie test above: that one pins WHICH clustered point
// wins M, this one pins that only clustered points are eligible. Noise points
// never have point_total written -- the intra pass writes it for cluster
// members and the cross pass adds to it for cluster members -- so a noise
// point keeps its 0.0 initialiser and wins the minimum outright if the
// owner[i] != NO_OWNER guard is dropped.
//
// Two properties of sample 5's row are load-bearing rather than tidy. Every
// distance on it is finite, which denies the dropped guard the refusal escape
// hatch: a NaN there would make the mutant throw instead of reporting a
// number, and that is exactly how this guard was incidentally covered before
// this test existed -- by an unrelated coverage-threshold fixture that sets
// d(0, 5) to NaN, whose failure message named finiteness and pointed at the
// wrong part of the file. And each one is set to 0.5 rather than left on
// DenseStorage's 0.0 fill, so the mutant reports a thoroughly plausible 20.0
// rather than a degenerate 0.0 that no reader would mistake for an answer.
//
// Those same two properties make this the only fixture here that can pin
// point-biserial's pair count to clustered pairs, so the test carries a second
// assertion on that field; both pin the one contract its name states.
TEST(ClusterReportTest, GlobalMedoidIsChosenAmongClusteredPointsOnly) {
    DenseStorage storage(6);
    storage.Set(0, 1, 0.25);
    storage.Set(0, 2, 0.25);
    storage.Set(1, 2, 0.5);
    storage.Set(3, 4, 0.25);
    for (size_t i = 0; i < 3; ++i) {
        storage.Set(i, 3, 0.5);
        storage.Set(i, 4, 0.5);
    }
    for (size_t i = 0; i < 5; ++i) {
        storage.Set(i, 5, 0.5);
    }

    const ClusteringResult result = MakeResult({0, 0, 0, 1, 1, -1});
    const ClusterReport r = cluster_report(result, storage, ClusterReportOptions());

    // Totals over clustered points: 0 -> 1.5, 1 -> 1.75, 2 -> 1.75,
    // 3 -> 1.75, 4 -> 1.75, and noise sample 5 -> 0.0. So M = 0.
    // Medoids are 0 and 3. Between-scatter = 3*d(0,0)^2 + 2*d(3,0)^2 = 0.5,
    // within-scatter = (0 + 0.0625 + 0.0625) + (0 + 0.0625) = 0.1875,
    // CH = (0.5 / 1) / (0.1875 / 3) = 8.
    //
    // Drop the guard and M becomes sample 5, whose 0.0 total is unbeatable:
    // between-scatter rises to 3*0.25 + 2*0.25 = 1.25 and CH reads 20. Both
    // numbers are entirely plausible, which is the point.
    EXPECT_DOUBLE_EQ(r.calinski_harabasz_medoid, 8.0);

    // The same contract one level down. Four within-pairs (0.25, 0.25, 0.5 and
    // 0.25, mean 0.3125) and six between-pairs (0.5 each) make ten clustered
    // pairs, against the fifteen of C(6, 2); the five pairs on noise sample 5
    // are exactly the difference. Pooled over the ten: mean 0.425, population
    // variance 0.013125, so the index is
    // (0.5 - 0.3125)/sqrt(0.013125) * sqrt(4*6)/10 = 1.5*sqrt(2/7). Divide by
    // fifteen instead and it reads 0.535, an unremarkable correlation. Every
    // other fixture in this file that reaches this line is noise-free, where
    // the two denominators coincide, and every noise-bearing one returns NaN
    // before it -- this is the only place the two can be told apart.
    EXPECT_NEAR(r.point_biserial, 1.5 * std::sqrt(2.0 / 7.0), 1e-12);
}

// The third property of M, after "which point wins a tie" and "which points are
// eligible": the search runs over every clustered point, not over the K elected
// medoids. Restricting it to true_medoids reads as a cheap optimisation and no
// other fixture here refutes it, because on each of them the winner happens to
// be somebody's medoid as well.
TEST(ClusterReportTest, GlobalMedoidIsNotRestrictedToClusterMedoids) {
    const DenseStorage storage = MakeOffMedoidGlobalStorage();
    const ClusterReport r =
        cluster_report(MakeResult({0, 0, 1, 1}), storage, ClusterReportOptions());

    // Medoids 0 and 2, by the earliest-member rule on two internal ties.
    // Totals over the clustered set are 0 -> 9/4, 1 -> 3/4, 2 -> 3/2 and
    // 3 -> 3/2, so M = 1 outright: the one clustered point that is not a
    // medoid, and the reason this fixture exists.
    //
    // Between-scatter = 2*d(0,1)^2 + 2*d(2,1)^2 = 1/2 + 1/32 = 17/32;
    // within-scatter = 1/4 + 1/4 = 1/2. With N = 4 and K = 2 the index is
    // (17/32) / ((1/2)/2) = 17/8.
    //
    // Scan only the medoids and 2 wins with 3/2, which puts M on top of a
    // medoid: between-scatter becomes 2*d(0,2)^2 = 49/32 and the index reads
    // 49/8. Both are ordinary-looking Calinski-Harabasz scores.
    EXPECT_DOUBLE_EQ(r.calinski_harabasz_medoid, 17.0 / 8.0);
}

// The medoid Dunn variant divides the smallest distance between two MEDOIDS by
// the largest medoid spread. min_inter -- the smallest distance between two
// members of different clusters -- is in scope at that line and already feeds
// the other Dunn variant, so substituting it is a one-token slip that a flat
// cross block cannot detect, since there every member distance is also the
// medoid distance.
TEST(ClusterReportTest, MedoidDunnUsesMedoidSeparationNotMemberSeparation) {
    const DenseStorage storage = MakeOffMedoidGlobalStorage();
    const ClusterReport r =
        cluster_report(MakeResult({0, 0, 1, 1}), storage, ClusterReportOptions());

    // Both medoid scatters are (0 + 1/2)/2 = 1/4, so the spread denominator is
    // 2 * 1/4 = 1/2. The only medoid separation is d(0,2) = 7/8, giving
    // (7/8)/(1/2) = 7/4.
    //
    // The smallest cross-cluster member distance is d(1,2) = 1/8, so a
    // numerator taken from there reads 1/4 -- a seventh of the right answer,
    // and still a plausible Dunn score rather than a visible failure.
    EXPECT_DOUBLE_EQ(r.dunn_medoid_separation_medoid_spread, 7.0 / 4.0);
}

// point_biserial is a signed correlation, and the sign is the part of it that
// carries the verdict: a labelling that groups the far pairs together has to
// score negative. Every other fixture in this file separates its clusters in
// the expected direction, so none of them would notice the mean difference
// being wrapped in std::fabs.
TEST(ClusterReportTest, PointBiserialKeepsTheSignOfTheSeparation) {
    // The separation inverted: the two within-cluster pairs are the far ones
    // and all four cross pairs are the near ones.
    DenseStorage storage(4);
    storage.Set(0, 1, 0.75);
    storage.Set(2, 3, 0.75);
    storage.Set(0, 2, 0.25);
    storage.Set(0, 3, 0.25);
    storage.Set(1, 2, 0.25);
    storage.Set(1, 3, 0.25);
    const ClusterReport r =
        cluster_report(MakeResult({0, 0, 1, 1}), storage, ClusterReportOptions());

    // Two within-pairs at 3/4 against four between-pairs at 1/4 pool to a mean
    // of 5/12 and a population variance of 1/18 over the six clustered pairs,
    // so the index is ((1/4 - 3/4)/sqrt(1/18)) * sqrt(2*4)/6 = -1. Distance is
    // a perfect decreasing function of the within/between indicator here, so
    // the floor of the statistic is the correct answer. Take the absolute
    // value of the difference and the worst possible clustering reports +1,
    // the best possible score.
    //
    // EXPECT_NEAR rather than EXPECT_DOUBLE_EQ for the reason given above: the
    // implementation reaches this through Welford updates, merge_moments and a
    // square root rather than through the closed form.
    EXPECT_NEAR(r.point_biserial, -1.0, 1e-12);
}

// Section 7.1 item 11 and the section 5.5 table.
TEST(ClusterReportTest, InternalIndicesUndefinedCases) {
    const DenseStorage storage = MakeSixPointStorage();

    const ClusterReport single =
        cluster_report(MakeResult({0, 0, 0, -1, -1, -1}), storage, ClusterReportOptions());
    EXPECT_TRUE(std::isnan(single.calinski_harabasz_medoid));
    EXPECT_TRUE(std::isnan(single.davies_bouldin_medoid));
    EXPECT_TRUE(std::isnan(single.dunn_mean_separation_mean_diameter));
    EXPECT_TRUE(std::isnan(single.dunn_medoid_separation_medoid_spread));
    // The one scenario here with within-pairs but no between-pairs, and the
    // only assertion that isolates point-biserial's between-count predicate
    // from its within-count one. It is not a vacuous NaN check: with the
    // between predicate removed the mutant does not throw or NaN, because
    // DistanceMoments::mean is 0.0 at count 0 and merge_moments returns the
    // non-empty stream untouched, so sqrt(P_w * 0) zeroes the product and the
    // field reads a finite -0.0 -- a perfectly plausible "no correlation"
    // answer for a clustering that cannot have one.
    EXPECT_TRUE(std::isnan(single.point_biserial));

    const ClusterReport none = cluster_report(
        ClusteringResult(std::vector<ClusterLabel>{-1, -1, -1, -1, -1, -1}, Clusters{}),
        storage,
        ClusterReportOptions());
    EXPECT_TRUE(std::isnan(none.calinski_harabasz_medoid));
    EXPECT_TRUE(std::isnan(none.davies_bouldin_medoid));
    EXPECT_TRUE(std::isnan(none.dunn_mean_separation_mean_diameter));
    EXPECT_TRUE(std::isnan(none.dunn_medoid_separation_medoid_spread));
    EXPECT_TRUE(std::isnan(none.point_biserial));

    // All singletons: Nc == K, so CH's denominator has no degrees of freedom,
    // and there are no within-pairs for point-biserial.
    const ClusterReport singletons =
        cluster_report(MakeResult({0, 1, 2, 3, 4, 5}), storage, ClusterReportOptions());
    EXPECT_TRUE(std::isnan(singletons.calinski_harabasz_medoid));
    EXPECT_TRUE(std::isnan(singletons.point_biserial));

    // Constant distances leave point-biserial's s_d at zero.
    DenseStorage flat(4);
    for (size_t i = 0; i < 4; ++i) {
        for (size_t j = i + 1; j < 4; ++j) {
            flat.Set(i, j, 0.5);
        }
    }
    const ClusterReport constant =
        cluster_report(MakeResult({0, 0, 1, 1}), flat, ClusterReportOptions());
    EXPECT_TRUE(std::isnan(constant.point_biserial));
}

// The three zero-denominator guards that no scenario above reaches. The
// all-singletons case fails the clustered_count > cluster_count test and takes
// the else-NaN arm before CH's denominator is ever formed, and the flat-0.5
// case has a medoid scatter of 0.25 per cluster, so neither Dunn denominator
// is zero there. Two perfectly tight clusters at a positive distance is the
// shape that drives all three denominators to zero at once while still
// entering the guarded branch.
//
// The two zero distances are set explicitly rather than left to
// DenseStorage's zero fill because they are the mechanism of the test, not
// incidental to it: deleting them as redundant would leave every assertion
// passing and the test toothless about what it is pinning.
//
// The two finite values carry as much weight as the three NaNs. They are what
// show this is a well-formed input the report answers rather than a
// degenerate one it gives up on, which is what makes the NaNs a deliberate
// "undefined" rather than a side effect of arithmetic falling over.
TEST(ClusterReportTest, ZeroWithinScatterLeavesCalinskiHarabaszAndDunnUndefined) {
    DenseStorage storage(4);
    storage.Set(0, 1, 0.0);
    storage.Set(2, 3, 0.0);
    storage.Set(0, 2, 0.5);
    storage.Set(0, 3, 0.5);
    storage.Set(1, 2, 0.5);
    storage.Set(1, 3, 0.5);
    const ClusterReport r =
        cluster_report(MakeResult({0, 0, 1, 1}), storage, ClusterReportOptions());

    // Medoids 0 and 2; every point total is 1.0, so M = 0 by the tiebreak.
    // Between-scatter = 2 * 0.5^2 = 0.5 and within-scatter = 0, and
    // clustered_count (4) exceeds cluster_count (2), so CH forms the quotient
    // and finds 0/(4-2) == 0 underneath it. Unguarded it would read +inf.
    EXPECT_TRUE(std::isnan(r.calinski_harabasz_medoid));
    // Both mean within-pair distances are 0, so the maximum is 0 against a
    // minimum mean separation of 0.5. Unguarded: +inf.
    EXPECT_TRUE(std::isnan(r.dunn_mean_separation_mean_diameter));
    // Both medoid spreads are 0 against a minimum medoid separation of 0.5.
    // Unguarded: +inf, and by a different guard than the line above.
    EXPECT_TRUE(std::isnan(r.dunn_medoid_separation_medoid_spread));

    // Defined, and genuinely 0: the scatters are zero and the separation is
    // not, so the ratio is a real best-possible score rather than a missing
    // one. Contrast CoincidentZeroScatterSingletonsGiveInfiniteDaviesBouldin,
    // where the separation is zero too and the answer is +inf.
    EXPECT_DOUBLE_EQ(r.davies_bouldin_medoid, 0.0);
    // Also defined, and at its ceiling: distance is an exact linear function
    // of the within/between indicator here, 0.0 inside and 0.5 across, so the
    // correlation is perfect.
    EXPECT_DOUBLE_EQ(r.point_biserial, 1.0);
}

// Coincident medoids give inf, not NaN: "two clusters share a representative"
// is a real answer, and NaN would hide it.
TEST(ClusterReportTest, CoincidentMedoidsGiveInfiniteDaviesBouldin) {
    DenseStorage storage(4);
    storage.Set(0, 1, 0.4);
    storage.Set(2, 3, 0.4);
    storage.Set(0, 2, 0.0);
    storage.Set(0, 3, 0.8);
    storage.Set(1, 2, 0.8);
    storage.Set(1, 3, 0.8);
    const ClusterReport r =
        cluster_report(MakeResult({0, 0, 1, 1}), storage, ClusterReportOptions());
    EXPECT_TRUE(std::isinf(r.davies_bouldin_medoid));
    EXPECT_FALSE(std::isnan(r.davies_bouldin_medoid));
}

// The case the explicit separation == 0.0 branch exists for. Both scatters are
// zero as well, so the arithmetic reads 0.0/0.0 == NaN, and std::max(0.0, NaN)
// keeps the 0.0 -- the most degenerate clustering possible would otherwise
// report a perfect Davies-Bouldin. The test above cannot catch this: its
// scatters are 0.2, so the division already yields inf.
TEST(ClusterReportTest, CoincidentZeroScatterSingletonsGiveInfiniteDaviesBouldin) {
    DenseStorage storage(2);
    storage.Set(0, 1, 0.0);
    const ClusterReport r =
        cluster_report(MakeResult({0, 1}), storage, ClusterReportOptions());
    ASSERT_FALSE(std::isnan(r.davies_bouldin_medoid));
    EXPECT_TRUE(std::isinf(r.davies_bouldin_medoid));
    EXPECT_NE(r.davies_bouldin_medoid, 0.0);
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

    // The first ordinal in that message is read out of owner, which is indexed
    // by sample rather than by cluster. The case above cannot separate the two
    // index spaces: the duplicated sample is 2 and the scanning cluster's
    // ordinal is 1, and owner holds 0 at both subscripts, so a lookup keyed on
    // the ordinal prints the same sentence. Three clusters pull them apart --
    // sample 1 is the duplicate and belongs to cluster 0, while sample 2, the
    // index the scanning ordinal would supply, belongs to cluster 1.
    const ClusteringResult across_three(
        std::vector<ClusterLabel>{0, 0, 1, 1, 2, 2}, Clusters{{0, 1}, {2, 3}, {1, 4}});
    try {
        cluster_report(across_three, storage, ClusterReportOptions());
        FAIL() << "expected std::invalid_argument";
    } catch (const std::invalid_argument& e) {
        const std::string message(e.what());
        EXPECT_NE(message.find("sample 1 appears in clusters 0 and 2"),
                  std::string::npos)
            << message;
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

    // The cluster this message names also comes out of owner at the sample's
    // own subscript, and the label that triggered the refusal is a cluster
    // ordinal sitting in scope beside it. Above, the two coincide: sample 2
    // carries label 1, and owner[1] is 0, exactly what owner[2] holds. The
    // mismatch has to sit on a sample whose label does not point back at its
    // own cluster for the reads to diverge -- sample 3 is in cluster 1 and
    // carries label 0, and owner[0] is 0 rather than 1.
    const ClusteringResult label_points_elsewhere(
        std::vector<ClusterLabel>{0, 0, 0, 0, 1, 1}, Clusters{{0, 1, 2}, {3, 4, 5}});
    try {
        cluster_report(label_points_elsewhere, storage, ClusterReportOptions());
        FAIL() << "expected std::invalid_argument";
    } catch (const std::invalid_argument& e) {
        const std::string message(e.what());
        EXPECT_NE(message.find("sample 3 has label 0 but appears in cluster 1"),
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

// Section 5.2 precondition 2. Previously a NaN reached std::sort through
// median_distance(intra_pairs), which is undefined behaviour, so this removes a
// hazard rather than a defined result.
TEST(ClusterReportTest, NonFiniteDistanceIsRefused) {
    const ClusteringResult result = MakeResult({0, 0, 0, 1, 1, 1});

    // Both cases assert the offending pair, not just "not finite". A refusal
    // that names the wrong pair, or names none, is the failure this precondition
    // exists to prevent -- a caller who cannot find the bad cell has been
    // refused without being helped. INVARIANT 1.
    DenseStorage nan_storage = MakeSixPointStorage();
    nan_storage.Set(1, 2, std::numeric_limits<double>::quiet_NaN());
    try {
        cluster_report(result, nan_storage, ClusterReportOptions());
        FAIL() << "expected std::invalid_argument";
    } catch (const std::invalid_argument& e) {
        const std::string message(e.what());
        EXPECT_NE(message.find("not finite"), std::string::npos) << message;
        // Read by the intra pass, which walks cluster {0,1,2} in ascending
        // member order, so 1 is named before 2.
        EXPECT_NE(message.find("samples 1 and 2"), std::string::npos) << message;
    }

    // A cross-cluster pair, reached by the cross pass rather than the intra one.
    DenseStorage inf_storage = MakeSixPointStorage();
    inf_storage.Set(0, 4, std::numeric_limits<double>::infinity());
    try {
        cluster_report(result, inf_storage, ClusterReportOptions());
        FAIL() << "expected std::invalid_argument";
    } catch (const std::invalid_argument& e) {
        const std::string message(e.what());
        EXPECT_NE(message.find("not finite"), std::string::npos) << message;
        EXPECT_NE(message.find("samples 0 and 4"), std::string::npos) << message;
    }
}

// A noise point's distance to a representative is not an intra pair and not a
// cross pair, but the coverage scan below reads it, so it is covered too.
// Cluster {0,1,2} has representative 0 and cluster {3,4} has representative 3,
// so (5, 0) is one of the four reads that scan makes for sample 5.
TEST(ClusterReportTest, NonFiniteNoiseToRepresentativeDistanceIsRefused) {
    DenseStorage storage = MakeSixPointStorage();
    storage.Set(0, 5, std::numeric_limits<double>::quiet_NaN());
    const ClusteringResult result = MakeResult({0, 0, 0, 1, 1, -1});
    try {
        cluster_report(result, storage, ClusterReportOptions());
        FAIL() << "expected std::invalid_argument";
    } catch (const std::invalid_argument& e) {
        const std::string message(e.what());
        EXPECT_NE(message.find("not finite"), std::string::npos) << message;
        // checked_distance(storage, point, representative), so the noise point
        // is named first.
        EXPECT_NE(message.find("samples 5 and 0"), std::string::npos) << message;
    }
}

// The complement of the test above, and the boundary of the coverage
// precondition. The poisoned cell is the same one -- d(5, 0), a read the
// coverage scan makes when thresholds exist -- but with the threshold list
// cleared, coverage_at is empty and no reported value reads it. Refusing here
// would refuse a report every field of which is defined. INVARIANT 3.
TEST(ClusterReportTest, EmptyCoverageThresholdsDoNotRefuseUnreadDistances) {
    DenseStorage storage = MakeSixPointStorage();
    storage.Set(0, 5, std::numeric_limits<double>::quiet_NaN());
    const ClusteringResult result = MakeResult({0, 0, 0, 1, 1, -1});
    ClusterReportOptions options;
    options.coverage_thresholds.clear();

    const ClusterReport report = cluster_report(result, storage, options);

    EXPECT_TRUE(report.coverage_at.empty());
    EXPECT_DOUBLE_EQ(report.median_medoid_member_distance, 0.2);
    EXPECT_DOUBLE_EQ(report.silhouette, 0.7);
    EXPECT_DOUBLE_EQ(report.dunn_index, 2.0);
}

// The complement, and the boundary of the precondition: sample 1 is a member of
// cluster {0,1,2} but not its representative, and sample 5 is noise, so no
// reported value reads d(1, 5). The intra pass covers only within-cluster
// pairs, the cross pass only clustered-to-clustered, and the coverage scan only
// sample-to-representative.
//
// cluster_representative does read d(1, 5) -- its external-neighbour scan
// (Representative.cpp:81-92) touches every item outside the candidate's cluster
// through a raw storage.Get -- but the value it produces cannot reach the
// report: nearest_external_distance and silhouette_like_score are written into
// RepresentativeMetrics and then read by nothing. representative_score
// (Representative.cpp:179-198) consults neither for any of the four methods,
// rank_representatives sorts on score alone, and cluster_representative returns
// ranked.front().member and discards the metrics. So the report below is fully
// determined, and refusing it would be over-refusal. INVARIANT 3.
//
// median_medoid_member_distance is the assertion because it is the field a
// changed representative would move: 0.2 for representatives 0 and 3, but 0.25
// if the NaN pushed cluster {0,1,2} onto member 1.
TEST(ClusterReportTest, NonFiniteDistanceNoReportedValueReadsIsAccepted) {
    DenseStorage storage = MakeSixPointStorage();
    storage.Set(1, 5, std::numeric_limits<double>::quiet_NaN());
    const ClusteringResult result = MakeResult({0, 0, 0, 1, 1, -1});
    const ClusterReport report =
        cluster_report(result, storage, ClusterReportOptions());
    EXPECT_EQ(report.num_noise, 1u);
    EXPECT_DOUBLE_EQ(report.median_medoid_member_distance, 0.2);
}

// The structural precondition is reported first when a call is wrong for both
// reasons: a caller fixes the partition before the matrix means anything.
// INVARIANT 1.
TEST(ClusterReportTest, PartitionErrorOutranksNonFiniteDistance) {
    DenseStorage storage = MakeSixPointStorage();
    storage.Set(1, 2, std::numeric_limits<double>::quiet_NaN());
    const ClusteringResult result(
        std::vector<ClusterLabel>{0, 0, 0, 1, 1, 1}, Clusters{{0, 1, 2}, {2, 3, 4}});
    try {
        cluster_report(result, storage, ClusterReportOptions());
        FAIL() << "expected std::invalid_argument";
    } catch (const std::invalid_argument& e) {
        EXPECT_NE(std::string(e.what()).find("sample 2"), std::string::npos)
            << e.what();
        EXPECT_EQ(std::string(e.what()).find("not finite"), std::string::npos)
            << e.what();
    }
}

// A representative method ClusterReportOptions cannot configure is a
// configuration error, and no repair of the distance matrix rescues it, so it
// is named ahead of the finiteness refusal. INVARIANT 1, same precedence
// argument as the partition check above. The clean case pins the refusal the
// header documents at ClusterReport.h:198-202, which no other test covers.
TEST(ClusterReportTest, UnsupportedRepresentativeMethodOutranksNonFiniteDistance) {
    const ClusteringResult result = MakeResult({0, 0, 0, 1, 1, 1});
    ClusterReportOptions options;
    options.representative_method = RepresentativeMethod::HighestNeighborhood;

    const DenseStorage clean_storage = MakeSixPointStorage();
    EXPECT_THROW(cluster_report(result, clean_storage, options), std::invalid_argument);

    DenseStorage nan_storage = MakeSixPointStorage();
    nan_storage.Set(1, 2, std::numeric_limits<double>::quiet_NaN());
    try {
        cluster_report(result, nan_storage, options);
        FAIL() << "expected std::invalid_argument";
    } catch (const std::invalid_argument& e) {
        const std::string message(e.what());
        // Deliberately not "neighbor threshold": Representative.cpp's
        // validate_options says that too, and matching on it could not tell the
        // two refusals apart. cluster_report must name the option the caller
        // can actually act on.
        EXPECT_NE(
            message.find(
                "cluster_report: representative_method HighestNeighborhood is unsupported"),
            std::string::npos)
            << message;
        EXPECT_EQ(message.find("not finite"), std::string::npos) << message;
    }
}

// A NaN boundary_threshold is not an undefined answer, it is a confident wrong
// one: `distance <= threshold` is false for every pair against a NaN, so the
// unguarded call returned boundary_violations = 0 -- the same number a
// perfectly separated clustering reports, with nothing marking it suspect.
// Passing an explicit NaN must not read as passing nothing. INVARIANT 2.
TEST(ClusterReportTest, NanBoundaryThresholdIsRefused) {
    const DenseStorage storage = MakeThreeClusterStorage();
    const ClusteringResult result = MakeResult({0, 0, 1, 1, 1, 2, 2, 2, 2});
    ClusterReportOptions options;
    options.boundary_threshold = std::numeric_limits<double>::quiet_NaN();

    try {
        cluster_report(result, storage, options);
        FAIL() << "expected std::invalid_argument";
    } catch (const std::invalid_argument& e) {
        // The whole message, not the shared "must not be NaN" tail. The
        // coverage refusal below ends in those same words, and a test matching
        // only the tail would stay green if the two throw sites exchanged
        // messages -- each would then name the option the caller did not set.
        EXPECT_EQ(
            std::string(e.what()),
            "cluster_report: boundary_threshold must not be NaN");
    }
}

// The coverage list is a list, so the refusal has to say which entry is bad --
// the caller cannot act on "one of them is NaN". The NaN sits at index 2 of
// five, so a message that hardcoded either end of the list names the wrong
// entry and fails here.
//
// The second NaN at index 4 is what makes the refusal name the FIRST offender
// rather than merely an offender. With one NaN in the list, a guard rewritten
// to record the offending index and throw after the loop -- naming the last --
// is indistinguishable from the shipped one. A caller told to fix threshold 4
// when threshold 2 is also NaN just earns a second refusal on the next call.
TEST(ClusterReportTest, NanCoverageThresholdIsRefusedByIndex) {
    const DenseStorage storage = MakeThreeClusterStorage();
    const ClusteringResult result = MakeResult({0, 0, 1, 1, 1, 2, 2, 2, 2});
    ClusterReportOptions options;
    options.coverage_thresholds = {
        0.05,
        0.10,
        std::numeric_limits<double>::quiet_NaN(),
        0.15,
        std::numeric_limits<double>::quiet_NaN()};

    try {
        cluster_report(result, storage, options);
        FAIL() << "expected std::invalid_argument";
    } catch (const std::invalid_argument& e) {
        EXPECT_EQ(
            std::string(e.what()),
            "cluster_report: coverage threshold 2 must not be NaN");
    }
}

// A malformed threshold is a configuration error of the same class as an
// unsupported representative_method: no repair of the distance matrix rescues
// it, so it is named ahead of the finiteness refusal. INVARIANT 1, and the same
// precedence argument as UnsupportedRepresentativeMethodOutranksNonFiniteDistance
// above.
//
// This is the only test that combines a NaN threshold with a holey matrix, and
// so the only one that pins WHERE the two threshold guards sit rather than
// merely that they exist. Moving them down beside the code that reads the
// options -- a plausible tidy-up -- puts them after the intra pass, and then
// detail::checked_distance answers first and the caller is told to fix the
// matrix instead of the option they can actually act on.
TEST(ClusterReportTest, NanThresholdOutranksNonFiniteDistance) {
    DenseStorage storage = MakeSixPointStorage();
    // Inside cluster {0,1,2}, so the intra pass reads it and throws.
    storage.Set(1, 2, std::numeric_limits<double>::quiet_NaN());
    const ClusteringResult result = MakeResult({0, 0, 0, 1, 1, 1});
    ClusterReportOptions options;
    options.boundary_threshold = std::numeric_limits<double>::quiet_NaN();

    try {
        cluster_report(result, storage, options);
        FAIL() << "expected std::invalid_argument";
    } catch (const std::invalid_argument& e) {
        const std::string message(e.what());
        EXPECT_EQ(message, "cluster_report: boundary_threshold must not be NaN");
        EXPECT_EQ(message.find("not finite"), std::string::npos) << message;
    }
}

// The over-refusal guard for boundary_threshold. An infinite threshold is a
// well-formed "count every cross pair" request that the report has always
// answered, so a guard written as !std::isfinite would take a defined result
// away. INVARIANT 3.
//
// The assertion is the count. "Did not throw" would also hold if the threshold
// stopped being read at all, and the narrow case below is what separates the 26
// from a cross pass that counts every pair regardless of the option.
TEST(ClusterReportTest, InfiniteBoundaryThresholdIsAcceptedAndCountsEveryCrossPair) {
    const DenseStorage storage = MakeThreeClusterStorage();
    const ClusteringResult result = MakeResult({0, 0, 1, 1, 1, 2, 2, 2, 2});

    ClusterReportOptions wide;
    wide.boundary_threshold = std::numeric_limits<double>::infinity();
    // A-B 6 pairs at 0.80, A-C 8 at 0.90, B-C 12 at 0.60. Twenty-six is none of
    // the neighbouring counts on this fixture: nine samples, three clusters,
    // ten intra pairs, thirty-six pairs in total.
    EXPECT_EQ(cluster_report(result, storage, wide).boundary_violations, 26u);

    ClusterReportOptions narrow;
    narrow.boundary_threshold = 0.60;
    EXPECT_EQ(cluster_report(result, storage, narrow).boundary_violations, 12u);
}

// The same guard for the coverage list. Both a threshold far past the largest
// distance and an infinite one saturate, and the 0.05 entry sharing the list
// with them does not, so the two 1.0s are the scan answering rather than the
// scan having stopped reading. INVARIANT 3.
TEST(ClusterReportTest, LargeAndInfiniteCoverageThresholdsSaturateCoverage) {
    const DenseStorage storage = MakeThreeClusterStorage();
    const ClusteringResult result = MakeResult({0, 0, 1, 1, 1, 2, 2, 2, 2});

    ClusterReportOptions options;
    options.coverage_thresholds = {
        0.05, 1.0e6, std::numeric_limits<double>::infinity()};
    const ClusterReport r = cluster_report(result, storage, options);

    // Distances to the nearest representative (0, 2 or 5), in sample order:
    // 0, 0.32, 0, 0.12, 0.24, 0, 0.04, 0.06, 0.08. Four of the nine are within
    // 0.05, and 0.32 is the farthest, so both wide thresholds admit all nine.
    ASSERT_EQ(r.coverage_at.size(), 3u);
    EXPECT_DOUBLE_EQ(r.coverage_at[0], 4.0 / 9.0);
    EXPECT_DOUBLE_EQ(r.coverage_at[1], 1.0);
    EXPECT_DOUBLE_EQ(r.coverage_at[2], 1.0);
}

// A call wrong in both ways names the fault the caller must fix first. No
// threshold value rescues a representative_method ClusterReportOptions cannot
// configure, so the method outranks the malformed threshold. INVARIANT 1, and
// the only test pinning that the two threshold guards sit behind the method
// refusal rather than ahead of it.
TEST(ClusterReportTest, UnsupportedRepresentativeMethodOutranksNanThreshold) {
    const DenseStorage storage = MakeThreeClusterStorage();
    const ClusteringResult result = MakeResult({0, 0, 1, 1, 1, 2, 2, 2, 2});
    ClusterReportOptions options;
    options.representative_method = RepresentativeMethod::HighestNeighborhood;
    options.boundary_threshold = std::numeric_limits<double>::quiet_NaN();
    options.coverage_thresholds = {
        0.05, std::numeric_limits<double>::quiet_NaN()};

    try {
        cluster_report(result, storage, options);
        FAIL() << "expected std::invalid_argument";
    } catch (const std::invalid_argument& e) {
        const std::string message(e.what());
        EXPECT_NE(
            message.find(
                "cluster_report: representative_method HighestNeighborhood is unsupported"),
            std::string::npos)
            << message;
        // Neither threshold refusal may surface here, and "NaN" appears in both
        // of their messages and in neither of the method's.
        EXPECT_EQ(message.find("NaN"), std::string::npos) << message;
    }
}

// Both NaN refusals sit inside the members-non-empty guard, so a partition with
// no clusters at all still accepts a NaN threshold. That is deliberate: no
// comparison in a K == 0 report consumes either option -- coverage_at stays
// empty and boundary_violations is zero for want of a pair to count -- so there
// is no wrong number for the NaN to hide behind and refusing it would be
// over-refusal. The report is not option-free, though: the NaN coverage
// threshold passed below is echoed verbatim into report.coverage_thresholds, so
// it does come back out. That is the caller's own value returned to them rather
// than a number the report computed, which is why it does not argue for a
// refusal. Nothing else pins the placement -- hoisting the guard out of that
// block is a one-line change no other test in this file would notice.
//
// The rule the guards implement is a placement rule, "at least one cluster", and
// not "would this call have read the threshold". SingleClusterRefusesANanThreshold
// below pins the other side of it: at K == 1 the boundary threshold is provably
// never read and a NaN is still refused. Uniform placement is the more
// predictable contract, and the trade is a real one rather than a free choice:
// at K == 1 the refusal withholds a report that was already fully determined --
// median_radius, p95_diameter, the intra-distance and size statistics and
// coverage_at are all computed there without reading boundary_threshold. The
// production comment beside the guards records the same concession.
TEST(ClusterReportTest, EmptyPartitionAcceptsANanThreshold) {
    const DenseStorage storage = MakeTwoClusterStorage();
    // labels_to_clusters returns an empty Clusters as soon as the largest label
    // is negative, so an all-noise label vector is what empties the member list.
    // treat_noise_as_singletons plays no part -- it only picks the denominator
    // for singleton_fraction, which this test does not assert.
    const ClusteringResult result = MakeResult({-1, -1, -1, -1});
    ClusterReportOptions options;
    options.boundary_threshold = std::numeric_limits<double>::quiet_NaN();
    options.coverage_thresholds = {std::numeric_limits<double>::quiet_NaN()};

    const ClusterReport r = cluster_report(result, storage, options);

    ASSERT_EQ(r.num_clusters, 0u);
    EXPECT_EQ(r.num_noise, 4u);
    // Zero for want of a pair to count, not for want of a comparison that held.
    EXPECT_EQ(r.boundary_violations, 0u);
    EXPECT_TRUE(r.coverage_at.empty());
}

// K == 1, the cluster count the other NaN tests skip: every one of them uses a
// K = 2 or K = 3 partition and the acceptance test above uses K = 0. Both
// guards are therefore pinned only at the two ends, and scoping either of them
// under a cluster count -- wrapping it in `if (cluster_count >= 2)`, which the
// guards' own rationale comment can be read as inviting -- keeps the rest of
// the suite green while restoring the defect the guards exist to close.
//
// It restores it for real, not in principle. The coverage scan is not gated on
// the cluster count; it needs only a non-empty sample count, a non-empty
// representative list and a non-empty threshold list, all of which hold at
// K == 1. So a scoped-away guard lets `nearest <= NaN` come back false for every
// point and publishes coverage_at[0] = 0.0, a plausible "nothing is covered" for
// a question that was never answerable.
//
// The boundary half is here for the opposite reason. At K == 1 the threshold is
// provably never read -- it appears only inside the double cluster loop, which
// no single cluster enters -- and it is refused anyway. That is the placement
// rule stated plainly: the guards refuse from the first cluster onwards, whether
// or not this particular call would have consumed the value.
//
// Both halves share one fixture and one cluster count, so they are one test:
// each is the other's context, and a future reader relaxing either guard needs
// to see both consequences together.
TEST(ClusterReportTest, SingleClusterRefusesANanThreshold) {
    const DenseStorage storage = MakeTwoClusterStorage();
    const ClusteringResult result = MakeResult({0, 0, 0, 0});

    ClusterReportOptions coverage_nan;
    coverage_nan.coverage_thresholds = {std::numeric_limits<double>::quiet_NaN()};
    try {
        cluster_report(result, storage, coverage_nan);
        FAIL() << "expected std::invalid_argument";
    } catch (const std::invalid_argument& e) {
        EXPECT_EQ(
            std::string(e.what()),
            "cluster_report: coverage threshold 0 must not be NaN");
    }

    ClusterReportOptions boundary_nan;
    boundary_nan.boundary_threshold = std::numeric_limits<double>::quiet_NaN();
    try {
        cluster_report(result, storage, boundary_nan);
        FAIL() << "expected std::invalid_argument";
    } catch (const std::invalid_argument& e) {
        EXPECT_EQ(
            std::string(e.what()),
            "cluster_report: boundary_threshold must not be NaN");
    }
}

// The shape that makes the selector's comparator genuinely inconsistent, rather
// than merely wrong. cluster_representative sorts each candidate's intra
// distances through median_distance and then stable_sorts the candidate scores;
// a NaN alongside two DISTINCT finite values makes operator< a non-strict-weak
// ordering, which is undefined behaviour. Four members with unequal finite
// companions is the smallest cluster that produces that, which three members
// cannot -- so this pins that such an input is refused at all.
//
// It does NOT pin that the refusal precedes the selection: it asserts only the
// eventual message, and stays green under the hoisted ordering too.
// RepresentativeSelectionRunsAfterTheFinitenessCheck is the test that pins the
// order. Both are needed, and neither substitutes for the other.
TEST(ClusterReportTest, NonFiniteIntraDistanceWithInconsistentComparatorIsRefused) {
    DenseStorage storage(4);
    storage.Set(0, 1, 0.10);
    storage.Set(0, 2, 0.30);
    storage.Set(0, 3, std::numeric_limits<double>::quiet_NaN());
    storage.Set(1, 2, 0.50);
    storage.Set(1, 3, 0.20);
    storage.Set(2, 3, 0.40);

    const ClusteringResult result = MakeResult({0, 0, 0, 0});
    try {
        cluster_report(result, storage, ClusterReportOptions());
        FAIL() << "expected std::invalid_argument";
    } catch (const std::invalid_argument& e) {
        const std::string message(e.what());
        EXPECT_EQ(
            message,
            "cluster_report: distance between samples 0 and 3 is not finite")
            << message;
    }
}

namespace {

// Finite everywhere except the FIRST read of edge (0,1), which yields NaN. The
// poison is consumed by whichever caller reads that edge first, which is what
// makes the read order observable without a sanitizer: StorageBackend::Get is
// virtual and detail::checked_distance reads through it, so the checked loop and
// cluster_representative's raw reads both arrive here.
class FirstReadPoisonedStorage : public DenseStorage {
public:
    explicit FirstReadPoisonedStorage(size_t n) : DenseStorage(n) {}

    double Get(size_t i, size_t j) const override {
        const bool poisoned_edge = (i == 0 && j == 1) || (i == 1 && j == 0);
        if (poisoned_edge && !consumed_) {
            consumed_ = true;
            return std::numeric_limits<double>::quiet_NaN();
        }
        return DenseStorage::Get(i, j);
    }

    bool PoisonConsumed() const { return consumed_; }

private:
    mutable bool consumed_ = false;
};

}  // namespace

// The read ORDER is what this pins, not the refusal itself. Both NaN tests above
// assert only the eventual message, and both stay green if the representative
// selection is hoisted above the checked loop -- the selector absorbs the NaN,
// produces some answer, and the refusal still fires later or not at all. Here
// the poison exists for exactly one read, so whichever caller reads d(0,1) first
// is the one that consumes it: with the selection in its current position the
// checked loop takes the poison and refuses, and with the selection hoisted the
// call returns normally. Throw versus return is the discriminator, which is why
// this needs no sanitizer to observe an ordering whose only other symptom is
// undefined behaviour.
TEST(ClusterReportTest, RepresentativeSelectionRunsAfterTheFinitenessCheck) {
    FirstReadPoisonedStorage storage(3);
    storage.Set(0, 1, 0.10);
    storage.Set(0, 2, 0.30);
    storage.Set(1, 2, 0.50);

    const ClusteringResult result = MakeResult({0, 0, 0});
    try {
        cluster_report(result, storage, ClusterReportOptions());
        FAIL() << "expected std::invalid_argument";
    } catch (const std::invalid_argument& e) {
        EXPECT_EQ(
            std::string(e.what()),
            "cluster_report: distance between samples 0 and 1 is not finite");
    }
    // Not the discriminator -- the hoisted ordering consumes the poison too.
    // This is here so the test cannot pass vacuously if a later change stops
    // reading that edge at all.
    EXPECT_TRUE(storage.PoisonConsumed());
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

// Section 7.1 item 14 and the section 5.5 shape rules. The invariant asserted
// in every case is noise_coverage_at.size() == coverage_at.size() -- not
// "always coverage_thresholds.size()", which is false for K == 0.
TEST(ClusterReportTest, NoiseCoverageCurveShapes) {
    // A dedicated fixture, not MakeSixPointStorage. There every point sits
    // within 0.2 of a representative while the smallest default threshold is
    // 0.25, so both curves would be [1, 1, 1] and the EXPECT_NE below could
    // never fire. Here the noise point is deliberately placed exactly on the
    // third threshold (0.45) to make the <= comparison observable and to create
    // a visible step in the curve.
    DenseStorage storage(6);
    storage.Set(0, 1, 0.2);
    storage.Set(0, 2, 0.2);
    storage.Set(1, 2, 0.4);
    storage.Set(3, 4, 0.2);
    storage.Set(3, 5, 0.45);
    storage.Set(4, 5, 0.5);
    for (size_t i = 0; i < 3; ++i) {
        storage.Set(i, 3, 0.8);
        storage.Set(i, 4, 0.8);
        storage.Set(i, 5, 0.9);
    }

    // Clusters plus noise. Representatives are 0 and 3; sample 5 is noise and
    // sits exactly 0.45 from representative 3. The equality at the boundary is
    // exact (same decimal literal, no arithmetic) and pins the <= comparison.
    const ClusterReport with_noise =
        cluster_report(MakeResult({0, 0, 0, 1, 1, -1}), storage, ClusterReportOptions());
    ASSERT_EQ(with_noise.coverage_thresholds,
              (std::vector<double>{0.25, 0.35, 0.45}));
    EXPECT_EQ(with_noise.noise_coverage_at.size(), with_noise.coverage_at.size());
    EXPECT_EQ(with_noise.noise_coverage_at.size(),
              with_noise.coverage_thresholds.size());

    // coverage_at counts all six samples: five are within 0.2 of a
    // representative, and the noise point joins them only at 0.45.
    EXPECT_DOUBLE_EQ(with_noise.coverage_at[0], 5.0 / 6.0);
    EXPECT_DOUBLE_EQ(with_noise.coverage_at[1], 5.0 / 6.0);
    EXPECT_DOUBLE_EQ(with_noise.coverage_at[2], 1.0);
    // noise_coverage_at counts only sample 5.
    EXPECT_DOUBLE_EQ(with_noise.noise_coverage_at[0], 0.0);
    EXPECT_DOUBLE_EQ(with_noise.noise_coverage_at[1], 0.0);
    EXPECT_DOUBLE_EQ(with_noise.noise_coverage_at[2], 1.0);
    EXPECT_NE(with_noise.noise_coverage_at, with_noise.coverage_at);

    // Two noise points at different distances. Representatives are 0 and 3;
    // noise samples 2 and 5 are at 0.2 and 0.45 from their nearest reps. This
    // pins the denominator (num_noise, not num_samples or clustered_count), the
    // per-threshold reset (no accumulation across thresholds), and the inclusive
    // comparison (<= not <, with sample 5 exactly on the 0.45 boundary).
    const ClusterReport two_noise =
        cluster_report(MakeResult({0, 0, -1, 1, 1, -1}), storage, ClusterReportOptions());
    EXPECT_EQ(two_noise.noise_coverage_at.size(), two_noise.coverage_at.size());
    EXPECT_DOUBLE_EQ(two_noise.coverage_at[0], 5.0 / 6.0);
    EXPECT_DOUBLE_EQ(two_noise.coverage_at[1], 5.0 / 6.0);
    EXPECT_DOUBLE_EQ(two_noise.coverage_at[2], 1.0);
    // noise_coverage_at[0] == 0.5 pins the denominator: num_samples would give
    // 1/6, clustered_count would give 1/4. noise_coverage_at[1] == 0.5 pins the
    // per-threshold reset: hoisting covered_noise out of the loop yields 1.0.
    EXPECT_DOUBLE_EQ(two_noise.noise_coverage_at[0], 0.5);
    EXPECT_DOUBLE_EQ(two_noise.noise_coverage_at[1], 0.5);
    EXPECT_DOUBLE_EQ(two_noise.noise_coverage_at[2], 1.0);

    // No noise: full length, every entry NaN. Not 0.0, which would read as
    // "no noise point is covered" rather than "the question does not apply".
    const ClusterReport noise_free =
        cluster_report(MakeResult({0, 0, 0, 1, 1, 1}), storage, ClusterReportOptions());
    EXPECT_EQ(noise_free.noise_coverage_at.size(), noise_free.coverage_at.size());
    ASSERT_EQ(noise_free.noise_coverage_at.size(), 3u);
    for (const double value : noise_free.noise_coverage_at) {
        EXPECT_TRUE(std::isnan(value));
    }

    // All noise is the K == 0 case: there is no representative to measure
    // against, so both curves are empty.
    const ClusterReport all_noise = cluster_report(
        ClusteringResult(std::vector<ClusterLabel>{-1, -1, -1, -1, -1, -1}, Clusters{}),
        storage,
        ClusterReportOptions());
    EXPECT_TRUE(all_noise.coverage_at.empty());
    EXPECT_TRUE(all_noise.noise_coverage_at.empty());
    EXPECT_EQ(all_noise.noise_coverage_at.size(), all_noise.coverage_at.size());
    // The thresholds are still reported, whatever K is.
    EXPECT_EQ(all_noise.coverage_thresholds.size(), 3u);

    const DenseStorage empty_storage(0);
    const ClusterReport zero_samples = cluster_report(
        ClusteringResult(std::vector<ClusterLabel>{}, Clusters{}),
        empty_storage,
        ClusterReportOptions());
    EXPECT_TRUE(zero_samples.coverage_at.empty());
    EXPECT_TRUE(zero_samples.noise_coverage_at.empty());
    EXPECT_EQ(zero_samples.noise_coverage_at.size(), zero_samples.coverage_at.size());

    ClusterReportOptions no_thresholds;
    no_thresholds.coverage_thresholds.clear();
    const ClusterReport unthresholded =
        cluster_report(MakeResult({0, 0, 0, 1, 1, -1}), storage, no_thresholds);
    EXPECT_TRUE(unthresholded.coverage_at.empty());
    EXPECT_TRUE(unthresholded.noise_coverage_at.empty());
    EXPECT_EQ(unthresholded.noise_coverage_at.size(),
              unthresholded.coverage_at.size());
}

// Section 7.1 items 2 and 4. The brute-force reference is the whole point: it
// recomputes gamma by the O(P_w * P_b) definition and c_index by explicitly
// sorting all pooled pairs, so it shares no code path with the run-at-a-time
// walk it checks. Item 4 asks for exactly this on c_index. The oracle shares
// the naive three-sum formulation, so it cannot be used on near-tied inputs.
TEST(InternalIndicesTest, PairRankMatchesBruteForce) {
    const std::vector<double> within{0.2, 0.35, 0.5, 0.1};
    const std::vector<double> between{0.6, 0.15, 0.9, 0.45, 0.7};
    const detail::PairRankIndices got = detail::pair_rank_indices(within, between);
    const detail::PairRankIndices want = BruteForcePairRank(within, between);
    // Gamma is a ratio of two exactly-represented integer counts, so both
    // routes land on the identical double and exact equality is honest.
    EXPECT_DOUBLE_EQ(got.baker_hubert_gamma, want.baker_hubert_gamma);
    // C-index is not, and must not be compared that way. The two routes
    // accumulate different sequences of differences, and the divergence is a
    // property of this particular fixture's summation orders rather than a
    // bound that holds for every input. Measured: 0.18421052631578946 against
    // 0.18421052631578952, two ULP apart and comfortably inside the tolerance.
    // The tolerance is absolute, and the nearest real disagreement is nowhere
    // near it: the smallest single-element misselection this fixture admits --
    // taking 0.45 rather than 0.5 into S_max -- moves c_index to
    // 0.18918918918918914, a shift of 5e-3, some ten orders of magnitude above
    // 1e-12.
    EXPECT_NEAR(got.c_index, want.c_index, 1e-12);
}

// Section 5.4. With no between-pairs the pooled set is `within` itself, so
// S_min and S_max are the same sum in two different orders. They agree only in
// exact arithmetic: the `sum_max != sum_min` test alone lets a rounding
// difference through and reports 0.0 -- the best possible compactness -- for an
// index that is not defined at all.
TEST(InternalIndicesTest, PairRankCIndexIsNaNWithoutBetweenPairs) {
    const std::vector<double> within{0.2, 0.35, 0.5, 0.1};
    const std::vector<double> between;
    const detail::PairRankIndices got = detail::pair_rank_indices(within, between);
    EXPECT_TRUE(std::isnan(got.c_index));
    // No couples exist either, so Gamma is undefined for the same input.
    EXPECT_TRUE(std::isnan(got.baker_hubert_gamma));
}

// Section 7.1 item 2, second fixture. The first one is tie-free, so every run
// has length one and the per-element accumulation is indistinguishable from a
// single addition per run. Here `within` holds a duplicate at 0.3 and `between`
// a duplicate at 0.5, and both runs fall where the opposite array's seen-count
// is already non-zero -- which is what makes the repetition observable.
TEST(InternalIndicesTest, PairRankMatchesBruteForceWithTies) {
    const std::vector<double> within{0.3, 0.3, 0.9};
    const std::vector<double> between{0.1, 0.5, 0.5};
    const detail::PairRankIndices got = detail::pair_rank_indices(within, between);
    const detail::PairRankIndices want = BruteForcePairRank(within, between);
    // Four concordant couples against five discordant.
    EXPECT_DOUBLE_EQ(got.baker_hubert_gamma, -1.0 / 9.0);
    EXPECT_DOUBLE_EQ(got.baker_hubert_gamma, want.baker_hubert_gamma);
    EXPECT_NEAR(got.c_index, want.c_index, 1e-12);
}

// Section 7.1 item 4, exhaustion arm. Both merge walks take within_count
// elements from the pooled set, and each needs a guard for the case where
// `between` runs out first -- without the descending one, the comparison
// `between[high_between - 1]` at zero underflows to SIZE_MAX and reads out of
// bounds. No other fixture reaches either guard, but not for one reason: those
// that assert a c_index all have at least as many between-pairs as
// within-pairs, and the one with no between-pairs at all is turned away by the
// between_count > 0 gate before either walk runs. Here the single
// between-distance sits in the middle of the pooled range with `within` the
// larger array, so the ascending walk exhausts `between` from the bottom and
// the descending walk exhausts it from the top.
TEST(InternalIndicesTest, PairRankCIndexExhaustsBetweenInBothWalks) {
    const std::vector<double> within{0.1, 0.2, 0.8, 0.9};
    const std::vector<double> between{0.5};
    const detail::PairRankIndices got = detail::pair_rank_indices(within, between);
    const detail::PairRankIndices want = BruteForcePairRank(within, between);
    // S_w = 2.0, S_min = 1.6, S_max = 2.4, so the index is one half.
    EXPECT_NEAR(got.c_index, 0.5, 1e-12);
    EXPECT_NEAR(got.c_index, want.c_index, 1e-12);
}

// Section 5.4, near-tied arm. S_w, S_min and S_max agree to within one ULP
// here, so a c_index built by differencing them afterwards loses the entire
// quantity it is trying to measure: exact S_min is 2 - 2^-53, which is the
// midpoint between two doubles and rounds to the same 2.0 as S_max, and the
// index came out NaN instead of 1. Fingerprint distances are frequently
// near-constant, so this is the shape of a real dataset and not only of an
// adversarial one. BruteForcePairRank is deliberately not the oracle: it sums
// the same way and returns NaN on this input too.
TEST(InternalIndicesTest, PairRankCIndexSurvivesNearTiedDistances) {
    const double just_below_one = std::nextafter(1.0, 0.0);
    const std::vector<double> within{1.0, 1.0};
    const std::vector<double> between{just_below_one, 1.0, 1.0, 1.0};
    const detail::PairRankIndices got = detail::pair_rank_indices(within, between);
    // Every within-distance is at the top of the pooled range, so the index is
    // at its worst-case value, and the gap formulation reaches it exactly.
    EXPECT_DOUBLE_EQ(got.c_index, 1.0);
}

// Section 5.4, lower endpoint. Every within-distance is below every
// between-distance, so the within-pairs already are the smallest pooled pairs
// and the index is 0 -- a real, defined, best-case score. It is the only
// fixture that separates "the denominator is zero" from "the numerator is
// zero", and reporting NaN for a perfectly separated clustering would be
// exactly the over-refusal this branch treats as a defect.
TEST(InternalIndicesTest, PairRankCIndexIsZeroForPerfectSeparation) {
    const std::vector<double> within{0.1, 0.2};
    const std::vector<double> between{0.8, 0.9};
    const detail::PairRankIndices got = detail::pair_rank_indices(within, between);
    EXPECT_FALSE(std::isnan(got.c_index));
    EXPECT_DOUBLE_EQ(got.c_index, 0.0);
}

// Section 7.1 item 3. A couple whose two distances are equal contributes to
// neither counter; the run-at-a-time walk is what makes that true. The fixture
// needs a discordant couple as well as a concordant one -- with none, scoring
// the tie as concordant leaves the ratio at 1.0 and the test observes nothing.
TEST(InternalIndicesTest, PairRankExcludesTiedCouples) {
    const std::vector<double> within{0.5, 0.2};
    const std::vector<double> between{0.5, 0.9, 0.1};
    const detail::PairRankIndices got = detail::pair_rank_indices(within, between);
    // Three concordant couples, two discordant, and one tie that scores as
    // neither, giving (3 - 2) / 5. Counting the tie as concordant would make it
    // (4 - 2) / 6 = 1/3 instead.
    EXPECT_DOUBLE_EQ(got.baker_hubert_gamma, 0.2);
    EXPECT_DOUBLE_EQ(got.baker_hubert_gamma,
                     BruteForcePairRank(within, between).baker_hubert_gamma);
}

// Section 7.1 item 16. An unsigned subtraction wraps here and returns a value
// near +1, scoring the worst clustering as the best one. The fixture is
// deliberately interior rather than saturated: at -1.0 the concordant count is
// zero, and a test sitting at the end of the range cannot see an error that
// pushes it further that way.
TEST(InternalIndicesTest, PairRankGammaGoesNegative) {
    const std::vector<double> within{0.9, 0.4};
    const std::vector<double> between{0.1, 0.1, 0.5, 0.8};
    const detail::PairRankIndices got = detail::pair_rank_indices(within, between);
    // Two concordant couples against six discordant.
    EXPECT_LT(got.baker_hubert_gamma, 0.0);
    EXPECT_DOUBLE_EQ(got.baker_hubert_gamma, -0.5);
    EXPECT_DOUBLE_EQ(got.baker_hubert_gamma,
                     BruteForcePairRank(within, between).baker_hubert_gamma);
}

TEST(InternalIndicesTest, PairRankAllTiedIsNaNNotARefusal) {
    const std::vector<double> within{0.5, 0.5};
    const std::vector<double> between{0.5, 0.5, 0.5, 0.5};
    detail::PairRankIndices got{0.0, 0.0};
    EXPECT_NO_THROW(got = detail::pair_rank_indices(within, between));
    EXPECT_TRUE(std::isnan(got.baker_hubert_gamma));
    EXPECT_TRUE(std::isnan(got.c_index));
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

// Section 5.3. Pins that a merged pair of streams and one accumulated stream
// reach the same count, mean and standard deviation, and that both match values
// computed independently.
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
    // Both routes are also checked against the values themselves, not only
    // against each other: a self-comparison stays green if the shared
    // implementation is wrong in the same way on both sides.
    EXPECT_DOUBLE_EQ(single.mean, 4.0 / 7.0);
    EXPECT_NEAR(detail::population_stddev(single), 0.27105237087157541, 1e-12);
}

TEST(InternalIndicesTest, PopulationStddevIsZeroBelowTwoValues) {
    const detail::DistanceMoments empty;
    EXPECT_DOUBLE_EQ(detail::population_stddev(empty), 0.0);

    detail::DistanceMoments one;
    one.Add(0.42);
    EXPECT_DOUBLE_EQ(detail::population_stddev(one), 0.0);
}

// Section 5.3, the reason DistanceMoments is Welford at all. On a
// near-constant stream -- the common case for fingerprint distances -- the
// naive sqrt(sum_sq/n - mean^2) subtracts two quantities that agree to more
// digits than a double carries. Measured on this fixture its radicand comes out
// at -6.25e-18 and the result is NaN, where Welford lands on the right answer.
// WelfordMergeMatchesSinglePass cannot see that: its fixture is well
// conditioned enough that both formulations agree to 5.6e-17.
TEST(InternalIndicesTest, PopulationStddevSurvivesNearConstantValues) {
    detail::DistanceMoments moments;
    for (const double value : {1.0 - 1e-8, 1.0, 1.0, 1.0}) {
        moments.Add(value);
    }
    // The closed form is sqrt(3)/4 * 1e-8 = 4.330127018922193e-09; the computed
    // value differs from it at 1.8e-17 because 1 - 1e-8 is not exact in binary,
    // so the assertion is against the measured double. The tolerance is far
    // above any reassociation or FMA-contraction difference and far below the
    // gap to any wrong answer -- the naive form does not miss by a little here,
    // it returns NaN.
    EXPECT_NEAR(detail::population_stddev(moments), 4.3301270006183164e-09, 1e-16);
}

// Section 5.3. The short-circuits exist because Chan's update divides by the
// combined count and would produce NaN on an empty operand, and cluster_report
// merges per-cluster streams of which some can legitimately be empty. Both
// orders are checked: returning the wrong operand is the natural way to write
// this wrong, and it is invisible unless one side is empty.
TEST(InternalIndicesTest, MergeMomentsPassesThroughEmptyStreams) {
    detail::DistanceMoments filled;
    filled.Add(0.25);
    filled.Add(0.75);
    const detail::DistanceMoments empty;

    const detail::DistanceMoments left = detail::merge_moments(empty, filled);
    EXPECT_EQ(left.count, 2u);
    EXPECT_DOUBLE_EQ(left.mean, 0.5);
    EXPECT_DOUBLE_EQ(left.m2, 0.125);

    const detail::DistanceMoments right = detail::merge_moments(filled, empty);
    EXPECT_EQ(right.count, 2u);
    EXPECT_DOUBLE_EQ(right.mean, 0.5);
    EXPECT_DOUBLE_EQ(right.m2, 0.125);

    const detail::DistanceMoments neither = detail::merge_moments(empty, empty);
    EXPECT_EQ(neither.count, 0u);
    EXPECT_DOUBLE_EQ(detail::population_stddev(neither), 0.0);
}

// Section 5.2. The refusal is the precondition every later index depends on:
// a NaN distance gives std::sort no strict weak ordering, so the pair-rank
// sorts would be undefined behaviour rather than merely wrong.
TEST(InternalIndicesTest, CheckedDistanceRefusesNonFinite) {
    DenseStorage storage(3);
    storage.Set(0, 1, 0.25);
    storage.Set(0, 2, std::numeric_limits<double>::quiet_NaN());
    storage.Set(1, 2, std::numeric_limits<double>::infinity());

    EXPECT_DOUBLE_EQ(detail::checked_distance(storage, 0, 1), 0.25);
    EXPECT_THROW(detail::checked_distance(storage, 0, 2), std::invalid_argument);
    EXPECT_THROW(detail::checked_distance(storage, 1, 2), std::invalid_argument);
}

// Section 7.1 item 10. The records must decompose the scalars, not merely
// resemble them.
//
// The three-cluster split and the raised boundary_threshold are both
// deliberate. At the default 0.30 the per-cluster counts are 1/2/1 (d(0,2) and
// d(3,4) both trip at 0.2), so clusters A and C tie and summing the wrong
// cluster's total could still pass. At 0.5 the third pair d(1,2) = 0.4 trips
// as well, giving 2/3/1, all distinct, so summing the wrong cluster's total is
// caught.
TEST(ClusterReportTest, RecordsRecomputeTheAggregates) {
    const DenseStorage storage = MakeSixPointStorage();
    ClusterReportOptions options;
    options.compute_per_cluster_records = true;
    options.boundary_threshold = 0.5;
    const ClusterReport r =
        cluster_report(MakeResult({0, 0, 1, 1, 2, -1}), storage, options);

    ASSERT_EQ(r.records.size(), r.num_clusters);
    ASSERT_EQ(r.records.size(), 3u);
    EXPECT_TRUE(r.requested.per_cluster_records);

    std::vector<double> radii;
    std::vector<double> representative_means;
    size_t record_violations = 0;
    for (size_t k = 0; k < r.records.size(); ++k) {
        EXPECT_EQ(r.records[k].label, static_cast<ClusterLabel>(k));
        radii.push_back(r.records[k].radius);
        representative_means.push_back(r.records[k].mean_representative_distance);
        record_violations += r.records[k].boundary_violations;
    }

    std::sort(radii.begin(), radii.end());
    std::sort(representative_means.begin(), representative_means.end());
    EXPECT_DOUBLE_EQ(detail::median_distance(radii), r.median_radius);
    EXPECT_DOUBLE_EQ(detail::median_distance(representative_means),
                     r.median_medoid_member_distance);

    // Pinned, not just self-consistent: {0,1}-{2,3} contributes d(0,2)=0.2 and
    // d(1,2)=0.4; {2,3}-{4} contributes d(3,4)=0.2; {0,1}-{4} contributes
    // nothing. Three violating pairs in total.
    ASSERT_EQ(r.boundary_violations, 3u);
    EXPECT_EQ(r.records[0].boundary_violations, 2u);
    EXPECT_EQ(r.records[1].boundary_violations, 3u);
    EXPECT_EQ(r.records[2].boundary_violations, 1u);

    // Each violating pair is counted by both of its endpoints on the records,
    // and once on the scorecard.
    EXPECT_EQ(record_violations, 2u * r.boundary_violations);
}

// The per-record violation count of cluster b takes one contribution per
// lower-ordinal partner, so only a cluster reached by two POSITIVE
// contributions separates an accumulation from an overwrite.
// RecordsRecomputeTheAggregates cannot do that: at its 0.5 threshold cluster C
// is handed 0 and then 1, and discarding a zero changes nothing. Raising the
// threshold to 0.9 makes every cross pair trip, so C is handed 2 and 2.
TEST(ClusterReportTest, RecordBoundaryViolationsAccumulateAcrossPartners) {
    const DenseStorage storage = MakeSixPointStorage();
    ClusterReportOptions options;
    options.compute_per_cluster_records = true;
    options.boundary_threshold = 0.9;
    const ClusterReport r =
        cluster_report(MakeResult({0, 0, 1, 1, 2, -1}), storage, options);

    ASSERT_EQ(r.records.size(), 3u);

    // Every cross pair among the clustered points trips at 0.9: A-B has four
    // (0.2, 0.8, 0.4, 0.8), A-C two (0.8, 0.8) and B-C two (0.8, 0.2).
    ASSERT_EQ(r.boundary_violations, 8u);

    // C is the only discriminator. It takes 2 from A-C and 2 from B-C, so an
    // overwrite would leave it at 2. A is never a later partner, and B's two
    // contributions reach it through different slots, so both land on 6
    // either way.
    EXPECT_EQ(r.records[0].boundary_violations, 6u);
    EXPECT_EQ(r.records[1].boundary_violations, 6u);
    EXPECT_EQ(r.records[2].boundary_violations, 4u);

    size_t record_violations = 0;
    for (const ClusterRecord& record : r.records) {
        record_violations += record.boundary_violations;
    }
    EXPECT_EQ(record_violations, 2u * r.boundary_violations);
}

// Section 7.1 item 10 also names p95_diameter. percentile() is file-local to
// ClusterReport.cpp and cannot be called from here, so the expected value is
// worked out by hand instead -- which means the three diameters must differ.
// A fixture where they are all equal would pass against first, min, median and
// max alike and would test nothing about the percentile.
TEST(ClusterReportTest, RecordDiametersRecomputeP95) {
    const DenseStorage storage = MakeSixPointStorage();
    ClusterReportOptions options;
    options.compute_per_cluster_records = true;
    const ClusterReport r =
        cluster_report(MakeResult({0, 0, 1, 1, 2, -1}), storage, options);

    ASSERT_EQ(r.records.size(), 3u);
    EXPECT_DOUBLE_EQ(r.records[0].diameter, 0.2);  // d(0,1)
    EXPECT_DOUBLE_EQ(r.records[1].diameter, 0.8);  // d(2,3), a cross-group pair
    EXPECT_DOUBLE_EQ(r.records[2].diameter, 0.0);  // singleton

    // Also the only fixture that can catch a per-cluster median sourced from
    // the cumulative pool. intra_pairs accumulates across clusters while
    // cluster_distances is cleared per cluster (ClusterReport.cpp:331), so the
    // two agree on the FIRST cluster and diverge afterwards -- every other
    // record-median assertion in this file reads records[0] or a singleton.
    // Reading B's median is what discriminates: its own single pair is 0.8,
    // while the pool at that point is {0.2, 0.8} and medians to 0.5.
    EXPECT_DOUBLE_EQ(r.records[0].median_intra_distance, 0.2);  // d(0,1)
    EXPECT_DOUBLE_EQ(r.records[1].median_intra_distance, 0.8);  // d(2,3)
    // The guard is !cluster_distances.empty(), independent of the value's
    // source, so this stays NaN under that mutation -- pinned so a future
    // change that folds guard and source together cannot pass silently.
    EXPECT_TRUE(std::isnan(r.records[2].median_intra_distance));
    // The same discrimination for the mean. A guard-preserving mutation --
    // replacing only the ternary's defined branch with report.mean_intra_distance
    // -- keeps every singleton's NaN, so the isnan pin at :2359 cannot see it.
    // Only a populated cluster whose own mean differs from the global mean can:
    // here the locals are 0.2 and 0.8 while the pool means 0.5.
    EXPECT_DOUBLE_EQ(r.records[0].mean_intra_distance, 0.2);  // d(0,1)
    EXPECT_DOUBLE_EQ(r.records[1].mean_intra_distance, 0.8);  // d(2,3)
    EXPECT_TRUE(std::isnan(r.records[2].mean_intra_distance));

    // Sorted diameters [0.0, 0.2, 0.8]; fractional rank 0.95 * 2 = 1.9, so the
    // result interpolates 90% of the way from 0.2 to 0.8: 0.2 + 0.9 * 0.6.
    // Distinct from the first (0.2), the min (0.0), the median (0.2) and the
    // max (0.8), so a wrong reduction cannot pass.
    EXPECT_NEAR(r.p95_diameter, 0.74, 1e-12);
}

// The record-level half of Task 6's MedoidNamedIndicesIgnoreRepresentativeMethod:
// ClusterRecord::representative is the CONFIGURED representative, section 3.7,
// while the medoid-named scalars are not. Task 6 could not assert this because
// the table did not exist yet.
TEST(ClusterReportTest, RecordRepresentativeFollowsTheConfiguredMethod) {
    const DenseStorage storage = MakeDivergentRepresentativeStorage();
    const ClusteringResult result = MakeResult({0, 0, 0, 0, 1, 1, 1, 1});

    ClusterReportOptions medoid_options;
    medoid_options.representative_method = RepresentativeMethod::Medoid;
    medoid_options.compute_per_cluster_records = true;
    ClusterReportOptions minimax_options;
    minimax_options.representative_method = RepresentativeMethod::Minimax;
    minimax_options.compute_per_cluster_records = true;

    const ClusterReport by_medoid = cluster_report(result, storage, medoid_options);
    const ClusterReport by_minimax = cluster_report(result, storage, minimax_options);

    ASSERT_EQ(by_medoid.records.size(), 2u);
    ASSERT_EQ(by_minimax.records.size(), 2u);
    EXPECT_EQ(by_medoid.records[0].representative, 0u);
    EXPECT_EQ(by_medoid.records[1].representative, 4u);
    EXPECT_EQ(by_minimax.records[0].representative, 1u);
    EXPECT_EQ(by_minimax.records[1].representative, 5u);

    EXPECT_NEAR(by_medoid.records[0].mean_representative_distance, 0.3, 1e-12);
    EXPECT_DOUBLE_EQ(by_medoid.records[0].radius, 0.7);
    EXPECT_NEAR(by_minimax.records[0].mean_representative_distance, 1.1 / 3.0, 1e-12);
    EXPECT_DOUBLE_EQ(by_minimax.records[0].radius, 0.5);

    // The two fields that must NOT move with the configured representative.
    EXPECT_NEAR(by_medoid.records[0].mean_intra_distance, 0.4, 1e-12);
    EXPECT_DOUBLE_EQ(by_medoid.records[0].median_intra_distance, 0.5);
    EXPECT_NEAR(by_minimax.records[0].mean_intra_distance, 0.4, 1e-12);
    EXPECT_DOUBLE_EQ(by_minimax.records[0].median_intra_distance, 0.5);

    // And the medoid-named scalars still did not move with them.
    EXPECT_DOUBLE_EQ(by_medoid.davies_bouldin_medoid,
                     by_minimax.davies_bouldin_medoid);
}

// Section 7.1 item 12, against the section 4.4 table.
TEST(ClusterReportTest, RecordUndefinedCases) {
    const DenseStorage storage = MakeSixPointStorage();
    ClusterReportOptions options;
    options.compute_per_cluster_records = true;

    const ClusterReport single =
        cluster_report(MakeResult({0, 0, 0, -1, -1, -1}), storage, options);
    ASSERT_EQ(single.records.size(), 1u);
    EXPECT_EQ(single.records[0].nearest_cluster, NO_NEAREST_CLUSTER);
    EXPECT_TRUE(std::isnan(single.records[0].nearest_cluster_distance));
    EXPECT_NE(single.records[0].nearest_cluster_distance, 0.0);
    EXPECT_TRUE(std::isnan(single.records[0].silhouette));

    const ClusterReport with_singleton =
        cluster_report(MakeResult({0, 0, 0, 1, -1, -1}), storage, options);
    ASSERT_EQ(with_singleton.records.size(), 2u);
    const ClusterRecord& lone = with_singleton.records[1];
    EXPECT_EQ(lone.size, 1u);
    EXPECT_TRUE(std::isnan(lone.mean_intra_distance));
    EXPECT_TRUE(std::isnan(lone.median_intra_distance));
    EXPECT_DOUBLE_EQ(lone.radius, 0.0);
    EXPECT_DOUBLE_EQ(lone.diameter, 0.0);
    EXPECT_DOUBLE_EQ(lone.mean_representative_distance, 0.0);
    EXPECT_EQ(lone.nearest_cluster, 0);
    EXPECT_DOUBLE_EQ(lone.nearest_cluster_distance, 0.8);
    EXPECT_NEAR(lone.silhouette, 1.0, 1e-12);

    // The defined cases, which are what make the NaN assertions above mean
    // anything. Since Task 2's gate these four fields DEFAULT to NaN, so an
    // isnan() assertion alone passes just as happily when the population path
    // never touched the field. Asserting a real value on a populated record is
    // the half that catches a missing assignment.
    const ClusterRecord& populated = with_singleton.records[0];
    EXPECT_EQ(populated.size, 3u);
    EXPECT_NEAR(populated.mean_intra_distance, 0.8 / 3.0, 1e-12);
    EXPECT_DOUBLE_EQ(populated.median_intra_distance, 0.2);
    EXPECT_DOUBLE_EQ(populated.radius, 0.2);
    EXPECT_DOUBLE_EQ(populated.diameter, 0.4);
    EXPECT_DOUBLE_EQ(populated.mean_representative_distance, 0.2);
    EXPECT_DOUBLE_EQ(populated.nearest_cluster_distance, 0.8);
    EXPECT_NEAR(populated.silhouette, 2.0 / 3.0, 1e-12);
    EXPECT_EQ(populated.nearest_cluster, 1);
}

// Section 7.1 item 13, the per_cluster_records half. This is what the requested
// struct exists for: an empty table means two different things.
TEST(ClusterReportTest, RecordsHonestyProperty) {
    const DenseStorage storage = MakeSixPointStorage();

    const ClusterReport unrequested =
        cluster_report(MakeResult({0, 0, 0, 1, 1, 1}), storage, ClusterReportOptions());
    EXPECT_TRUE(unrequested.records.empty());
    EXPECT_FALSE(unrequested.requested.per_cluster_records);

    ClusterReportOptions options;
    options.compute_per_cluster_records = true;
    const ClusterReport requested_but_empty = cluster_report(
        ClusteringResult(std::vector<ClusterLabel>{-1, -1, -1, -1, -1, -1}, Clusters{}),
        storage,
        options);
    EXPECT_TRUE(requested_but_empty.records.empty());
    EXPECT_TRUE(requested_but_empty.requested.per_cluster_records);

    ClusterReportOptions pair_rank_options;
    pair_rank_options.compute_pair_rank_indices = true;
    const ClusterReport with_pair_rank =
        cluster_report(MakeResult({0, 0, 0, 1, 1, 1}), storage, pair_rank_options);
    EXPECT_TRUE(with_pair_rank.requested.pair_rank_indices);
}

// The nearest_cluster tiebreak rule ("lowest ordinal", ClusterReport.h:76) and
// the per-cluster nearest_cluster_distance field (distinct from the global
// min_inter). MakeSixPointStorage with three clusters A={0,1,2}, B={3,4}, C={5}.
TEST(ClusterReportTest, RecordNearestClusterTiesAndPerClusterDistance) {
    const DenseStorage storage = MakeSixPointStorage();
    ClusterReportOptions options;
    options.compute_per_cluster_records = true;
    const ClusterReport r =
        cluster_report(MakeResult({0, 0, 0, 1, 1, 2}), storage, options);

    ASSERT_EQ(r.records.size(), 3u);

    // A is tied: minima to B and to C are both 0.8, and B is visited first in
    // the cross loop, so the strict < keeps ordinal 1. Under <= it becomes 2.
    EXPECT_EQ(r.records[0].nearest_cluster, 1);
    // A's distance differs from min_inter (0.2). Under the global-min mutation
    // record 0 would report 0.2 while its own nearest_cluster field names B.
    EXPECT_DOUBLE_EQ(r.records[0].nearest_cluster_distance, 0.8);

    EXPECT_EQ(r.records[1].nearest_cluster, 2);
    EXPECT_DOUBLE_EQ(r.records[1].nearest_cluster_distance, 0.2);

    // C is a singleton inside a K == 3 clustering. Under the size-guard mutation
    // its three guard-controlled fields go undefined.
    EXPECT_EQ(r.records[2].nearest_cluster, 1);
    EXPECT_DOUBLE_EQ(r.records[2].nearest_cluster_distance, 0.2);
    EXPECT_NEAR(r.records[2].silhouette, 1.0, 1e-12);
}

// The a-side nearest-cluster write is a running minimum, and the tie fixture
// above cannot see that: A's two candidates are both 0.8, so a guard that
// accepts any different value rejects the second one for the same reason the
// shipped strict < does. Under {0, 0, 1, 1, 2, -1} the candidates rise --
// cluster A is offered 0.2 from B and then 0.8 from C -- so a != in place of
// the < keeps the last distinct partner and names C at 0.8 as A's nearest
// cluster. That answer is in range, is not a NaN, and nothing else in the
// report contradicts it.
TEST(ClusterReportTest, RecordNearestClusterOnTheASideKeepsTheSmallestNotTheLast) {
    const DenseStorage storage = MakeSixPointStorage();
    ClusterReportOptions options;
    options.compute_per_cluster_records = true;
    const ClusterReport r =
        cluster_report(MakeResult({0, 0, 1, 1, 2, -1}), storage, options);

    ASSERT_EQ(r.records.size(), 3u);

    // The discriminating pair: A's later candidate is the farther one.
    EXPECT_EQ(r.records[0].nearest_cluster, 1);
    EXPECT_DOUBLE_EQ(r.records[0].nearest_cluster_distance, 0.2);

    // Asserted so a failure separates the rising-sequence defect from a
    // wholesale change: B ties at 0.2 between A and C and keeps the lowest
    // ordinal, and C is decided outright.
    EXPECT_EQ(r.records[1].nearest_cluster, 0);
    EXPECT_DOUBLE_EQ(r.records[1].nearest_cluster_distance, 0.2);
    EXPECT_EQ(r.records[2].nearest_cluster, 1);
    EXPECT_DOUBLE_EQ(r.records[2].nearest_cluster_distance, 0.2);
}

// The section 7.1 item 1 fixture again: a perfect clustering has C-index 0 and
// Gamma 1, both exactly.
TEST(ClusterReportTest, PairRankIndicesOnTheHandComputedFixture) {
    const DenseStorage storage = MakeSixPointStorage();
    ClusterReportOptions options;
    options.compute_pair_rank_indices = true;
    const ClusterReport r =
        cluster_report(MakeResult({0, 0, 0, 1, 1, 1}), storage, options);
    EXPECT_TRUE(r.requested.pair_rank_indices);
    EXPECT_NEAR(r.c_index, 0.0, 1e-12);
    EXPECT_NEAR(r.baker_hubert_gamma, 1.0, 1e-12);
}

// Section 7.1 item 13, the pair_rank half.
TEST(ClusterReportTest, PairRankHonestyProperty) {
    const DenseStorage storage = MakeSixPointStorage();

    const ClusterReport unrequested =
        cluster_report(MakeResult({0, 0, 0, 1, 1, 1}), storage, ClusterReportOptions());
    EXPECT_TRUE(std::isnan(unrequested.c_index));
    EXPECT_TRUE(std::isnan(unrequested.baker_hubert_gamma));
    EXPECT_FALSE(unrequested.requested.pair_rank_indices);

    ClusterReportOptions options;
    options.compute_pair_rank_indices = true;
    const ClusterReport single =
        cluster_report(MakeResult({0, 0, 0, -1, -1, -1}), storage, options);
    EXPECT_TRUE(std::isnan(single.c_index));
    EXPECT_TRUE(std::isnan(single.baker_hubert_gamma));
    EXPECT_TRUE(single.requested.pair_rank_indices);
}

// Section 7.1 item 18. A tie-heavy clustering is answered, not refused: s+ and
// s- both stay 0 and Gamma is the section 5.5 NaN. The discarded pre-walk guard
// on P_w * P_b would have refused this shape at scale, but this fixture cannot
// discriminate between the two guards -- at any size a test can reach,
// P_w * P_b sits many orders of magnitude below 2^64.
TEST(ClusterReportTest, TieHeavyClusteringIsAnsweredNotRefused) {
    DenseStorage storage(4);
    for (size_t i = 0; i < 4; ++i) {
        for (size_t j = i + 1; j < 4; ++j) {
            storage.Set(i, j, 0.5);
        }
    }
    ClusterReportOptions options;
    options.compute_pair_rank_indices = true;
    ClusterReport r;
    EXPECT_NO_THROW(r = cluster_report(MakeResult({0, 0, 1, 1}), storage, options));
    EXPECT_TRUE(std::isnan(r.baker_hubert_gamma));
    // All ties also put S_max == S_min, so the C-index is the same 0/0.
    EXPECT_TRUE(std::isnan(r.c_index));
}

// Section 7.1 item 16 at the report level.
TEST(ClusterReportTest, ReportGammaGoesNegativeOnABadClustering) {
    DenseStorage storage(4);
    storage.Set(0, 1, 0.9);
    storage.Set(2, 3, 0.9);
    storage.Set(0, 2, 0.1);
    storage.Set(0, 3, 0.1);
    storage.Set(1, 2, 0.1);
    storage.Set(1, 3, 0.1);
    ClusterReportOptions options;
    options.compute_pair_rank_indices = true;
    const ClusterReport r = cluster_report(MakeResult({0, 0, 1, 1}), storage, options);
    EXPECT_NEAR(r.baker_hubert_gamma, -1.0, 1e-12);
}

// Every other report-level c_index assertion sits at 0.0, the value a perfect
// clustering takes, and zero is a fixed point of any distortion that maps zero
// to zero -- squaring the result, scaling it, or clamping it low all leave the
// perfect fixture and the NaN cases green. Only a fixture whose correct index
// falls strictly inside (0, 1) pins the wiring to the value the helper
// returned. This clustering is deliberately a poor one: both within-distances
// sit inside the between-range, which is what lifts the index off its
// endpoint.
TEST(ClusterReportTest, ReportCIndexIsPinnedStrictlyInsideItsRange) {
    DenseStorage storage(4);
    storage.Set(0, 1, 0.4);
    storage.Set(2, 3, 0.6);
    storage.Set(0, 2, 0.1);
    storage.Set(0, 3, 0.2);
    storage.Set(1, 2, 0.8);
    storage.Set(1, 3, 0.9);

    ClusterReportOptions options;
    options.compute_pair_rank_indices = true;
    const ClusterReport r =
        cluster_report(MakeResult({0, 0, 1, 1}), storage, options);

    // Pooled ascending: 0.1 0.2 0.4 0.6 0.8 0.9. The two within-pairs total
    // S_w = 1.0, the two smallest pooled pairs S_min = 0.3, the two largest
    // S_max = 1.7, so the index is (1.0 - 0.3) / (1.7 - 0.3) = 1/2.
    EXPECT_NEAR(r.c_index, 0.5, 1e-12);

    // Gamma is 0 here for a reason worth separating from the tie-heavy case:
    // each within-distance beats two between-distances and loses to two, so
    // s+ and s- are both 4. That is a defined 0.0, not the 0/0 that yields NaN.
    EXPECT_FALSE(std::isnan(r.baker_hubert_gamma));
    EXPECT_NEAR(r.baker_hubert_gamma, 0.0, 1e-12);
}

// The report-level Gamma assertions elsewhere all sit on -1, 0 or 1, and each
// of those is a fixed point of any rounding or sign-preserving snap applied to
// the wiring at the report boundary. A fixture whose correct Gamma falls
// strictly between an endpoint and zero is what makes that class of distortion
// observable. The C-index here is also deliberately count-sensitive: the
// between-set is asymmetric, so duplicating every between-pair moves S_min and
// S_max and changes the answer, which a symmetric fixture does not.
TEST(ClusterReportTest, ReportPairRankPinsInteriorGammaAndPairCounts) {
    DenseStorage storage(4);
    storage.Set(0, 1, 0.3);
    storage.Set(2, 3, 0.4);
    storage.Set(0, 2, 0.1);
    storage.Set(0, 3, 0.5);
    storage.Set(1, 2, 0.6);
    storage.Set(1, 3, 0.7);

    ClusterReportOptions options;
    options.compute_pair_rank_indices = true;
    const ClusterReport r =
        cluster_report(MakeResult({0, 0, 1, 1}), storage, options);

    // Within {0.3, 0.4}, between {0.1, 0.5, 0.6, 0.7}. Pooled ascending:
    // 0.1 0.3 0.4 0.5 0.6 0.7. S_w = 0.7, the two smallest pooled sum to
    // S_min = 0.4 and the two largest to S_max = 1.3, so the index is
    // (0.7 - 0.4) / (1.3 - 0.4) = 1/3.
    EXPECT_NEAR(r.c_index, 1.0 / 3.0, 1e-12);

    // Each within-distance beats one between-distance (0.1) and loses to the
    // other three, so s+ = 6 and s- = 2, giving (6 - 2) / 8 = 0.5. Strictly
    // interior on both sides: not an endpoint, and not zero.
    EXPECT_NEAR(r.baker_hubert_gamma, 0.5, 1e-12);
}

// Section 5.4 requires the within and between arrays to partition every one of
// the Nc(Nc-1)/2 clustered pairs, and a singleton cluster contributes no
// within-pair but every one of its cross pairs. No other pair-rank fixture has
// a cluster smaller than two, so a rule that skipped singletons when collecting
// between-distances would leave the whole suite green. Both expected values
// move under such a rule, which is what makes this fixture discriminating
// rather than merely present.
TEST(ClusterReportTest, ReportPairRankCountsSingletonClusterPairs) {
    DenseStorage storage(5);
    storage.Set(0, 1, 0.3);
    storage.Set(0, 2, 0.5);
    storage.Set(0, 3, 0.6);
    storage.Set(0, 4, 0.1);
    storage.Set(1, 2, 0.7);
    storage.Set(1, 3, 0.8);
    storage.Set(1, 4, 0.2);
    storage.Set(2, 3, 0.4);
    storage.Set(2, 4, 0.9);
    storage.Set(3, 4, 1.0);

    ClusterReportOptions options;
    options.compute_pair_rank_indices = true;
    const ClusterReport r =
        cluster_report(MakeResult({0, 0, 1, 1, 2}), storage, options);

    // Clusters {0,1}, {2,3} and the singleton {4}. Within {0.3, 0.4}; the
    // eight between-distances are every remaining pair. Pooled ascending:
    // 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0. S_w = 0.7, S_min = 0.3,
    // S_max = 1.9, so the index is (0.7 - 0.3) / (1.9 - 0.3) = 0.25.
    EXPECT_NEAR(r.c_index, 0.25, 1e-12);

    // Each within-distance beats 0.1 and 0.2 and loses to the other six, so
    // s+ = 12 and s- = 4, giving (12 - 4) / 16 = 0.5.
    EXPECT_NEAR(r.baker_hubert_gamma, 0.5, 1e-12);
}

// The singleton fixture above has a single singleton, so every one of its
// cross pairs touches a cluster of size two or more. That leaves the
// singleton-to-singleton pair unreached: a rule collecting a cross pair when
// either side is a non-singleton would drop only those pairs, satisfy every
// other pair-rank assertion, and still violate the section 5.4 partition. A
// {2,1,1} shape is the smallest clustering that contains such a pair, and both
// metrics here move when it goes missing.
TEST(ClusterReportTest, ReportPairRankCountsPairsBetweenTwoSingletons) {
    DenseStorage storage(4);
    storage.Set(0, 1, 0.4);
    storage.Set(0, 2, 0.1);
    storage.Set(1, 2, 0.5);
    storage.Set(0, 3, 0.6);
    storage.Set(1, 3, 0.7);
    storage.Set(2, 3, 0.9);

    ClusterReportOptions options;
    options.compute_pair_rank_indices = true;
    const ClusterReport r =
        cluster_report(MakeResult({0, 0, 1, 2}), storage, options);

    // Clusters {0,1}, {2} and {3}. One within-distance, 0.4; five between,
    // {0.1, 0.5, 0.6, 0.7, 0.9}, of which 0.9 is the pair joining the two
    // singletons. Pooled ascending: 0.1 0.4 0.5 0.6 0.7 0.9. With one
    // within-pair, S_w = 0.4, S_min = 0.1 and S_max = 0.9, so the index is
    // (0.4 - 0.1) / (0.9 - 0.1) = 0.375. Dropping 0.9 would make 0.7 the
    // largest pooled distance and the index 0.5.
    EXPECT_NEAR(r.c_index, 0.375, 1e-12);

    // The single within-distance loses to 0.1 and beats the other four, so
    // s+ = 4 and s- = 1, giving (4 - 1) / 5 = 0.6. Dropping 0.9 would leave
    // s+ = 3 against s- = 1, and 0.5.
    EXPECT_NEAR(r.baker_hubert_gamma, 0.6, 1e-12);
}

// Cluster labels are arbitrary: section 5.4 constrains the partition, not the
// numbering, so relabelling a clustering must leave every pair-rank value
// unchanged. This is the fixture above with its clusters renumbered, which
// makes the member-list sizes ascend rather than descend. Every other pair-rank
// fixture numbers its clusters in nonincreasing size order, so a collection
// rule comparing the two sides' sizes would be satisfied everywhere in the
// suite and would still drop cross pairs here. The expected values are the
// ones the fixture above asserts, because it is the same partition.
TEST(ClusterReportTest, ReportPairRankIsUnchangedByClusterRelabelling) {
    DenseStorage storage(4);
    storage.Set(0, 1, 0.4);
    storage.Set(0, 2, 0.1);
    storage.Set(1, 2, 0.5);
    storage.Set(0, 3, 0.6);
    storage.Set(1, 3, 0.7);
    storage.Set(2, 3, 0.9);

    ClusterReportOptions options;
    options.compute_pair_rank_indices = true;
    // Clusters {2}, {3} and {0,1} in that order: sizes 1, 1, 2.
    const ClusterReport r =
        cluster_report(MakeResult({2, 2, 0, 1}), storage, options);

    // The partition is unchanged, so the within-distance is still 0.4 and the
    // five between-distances are still {0.1, 0.5, 0.6, 0.7, 0.9}. Pooled
    // ascending: 0.1 0.4 0.5 0.6 0.7 0.9, giving S_w = 0.4, S_min = 0.1 and
    // S_max = 0.9, so the index is (0.4 - 0.1) / (0.9 - 0.1) = 0.375.
    EXPECT_NEAR(r.c_index, 0.375, 1e-12);

    // Likewise s+ = 4 and s- = 1, giving (4 - 1) / 5 = 0.6.
    EXPECT_NEAR(r.baker_hubert_gamma, 0.6, 1e-12);
}

namespace {

// Six points in three clusters of two -- C0 = {0,1}, C1 = {2,3}, C2 = {4,5} --
// every intra pair at 0.1, and the three cross blocks flat but unequal:
// C0-C1 0.9, C0-C2 0.4, C1-C2 0.9.
//
// Three of the fused cross loop's accumulators are written through both loop
// indices, and a cluster is the inner index b only for partners of lower
// ordinal. C2 is therefore never the outer a, so its silhouette b term is
// written on the b-side alone -- and its two writes disagree: 0.4 from the pair
// (C0, C2), then 0.9 from (C1, C2). A b-side write that kept the last value
// instead of the smallest would report C2's b term as 0.9. Every other
// multi-cluster fixture in this file has cross blocks that are flat across all
// partners, where the two candidates coincide.
//
// The same asymmetry gives the Dunn numerator a unique minimum, 0.4, which is
// not the separation of the last cluster pair the loop visits.
DenseStorage MakeAsymmetricCrossBlockStorage() {
    DenseStorage storage(6);
    storage.Set(0, 1, 0.1);
    storage.Set(2, 3, 0.1);
    storage.Set(4, 5, 0.1);
    for (size_t i = 0; i < 2; ++i) {
        for (size_t j = 2; j < 4; ++j) {
            storage.Set(i, j, 0.9);
        }
        for (size_t j = 4; j < 6; ++j) {
            storage.Set(i, j, 0.4);
        }
    }
    for (size_t i = 2; i < 4; ++i) {
        for (size_t j = 4; j < 6; ++j) {
            storage.Set(i, j, 0.9);
        }
    }
    return storage;
}

// The same three-cluster shape with C2 equidistant from both partners at 0.5,
// so its nearest-cluster answer is a genuine tie and the documented
// lowest-ordinal rule has something to decide. C2 is again never the outer a,
// so the tie is resolved entirely by the b-side comparison.
DenseStorage MakeTiedThirdClusterStorage() {
    DenseStorage storage(6);
    storage.Set(0, 1, 0.1);
    storage.Set(2, 3, 0.1);
    storage.Set(4, 5, 0.1);
    for (size_t i = 0; i < 2; ++i) {
        for (size_t j = 2; j < 4; ++j) {
            storage.Set(i, j, 0.9);
        }
        for (size_t j = 4; j < 6; ++j) {
            storage.Set(i, j, 0.5);
        }
    }
    for (size_t i = 2; i < 4; ++i) {
        for (size_t j = 4; j < 6; ++j) {
            storage.Set(i, j, 0.5);
        }
    }
    return storage;
}

}  // namespace

TEST(ClusterReportTest, SilhouetteBTermKeepsTheNearestClusterOnTheBSide) {
    const DenseStorage storage = MakeAsymmetricCrossBlockStorage();
    ClusterReportOptions options;
    options.compute_per_cluster_records = true;
    const ClusterReport r =
        cluster_report(MakeResult({0, 0, 1, 1, 2, 2}), storage, options);

    ASSERT_EQ(r.records.size(), 3u);

    // Every a term is 0.1, the cluster's single intra distance. C0 is the outer
    // a in both its pairs and sees 0.9 then 0.4; C1 sees 0.9 on both sides.
    EXPECT_DOUBLE_EQ(r.records[0].silhouette, (0.4 - 0.1) / 0.4);
    EXPECT_DOUBLE_EQ(r.records[1].silhouette, (0.9 - 0.1) / 0.9);
    // The discriminating assertion: C2 sees 0.4 and then 0.9, both on the
    // b-side. Keeping the last write reports 0.888..., which is larger, still
    // inside the documented [-1, 1], and wrong.
    EXPECT_DOUBLE_EQ(r.records[2].silhouette, (0.4 - 0.1) / 0.4);

    // The scalar averages the same six terms, so it moves as well: 0.7963
    // against 0.8426.
    EXPECT_NEAR(r.silhouette,
                (4.0 * ((0.4 - 0.1) / 0.4) + 2.0 * ((0.9 - 0.1) / 0.9)) / 6.0,
                1e-12);
}

TEST(ClusterReportTest, DunnNumeratorIsTheGlobalMinimumSeparation) {
    const DenseStorage storage = MakeAsymmetricCrossBlockStorage();
    const ClusterReport r = cluster_report(
        MakeResult({0, 0, 1, 1, 2, 2}), storage, ClusterReportOptions());

    // The cross loop visits the pairs (C0,C1), (C0,C2) and (C1,C2), separating
    // at 0.9, 0.4 and 0.9, and every diameter is 0.1. Assigning each pair's
    // minimum rather than folding it leaves the last pair's 0.9 and reports
    // 9.0 -- a better-looking index than the true 4.0.
    EXPECT_NEAR(r.dunn_index, 0.4 / 0.1, 1e-12);
}

TEST(ClusterReportTest, RecordNearestClusterTieOnTheBSideKeepsTheLowestOrdinal) {
    const DenseStorage storage = MakeTiedThirdClusterStorage();
    ClusterReportOptions options;
    options.compute_per_cluster_records = true;
    const ClusterReport r =
        cluster_report(MakeResult({0, 0, 1, 1, 2, 2}), storage, options);

    ASSERT_EQ(r.records.size(), 3u);

    // C2's two candidates are both at 0.5 and are compared on the b-side, in
    // the order C0 then C1. ClusterReport.h documents the lowest ordinal as the
    // winner, so a <= there -- which the a-side tie test cannot see -- names
    // C1 instead.
    EXPECT_EQ(r.records[2].nearest_cluster, 0);
    EXPECT_DOUBLE_EQ(r.records[2].nearest_cluster_distance, 0.5);

    // The other two clusters are decided outright at 0.5 against 0.9, and are
    // asserted so a failure separates the tie from a wholesale change.
    EXPECT_EQ(r.records[0].nearest_cluster, 2);
    EXPECT_EQ(r.records[1].nearest_cluster, 2);
}

// The b-side counterpart. The tie test above pins the <= direction, where
// both of C2's candidates are equal; it cannot see a guard that accepts any
// different value. Here C2 is offered 0.4 from C0 and then 0.9 from C1, both
// on the b-side, so a != in place of the < names C1 at 0.9 -- the farther of
// the two clusters reported as the nearest.
// SilhouetteBTermKeepsTheNearestClusterOnTheBSide uses this same fixture but
// reads best_other_mean, a different accumulator, so it leaves these two
// fields unpinned.
TEST(ClusterReportTest, RecordNearestClusterOnTheBSideKeepsTheSmallestNotTheLast) {
    const DenseStorage storage = MakeAsymmetricCrossBlockStorage();
    ClusterReportOptions options;
    options.compute_per_cluster_records = true;
    const ClusterReport r =
        cluster_report(MakeResult({0, 0, 1, 1, 2, 2}), storage, options);

    ASSERT_EQ(r.records.size(), 3u);

    // The discriminating pair: C2's later candidate is the farther one.
    EXPECT_EQ(r.records[2].nearest_cluster, 0);
    EXPECT_DOUBLE_EQ(r.records[2].nearest_cluster_distance, 0.4);

    // C0 falls 0.9 then 0.4 on the a-side, so it is decided outright. C1 is
    // not: it ties at 0.9 between C0 on the b-side and C2 on the a-side, and
    // keeps the lowest ordinal. The line below therefore pins more than it
    // looks like -- a <= at the a-side guard names C2 there -- so it is not a
    // canary to be dropped when this test is next simplified.
    EXPECT_EQ(r.records[0].nearest_cluster, 2);
    EXPECT_DOUBLE_EQ(r.records[0].nearest_cluster_distance, 0.4);
    EXPECT_EQ(r.records[1].nearest_cluster, 0);
    EXPECT_DOUBLE_EQ(r.records[1].nearest_cluster_distance, 0.9);
}

TEST(ClusterReportTest, SilhouetteDenominatorIsTheLargerOfTheTwoTerms) {
    // Inverted separation: the within-distance 0.75 exceeds every cross
    // distance at 0.25, so a > b and the silhouette is negative. Every other
    // silhouette fixture in this file separates cleanly, and there a < b makes
    // the larger term and the b term the same number.
    DenseStorage storage(4);
    storage.Set(0, 1, 0.75);
    storage.Set(2, 3, 0.75);
    storage.Set(0, 2, 0.25);
    storage.Set(0, 3, 0.25);
    storage.Set(1, 2, 0.25);
    storage.Set(1, 3, 0.25);

    const ClusterReport r = cluster_report(
        MakeResult({0, 0, 1, 1}), storage, ClusterReportOptions());

    // Pinning the exact negative value pins the sign convention with it.
    // Dividing by the b term alone gives -2.0.
    EXPECT_DOUBLE_EQ(r.silhouette, (0.25 - 0.75) / 0.75);
    // The range is published, so it is asserted on its own lines: a term that
    // leaves [-1, 1] then fails as a range violation rather than only as a
    // changed number.
    EXPECT_GE(r.silhouette, -1.0);
    EXPECT_LE(r.silhouette, 1.0);
}

TEST(ClusterReportTest, RecordRadiusIsTheFarthestMemberNotTheLastOne) {
    // C0 = {0,1,2} elects sample 0: the member means are 0.25, 0.4 and 0.25,
    // and the earliest-member rule breaks the tie. Its farthest member is 1 at
    // 0.4, while the last member the fill walks is 2 at 0.1. The file's other
    // radius fixtures all happen to iterate the farthest member last, where
    // "largest so far" and "last seen" agree.
    DenseStorage storage(5);
    storage.Set(0, 1, 0.4);
    storage.Set(0, 2, 0.1);
    storage.Set(1, 2, 0.4);
    storage.Set(3, 4, 0.2);
    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 3; j < 5; ++j) {
            storage.Set(i, j, 0.9);
        }
    }

    ClusterReportOptions options;
    options.compute_per_cluster_records = true;
    const ClusterReport r =
        cluster_report(MakeResult({0, 0, 0, 1, 1}), storage, options);

    ASSERT_EQ(r.records.size(), 2u);
    EXPECT_EQ(r.records[0].representative, 0u);
    // Assigning the last distance instead of folding the maximum reports 0.1.
    EXPECT_DOUBLE_EQ(r.records[0].radius, 0.4);
    EXPECT_DOUBLE_EQ(r.records[1].radius, 0.2);
    // The scalar reads the same vector: the median of {0.4, 0.2} against the
    // median of {0.1, 0.2}.
    EXPECT_DOUBLE_EQ(r.median_radius, 0.3);
}

TEST(ClusterReportTest, GlobalMedoidCountsWithinClusterDistancesToo) {
    // Five samples, A = {0,1,2} and B = {3,4}. Sample 0 sits far from its own
    // mates at 0.9 and close to B at 0.2, which is what makes "least total
    // distance to all clustered points" and "least total distance to points in
    // other clusters" disagree. The file's three global-medoid tests pin the
    // tiebreak and the clustered-only filter but not the summand, and the one
    // with three clusters uses size-2 clusters throughout, where the
    // within-cluster term is a constant offset that cannot discriminate.
    DenseStorage storage(5);
    storage.Set(0, 1, 0.9);
    storage.Set(0, 2, 0.9);
    storage.Set(1, 2, 0.1);
    storage.Set(3, 4, 0.1);
    storage.Set(0, 3, 0.2);
    storage.Set(0, 4, 0.2);
    storage.Set(1, 3, 0.6);
    storage.Set(1, 4, 0.6);
    storage.Set(2, 3, 0.6);
    storage.Set(2, 4, 0.6);

    const ClusterReport r = cluster_report(
        MakeResult({0, 0, 0, 1, 1}), storage, ClusterReportOptions());

    // Totals over the other four clustered points are 2.2, 2.2, 2.2, 1.5, 1.5,
    // so M = 3. Counting only cross-cluster distances they are 0.4, 1.2, 1.2,
    // 1.4, 1.4, and M = 0. The cluster medoids are 1 and 3 either way, so the
    // within-scatter stays 0.9^2 + 0.1^2 + 0.1^2 = 0.83 and only the between
    // term moves: 3 * d(1,3)^2 = 1.08 at M = 3, against
    // 3 * d(1,0)^2 + 2 * d(3,0)^2 = 2.51 at M = 0, i.e. 3.904 against 9.072.
    EXPECT_NEAR(r.calinski_harabasz_medoid, 1.08 / (0.83 / 3.0), 1e-12);
}
