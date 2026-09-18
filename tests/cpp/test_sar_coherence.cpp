#include <gtest/gtest.h>

#include <cfloat>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "oecluster/StorageBackend.h"
#include "oecluster/clustering/ClusterTypes.h"
#include "oecluster/clustering/PartitionAgreement.h"
#include "oecluster/clustering/SARCoherence.h"
#include "../../src/clustering/ActivityMetrics.h"
#include "../../src/clustering/ContingencyTable.h"

namespace {

using OECluster::ActivityLandscape;
using OECluster::ActivityLandscapeOptions;
using OECluster::ClusterLabel;
using OECluster::Modelability;
using OECluster::ModelabilityOptions;
using OECluster::NoiseHandling;
using OECluster::SARCoherence;
using OECluster::SARCoherenceOptions;

constexpr double NOT_A_NUMBER = std::numeric_limits<double>::quiet_NaN();

/// Fills the upper triangle in condensed (i, j) order. The list must hold
/// exactly n * (n - 1) / 2 entries.
void FillStorage(OECluster::DenseStorage& storage,
                 const std::vector<double>& condensed) {
    const std::size_t n = storage.NumSamples();
    std::size_t k = 0;
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = i + 1; j < n; ++j) {
            storage.Set(i, j, condensed[k]);
            ++k;
        }
    }
}

OECluster::ClusteringResult MakeResult(std::vector<ClusterLabel> labels) {
    OECluster::Clusters members = OECluster::labels_to_clusters(labels);
    return OECluster::ClusteringResult(std::move(labels), std::move(members));
}

SARCoherence Coherence(const std::vector<ClusterLabel>& labels,
                       const std::vector<double>& activity,
                       NoiseHandling noise = NoiseHandling::Excluded) {
    SARCoherenceOptions options;
    options.noise_handling = noise;
    return OECluster::sar_coherence(labels, activity, options);
}

}  // namespace

TEST(ActivityMetricsTest, MarkMissingActivityLeavesOtherBitsAsFound) {
    const std::vector<double> activity = {1.0, NOT_A_NUMBER, 3.0};
    std::vector<bool> drop = {true, false, false};

    OECluster::detail::mark_missing_activity(activity, "caller", drop);

    EXPECT_TRUE(drop[0]);
    EXPECT_TRUE(drop[1]);
    EXPECT_FALSE(drop[2]);
}

TEST(ActivityMetricsTest, MarkMissingActivityRejectsAnInfiniteValue) {
    const std::vector<double> activity = {1.0,
                                          std::numeric_limits<double>::infinity()};
    std::vector<bool> drop(activity.size(), false);

    try {
        OECluster::detail::mark_missing_activity(activity, "sar_coherence", drop);
        FAIL() << "expected an infinite activity to be rejected";
    } catch (const std::invalid_argument& error) {
        const std::string message = error.what();
        EXPECT_NE(message.find("sar_coherence"), std::string::npos);
        EXPECT_NE(message.find("activity[1]"), std::string::npos);
        EXPECT_NE(message.find("infinite"), std::string::npos);
        EXPECT_NE(message.find("NaN"), std::string::npos);
    }
}

TEST(ActivityMetricsTest, GatherScoredReadsAFinalisedMask) {
    const std::vector<double> activity = {1.0, 2.0, 3.0, 4.0};
    const std::vector<bool> drop = {false, true, false, true};

    const OECluster::detail::ScoredActivity scored =
        OECluster::detail::gather_scored(activity, drop);

    EXPECT_EQ(scored.indices, (std::vector<std::size_t>{0, 2}));
    EXPECT_EQ(scored.values, (std::vector<double>{1.0, 3.0}));
    EXPECT_EQ(OECluster::detail::gather_indices(drop),
              (std::vector<std::size_t>{0, 2}));
}

// Both exclusion reasons must be ORed into one mask before either the values
// or the group ids are gathered. Gathering after only one of them produces a
// value vector and an id vector of different lengths, which is the failure
// this fixture pins.
TEST(ActivityMetricsTest, ExclusionAndMissingnessShareOneMask) {
    const std::vector<ClusterLabel> labels = {-1, 0, 0, -1, 1};
    const std::vector<double> activity = {1.0, 2.0, NOT_A_NUMBER, 5.0, 4.0};

    std::vector<bool> drop(activity.size(), false);
    OECluster::detail::mark_excluded(labels, NoiseHandling::Excluded, drop);
    OECluster::detail::mark_missing_activity(activity, "sar_coherence", drop);

    std::uint32_t num_ids = 0;
    const std::vector<std::uint32_t> ids = OECluster::detail::intern_side(
        labels, NoiseHandling::Excluded, drop, num_ids);
    const OECluster::detail::ScoredActivity scored =
        OECluster::detail::gather_scored(activity, drop);

    ASSERT_EQ(scored.values.size(), ids.size());
    ASSERT_EQ(scored.indices.size(), ids.size());
    EXPECT_EQ(scored.indices, (std::vector<std::size_t>{1, 4}));
    EXPECT_EQ(scored.values, (std::vector<double>{2.0, 4.0}));
    EXPECT_EQ(num_ids, 2u);

    // The alignment itself: position p's group id must be the id of the label
    // that position p was gathered from. Phrasing it as "equal ids iff equal
    // labels" asserts the mapping without reimplementing the interning order,
    // which is the property whose absence makes every group assignment wrong.
    for (std::size_t p = 0; p < ids.size(); ++p) {
        for (std::size_t q = 0; q < ids.size(); ++q) {
            EXPECT_EQ(ids[p] == ids[q],
                      labels[scored.indices[p]] == labels[scored.indices[q]])
                << "positions " << p << " and " << q;
        }
    }
}

TEST(ActivityMetricsTest, PopulationStddevMatchesAHandValue) {
    EXPECT_NEAR(OECluster::detail::population_stddev({1.0, 3.0, 5.0, 7.0}),
                2.23606797749979, 1e-12);
}

TEST(ActivityMetricsTest, PopulationStddevIsUndefinedBelowTwoValues) {
    EXPECT_TRUE(std::isnan(OECluster::detail::population_stddev({})));
    EXPECT_TRUE(std::isnan(OECluster::detail::population_stddev({4.0})));
}

// Two DBL_MAX values have an exact population stddev of zero, but the mean
// cannot be formed in double. The helper reports the infinity rather than a
// plausible zero; the caller turns it into a throw.
TEST(ActivityMetricsTest, PopulationStddevReportsAnUnformableMean) {
    const double spread = OECluster::detail::population_stddev({DBL_MAX, DBL_MAX});
    EXPECT_TRUE(std::isinf(spread));
    EXPECT_GT(spread, 0.0);
}

// The scaled accumulation could return a finite 1e300 here -- the spread is
// representable even though the sum of squared deviations is not -- and that is
// exactly what must not happen: a sum of squared deviations that does not fit
// in a double is outside the supported domain, sums_of_squares throws on it,
// and the two helpers have to agree about where the domain ends. This is the
// unit-level half of ActivityLandscapeTest.RejectsAnOverflowingSpread;
// without it, a scaling change can silently widen the domain and only a much
// later integration test would notice.
TEST(ActivityMetricsTest, PopulationStddevReportsAnUnrepresentableSpread) {
    const double spread = OECluster::detail::population_stddev({-1e300, 1e300});
    EXPECT_TRUE(std::isinf(spread));
    EXPECT_GT(spread, 0.0);

    // The companion refusal, so the disagreement would fail here too.
    EXPECT_THROW(OECluster::detail::sums_of_squares({0, 1}, {-1e300, 1e300}, 2,
                                                    "sar_coherence"),
                 std::invalid_argument);
}

TEST(ActivityMetricsTest, SumsOfSquaresDecomposesAHandFixture) {
    const std::vector<std::uint32_t> ids = {0, 0, 1, 1, 2, 2};
    const std::vector<double> values = {1.0, 3.0, 5.0, 7.0, 9.0, 11.0};

    const OECluster::detail::SumsOfSquares ss =
        OECluster::detail::sums_of_squares(ids, values, 3, "sar_coherence");

    // The struct holds the sums in units of scale squared, so the hand values
    // are read back through the accessors. The scale is 5: the grand mean is
    // 6 and the widest deviation is the one at 1 and at 11.
    EXPECT_EQ(ss.scale, 5.0);
    EXPECT_NEAR(OECluster::detail::ss_total(ss), 70.0, 1e-12);
    EXPECT_NEAR(OECluster::detail::ss_within(ss), 6.0, 1e-12);
    EXPECT_NEAR(OECluster::detail::ss_between(ss), 64.0, 1e-12);
    EXPECT_NEAR(ss.total, 70.0 / 25.0, 1e-15);
    EXPECT_EQ(ss.num_scored, 6u);
    EXPECT_EQ(ss.num_groups, 3u);

    // The two general invariants, at the two different strengths
    // sums_of_squares documents. The bounds are exact because between is
    // clamped into them; the additivity is not, because between is derived by
    // a subtraction that rounds, so asserting bitwise equality here would be
    // asserting something the derivation does not deliver. Checking these on
    // the hand fixture is what makes them invariants rather than three lucky
    // numbers: the expected values above would still pass if the clamp were
    // dropped.
    EXPECT_GE(ss.between, 0.0);
    EXPECT_LE(ss.between, ss.total);
    EXPECT_NEAR(ss.between + ss.within, ss.total,
                4.0 * std::numeric_limits<double>::epsilon() * ss.total);
}

// SS_between is derived, not computed, so a perfect separation must land on
// the endpoints exactly rather than within a tolerance.
TEST(ActivityMetricsTest, SumsOfSquaresIsExactOnAPerfectSeparation) {
    const std::vector<std::uint32_t> ids = {0, 0, 1, 1};
    const std::vector<double> values = {1.0, 1.0, 5.0, 5.0};

    const OECluster::detail::SumsOfSquares ss =
        OECluster::detail::sums_of_squares(ids, values, 2, "sar_coherence");

    EXPECT_EQ(ss.within, 0.0);
    EXPECT_EQ(ss.between, ss.total);
    EXPECT_EQ(OECluster::detail::eta_squared(ss), 1.0);
}

TEST(ActivityMetricsTest, SumsOfSquaresRejectsMismatchedLengths) {
    const std::vector<std::uint32_t> ids = {0, 0, 1};
    const std::vector<double> values = {1.0, 2.0};

    EXPECT_THROW(OECluster::detail::sums_of_squares(ids, values, 2, "caller"),
                 std::invalid_argument);
}

TEST(ActivityMetricsTest, EffectSizesMatchAHandFixture) {
    const std::vector<std::uint32_t> ids = {0, 0, 1, 1, 2, 2};
    const std::vector<double> values = {1.0, 3.0, 5.0, 7.0, 9.0, 11.0};
    const OECluster::detail::SumsOfSquares ss =
        OECluster::detail::sums_of_squares(ids, values, 3, "sar_coherence");

    EXPECT_NEAR(OECluster::detail::eta_squared(ss), 0.9142857142857143, 1e-12);
    EXPECT_NEAR(OECluster::detail::omega_squared(ss), 0.8333333333333334, 1e-12);
}

// Chance correction may go below zero: a labelling that separates activity
// worse than a random one of the same shape earns a negative omega squared.
TEST(ActivityMetricsTest, OmegaSquaredGoesNegativeBelowChance) {
    const std::vector<std::uint32_t> ids = {0, 0, 1, 1};
    const std::vector<double> values = {1.0, 5.0, 1.0, 5.0};
    const OECluster::detail::SumsOfSquares ss =
        OECluster::detail::sums_of_squares(ids, values, 2, "sar_coherence");

    EXPECT_EQ(OECluster::detail::eta_squared(ss), 0.0);
    EXPECT_NEAR(OECluster::detail::omega_squared(ss), -0.3333333333333333,
                1e-12);
}

// One group explains nothing, however wide its spread. That is 0.0, not the
// NaN that a naive zero-check on SS_between would produce.
TEST(ActivityMetricsTest, ASingleGroupExplainsNothing) {
    const std::vector<std::uint32_t> ids = {0, 0, 0};
    const std::vector<double> values = {1.0, 2.0, 3.0};
    const OECluster::detail::SumsOfSquares ss =
        OECluster::detail::sums_of_squares(ids, values, 1, "sar_coherence");

    EXPECT_GT(ss.total, 0.0);
    EXPECT_EQ(OECluster::detail::eta_squared(ss), 0.0);
    EXPECT_EQ(OECluster::detail::omega_squared(ss), 0.0);
}

TEST(ActivityMetricsTest, SumsOfSquaresRejectsAnUnformableMean) {
    const std::vector<std::uint32_t> ids = {0, 0};
    const std::vector<double> values = {DBL_MAX, DBL_MAX};

    try {
        OECluster::detail::sums_of_squares(ids, values, 1, "sar_coherence");
        FAIL() << "expected an unformable mean to be rejected";
    } catch (const std::invalid_argument& error) {
        const std::string message = error.what();
        EXPECT_NE(message.find("mean"), std::string::npos);
        EXPECT_EQ(message.find("SS_total"), std::string::npos);
    }
}

TEST(ActivityMetricsTest, SumsOfSquaresRejectsAnUnformableTotal) {
    const double big = std::sqrt(DBL_MAX);
    const std::vector<std::uint32_t> ids = {0, 0, 1, 1};
    const std::vector<double> values = {-big, big, -big, big};

    try {
        OECluster::detail::sums_of_squares(ids, values, 2, "sar_coherence");
        FAIL() << "expected an unformable SS_total to be rejected";
    } catch (const std::invalid_argument& error) {
        const std::string message = error.what();
        EXPECT_NE(message.find("SS_total"), std::string::npos);
        EXPECT_EQ(message.find("mean"), std::string::npos);
    }
}

// The specified evaluation order divides by SS_total first, which cancels the
// input's scale. These two fixtures sit at the extreme ends of the supported
// domain and both have an exact omega squared of -0.2; the two rejected
// orderings lose one end to overflow and the other to underflow.
TEST(ActivityMetricsTest, OmegaSquaredSurvivesTheOverflowEnd) {
    const std::vector<std::uint32_t> ids = {0, 0, 1};
    const std::vector<double> values = {-std::sqrt(5e307), std::sqrt(5e307),
                                        std::sqrt(7.5e307)};
    const OECluster::detail::SumsOfSquares ss =
        OECluster::detail::sums_of_squares(ids, values, 2, "sar_coherence");

    EXPECT_NEAR(OECluster::detail::omega_squared(ss), -0.2, 1e-12);
}

TEST(ActivityMetricsTest, OmegaSquaredSurvivesTheUnderflowEnd) {
    const double tiny = std::ldexp(1.0, -537);
    const std::vector<std::uint32_t> ids = {0, 0, 0, 1, 1, 1};
    const std::vector<double> values = {-tiny, 0.0, tiny, 0.0, 0.0, 0.0};
    const OECluster::detail::SumsOfSquares ss =
        OECluster::detail::sums_of_squares(ids, values, 2, "sar_coherence");

    EXPECT_NEAR(OECluster::detail::omega_squared(ss), -0.2, 1e-12);
}

// The two above pin the evaluation order of omega squared. These two pin the
// scaling of the sums it reads, which is a separate failure and a worse one:
// at 2^-539 a squared deviation is below the smallest subnormal and flushes to
// zero, and it does so for SS_within before SS_total, so the unscaled
// arithmetic reports eta squared as exactly 1.0 -- perfect separation -- on a
// fixture whose true value is 18/31. Both fixtures are inside the supported
// domain: every value, difference and sum here is a finite double. The expected
// numbers are exact rationals, so a tolerance would hide a partial fix.
TEST(ActivityMetricsTest, EffectSizesSurviveSquaredDeviationUnderflow) {
    const double x = std::ldexp(1.0, -539);
    const std::vector<std::uint32_t> ids = {0, 0, 0, 1, 1, 1};
    const std::vector<double> values = {-3.0 * x, -2.0 * x, x,
                                        x, 2.0 * x, 5.0 * x};
    const OECluster::detail::SumsOfSquares ss =
        OECluster::detail::sums_of_squares(ids, values, 2, "sar_coherence");

    // Scaled, the sums are ordinary numbers: the largest term is 1.0 and the
    // total cannot leave [1, num_scored].
    EXPECT_GT(ss.total, 1.0);
    EXPECT_GT(ss.within, 0.0);
    EXPECT_NEAR(OECluster::detail::eta_squared(ss), 18.0 / 31.0, 1e-15);
    EXPECT_NEAR(OECluster::detail::omega_squared(ss), 59.0 / 137.0, 1e-15);
}

TEST(ActivityMetricsTest, PopulationStddevSurvivesSquaredDeviationUnderflow) {
    const double x = std::ldexp(1.0, -539);
    const double spread = OECluster::detail::population_stddev(
        {-3.0 * x, -2.0 * x, x, x, 2.0 * x, 5.0 * x});

    // Unscaled this is 0.0, which would put every sample inside the RMODI
    // band and report a landscape with no structure to contradict.
    EXPECT_GT(spread, 0.0);
    EXPECT_NEAR(spread, 1.4585016579561791e-162, 1e-177);
}

// A naive mean accumulation's error scales with magnitude, while deviations
// scale with spread. When the spread is near the values' rounding floor the
// naive mean can land a full ulp away, reversing the effect sizes. This fixture
// has values {1.5 + 9u, 1.5 + 10u, 1.5 + 10u, 1.5 + 11u} where u = 2^-52; a
// naive grand mean lands at 1.5 + 11u instead of 1.5 + 10u, driving eta_squared
// to exactly 0.0 and omega_squared below zero. The correction pass recovers
// positive effect sizes.
TEST(ActivityMetricsTest, MeanCorrectionRescuesSpreadAtTheRoundingFloor) {
    const std::vector<std::uint32_t> ids = {0, 0, 0, 1};
    const std::vector<double> values = {0x1.8000000000009p+0, 0x1.800000000000ap+0,
                                        0x1.800000000000ap+0, 0x1.800000000000bp+0};
    const OECluster::detail::SumsOfSquares ss =
        OECluster::detail::sums_of_squares(ids, values, 2, "test");

    const double eta2 = OECluster::detail::eta_squared(ss);
    const double omega2 = OECluster::detail::omega_squared(ss);
    const double stddev = OECluster::detail::population_stddev(values);

    // The exact answers are eta_squared = 2/3, omega_squared = 3/7, and
    // stddev = 1.5700924586837752e-16. The true group mean 1.5 + (29/3)u is not
    // representable, so the correction cannot recover the exact 2/3; it does
    // recover the exact stddev.
    EXPECT_DOUBLE_EQ(eta2, 0.5);
    EXPECT_DOUBLE_EQ(omega2, 0.2);
    EXPECT_DOUBLE_EQ(stddev, 1.5700924586837752e-16);

    // Before the fix these were exactly 0.0, -1/3, and slightly high.
    EXPECT_GT(eta2, 0.0);
    EXPECT_GT(omega2, 0.0);
}

// The ss_between clamp prevents negative effect sizes when within exceeds
// total by rounding. Mathematically within <= total because each value is at
// least as close to its group mean as to the grand mean, but the two sums
// accumulate in separate loops and round independently. On tightly spaced input
// where the group means straddle the grand mean, within can finish one ulp above
// total. This fixture — {1.5, 1.5, 1.5 + 10u, 1.5 + 11u} with u = 2^-52,
// interleaved across two groups — produces total = 3.0833333333333335 and
// within = 3.0833333333333339, so total - within is strictly negative. Without
// the clamp, eta_squared would report a negative value.
TEST(ActivityMetricsTest, BetweenClampsToZeroWhenWithinExceedsTotal) {
    const std::vector<std::uint32_t> ids = {0, 1, 0, 1};
    const std::vector<double> values = {0x1.8p+0, 0x1.8p+0,
                                        0x1.800000000000ap+0, 0x1.800000000000bp+0};
    const OECluster::detail::SumsOfSquares ss =
        OECluster::detail::sums_of_squares(ids, values, 2, "test");

    // The clamp drives between and eta_squared to exactly zero where the
    // unclamped total - within is strictly negative.
    EXPECT_DOUBLE_EQ(ss.between, 0.0);
    EXPECT_DOUBLE_EQ(OECluster::detail::eta_squared(ss), 0.0);
}

TEST(SARCoherenceTest, DecomposesAThreeClusterFixture) {
    const SARCoherence coherence =
        Coherence({0, 0, 1, 1, 2, 2}, {1.0, 3.0, 5.0, 7.0, 9.0, 11.0});

    EXPECT_EQ(coherence.num_samples, 6u);
    EXPECT_EQ(coherence.num_scored, 6u);
    EXPECT_EQ(coherence.num_clusters, 3u);
    EXPECT_NEAR(coherence.eta_squared, 0.9142857142857143, 1e-12);
    EXPECT_NEAR(coherence.omega_squared, 0.8333333333333334, 1e-12);

    ASSERT_EQ(coherence.clusters.size(), 3u);
    EXPECT_EQ(coherence.clusters[0].label, 0);
    EXPECT_EQ(coherence.clusters[0].num_scored, 2u);
    EXPECT_NEAR(coherence.clusters[0].mean_activity, 2.0, 1e-12);
    EXPECT_NEAR(coherence.clusters[0].stddev_activity, 1.0, 1e-12);
    EXPECT_NEAR(coherence.clusters[1].mean_activity, 6.0, 1e-12);
    EXPECT_NEAR(coherence.clusters[2].mean_activity, 10.0, 1e-12);
}

TEST(SARCoherenceTest, ReportsAPerfectSeparationExactly) {
    const SARCoherence coherence = Coherence({0, 0, 1, 1}, {1.0, 1.0, 5.0, 5.0});

    EXPECT_EQ(coherence.eta_squared, 1.0);
    EXPECT_EQ(coherence.omega_squared, 1.0);
}

TEST(SARCoherenceTest, ReportsNoSeparationExactly) {
    const SARCoherence coherence = Coherence({0, 0, 1, 1}, {1.0, 5.0, 1.0, 5.0});

    EXPECT_EQ(coherence.eta_squared, 0.0);
    EXPECT_NEAR(coherence.omega_squared, -0.3333333333333333, 1e-12);
}

// Promoting noise to singletons gives every noise point a group mean equal to
// its own value, which reads as perfectly explained variance. That is why the
// default is Excluded, and this fixture is the evidence.
TEST(SARCoherenceTest, NoiseHandlingChangesTheEffectSize) {
    const std::vector<ClusterLabel> labels = {-1, -2, 0, 0, 1, 1};
    const std::vector<double> activity = {10.0, 0.0, 2.0, 4.0, 6.0, 8.0};

    const SARCoherence excluded =
        Coherence(labels, activity, NoiseHandling::Excluded);
    EXPECT_EQ(excluded.num_samples, 6u);
    EXPECT_EQ(excluded.num_scored, 4u);
    EXPECT_EQ(excluded.num_clusters, 2u);
    EXPECT_NEAR(excluded.eta_squared, 0.8, 1e-12);
    EXPECT_NEAR(excluded.omega_squared, 0.6363636363636364, 1e-12);

    const SARCoherence grouped =
        Coherence(labels, activity, NoiseHandling::Grouped);
    EXPECT_EQ(grouped.num_samples, 6u);
    EXPECT_EQ(grouped.num_scored, 6u);
    EXPECT_EQ(grouped.num_clusters, 3u);
    EXPECT_NEAR(grouped.eta_squared, 0.22857142857142856, 1e-12);
    EXPECT_NEAR(grouped.omega_squared, -0.22727272727272727, 1e-12);
    EXPECT_EQ(grouped.clusters[0].label, -1);
    EXPECT_EQ(grouped.clusters[0].num_scored, 2u);

    const SARCoherence singletons =
        Coherence(labels, activity, NoiseHandling::Singletons);
    EXPECT_EQ(singletons.num_samples, 6u);
    EXPECT_EQ(singletons.num_clusters, 4u);
    EXPECT_NEAR(singletons.eta_squared, 0.9428571428571428, 1e-12);
    EXPECT_NEAR(singletons.omega_squared, 0.8333333333333334, 1e-12);
    EXPECT_GT(singletons.eta_squared, excluded.eta_squared);
    EXPECT_EQ(singletons.clusters[0].label, -1);
    EXPECT_EQ(singletons.clusters[1].label, -2);
}

// num_samples is the length of the input whatever noise handling drops, so a
// caller can always tell how many measurements were missing.
TEST(SARCoherenceTest, MissingActivitiesAreExcludedWithoutShrinkingNumSamples) {
    const SARCoherence coherence =
        Coherence({0, 0, 1, 1, 1}, {1.0, 3.0, NOT_A_NUMBER, 5.0, 7.0});

    EXPECT_EQ(coherence.num_samples, 5u);
    EXPECT_EQ(coherence.num_scored, 4u);
    EXPECT_EQ(coherence.num_clusters, 2u);
    EXPECT_NEAR(coherence.eta_squared, 0.8, 1e-12);
    EXPECT_EQ(coherence.clusters[1].num_scored, 2u);

    // Dropping a NaN must give the same answer as never having supplied it:
    // exactly, not to a tolerance, since the surviving values reach the
    // decomposition in the same order either way. Only num_samples differs,
    // which is the whole point of reporting it separately.
    const SARCoherence prefiltered = Coherence({0, 0, 1, 1}, {1.0, 3.0, 5.0, 7.0});
    EXPECT_EQ(coherence.eta_squared, prefiltered.eta_squared);
    EXPECT_EQ(coherence.omega_squared, prefiltered.omega_squared);
    EXPECT_EQ(coherence.num_scored, prefiltered.num_scored);
    EXPECT_EQ(prefiltered.num_samples, 4u);
}

// Row order is a documented part of the contract, and the interning runs over
// the finalized drop mask rather than the raw labels -- so a cluster's position
// comes from its earliest *scored* member, not its earliest member. Here label
// 1 appears first in the input but its first sample has no measurement, so
// label 0 takes the leading row. Writing the contract the other way round is
// the easy mistake; this fixture is the one that tells the two apart.
TEST(SARCoherenceTest, RowOrderFollowsTheFirstScoredOccurrence) {
    const SARCoherence coherence = Coherence({1, 0, 1}, {NOT_A_NUMBER, 2.0, 3.0});

    ASSERT_EQ(coherence.clusters.size(), 2u);
    EXPECT_EQ(coherence.clusters[0].label, 0);
    EXPECT_EQ(coherence.clusters[1].label, 1);
    EXPECT_EQ(coherence.clusters[0].num_scored, 1u);
    EXPECT_EQ(coherence.clusters[1].num_scored, 1u);
}

TEST(SARCoherenceTest, BothOverloadsAgree) {
    const std::vector<ClusterLabel> labels = {-1, 0, 0, 1, 1};
    const std::vector<double> activity = {1.0, 2.0, NOT_A_NUMBER, 3.0, 4.0};

    const SARCoherence from_labels = OECluster::sar_coherence(labels, activity);
    const SARCoherence from_result =
        OECluster::sar_coherence(MakeResult(labels), activity);

    EXPECT_EQ(from_result.num_samples, from_labels.num_samples);
    EXPECT_EQ(from_result.num_scored, from_labels.num_scored);
    EXPECT_EQ(from_result.num_clusters, from_labels.num_clusters);
    EXPECT_EQ(from_result.eta_squared, from_labels.eta_squared);
    EXPECT_EQ(from_result.omega_squared, from_labels.omega_squared);
    ASSERT_EQ(from_result.clusters.size(), from_labels.clusters.size());
    for (std::size_t i = 0; i < from_result.clusters.size(); ++i) {
        EXPECT_EQ(from_result.clusters[i].label, from_labels.clusters[i].label);
        EXPECT_EQ(from_result.clusters[i].num_scored,
                  from_labels.clusters[i].num_scored);
        EXPECT_EQ(from_result.clusters[i].mean_activity,
                  from_labels.clusters[i].mean_activity);
        const double result_stddev = from_result.clusters[i].stddev_activity;
        const double labels_stddev = from_labels.clusters[i].stddev_activity;
        EXPECT_EQ(std::isnan(result_stddev), std::isnan(labels_stddev));
        if (!std::isnan(result_stddev)) {
            EXPECT_EQ(result_stddev, labels_stddev);
        }
    }
}

TEST(SARCoherenceTest, UndefinedValuesFollowTheDocumentedTable) {
    const SARCoherence one_scored = Coherence({0, 0}, {1.0, NOT_A_NUMBER});
    EXPECT_EQ(one_scored.num_scored, 1u);
    EXPECT_TRUE(std::isnan(one_scored.eta_squared));
    EXPECT_TRUE(std::isnan(one_scored.omega_squared));
    ASSERT_EQ(one_scored.clusters.size(), 1u);
    EXPECT_NEAR(one_scored.clusters[0].mean_activity, 1.0, 1e-12);
    EXPECT_TRUE(std::isnan(one_scored.clusters[0].stddev_activity));

    const SARCoherence no_variance = Coherence({0, 0, 1, 1}, {3.0, 3.0, 3.0, 3.0});
    EXPECT_TRUE(std::isnan(no_variance.eta_squared));
    EXPECT_TRUE(std::isnan(no_variance.omega_squared));
    EXPECT_NEAR(no_variance.clusters[0].mean_activity, 3.0, 1e-12);

    const SARCoherence all_singletons = Coherence({0, 1, 2}, {1.0, 2.0, 3.0});
    EXPECT_EQ(all_singletons.eta_squared, 1.0);
    EXPECT_TRUE(std::isnan(all_singletons.omega_squared));

    const SARCoherence all_missing =
        Coherence({0, 0}, {NOT_A_NUMBER, NOT_A_NUMBER});
    EXPECT_EQ(all_missing.num_samples, 2u);
    EXPECT_EQ(all_missing.num_scored, 0u);
    EXPECT_EQ(all_missing.num_clusters, 0u);
    EXPECT_TRUE(all_missing.clusters.empty());
    EXPECT_TRUE(std::isnan(all_missing.eta_squared));
    EXPECT_TRUE(std::isnan(all_missing.omega_squared));

    // One cluster with a real spread explains none of it. That is zero, not
    // undefined: num_clusters < 2 is deliberately absent from the NaN table.
    const SARCoherence one_cluster = Coherence({0, 0, 0}, {1.0, 2.0, 3.0});
    EXPECT_EQ(one_cluster.num_clusters, 1u);
    EXPECT_EQ(one_cluster.eta_squared, 0.0);
    EXPECT_EQ(one_cluster.omega_squared, 0.0);
}

TEST(SARCoherenceTest, ReportsDistinctStandardDeviationsPerCluster) {
    const SARCoherence coherence = Coherence({0, 0, 1, 1}, {3.0, 3.0, 1.0, 5.0});

    ASSERT_EQ(coherence.clusters.size(), 2u);
    EXPECT_NEAR(coherence.clusters[0].mean_activity, 3.0, 1e-12);
    EXPECT_EQ(coherence.clusters[0].stddev_activity, 0.0);
    EXPECT_NEAR(coherence.clusters[1].mean_activity, 3.0, 1e-12);
    EXPECT_NEAR(coherence.clusters[1].stddev_activity, 2.0, 1e-12);
}

TEST(SARCoherenceTest, RowMeansMatchTheDecompositionAtTheRoundingFloor) {
    const SARCoherence coherence = Coherence(
        {0, 0, 0, 1},
        {0x1.8000000000009p+0, 0x1.800000000000ap+0, 0x1.800000000000ap+0,
         0x1.800000000000bp+0});

    ASSERT_EQ(coherence.clusters.size(), 2u);
    EXPECT_EQ(coherence.clusters[0].num_scored, 3u);
    EXPECT_EQ(coherence.clusters[0].mean_activity, 0x1.800000000000ap+0);
    EXPECT_EQ(coherence.clusters[1].num_scored, 1u);
    EXPECT_EQ(coherence.clusters[1].mean_activity, 0x1.800000000000bp+0);
}

// Two tight groups a hair apart on a large offset. Catastrophic cancellation
// in a one-pass decomposition reports 1.0 or a negative total here; the
// two-pass form reports a plain interior value.
TEST(SARCoherenceTest, SurvivesALargeOffsetWithATinySpread) {
    const SARCoherence coherence = Coherence(
        {0, 0, 1, 1}, {1e8, 1e8 + 1e-4, 1e8 + 1e-4, 1e8 + 2e-4});

    EXPECT_GT(coherence.eta_squared, 0.0);
    EXPECT_LT(coherence.eta_squared, 1.0);
    EXPECT_NEAR(coherence.eta_squared, 0.5, 1e-6);
}

TEST(SARCoherenceTest, RejectsActivitiesOutsideTheSupportedRange) {
    const double big = std::sqrt(DBL_MAX);
    try {
        Coherence({0, 0, 1, 1}, {-big, big, -big, big});
        FAIL() << "expected an unformable SS_total to be rejected";
    } catch (const std::invalid_argument& error) {
        EXPECT_NE(std::string(error.what()).find("SS_total"), std::string::npos);
    }

    try {
        Coherence({0, 0}, {DBL_MAX, DBL_MAX});
        FAIL() << "expected an unformable mean to be rejected";
    } catch (const std::invalid_argument& error) {
        EXPECT_NE(std::string(error.what()).find("mean"), std::string::npos);
    }
}

TEST(SARCoherenceTest, RejectsAnInfiniteActivity) {
    try {
        Coherence({0, 0}, {1.0, std::numeric_limits<double>::infinity()});
        FAIL() << "expected an infinite activity to be rejected";
    } catch (const std::invalid_argument& error) {
        const std::string message = error.what();
        EXPECT_NE(message.find("sar_coherence"), std::string::npos);
        EXPECT_NE(message.find("activity[1]"), std::string::npos);
    }
}

TEST(SARCoherenceTest, RejectsAnEmptyOrMismatchedActivity) {
    try {
        Coherence({0, 1}, {});
        FAIL() << "expected an empty activity to be rejected";
    } catch (const std::invalid_argument& error) {
        EXPECT_NE(std::string(error.what()).find("non-empty"),
                  std::string::npos);
    }

    try {
        Coherence({0, 1, 2}, {1.0, 2.0});
        FAIL() << "expected a length mismatch to be rejected";
    } catch (const std::invalid_argument& error) {
        const std::string message = error.what();
        EXPECT_NE(message.find("activity has 2 entries"), std::string::npos);
        EXPECT_NE(message.find("the clustering has 3 samples"),
                  std::string::npos);
    }
}

TEST(ActivityLandscapeTest, CountsCliffsOnBothBoundariesInclusively) {
    OECluster::DenseStorage storage(4);
    FillStorage(storage, {0.30, 0.31, 0.20, 0.50, 0.10, 0.40});

    const ActivityLandscape landscape =
        OECluster::activity_landscape(storage, {0.0, 1.0, 5.0, 0.5});

    EXPECT_EQ(landscape.num_samples, 4u);
    EXPECT_EQ(landscape.num_scored, 4u);
    EXPECT_EQ(landscape.num_pairs_scored, 6u);
    EXPECT_EQ(landscape.num_cliffs, 1u);
    EXPECT_NEAR(landscape.cliff_density, 1.0 / 6.0, 1e-12);
}

// The boundary fixture above pins the default thresholds. This one pins that
// they are read from the options at all: hard-coding 0.30 and 1.0 in the sweep
// would leave validate_landscape_options intact and every other test passing.
TEST(ActivityLandscapeTest, CliffThresholdsChangeWhichPairsCount) {
    OECluster::DenseStorage storage(4);
    FillStorage(storage, {0.30, 0.31, 0.20, 0.50, 0.10, 0.40});
    const std::vector<double> activity = {0.0, 1.0, 5.0, 0.5};

    ActivityLandscapeOptions loose_distance;
    loose_distance.distance_threshold = 0.45;
    const ActivityLandscape wider =
        OECluster::activity_landscape(storage, activity, loose_distance);
    EXPECT_EQ(wider.num_cliffs, 3u);
    EXPECT_DOUBLE_EQ(wider.cliff_density, 0.5);

    ActivityLandscapeOptions strict_activity;
    strict_activity.activity_threshold = 4.5;
    const ActivityLandscape steeper =
        OECluster::activity_landscape(storage, activity, strict_activity);
    EXPECT_EQ(steeper.num_cliffs, 0u);
    EXPECT_DOUBLE_EQ(steeper.cliff_density, 0.0);
}

TEST(ActivityLandscapeTest, ScoresSaliAgainstAHandFixture) {
    OECluster::DenseStorage storage(3);
    FillStorage(storage, {0.5, 0.25, 0.125});

    const ActivityLandscape landscape =
        OECluster::activity_landscape(storage, {0.0, 1.0, 3.0});

    EXPECT_EQ(landscape.num_zero_distance_pairs, 0u);
    EXPECT_DOUBLE_EQ(landscape.max_sali, 16.0);
    EXPECT_DOUBLE_EQ(landscape.mean_sali, 10.0);
    EXPECT_EQ(landscape.num_cliffs, 2u);
    EXPECT_NEAR(landscape.cliff_density, 2.0 / 3.0, 1e-12);
    EXPECT_NEAR(landscape.activity_stddev, 1.2472191289246473, 1e-12);
    EXPECT_EQ(landscape.rmodi, 0.0);
}

// The mirror of the hand fixture above, with the activity descending instead of
// ascending. Every other value-asserting fixture is nondecreasing in scored
// order, where fabs(a_p - a_q) and a_q - a_p agree on every pair; this one is
// the only place the absolute value is load-bearing.
TEST(ActivityLandscapeTest, ScoresSaliOnADescendingActivity) {
    OECluster::DenseStorage storage(3);
    FillStorage(storage, {0.5, 0.25, 0.125});

    const ActivityLandscape landscape =
        OECluster::activity_landscape(storage, {3.0, 1.0, 0.0});

    EXPECT_DOUBLE_EQ(landscape.max_sali, 12.0);
    EXPECT_DOUBLE_EQ(landscape.mean_sali, 8.0);
    EXPECT_EQ(landscape.num_cliffs, 2u);
    EXPECT_NEAR(landscape.cliff_density, 2.0 / 3.0, 1e-12);
    EXPECT_NEAR(landscape.activity_stddev, 1.2472191289246473, 1e-12);
    EXPECT_EQ(landscape.rmodi, 0.0);
}

TEST(ActivityLandscapeTest, ReportsZeroDistancePairsRatherThanScoringThem) {
    OECluster::DenseStorage storage(3);
    FillStorage(storage, {0.0, 0.5, 0.25});

    const ActivityLandscape landscape =
        OECluster::activity_landscape(storage, {0.0, 1.0, 2.0});

    EXPECT_EQ(landscape.num_zero_distance_pairs, 1u);
    EXPECT_DOUBLE_EQ(landscape.max_sali, 4.0);
    EXPECT_DOUBLE_EQ(landscape.mean_sali, 4.0);
}

TEST(ActivityLandscapeTest, LeavesSaliUndefinedWhenEveryPairIsCoincident) {
    OECluster::DenseStorage storage(3);
    FillStorage(storage, {0.0, 0.0, 0.0});

    const ActivityLandscape landscape =
        OECluster::activity_landscape(storage, {0.0, 1.0, 2.0});

    EXPECT_EQ(landscape.num_zero_distance_pairs, 3u);
    EXPECT_TRUE(std::isnan(landscape.max_sali));
    EXPECT_TRUE(std::isnan(landscape.mean_sali));
    EXPECT_EQ(landscape.num_cliffs, 3u);
    EXPECT_DOUBLE_EQ(landscape.cliff_density, 1.0);
    EXPECT_EQ(landscape.rmodi, 0.0);
}

// A zero spread makes the band zero, so every neighbour is in the same band
// and every molecule is concordant. Nothing here is undefined.
TEST(ActivityLandscapeTest, HandlesAFlatActivity) {
    OECluster::DenseStorage storage(3);
    FillStorage(storage, {0.5, 0.25, 0.125});

    const ActivityLandscape landscape =
        OECluster::activity_landscape(storage, {5.0, 5.0, 5.0});

    EXPECT_EQ(landscape.activity_stddev, 0.0);
    EXPECT_EQ(landscape.rmodi, 1.0);
    EXPECT_EQ(landscape.num_cliffs, 0u);
    EXPECT_EQ(landscape.max_sali, 0.0);
    EXPECT_EQ(landscape.mean_sali, 0.0);
}

// A molecule with no same-band neighbour keeps same_min at infinity and
// cannot be concordant.
TEST(ActivityLandscapeTest, AMoleculeWithNoSameBandNeighbourIsDiscordant) {
    OECluster::DenseStorage storage(3);
    FillStorage(storage, {0.5, 0.25, 0.125});

    const ActivityLandscape landscape =
        OECluster::activity_landscape(storage, {0.0, 10.0, 10.5});

    EXPECT_NEAR(landscape.rmodi, 2.0 / 3.0, 1e-12);
}

// A molecule with no different-band neighbour keeps diff_min at infinity and
// is concordant, which is the opposite convention and easy to get backwards.
TEST(ActivityLandscapeTest, AMoleculeWithNoDifferentBandNeighbourIsConcordant) {
    OECluster::DenseStorage storage(4);
    FillStorage(storage, {0.8, 0.9, 0.1, 0.2, 0.3, 0.4});
    ActivityLandscapeOptions options;
    options.rmodi_delta = 2.0;

    const ActivityLandscape landscape =
        OECluster::activity_landscape(storage, {0.0, 1.0, 2.0, 3.0}, options);

    EXPECT_NEAR(landscape.activity_stddev, 1.118033988749895, 1e-12);
    EXPECT_DOUBLE_EQ(landscape.rmodi, 0.5);
}

// Equal band minima do not count: the comparison is strict, so a molecule
// whose nearest same-band and nearest different-band neighbours are both at
// distance zero is discordant.
TEST(ActivityLandscapeTest, EqualBandMinimaAtZeroDoNotCount) {
    OECluster::DenseStorage storage(4);
    FillStorage(storage, {0.0, 0.0, 0.7, 0.6, 0.5, 0.4});

    const ActivityLandscape landscape =
        OECluster::activity_landscape(storage, {0.0, 0.0, 5.0, 5.0});

    EXPECT_EQ(landscape.num_zero_distance_pairs, 2u);
    EXPECT_DOUBLE_EQ(landscape.rmodi, 0.5);
}

TEST(ActivityLandscapeTest, EqualBandMinimaAtANonzeroDistanceDoNotCount) {
    OECluster::DenseStorage storage(4);
    FillStorage(storage, {0.3, 0.3, 0.9, 0.8, 0.7, 0.6});

    const ActivityLandscape landscape =
        OECluster::activity_landscape(storage, {0.0, 0.0, 5.0, 5.0});

    // A non-strict comparison would report 0.75 here.
    EXPECT_DOUBLE_EQ(landscape.rmodi, 0.5);
}

TEST(ActivityLandscapeTest, MissingActivitiesLeaveTheirPairsUnscored) {
    OECluster::DenseStorage storage(4);
    // Sample 1 is the dropped one, so scored position p no longer equals
    // original sample p. The distances involving sample 1 are all 0.9 and
    // differ from every scored distance, so reading the matrix at the scored
    // position instead of the original index changes the answer.
    FillStorage(storage, {0.9, 0.5, 0.25, 0.9, 0.9, 0.125});

    const ActivityLandscape landscape = OECluster::activity_landscape(
        storage, {0.0, NOT_A_NUMBER, 1.0, 3.0});

    EXPECT_EQ(landscape.num_samples, 4u);
    EXPECT_EQ(landscape.num_scored, 3u);
    EXPECT_EQ(landscape.num_pairs_scored, 3u);
    EXPECT_DOUBLE_EQ(landscape.max_sali, 16.0);
    EXPECT_DOUBLE_EQ(landscape.mean_sali, 10.0);
}

// Guards 8 and 9 must not reject a degenerate but legal input. Each of these
// returns the documented undefined values rather than throwing.
TEST(ActivityLandscapeTest, DegenerateInputsReturnUndefinedValues) {
    OECluster::DenseStorage storage(2);
    FillStorage(storage, {0.5});

    const ActivityLandscape all_missing = OECluster::activity_landscape(
        storage, {NOT_A_NUMBER, NOT_A_NUMBER});
    EXPECT_EQ(all_missing.num_scored, 0u);
    EXPECT_EQ(all_missing.num_pairs_scored, 0u);
    EXPECT_TRUE(std::isnan(all_missing.cliff_density));
    EXPECT_TRUE(std::isnan(all_missing.rmodi));
    EXPECT_TRUE(std::isnan(all_missing.activity_stddev));

    const ActivityLandscape one_scored =
        OECluster::activity_landscape(storage, {1.0, NOT_A_NUMBER});
    EXPECT_EQ(one_scored.num_scored, 1u);
    EXPECT_EQ(one_scored.num_pairs_scored, 0u);
    EXPECT_TRUE(std::isnan(one_scored.rmodi));
    EXPECT_TRUE(std::isnan(one_scored.activity_stddev));
}

TEST(ActivityLandscapeTest, RejectsANegativeOrNonFiniteDistance) {
    OECluster::DenseStorage negative(3);
    FillStorage(negative, {0.5, -0.25, 0.125});
    try {
        OECluster::activity_landscape(negative, {0.0, 1.0, 3.0});
        FAIL() << "expected a negative distance to be rejected";
    } catch (const std::invalid_argument& error) {
        const std::string message = error.what();
        EXPECT_NE(message.find("activity_landscape"), std::string::npos);
        EXPECT_NE(message.find("0 and 2"), std::string::npos);
    }

    OECluster::DenseStorage infinite(3);
    FillStorage(infinite,
                {0.5, std::numeric_limits<double>::infinity(), 0.125});
    EXPECT_THROW(OECluster::activity_landscape(infinite, {0.0, 1.0, 3.0}),
                 std::invalid_argument);

    // isinf alone would pass the two cases above while letting a NaN distance
    // reach the arithmetic, where every comparison reads false and the sweep
    // reports a silently wrong landscape instead of throwing.
    OECluster::DenseStorage not_a_number(3);
    FillStorage(not_a_number, {0.5, NOT_A_NUMBER, 0.125});
    try {
        OECluster::activity_landscape(not_a_number, {0.0, 1.0, 3.0});
        FAIL() << "expected a NaN distance to be rejected";
    } catch (const std::invalid_argument& error) {
        const std::string message = error.what();
        EXPECT_NE(message.find("activity_landscape"), std::string::npos);
        EXPECT_NE(message.find("0 and 2"), std::string::npos);
    }
}

TEST(ActivityLandscapeTest, RejectsAnOverflowingSali) {
    OECluster::DenseStorage storage(2);
    FillStorage(storage, {1e-200});

    try {
        OECluster::activity_landscape(storage, {0.0, 1e154});
        FAIL() << "expected an overflowing SALI accumulator to be rejected";
    } catch (const std::invalid_argument& error) {
        EXPECT_NE(std::string(error.what()).find("mean_sali"),
                  std::string::npos);
    }
}

TEST(ActivityLandscapeTest, RejectsAnOverflowingSpread) {
    OECluster::DenseStorage storage(2);
    FillStorage(storage, {0.5});

    for (const std::vector<double>& activity :
         {std::vector<double>{-1e300, 1e300}, std::vector<double>{DBL_MAX, DBL_MAX}}) {
        try {
            OECluster::activity_landscape(storage, activity);
            FAIL() << "expected an overflowing spread to be rejected";
        } catch (const std::invalid_argument& error) {
            EXPECT_NE(std::string(error.what()).find("activity_stddev"),
                      std::string::npos);
        }
    }

    // The arrangement that would invert RMODI if the spread guard ever stopped
    // firing: |a - b| overflows to +inf, the band overflows to +inf too, and
    // `delta <= band` reads as true even though 2*DBL_MAX exceeds
    // 1.5*DBL_MAX. A zero distance is what makes it dangerous, because it skips
    // the SALI division where the infinity would otherwise be caught -- the
    // sweep would report a perfectly concordant landscape. The sweep carries no
    // check of its own against this, deliberately: an overflowing difference
    // forces a scale above DBL_MAX/2, whose square cannot be represented, so
    // the spread guard necessarily precedes it. That argument is only worth
    // making if something pins it, which is this block.
    OECluster::DenseStorage coincident(2);
    FillStorage(coincident, {0.0});
    ActivityLandscapeOptions wide_band;
    wide_band.rmodi_delta = 1.5;
    try {
        OECluster::activity_landscape(coincident, {-DBL_MAX, DBL_MAX},
                                      wide_band);
        FAIL() << "expected an overflowing difference to be rejected";
    } catch (const std::invalid_argument& error) {
        EXPECT_NE(std::string(error.what()).find("activity_stddev"),
                  std::string::npos);
    }
}

TEST(ActivityLandscapeTest, RejectsNegativeOrNonFiniteOptions) {
    OECluster::DenseStorage storage(3);
    FillStorage(storage, {0.5, 0.25, 0.125});
    const std::vector<double> activity = {0.0, 1.0, 3.0};

    ActivityLandscapeOptions negative_distance;
    negative_distance.distance_threshold = -0.1;
    try {
        OECluster::activity_landscape(storage, activity, negative_distance);
        FAIL() << "expected a negative distance_threshold to be rejected";
    } catch (const std::invalid_argument& error) {
        EXPECT_NE(std::string(error.what()).find("distance_threshold"),
                  std::string::npos);
    }

    ActivityLandscapeOptions negative_activity;
    negative_activity.activity_threshold = -1.0;
    EXPECT_THROW(
        OECluster::activity_landscape(storage, activity, negative_activity),
        std::invalid_argument);

    ActivityLandscapeOptions nonfinite_delta;
    nonfinite_delta.rmodi_delta = std::numeric_limits<double>::infinity();
    try {
        OECluster::activity_landscape(storage, activity, nonfinite_delta);
        FAIL() << "expected a non-finite rmodi_delta to be rejected";
    } catch (const std::invalid_argument& error) {
        EXPECT_NE(std::string(error.what()).find("rmodi_delta"),
                  std::string::npos);
    }

    // Zero is legal for both thresholds that can meaningfully be zero.
    ActivityLandscapeOptions zeroed;
    zeroed.rmodi_delta = 0.0;
    zeroed.activity_threshold = 0.0;
    EXPECT_NO_THROW(OECluster::activity_landscape(storage, activity, zeroed));
}

TEST(ActivityLandscapeTest, RejectsAnEmptyOrMismatchedActivity) {
    OECluster::DenseStorage storage(3);
    FillStorage(storage, {0.5, 0.25, 0.125});

    try {
        OECluster::activity_landscape(storage, {});
        FAIL() << "expected an empty activity to be rejected";
    } catch (const std::invalid_argument& error) {
        EXPECT_NE(std::string(error.what()).find("non-empty"),
                  std::string::npos);
    }

    try {
        OECluster::activity_landscape(storage, {1.0, 2.0});
        FAIL() << "expected a length mismatch to be rejected";
    } catch (const std::invalid_argument& error) {
        const std::string message = error.what();
        EXPECT_NE(message.find("activity has 2 entries"), std::string::npos);
        EXPECT_NE(message.find("the storage has 3 samples"), std::string::npos);
    }

    // Both directions of the mismatch, because the check is a comparison and a
    // one-sided test would pass against `activity.size() < expected`.
    try {
        OECluster::activity_landscape(storage, {1.0, 2.0, 3.0, 4.0});
        FAIL() << "expected an overlong activity to be rejected";
    } catch (const std::invalid_argument& error) {
        const std::string message = error.what();
        EXPECT_NE(message.find("activity has 4 entries"), std::string::npos);
        EXPECT_NE(message.find("the storage has 3 samples"), std::string::npos);
    }

    // An empty activity against empty storage takes the non-empty refusal, not
    // the cardinality one: the emptiness check runs first, so "0 entries and 0
    // samples" never reads as agreement.
    OECluster::DenseStorage empty(0);
    try {
        OECluster::activity_landscape(empty, {});
        FAIL() << "expected an empty activity to be rejected";
    } catch (const std::invalid_argument& error) {
        EXPECT_NE(std::string(error.what()).find("non-empty"),
                  std::string::npos);
    }
}

TEST(ActivityLandscapeTest, RefusesSparseStorage) {
    OECluster::SparseStorage sparse(3, 0.5);

    try {
        OECluster::activity_landscape(sparse, {0.0, 1.0, 3.0});
        FAIL() << "expected SparseStorage to be refused";
    } catch (const std::invalid_argument& error) {
        const std::string message = error.what();
        EXPECT_NE(message.find("activity_landscape"), std::string::npos);
        EXPECT_NE(message.find("SparseStorage"), std::string::npos);
    }
}

// Every reported value must be identical bit for bit across thread counts.
// The floating-point sums are the fragile part: they are only invariant
// because the per-row partials are combined in ascending row order after the
// join, never inside a worker.
TEST(ActivityLandscapeTest, IsInvariantUnderTheThreadCount) {
    constexpr std::size_t SAMPLES = 200;
    OECluster::DenseStorage storage(SAMPLES);
    for (std::size_t i = 0; i < SAMPLES; ++i) {
        for (std::size_t j = i + 1; j < SAMPLES; ++j) {
            storage.Set(i, j, static_cast<double>((i * 37 + j * 11) % 100) / 100.0);
        }
    }
    std::vector<double> activity(SAMPLES);
    for (std::size_t i = 0; i < SAMPLES; ++i) {
        activity[i] = static_cast<double>((i * 13) % 29) / 7.0;
    }

    ActivityLandscapeOptions single;
    single.num_threads = 1;
    const ActivityLandscape reference =
        OECluster::activity_landscape(storage, activity, single);

    for (const std::size_t threads : {2u, 4u, 8u}) {
        ActivityLandscapeOptions options;
        options.num_threads = threads;
        const ActivityLandscape landscape =
            OECluster::activity_landscape(storage, activity, options);

        EXPECT_EQ(landscape.num_cliffs, reference.num_cliffs);
        EXPECT_EQ(landscape.num_zero_distance_pairs,
                  reference.num_zero_distance_pairs);
        EXPECT_EQ(landscape.max_sali, reference.max_sali);
        EXPECT_EQ(landscape.mean_sali, reference.mean_sali);
        EXPECT_EQ(landscape.rmodi, reference.rmodi);
        EXPECT_EQ(landscape.cliff_density, reference.cliff_density);
    }
}

// num_threads is a size_t, so "more threads than the machine could ever run"
// is a value a caller can pass. The cap at the row count is what keeps
// ThreadPool from being asked to create that many OS threads, and the answer
// must still be the single-threaded one. The chunk-size division is safe for a
// separate reason: the thread count it divides by is itself bounded by the row
// count, not by num_threads.
TEST(ActivityLandscapeTest, CapsAnAbsurdThreadCount) {
    OECluster::DenseStorage storage(4);
    FillStorage(storage, {0.1, 0.2, 0.3, 0.15, 0.25, 0.35});
    const std::vector<double> activity = {1.0, 2.0, 4.0, 8.0};

    ActivityLandscapeOptions single;
    single.num_threads = 1;
    const ActivityLandscape reference =
        OECluster::activity_landscape(storage, activity, single);

    ActivityLandscapeOptions absurd;
    absurd.num_threads = std::size_t{1} << 61;
    const ActivityLandscape landscape =
        OECluster::activity_landscape(storage, activity, absurd);

    EXPECT_EQ(landscape.num_cliffs, reference.num_cliffs);
    EXPECT_EQ(landscape.max_sali, reference.max_sali);
    EXPECT_EQ(landscape.mean_sali, reference.mean_sali);
    EXPECT_EQ(landscape.rmodi, reference.rmodi);
}

TEST(ModelabilityTest, ScoresATwoClassHandFixture) {
    OECluster::DenseStorage storage(4);
    FillStorage(storage, {0.5, 0.1, 0.6, 0.7, 0.8, 0.2});

    const Modelability model =
        OECluster::modelability(storage, {"A", "A", "B", "B"});

    EXPECT_EQ(model.num_samples, 4u);
    EXPECT_EQ(model.num_scored, 4u);
    EXPECT_EQ(model.num_classes, 2u);
    EXPECT_DOUBLE_EQ(model.modi, 0.5);

    ASSERT_EQ(model.classes.size(), 2u);
    EXPECT_EQ(model.classes[0].label, "A");
    EXPECT_EQ(model.classes[0].num_members, 2u);
    EXPECT_DOUBLE_EQ(model.classes[0].fraction_same_class, 0.5);
    EXPECT_EQ(model.classes[1].label, "B");
    EXPECT_DOUBLE_EQ(model.classes[1].fraction_same_class, 0.5);
}

// The only fixture that pins exact per-class fractions across more than two
// classes. It does not pin MODI as the unweighted mean over classes, and
// could not: every class here has exactly two members, so weighting by
// membership gives the same 0.75. ResolvesTiesToTheLowestScoredIndex
// (unweighted 0.25 against a weighted 1/3) and ExcludesSamplesWithNoClass
// (0.5 against 2/3) are where that is pinned.
TEST(ModelabilityTest, ScoresAFourClassHandFixture) {
    constexpr std::size_t SAMPLES = 8;
    OECluster::DenseStorage storage(SAMPLES);
    for (std::size_t i = 0; i < SAMPLES; ++i) {
        for (std::size_t j = i + 1; j < SAMPLES; ++j) {
            storage.Set(i, j, 0.9);
        }
    }
    storage.Set(0, 1, 0.1);
    storage.Set(2, 3, 0.3);
    storage.Set(2, 4, 0.2);
    storage.Set(4, 5, 0.4);
    storage.Set(6, 7, 0.15);

    const Modelability model = OECluster::modelability(
        storage, {"A", "A", "B", "B", "C", "C", "D", "D"});

    EXPECT_EQ(model.num_classes, 4u);
    ASSERT_EQ(model.classes.size(), 4u);
    EXPECT_DOUBLE_EQ(model.classes[0].fraction_same_class, 1.0);
    EXPECT_DOUBLE_EQ(model.classes[1].fraction_same_class, 0.5);
    EXPECT_DOUBLE_EQ(model.classes[2].fraction_same_class, 0.5);
    EXPECT_DOUBLE_EQ(model.classes[3].fraction_same_class, 1.0);
    EXPECT_DOUBLE_EQ(model.modi, 0.75);
}

TEST(ModelabilityTest, ReportsOneOnAPerfectSeparation) {
    OECluster::DenseStorage storage(4);
    FillStorage(storage, {0.1, 0.9, 0.9, 0.9, 0.9, 0.1});

    const Modelability model =
        OECluster::modelability(storage, {"A", "A", "B", "B"});

    EXPECT_DOUBLE_EQ(model.modi, 1.0);
}

// A molecule is never its own nearest neighbour. With every off-diagonal
// distance equal, the diagonal's zero would otherwise make every molecule
// concordant and report 1.0 instead of 0.0.
TEST(ModelabilityTest, NeverReadsTheDiagonal) {
    OECluster::DenseStorage storage(3);
    FillStorage(storage, {1.0, 1.0, 1.0});

    const Modelability model = OECluster::modelability(storage, {"A", "B", "B"});

    EXPECT_DOUBLE_EQ(model.modi, 0.0);
    EXPECT_DOUBLE_EQ(model.classes[0].fraction_same_class, 0.0);
    EXPECT_DOUBLE_EQ(model.classes[1].fraction_same_class, 0.0);
}

TEST(ModelabilityTest, ResolvesTiesToTheLowestScoredIndex) {
    OECluster::DenseStorage storage(3);
    FillStorage(storage, {0.5, 0.5, 0.9});

    const Modelability model = OECluster::modelability(storage, {"A", "B", "A"});

    // Breaking the tie towards the highest index would report 0.5.
    EXPECT_DOUBLE_EQ(model.modi, 0.25);
    EXPECT_DOUBLE_EQ(model.classes[0].fraction_same_class, 0.5);
    EXPECT_DOUBLE_EQ(model.classes[1].fraction_same_class, 0.0);
}

TEST(ModelabilityTest, ExcludesSamplesWithNoClass) {
    OECluster::DenseStorage storage(4);
    FillStorage(storage, {0.1, 0.2, 0.9, 0.9, 0.9, 0.3});

    const Modelability model =
        OECluster::modelability(storage, {"A", "", "A", "B"});

    EXPECT_EQ(model.num_samples, 4u);
    EXPECT_EQ(model.num_scored, 3u);
    EXPECT_EQ(model.num_classes, 2u);
    ASSERT_EQ(model.classes.size(), 2u);
    EXPECT_EQ(model.classes[0].num_members, 2u);
    EXPECT_DOUBLE_EQ(model.classes[0].fraction_same_class, 1.0);
    // The B row's label is read from the original sample index, which is 3
    // rather than the scored position 2. Reading it at the scored position
    // reports original sample 2's class, "A".
    EXPECT_EQ(model.classes[1].label, "B");
    EXPECT_DOUBLE_EQ(model.classes[1].fraction_same_class, 0.0);
    EXPECT_DOUBLE_EQ(model.modi, 0.5);
}

// indices[p] equals p unless something was dropped, so on every fixture with
// no empty class string the sweep reads the same matrix whether it indexes by
// scored position or by original sample. ExcludesSamplesWithNoClass has the
// drop but cannot see the difference: both readings hand every scored sample a
// nearest neighbour of the same class, so none of its numbers move. Dropping a
// middle sample here separates the two readings, 0.75 against 0.25.
TEST(ModelabilityTest, MapsScoredPositionsToSampleIndices) {
    OECluster::DenseStorage storage(5);
    FillStorage(storage, {0.5, 0.9, 0.2, 0.3, 0.1, 0.5, 0.9, 0.1, 0.5, 0.1});

    const Modelability model =
        OECluster::modelability(storage, {"A", "A", "", "B", "B"});

    EXPECT_EQ(model.num_samples, 5u);
    EXPECT_EQ(model.num_scored, 4u);
    EXPECT_EQ(model.num_classes, 2u);

    ASSERT_EQ(model.classes.size(), 2u);
    EXPECT_EQ(model.classes[0].label, "A");
    EXPECT_EQ(model.classes[0].num_members, 2u);
    EXPECT_DOUBLE_EQ(model.classes[0].fraction_same_class, 0.5);
    // Sample 2 is the dropped one, so a label read at the scored position
    // reports "" here instead of "B".
    EXPECT_EQ(model.classes[1].label, "B");
    EXPECT_EQ(model.classes[1].num_members, 2u);
    EXPECT_DOUBLE_EQ(model.classes[1].fraction_same_class, 1.0);
    EXPECT_DOUBLE_EQ(model.modi, 0.75);
}

// The same mapping, witnessed from a drop that lands before the first scored
// sample. That is the sharper case for the label lookup: a leading drop shifts
// every scored position, so reading the labels at the scored position corrupts
// both rows rather than the one the middle-drop fixture above corrupts.
//
// classes[0] is "B", not "A", because the order is first appearance among the
// scored samples and the earliest of those is original sample 1.
TEST(ModelabilityTest, MapsScoredPositionsWhenTheFirstSampleIsDropped) {
    OECluster::DenseStorage storage(4);
    FillStorage(storage, {0.1, 0.8, 0.7, 0.5, 0.9, 0.2});

    const Modelability model =
        OECluster::modelability(storage, {"", "B", "A", "A"});

    EXPECT_EQ(model.num_samples, 4u);
    EXPECT_EQ(model.num_scored, 3u);
    EXPECT_EQ(model.num_classes, 2u);

    ASSERT_EQ(model.classes.size(), 2u);
    EXPECT_EQ(model.classes[0].label, "B");
    EXPECT_EQ(model.classes[0].num_members, 1u);
    EXPECT_DOUBLE_EQ(model.classes[0].fraction_same_class, 0.0);
    EXPECT_EQ(model.classes[1].label, "A");
    EXPECT_EQ(model.classes[1].num_members, 2u);
    EXPECT_DOUBLE_EQ(model.classes[1].fraction_same_class, 1.0);
    EXPECT_DOUBLE_EQ(model.modi, 0.5);
}

// With one class nobody has a neighbour that could differ, so concordance is
// not a question that has an answer.
TEST(ModelabilityTest, IsUndefinedBelowTwoClasses) {
    OECluster::DenseStorage storage(3);
    FillStorage(storage, {0.5, 0.25, 0.125});

    const Modelability single =
        OECluster::modelability(storage, {"A", "A", "A"});
    EXPECT_EQ(single.num_classes, 1u);
    EXPECT_TRUE(std::isnan(single.modi));
    ASSERT_EQ(single.classes.size(), 1u);
    EXPECT_EQ(single.classes[0].num_members, 3u);
    EXPECT_TRUE(std::isnan(single.classes[0].fraction_same_class));

    const Modelability none = OECluster::modelability(storage, {"", "", ""});
    EXPECT_EQ(none.num_scored, 0u);
    EXPECT_EQ(none.num_classes, 0u);
    EXPECT_TRUE(none.classes.empty());
    EXPECT_TRUE(std::isnan(none.modi));
}

// The pair indices are sorted before they are formatted, so the message reads
// the same whichever of the two rows reached the bad entry first. Without that
// the assertion below is a coin flip on the thread schedule.
TEST(ModelabilityTest, RejectsANegativeOrNonFiniteDistance) {
    OECluster::DenseStorage infinite(3);
    FillStorage(infinite,
                {0.5, std::numeric_limits<double>::infinity(), 0.125});

    try {
        OECluster::modelability(infinite, {"A", "B", "B"});
        FAIL() << "expected a non-finite distance to be rejected";
    } catch (const std::invalid_argument& error) {
        const std::string message = error.what();
        EXPECT_NE(message.find("modelability"), std::string::npos);
        EXPECT_NE(message.find("0 and 2"), std::string::npos);
        EXPECT_EQ(message.find("2 and 0"), std::string::npos);
    }

    OECluster::DenseStorage negative(3);
    FillStorage(negative, {0.5, -0.25, 0.125});

    try {
        OECluster::modelability(negative, {"A", "B", "B"});
        FAIL() << "expected a negative distance to be rejected";
    } catch (const std::invalid_argument& error) {
        const std::string message = error.what();
        EXPECT_NE(message.find("0 and 2"), std::string::npos);
        EXPECT_NE(message.find("finite and non-negative"), std::string::npos);
    }

    // A NaN is the case a guard written as isinf would let through, and it is
    // the one that fails silently: every comparison against a NaN reads false,
    // so the sweep would keep an earlier neighbour and return a plausible
    // number over a corrupt matrix instead of refusing it.
    OECluster::DenseStorage not_a_number(3);
    FillStorage(not_a_number, {0.5, NOT_A_NUMBER, 0.125});

    try {
        OECluster::modelability(not_a_number, {"A", "B", "B"});
        FAIL() << "expected a NaN distance to be rejected";
    } catch (const std::invalid_argument& error) {
        const std::string message = error.what();
        EXPECT_NE(message.find("modelability"), std::string::npos);
        EXPECT_NE(message.find("0 and 2"), std::string::npos);
    }

    // The message names original sample indices, which is what a caller needs
    // to locate the entry in their own matrix. A dropped sample is the only
    // thing that makes those differ from the scored positions, and no fixture
    // above has one: sample 0 is unannotated here, so the corrupt entry
    // between samples 1 and 2 sits at scored positions 0 and 1. Both surviving
    // rows read it and both throw, but bad_distance sorts its pair, so the
    // string is the same whichever worker wins.
    OECluster::DenseStorage dropped(3);
    FillStorage(dropped, {0.5, 0.5, std::numeric_limits<double>::infinity()});

    try {
        OECluster::modelability(dropped, {"", "A", "B"});
        FAIL() << "expected a non-finite distance to be rejected";
    } catch (const std::invalid_argument& error) {
        const std::string message = error.what();
        EXPECT_NE(message.find("1 and 2"), std::string::npos);
        EXPECT_EQ(message.find("0 and 1"), std::string::npos);
    }
}

// The single-class path returns before the concordance sweep, so it is the one
// shape that could skip validation entirely. It must refuse on the same terms
// as every other shape rather than hand back NaNs over a corrupt matrix.
TEST(ModelabilityTest, RejectsABadDistanceEvenWithOneClass) {
    OECluster::DenseStorage infinite(3);
    FillStorage(infinite,
                {0.5, std::numeric_limits<double>::infinity(), 0.125});

    try {
        OECluster::modelability(infinite, {"A", "A", "A"});
        FAIL() << "expected a non-finite distance to be rejected";
    } catch (const std::invalid_argument& error) {
        EXPECT_NE(std::string(error.what()).find("0 and 2"), std::string::npos);
    }

    // The serial path needs the NaN case for the same reason the parallel one
    // does: an isinf-only guard would return NaNs over a corrupt matrix rather
    // than refuse it, and this is the shape where nothing else reads a
    // distance.
    OECluster::DenseStorage not_a_number(3);
    FillStorage(not_a_number, {0.5, NOT_A_NUMBER, 0.125});

    try {
        OECluster::modelability(not_a_number, {"A", "A", "A"});
        FAIL() << "expected a NaN distance to be rejected";
    } catch (const std::invalid_argument& error) {
        EXPECT_NE(std::string(error.what()).find("0 and 2"), std::string::npos);
    }

    // A finite negative distance is the guard's other clause, and nothing else
    // reaches it here: the only other negative fixture for modelability is
    // annotated with two classes, so the concordance sweep refuses it before
    // the serial scan is ever entered. Without this case the serial scan could
    // lose its negative test and a one-class call over a negative matrix would
    // return the same NaN that a legitimate single-class call returns.
    OECluster::DenseStorage negative(3);
    FillStorage(negative, {0.5, -0.25, 0.125});

    try {
        OECluster::modelability(negative, {"A", "A", "A"});
        FAIL() << "expected a negative distance to be rejected";
    } catch (const std::invalid_argument& error) {
        const std::string message = error.what();
        EXPECT_NE(message.find("0 and 2"), std::string::npos);
        EXPECT_NE(message.find("finite and non-negative"), std::string::npos);
    }

    // The dropped sample takes its distances out of the scan with it: sample 1
    // is unannotated, so the infinity at (0, 2) is still the pair that bites,
    // while a corrupt entry touching only sample 1 would not be read at all.
    OECluster::DenseStorage unreachable(3);
    FillStorage(unreachable,
                {std::numeric_limits<double>::infinity(), 0.25, 0.125});
    const Modelability dropped =
        OECluster::modelability(unreachable, {"A", "", "A"});
    EXPECT_EQ(dropped.num_scored, 2u);
    EXPECT_TRUE(std::isnan(dropped.modi));

    // The same index mapping as the parallel path, on the serial scan. A drop
    // is the only shape where an original sample index and a scored position
    // can differ, and the throwing cases above have none: with sample 0
    // unannotated, the corrupt entry is between samples 1 and 2 and must be
    // named that way rather than as scored positions 0 and 1.
    OECluster::DenseStorage shifted(3);
    FillStorage(shifted, {0.5, 0.5, std::numeric_limits<double>::infinity()});

    try {
        OECluster::modelability(shifted, {"", "A", "A"});
        FAIL() << "expected a non-finite distance to be rejected";
    } catch (const std::invalid_argument& error) {
        const std::string message = error.what();
        EXPECT_NE(message.find("1 and 2"), std::string::npos);
        EXPECT_EQ(message.find("0 and 1"), std::string::npos);
    }
}

TEST(ModelabilityTest, RejectsAnEmptyOrMismatchedAnnotation) {
    OECluster::DenseStorage storage(3);
    FillStorage(storage, {0.5, 0.25, 0.125});

    try {
        OECluster::modelability(storage, {});
        FAIL() << "expected an empty annotation to be rejected";
    } catch (const std::invalid_argument& error) {
        EXPECT_NE(std::string(error.what()).find("non-empty"),
                  std::string::npos);
    }

    // Both directions of the cardinality check, since a size comparison
    // written with the wrong operator passes one of them.
    try {
        OECluster::modelability(storage, {"A", "B"});
        FAIL() << "expected a short annotation to be rejected";
    } catch (const std::invalid_argument& error) {
        const std::string message = error.what();
        EXPECT_NE(message.find("activity_classes has 2 entries"),
                  std::string::npos);
        EXPECT_NE(message.find("the storage has 3 samples"), std::string::npos);
    }

    try {
        OECluster::modelability(storage, {"A", "B", "A", "B"});
        FAIL() << "expected an overlong annotation to be rejected";
    } catch (const std::invalid_argument& error) {
        const std::string message = error.what();
        EXPECT_NE(message.find("activity_classes has 4 entries"),
                  std::string::npos);
        EXPECT_NE(message.find("the storage has 3 samples"), std::string::npos);
    }

    // An empty annotation against empty storage takes the non-empty refusal,
    // not the cardinality one: the emptiness check runs first, so "0 entries
    // and 0 samples" never reads as agreement.
    OECluster::DenseStorage empty(0);
    try {
        OECluster::modelability(empty, {});
        FAIL() << "expected an empty annotation to be rejected";
    } catch (const std::invalid_argument& error) {
        EXPECT_NE(std::string(error.what()).find("non-empty"),
                  std::string::npos);
    }
}

TEST(ModelabilityTest, RefusesSparseStorage) {
    OECluster::SparseStorage sparse(3, 0.5);

    try {
        OECluster::modelability(sparse, {"A", "B", "B"});
        FAIL() << "expected SparseStorage to be refused";
    } catch (const std::invalid_argument& error) {
        const std::string message = error.what();
        EXPECT_NE(message.find("modelability"), std::string::npos);
        EXPECT_NE(message.find("SparseStorage"), std::string::npos);
    }
}

TEST(ModelabilityTest, IsInvariantUnderTheThreadCount) {
    constexpr std::size_t SAMPLES = 200;
    OECluster::DenseStorage storage(SAMPLES);
    for (std::size_t i = 0; i < SAMPLES; ++i) {
        for (std::size_t j = i + 1; j < SAMPLES; ++j) {
            storage.Set(i, j, static_cast<double>((i * 37 + j * 11) % 100) / 100.0);
        }
    }
    // The class must not be a function of i mod 4. A stored distance is zero
    // when 37p + 11q is a multiple of 100 for the ordered pair (p, q) the loop
    // above filled, which puts row i's zero-distance partner at 33i modulo 100
    // when the partner is the larger index and at 97i modulo 100 when it is the
    // smaller. Both multipliers are congruent to 1 modulo 4 and 100 is a
    // multiple of 4, so either way the partner falls in row i's own residue
    // class modulo 4. A class read off i mod 4 therefore makes every row
    // concordant and drives modi and all four fractions to exactly 1.0, which
    // no longer tells a correct threaded sweep apart from one that never
    // compared a neighbour's label. Dividing by 7 before the stride breaks
    // that alignment.
    std::vector<std::string> classes(SAMPLES);
    for (std::size_t i = 0; i < SAMPLES; ++i) {
        classes[i] = "class" + std::to_string(((i / 7) * 3) % 4);
    }

    ModelabilityOptions single;
    single.num_threads = 1;
    const Modelability reference =
        OECluster::modelability(storage, classes, single);

    // Pin the single-threaded answer before comparing the others against it.
    // Without this the loop below only establishes that the threaded paths
    // agree with each other, which a branch that never looked at the
    // neighbour's label would also satisfy.
    EXPECT_EQ(reference.modi, 0.24894108586830957);
    ASSERT_EQ(reference.classes.size(), 4u);
    EXPECT_EQ(reference.classes[0].label, "class0");
    EXPECT_EQ(reference.classes[0].num_members, 53u);
    EXPECT_EQ(reference.classes[0].fraction_same_class, 0.30188679245283018);
    EXPECT_EQ(reference.classes[1].label, "class3");
    EXPECT_EQ(reference.classes[1].num_members, 49u);
    EXPECT_EQ(reference.classes[1].fraction_same_class, 0.32653061224489793);
    EXPECT_EQ(reference.classes[2].label, "class2");
    EXPECT_EQ(reference.classes[2].num_members, 49u);
    EXPECT_EQ(reference.classes[2].fraction_same_class, 0.16326530612244897);
    EXPECT_EQ(reference.classes[3].label, "class1");
    EXPECT_EQ(reference.classes[3].num_members, 49u);
    EXPECT_EQ(reference.classes[3].fraction_same_class, 0.20408163265306123);

    for (const std::size_t threads : {2u, 4u, 8u}) {
        ModelabilityOptions options;
        options.num_threads = threads;
        const Modelability model =
            OECluster::modelability(storage, classes, options);

        EXPECT_EQ(model.modi, reference.modi);
        ASSERT_EQ(model.classes.size(), reference.classes.size());
        for (std::size_t k = 0; k < model.classes.size(); ++k) {
            EXPECT_EQ(model.classes[k].label, reference.classes[k].label);
            EXPECT_EQ(model.classes[k].fraction_same_class,
                      reference.classes[k].fraction_same_class);
        }
    }
}

// The same cap as `activity_landscape`, tested the same way: an absurd
// num_threads must not be handed to ThreadPool as a thread count, and the
// answer must not change.
TEST(ModelabilityTest, CapsAnAbsurdThreadCount) {
    OECluster::DenseStorage storage(4);
    FillStorage(storage, {0.1, 0.8, 0.9, 0.7, 0.6, 0.2});
    const std::vector<std::string> classes = {"A", "A", "B", "B"};

    ModelabilityOptions single;
    single.num_threads = 1;
    const Modelability reference =
        OECluster::modelability(storage, classes, single);

    ModelabilityOptions absurd;
    absurd.num_threads = std::size_t{1} << 61;
    const Modelability model =
        OECluster::modelability(storage, classes, absurd);

    EXPECT_EQ(model.modi, reference.modi);
    ASSERT_EQ(model.classes.size(), reference.classes.size());
    for (std::size_t k = 0; k < model.classes.size(); ++k) {
        EXPECT_EQ(model.classes[k].fraction_same_class,
                  reference.classes[k].fraction_same_class);
    }
}
