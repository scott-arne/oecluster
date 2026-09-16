#include <gtest/gtest.h>

#include <cfloat>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "oecluster/clustering/ClusterTypes.h"
#include "oecluster/clustering/PartitionAgreement.h"
#include "../../src/clustering/ActivityMetrics.h"
#include "../../src/clustering/ContingencyTable.h"

namespace {

using OECluster::ClusterLabel;
using OECluster::NoiseHandling;

constexpr double NOT_A_NUMBER = std::numeric_limits<double>::quiet_NaN();

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
// exactly what must not happen: §3.3 puts this input outside the supported
// domain, sums_of_squares throws on it, and the two helpers have to agree about
// where the domain ends. This is the unit-level half of
// ActivityLandscapeTest.RejectsAnOverflowingSpread; without it, a scaling
// change can silently widen the domain and only a much later integration test
// would notice.
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

    // The two general invariants, at the two different strengths §5.2
    // guarantees them. The bounds are exact because between is clamped into
    // them; the additivity is not, because between is derived by a subtraction
    // that rounds, so asserting bitwise equality here would be asserting
    // something the derivation does not deliver. Checking these on the hand
    // fixture is what makes them invariants rather than three lucky numbers:
    // the expected values above would still pass if the clamp were dropped.
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
// fixture whose true value is 18/31. Both fixtures are inside §3.3's domain:
// every value, difference and sum here is a finite double. The expected
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
