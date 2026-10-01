#include <gtest/gtest.h>

#include <cfloat>
#include <cmath>
#include <cstdint>
#include <limits>
#include <random>
#include <stdexcept>
#include <vector>

#include "../../src/clustering/ClusterMetrics.h"
#include "../../src/clustering/ExactMedian.h"

using OECluster::detail::ExactMedian;
using OECluster::detail::median_key;
using OECluster::detail::median_value;

namespace {

// The reference: zeros normalized as the selector promises, then the
// sort-based median the report used before.
double Oracle(std::vector<double> values) {
    if (values.empty()) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    for (double& v : values) {
        if (v == 0.0) {
            v = 0.0;
        }
    }
    const double median = OECluster::detail::median_distance(values);
    return median == 0.0 ? 0.0 : median;
}

double RunMedian(const std::vector<double>& values, size_t budget,
                 size_t* passes = nullptr) {
    ExactMedian median(values.size(), budget);
    size_t count = 0;
    while (!median.Done()) {
        for (const double v : values) {
            median.Visit(v);
        }
        median.EndPass();
        ++count;
    }
    if (passes != nullptr) {
        *passes = count;
    }
    return median.Result();
}

bool SameBits(double a, double b) {
    if (std::isnan(a) || std::isnan(b)) {
        return std::isnan(a) && std::isnan(b);
    }
    return a == b && std::signbit(a) == std::signbit(b);
}

void ExpectAllBudgets(const std::vector<double>& values) {
    const size_t n = values.size();
    const double expected = Oracle(values);
    const size_t budgets[] = {0, n == 0 ? 0 : n - 1, n, n + 1, SIZE_MAX};
    for (const size_t budget : budgets) {
        const double got = RunMedian(values, budget);
        EXPECT_TRUE(SameBits(got, expected))
            << "budget " << budget << ": got " << got << ", expected " << expected;
    }
}

}  // namespace

TEST(ExactMedianTest, KeysOrderLikeValuesAndRoundTrip) {
    const std::vector<double> ascending = {
        -std::numeric_limits<double>::infinity(), -DBL_MAX, -1.5, -DBL_MIN,
        -std::numeric_limits<double>::denorm_min(), 0.0,
        std::numeric_limits<double>::denorm_min(), DBL_MIN, 1.5, DBL_MAX,
        std::numeric_limits<double>::infinity()};
    for (size_t i = 0; i < ascending.size(); ++i) {
        EXPECT_EQ(median_value(median_key(ascending[i])), ascending[i]);
        if (i > 0) {
            EXPECT_LT(median_key(ascending[i - 1]), median_key(ascending[i]));
        }
    }
}

TEST(ExactMedianTest, OddAndEvenCounts) {
    ExpectAllBudgets({3.0, 1.0, 2.0});
    ExpectAllBudgets({4.0, 1.0, 3.0, 2.0});
    ExpectAllBudgets({0.25});
}

TEST(ExactMedianTest, RandomSetsMatchTheOracle) {
    std::mt19937 generator(20261001);
    std::uniform_real_distribution<double> uniform(-1.0, 1.0);
    for (const size_t n : {size_t{1000}, size_t{1001}}) {
        std::vector<double> values(n);
        for (double& v : values) {
            v = uniform(generator);
        }
        ExpectAllBudgets(values);
    }
}

TEST(ExactMedianTest, HeavyTiesAndAllEqual) {
    ExpectAllBudgets({0.5, 0.5, 0.5, 0.25, 0.5, 0.75, 0.5});
    ExpectAllBudgets({0.5, 0.5, 0.25, 0.25});
    ExpectAllBudgets(std::vector<double>(64, 0.375));
}

TEST(ExactMedianTest, EvenRanksSplitAcrossEachDigit) {
    // lo and hi share every digit above d and differ at digit d, so the two
    // target ranks follow different prefixes from that digit on.
    const double base = 0.3;
    for (int d = 0; d < 4; ++d) {
        const double lo = base;
        const double hi = median_value(median_key(base) + (uint64_t{1} << (48 - 16 * d)));
        ASSERT_LT(lo, hi);
        ExpectAllBudgets({lo, lo, hi, hi});
        ExpectAllBudgets({hi, lo, hi, lo, lo, hi});
    }
}

TEST(ExactMedianTest, NegativesAndMixedSigns) {
    ExpectAllBudgets({-0.5, -0.25, -0.75});
    ExpectAllBudgets({-0.5, 0.25, -0.75, 0.5});
    ExpectAllBudgets({-1e-300, 1e-300, -2.0, 2.0, 0.0});
}

TEST(ExactMedianTest, NegativeZeroBecomesPositiveZero) {
    const std::vector<double> values = {-0.0, -0.0, 1.0};
    for (const size_t budget : {size_t{0}, SIZE_MAX}) {
        const double got = RunMedian(values, budget);
        EXPECT_EQ(got, 0.0);
        EXPECT_FALSE(std::signbit(got));
    }
    ExpectAllBudgets({-0.0, 0.0, -0.0, 0.0});
}

TEST(ExactMedianTest, AnAverageThatUnderflowsToZeroIsPositive) {
    // (-denorm_min + 0.0) / 2 rounds to -0.0, so normalizing the inputs alone
    // would not keep the promise.
    const double tiny = std::numeric_limits<double>::denorm_min();
    for (const size_t budget : {size_t{0}, SIZE_MAX}) {
        const double got = RunMedian({-tiny, 0.0}, budget);
        EXPECT_EQ(got, 0.0);
        EXPECT_FALSE(std::signbit(got)) << "budget " << budget;
    }
    std::vector<double> values = {0.0, -tiny};
    EXPECT_FALSE(std::signbit(ExactMedian::OfInPlace(values)));
}

TEST(ExactMedianTest, ExtremeMagnitudes) {
    ExpectAllBudgets({DBL_MAX, std::numeric_limits<double>::denorm_min(),
                      std::numeric_limits<double>::infinity(), -DBL_MAX});
}

TEST(ExactMedianTest, EmptyIsNaNAndAlreadyDone) {
    ExactMedian median(0);
    EXPECT_TRUE(median.Done());
    EXPECT_TRUE(std::isnan(median.Result()));
    EXPECT_THROW(median.Visit(1.0), std::logic_error);
    EXPECT_THROW(median.EndPass(), std::logic_error);
}

TEST(ExactMedianTest, PassCounts) {
    const std::vector<double> values = {0.1, 0.4, 0.2, 0.3, 0.5};
    size_t passes = 0;
    RunMedian(values, values.size(), &passes);
    EXPECT_EQ(passes, 1u);
    RunMedian(values, values.size() - 1, &passes);
    EXPECT_EQ(passes, 4u);
}

TEST(ExactMedianTest, ProtocolViolationsAreLogicErrors) {
    for (const size_t budget : {size_t{0}, SIZE_MAX}) {
        ExactMedian short_pass(3, budget);
        short_pass.Visit(1.0);
        EXPECT_THROW(short_pass.EndPass(), std::logic_error);

        ExactMedian extra(1, budget);
        extra.Visit(1.0);
        EXPECT_THROW(extra.Visit(2.0), std::logic_error);

        ExactMedian early(1, budget);
        EXPECT_THROW(early.Result(), std::logic_error);
    }

    ExactMedian after(1, SIZE_MAX);
    after.Visit(1.0);
    after.EndPass();
    EXPECT_THROW(after.Visit(1.0), std::logic_error);

    // Values that change between radix passes leave a rank with no bucket.
    ExactMedian changed(2, 0);
    changed.Visit(1.0);
    changed.Visit(1.0);
    changed.EndPass();
    changed.Visit(-1.0);
    changed.Visit(-1.0);
    EXPECT_THROW(changed.EndPass(), std::logic_error);
}

TEST(ExactMedianTest, OfInPlaceMatchesTheOracle) {
    std::vector<double> empty;
    EXPECT_TRUE(std::isnan(ExactMedian::OfInPlace(empty)));
    std::vector<double> values = {0.4, -0.0, 0.1, -0.0};
    const double expected = Oracle(values);
    const double got = ExactMedian::OfInPlace(values);
    EXPECT_TRUE(SameBits(got, expected));
    EXPECT_FALSE(std::signbit(got));
}
