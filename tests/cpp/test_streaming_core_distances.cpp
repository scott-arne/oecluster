/**
 * @file test_streaming_core_distances.cpp
 * @brief HDBSCAN core distances from one pass over a comparison.
 */

#include <gtest/gtest.h>

#include <cmath>
#include <cstddef>
#include <limits>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "oecluster/StorageBackend.h"
#include "oecluster/clustering/HDBSCAN.h"

#include "../../src/clustering/StreamingCoreDistances.h"
#include "diversity_test_support.h"
#include "mst_test_support.h"

using namespace OECluster;
using namespace diversity_test;
using namespace mst_test;

namespace {

const double NaN = std::numeric_limits<double>::quiet_NaN();
const double INF = std::numeric_limits<double>::infinity();

// Every (block, block) tile the schedule runs, so pair coverage is checked by
// expanding tiles into item pairs.
std::vector<size_t> PairCoverage(size_t n, const detail::CoreBlockSchedule& schedule) {
    std::vector<size_t> covered(n * (n - 1) / 2, 0);
    for (const auto& tiles : schedule.rounds) {
        std::set<size_t> blocks_in_round;
        for (const auto& [a, b] : tiles) {
            EXPECT_LE(a, b);
            EXPECT_TRUE(blocks_in_round.insert(a).second) << "block " << a << " twice";
            if (b != a) {
                EXPECT_TRUE(blocks_in_round.insert(b).second) << "block " << b << " twice";
            }
            for (size_t i = schedule.bounds[a]; i < schedule.bounds[a + 1]; ++i) {
                const size_t first = (a == b) ? i + 1 : schedule.bounds[b];
                for (size_t j = first; j < schedule.bounds[b + 1]; ++j) {
                    ++covered[n * i - i * (i + 1) / 2 + j - i - 1];
                }
            }
        }
    }
    return covered;
}

}  // namespace

TEST(StreamingCoreDistancesTest, TheScheduleCoversEveryPairOnceWithDisjointRounds) {
    std::vector<size_t> sizes;
    for (size_t n = 1; n <= 70; ++n) {
        sizes.push_back(n);
    }
    for (size_t n : {127, 128, 257, 300, 511}) {
        sizes.push_back(n);
    }
    for (size_t n : sizes) {
        for (size_t participants : {1, 2, 3, 4, 7, 16}) {
            const detail::CoreBlockSchedule schedule =
                detail::core_block_schedule(n, participants);
            const size_t blocks = std::min(n, 4 * participants);
            ASSERT_EQ(schedule.bounds.size(), blocks + 1);
            EXPECT_EQ(schedule.bounds.front(), 0u);
            EXPECT_EQ(schedule.bounds.back(), n);
            if (n >= 2) {
                EXPECT_EQ(PairCoverage(n, schedule),
                          std::vector<size_t>(n * (n - 1) / 2, 1))
                    << "n=" << n << " participants=" << participants;
            }
        }
    }
    EXPECT_TRUE(detail::core_block_schedule(0, 4).rounds.empty());
}

// The streaming parity test and the legacy HDBSCAN oracle both take their core
// distances from this function, so only literals worked out by hand can catch a
// change common to all of them. Items on a line at 0, 1, 3, 7 and 15; a core
// distance is the distance to the (min_samples - 1)-th nearest neighbor.
TEST(StreamingCoreDistancesTest, TheMatrixPassReproducesHandComputedValues) {
    const std::vector<double> condensed = Positions({0.0, 1.0, 3.0, 7.0, 15.0});
    const DenseStorage storage = MakeStorage(5, condensed);
    const std::vector<std::pair<size_t, std::vector<double>>> expected{
        {1, {0.0, 0.0, 0.0, 0.0, 0.0}},
        {2, {1.0, 1.0, 2.0, 4.0, 8.0}},
        {3, {3.0, 2.0, 3.0, 6.0, 12.0}},
        {5, {15.0, 14.0, 12.0, 8.0, 15.0}}};
    for (const auto& [min_samples, values] : expected) {
        const detail::CoreDistances core =
            detail::matrix_core_distances(storage, min_samples, 1, 0, "test");
        EXPECT_EQ(core.values, values) << "min_samples=" << min_samples;
        EXPECT_EQ(core.max_distance, min_samples == 1 ? 0.0 : 15.0)
            << "min_samples=" << min_samples;
    }
}

TEST(StreamingCoreDistancesTest, MatchesTheMatrixPassBitForBit) {
    struct Case {
        size_t n;
        std::vector<double> condensed;
    };
    const std::vector<Case> cases{{2, Scrambled(2)},
                                  {9, ScrambledSixths(9)},
                                  {25, Hashed(25)},
                                  {61, Quantized(61, 3, 4)},
                                  {130, Quantized(130, 9, 1000)}};
    for (const Case& c : cases) {
        const DenseStorage storage = MakeStorage(c.n, c.condensed);
        for (size_t min_samples = 1; min_samples <= c.n; ++min_samples) {
            const detail::CoreDistances expected =
                detail::matrix_core_distances(storage, min_samples, 1, 0, "test");
            for (size_t threads : {1, 4, 8}) {
                TableComparison comparison(c.n, c.condensed);
                const detail::CoreDistances observed = detail::streaming_core_distances(
                    comparison, min_samples, threads, "test");
                ASSERT_EQ(observed.values, expected.values)
                    << "n=" << c.n << " min_samples=" << min_samples << " threads=" << threads;
                EXPECT_EQ(observed.max_distance, expected.max_distance);
            }
        }
    }
}

TEST(StreamingCoreDistancesTest, ComparesEveryPairExactlyOnce) {
    const size_t n = 53;
    const std::vector<double> condensed = Quantized(n, 4, 7);
    for (size_t threads : {1, 4, 16}) {
        PairCountingComparison comparison(n, condensed);
        detail::streaming_core_distances(comparison, 5, threads, "test");
        EXPECT_EQ(comparison.Counts(), std::vector<size_t>(condensed.size(), 1))
            << "threads=" << threads;
    }
}

TEST(StreamingCoreDistancesTest, AThrowingCompareStopsThePassAndTheNextPassSucceeds) {
    const size_t n = 150;
    const std::vector<double> condensed = Quantized(n, 17, 40);
    const detail::CoreDistances expected = detail::matrix_core_distances(
        MakeStorage(n, condensed), 4, 1, 0, "test");
    for (size_t threads : {1, 4}) {
        ThrowingComparison failing(n, condensed, 20, 130);
        try {
            detail::streaming_core_distances(failing, 4, threads, "test");
            FAIL() << "the exception was lost at threads=" << threads;
        } catch (const std::runtime_error& error) {
            EXPECT_STREQ(error.what(), "ThrowingComparison refused (20, 130)");
        }
        PairCountingComparison healthy(n, condensed);
        EXPECT_EQ(detail::streaming_core_distances(healthy, 4, threads, "test").values,
                  expected.values);
    }
}

TEST(StreamingCoreDistancesTest, MinSamplesOneComparesNothing) {
    const size_t n = 12;
    PairCountingComparison comparison(n, Quantized(n, 4, 7));
    const detail::CoreDistances core =
        detail::streaming_core_distances(comparison, 1, 4, "test");
    EXPECT_EQ(core.values, std::vector<double>(n, 0.0));
    EXPECT_EQ(core.max_distance, 0.0);
    EXPECT_EQ(comparison.Total(), 0u);
}

TEST(StreamingCoreDistancesTest, RefusesAnOutOfRangeMinSamplesBeforeReading) {
    PairCountingComparison comparison(5, Line(5));
    EXPECT_THROW(detail::streaming_core_distances(comparison, 0, 1, "test"),
                 std::invalid_argument);
    EXPECT_THROW(detail::streaming_core_distances(comparison, 6, 1, "test"),
                 std::invalid_argument);
    EXPECT_EQ(comparison.Total(), 0u);
}

TEST(StreamingCoreDistancesTest, RefusesValuesOutsideTheDomain) {
    for (double bad : {NaN, INF, -INF, -0.25}) {
        std::vector<double> condensed = Line(6);
        condensed[7] = bad;  // the pair (1, 4)
        const std::string expected =
            bad < 0.0 && std::isfinite(bad)
                ? "test read a negative distance between items 1 and 4"
                : "test read a non-finite distance between items 1 and 4";
        for (size_t threads : {1, 4}) {
            TableComparison comparison(6, condensed);
            try {
                detail::streaming_core_distances(comparison, 2, threads, "test");
                FAIL() << "accepted " << bad;
            } catch (const std::runtime_error& error) {
                EXPECT_EQ(std::string(error.what()), expected);
            }
        }
        const DenseStorage storage = MakeStorage(6, condensed);
        try {
            detail::matrix_core_distances(storage, 2, 1, 0, "test");
            FAIL() << "matrix accepted " << bad;
        } catch (const std::runtime_error& error) {
            EXPECT_EQ(std::string(error.what()), expected);
        }
    }
}

TEST(StreamingCoreDistancesTest, NegativeZeroBecomesPositiveZero) {
    std::vector<double> condensed = Line(4);
    condensed[0] = -0.0;  // the pair (0, 1)
    TableComparison comparison(4, condensed);
    const detail::CoreDistances streamed =
        detail::streaming_core_distances(comparison, 2, 1, "test");
    const detail::CoreDistances matrix =
        detail::matrix_core_distances(MakeStorage(4, condensed), 2, 1, 0, "test");
    for (const auto* core : {&streamed, &matrix}) {
        EXPECT_EQ(core->values[0], 0.0);
        EXPECT_FALSE(std::signbit(core->values[0]));
        EXPECT_FALSE(std::signbit(core->values[1]));
    }
}

TEST(StreamingCoreDistancesTest, AnOverflowingHeapSizeIsALengthError) {
    HugeComparison huge(std::numeric_limits<size_t>::max() / 4);
    EXPECT_THROW(detail::streaming_core_distances(huge, 5, 1, "test"), std::length_error);
}

TEST(StreamingCoreDistancesTest, TheMatrixChunkSizeChangesNothing) {
    const size_t n = 90;
    const std::vector<double> condensed = Quantized(n, 21, 6);
    const DenseStorage storage = MakeStorage(n, condensed);
    const detail::CoreDistances expected =
        detail::matrix_core_distances(storage, 4, 1, 0, "test");
    for (size_t chunk : {0, 1, 7, 64, 100000}) {
        for (size_t threads : {1, 4}) {
            EXPECT_EQ(detail::matrix_core_distances(storage, 4, threads, chunk, "test").values,
                      expected.values);
        }
    }
}
