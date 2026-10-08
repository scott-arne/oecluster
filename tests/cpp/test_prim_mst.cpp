/**
 * @file test_prim_mst.cpp
 * @brief The shared Prim kernel against the 5.19.0 sequential scan.
 */

#include <gtest/gtest.h>

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <functional>
#include <limits>
#include <optional>
#include <thread>
#include <stdexcept>
#include <string>
#include <vector>

#include "oecluster/StorageBackend.h"
#include "oecluster/clustering/HDBSCAN.h"

#include "../../src/clustering/DistanceAccess.h"
#include "../../src/clustering/PrimMST.h"
#include "../../src/clustering/StepTeam.h"
#include "diversity_test_support.h"
#include "mst_test_support.h"

using namespace OECluster;
using namespace diversity_test;
using namespace mst_test;

namespace {

const double NaN = std::numeric_limits<double>::quiet_NaN();
const double INF = std::numeric_limits<double>::infinity();
const size_t ALL_TEAM = 0;
const size_t ALL_SERIAL = std::numeric_limits<size_t>::max();

struct Fixture {
    size_t n;
    std::vector<double> condensed;
};

// Tie-heavy tables, where the lowest-index rule decides most steps, and
// seeded random levels.
std::vector<Fixture> Fixtures() {
    return {{2, Scrambled(2)},
            {3, Scrambled(3)},
            {7, Scrambled(7)},
            {12, ScrambledSixths(12)},
            {25, Hashed(25)},
            {40, Quantized(40, 3, 5)},
            {97, Quantized(97, 11, 9)},
            {150, Quantized(150, 5, 1000)}};
}

detail::PrimOptions Options(size_t threads, std::optional<size_t> cutoff) {
    detail::PrimOptions options;
    options.num_threads = threads;
    options.serial_cutoff = cutoff;
    options.caller = "test";
    return options;
}

std::vector<std::optional<size_t>> Cutoffs() {
    return {ALL_TEAM, ALL_SERIAL, std::nullopt};
}

}  // namespace

TEST(PrimMSTTest, SingleLinkageMatchesTheLegacyScan) {
    for (const Fixture& fixture : Fixtures()) {
        const DenseStorage storage = MakeStorage(fixture.n, fixture.condensed);
        const auto expected =
            LegacyPrim(fixture.n, fixture.condensed,
                       std::vector<double>(fixture.n, 0.0), 1.0);
        for (size_t threads : {1, 4, 8}) {
            for (const auto& cutoff : Cutoffs()) {
                const auto observed =
                    detail::prim_mst(storage, detail::PrimWeights(), Options(threads, cutoff));
                EXPECT_TRUE(SameEdges(observed, expected))
                    << "n=" << fixture.n << " threads=" << threads;
            }
        }
    }
}

TEST(PrimMSTTest, MutualReachabilityMatchesTheLegacyScan) {
    for (const Fixture& fixture : Fixtures()) {
        const DenseStorage storage = MakeStorage(fixture.n, fixture.condensed);
        for (size_t min_samples : {size_t{1}, size_t{2}, std::min<size_t>(5, fixture.n)}) {
            const std::vector<double> core =
                detail::compute_core_distances(storage, min_samples, 1);
            for (double alpha : {1.0, 0.5, 1.7}) {
                const auto expected = LegacyPrim(fixture.n, fixture.condensed, core, alpha);
                detail::PrimWeights weights;
                weights.core = core;
                weights.alpha = alpha;
                for (bool prune : {false, min_samples > 1}) {
                    weights.prune = prune;
                    for (size_t threads : {1, 4, 8}) {
                        for (const auto& cutoff : Cutoffs()) {
                            const auto observed = detail::prim_mst(
                                storage, weights, Options(threads, cutoff));
                            EXPECT_TRUE(SameEdges(observed, expected))
                                << "n=" << fixture.n << " min_samples=" << min_samples
                                << " alpha=" << alpha << " prune=" << prune
                                << " threads=" << threads;
                        }
                    }
                }
            }
        }
    }
}

TEST(PrimMSTTest, TheComparisonProviderMatchesTheMatrixProvider) {
    for (const Fixture& fixture : Fixtures()) {
        const DenseStorage storage = MakeStorage(fixture.n, fixture.condensed);
        detail::PrimWeights weights;
        weights.core = detail::compute_core_distances(
            storage, std::min<size_t>(3, fixture.n), 1);
        weights.prune = fixture.n >= 3;
        for (bool mutual : {false, true}) {
            const detail::PrimWeights& used = mutual ? weights : detail::PrimWeights();
            const auto expected = detail::prim_mst(storage, used, Options(1, std::nullopt));
            for (size_t threads : {1, 4, 8}) {
                for (const auto& cutoff : Cutoffs()) {
                    TableComparison comparison(fixture.n, fixture.condensed);
                    EXPECT_TRUE(SameEdges(
                        detail::prim_mst(comparison, used, Options(threads, cutoff)),
                        expected))
                        << "n=" << fixture.n << " mutual=" << mutual
                        << " threads=" << threads;
                }
            }
        }
    }
}

TEST(PrimMSTTest, TheComparisonPathClonesOncePerParticipant) {
    const Fixture fixture{40, Quantized(40, 3, 5)};
    TableComparison comparison(fixture.n, fixture.condensed);
    detail::prim_mst(comparison, detail::PrimWeights(), Options(4, ALL_TEAM));
    EXPECT_EQ(comparison.NumClones(), 4u);
}

// Counting the clones does not say they were used: a pass that called one
// clone from every participant would count the same. The hazard the serial
// clone build exists to prevent is a shared instance, so the witness has to be
// which instance served each call.
TEST(PrimMSTTest, ParticipantsReadThroughDifferentClones) {
    const Fixture fixture{150, Quantized(150, 5, 1000)};
    CloneWitnessComparison comparison(fixture.n, fixture.condensed);
    detail::prim_mst(comparison, detail::PrimWeights(), Options(4, ALL_TEAM));
    EXPECT_GE(comparison.ServingClones(), 2u);
    EXPECT_FALSE(comparison.PrototypeServed());
}

TEST(PrimMSTTest, SingleLinkageComparesEveryPairExactlyOnce) {
    const Fixture fixture{60, Quantized(60, 7, 6)};
    for (size_t threads : {1, 4}) {
        for (const auto& cutoff : Cutoffs()) {
            PairCountingComparison comparison(fixture.n, fixture.condensed);
            detail::prim_mst(comparison, detail::PrimWeights(), Options(threads, cutoff));
            EXPECT_EQ(comparison.Counts(), std::vector<size_t>(fixture.condensed.size(), 1));
        }
    }
}

TEST(PrimMSTTest, UnprunedMutualReachabilityComparesEveryPairExactlyOnce) {
    const Fixture fixture{60, Quantized(60, 7, 6)};
    detail::PrimWeights weights;
    weights.core.assign(fixture.n, 0.0);
    PairCountingComparison comparison(fixture.n, fixture.condensed);
    detail::prim_mst(comparison, weights, Options(4, ALL_TEAM));
    EXPECT_EQ(comparison.Counts(), std::vector<size_t>(fixture.condensed.size(), 1));
}

// Two well-separated groups of evenly spaced points: once a candidate's reach
// equals its own core distance it is never compared again, so pruning skips
// pairs without changing an edge. Candidates across the gap keep a large reach
// and are compared at every step, which is why the saving is modest.
TEST(PrimMSTTest, PruningSkipsPairsWithoutChangingTheTree) {
    std::vector<double> x;
    for (int i = 0; i < 40; ++i) {
        x.push_back(i * 0.01);
        x.push_back(100.0 + i * 0.01);
    }
    const std::vector<double> condensed = Positions(x);
    const size_t n = x.size();
    const DenseStorage storage = MakeStorage(n, condensed);
    detail::PrimWeights weights;
    weights.core = detail::compute_core_distances(storage, 5, 1);

    PairCountingComparison unpruned(n, condensed);
    const auto expected = detail::prim_mst(unpruned, weights, Options(4, std::nullopt));
    EXPECT_EQ(unpruned.Total(), condensed.size());

    weights.prune = true;
    PairCountingComparison pruned(n, condensed);
    EXPECT_TRUE(SameEdges(detail::prim_mst(pruned, weights, Options(4, std::nullopt)),
                          expected));
    for (const size_t count : pruned.Counts()) {
        EXPECT_LE(count, 1u);
    }
    EXPECT_LT(pruned.Total(), unpruned.Total());
}

TEST(PrimMSTTest, AThrowingCompareStopsTheRunAndTheNextRunSucceeds) {
    const size_t n = 200;
    const std::vector<double> condensed = Quantized(n, 13, 50);
    const DenseStorage storage = MakeStorage(n, condensed);
    const auto expected = detail::prim_mst(storage, detail::PrimWeights(),
                                           Options(1, std::nullopt));
    for (size_t threads : {1, 4, 8}) {
        ThrowingComparison failing(n, condensed, 70, 140);
        try {
            detail::prim_mst(failing, detail::PrimWeights(), Options(threads, ALL_TEAM));
            FAIL() << "the exception was lost at threads=" << threads;
        } catch (const std::runtime_error& error) {
            EXPECT_STREQ(error.what(), "ThrowingComparison refused (70, 140)");
        }
        PairCountingComparison healthy(n, condensed);
        EXPECT_TRUE(SameEdges(
            detail::prim_mst(healthy, detail::PrimWeights(), Options(threads, ALL_TEAM)),
            expected));
    }
}

TEST(PrimMSTTest, FewerThanTwoItemsHaveNoEdges) {
    for (size_t n : {0, 1}) {
        DenseStorage storage(n);
        EXPECT_TRUE(detail::prim_mst(storage, detail::PrimWeights(), Options(4, ALL_TEAM))
                        .empty());
    }
}

TEST(PrimMSTTest, ANonFiniteDistanceIsRefusedNamingThePair) {
    for (double bad : {NaN, INF, -INF}) {
        std::vector<double> condensed = Line(6);
        condensed[7] = bad;  // the pair (1, 4)
        for (size_t threads : {1, 4}) {
            DenseStorage storage = MakeStorage(6, condensed);
            TableComparison comparison(6, condensed);
            for (bool mutual : {false, true}) {
                detail::PrimWeights weights;
                if (mutual) {
                    weights.core.assign(6, 0.0);
                }
                try {
                    detail::prim_mst(storage, weights, Options(threads, ALL_TEAM));
                    FAIL() << "matrix accepted " << bad;
                } catch (const std::runtime_error& error) {
                    EXPECT_STREQ(error.what(),
                                 "test read a non-finite distance between items 1 and 4");
                }
                EXPECT_THROW(detail::prim_mst(comparison, weights, Options(threads, ALL_TEAM)),
                             std::runtime_error);
            }
        }
    }
}

TEST(PrimMSTTest, MutualReachabilityRefusesANegativeDistance) {
    std::vector<double> condensed = Line(6);
    condensed[7] = -0.5;
    DenseStorage storage = MakeStorage(6, condensed);
    detail::PrimWeights weights;
    weights.core.assign(6, 0.0);
    try {
        detail::prim_mst(storage, weights, Options(1, std::nullopt));
        FAIL() << "accepted a negative distance";
    } catch (const std::runtime_error& error) {
        EXPECT_STREQ(error.what(), "test read a negative distance between items 1 and 4");
    }
    // Single linkage has no self-distance to undercut and accepts it.
    EXPECT_EQ(detail::prim_mst(storage, detail::PrimWeights(), Options(1, std::nullopt))
                  .size(),
              5u);
}

TEST(PrimMSTTest, AnOverflowingQuotientRaisesTheAlphaError) {
    DenseStorage storage = MakeStorage(4, Line(4));
    detail::PrimWeights weights;
    weights.core.assign(4, 0.0);
    weights.alpha = 1e-310;
    try {
        detail::prim_mst(storage, weights, Options(1, std::nullopt));
        FAIL() << "accepted an overflowing quotient";
    } catch (const std::invalid_argument& error) {
        EXPECT_NE(std::string(error.what()).find("test alpha=1e-310"), std::string::npos)
            << error.what();
    }
}

TEST(PrimMSTTest, NegativeZeroIsReadAsPositiveZero) {
    std::vector<double> condensed = Line(5);
    condensed[0] = -0.0;  // the pair (0, 1)
    DenseStorage storage = MakeStorage(5, condensed);
    const auto mst = detail::prim_mst(storage, detail::PrimWeights(), Options(1, std::nullopt));
    ASSERT_EQ(mst[0].next_node, 1u);
    EXPECT_EQ(mst[0].distance, 0.0);
    EXPECT_FALSE(std::signbit(mst[0].distance));
}

TEST(PrimMSTTest, TheDefaultParticipantCountIsCappedAtEight) {
    EXPECT_EQ(detail::prim_participants(0, 100, 14), 8u);
    EXPECT_EQ(detail::prim_participants(0, 100, 4), 4u);
    EXPECT_EQ(detail::prim_participants(0, 100, 0), 1u);
    EXPECT_EQ(detail::prim_participants(0, 3, 14), 3u);
    // An explicit request is honored, clamped only to the item count.
    EXPECT_EQ(detail::prim_participants(14, 100, 4), 14u);
    EXPECT_EQ(detail::prim_participants(14, 10, 4), 10u);
}

// Not a check: the measurement that sets PRIM_SERIAL_CUTOFF_MATRIX. Run it with
//   oecluster_tests --gtest_also_run_disabled_tests
//       --gtest_filter=PrimMSTTest.DISABLED_MeasureMatrixSerialCutoff
// only while the 1-minute load average is below half the core count. Each row
// times one serial step and one team step over R candidates at the default
// participant count, the median of 5, reading a 512 MB buffer through the
// condensed index of a (2^17 + 1)-item matrix so the reads keep a real matrix's
// index arithmetic and cache behavior. The constant is the smallest R from which
// the team is faster at every larger R in the grid, or 2^18 if none.
TEST(PrimMSTTest, DISABLED_MeasureMatrixSerialCutoff) {
    const size_t virtual_n = (size_t{1} << 17) + 1;
    const size_t buffer_size = size_t{1} << 26;
    std::vector<double> buffer(buffer_size);
    for (size_t k = 0; k < buffer_size; ++k) {
        buffer[k] = static_cast<double>((k * 2654435761u) % 1000) / 1000.0;
    }
    const size_t participants =
        detail::prim_participants(0, virtual_n, std::thread::hardware_concurrency());
    detail::StepTeam team(participants);
    const size_t current = virtual_n / 2;
    std::printf("participants=%zu\n%8s %12s %12s\n", participants, "R", "serial_us",
                "team_us");
    for (size_t r = size_t{1} << 6; r <= (size_t{1} << 17); r <<= 1) {
        // Distinct candidates spread over the virtual matrix, skipping the
        // current node: a duplicate would let two chunks write one reach.
        std::vector<size_t> remaining(r);
        const size_t stride = (virtual_n - 1) / r;
        for (size_t k = 0; k < r; ++k) {
            const size_t index = k * stride;
            remaining[k] = index >= current ? index + 1 : index;
        }
        ASSERT_TRUE(std::adjacent_find(remaining.begin(), remaining.end()) ==
                    remaining.end());
        ASSERT_TRUE(std::find(remaining.begin(), remaining.end(), current) ==
                    remaining.end());
        std::vector<double> reach(virtual_n);
        auto scan = [&](size_t begin, size_t end, double& best) {
            for (size_t position = begin; position < end; ++position) {
                const size_t c = remaining[position];
                const size_t index =
                    detail::condensed_index(virtual_n, current, c) % buffer_size;
                const double weight = buffer[index];
                if (weight < reach[c]) {
                    reach[c] = weight;
                }
                best = std::min(best, reach[c]);
            }
        };
        auto median_us = [&](auto&& step) {
            std::vector<double> times;
            for (int rep = 0; rep < 5; ++rep) {
                std::fill(reach.begin(), reach.end(),
                          std::numeric_limits<double>::infinity());
                const auto start = std::chrono::steady_clock::now();
                step();
                times.push_back(std::chrono::duration<double, std::micro>(
                                    std::chrono::steady_clock::now() - start)
                                    .count());
            }
            std::sort(times.begin(), times.end());
            return times[2];
        };
        // Every timed step's best feeds a checksum printed below, so no
        // optimizer can treat the work as dead.
        double checksum = 0.0;
        const double serial = median_us([&] {
            double best = std::numeric_limits<double>::infinity();
            scan(0, r, best);
            checksum += best;
        });
        std::atomic<size_t> next{0};
        const size_t unit = detail::work_unit(r, participants);
        std::vector<double> published(participants);
        const std::function<void(size_t)> body = [&](size_t participant) {
            double best = std::numeric_limits<double>::infinity();
            while (true) {
                const size_t begin = next.fetch_add(unit);
                if (begin >= r) {
                    break;
                }
                scan(begin, std::min(begin + unit, r), best);
            }
            published[participant] = best;
        };
        const double teamed = median_us([&] {
            next.store(0);
            team.Run(body);
            checksum += *std::min_element(published.begin(), published.end());
        });
        std::printf("%8zu %12.1f %12.1f   checksum %.3f\n", r, serial, teamed, checksum);
    }
}
