/**
 * @file test_agglomerative_row_cache.cpp
 * @brief The row-cache kernel against the 5.20.0 heap, compared bitwise.
 */

#include <gtest/gtest.h>

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <limits>
#include <optional>
#include <random>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#include "oecluster/StorageBackend.h"
#include "oecluster/clustering/Agglomerative.h"

#include "../../src/clustering/AgglomerativeRowCache.h"
#include "../../src/clustering/StepTeam.h"
#include "agglomerative_oracle.h"

using namespace OECluster;

namespace {

const double INF = std::numeric_limits<double>::infinity();
const double NaN = std::numeric_limits<double>::quiet_NaN();
const size_t ALL_TEAM = 0;
const size_t ALL_SERIAL = std::numeric_limits<size_t>::max();
const size_t OVERFLOW_CHUNK = std::numeric_limits<size_t>::max();

size_t gTotal = 0;
size_t gMatched = 0;
std::vector<std::string> gMismatches;

bool SameBits(double lhs, double rhs) {
    return std::memcmp(&lhs, &rhs, sizeof(double)) == 0;
}

std::vector<double> Condensed(size_t n,
                              const std::function<double(size_t, size_t)>& distance) {
    std::vector<double> condensed;
    condensed.reserve(n * (n - 1) / 2);
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i + 1; j < n; ++j) {
            condensed.push_back(distance(i, j));
        }
    }
    return condensed;
}

DenseStorage MakeStorage(size_t n, const std::vector<double>& condensed) {
    DenseStorage storage(n);
    size_t k = 0;
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i + 1; j < n; ++j) {
            storage.Set(i, j, condensed[k++]);
        }
    }
    return storage;
}

struct Fixture {
    std::string name;
    size_t n;
    std::vector<double> condensed;
};

std::vector<double> Random(size_t n, unsigned seed) {
    std::mt19937 generator(seed);
    std::uniform_real_distribution<double> value(0.0, 1.0);
    return Condensed(n, [&](size_t, size_t) { return value(generator); });
}

// Random values rounded to two decimals: no visible pattern, exact ties often.
std::vector<double> RandomTied(size_t n, unsigned seed) {
    std::mt19937 generator(seed);
    std::uniform_int_distribution<int> value(0, 99);
    return Condensed(n, [&](size_t, size_t) {
        return static_cast<double>(value(generator)) / 100.0;
    });
}

std::vector<double> Scrambled(size_t n) {
    return Condensed(n, [](size_t i, size_t j) {
        return static_cast<double>((i * 7 + j * 13) % 6 + 1);
    });
}

std::vector<double> Quantized(size_t n, unsigned seed, int levels) {
    std::mt19937 generator(seed);
    std::uniform_int_distribution<int> level(0, levels);
    return Condensed(n, [&](size_t, size_t) {
        return static_cast<double>(level(generator)) / levels;
    });
}

// Equal-sized blocks: one distance inside a block, another between any two.
// Every within-block pair ties with every other, and so does every across.
std::vector<double> Blocked(size_t n, size_t blocks) {
    const size_t width = (n + blocks - 1) / blocks;
    return Condensed(n, [&](size_t i, size_t j) {
        return i / width == j / width ? 0.25 : 0.75;
    });
}

std::vector<double> AllEqual(size_t n, double value) {
    return Condensed(n, [&](size_t, size_t) { return value; });
}

// A library of copies: items inside a copy group are identical, and the groups
// sit at quantized distances from each other.
std::vector<double> DuplicateHeavy(size_t n, size_t group, unsigned seed) {
    std::mt19937 generator(seed);
    std::uniform_int_distribution<int> level(1, 9);
    const size_t groups = (n + group - 1) / group;
    std::vector<double> between(groups * groups, 0.0);
    for (size_t a = 0; a < groups; ++a) {
        for (size_t b = a + 1; b < groups; ++b) {
            const double value = static_cast<double>(level(generator)) / 10.0;
            between[a * groups + b] = value;
            between[b * groups + a] = value;
        }
    }
    return Condensed(n, [&](size_t i, size_t j) {
        const size_t a = i / group;
        const size_t b = j / group;
        return a == b ? 0.0 : between[a * groups + b];
    });
}

std::vector<Fixture> Fixtures() {
    std::vector<Fixture> fixtures;
    for (size_t n : {size_t{2}, size_t{3}, size_t{4}, size_t{5}, size_t{6}, size_t{7},
                     size_t{8}, size_t{9}, size_t{10}, size_t{11}, size_t{12}}) {
        fixtures.push_back({"scrambled" + std::to_string(n), n, Scrambled(n)});
        fixtures.push_back({"random" + std::to_string(n), n, Random(n, 7u)});
    }
    for (size_t n : {size_t{17}, size_t{40}, size_t{97}, size_t{200}, size_t{333}}) {
        fixtures.push_back({"random" + std::to_string(n), n, Random(n, 11u)});
        fixtures.push_back({"randomtied" + std::to_string(n), n, RandomTied(n, 13u)});
        fixtures.push_back({"scrambled" + std::to_string(n), n, Scrambled(n)});
        fixtures.push_back({"quantized3_" + std::to_string(n), n, Quantized(n, 5u, 3)});
        fixtures.push_back({"blocked" + std::to_string(n), n, Blocked(n, 5)});
        fixtures.push_back({"allequal" + std::to_string(n), n, AllEqual(n, 0.5)});
        fixtures.push_back({"duplicates" + std::to_string(n), n,
                            DuplicateHeavy(n, 4, 17u)});
    }
    // Derived positive overflow: average and weighted reach +inf from a finite
    // input, which 5.20.0 orders normally and the new path must reproduce.
    for (size_t n : {size_t{3}, size_t{5}, size_t{10}, size_t{64}}) {
        fixtures.push_back({"dblmax" + std::to_string(n), n,
                            AllEqual(n, std::numeric_limits<double>::max())});
    }
    fixtures.push_back({"zeros8", 8, AllEqual(8, 0.0)});
    // The one input where the operand order of a complete-linkage update is
    // observable: std::max(-0.0, +0.0) keeps its first argument, so swapping
    // the children changes the sign bit of the second merge height.
    fixtures.push_back({"signedzero3", 3, {0.0, -0.0, 0.0}});
    return fixtures;
}

double CutThreshold(const Fixture& fixture) {
    std::vector<double> sorted = fixture.condensed;
    std::sort(sorted.begin(), sorted.end());
    return sorted[sorted.size() * 2 / 5];
}

struct Config {
    size_t num_threads = 0;
    size_t chunk_size = 4096;
    std::optional<size_t> serial_cutoff;
};

std::string Describe(const Fixture& fixture, const AgglomerativeOptions& options,
                     const Config& config) {
    const char* linkage = options.linkage == AgglomerativeLinkageMethod::Complete
                              ? "complete"
                              : (options.linkage == AgglomerativeLinkageMethod::Average
                                     ? "average"
                                     : "weighted");
    std::string text = fixture.name;
    text += " ";
    text += linkage;
    text += options.distance_threshold >= 0.0 ? " threshold" : " n_clusters";
    text += options.compute_full_tree ? " full" : " early";
    text += " threads=" + std::to_string(config.num_threads);
    text += " chunk=" + std::to_string(config.chunk_size);
    text += " cutoff=" +
            (config.serial_cutoff ? std::to_string(*config.serial_cutoff)
                                  : std::string("default"));
    return text;
}

// Bitwise on the heights, exactly equal on everything else. A tolerance would
// pass the float-reassociation class the design was chosen to avoid, so the
// heights are compared as bit patterns.
bool SameResult(const AgglomerativeResult& expected, const AgglomerativeResult& actual,
                std::string& diagnosis) {
    if (expected.Labels() != actual.Labels()) {
        diagnosis = "labels";
        return false;
    }
    if (expected.Members() != actual.Members()) {
        diagnosis = "members";
        return false;
    }
    if (expected.ChildrenLeft() != actual.ChildrenLeft()) {
        diagnosis = "children_left";
        return false;
    }
    if (expected.ChildrenRight() != actual.ChildrenRight()) {
        diagnosis = "children_right";
        return false;
    }
    if (expected.ClusterSizes() != actual.ClusterSizes()) {
        diagnosis = "cluster_sizes";
        return false;
    }
    if (expected.Distances().size() != actual.Distances().size()) {
        diagnosis = "height count";
        return false;
    }
    for (size_t k = 0; k < expected.Distances().size(); ++k) {
        if (!SameBits(expected.Distances()[k], actual.Distances()[k])) {
            diagnosis = "height " + std::to_string(k);
            return false;
        }
    }
    return true;
}

void Compare(const Fixture& fixture, const AgglomerativeResult& expected,
             const AgglomerativeOptions& options, const Config& config) {
    AgglomerativeOptions actual_options = options;
    actual_options.num_threads = config.num_threads;
    actual_options.chunk_size = config.chunk_size;
    const DenseStorage storage = MakeStorage(fixture.n, fixture.condensed);
    const AgglomerativeResult actual = detail::row_cache_result(
        storage, actual_options, config.serial_cutoff);
    ++gTotal;
    std::string diagnosis;
    if (SameResult(expected, actual, diagnosis)) {
        ++gMatched;
        return;
    }
    gMismatches.push_back(Describe(fixture, options, config) + ": " + diagnosis);
}

std::vector<AgglomerativeOptions> CutModes(const Fixture& fixture) {
    std::vector<AgglomerativeOptions> modes;
    for (const AgglomerativeLinkageMethod linkage :
         {AgglomerativeLinkageMethod::Complete, AgglomerativeLinkageMethod::Average,
          AgglomerativeLinkageMethod::Weighted}) {
        for (const bool full_tree : {true, false}) {
            for (const size_t clusters :
                 {size_t{1}, size_t{2}, std::max<size_t>(1, fixture.n / 3)}) {
                AgglomerativeOptions options;
                options.linkage = linkage;
                options.compute_full_tree = full_tree;
                options.n_clusters = std::min(clusters, fixture.n);
                options.distance_threshold = -1.0;
                modes.push_back(options);
            }
            AgglomerativeOptions options;
            options.linkage = linkage;
            options.compute_full_tree = full_tree;
            options.distance_threshold = CutThreshold(fixture);
            modes.push_back(options);
        }
    }
    return modes;
}

}  // namespace

// The centrepiece: every fixture, linkage, cut mode and compute_full_tree,
// against the shipped 5.20.0 heap at the library's default options.
TEST(AgglomerativeRowCacheTest, MatchesTheHeapBitwise) {
    const Config config;
    for (const Fixture& fixture : Fixtures()) {
        const DenseStorage storage = MakeStorage(fixture.n, fixture.condensed);
        for (const AgglomerativeOptions& options : CutModes(fixture)) {
            const AgglomerativeResult expected =
                agglomerative_oracle::heap_cluster(storage, options);
            Compare(fixture, expected, options, config);
        }
    }
    for (const std::string& mismatch : gMismatches) {
        ADD_FAILURE() << mismatch;
    }
    std::printf("[default config] matched %zu of %zu\n", gMatched, gTotal);
}

// Thread counts, chunk sizes and both sides of the serial cutoff, on the
// fixtures whose ties the key has to break.
TEST(AgglomerativeRowCacheTest, MatchesTheHeapUnderEveryConfiguration) {
    const size_t before = gTotal;
    const size_t matched_before = gMatched;
    const std::vector<std::string> interesting = {
        "scrambled7",  "scrambled12", "randomtied40",  "quantized3_40",
        "blocked40",   "allequal40",  "duplicates40",  "random97",
        "blocked97",   "allequal97",  "duplicates200", "dblmax10",
        "dblmax64",    "zeros8",      "signedzero3"};
    for (const Fixture& fixture : Fixtures()) {
        if (std::find(interesting.begin(), interesting.end(), fixture.name) ==
            interesting.end()) {
            continue;
        }
        const DenseStorage storage = MakeStorage(fixture.n, fixture.condensed);
        for (const AgglomerativeOptions& options : CutModes(fixture)) {
            const AgglomerativeResult expected =
                agglomerative_oracle::heap_cluster(storage, options);
            for (const size_t threads : {size_t{1}, size_t{2}, size_t{8}, size_t{0}}) {
                for (const size_t chunk :
                     {size_t{1}, size_t{7}, fixture.n, size_t{4096}}) {
                    for (const std::optional<size_t> cutoff :
                         {std::optional<size_t>(ALL_TEAM),
                          std::optional<size_t>(ALL_SERIAL),
                          std::optional<size_t>()}) {
                        Compare(fixture, expected, options,
                                Config{threads, chunk, cutoff});
                    }
                }
            }
        }
    }
    for (size_t k = 0; k < gMismatches.size(); ++k) {
        ADD_FAILURE() << gMismatches[k];
    }
    std::printf("[config matrix] matched %zu of %zu\n", gMatched - matched_before,
                gTotal - before);
}

// 5.20.0 silently skipped the copy at this chunk size and returned a tree built
// from an all-infinity table, so the overflow chunk is pinned against the new
// path's own defined result rather than against the heap's.
TEST(AgglomerativeRowCacheTest, AnOverflowingChunkDoesNotSkipTheCopy) {
    size_t cases = 0;
    size_t matched = 0;
    size_t heap_differed = 0;
    for (const Fixture& fixture : Fixtures()) {
        const DenseStorage storage = MakeStorage(fixture.n, fixture.condensed);
        for (const AgglomerativeOptions& options : CutModes(fixture)) {
            AgglomerativeOptions sane = options;
            sane.chunk_size = 4096;
            const AgglomerativeResult reference =
                detail::row_cache_result(storage, sane);
            AgglomerativeOptions overflowing = options;
            overflowing.chunk_size = OVERFLOW_CHUNK;
            const AgglomerativeResult actual =
                detail::row_cache_result(storage, overflowing);
            std::string diagnosis;
            ++cases;
            if (SameResult(reference, actual, diagnosis)) {
                ++matched;
            } else {
                ADD_FAILURE() << fixture.name << " " << diagnosis;
            }
            const AgglomerativeResult heap =
                agglomerative_oracle::heap_cluster(storage, overflowing);
            if (!SameResult(reference, heap, diagnosis)) {
                ++heap_differed;
            }
        }
    }
    std::printf("[overflow chunk] matched %zu of %zu; 5.20.0 differed in %zu\n",
                matched, cases, heap_differed);
}

TEST(AgglomerativeRowCacheTest, RefusesANonFiniteInputValue) {
    for (const double bad : {INF, -INF, NaN}) {
        for (const AgglomerativeLinkageMethod linkage :
             {AgglomerativeLinkageMethod::Complete, AgglomerativeLinkageMethod::Average,
              AgglomerativeLinkageMethod::Weighted}) {
            std::vector<double> condensed = Random(6, 3u);
            condensed[7] = bad;  // the pair (1, 4)
            const DenseStorage storage = MakeStorage(6, condensed);
            AgglomerativeOptions options;
            options.linkage = linkage;
            options.n_clusters = 2;
            try {
                detail::row_cache_result(storage, options);
                FAIL() << "accepted a non-finite input distance";
            } catch (const std::runtime_error& error) {
                EXPECT_NE(std::string(error.what())
                              .find("read a non-finite distance between items 1 and 4"),
                          std::string::npos)
                    << error.what();
            }
        }
    }
}

// The node-id rule, isolated: slots are recycled and node ids are not, so a
// comparison on slot indices reorders a tie the moment a merged cluster takes
// a retired slot. The fixture below reaches that state, and the heap's own
// answer is the thing being matched.
TEST(AgglomerativeRowCacheTest, BreaksTiesOnNodeIdsNotSlots) {
    for (size_t n = 4; n <= 40; ++n) {
        const std::vector<double> condensed = Blocked(n, 3);
        const DenseStorage storage = MakeStorage(n, condensed);
        for (const AgglomerativeLinkageMethod linkage :
             {AgglomerativeLinkageMethod::Complete, AgglomerativeLinkageMethod::Average,
              AgglomerativeLinkageMethod::Weighted}) {
            AgglomerativeOptions options;
            options.linkage = linkage;
            options.n_clusters = 1;
            const AgglomerativeResult expected =
                agglomerative_oracle::heap_cluster(storage, options);
            const AgglomerativeResult actual =
                detail::row_cache_result(storage, options);
            std::string diagnosis;
            ASSERT_TRUE(SameResult(expected, actual, diagnosis))
                << "n=" << n << " " << diagnosis;
        }
    }
}

// Not a check: the measurement that sets ROW_CACHE_SERIAL_CUTOFF. Run it with
//   proto_tests --gtest_also_run_disabled_tests
//       --gtest_filter=AgglomerativeRowCacheTest.DISABLED_MeasureMergeStepSerialCutoff
// only while the 1-minute load average is below half the core count. Each row
// times one serial merge step and one team step over R active slots, the
// median of 5, reading a 512 MB buffer through the condensed index of a
// (2^17 + 1)-slot matrix so the reads keep a real workspace's index arithmetic
// and cache behavior. Every row gets the two reads and the write of the
// linkage update, and RESCANS_PER_STEP of them additionally rescan the whole
// row; that rate is not a guess but the measured one, between 2.5 rescans per
// step at N = 50 and 5.8 at N = 10000 on random distances, and it is what the
// step actually costs -- the rescans read about 4.6 times as many doubles as
// the updates do, so modelling the update alone puts the crossover several
// times too high. The constant is the smallest R from which the team is faster
// at every larger R in the grid, or 2^18 if none.
TEST(AgglomerativeRowCacheTest, DISABLED_MeasureMergeStepSerialCutoff) {
    const size_t RESCANS_PER_STEP = 5;
    const size_t virtual_n = (size_t{1} << 17) + 1;
    const size_t buffer_size = size_t{1} << 26;
    std::vector<double> buffer(buffer_size);
    for (size_t k = 0; k < buffer_size; ++k) {
        buffer[k] = static_cast<double>((k * 2654435761u) % 1000) / 1000.0;
    }
    std::vector<size_t> base(virtual_n);
    for (size_t i = 0; i < virtual_n; ++i) {
        base[i] = virtual_n * i - i * (i + 1) / 2;
    }
    auto pair_index = [&](size_t i, size_t j) {
        const size_t low = i < j ? i : j;
        const size_t high = i < j ? j : i;
        return (base[low] + high - low - 1) % buffer_size;
    };
    const size_t participants =
        detail::row_cache_participants(0, virtual_n, std::thread::hardware_concurrency());
    detail::StepTeam team(participants);
    const size_t slot_low = virtual_n / 3;
    const size_t slot_high = virtual_n / 2;
    const size_t kept = slot_low;
    std::printf("participants=%zu\n%8s %12s %12s\n", participants, "R", "serial_us",
                "team_us");
    for (size_t r = size_t{1} << 6; r <= (size_t{1} << 17); r <<= 1) {
        std::vector<size_t> live(r);
        const size_t stride = (virtual_n - 1) / r;
        for (size_t k = 0; k < r; ++k) {
            live[k] = 1 + k * stride;
        }
        // One row in every `rescan_period` rescans, so a whole step rescans
        // RESCANS_PER_STEP rows whatever R is -- the shape the counters show.
        const size_t rescan_period = std::max<size_t>(1, r / RESCANS_PER_STEP);
        auto run = [&](size_t begin, size_t end, double& best) {
            for (size_t p = begin; p < end; ++p) {
                const size_t x = live[p];
                const double distance_low = buffer[pair_index(x, slot_low)];
                const double distance_high = buffer[pair_index(x, slot_high)];
                const double merged = ((3.0 * distance_low) + (5.0 * distance_high)) / 8.0;
                buffer[pair_index(x, kept)] = merged;
                best = std::min(best, merged);
                if (p % rescan_period == 0) {
                    for (size_t q = 0; q < r; ++q) {
                        best = std::min(best, buffer[pair_index(x, live[q])]);
                    }
                }
            }
        };
        auto median_us = [&](auto&& step) {
            std::vector<double> times;
            for (int rep = 0; rep < 5; ++rep) {
                const auto start = std::chrono::steady_clock::now();
                step();
                times.push_back(std::chrono::duration<double, std::micro>(
                                    std::chrono::steady_clock::now() - start)
                                    .count());
            }
            std::sort(times.begin(), times.end());
            return times[2];
        };
        double checksum = 0.0;
        const double serial = median_us([&] {
            double best = std::numeric_limits<double>::infinity();
            run(0, r, best);
            checksum += best;
        });
        std::atomic<size_t> next{0};
        const size_t unit = detail::work_unit(r, participants);
        std::vector<double> published(participants,
                                      std::numeric_limits<double>::infinity());
        const double parallel = median_us([&] {
            next.store(0, std::memory_order_relaxed);
            team.Run([&](size_t participant) {
                double best = std::numeric_limits<double>::infinity();
                while (true) {
                    const size_t begin = next.fetch_add(unit, std::memory_order_relaxed);
                    if (begin >= r) {
                        break;
                    }
                    run(begin, std::min(begin + unit, r), best);
                }
                published[participant] = best;
            });
            for (const double value : published) {
                // A participant that claimed no unit publishes the identity.
                if (std::isfinite(value)) {
                    checksum += value;
                }
            }
        });
        std::printf("%8zu %12.2f %12.2f\n", r, serial, parallel);
        EXPECT_TRUE(std::isfinite(checksum));
    }
}

TEST(AgglomerativeRowCacheTest, TheDefaultParticipantCountIsCappedAtEight) {
    EXPECT_EQ(detail::row_cache_participants(0, 100, 14), 8u);
    EXPECT_EQ(detail::row_cache_participants(0, 100, 4), 4u);
    EXPECT_EQ(detail::row_cache_participants(0, 100, 0), 1u);
    EXPECT_EQ(detail::row_cache_participants(0, 3, 14), 3u);
    EXPECT_EQ(detail::row_cache_participants(14, 100, 4), 14u);
    EXPECT_EQ(detail::row_cache_participants(14, 10, 4), 10u);
}

TEST(AgglomerativeRowCacheTest, ATeamIsBuiltOnlyWhenAStepCanReachTheCutoff) {
    // One participant never needs a team, whatever the cutoff.
    EXPECT_FALSE(detail::row_cache_team_is_useful(1, 100000, 0));

    // The largest step updates n - 2 rows, so n - 2 is the boundary.
    EXPECT_FALSE(detail::row_cache_team_is_useful(8, 513, 512));
    EXPECT_TRUE(detail::row_cache_team_is_useful(8, 514, 512));

    // A forced cutoff of 0 always wants the team; that is how the differential
    // exercises the team path on small fixtures.
    EXPECT_TRUE(detail::row_cache_team_is_useful(8, 2, 0));

    // n below 2 must not wrap through the unsigned subtraction.
    EXPECT_FALSE(detail::row_cache_team_is_useful(8, 0, 0));
    EXPECT_FALSE(detail::row_cache_team_is_useful(8, 1, 0));
}
