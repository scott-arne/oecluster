/**
 * @file test_maxmin_select.cpp
 * @brief Farthest-first (MaxMin) selection over distance matrices and
 * comparisons.
 */

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "oecluster/Error.h"
#include "oecluster/GateFacts.h"
#include "oecluster/StorageBackend.h"
#include "oecluster/clustering/DiversitySelection.h"

#include "diversity_test_support.h"

using namespace OECluster;
using namespace diversity_test;

namespace {

const double NaN = std::numeric_limits<double>::quiet_NaN();

MaxMinOptions CountOptions(size_t count, size_t seed = 0) {
    MaxMinOptions options;
    options.count = count;
    options.seed = seed;
    return options;
}

// Brute force from the definition: each pick maximizes the minimum distance
// to everything already selected, ties to the smaller index. It shares no
// code with the kernel, so agreement is evidence rather than tautology.
MaxMinSelection OracleSelect(size_t n, const std::vector<double>& condensed,
                             size_t count, double threshold,
                             std::vector<size_t> initial) {
    const auto distance = [&](size_t a, size_t b) {
        const size_t i = std::min(a, b);
        const size_t j = std::max(a, b);
        return condensed[n * i - i * (i + 1) / 2 + j - i - 1];
    };
    MaxMinSelection selection;
    selection.indices = initial;
    selection.pick_distances.assign(initial.size(), NaN);
    while (true) {
        if (count != 0 && selection.indices.size() == count) {
            selection.stop = MaxMinStop::Count;
            return selection;
        }
        size_t best = n;
        double best_distance = 0.0;
        for (size_t candidate = 0; candidate < n; ++candidate) {
            bool taken = false;
            double nearest = std::numeric_limits<double>::infinity();
            for (const size_t member : selection.indices) {
                taken = taken || member == candidate;
                if (member != candidate) {
                    nearest = std::min(nearest, distance(member, candidate));
                }
            }
            if (!taken && (best == n || nearest > best_distance)) {
                best = candidate;
                best_distance = nearest;
            }
        }
        if (best == n) {
            selection.stop = MaxMinStop::Exhausted;
            return selection;
        }
        if (!std::isnan(threshold) && best_distance <= threshold) {
            selection.stop = MaxMinStop::Threshold;
            return selection;
        }
        selection.indices.push_back(best);
        selection.pick_distances.push_back(best_distance);
    }
}

void ExpectSameSelection(const MaxMinSelection& actual,
                         const MaxMinSelection& expected) {
    EXPECT_EQ(actual.indices, expected.indices);
    EXPECT_EQ(actual.stop, expected.stop);
    ASSERT_EQ(actual.pick_distances.size(), expected.pick_distances.size());
    for (size_t i = 0; i < expected.pick_distances.size(); ++i) {
        if (std::isnan(expected.pick_distances[i])) {
            EXPECT_TRUE(std::isnan(actual.pick_distances[i])) << "position " << i;
        } else {
            // Bit-identical: both sides pick the same stored double.
            EXPECT_EQ(actual.pick_distances[i], expected.pick_distances[i])
                << "position " << i;
        }
    }
}

void ExpectInvalidArgument(const std::function<void()>& call,
                           const std::string& message) {
    try {
        call();
        FAIL() << "expected std::invalid_argument: " << message;
    } catch (const std::invalid_argument& error) {
        EXPECT_EQ(std::string(error.what()), message);
    }
}

}  // namespace

TEST(MaxMinSelectTest, MatchesABruteForceOracle) {
    for (const size_t n : {size_t{1}, size_t{2}, size_t{5}, size_t{9}}) {
        const std::vector<double> condensed = Scrambled(n);
        const DenseStorage storage = MakeStorage(n, condensed);
        for (size_t seed = 0; seed < n; ++seed) {
            for (size_t count = 1; count <= n; ++count) {
                SCOPED_TRACE("n " + std::to_string(n) + ", seed " +
                             std::to_string(seed) + ", count " +
                             std::to_string(count));
                ExpectSameSelection(
                    maxmin_select(storage, CountOptions(count, seed)),
                    OracleSelect(n, condensed, count, NaN, {seed}));
            }
            for (const double threshold : {0.0, 1.0, 2.5, 6.0}) {
                MaxMinOptions options;
                options.threshold = threshold;
                options.seed = seed;
                SCOPED_TRACE("n " + std::to_string(n) + ", seed " +
                             std::to_string(seed) + ", threshold " +
                             std::to_string(threshold));
                ExpectSameSelection(
                    maxmin_select(storage, options),
                    OracleSelect(n, condensed, 0, threshold, {seed}));
            }
        }
    }
}

TEST(MaxMinSelectTest, PickDistancesNeverIncreaseAfterTheSeed) {
    const size_t n = 17;
    const DenseStorage storage = MakeStorage(n, Scrambled(n));

    const MaxMinSelection selection = maxmin_select(storage, CountOptions(n));

    EXPECT_TRUE(std::isnan(selection.pick_distances[0]));
    for (size_t i = 2; i < selection.pick_distances.size(); ++i) {
        EXPECT_LE(selection.pick_distances[i], selection.pick_distances[i - 1]);
    }
}

TEST(MaxMinSelectTest, TheCountStopWinsWhenEveryItemIsSelected) {
    const DenseStorage storage = MakeStorage(4, Line(4));

    const MaxMinSelection selection = maxmin_select(storage, CountOptions(4));

    EXPECT_EQ(selection.indices, std::vector<size_t>({0, 3, 1, 2}));
    EXPECT_EQ(selection.stop, MaxMinStop::Count);
}

TEST(MaxMinSelectTest, AThresholdAloneRunsToExhaustionOnDistinctItems) {
    const DenseStorage storage = MakeStorage(4, Line(4));
    MaxMinOptions options;
    options.threshold = 0.0;

    const MaxMinSelection selection = maxmin_select(storage, options);

    EXPECT_EQ(selection.indices, std::vector<size_t>({0, 3, 1, 2}));
    EXPECT_EQ(selection.stop, MaxMinStop::Exhausted);
}

TEST(MaxMinSelectTest, TheFirstOfCountAndThresholdToTriggerWins) {
    const DenseStorage storage = MakeStorage(6, Line(6));
    MaxMinOptions options;
    options.count = 5;
    options.threshold = 2.0;

    const MaxMinSelection selection = maxmin_select(storage, options);

    EXPECT_EQ(selection.indices, std::vector<size_t>({0, 5}));
    EXPECT_EQ(selection.stop, MaxMinStop::Threshold);
}

TEST(MaxMinSelectTest, ExtendsAnInitialSelection) {
    const DenseStorage storage = MakeStorage(6, Line(6));
    MaxMinOptions options;
    options.count = 4;
    options.initial = {2, 3};

    const MaxMinSelection selection = maxmin_select(storage, options);

    ExpectSameSelection(selection, OracleSelect(6, Line(6), 4, NaN, {2, 3}));
    EXPECT_EQ(selection.indices, std::vector<size_t>({2, 3, 0, 5}));
}

TEST(MaxMinSelectTest, AnInitialSelectionThatMeetsTheCountIsReturnedAsIs) {
    const DenseStorage storage = MakeStorage(6, Line(6));
    MaxMinOptions options;
    options.count = 2;
    options.initial = {4, 1};

    const MaxMinSelection selection = maxmin_select(storage, options);

    EXPECT_EQ(selection.indices, std::vector<size_t>({4, 1}));
    EXPECT_EQ(selection.stop, MaxMinStop::Count);
}

TEST(MaxMinSelectTest, TheFarthestSeedIsFarthestFromItemZero) {
    // Items 3 and 4 tie at 9 from item 0; the smaller index wins.
    const DenseStorage storage = MakeStorage(5, Positions({0, 4, 1, 9, 9}));
    MaxMinOptions options = CountOptions(1);
    options.seed_mode = MaxMinSeed::Farthest;

    EXPECT_EQ(maxmin_select(storage, options).indices, std::vector<size_t>({3}));
}

TEST(MaxMinSelectTest, TheFarthestSeedOfOneItemIsItemZero) {
    const DenseStorage storage = MakeStorage(1, {});
    MaxMinOptions options = CountOptions(1);
    options.seed_mode = MaxMinSeed::Farthest;

    EXPECT_EQ(maxmin_select(storage, options).indices, std::vector<size_t>({0}));
}

TEST(MaxMinSelectTest, TheMedoidSeedHasTheSmallestDistanceSum) {
    // Sums: 24, 21, 20, 28, 31.
    const DenseStorage storage = MakeStorage(5, Positions({0, 1, 2, 10, 11}));
    MaxMinOptions options = CountOptions(1);
    options.seed_mode = MaxMinSeed::Medoid;

    EXPECT_EQ(maxmin_select(storage, options).indices, std::vector<size_t>({2}));
}

TEST(MaxMinSelectTest, TheMedoidSeedIgnoresThreadCountAndChunkSize) {
    const size_t n = 23;
    const DenseStorage storage = MakeStorage(n, Scrambled(n));
    MaxMinOptions options = CountOptions(6);
    options.seed_mode = MaxMinSeed::Medoid;
    options.num_threads = 1;
    const MaxMinSelection expected = maxmin_select(storage, options);

    for (const size_t threads : {size_t{2}, size_t{8}, size_t{1000}}) {
        for (const size_t chunk : {size_t{1}, size_t{2},
                                   std::numeric_limits<size_t>::max()}) {
            options.num_threads = threads;
            options.chunk_size = chunk;
            ExpectSameSelection(maxmin_select(storage, options), expected);
        }
    }
}

TEST(MaxMinSelectTest, RefusesANonFiniteDistanceItReads) {
    std::vector<double> condensed = Line(6);
    condensed[2] = NaN;  // (0, 3)
    const DenseStorage storage = MakeStorage(6, condensed);

    ExpectInvalidArgument(
        [&] { maxmin_select(storage, CountOptions(3)); },
        "Diversity selection read a non-finite distance between items 0 and 3");
}

// A non-finite entry the selection itself never reads is still refused by the
// Medoid seed, which sums every row; the first pair in condensed order is
// named whatever the thread count.
TEST(MaxMinSelectTest, TheMedoidSeedRefusesAnyNonFiniteEntry) {
    std::vector<double> condensed = Line(6);
    condensed[9] = std::numeric_limits<double>::infinity();  // (2, 3)
    condensed[14] = NaN;                                      // (4, 5)
    const DenseStorage storage = MakeStorage(6, condensed);
    MaxMinOptions options = CountOptions(1);
    options.seed_mode = MaxMinSeed::Medoid;
    options.num_threads = 4;
    options.chunk_size = 1;

    ExpectInvalidArgument(
        [&] { maxmin_select(storage, options); },
        "Diversity selection read a non-finite distance between items 2 and 3");
}

TEST(MaxMinSelectTest, TheMedoidSeedRefusesADistanceSumThatOverflows) {
    const DenseStorage storage = MakeStorage(3, {1e308, 1e308, 1e308});
    MaxMinOptions options = CountOptions(1);
    options.seed_mode = MaxMinSeed::Medoid;

    ExpectInvalidArgument(
        [&] { maxmin_select(storage, options); },
        "MaxMin selection medoid seed: the distance sum of item 0 overflows");
}

TEST(MaxMinSelectValidationTest, RefusesEachInvalidRequest) {
    const DenseStorage storage = MakeStorage(6, Line(6));
    const auto expect = [&](MaxMinOptions options, const std::string& message) {
        SCOPED_TRACE(message);
        ExpectInvalidArgument([&] { maxmin_select(storage, options); }, message);
    };

    MaxMinOptions options = CountOptions(2);
    options.seed_mode = static_cast<MaxMinSeed>(7);
    expect(options, "Unknown MaxMin selection seed mode");

    options = CountOptions(2);
    options.chunk_size = 0;
    expect(options, "MaxMin selection chunk_size must be at least one");

    expect(MaxMinOptions(),
           "MaxMin selection requires a count, a threshold, or both");

    expect(CountOptions(7),
           "MaxMin selection count must be at most the item count (6)");

    for (const double threshold :
         {std::numeric_limits<double>::infinity(), -1.0}) {
        options = MaxMinOptions();
        options.threshold = threshold;
        expect(options,
               "MaxMin selection threshold must be finite and non-negative");
    }

    expect(CountOptions(2, 6), "MaxMin selection seed is outside the item range");

    options = CountOptions(3, 1);
    options.initial = {0};
    expect(options, "MaxMin selection initial cannot be combined with a seed");

    options = CountOptions(3);
    options.seed_mode = MaxMinSeed::Farthest;
    options.initial = {0};
    expect(options, "MaxMin selection initial cannot be combined with a seed");

    options = CountOptions(1);
    options.initial = {0, 1};
    expect(options, "MaxMin selection initial holds more entries than count");

    options = CountOptions(3);
    options.initial = {0, 6};
    expect(options, "MaxMin selection initial index is outside the item range");

    options = CountOptions(3);
    options.initial = {2, 2};
    expect(options, "MaxMin selection initial entries must be unique");
}

TEST(MaxMinSelectValidationTest, RefusesStorageItCannotRead) {
    ExpectInvalidArgument(
        [] { maxmin_select(SparseStorage(4, 0.5), CountOptions(2)); },
        "MaxMin selection requires complete pairwise distances; SparseStorage "
        "is not supported");
    ExpectInvalidArgument(
        [] { maxmin_select(NullDataStorage(4), CountOptions(2)); },
        "MaxMin selection requires contiguous dense or memory-mapped storage");
    ExpectInvalidArgument(
        [] { maxmin_select(DenseStorage(0), CountOptions(1)); },
        "MaxMin selection requires at least one item");
}

// The enum is checked before anything else, so a caller with two mistakes
// hears about the one that does not depend on the data.
TEST(MaxMinSelectValidationTest, ChecksTheSeedModeFirst) {
    MaxMinOptions options = CountOptions(2);
    options.seed_mode = static_cast<MaxMinSeed>(7);

    ExpectInvalidArgument([&] { maxmin_select(SparseStorage(4, 0.5), options); },
                          "Unknown MaxMin selection seed mode");
}

TEST(MaxMinSelectTest, MaskingYieldsDistinctIndicesWhenEveryDistanceTies) {
    const DenseStorage storage = MakeStorage(5, std::vector<double>(10, 0.0));

    const MaxMinSelection selection = maxmin_select(storage, CountOptions(5, 2));

    EXPECT_EQ(selection.indices, std::vector<size_t>({2, 0, 1, 3, 4}));
    EXPECT_EQ(selection.stop, MaxMinStop::Count);
}

TEST(MaxMinSelectTest, AThresholdOnlySelectionExtendsAnInitialSelection) {
    const DenseStorage storage = MakeStorage(6, Line(6));
    MaxMinOptions options;
    options.threshold = 1.0;
    options.initial = {2, 3};

    const MaxMinSelection selection = maxmin_select(storage, options);

    // After {2, 3}: 0 and 5 at 2 are beyond the threshold; then every
    // remaining item sits at 1, which is not.
    EXPECT_EQ(selection.indices, std::vector<size_t>({2, 3, 0, 5}));
    EXPECT_EQ(selection.stop, MaxMinStop::Threshold);
    ExpectSameSelection(selection, OracleSelect(6, Line(6), 0, 1.0, {2, 3}));
}

namespace {

GateFacts FactsWith(Capability is_distance, Capability zero_self,
                    DataIntegrity data_integrity) {
    GateFacts facts;
    facts.is_distance = is_distance;
    facts.zero_self = zero_self;
    facts.data_integrity = data_integrity;
    return facts;
}

void ExpectComparisonError(const std::function<void()>& call,
                           const std::string& message) {
    try {
        call();
        FAIL() << "expected ComparisonError: " << message;
    } catch (const ComparisonError& error) {
        EXPECT_EQ(std::string(error.what()), message);
    }
}

}  // namespace

// TableComparison throws on a self-pair or a reversed pair, so a clean run is
// itself the proof that the lazy path read only (min, max) off-diagonal pairs.
TEST(MaxMinSelectComparisonTest, TheTableEnforcesItsReadingContract) {
    TableComparison table(3, Line(3));

    EXPECT_EQ(table.Compare(0, 2), 2.0);
    EXPECT_THROW(table.Compare(1, 1), std::logic_error);
    EXPECT_THROW(table.Compare(2, 0), std::logic_error);
}

TEST(MaxMinSelectComparisonTest, MatchesTheMatrixAtEveryThreadCountAndChunkSize) {
    for (const size_t n : {size_t{1}, size_t{2}, size_t{7}, size_t{12}}) {
        const std::vector<double> condensed = Scrambled(n);
        const DenseStorage storage = MakeStorage(n, condensed);
        std::vector<MaxMinOptions> requests;
        for (const size_t count : {size_t{1}, (n + 1) / 2, n}) {
            requests.push_back(CountOptions(count, n / 2));
        }
        MaxMinOptions threshold_request;
        threshold_request.threshold = 2.0;
        requests.push_back(threshold_request);
        MaxMinOptions farthest = CountOptions(n);
        farthest.seed_mode = MaxMinSeed::Farthest;
        requests.push_back(farthest);
        if (n >= 3) {
            MaxMinOptions extend = CountOptions(n);
            extend.initial = {n - 1, 1};
            requests.push_back(extend);
        }

        for (size_t r = 0; r < requests.size(); ++r) {
            const MaxMinSelection expected = maxmin_select(storage, requests[r]);
            for (const size_t threads : {size_t{1}, size_t{2}, size_t{8}}) {
                for (const size_t chunk : {size_t{1}, size_t{2}, size_t{256}}) {
                    SCOPED_TRACE("n " + std::to_string(n) + ", request " +
                                 std::to_string(r) + ", threads " +
                                 std::to_string(threads) + ", chunk " +
                                 std::to_string(chunk));
                    MaxMinOptions options = requests[r];
                    options.num_threads = threads;
                    options.chunk_size = chunk;
                    TableComparison table(n, condensed);
                    ExpectSameSelection(maxmin_select(table, options), expected);
                }
            }
        }
    }
}

TEST(MaxMinSelectComparisonTest, AcceptsAnExtremeChunkSize) {
    const size_t n = 9;
    const std::vector<double> condensed = Scrambled(n);
    const MaxMinSelection expected =
        maxmin_select(MakeStorage(n, condensed), CountOptions(n));

    MaxMinOptions options = CountOptions(n);
    options.num_threads = 4;
    options.chunk_size = std::numeric_limits<size_t>::max();
    TableComparison table(n, condensed);

    ExpectSameSelection(maxmin_select(table, options), expected);
}

// chunk_size 1 keeps every fold on the ThreadPool path; a larger chunk would
// take the single-chunk serial shortcut and never construct the workers the
// cap exists to bound. Uncapped, the pool's reserve for this count terminates
// the process, as KMedoidsThreadCapTest documents.
TEST(MaxMinSelectComparisonTest, CapsAnAbsurdThreadCount) {
    const size_t n = 9;
    const std::vector<double> condensed = Scrambled(n);
    const MaxMinSelection expected =
        maxmin_select(MakeStorage(n, condensed), CountOptions(n));

    MaxMinOptions options = CountOptions(n);
    options.num_threads = std::size_t{1} << 61;
    options.chunk_size = 1;
    TableComparison table(n, condensed);

    ExpectSameSelection(maxmin_select(table, options), expected);
}

// Clones circulate through a free list, so their number is bounded by the
// workers that ran at once, not by the rows folded.
TEST(MaxMinSelectComparisonTest, ClonesAtMostOncePerWorker) {
    const size_t n = 40;
    TableComparison table(n, Scrambled(n));
    MaxMinOptions options = CountOptions(n);
    options.num_threads = 3;
    options.chunk_size = 1;

    maxmin_select(table, options);

    EXPECT_GE(table.NumClones(), 1u);
    EXPECT_LE(table.NumClones(), 3u);
}

// Masking item 0 keeps the diagonal unread, so a comparison whose d(x, x) is
// not zero -- and does not claim it is -- still picks the matrix's seed.
TEST(MaxMinSelectComparisonTest, TheFarthestSeedNeverReadsTheDiagonal) {
    const std::vector<double> condensed = Positions({0, 4, 1, 9, 9});
    TableComparison table(5, condensed, GateFacts(), 100.0);
    MaxMinOptions options = CountOptions(3);
    options.seed_mode = MaxMinSeed::Farthest;

    ExpectSameSelection(maxmin_select(table, options),
                        maxmin_select(MakeStorage(5, condensed), options));
}

TEST(MaxMinSelectComparisonTest, RefusesANonFiniteDistanceItReads) {
    std::vector<double> condensed = Line(6);
    condensed[2] = NaN;  // (0, 3)
    for (const size_t threads : {size_t{1}, size_t{4}}) {
        MaxMinOptions options = CountOptions(3);
        options.num_threads = threads;
        options.chunk_size = 1;
        TableComparison table(6, condensed);

        ExpectInvalidArgument(
            [&] { maxmin_select(table, options); },
            "Diversity selection read a non-finite distance between items 0 "
            "and 3");
    }
}

TEST(MaxMinSelectComparisonTest, RefusesComparisonsItsFactsRuleOut) {
    const auto expect = [](GateFacts facts, const std::string& message) {
        SCOPED_TRACE(message);
        // count=2 would call Compare if the facts gate didn't refuse first.
        CountingComparison counter(4, facts);
        ExpectComparisonError([&] { maxmin_select(counter, CountOptions(2)); },
                              message);
        EXPECT_EQ(counter.Count(), 0u) << "Compare called before facts refusal";
    };

    expect(FactsWith(Capability::No, Capability::Unknown, DataIntegrity::Complete),
           "MaxMin selection requires distances, but the comparison reports "
           "similarities");
    expect(FactsWith(Capability::Yes, Capability::No, DataIntegrity::Complete),
           "MaxMin selection requires a zero self-distance, but the comparison "
           "reports that d(x, x) is not zero");
    expect(FactsWith(Capability::Yes, Capability::Yes, DataIntegrity::NaNPresent),
           "MaxMin selection cannot rank distances the comparison declares may "
           "be non-finite (missing='propagate')");
    expect(FactsWith(Capability::Yes, Capability::Yes,
                     DataIntegrity::SubsetScored),
           "MaxMin selection cannot rank distances scored on per-pair feature "
           "subsets (missing='ignore'); they are not mutually comparable");
}

TEST(MaxMinSelectComparisonTest, AcceptsFactsThatAreKnownGoodOrUnknown) {
    TableComparison known(4, Line(4),
                          FactsWith(Capability::Yes, Capability::Yes,
                                    DataIntegrity::Complete));
    TableComparison unknown(4, Line(4));

    EXPECT_EQ(maxmin_select(known, CountOptions(2)).indices,
              std::vector<size_t>({0, 3}));
    EXPECT_EQ(maxmin_select(unknown, CountOptions(2)).indices,
              std::vector<size_t>({0, 3}));
}

TEST(MaxMinSelectComparisonTest, RefusesTheMedoidSeed) {
    TableComparison table(4, Line(4));
    MaxMinOptions options = CountOptions(2);
    options.seed_mode = MaxMinSeed::Medoid;

    ExpectInvalidArgument(
        [&] { maxmin_select(table, options); },
        "MaxMin selection seed mode Medoid requires a distance matrix");

    // Option validation comes first, in the matrix overload's order.
    options.count = 5;
    ExpectInvalidArgument(
        [&] { maxmin_select(table, options); },
        "MaxMin selection count must be at most the item count (4)");
}

TEST(MaxMinSelectComparisonTest, SharesTheMatrixValidation) {
    TableComparison table(6, Line(6));
    ExpectInvalidArgument([&] { maxmin_select(table, MaxMinOptions()); },
                          "MaxMin selection requires a count, a threshold, or "
                          "both");

    TableComparison empty(0, {});
    ExpectInvalidArgument([&] { maxmin_select(empty, CountOptions(1)); },
                          "MaxMin selection requires at least one item");
}
