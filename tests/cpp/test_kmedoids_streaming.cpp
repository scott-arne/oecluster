/**
 * @file test_kmedoids_streaming.cpp
 * @brief The streaming distance source, the budget arithmetic and the entry point.
 */
#include <gtest/gtest.h>

#include <algorithm>
#include <cstddef>
#include <limits>
#include <mutex>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "oecluster/Error.h"
#include "oecluster/GateFacts.h"
#include "oecluster/StorageBackend.h"
#include "oecluster/clustering/KMedoids.h"

#include "../../src/clustering/KMedoidsBudget.h"
#include "../../src/clustering/KMedoidsDistanceSource.h"
#include "diversity_test_support.h"
#include "mst_test_support.h"
#include "streaming_test_support.h"

using namespace OECluster;
using namespace diversity_test;
using namespace mst_test;
using namespace streaming_test;

// diversity_test::Line(n) is the fixture here: items on a line with
// d(i, j) = j - i, condensed for pairs i < j.

TEST(KMedoidsBudgetTest, ASweepIsEveryOrderedPairExceptTheDiagonal) {
    EXPECT_EQ(detail::sweep_comparisons(1), 0u);
    EXPECT_EQ(detail::sweep_comparisons(2), 2u);
    EXPECT_EQ(detail::sweep_comparisons(100), 9900u);
}

// A verification candidate is rescored against every item, for every slot, so
// the pass is k sweeps and not one. Pricing it as one is what let a budgeted
// run overrun its limit by a factor of k.
TEST(KMedoidsBudgetTest, AVerificationPassCostsOneSweepPerMedoid) {
    EXPECT_EQ(detail::verification_comparisons(100, 5),
              5 * detail::sweep_comparisons(100));
    EXPECT_EQ(detail::verification_comparisons(100, 0), 0u);
    // Saturates for the reason the counts below do: a wrapped cost would read
    // as a pass cheap enough to start.
    const size_t widest = size_t{1} << 32;
    EXPECT_EQ(detail::verification_comparisons(widest, 2),
              std::numeric_limits<size_t>::max());
}

TEST(KMedoidsBudgetTest, UnavoidableWorkCountsInitAndVerificationAndTheCache) {
    const size_t n = 100;
    const size_t k = 5;
    const size_t sweep = detail::sweep_comparisons(n);
    // Build: k init sweeps -- the global-medoid pass that picks the first
    // seed, then k - 1 gain passes -- plus k verification sweeps and one
    // cache fill.
    EXPECT_EQ(detail::unavoidable_comparisons(n, k, KMedoidsInit::Build),
              k * sweep + k * sweep + n * k);
    // FarthestFirst: about one init sweep.
    EXPECT_EQ(detail::unavoidable_comparisons(n, k, KMedoidsInit::FarthestFirst),
              sweep + k * sweep + n * k);
    // Explicit: no init sweep at all, so the verification is the whole of it.
    EXPECT_EQ(detail::unavoidable_comparisons(n, k, KMedoidsInit::Explicit),
              k * sweep + n * k);
}

// The products below are the whole reason both budget functions saturate: a
// wrapped count reads as *cheap* work, so the ceiling it feeds would wave
// through exactly the runs it exists to refuse. Each case names what the
// unguarded arithmetic would have answered instead.
TEST(KMedoidsBudgetTest, ASweepSaturatesRatherThanWrapping) {
    constexpr size_t LIMIT = std::numeric_limits<size_t>::max();
    // The widest item count whose sweep still fits: n * (n - 1) is one 2^32
    // short of 2^64.
    const size_t widest = size_t{1} << 32;
    EXPECT_EQ(detail::sweep_comparisons(widest), widest * (widest - 1));
    // One more item and the product passes 2^64, wrapping to a mere 2^32.
    EXPECT_EQ(detail::sweep_comparisons(widest + 1), LIMIT);
}

TEST(KMedoidsBudgetTest, UnavoidableWorkSaturatesOnEitherProduct) {
    constexpr size_t LIMIT = std::numeric_limits<size_t>::max();
    const size_t widest = size_t{1} << 32;

    // Sweeps product: one init sweep plus one verification sweep is twice a
    // sweep that already fills most of the range. Unguarded it would wrap to
    // roughly 2^64 - 2^33 and read as a run worth attempting.
    EXPECT_EQ(detail::unavoidable_comparisons(widest, 1, KMedoidsInit::FarthestFirst),
              LIMIT);

    // Cache-fill product: Explicit takes no init sweep, so one verification
    // sweep is the whole sweeps total, and the n * k cache fill is what tips
    // it over. Unguarded the sum lands exactly on 2^64 and wraps to zero --
    // free work.
    EXPECT_EQ(detail::unavoidable_comparisons(widest, 1, KMedoidsInit::Explicit),
              LIMIT);

    // The sweep count is a sum before it is a product, so it is guarded in
    // its own right: Build's k init sweeps plus k verification sweeps are 2k,
    // which wraps from k = 2^63 up. At exactly that boundary the wrapped sum
    // is zero, and with a single item no later guard would notice -- the
    // answer would be a plain 2^63 rather than a refusal.
    const size_t half = size_t{1} << 63;
    EXPECT_EQ(detail::unavoidable_comparisons(1, half, KMedoidsInit::Build),
              LIMIT);
    // One medoid below the boundary 2k still fits, so the guard fires where
    // the arithmetic breaks and not a step early.
    EXPECT_EQ(detail::unavoidable_comparisons(1, half - 1, KMedoidsInit::Build),
              half - 1);

    // A control on the same scale: still saturating would make the cases
    // above prove nothing, so check the guard is not simply always firing.
    // Zero medoids buys no verification, leaving the lone init sweep.
    EXPECT_EQ(detail::unavoidable_comparisons(widest, 0, KMedoidsInit::FarthestFirst),
              widest * (widest - 1));
}

TEST(KMedoidsStreamingSourceTest, TheDiagonalIsZeroAndCostsNoComparison) {
    CountingComparison comparison(4, GateFacts());
    detail::StreamingDistanceRunner distances(comparison, "k_medoids",
                                              /*num_threads=*/1, /*chunk_size=*/1);
    distances.WithAccessor([](const auto& local) {
        EXPECT_EQ(local(2, 2), 0.0);
    });
    EXPECT_EQ(comparison.Count(), 0u);
}

TEST(KMedoidsStreamingSourceTest, ReadsTheSamePairInBothDirections) {
    // PairCountingComparison throws on any read with i >= j, so reaching the
    // second assertion at all is the proof that the source canonicalises.
    PairCountingComparison comparison(4, Line(4));
    detail::StreamingDistanceRunner distances(comparison, "k_medoids",
                                              /*num_threads=*/1, /*chunk_size=*/1);
    distances.WithAccessor([](const auto& local) {
        EXPECT_EQ(local(1, 3), 2.0);
        EXPECT_EQ(local(3, 1), 2.0);
    });
    EXPECT_EQ(comparison.Total(), 2u);
}

TEST(KMedoidsStreamingSourceTest, RefusesANonFiniteValueNamingThePair) {
    const double nan_value = std::numeric_limits<double>::quiet_NaN();
    ScriptedPairComparison comparison(4, Line(4), 1, 2, {nan_value});
    detail::StreamingDistanceRunner distances(comparison, "k_medoids",
                                              /*num_threads=*/1, /*chunk_size=*/1);
    distances.WithAccessor([](const auto& local) {
        try {
            local(1, 2);
            FAIL() << "a non-finite distance must be refused";
        } catch (const std::runtime_error& error) {
            const std::string message = error.what();
            EXPECT_NE(message.find("1"), std::string::npos);
            EXPECT_NE(message.find("2"), std::string::npos);
        }
    });
}

TEST(KMedoidsStreamingSourceTest, RefusesAnotherPassAfterOneFailed) {
    ConstantComparison comparison(8, 1.0);
    detail::StreamingDistanceRunner distances(comparison, "k_medoids",
                                              /*num_threads=*/2, /*chunk_size=*/2);
    EXPECT_THROW(distances.Run(8,
                               [](const auto&, size_t, size_t) {
                                   throw std::runtime_error("the body failed");
                               }),
                 std::runtime_error);

    // The pool is owned for the runner's lifetime and its cancelled flag is
    // sticky, so without the guard this pass would run no chunk at all and the
    // kernel would read its untouched output as a converged answer.
    bool ran = false;
    EXPECT_THROW(distances.Run(8, [&](const auto&, size_t, size_t) { ran = true; }),
                 std::logic_error);
    EXPECT_FALSE(ran);
}

// A chunk size of zero is rejected by validate_options today, so this is not a
// reachable input -- but the matrix runner hands its accessor straight to the
// body whatever the chunk size, and a serial pass on the calling thread has no
// chunks to size anyway. Leaving the streaming runner to answer differently
// would make a silent no-op wait for the next entry point to reach it, and a
// skipped body is the one failure RunGuarded cannot catch: nothing throws.
TEST(KMedoidsStreamingSourceTest, RunsSerialWorkEvenWithoutAUsableChunkSize) {
    ConstantComparison comparison(4, 1.0);
    detail::StreamingDistanceRunner distances(comparison, "k_medoids",
                                              /*num_threads=*/1, /*chunk_size=*/0);
    bool ran = false;
    distances.WithAccessor([&](const auto& local) {
        ran = true;
        EXPECT_EQ(local(0, 1), 1.0);
    });
    EXPECT_TRUE(ran);
}

// Bit-identical results across num_threads and chunk_size is the guarantee the
// streaming path exists to keep, and it rests entirely on the two runners
// splitting a pass the same way: a kernel sums a candidate's row within one
// chunk, so a different chunk boundary is a different summation order. Nothing
// else in the suite would notice that, because every other test asks whether a
// run is correct, not whether it is partitioned identically to the other input.
TEST(KMedoidsStreamingSourceTest, PartitionsAPassExactlyAsTheMatrixRunnerDoes) {
    const std::vector<size_t> sizes = {1, 2, 3, 5, 8, 13, 40};
    const std::vector<size_t> threads = {0, 1, 2, 3, 8};
    const std::vector<size_t> chunks = {1, 2, 3, 7, 4096,
                                        std::numeric_limits<size_t>::max()};
    size_t configs = 0;
    for (const size_t n : sizes) {
        const std::vector<double> condensed = Line(n);
        for (const size_t num_threads : threads) {
            for (const size_t chunk_size : chunks) {
                for (const size_t range : {n, n / 2 + 1}) {
                    SCOPED_TRACE(testing::Message()
                                 << "n=" << n << " num_threads=" << num_threads
                                 << " chunk_size=" << chunk_size
                                 << " range=" << range);
                    std::mutex chunk_mutex;
                    std::vector<std::pair<size_t, size_t>> matrix_chunks;
                    std::vector<std::pair<size_t, size_t>> stream_chunks;

                    detail::MatrixDistanceRunner matrix(
                        condensed.empty() ? nullptr : condensed.data(), n,
                        num_threads, chunk_size);
                    matrix.Run(range, [&](const auto&, size_t begin, size_t end) {
                        std::lock_guard<std::mutex> lock(chunk_mutex);
                        matrix_chunks.emplace_back(begin, end);
                    });

                    ConstantComparison comparison(n, 1.0);
                    detail::StreamingDistanceRunner stream(comparison, "probe",
                                                           num_threads, chunk_size);
                    stream.Run(range, [&](const auto&, size_t begin, size_t end) {
                        std::lock_guard<std::mutex> lock(chunk_mutex);
                        stream_chunks.emplace_back(begin, end);
                    });

                    std::sort(matrix_chunks.begin(), matrix_chunks.end());
                    std::sort(stream_chunks.begin(), stream_chunks.end());
                    EXPECT_EQ(matrix_chunks, stream_chunks);
                    EXPECT_EQ(matrix.ChunkSize(range), stream.ChunkSize(range));
                    ++configs;
                }
            }
        }
    }
    // A mistyped loop bound would make every assertion above vacuous.
    EXPECT_EQ(configs, 420u);
}

// The central claim: streamed equals the matrix filled through the same
// Compare, bit for bit, at every thread count and chunk size.
TEST(KMedoidsStreamingTest, MatchesTheMatrixPathBitwise) {
    for (size_t n : {size_t{6}, size_t{17}, size_t{40}}) {
        // PairCountingComparison additionally throws on any i >= j read, so
        // this sweep doubles as the orientation proof.
        PairCountingComparison comparison(n, Line(n));
        const DenseStorage storage = CompareFilled(comparison);
        for (size_t k : {size_t{1}, size_t{2}, size_t{5}}) {
            if (k >= n) {
                continue;
            }
            for (auto init : {KMedoidsInit::Build, KMedoidsInit::FarthestFirst}) {
                for (size_t threads : {size_t{1}, size_t{4}}) {
                    for (size_t chunk : {size_t{1}, size_t{4096}}) {
                        SCOPED_TRACE(testing::Message()
                                     << "n=" << n << " k=" << k << " init="
                                     << static_cast<int>(init)
                                     << " threads=" << threads
                                     << " chunk=" << chunk);
                        KMedoidsOptions options;
                        options.n_clusters = k;
                        options.init = init;
                        options.num_threads = threads;
                        options.chunk_size = chunk;

                        const KMedoidsResult expected =
                            k_medoids_cluster(storage, options);
                        const KMedoidsResult actual =
                            k_medoids_cluster(comparison, options);

                        EXPECT_EQ(actual.Labels(), expected.Labels());
                        EXPECT_EQ(actual.Medoids(), expected.Medoids());
                        EXPECT_EQ(actual.Cost(), expected.Cost());
                        EXPECT_EQ(actual.NumIterations(), expected.NumIterations());
                        EXPECT_EQ(actual.Converged(), expected.Converged());
                    }
                }
            }
        }
    }
}

TEST(KMedoidsStreamingTest, EveryItemItsOwnMedoidCostsNoComparison) {
    CountingComparison comparison(5, GateFacts());
    KMedoidsOptions options;
    options.n_clusters = 5;
    const KMedoidsResult result = k_medoids_cluster(comparison, options);
    EXPECT_EQ(result.NumClusters(), 5u);
    EXPECT_EQ(comparison.Count(), 0u);
}

TEST(KMedoidsStreamingTest, RefusesWhenUnavoidableWorkExceedsTheBudget) {
    ConstantComparison comparison(50, 1.0);
    KMedoidsOptions options;
    options.n_clusters = 4;
    options.init = KMedoidsInit::Build;
    // Unavoidable work is (k + k) sweeps + n * k; one below that refuses.
    options.max_comparisons =
        detail::unavoidable_comparisons(50, 4, KMedoidsInit::Build) - 1;
    EXPECT_THROW(k_medoids_cluster(comparison, options), std::invalid_argument);

    options.max_comparisons =
        detail::unavoidable_comparisons(50, 4, KMedoidsInit::Build);
    EXPECT_NO_THROW(k_medoids_cluster(comparison, options));
}

// The refusal's advice is release-note text, and no rule of thumb about it is
// true across the range: Build's floor is 2k passes over every pair against
// farthest-first's k + 1, so switching saves nothing at k = 1, a quarter at
// k = 2, and under half however large k grows. Pricing both floors for the
// k actually asked for is what makes the advice answerable; suppressing it
// where it saves nothing is what stops it sending a caller in a circle.
TEST(KMedoidsStreamingTest, ARefusalPricesTheAlternativeForTheMedoidCountAsked) {
    const size_t n = 50;
    auto refusal = [n](size_t k, size_t limit) {
        ConstantComparison comparison(n, 1.0);
        KMedoidsOptions options;
        options.n_clusters = k;
        options.init = KMedoidsInit::Build;
        options.max_comparisons = limit;
        try {
            k_medoids_cluster(comparison, options);
            return std::string("no refusal");
        } catch (const std::invalid_argument& error) {
            return std::string(error.what());
        }
    };

    // One medoid: the two floors are equal, so the switch would return this
    // very message unchanged and is not offered.
    const std::string single = refusal(1, 100);
    EXPECT_NE(single.find(std::to_string(detail::unavoidable_comparisons(
                  n, 1, KMedoidsInit::Build))),
              std::string::npos);
    EXPECT_EQ(single.find("farthest_first"), std::string::npos);

    // Four medoids, and a limit neither floor clears: the alternative is
    // named with its own number and reported as still out of reach.
    const size_t cheaper =
        detail::unavoidable_comparisons(n, 4, KMedoidsInit::FarthestFirst);
    const std::string unreachable = refusal(4, 100);
    EXPECT_NE(unreachable.find(std::to_string(cheaper)), std::string::npos);
    EXPECT_NE(unreachable.find("still over it"), std::string::npos);

    // The same four medoids against a limit the alternative does clear, which
    // is the one case where the advice is worth acting on.
    const std::string reachable = refusal(4, cheaper + 1);
    EXPECT_NE(reachable.find(std::to_string(cheaper)), std::string::npos);
    EXPECT_NE(reachable.find("fits within this limit"), std::string::npos);
}

TEST(KMedoidsStreamingTest, AnExhaustedBudgetStopsTheSameWayAtEveryThreadCount) {
    KMedoidsResult first;
    for (size_t threads : {size_t{1}, size_t{4}}) {
        PairCountingComparison comparison(40, Line(40));
        KMedoidsOptions options;
        options.n_clusters = 3;
        options.init = KMedoidsInit::FarthestFirst;
        options.num_threads = threads;
        // Enough for the unavoidable part plus exactly one swap sweep.
        options.max_comparisons =
            detail::unavoidable_comparisons(40, 3, KMedoidsInit::FarthestFirst) +
            detail::sweep_comparisons(40);
        const KMedoidsResult result = k_medoids_cluster(comparison, options);
        EXPECT_FALSE(result.Converged());
        EXPECT_LT(result.NumIterations(), options.max_iterations);
        if (threads == 1) {
            first = result;
        } else {
            EXPECT_EQ(result.Labels(), first.Labels());
            EXPECT_EQ(result.Medoids(), first.Medoids());
            EXPECT_EQ(result.Cost(), first.Cost());
        }
    }
}

// The budget's job is to refuse a run it cannot afford, so it must never
// price the floor *below* what the floor costs. For Explicit init it prices
// k verification sweeps plus the n * k cache fill, while the run that just
// fits spends the cache fill alone: the between-passes check stops the loop
// before best_predicted_swap, so no verification pass ever runs. The gap is
// the safety margin, and this pins its direction -- a formula that drifted
// under the real floor would wave through runs it exists to refuse.
TEST(KMedoidsStreamingTest, AnExplicitRunAtItsPricedFloorOnlyFillsTheCache) {
    const size_t n = 50;
    const size_t k = 4;
    CountingComparison comparison(n, GateFacts());
    KMedoidsOptions options;
    options.n_clusters = k;
    options.init = KMedoidsInit::Explicit;
    options.initial_medoids = {0, 1, 2, 3};
    options.max_comparisons =
        detail::unavoidable_comparisons(n, k, KMedoidsInit::Explicit);

    const KMedoidsResult result = k_medoids_cluster(comparison, options);
    EXPECT_FALSE(result.Converged());
    EXPECT_EQ(result.NumIterations(), 0u);
    // One column per slot, less the k diagonal entries the source answers
    // without scoring.
    EXPECT_EQ(comparison.Count(), n * k - k);
    EXPECT_LT(comparison.Count(), options.max_comparisons);
}

// The defect this guards, and the reason the loop checks twice: the
// verification pass rescans every item for each of the k slots, so it costs k
// sweeps. Charged as one, a budget that admitted it was overrun by a factor
// of k -- at exactly this shape, 11466 comparisons spent against a limit of
// 5100. Priced properly the pass does not start unless all k sweeps fit, and
// that limit is refused outright instead.
TEST(KMedoidsStreamingTest, TheVerificationPassOnlyRunsWhenItsKSweepsFit) {
    const size_t n = 50;
    const size_t k = 4;
    const size_t floor =
        detail::unavoidable_comparisons(n, k, KMedoidsInit::Explicit);
    const size_t sweep = detail::sweep_comparisons(n);
    const size_t verification = detail::verification_comparisons(n, k);

    KMedoidsOptions options;
    options.n_clusters = k;
    options.init = KMedoidsInit::Explicit;
    options.initial_medoids = {0, 1, 2, 3};

    // One predicted-delta scan fits, the verification pass behind it does not,
    // so the run stops between the two rather than starting a pass it cannot
    // pay for.
    CountingComparison stopped(n, GateFacts());
    options.max_comparisons = floor + sweep;
    const KMedoidsResult early = k_medoids_cluster(stopped, options);
    EXPECT_FALSE(early.Converged());
    EXPECT_EQ(early.NumIterations(), 0u);
    EXPECT_LE(stopped.Count(), options.max_comparisons);

    // k sweeps more and the pass runs. Convergence is the witness that it did:
    // only verification_pass can set it.
    CountingComparison verified(n, GateFacts());
    options.max_comparisons = floor + sweep + verification;
    const KMedoidsResult full = k_medoids_cluster(verified, options);
    EXPECT_TRUE(full.Converged());
    EXPECT_LE(verified.Count(), options.max_comparisons);
    EXPECT_GT(verified.Count(), stopped.Count());
}

// Both halves of what the floor does and does not pay for, at the shape that
// used to expose it: k == 1 with BUILD, where no verification sweep runs and
// so no over-pricing is left over to absorb a mistake. The floor now covers
// the init passes -- it prices k of them, and BUILD really makes k, the
// global-medoid pass included. What it still leaves out is the O(n * k) row
// bookkeeping around those passes: BUILD refreshes its nearest-distance table
// once per step. That remainder is deliberately unpriced, so it is pinned at
// its exact size rather than left to drift.
TEST(KMedoidsStreamingTest, ABuildFloorCoversItsInitPassesButNotItsRowBookkeeping) {
    const size_t n = 50;
    const size_t k = 1;
    CountingComparison comparison(n, GateFacts());
    KMedoidsOptions options;
    options.n_clusters = k;
    options.init = KMedoidsInit::Build;
    options.max_comparisons =
        detail::unavoidable_comparisons(n, k, KMedoidsInit::Build);

    const KMedoidsResult result = k_medoids_cluster(comparison, options);
    EXPECT_FALSE(result.Converged());
    // The global-medoid pass, the nearest-table fill, and the cache column.
    EXPECT_EQ(comparison.Count(), detail::sweep_comparisons(n) + 2 * (n - 1));
    // Covered, where an init priced at k - 1 passes left it 48 over.
    EXPECT_LE(comparison.Count(), options.max_comparisons);
    // And the unpriced remainder, over the one init pass and the cache fill
    // the formula does cover, is exactly the row bookkeeping.
    const size_t priced_before_verification =
        k * detail::sweep_comparisons(n) + n * k;
    EXPECT_EQ(comparison.Count() - priced_before_verification, n - 2);
}

// The ledger on the path it actually has to hold on. Every other counting
// test here scores zero for every pair, which converges without accepting a
// single swap, so the accepted-swap branch -- the one that refreshes a medoid
// column, rebuilds the assignment cache and comes round for another sweep --
// has never had its spend checked. On a line fixture the run accepts two
// swaps before it verifies, and the spend stays inside the limit at every
// budget from the floor up.
TEST(KMedoidsStreamingTest, AnAcceptedSwapRunStaysWithinEveryBudget) {
    const size_t n = 40;
    const size_t k = 3;
    const size_t floor =
        detail::unavoidable_comparisons(n, k, KMedoidsInit::FarthestFirst);
    const size_t sweep = detail::sweep_comparisons(n);

    size_t swapping_runs = 0;
    size_t stopped_after_swapping = 0;
    size_t converged_runs = 0;
    for (size_t passes = 0; passes <= 8; ++passes) {
        SCOPED_TRACE(testing::Message() << "extra passes=" << passes);
        // Throws on any i >= j read, so the whole sweep is also the
        // orientation proof on a multi-iteration run.
        PairCountingComparison comparison(n, Line(n));
        KMedoidsOptions options;
        options.n_clusters = k;
        options.init = KMedoidsInit::FarthestFirst;
        options.max_comparisons = floor + passes * sweep;

        const KMedoidsResult result = k_medoids_cluster(comparison, options);
        EXPECT_LE(comparison.Total(), options.max_comparisons);
        if (result.NumIterations() > 0) {
            ++swapping_runs;
            if (!result.Converged()) {
                ++stopped_after_swapping;
            }
        }
        if (result.Converged()) {
            ++converged_runs;
        }
    }

    // Without these the assertions above could all be covering runs that
    // stopped at the floor, which is what the zero-distance tests already
    // cover. The middle count is the case this test exists for: swaps were
    // accepted and the budget then stopped the run anyway.
    EXPECT_GT(swapping_runs, 0u);
    EXPECT_GT(stopped_after_swapping, 0u);
    EXPECT_GT(converged_runs, 0u);
}

TEST(KMedoidsStreamingTest, RefusesAnUnrepeatableComparison) {
    NamedRocsComparison rocs(6, Line(6));
    KMedoidsOptions options;
    options.n_clusters = 2;
    EXPECT_THROW(k_medoids_cluster(rocs, options), ComparisonError);
    // Refused on its declared name, before any pair is scored.
    EXPECT_EQ(rocs.Total(), 0u);
}

TEST(KMedoidsStreamingTest, RefusesAnExplicitMedoidOutsideTheItems) {
    PairCountingComparison comparison(6, Line(6));
    KMedoidsOptions options;
    options.n_clusters = 2;
    options.init = KMedoidsInit::Explicit;
    options.initial_medoids = {0, 99};
    EXPECT_THROW(k_medoids_cluster(comparison, options), std::out_of_range);
}
