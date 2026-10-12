/**
 * @file test_kmedoids_streaming.cpp
 * @brief The streaming distance source and the budget arithmetic.
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

#include "oecluster/GateFacts.h"
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

TEST(KMedoidsBudgetTest, UnavoidableWorkCountsInitAndVerificationAndTheCache) {
    const size_t n = 100;
    const size_t k = 5;
    const size_t sweep = detail::sweep_comparisons(n);
    // Build: k - 1 init sweeps, one verification sweep, one cache fill.
    EXPECT_EQ(detail::unavoidable_comparisons(n, k, KMedoidsInit::Build),
              (k - 1) * sweep + sweep + n * k);
    // FarthestFirst: about one init sweep.
    EXPECT_EQ(detail::unavoidable_comparisons(n, k, KMedoidsInit::FarthestFirst),
              sweep + sweep + n * k);
    // Explicit: no init sweep at all.
    EXPECT_EQ(detail::unavoidable_comparisons(n, k, KMedoidsInit::Explicit),
              sweep + n * k);
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

    // Sweeps product: one init sweep plus the verification sweep is twice a
    // sweep that already fills most of the range. Unguarded it would wrap to
    // roughly 2^64 - 2^33 and read as a run worth attempting.
    EXPECT_EQ(detail::unavoidable_comparisons(widest, 1, KMedoidsInit::FarthestFirst),
              LIMIT);

    // Cache-fill product: Explicit takes no init sweep, so the sweeps total
    // survives, and the n * k cache fill is what tips it over. Unguarded the
    // sum lands exactly on 2^64 and wraps to zero -- free work.
    EXPECT_EQ(detail::unavoidable_comparisons(widest, 1, KMedoidsInit::Explicit),
              LIMIT);

    // A control on the same scale: still saturating would make the two cases
    // above prove nothing, so check the guard is not simply always firing.
    EXPECT_EQ(detail::unavoidable_comparisons(widest, 0, KMedoidsInit::Explicit),
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
