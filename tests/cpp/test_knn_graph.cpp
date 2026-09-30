/**
 * @file test_knn_graph.cpp
 * @brief knn_graph over distance matrices and comparisons, and the KNNGraph
 * constructor's validation.
 */

#include <gtest/gtest.h>

#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "oecluster/Error.h"
#include "oecluster/GateFacts.h"
#include "oecluster/StorageBackend.h"
#include "oecluster/clustering/KNNGraph.h"

#include "diversity_test_support.h"
#include "knn_test_support.h"

using namespace OECluster;
using namespace diversity_test;
using namespace knn_test;

namespace {

KNNGraphOptions Options(size_t k, size_t num_threads = 0,
                        size_t chunk_size = 4096) {
    KNNGraphOptions options;
    options.k = k;
    options.num_threads = num_threads;
    options.chunk_size = chunk_size;
    return options;
}

const size_t CHUNK_SIZES[] = {1, 7, 4096, std::numeric_limits<size_t>::max()};
const size_t THREAD_COUNTS[] = {1, 4};

}  // namespace

TEST(KNNGraphTest, DenseMatchesTheBruteForceOracle) {
    for (const size_t n : {size_t{2}, size_t{3}, size_t{7}, size_t{16}}) {
        for (const unsigned seed : {1u, 2u, 3u}) {
            const std::vector<double> condensed = Quantized(n, seed, 4);
            const DenseStorage storage = MakeStorage(n, condensed);
            for (size_t k = 1; k < n; ++k) {
                ExpectSameGraph(knn_graph(storage, Options(k)),
                                Oracle(n, condensed, k));
            }
        }
    }
}

TEST(KNNGraphTest, EqualDistancesBreakTiesByTheLowerIndex) {
    const size_t n = 5;
    const std::vector<double> condensed(n * (n - 1) / 2, 1.0);
    const KNNGraph graph = knn_graph(MakeStorage(n, condensed), Options(2));
    EXPECT_EQ(graph.Indices(),
              (std::vector<size_t>{1, 2, 0, 2, 0, 1, 0, 1, 0, 1}));
    EXPECT_EQ(graph.Distances(), std::vector<double>(n * 2, 1.0));
}

TEST(KNNGraphTest, RowsAreOrderedByDistanceThenIndex) {
    // Items on a line at 0, 1, 2, 4, 7, 8: row 2 ties 0 and 3 at distance 2.
    const std::vector<double> condensed = Positions({0, 1, 2, 4, 7, 8});
    const KNNGraph graph = knn_graph(MakeStorage(6, condensed), Options(2));
    EXPECT_EQ(graph.NumItems(), 6u);
    EXPECT_EQ(graph.K(), 2u);
    EXPECT_EQ(graph.Indices(),
              (std::vector<size_t>{1, 2, 0, 2, 1, 0, 2, 1, 5, 3, 4, 3}));
    EXPECT_EQ(graph.Distances(),
              (std::vector<double>{1, 2, 1, 1, 1, 2, 2, 3, 1, 3, 1, 4}));
}

TEST(KNNGraphTest, ResultIsIndependentOfThreadsAndChunkSize) {
    const size_t n = 23;
    const std::vector<double> condensed = Quantized(n, 7, 5);
    const DenseStorage storage = MakeStorage(n, condensed);
    const KNNGraph expected = Oracle(n, condensed, 4);
    for (const size_t threads : THREAD_COUNTS) {
        for (const size_t chunk : CHUNK_SIZES) {
            ExpectSameGraph(knn_graph(storage, Options(4, threads, chunk)),
                            expected);
        }
    }
}

TEST(KNNGraphTest, MemoryMappedStorageMatchesDense) {
    const size_t n = 19;
    const std::vector<double> condensed = Quantized(n, 11, 5);
    const TempMMap mmap("oecluster_knn_graph_mmap.bin", n, condensed);
    const KNNGraph expected = Oracle(n, condensed, 3);
    for (const size_t threads : THREAD_COUNTS) {
        for (const size_t chunk : CHUNK_SIZES) {
            ExpectSameGraph(knn_graph(mmap.Storage(), Options(3, threads, chunk)),
                            expected);
        }
    }
}

TEST(KNNGraphTest, KOutsideOneToNMinusOneIsRefused) {
    const DenseStorage storage = MakeStorage(4, Line(4));
    ExpectInvalidArgument([&] { knn_graph(storage, Options(0)); },
                          "knn_graph k must be between 1 and 3 for 4 items, got 0");
    ExpectInvalidArgument([&] { knn_graph(storage, Options(4)); },
                          "knn_graph k must be between 1 and 3 for 4 items, got 4");
    ExpectInvalidArgument(
        [] { knn_graph(DenseStorage(1), Options(1)); },
        "knn_graph needs at least two items: a single item has no neighbors");
}

TEST(KNNGraphTest, ZeroItemsGiveAnEmptyGraphForAnyK) {
    for (const size_t k : {size_t{0}, size_t{1}, size_t{5}}) {
        const KNNGraph graph = knn_graph(DenseStorage(0), Options(k));
        EXPECT_EQ(graph.NumItems(), 0u);
        EXPECT_EQ(graph.K(), k);
        EXPECT_TRUE(graph.Indices().empty());
        EXPECT_TRUE(graph.Distances().empty());
    }
}

TEST(KNNGraphTest, ZeroChunkSizeIsRefusedEvenWithoutItems) {
    ExpectInvalidArgument(
        [] { knn_graph(MakeStorage(4, Line(4)), Options(1, 0, 0)); },
        "knn_graph chunk_size must be at least one");
    ExpectInvalidArgument([] { knn_graph(DenseStorage(0), Options(1, 0, 0)); },
                          "knn_graph chunk_size must be at least one");
}

TEST(KNNGraphTest, DataLessStorageIsRefused) {
    ExpectInvalidArgument(
        [] { knn_graph(NullDataStorage(4), Options(1)); },
        "knn_graph requires dense, memory-mapped or sparse storage; this "
        "storage has no data array");
}

TEST(KNNGraphTest, NonFiniteDistancesAreRefused) {
    for (const double bad : {NaN, INF, -INF}) {
        std::vector<double> condensed = Line(5);
        condensed[5] = bad;  // the pair (1, 3)
        const DenseStorage storage = MakeStorage(5, condensed);
        ExpectRuntimeError(
            [&] { knn_graph(storage, Options(2, 1)); },
            "knn_graph read a non-finite distance between items 1 and 3");
    }
}

TEST(KNNGraphConstructorTest, AValidGraphRoundTripsItsAccessors) {
    const KNNGraph graph(3, 1, {1, 0, 1}, {0.5, 0.5, 2.0});
    EXPECT_EQ(graph.NumItems(), 3u);
    EXPECT_EQ(graph.K(), 1u);
    EXPECT_EQ(graph.Indices(), (std::vector<size_t>{1, 0, 1}));
    EXPECT_EQ(graph.Distances(), (std::vector<double>{0.5, 0.5, 2.0}));
}

TEST(KNNGraphConstructorTest, TheDefaultGraphIsEmpty) {
    const KNNGraph graph;
    EXPECT_EQ(graph.NumItems(), 0u);
    EXPECT_EQ(graph.K(), 0u);
    EXPECT_TRUE(graph.Indices().empty());
    EXPECT_TRUE(graph.Distances().empty());
}

TEST(KNNGraphConstructorTest, AZeroItemGraphKeepsItsKAndNeedsEmptyArrays) {
    EXPECT_EQ(KNNGraph(0, 7, {}, {}).K(), 7u);
    ExpectInvalidArgument(
        [] { KNNGraph(0, 7, {1}, {1.0}); },
        "KNNGraph needs 0 indices and distances for 0 items at k = 7, got 1 "
        "indices and 1 distances");
}

TEST(KNNGraphConstructorTest, EveryRuleHasARefusingCase) {
    const size_t max = std::numeric_limits<size_t>::max();
    ExpectInvalidArgument([&] { KNNGraph(max, 2, {}, {}); },
                          "KNNGraph num_items * k overflows size_t");
    ExpectInvalidArgument(
        [] { KNNGraph(3, 1, {1, 0}, {1.0, 1.0, 1.0}); },
        "KNNGraph needs 3 indices and distances for 3 items at k = 1, got 2 "
        "indices and 3 distances");
    ExpectInvalidArgument(
        [] { KNNGraph(3, 1, {1, 0, 0}, {1.0, 1.0}); },
        "KNNGraph needs 3 indices and distances for 3 items at k = 1, got 3 "
        "indices and 2 distances");
    ExpectInvalidArgument([] { KNNGraph(3, 0, {}, {}); },
                          "KNNGraph k must be between 1 and 2 for 3 items, got 0");
    ExpectInvalidArgument(
        [] { KNNGraph(3, 3, std::vector<size_t>(9, 0), std::vector<double>(9, 1.0)); },
        "KNNGraph k must be between 1 and 2 for 3 items, got 3");
    ExpectInvalidArgument(
        [] { KNNGraph(1, 1, {0}, {1.0}); },
        "KNNGraph needs at least two items: a single item has no neighbors");
    ExpectInvalidArgument([] { KNNGraph(3, 1, {1, 7, 0}, {1.0, 1.0, 1.0}); },
                          "KNNGraph row 1 names item 7, outside the 3 items");
    ExpectInvalidArgument([] { KNNGraph(3, 1, {1, 1, 0}, {1.0, 1.0, 1.0}); },
                          "KNNGraph row 1 contains its own item");
    ExpectInvalidArgument(
        [] { KNNGraph(3, 2, {1, 2, 0, 0, 0, 1}, {1.0, 2.0, 1.0, 2.0, 1.0, 2.0}); },
        "KNNGraph row 1 repeats item 0");
    ExpectInvalidArgument(
        [] { KNNGraph(3, 1, {1, 0, 0}, {1.0, NaN, 1.0}); },
        "KNNGraph row 1 has a non-finite distance");
    ExpectInvalidArgument(
        [] { KNNGraph(3, 2, {1, 2, 2, 0, 0, 1}, {1.0, 2.0, 1.0, 1.0, 1.0, 2.0}); },
        "KNNGraph row 1 is not ordered by ascending (distance, index)");
    ExpectInvalidArgument(
        [] { KNNGraph(3, 2, {1, 2, 0, 2, 0, 1}, {1.0, 2.0, 2.0, 1.0, 1.0, 2.0}); },
        "KNNGraph row 1 is not ordered by ascending (distance, index)");
}

TEST(KNNGraphTest, EveryPathAgreesForEveryThreadCountAndChunkSize) {
    const size_t n = 17;
    const std::vector<double> condensed = Quantized(n, 5, 4);
    const DenseStorage dense = MakeStorage(n, condensed);
    const TempMMap mmap("oecluster_knn_graph_paths.bin", n, condensed);
    // Quantized distances lie in [0, 1]; a cutoff above that keeps every pair.
    const auto sparse = MakeSparse(n, condensed, 2.0);
    TableComparison comparison(n, condensed);
    const KNNGraph expected = Oracle(n, condensed, 5);
    for (const size_t threads : THREAD_COUNTS) {
        for (const size_t chunk : CHUNK_SIZES) {
            const KNNGraphOptions options = Options(5, threads, chunk);
            ExpectSameGraph(knn_graph(dense, options), expected);
            ExpectSameGraph(knn_graph(mmap.Storage(), options), expected);
            ExpectSameGraph(knn_graph(*sparse, options), expected);
            ExpectSameGraph(knn_graph(comparison, options), expected);
        }
    }
}

TEST(KNNGraphSparseTest, ACutoffBelowSomeDistancesMatchesDenseWhileNoItemIsShort) {
    // d(i, j) = j - i; at cutoff 2 every item keeps at least two neighbors.
    const std::vector<double> condensed = Line(8);
    const auto sparse = MakeSparse(8, condensed, 2.0);
    ExpectSameGraph(knn_graph(*sparse, Options(2)), Oracle(8, condensed, 2));
}

TEST(KNNGraphSparseTest, AShortItemIsRefused) {
    const auto sparse = MakeSparse(8, Line(8), 1.0);
    ExpectInvalidArgument(
        [&] { knn_graph(*sparse, Options(2)); },
        "knn_graph item 0 has only 1 of the k = 2 neighbors it needs within the "
        "sparse cutoff 1; raise the cutoff or lower k");
    const auto half = MakeSparse(8, Line(8), 0.5);
    ExpectInvalidArgument(
        [&] { knn_graph(*half, Options(1)); },
        "knn_graph item 0 has only 0 of the k = 1 neighbors it needs within the "
        "sparse cutoff 0.5; raise the cutoff or lower k");
}

TEST(KNNGraphSparseTest, UnfinalizedStorageIsRefusedAsShort) {
    SparseStorage storage(4, 10.0);
    storage.Set(0, 1, 1.0);
    ExpectInvalidArgument(
        [&] { knn_graph(storage, Options(1)); },
        "knn_graph item 0 has only 0 of the k = 1 neighbors it needs within the "
        "sparse cutoff 10; raise the cutoff or lower k");
}

TEST(KNNGraphSparseTest, DuplicatePairsCountOnceWithTheValueGetReports) {
    for (const double second : {0.5, 3.0}) {
        SparseStorage storage(4, 10.0);
        storage.Set(0, 1, 0.5);
        storage.Set(0, 1, second);
        storage.Set(0, 2, 1.0);
        storage.Set(0, 3, 2.0);
        storage.Set(1, 2, 1.5);
        storage.Set(1, 3, 2.5);
        storage.Set(2, 3, 0.25);
        storage.Finalize();
        const KNNGraph graph = knn_graph(storage, Options(3));
        const double reported = storage.Get(0, 1);
        // Row 0 holds items 1, 2 and 3 exactly once, and pair (0, 1) carries
        // Get()'s value.
        std::vector<size_t> row0(graph.Indices().begin(), graph.Indices().begin() + 3);
        std::sort(row0.begin(), row0.end());
        EXPECT_EQ(row0, (std::vector<size_t>{1, 2, 3}));
        for (size_t m = 0; m < 3; ++m) {
            if (graph.Indices()[m] == 1) {
                EXPECT_EQ(graph.Distances()[m], reported);
            }
            if (graph.Indices()[3 + m] == 0) {
                EXPECT_EQ(graph.Distances()[3 + m], reported);
            }
        }
    }
}

TEST(KNNGraphSparseTest, StoredNonFiniteEntriesAreRefused) {
    for (const double bad : {NaN, -INF}) {
        SparseStorage storage(4, 10.0);
        const std::vector<double> condensed = Line(4);
        size_t k = 0;
        for (size_t i = 0; i < 4; ++i) {
            for (size_t j = i + 1; j < 4; ++j) {
                storage.Set(i, j, (i == 1 && j == 3) ? bad : condensed[k]);
                ++k;
            }
        }
        storage.Finalize();
        ExpectRuntimeError(
            [&] { knn_graph(storage, Options(2, 1)); },
            "knn_graph read a non-finite distance between items 1 and 3");
    }
}

TEST(KNNGraphComparisonTest, ReadsEveryOrderedPairOnceAsCompareMinMax) {
    // TableComparison throws on a reversed pair or the diagonal, so passing
    // pins Compare(min, max) with no self read; the count pins N(N-1) calls.
    const size_t n = 9;
    TableComparison table(n, Quantized(n, 3, 4));
    EXPECT_NO_THROW(knn_graph(table, Options(3, 4, 7)));
    CountingComparison counting(n, GateFacts());
    knn_graph(counting, Options(3, 4, 7));
    EXPECT_EQ(counting.Count(), n * (n - 1));
}

TEST(KNNGraphComparisonTest, ClonesAreIsolatedAndBoundedByTheThreadCount) {
    // Sized like SphereExclusionComparisonTest.ProvesCloneIsolation: one row
    // per unit and 40 units give four workers room to overlap.
    const size_t n = 40;
    const std::vector<double> condensed = Quantized(n, 9, 5);
    IsolationComparison isolation(n, condensed);
    ExpectSameGraph(knn_graph(isolation, Options(4, 4, 1)), Oracle(n, condensed, 4));
    EXPECT_EQ(isolation.Violations(), 0u);
    EXPECT_TRUE(isolation.OverlapObserved())
        << "Overlap not observed; test may be flaky on this machine";
    TableComparison table(n, condensed);
    knn_graph(table, Options(4, 4, 1));
    EXPECT_GE(table.NumClones(), 1u);
    EXPECT_LE(table.NumClones(), 4u);
}

TEST(KNNGraphComparisonTest, FactsThatRuleOutRankingAreRefused) {
    GateFacts similarity;
    similarity.is_distance = Capability::No;
    CountingComparison comparison(5, similarity);
    try {
        knn_graph(comparison, Options(2));
        FAIL() << "expected ComparisonError";
    } catch (const ComparisonError& error) {
        EXPECT_EQ(std::string(error.what()),
                  "knn_graph requires distances, but the comparison reports "
                  "similarities");
    }
    EXPECT_EQ(comparison.Count(), 0u);
}

TEST(KNNGraphComparisonTest, BoundsMatchTheStorageOverload) {
    CountingComparison empty(0, GateFacts());
    EXPECT_EQ(knn_graph(empty, Options(3)).NumItems(), 0u);
    EXPECT_EQ(knn_graph(empty, Options(3)).K(), 3u);
    ExpectInvalidArgument([&] { knn_graph(empty, Options(1, 0, 0)); },
                          "knn_graph chunk_size must be at least one");
    CountingComparison one(1, GateFacts());
    ExpectInvalidArgument(
        [&] { knn_graph(one, Options(1)); },
        "knn_graph needs at least two items: a single item has no neighbors");
    CountingComparison four(4, GateFacts());
    ExpectInvalidArgument([&] { knn_graph(four, Options(4)); },
                          "knn_graph k must be between 1 and 3 for 4 items, got 4");
    EXPECT_EQ(four.Count(), 0u);
}

// Each case breaks two adjacent rules of the validation order and pins that
// the earlier rule's message wins.
TEST(KNNGraphTest, ValidationOrderHoldsWhenSeveralInputsAreInvalid) {
    // chunk_size before the item count.
    ExpectInvalidArgument([] { knn_graph(DenseStorage(1), Options(1, 0, 0)); },
                          "knn_graph chunk_size must be at least one");
    // k before the data-array check.
    ExpectInvalidArgument([] { knn_graph(NullDataStorage(4), Options(0)); },
                          "knn_graph k must be between 1 and 3 for 4 items, got 0");
    // k before the sparse short-item scan.
    const auto sparse = MakeSparse(8, Line(8), 0.5);
    ExpectInvalidArgument([&] { knn_graph(*sparse, Options(8)); },
                          "knn_graph k must be between 1 and 7 for 8 items, got 8");
    // The short-item scan before any distance read.
    SparseStorage short_and_nan(4, 10.0);
    short_and_nan.Set(0, 1, NaN);
    short_and_nan.Finalize();
    ExpectInvalidArgument(
        [&] { knn_graph(short_and_nan, Options(2)); },
        "knn_graph item 0 has only 1 of the k = 2 neighbors it needs within the "
        "sparse cutoff 10; raise the cutoff or lower k");
    // k before the comparison facts.
    GateFacts similarity;
    similarity.is_distance = Capability::No;
    CountingComparison counting(4, similarity);
    ExpectInvalidArgument([&] { knn_graph(counting, Options(4)); },
                          "knn_graph k must be between 1 and 3 for 4 items, got 4");
    EXPECT_EQ(counting.Count(), 0u);
    // The comparison facts before any distance read.
    std::vector<double> condensed = Line(5);
    condensed[5] = NaN;
    TableComparison table(5, condensed, similarity);
    EXPECT_THROW(knn_graph(table, Options(2)), ComparisonError);
}

TEST(KNNGraphComparisonTest, NonFiniteDistancesAreRefused) {
    for (const double bad : {NaN, INF}) {
        std::vector<double> condensed = Line(5);
        condensed[5] = bad;  // the pair (1, 3)
        TableComparison comparison(5, condensed);
        ExpectRuntimeError(
            [&] { knn_graph(comparison, Options(2, 1)); },
            "knn_graph read a non-finite distance between items 1 and 3");
    }
}
