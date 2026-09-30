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
