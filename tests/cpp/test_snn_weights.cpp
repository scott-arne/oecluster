/**
 * @file test_snn_weights.cpp
 * @brief Shared-nearest-neighbor Jaccard weights over a KNNGraph.
 */

#include <gtest/gtest.h>

#include <climits>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "oecluster/clustering/KNNGraph.h"

#include "../../src/clustering/SNNWeights.h"
#include "diversity_test_support.h"
#include "knn_test_support.h"

using namespace OECluster;
using namespace diversity_test;
using namespace knn_test;

namespace {

// Rows 0 -> {1, 2}, 1 -> {0, 2}, 2 -> {1, 3}, 3 -> {2, 1}. With N+(i) the row
// plus i: N+(0) = {0,1,2}, N+(1) = {0,1,2}, N+(2) = {1,2,3}, N+(3) = {1,2,3}.
// Pairs 0-1 and 2-3 share 3 of 6 slots, weight 3 / (6 - 3) = 1; pairs 0-2,
// 1-2 and 1-3 share 2, weight 2 / (6 - 2) = 0.5. Arcs 0-1, 1-2 and 2-3
// are mutual; 0 -> 2 and 3 -> 1 are one-way.
KNNGraph HandGraph() {
    return KNNGraph(4, 2, {1, 2, 0, 2, 1, 3, 2, 1},
                    {1.0, 2.0, 1.0, 2.0, 1.0, 2.0, 1.0, 2.0});
}

std::map<std::pair<uint32_t, uint32_t>, double> Edges(
    const detail::WeightedGraph& graph) {
    std::map<std::pair<uint32_t, uint32_t>, double> edges;
    for (size_t i = 0; i < graph.num_nodes; ++i) {
        for (size_t e = graph.offsets[i]; e < graph.offsets[i + 1]; ++e) {
            edges[{static_cast<uint32_t>(i), graph.neighbors[e]}] = graph.weights[e];
        }
    }
    return edges;
}

void ExpectSameWeights(const detail::WeightedGraph& actual,
                       const detail::WeightedGraph& expected) {
    EXPECT_EQ(actual.num_nodes, expected.num_nodes);
    EXPECT_EQ(actual.offsets, expected.offsets);
    EXPECT_EQ(actual.neighbors, expected.neighbors);
    EXPECT_EQ(actual.weights, expected.weights);
}

}  // namespace

TEST(SNNWeightsTest, HandWorkedWeightsOverTheUnionEdgeSet) {
    const detail::WeightedGraph graph = detail::snn_weights(HandGraph(), 0.5, 1);
    EXPECT_EQ(graph.num_nodes, 4u);
    EXPECT_EQ(graph.offsets, (std::vector<size_t>{0, 2, 5, 8, 10}));
    EXPECT_EQ(graph.neighbors,
              (std::vector<uint32_t>{1, 2, 0, 2, 3, 0, 1, 3, 1, 2}));
    EXPECT_EQ(graph.weights, (std::vector<double>{1.0, 0.5, 1.0, 0.5, 0.5, 0.5,
                                                  0.5, 1.0, 0.5, 1.0}));
}

TEST(SNNWeightsTest, AWeightAtPruneIsKeptAndOneBelowIsDropped) {
    const detail::WeightedGraph kept = detail::snn_weights(HandGraph(), 0.5, 1);
    EXPECT_EQ(kept.offsets[4], 10u);
    const detail::WeightedGraph pruned =
        detail::snn_weights(HandGraph(), std::nextafter(0.5, 1.0), 1);
    EXPECT_EQ(pruned.offsets, (std::vector<size_t>{0, 1, 2, 3, 4}));
    EXPECT_EQ(pruned.neighbors, (std::vector<uint32_t>{1, 0, 3, 2}));
    EXPECT_EQ(pruned.weights, (std::vector<double>{1.0, 1.0, 1.0, 1.0}));
}

TEST(SNNWeightsTest, MatchesABruteForceJaccardOverTheUnionEdgeSet) {
    for (const unsigned seed : {1u, 2u, 3u}) {
        const size_t n = 40;
        const size_t k = 6;
        const KNNGraph graph = Oracle(n, Quantized(n, seed, 9), k);
        std::vector<std::set<size_t>> hood(n);
        for (size_t i = 0; i < n; ++i) {
            hood[i].insert(i);
            for (size_t m = 0; m < k; ++m) {
                hood[i].insert(graph.Indices()[i * k + m]);
            }
        }
        std::map<std::pair<uint32_t, uint32_t>, double> expected;
        for (size_t i = 0; i < n; ++i) {
            for (size_t m = 0; m < k; ++m) {
                const size_t j = graph.Indices()[i * k + m];
                size_t shared = 0;
                for (const size_t x : hood[i]) {
                    shared += hood[j].count(x);
                }
                const double weight = static_cast<double>(shared) /
                                      static_cast<double>(2 * (k + 1) - shared);
                if (weight >= 0.2) {
                    expected[{static_cast<uint32_t>(i), static_cast<uint32_t>(j)}] =
                        weight;
                    expected[{static_cast<uint32_t>(j), static_cast<uint32_t>(i)}] =
                        weight;
                }
            }
        }
        EXPECT_EQ(Edges(detail::snn_weights(graph, 0.2, 2)), expected)
            << "seed " << seed;
    }
}

TEST(SNNWeightsTest, TheGraphIsSymmetricWithSortedRowsAndNoSelfLoops) {
    const size_t n = 50;
    const detail::WeightedGraph graph =
        detail::snn_weights(Oracle(n, Quantized(n, 4, 7), 5), 0.0, 3);
    const auto edges = Edges(graph);
    for (const auto& [pair, weight] : edges) {
        EXPECT_NE(pair.first, pair.second);
        const auto mirror = edges.find({pair.second, pair.first});
        ASSERT_NE(mirror, edges.end());
        EXPECT_EQ(mirror->second, weight);
    }
    for (size_t i = 0; i < n; ++i) {
        for (size_t e = graph.offsets[i] + 1; e < graph.offsets[i + 1]; ++e) {
            EXPECT_LT(graph.neighbors[e - 1], graph.neighbors[e]);
        }
    }
}

// With k = 30 a unit holds 65536 / (30 * 31) = 70 rows, so 420 rows make six
// units and up to six workers write slots at once; this is also the case the
// TSAN run relies on for real concurrency.
TEST(SNNWeightsTest, TheOutputDoesNotDependOnTheThreadCount) {
    const size_t n = 420;
    const KNNGraph graph = Oracle(n, Quantized(n, 5, 11), 30);
    const detail::WeightedGraph serial = detail::snn_weights(graph, 1.0 / 15.0, 1);
    for (const size_t threads : {size_t{2}, size_t{4}, size_t{0}}) {
        ExpectSameWeights(detail::snn_weights(graph, 1.0 / 15.0, threads), serial);
    }
}

TEST(SNNWeightsTest, ZeroItemsGiveAZeroNodeGraph) {
    const detail::WeightedGraph graph = detail::snn_weights(KNNGraph(), 0.1, 4);
    EXPECT_EQ(graph.num_nodes, 0u);
    EXPECT_EQ(graph.offsets, (std::vector<size_t>{0}));
    EXPECT_TRUE(graph.neighbors.empty());
    EXPECT_TRUE(graph.weights.empty());
}

TEST(SNNWeightsTest, MoreThanIntMaxItemsAreRefused) {
    EXPECT_NO_THROW(
        detail::validate_leiden_item_count(static_cast<size_t>(INT_MAX), "leiden"));
    ExpectInvalidArgument(
        [] {
            detail::validate_leiden_item_count(static_cast<size_t>(INT_MAX) + 1,
                                               "leiden");
        },
        "leiden supports at most 2147483647 items, got 2147483648");
}
