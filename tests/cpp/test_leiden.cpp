/**
 * @file test_leiden.cpp
 * @brief The Leiden engine phase by phase, and leiden's public overloads.
 */

#include <gtest/gtest.h>

#include <algorithm>
#include <climits>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numeric>
#include <queue>
#include <set>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "oecluster/GateFacts.h"
#include "oecluster/StorageBackend.h"
#include "oecluster/clustering/KNNGraph.h"
#include "oecluster/clustering/Leiden.h"

#include "../../src/clustering/LeidenEngine.h"
#include "../../src/clustering/SNNWeights.h"
#include "diversity_test_support.h"
#include "knn_test_support.h"

using namespace OECluster;
using namespace OECluster::detail;
using namespace diversity_test;
using namespace knn_test;

namespace {

using Edge = std::tuple<uint32_t, uint32_t, double>;

WeightedGraph Graph(size_t n, const std::vector<Edge>& edges) {
    std::vector<std::vector<std::pair<uint32_t, double>>> rows(n);
    for (const auto& [a, b, w] : edges) {
        rows[a].emplace_back(b, w);
        rows[b].emplace_back(a, w);
    }
    WeightedGraph graph;
    graph.num_nodes = n;
    graph.offsets.push_back(0);
    for (auto& row : rows) {
        std::sort(row.begin(), row.end());
        for (const auto& [neighbor, weight] : row) {
            graph.neighbors.push_back(neighbor);
            graph.weights.push_back(weight);
        }
        graph.offsets.push_back(graph.neighbors.size());
    }
    return graph;
}

LeidenLevel Level(size_t n, const std::vector<Edge>& edges) {
    return make_base_level(Graph(n, edges));
}

void AddClique(std::vector<Edge>& edges, uint32_t first, uint32_t count,
               double weight = 1.0) {
    for (uint32_t a = first; a < first + count; ++a) {
        for (uint32_t b = a + 1; b < first + count; ++b) {
            edges.emplace_back(a, b, weight);
        }
    }
}

LeidenParams Params(LeidenObjective objective, double resolution,
                    double theta = 0.01) {
    LeidenParams params;
    params.objective = objective;
    params.resolution = resolution;
    params.theta = theta;
    return params;
}

std::vector<uint32_t> Iota(size_t n) {
    std::vector<uint32_t> values(n);
    std::iota(values.begin(), values.end(), uint32_t{0});
    return values;
}

std::vector<ClusterLabel> Blocks(const std::vector<size_t>& sizes) {
    std::vector<ClusterLabel> labels;
    for (size_t block = 0; block < sizes.size(); ++block) {
        labels.insert(labels.end(), sizes[block], static_cast<ClusterLabel>(block));
    }
    return labels;
}

// Whether the nodes labeled `label` induce a connected subgraph.
template <typename Label>
bool InducesConnected(const WeightedGraph& graph, const std::vector<Label>& labels,
                      Label label) {
    std::vector<char> seen(graph.num_nodes, 0);
    size_t start = graph.num_nodes;
    size_t members = 0;
    for (size_t v = 0; v < graph.num_nodes; ++v) {
        if (labels[v] == label) {
            ++members;
            start = std::min(start, v);
        }
    }
    if (members == 0) {
        return true;
    }
    std::queue<size_t> frontier;
    frontier.push(start);
    seen[start] = 1;
    size_t reached = 1;
    while (!frontier.empty()) {
        const size_t v = frontier.front();
        frontier.pop();
        for (size_t e = graph.offsets[v]; e < graph.offsets[v + 1]; ++e) {
            const uint32_t u = graph.neighbors[e];
            if (labels[u] == label && !seen[u]) {
                seen[u] = 1;
                ++reached;
                frontier.push(u);
            }
        }
    }
    return reached == members;
}

void ExpectAllClustersConnected(const WeightedGraph& graph,
                                const std::vector<ClusterLabel>& labels) {
    const std::set<ClusterLabel> distinct(labels.begin(), labels.end());
    for (const ClusterLabel label : distinct) {
        EXPECT_TRUE(InducesConnected(graph, labels, label)) << "cluster " << label;
    }
}

// Two 5-cliques, 0-4 and 5-9, joined by the bridge 4-5.
WeightedGraph Barbell() {
    std::vector<Edge> edges;
    AddClique(edges, 0, 5);
    AddClique(edges, 5, 5);
    edges.emplace_back(4, 5, 1.0);
    return Graph(10, edges);
}

// Six 5-cliques in a ring, clique c holding 5c..5c+4, each joined to the
// next by one edge.
WeightedGraph RingOfCliques() {
    std::vector<Edge> edges;
    for (uint32_t c = 0; c < 6; ++c) {
        AddClique(edges, 5 * c, 5);
        edges.emplace_back(5 * c + 4, 5 * ((c + 1) % 6), 1.0);
    }
    return Graph(30, edges);
}

// Four unit triangles A, B, C, D (0-2, 3-5, 6-8, 9-11). A-B and C-D are each
// joined by all nine cross edges at 0.25; the two halves share no edge.
WeightedGraph CliquesOfCliques() {
    std::vector<Edge> edges;
    for (uint32_t t = 0; t < 4; ++t) {
        AddClique(edges, 3 * t, 3);
    }
    for (const uint32_t pair : {0u, 6u}) {
        for (uint32_t a = pair; a < pair + 3; ++a) {
            for (uint32_t b = pair + 3; b < pair + 6; ++b) {
                edges.emplace_back(a, b, 0.25);
            }
        }
    }
    return Graph(12, edges);
}

WeightedGraph RandomSNN(size_t n, unsigned seed, size_t k) {
    return snn_weights(Oracle(n, Quantized(n, seed, 9), k), 1.0 / 15.0, 1);
}

LeidenOptions Options(size_t k) {
    LeidenOptions options;
    options.k = k;
    return options;
}

void ExpectSameResult(const LeidenResult& actual, const LeidenResult& expected) {
    EXPECT_EQ(actual.Labels(), expected.Labels());
    EXPECT_EQ(actual.Members(), expected.Members());
    EXPECT_EQ(actual.Quality(), expected.Quality());
    EXPECT_EQ(actual.Iterations(), expected.Iterations());
    EXPECT_EQ(actual.Objective(), expected.Objective());
    EXPECT_EQ(actual.Resolution(), expected.Resolution());
    EXPECT_EQ(actual.K(), expected.K());
}

}  // namespace

// ---------------------------------------------------------------------------
// Level, RNG, quality, gain, well-connectedness and candidate selection.

TEST(LeidenEngineTest, TheBaseLevelHasUnitSizesAndNoSelfLoops) {
    const LeidenLevel level = Level(3, {{0, 1, 2.0}, {1, 2, 0.5}});
    EXPECT_EQ(level.num_nodes, 3u);
    EXPECT_EQ(level.self, (std::vector<double>{0.0, 0.0, 0.0}));
    EXPECT_EQ(level.strength, (std::vector<double>{2.0, 2.5, 0.5}));
    EXPECT_EQ(level.size, (std::vector<size_t>{1, 1, 1}));
    EXPECT_EQ(level.m, 2.5);
}

TEST(LeidenEngineTest, TheRngMatchesItsGoldenValues) {
    // std::mt19937_64 fixes the raw stream; the derived values follow from
    // it by the rules in LeidenEngine.h on every standard library.
    LeidenRng next(42);
    EXPECT_EQ(next.Next(), 13930160852258120406ull);
    EXPECT_EQ(next.Next(), 11788048577503494824ull);
    EXPECT_EQ(next.Next(), 13874630024467741450ull);
    LeidenRng below(42);
    std::vector<uint64_t> bounded;
    for (int i = 0; i < 5; ++i) {
        bounded.push_back(below.Below(10));
    }
    EXPECT_EQ(bounded, (std::vector<uint64_t>{6, 4, 0, 2, 1}));
    LeidenRng uniform(42);
    EXPECT_EQ(uniform.Uniform(), 0.75515553295453897);
    EXPECT_EQ(uniform.Uniform(), 0.63903139385469743);
    LeidenRng shuffle(42);
    std::vector<uint32_t> values = Iota(8);
    shuffle.Shuffle(values);
    EXPECT_EQ(values, (std::vector<uint32_t>{7, 0, 5, 1, 2, 4, 3, 6}));
}

TEST(LeidenEngineTest, BoundedIntegersStayInRange) {
    LeidenRng rng(3);
    for (int i = 0; i < 100; ++i) {
        EXPECT_EQ(rng.Below(1), 0u);
    }
    const uint64_t huge = std::numeric_limits<uint64_t>::max() - 1;
    for (int i = 0; i < 100; ++i) {
        EXPECT_LT(rng.Below(huge), huge);
    }
    for (int i = 0; i < 100; ++i) {
        const double u = rng.Uniform();
        EXPECT_GE(u, 0.0);
        EXPECT_LT(u, 1.0);
    }
}

TEST(LeidenEngineTest, QualityMatchesHandComputedValues) {
    // Two disjoint unit edges 0-1 and 2-3; m = 2, every strength 1.
    const LeidenLevel level = Level(4, {{0, 1, 1.0}, {2, 3, 1.0}});
    const std::vector<uint32_t> pairs = {0, 0, 2, 2};
    // Each pair: e = 1, K = 2: 1/2 - (2/4)^2 = 0.25, twice.
    EXPECT_DOUBLE_EQ(leiden_quality(level, pairs, Params(LeidenObjective::Modularity, 1.0)),
                     0.5);
    // Each pair: e = 1, N = 2: 1 - 0.5 * 1 = 0.5, twice.
    EXPECT_DOUBLE_EQ(leiden_quality(level, pairs, Params(LeidenObjective::CPM, 0.5)), 1.0);
    const LeidenLevel empty = Level(3, {});
    EXPECT_EQ(leiden_quality(empty, Iota(3), Params(LeidenObjective::Modularity, 1.0)),
              0.0);
    EXPECT_EQ(leiden_quality(empty, Iota(3), Params(LeidenObjective::CPM, 1.0)), 0.0);
}

TEST(LeidenEngineTest, TheGainOfAMoveEqualsTheChangeInQuality) {
    std::vector<Edge> edges = {{0, 1, 1.0}, {0, 2, 0.5}, {1, 2, 2.0},
                               {2, 3, 1.5}, {3, 4, 1.0}, {1, 4, 0.25}};
    const LeidenLevel level = Level(5, edges);
    // Node 0 alone; community 1 holds nodes 1 and 2.
    const std::vector<uint32_t> before = {0, 1, 1, 3, 3};
    const std::vector<uint32_t> after = {1, 1, 1, 3, 3};
    const double w = 1.0 + 0.5;
    for (const double r : {0.5, 1.0, 1.7}) {
        const LeidenParams modularity = Params(LeidenObjective::Modularity, r);
        const double community_strength = level.strength[1] + level.strength[2];
        const double gain = leiden_gain(w, level.strength[0], 1.0,
                                        community_strength, 2.0, level.m, modularity);
        const double change = leiden_quality(level, after, modularity) -
                              leiden_quality(level, before, modularity);
        EXPECT_NEAR(gain, level.m * change, 1e-12) << "modularity r " << r;
        const LeidenParams cpm = Params(LeidenObjective::CPM, r);
        const double cpm_gain =
            leiden_gain(w, level.strength[0], 1.0, community_strength, 2.0, level.m, cpm);
        const double cpm_change =
            leiden_quality(level, after, cpm) - leiden_quality(level, before, cpm);
        EXPECT_NEAR(cpm_gain, cpm_change, 1e-12) << "cpm r " << r;
    }
}

TEST(LeidenEngineTest, WellConnectednessAcceptsItsExactBoundary) {
    // Modularity: 1.0 * 2 * (6 - 2) / (2 * 4) = 1.
    const LeidenParams modularity = Params(LeidenObjective::Modularity, 1.0);
    EXPECT_TRUE(leiden_well_connected(1.0, 2.0, 1.0, 6.0, 3.0, 4.0, modularity));
    EXPECT_FALSE(leiden_well_connected(std::nextafter(1.0, 0.0), 2.0, 1.0, 6.0, 3.0,
                                       4.0, modularity));
    // CPM: 0.5 * 2 * (5 - 2) = 3.
    const LeidenParams cpm = Params(LeidenObjective::CPM, 0.5);
    EXPECT_TRUE(leiden_well_connected(3.0, 9.0, 2.0, 9.0, 5.0, 4.0, cpm));
    EXPECT_FALSE(
        leiden_well_connected(std::nextafter(3.0, 0.0), 9.0, 2.0, 9.0, 5.0, 4.0, cpm));
}

TEST(LeidenEngineTest, CandidateSelectionFollowsTheExponentialWeights) {
    const std::vector<double> us = {0.0, 0.25, 0.5, 0.75, 0.999999};
    // Negative gains are never chosen.
    for (const double u : us) {
        EXPECT_EQ(leiden_select_candidate({-1.0, 0.5, -0.2}, 0.01, u), 1u);
        EXPECT_EQ(leiden_select_candidate({-1.0, 0.0}, 0.01, u), 1u);
        EXPECT_EQ(leiden_select_candidate({-1.0, -0.5}, 0.01, u), 2u);
    }
    EXPECT_EQ(leiden_select_candidate({}, 0.01, 0.5), 0u);
    // Gains 0 and 0.01 at theta 0.01 weigh e^-1 and 1, so candidate 0 owns
    // u below e^-1 / (1 + e^-1).
    const double split = std::exp(-1.0) / (1.0 + std::exp(-1.0));
    EXPECT_EQ(leiden_select_candidate({0.0, 0.01}, 0.01, split * 0.999), 0u);
    EXPECT_EQ(leiden_select_candidate({0.0, 0.01}, 0.01, split * 1.001), 1u);
    // A tiny theta always picks the maximum.
    for (const double u : us) {
        EXPECT_EQ(leiden_select_candidate({0.1, 0.2, 0.15}, 1e-6, u), 1u);
    }
}

// ---------------------------------------------------------------------------
// Local moving.

TEST(LeidenEngineTest, LocalMovingBreaksAGainTieTowardTheLowestId) {
    // Node 0, alone in community 0, gains 1 toward community 3 (via node 1)
    // and toward community 2 (via node 2); the tie goes to 2. Every other
    // node is held in place by a weight-4 edge, whatever the visit order.
    const LeidenLevel level =
        Level(5, {{0, 1, 1.0}, {0, 2, 1.0}, {1, 3, 4.0}, {2, 4, 4.0}});
    for (const uint64_t seed : {0u, 1u, 2u, 3u, 4u}) {
        std::vector<uint32_t> community = {0, 3, 2, 3, 2};
        LeidenRng rng(seed);
        leiden_local_move(level, community, Params(LeidenObjective::CPM, 0.0), rng);
        EXPECT_EQ(community, (std::vector<uint32_t>{2, 3, 2, 3, 2})) << "seed " << seed;
    }
}

TEST(LeidenEngineTest, LocalMovingStaysWhenNoGainIsStrictlyBetter) {
    // Node 0 gains 1 by staying with node 1 in community 5 and 1 by joining
    // community 2; equal is not better, so it stays.
    const LeidenLevel level = Level(6, {{0, 1, 1.0}, {0, 2, 1.0}, {2, 3, 4.0}});
    for (const uint64_t seed : {0u, 1u, 2u, 3u, 4u}) {
        std::vector<uint32_t> community = {5, 5, 2, 2, 0, 1};
        LeidenRng rng(seed);
        leiden_local_move(level, community, Params(LeidenObjective::CPM, 0.0), rng);
        EXPECT_EQ(community, (std::vector<uint32_t>{5, 5, 2, 2, 0, 1})) << "seed " << seed;
    }
}

TEST(LeidenEngineTest, LocalMovingTakesTheSmallestEmptyId) {
    // At CPM resolution 2 the pair 0-1 (weight 1) loses 1 together, so the
    // first of them visited leaves for the smallest empty id, 0, and the
    // other, now its community's sole member, keeps 5. Nodes 2 and 3 gain
    // 4 - 2 = 2 together and stay.
    const LeidenLevel level = Level(6, {{0, 1, 1.0}, {2, 3, 4.0}});
    for (const uint64_t seed : {0u, 1u, 2u, 3u, 4u}) {
        std::vector<uint32_t> community = {5, 5, 1, 1, 4, 3};
        LeidenRng rng(seed);
        leiden_local_move(level, community, Params(LeidenObjective::CPM, 2.0), rng);
        EXPECT_EQ((std::set<uint32_t>{community[0], community[1]}),
                  (std::set<uint32_t>{0, 5}))
            << "seed " << seed;
        EXPECT_EQ(community[2], 1u);
        EXPECT_EQ(community[3], 1u);
        EXPECT_EQ(community[4], 4u);
        EXPECT_EQ(community[5], 3u);
    }
}

TEST(LeidenEngineTest, LocalMovingOnATieHeavyFixtureIsPinned) {
    // An 8-cycle of unit edges: every move is a tie between the two sides.
    std::vector<Edge> edges;
    for (uint32_t v = 0; v < 8; ++v) {
        edges.emplace_back(v, (v + 1) % 8, 1.0);
    }
    const LeidenLevel level = Level(8, edges);
    std::vector<uint32_t> community = Iota(8);
    LeidenRng rng(7);
    leiden_local_move(level, community, Params(LeidenObjective::Modularity, 1.0), rng);
    EXPECT_EQ(community, (std::vector<uint32_t>{7, 1, 1, 4, 4, 6, 6, 7}));
}

TEST(LeidenEngineTest, LocalMovingFromSingletonsFindsTheTriangles) {
    const LeidenLevel level = make_base_level(CliquesOfCliques());
    for (const uint64_t seed : {0u, 1u, 2u, 3u, 4u}) {
        std::vector<uint32_t> community = Iota(12);
        LeidenRng rng(seed);
        leiden_local_move(level, community, Params(LeidenObjective::CPM, 0.2), rng);
        for (uint32_t t = 0; t < 4; ++t) {
            EXPECT_EQ(community[3 * t], community[3 * t + 1]) << "seed " << seed;
            EXPECT_EQ(community[3 * t], community[3 * t + 2]) << "seed " << seed;
        }
        EXPECT_EQ((std::set<uint32_t>(community.begin(), community.end())).size(), 4u)
            << "seed " << seed;
    }
}

// ---------------------------------------------------------------------------
// Refinement.

TEST(LeidenEngineTest, RefinementLeavesAPoorlyConnectedNodeAlone) {
    // Node 0 hangs off triangle 1-2-3 by 0.5, below 0.4 * 1 * (4 - 1) = 1.2,
    // so {0} is neither a mover nor a target.
    const LeidenLevel level =
        Level(4, {{0, 1, 0.5}, {1, 2, 4.0}, {1, 3, 4.0}, {2, 3, 4.0}});
    for (const uint64_t seed : {0u, 1u, 2u, 3u, 4u}) {
        LeidenRng rng(seed);
        const std::vector<uint32_t> refined =
            leiden_refine(level, {0, 0, 0, 0}, Params(LeidenObjective::CPM, 0.4), rng);
        EXPECT_EQ(refined[0], 0u) << "seed " << seed;
        for (uint32_t v = 1; v < 4; ++v) {
            EXPECT_NE(refined[v], 0u) << "seed " << seed;
        }
        EXPECT_EQ(refined[1], refined[2]) << "seed " << seed;
        EXPECT_EQ(refined[2], refined[3]) << "seed " << seed;
    }
}

TEST(LeidenEngineTest, RefinementSkipsTargetsThatAreNotWellConnected) {
    // Path 0-1-2 at CPM 0.6: the ends' weight 1 to the rest is below
    // 0.6 * 1 * 2 = 1.2, so node 1 has no target and the ends do not move.
    const LeidenLevel level = Level(3, {{0, 1, 1.0}, {1, 2, 1.0}});
    for (const uint64_t seed : {0u, 1u, 2u, 3u, 4u}) {
        LeidenRng rng(seed);
        EXPECT_EQ(leiden_refine(level, {0, 0, 0}, Params(LeidenObjective::CPM, 0.6), rng),
                  Iota(3))
            << "seed " << seed;
    }
}

TEST(LeidenEngineTest, RefinementNeverJoinsTheDisconnectedPiecesOfACommunity) {
    // Community 0 holds two 4-cliques with no edge between them, as local
    // moving can leave a community after its bridge node moves away; node 8
    // is the departed bridge in community 1.
    std::vector<Edge> edges;
    AddClique(edges, 0, 4);
    AddClique(edges, 4, 4);
    edges.emplace_back(3, 8, 1.0);
    edges.emplace_back(4, 8, 1.0);
    const LeidenLevel level = Level(9, edges);
    const std::vector<uint32_t> community = {0, 0, 0, 0, 0, 0, 0, 0, 1};
    for (const LeidenObjective objective :
         {LeidenObjective::Modularity, LeidenObjective::CPM}) {
        for (const uint64_t seed : {0u, 1u, 2u, 3u, 4u, 5u, 6u, 7u}) {
            LeidenRng rng(seed);
            const std::vector<uint32_t> refined =
                leiden_refine(level, community, Params(objective, 0.05), rng);
            for (uint32_t a = 0; a < 4; ++a) {
                for (uint32_t b = 4; b < 8; ++b) {
                    EXPECT_NE(refined[a], refined[b]) << "seed " << seed;
                }
            }
            EXPECT_EQ(refined[8], 8u) << "seed " << seed;
        }
    }
}

TEST(LeidenEngineTest, RefinedCommunitiesNestAndAreConnected) {
    for (const unsigned seed : {1u, 2u, 3u}) {
        const WeightedGraph graph = RandomSNN(80, seed, 6);
        const LeidenLevel level = make_base_level(graph);
        const LeidenParams params = Params(LeidenObjective::Modularity, 1.0);
        std::vector<uint32_t> community = Iota(80);
        LeidenRng rng(seed);
        leiden_local_move(level, community, params, rng);
        const std::vector<uint32_t> refined = leiden_refine(level, community, params, rng);
        const size_t distinct_refined =
            std::set<uint32_t>(refined.begin(), refined.end()).size();
        EXPECT_LT(distinct_refined, 80u) << "seed " << seed;
        for (size_t v = 0; v < 80; ++v) {
            EXPECT_EQ(community[refined[v]], community[v]);
            EXPECT_TRUE(InducesConnected(graph, refined, refined[v]));
        }
    }
}

// ---------------------------------------------------------------------------
// Aggregation, one pass and the driver.

TEST(LeidenEngineTest, AggregationBuildsTheExactNextLevel) {
    LeidenLevel level =
        Level(4, {{0, 1, 1.0}, {0, 2, 2.0}, {1, 2, 0.5}, {2, 3, 4.0}, {1, 3, 0.25}});
    level.self = {0.5, 0.0, 0.25, 0.0};
    level.size = {1, 2, 1, 3};
    level.strength = {4.0, 1.75, 7.0, 4.25};
    level.m = 8.5;
    const std::vector<uint32_t> refined = {0, 1, 0, 3};
    const std::vector<uint32_t> community = {1, 0, 1, 0};
    const LeidenAggregate aggregate = leiden_aggregate(level, refined, community);
    EXPECT_EQ(aggregate.node_of, (std::vector<uint32_t>{0, 1, 0, 2}));
    EXPECT_EQ(aggregate.community, (std::vector<uint32_t>{0, 1, 1}));
    const LeidenLevel& next = aggregate.level;
    EXPECT_EQ(next.num_nodes, 3u);
    EXPECT_EQ(next.offsets, (std::vector<size_t>{0, 2, 4, 6}));
    EXPECT_EQ(next.neighbors, (std::vector<uint32_t>{1, 2, 0, 2, 0, 1}));
    EXPECT_EQ(next.weights, (std::vector<double>{1.5, 4.0, 1.5, 0.25, 4.0, 0.25}));
    EXPECT_EQ(next.self, (std::vector<double>{2.75, 0.0, 0.0}));
    EXPECT_EQ(next.strength, (std::vector<double>{11.0, 1.75, 4.25}));
    EXPECT_EQ(next.size, (std::vector<size_t>{2, 2, 3}));
    EXPECT_EQ(next.m, 8.5);
    for (const LeidenObjective objective :
         {LeidenObjective::Modularity, LeidenObjective::CPM}) {
        const LeidenParams params = Params(objective, 0.3);
        EXPECT_DOUBLE_EQ(leiden_quality(next, aggregate.community, params),
                         leiden_quality(level, community, params));
    }
}

TEST(LeidenEngineTest, APassProjectsAMultiLevelPartitionToTheBase) {
    const LeidenLevel base = make_base_level(CliquesOfCliques());
    for (const uint64_t seed : {0u, 1u, 2u, 3u, 4u}) {
        LeidenRng rng(seed);
        const std::vector<uint32_t> labels =
            leiden_pass(base, Iota(12), Params(LeidenObjective::CPM, 0.2), rng);
        for (uint32_t v = 1; v < 6; ++v) {
            EXPECT_EQ(labels[v], labels[0]) << "seed " << seed;
            EXPECT_EQ(labels[6 + v], labels[6]) << "seed " << seed;
        }
        EXPECT_NE(labels[0], labels[6]) << "seed " << seed;
    }
}

TEST(LeidenEngineTest, ABridgedPairOfCliquesSplitsUnderBothObjectives) {
    const std::vector<ClusterLabel> expected = Blocks({5, 5});
    for (const uint64_t seed : {0u, 1u, 2u}) {
        EXPECT_EQ(run_leiden(Barbell(), Params(LeidenObjective::Modularity, 1.0), -1,
                             seed)
                      .labels,
                  expected);
        EXPECT_EQ(run_leiden(Barbell(), Params(LeidenObjective::CPM, 0.1), -1, seed)
                      .labels,
                  expected);
    }
}

TEST(LeidenEngineTest, CPMResolutionSpansSingletonsToComponents) {
    std::vector<Edge> edges;
    AddClique(edges, 0, 3);
    AddClique(edges, 3, 3);
    const WeightedGraph graph = Graph(7, edges);
    EXPECT_EQ(run_leiden(graph, Params(LeidenObjective::CPM, 1.5), -1, 0).labels,
              (std::vector<ClusterLabel>{0, 1, 2, 3, 4, 5, 6}));
    EXPECT_EQ(run_leiden(graph, Params(LeidenObjective::CPM, 0.0), -1, 0).labels,
              (std::vector<ClusterLabel>{0, 0, 0, 1, 1, 1, 2}));
}

TEST(LeidenEngineTest, ARingOfCliquesSplitsIntoTheCliques) {
    const std::vector<ClusterLabel> expected = Blocks({5, 5, 5, 5, 5, 5});
    for (const uint64_t seed : {0u, 1u, 2u, 3u, 4u}) {
        const LeidenRun run =
            run_leiden(RingOfCliques(), Params(LeidenObjective::Modularity, 1.0), -1, seed);
        EXPECT_EQ(run.labels, expected) << "seed " << seed;
    }
}

TEST(LeidenEngineTest, TheSameSeedGivesTheSameRun) {
    const WeightedGraph graph = RandomSNN(120, 6, 8);
    const LeidenParams params = Params(LeidenObjective::Modularity, 1.0);
    const LeidenRun first = run_leiden(graph, params, -1, 11);
    const LeidenRun second = run_leiden(graph, params, -1, 11);
    EXPECT_EQ(first.labels, second.labels);
    EXPECT_EQ(first.quality, second.quality);
    EXPECT_EQ(first.iterations, second.iterations);
    for (const uint64_t seed : {0u, 1u, 2u, 3u}) {
        const LeidenRun run = run_leiden(graph, params, -1, seed);
        ASSERT_EQ(run.labels.size(), 120u);
        EXPECT_EQ(run.labels[0], 0);
        ClusterLabel highest = 0;
        for (const ClusterLabel label : run.labels) {
            EXPECT_LE(label, highest + 1);
            highest = std::max(highest, label);
        }
    }
}

TEST(LeidenEngineTest, AFixedIterationCountRunsThatManyPasses) {
    const WeightedGraph graph = RandomSNN(60, 7, 5);
    const LeidenParams params = Params(LeidenObjective::Modularity, 1.0);
    const LeidenRun none = run_leiden(graph, params, 0, 0);
    EXPECT_EQ(none.iterations, 0u);
    std::vector<ClusterLabel> singletons(60);
    std::iota(singletons.begin(), singletons.end(), 0);
    EXPECT_EQ(none.labels, singletons);
    EXPECT_EQ(run_leiden(graph, params, 3, 0).iterations, 3u);
}

TEST(LeidenEngineTest, IteratingToConvergenceStopsAfterTheFirstUnchangedPass) {
    for (const unsigned seed : {1u, 2u, 3u, 4u}) {
        const WeightedGraph graph = RandomSNN(60, seed, 5);
        const LeidenParams params = Params(LeidenObjective::Modularity, 1.0);
        const LeidenRun converged = run_leiden(graph, params, -1, seed);
        ASSERT_GE(converged.iterations, 1u);
        const int64_t passes = static_cast<int64_t>(converged.iterations);
        EXPECT_EQ(run_leiden(graph, params, passes - 1, seed).labels, converged.labels)
            << "seed " << seed;
        if (passes >= 2) {
            EXPECT_NE(run_leiden(graph, params, passes - 2, seed).labels,
                      converged.labels)
                << "seed " << seed;
        }
    }
}

TEST(LeidenEngineTest, AGraphWithoutEdgesGivesSingletonsAndZeroQuality) {
    const WeightedGraph graph = Graph(4, {});
    for (const LeidenObjective objective :
         {LeidenObjective::Modularity, LeidenObjective::CPM}) {
        const LeidenRun run = run_leiden(graph, Params(objective, 1.0), -1, 0);
        EXPECT_EQ(run.labels, (std::vector<ClusterLabel>{0, 1, 2, 3}));
        EXPECT_EQ(run.quality, 0.0);
    }
}

TEST(LeidenEngineTest, ZeroNodesGiveAnEmptyRun) {
    const LeidenRun run =
        run_leiden(Graph(0, {}), Params(LeidenObjective::Modularity, 1.0), -1, 0);
    EXPECT_TRUE(run.labels.empty());
    EXPECT_EQ(run.quality, 0.0);
}

TEST(LeidenEngineTest, EveryClusterIsConnected) {
    // A Fig. 2-style graph: node 0 bridges triangles 1-3 and 4-6 and is
    // pulled toward the 4-clique 7-10.
    std::vector<Edge> edges;
    AddClique(edges, 1, 3);
    AddClique(edges, 4, 3);
    AddClique(edges, 7, 4);
    for (const uint32_t v : {1u, 4u}) {
        edges.emplace_back(0, v, 1.0);
    }
    for (const uint32_t v : {7u, 8u, 9u}) {
        edges.emplace_back(0, v, 1.0);
    }
    edges.emplace_back(3, 6, 0.5);
    const WeightedGraph figure = Graph(11, edges);
    for (const uint64_t seed : {0u, 1u, 2u, 3u, 4u, 5u, 6u, 7u}) {
        for (const LeidenObjective objective :
             {LeidenObjective::Modularity, LeidenObjective::CPM}) {
            ExpectAllClustersConnected(
                figure, run_leiden(figure, Params(objective, 0.3), -1, seed).labels);
        }
    }
    for (const unsigned seed : {1u, 2u, 3u, 4u, 5u}) {
        const WeightedGraph graph = RandomSNN(150, seed, 7);
        for (const LeidenObjective objective :
             {LeidenObjective::Modularity, LeidenObjective::CPM}) {
            const double resolution =
                objective == LeidenObjective::Modularity ? 1.0 : 0.05;
            ExpectAllClustersConnected(
                graph, run_leiden(graph, Params(objective, resolution), -1, seed).labels);
        }
    }
}

// ---------------------------------------------------------------------------
// Public overloads and the result type.

TEST(LeidenTest, TheGraphAndRawInputOverloadsAgree) {
    const size_t n = 40;
    const std::vector<double> condensed = Quantized(n, 9, 7);
    const DenseStorage dense = MakeStorage(n, condensed);
    const auto sparse = MakeSparse(n, condensed, 2.0);
    TableComparison comparison(n, condensed);
    LeidenOptions options = Options(6);
    options.seed = 5;
    const LeidenResult expected = leiden(Oracle(n, condensed, 6), options);
    EXPECT_EQ(expected.K(), 6u);
    for (const size_t threads : {size_t{1}, size_t{4}}) {
        options.num_threads = threads;
        ExpectSameResult(leiden(dense, options), expected);
        ExpectSameResult(leiden(*sparse, options), expected);
        ExpectSameResult(leiden(comparison, options), expected);
    }
}

TEST(LeidenTest, TheResultCarriesTheEngineOutput) {
    const size_t n = 40;
    const std::vector<double> condensed = Quantized(n, 10, 7);
    LeidenOptions options = Options(5);
    options.objective = LeidenObjective::CPM;
    options.resolution = 0.05;
    options.seed = 3;
    const KNNGraph graph = Oracle(n, condensed, 5);
    const LeidenResult result = leiden(graph, options);
    LeidenParams params = Params(LeidenObjective::CPM, 0.05);
    const LeidenRun run = run_leiden(snn_weights(graph, options.prune, 1), params, -1, 3);
    EXPECT_EQ(result.Labels(), run.labels);
    EXPECT_EQ(result.Members(), labels_to_clusters(run.labels));
    EXPECT_EQ(result.Quality(), run.quality);
    EXPECT_EQ(result.Iterations(), run.iterations);
    EXPECT_EQ(result.Objective(), LeidenObjective::CPM);
    EXPECT_EQ(result.Resolution(), 0.05);
    EXPECT_EQ(result.K(), 5u);
    EXPECT_EQ(result.Method(), "leiden");
}

TEST(LeidenTest, InvalidOptionsAreRefusedBeforeAnyComparison) {
    const auto refused = [](LeidenOptions options, const std::string& message) {
        CountingComparison comparison(8, GateFacts());
        ExpectInvalidArgument([&] { leiden(comparison, options); }, message);
        EXPECT_EQ(comparison.Count(), 0u) << message;
        ExpectInvalidArgument([&] { leiden(Oracle(4, Line(4), 2), options); }, message);
    };
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double inf = std::numeric_limits<double>::infinity();
    LeidenOptions forged = Options(3);
    forged.objective = static_cast<LeidenObjective>(7);
    refused(forged, "leiden objective is not a known LeidenObjective");
    for (const double resolution : {-0.5, nan, inf}) {
        LeidenOptions options = Options(3);
        options.resolution = resolution;
        refused(options, "leiden resolution must be finite and non-negative");
    }
    for (const double prune : {-0.1, 1.0, nan, inf}) {
        LeidenOptions options = Options(3);
        options.prune = prune;
        refused(options, "leiden prune must be finite and in [0, 1)");
    }
    for (const double theta : {0.0, -1.0, nan, inf}) {
        LeidenOptions options = Options(3);
        options.theta = theta;
        refused(options, "leiden theta must be finite and positive");
    }
    LeidenOptions iterations = Options(3);
    iterations.n_iterations = -2;
    refused(iterations, "leiden n_iterations must be -1 or non-negative, got -2");
}

TEST(LeidenTest, RawInputBoundsNameLeiden) {
    CountingComparison comparison(4, GateFacts());
    LeidenOptions chunk = Options(1);
    chunk.chunk_size = 0;
    ExpectInvalidArgument([&] { leiden(comparison, chunk); },
                          "leiden chunk_size must be at least one");
    ExpectInvalidArgument([] { leiden(DenseStorage(1), Options(1)); },
                          "leiden needs at least two items: a single item has no "
                          "neighbors");
    ExpectInvalidArgument([&] { leiden(comparison, Options(4)); },
                          "leiden k must be between 1 and 3 for 4 items, got 4");
    ExpectInvalidArgument([&] { leiden(comparison, Options(0)); },
                          "leiden k must be between 1 and 3 for 4 items, got 0");
    EXPECT_EQ(comparison.Count(), 0u);
}

// Each case breaks two adjacent rules of the validation order and pins that
// the earlier rule's message wins.
TEST(LeidenTest, ValidationOrderHoldsWhenSeveralInputsAreInvalid) {
    LeidenOptions both = Options(9);
    both.chunk_size = 0;
    both.theta = 0.0;
    // chunk_size before the options.
    ExpectInvalidArgument([&] { leiden(MakeStorage(4, Line(4)), both); },
                          "leiden chunk_size must be at least one");
    // The options before zero items.
    LeidenOptions theta = Options(9);
    theta.theta = 0.0;
    ExpectInvalidArgument([&] { leiden(DenseStorage(0), theta); },
                          "leiden theta must be finite and positive");
    // The options before k.
    ExpectInvalidArgument([&] { leiden(MakeStorage(4, Line(4)), theta); },
                          "leiden theta must be finite and positive");
    // The options before the comparison facts.
    GateFacts similarity;
    similarity.is_distance = Capability::No;
    CountingComparison counting(4, similarity);
    ExpectInvalidArgument([&] { leiden(counting, theta); },
                          "leiden theta must be finite and positive");
    EXPECT_EQ(counting.Count(), 0u);
}

TEST(LeidenTest, ZeroItemsGiveAnEmptyResultOnEveryOverload) {
    LeidenOptions options = Options(4);
    options.objective = LeidenObjective::CPM;
    options.resolution = 0.25;
    CountingComparison empty(0, GateFacts());
    for (const LeidenResult& result :
         {leiden(KNNGraph(), options), leiden(DenseStorage(0), options),
          leiden(empty, options)}) {
        EXPECT_TRUE(result.Labels().empty());
        EXPECT_TRUE(result.Members().empty());
        EXPECT_EQ(result.Quality(), 0.0);
        EXPECT_EQ(result.Iterations(), 0u);
        EXPECT_EQ(result.Objective(), LeidenObjective::CPM);
        EXPECT_EQ(result.Resolution(), 0.25);
    }
    EXPECT_EQ(leiden(KNNGraph(), options).K(), 0u);
    EXPECT_EQ(leiden(DenseStorage(0), options).K(), 4u);
}

TEST(LeidenTest, TheDefaultResultHasZeroFields) {
    const LeidenResult result;
    EXPECT_TRUE(result.Labels().empty());
    EXPECT_EQ(result.Quality(), 0.0);
    EXPECT_EQ(result.Iterations(), 0u);
    EXPECT_EQ(result.Objective(), LeidenObjective::Modularity);
    EXPECT_EQ(result.Resolution(), 0.0);
    EXPECT_EQ(result.K(), 0u);
    EXPECT_EQ(result.Method(), "leiden");
}

TEST(LeidenTest, TheValueConstructorStoresEveryField) {
    const LeidenResult result({0, 1, 0}, {{0, 2}, {1}}, 0.75, 4, LeidenObjective::CPM,
                              0.2, 7);
    EXPECT_EQ(result.Labels(), (std::vector<ClusterLabel>{0, 1, 0}));
    EXPECT_EQ(result.Members(), (Clusters{{0, 2}, {1}}));
    EXPECT_EQ(result.Quality(), 0.75);
    EXPECT_EQ(result.Iterations(), 4u);
    EXPECT_EQ(result.Objective(), LeidenObjective::CPM);
    EXPECT_EQ(result.Resolution(), 0.2);
    EXPECT_EQ(result.K(), 7u);
}
