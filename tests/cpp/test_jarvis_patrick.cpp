/**
 * @file test_jarvis_patrick.cpp
 * @brief Jarvis-Patrick clustering over a KNNGraph and over raw input.
 */

#include <gtest/gtest.h>

#include <algorithm>
#include <cstddef>
#include <limits>
#include <map>
#include <set>
#include <string>
#include <vector>

#include "oecluster/GateFacts.h"
#include "oecluster/StorageBackend.h"
#include "oecluster/clustering/JarvisPatrick.h"
#include "oecluster/clustering/KNNGraph.h"

#include "diversity_test_support.h"
#include "knn_test_support.h"

using namespace OECluster;
using namespace diversity_test;
using namespace knn_test;

namespace {

JarvisPatrickOptions Options(size_t k, size_t kmin, size_t num_threads = 0,
                             size_t chunk_size = 4096) {
    JarvisPatrickOptions options;
    options.k = k;
    options.kmin = kmin;
    options.num_threads = num_threads;
    options.chunk_size = chunk_size;
    return options;
}

// Brute-force linkage over the oracle graph: every pair i < j is tested for
// mutual membership and a shared count, with no sorting or binary search,
// and components are found by repeated relabeling.
std::vector<ClusterLabel> OracleLabels(const KNNGraph& graph, size_t kmin) {
    const size_t n = graph.NumItems();
    const size_t k = graph.K();
    std::vector<std::set<size_t>> rows(n);
    for (size_t i = 0; i < n; ++i) {
        rows[i].insert(graph.Indices().begin() + i * k,
                       graph.Indices().begin() + (i + 1) * k);
    }
    std::vector<size_t> component(n);
    for (size_t i = 0; i < n; ++i) {
        component[i] = i;
    }
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i + 1; j < n; ++j) {
            if (rows[i].count(j) == 0 || rows[j].count(i) == 0) {
                continue;
            }
            size_t shared = 0;
            for (const size_t x : rows[i]) {
                shared += rows[j].count(x);
            }
            if (shared < kmin) {
                continue;
            }
            const size_t from = std::max(component[i], component[j]);
            const size_t to = std::min(component[i], component[j]);
            for (size_t& c : component) {
                if (c == from) {
                    c = to;
                }
            }
        }
    }
    // Components are named by their smallest member; number them in order.
    std::map<size_t, ClusterLabel> label_of;
    std::vector<ClusterLabel> labels(n);
    for (size_t i = 0; i < n; ++i) {
        const auto inserted = label_of.emplace(
            component[i], static_cast<ClusterLabel>(label_of.size()));
        labels[i] = inserted.first->second;
    }
    return labels;
}

void ExpectSameResult(const JarvisPatrickResult& actual,
                      const JarvisPatrickResult& expected) {
    EXPECT_EQ(actual.Labels(), expected.Labels());
    EXPECT_EQ(actual.Members(), expected.Members());
    EXPECT_EQ(actual.K(), expected.K());
    EXPECT_EQ(actual.KMin(), expected.KMin());
}

}  // namespace

TEST(JarvisPatrickTest, AHandWorkedExamplePinsLabelsAndMembers) {
    // Items at 0, 1, 2, 4, 7, 8 with k = 2 give the rows
    //   0: {1, 2}  1: {0, 2}  2: {1, 0}  3: {2, 1}  4: {5, 3}  5: {4, 3}.
    // 0, 1 and 2 are pairwise mutual and share one item; 4 and 5 are mutual
    // and share 3. Item 3 names 2 and 1, but neither names 3, so it links to
    // nothing and is a singleton.
    const DenseStorage storage = MakeStorage(6, Positions({0, 1, 2, 4, 7, 8}));
    const JarvisPatrickResult result = jarvis_patrick(storage, Options(2, 1));
    EXPECT_EQ(result.Labels(), (std::vector<ClusterLabel>{0, 0, 0, 1, 2, 2}));
    EXPECT_EQ(result.Members(), (Clusters{{0, 1, 2}, {3}, {4, 5}}));
    EXPECT_EQ(result.K(), 2u);
    EXPECT_EQ(result.KMin(), 1u);
    EXPECT_EQ(result.Method(), "jarvis_patrick");
}

TEST(JarvisPatrickTest, AOneDirectionalNeighborDoesNotLink) {
    // Row 2 names 0, but row 0 does not name 2.
    const KNNGraph graph(3, 1, {1, 0, 0}, {1.0, 1.0, 2.0});
    const JarvisPatrickResult result = jarvis_patrick(graph, 0);
    EXPECT_EQ(result.Members(), (Clusters{{0, 1}, {2}}));
}

TEST(JarvisPatrickTest, KMinMinusOneSharedItemsDoNotLinkAndKMinDo) {
    // Pair (0, 1) is mutual and shares {2, 3}; pairs (0, 2), (0, 3), (1, 2),
    // (1, 3), (2, 4) and (3, 4) are mutual and share exactly one item.
    const KNNGraph graph(5, 3,
                         {1, 2, 3, 0, 2, 3, 0, 1, 4, 0, 1, 4, 2, 3, 0},
                         {1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3});
    EXPECT_EQ(jarvis_patrick(graph, 2).Members(),
              (Clusters{{0, 1}, {2}, {3}, {4}}));
    EXPECT_EQ(jarvis_patrick(graph, 1).Members(), (Clusters{{0, 1, 2, 3, 4}}));
    EXPECT_EQ(jarvis_patrick(graph, 0).Members(), (Clusters{{0, 1, 2, 3, 4}}));
}

TEST(JarvisPatrickTest, MatchesABruteForceLinkageOracle) {
    for (const unsigned seed : {1u, 2u, 3u, 4u, 5u}) {
        const size_t n = 30;
        const std::vector<double> condensed = Quantized(n, seed, 6);
        for (const size_t k : {size_t{2}, size_t{4}, size_t{7}}) {
            const KNNGraph graph = Oracle(n, condensed, k);
            for (size_t kmin = 0; kmin < k; ++kmin) {
                const JarvisPatrickResult result = jarvis_patrick(graph, kmin);
                EXPECT_EQ(result.Labels(), OracleLabels(graph, kmin))
                    << "seed " << seed << " k " << k << " kmin " << kmin;
                EXPECT_EQ(result.Members(), labels_to_clusters(result.Labels()));
            }
        }
    }
}

TEST(JarvisPatrickTest, TheGraphAndRawInputOverloadsAgree) {
    const size_t n = 21;
    const std::vector<double> condensed = Quantized(n, 8, 5);
    const DenseStorage dense = MakeStorage(n, condensed);
    const auto sparse = MakeSparse(n, condensed, 2.0);
    TableComparison comparison(n, condensed);
    const JarvisPatrickResult expected = jarvis_patrick(Oracle(n, condensed, 5), 2);
    for (const size_t threads : {size_t{1}, size_t{4}}) {
        for (const size_t chunk :
             {size_t{1}, size_t{7}, std::numeric_limits<size_t>::max()}) {
            const JarvisPatrickOptions options = Options(5, 2, threads, chunk);
            ExpectSameResult(jarvis_patrick(dense, options), expected);
            ExpectSameResult(jarvis_patrick(*sparse, options), expected);
            ExpectSameResult(jarvis_patrick(comparison, options), expected);
        }
    }
}

TEST(JarvisPatrickTest, KMinAtOrAboveKIsRefusedBeforeAnyComparison) {
    CountingComparison comparison(6, GateFacts());
    ExpectInvalidArgument(
        [&] { jarvis_patrick(comparison, Options(3, 3)); },
        "jarvis_patrick kmin must be less than k = 3, got 3; a mutual pair "
        "shares at most k - 1 neighbors");
    EXPECT_EQ(comparison.Count(), 0u);
    ExpectInvalidArgument(
        [] { jarvis_patrick(MakeStorage(6, Line(6)), Options(2, 5)); },
        "jarvis_patrick kmin must be less than k = 2, got 5; a mutual pair "
        "shares at most k - 1 neighbors");
    ExpectInvalidArgument(
        [] { jarvis_patrick(Oracle(4, Line(4), 2), 2); },
        "jarvis_patrick kmin must be less than k = 2, got 2; a mutual pair "
        "shares at most k - 1 neighbors");
}

TEST(JarvisPatrickTest, ZeroItemsGiveAnEmptyResultOnEveryOverload) {
    for (const size_t kmin : {size_t{0}, size_t{9}}) {
        const JarvisPatrickResult from_graph = jarvis_patrick(KNNGraph(), kmin);
        EXPECT_TRUE(from_graph.Labels().empty());
        EXPECT_TRUE(from_graph.Members().empty());
        EXPECT_EQ(from_graph.KMin(), kmin);
        const JarvisPatrickResult from_storage =
            jarvis_patrick(DenseStorage(0), Options(4, kmin));
        EXPECT_TRUE(from_storage.Labels().empty());
        EXPECT_EQ(from_storage.K(), 4u);
        EXPECT_EQ(from_storage.KMin(), kmin);
        CountingComparison empty(0, GateFacts());
        const JarvisPatrickResult from_comparison =
            jarvis_patrick(empty, Options(4, kmin));
        EXPECT_TRUE(from_comparison.Labels().empty());
        EXPECT_EQ(from_comparison.K(), 4u);
        EXPECT_EQ(from_comparison.KMin(), kmin);
    }
}

TEST(JarvisPatrickTest, RawInputBoundsNameJarvisPatrick) {
    ExpectInvalidArgument(
        [] { jarvis_patrick(DenseStorage(0), Options(1, 0, 0, 0)); },
        "jarvis_patrick chunk_size must be at least one");
    ExpectInvalidArgument(
        [] { jarvis_patrick(DenseStorage(1), Options(1, 0)); },
        "jarvis_patrick needs at least two items: a single item has no neighbors");
    ExpectInvalidArgument(
        [] { jarvis_patrick(MakeStorage(4, Line(4)), Options(4, 0)); },
        "jarvis_patrick k must be between 1 and 3 for 4 items, got 4");
    ExpectInvalidArgument(
        [] { jarvis_patrick(MakeStorage(4, Line(4)), Options(0, 0)); },
        "jarvis_patrick k must be between 1 and 3 for 4 items, got 0");
}

// Each case breaks two adjacent rules of the validation order and pins that
// the earlier rule's message wins.
TEST(JarvisPatrickTest, ValidationOrderHoldsWhenSeveralInputsAreInvalid) {
    const std::string kmin_message =
        "jarvis_patrick kmin must be less than k = 2, got 2; a mutual pair "
        "shares at most k - 1 neighbors";
    // chunk_size before kmin.
    ExpectInvalidArgument(
        [] { jarvis_patrick(MakeStorage(4, Line(4)), Options(2, 5, 0, 0)); },
        "jarvis_patrick chunk_size must be at least one");
    // k before kmin.
    ExpectInvalidArgument(
        [] { jarvis_patrick(MakeStorage(4, Line(4)), Options(0, 5)); },
        "jarvis_patrick k must be between 1 and 3 for 4 items, got 0");
    // kmin before the data-array check.
    ExpectInvalidArgument([] { jarvis_patrick(NullDataStorage(4), Options(2, 2)); },
                          kmin_message);
    // kmin before the sparse short-item scan.
    const auto sparse = MakeSparse(8, Line(8), 0.5);
    ExpectInvalidArgument([&] { jarvis_patrick(*sparse, Options(2, 2)); },
                          kmin_message);
    // kmin before the comparison facts.
    GateFacts similarity;
    similarity.is_distance = Capability::No;
    CountingComparison counting(4, similarity);
    ExpectInvalidArgument([&] { jarvis_patrick(counting, Options(2, 2)); },
                          kmin_message);
    EXPECT_EQ(counting.Count(), 0u);
}

TEST(JarvisPatrickTest, TheDefaultResultHasZeroKAndKMin) {
    const JarvisPatrickResult result;
    EXPECT_EQ(result.K(), 0u);
    EXPECT_EQ(result.KMin(), 0u);
    EXPECT_EQ(result.Method(), "jarvis_patrick");
}
