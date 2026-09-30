/**
 * @file test_butina_golden.cpp
 * @brief Butina outputs pinned before the sphere-exclusion engine replaced
 * the loop that produced them.
 */

#include <gtest/gtest.h>

#include <algorithm>
#include <cstddef>
#include <functional>
#include <utility>
#include <vector>

#include "oecluster/StorageBackend.h"
#include "oecluster/clustering/Butina.h"
#include "../../src/clustering/SphereExclusionEngine.h"
#include "../../src/clustering/ThresholdGraph.h"

#include "diversity_test_support.h"

using namespace OECluster;
using namespace diversity_test;

namespace {

ButinaResult RunButina(const DenseStorage& storage, double threshold,
                       bool reordering, size_t threads = 1,
                       size_t chunk = 4096) {
    ButinaOptions options;
    options.distance_threshold = threshold;
    options.reordering = reordering;
    options.num_threads = threads;
    options.chunk_size = chunk;
    return butina_cluster(storage, options);
}

void ExpectClusters(const ButinaResult& result, const Clusters& expected) {
    EXPECT_EQ(result.Members(), expected);
    size_t n = 0;
    for (const auto& cluster : expected) {
        n += cluster.size();
    }
    ASSERT_EQ(result.Labels().size(), n);
    for (size_t c = 0; c < expected.size(); ++c) {
        for (const size_t member : expected[c]) {
            EXPECT_EQ(result.Labels()[member], static_cast<ClusterLabel>(c));
        }
    }
}

double At(size_t n, const std::vector<double>& condensed, size_t a,
          size_t b) {
    const size_t i = std::min(a, b);
    const size_t j = std::max(a, b);
    return condensed[n * i - i * (i + 1) / 2 + j - i - 1];
}

// A frozen copy of butina_cluster as it stood in 5.9.0, over brute-force
// neighbor lists instead of the threshold graph. The duplication is
// deliberate: this is the oracle that the engine-backed adapter must keep
// matching, so it must not share code with that adapter.
Clusters LegacyButina(size_t n, const std::vector<double>& condensed,
                      double threshold, bool reordering) {
    std::vector<std::vector<size_t>> neighbors(n);
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            if (i == j || At(n, condensed, i, j) <= threshold) {
                neighbors[i].push_back(j);
            }
        }
    }
    using Candidate = std::pair<size_t, size_t>;
    std::vector<Candidate> candidates;
    for (size_t i = 0; i < n; ++i) {
        candidates.emplace_back(neighbors[i].size(), i);
    }
    std::sort(candidates.begin(), candidates.end(), std::greater<Candidate>());
    std::vector<bool> seen(n, false);
    Clusters clusters;
    while (!candidates.empty() && candidates.front().first > 1) {
        const size_t idx = candidates.front().second;
        candidates.erase(candidates.begin());
        if (seen[idx]) {
            continue;
        }
        Cluster cluster{idx};
        seen[idx] = true;
        for (const size_t neighbor : neighbors[idx]) {
            if (!seen[neighbor]) {
                cluster.push_back(neighbor);
                seen[neighbor] = true;
            }
        }
        clusters.push_back(cluster);
        if (reordering) {
            std::vector<bool> affected(n, false);
            for (const size_t member : cluster) {
                for (const size_t neighbor : neighbors[member]) {
                    if (!seen[neighbor]) {
                        affected[neighbor] = true;
                    }
                }
            }
            for (auto& candidate : candidates) {
                if (!affected[candidate.second]) {
                    continue;
                }
                size_t unseen = 0;
                for (const size_t neighbor : neighbors[candidate.second]) {
                    if (!seen[neighbor]) {
                        ++unseen;
                    }
                }
                candidate.first = unseen;
            }
            std::sort(candidates.begin(), candidates.end(),
                      std::greater<Candidate>());
        }
    }
    while (!candidates.empty()) {
        const size_t idx = candidates.front().second;
        candidates.erase(candidates.begin());
        if (!seen[idx]) {
            clusters.push_back(Cluster{idx});
            seen[idx] = true;
        }
    }
    return clusters;
}

}  // namespace

// Expected values were captured from butina_cluster in 5.9.0.
TEST(ButinaGoldenTest, ScrambledSixthsWithoutReordering) {
    const DenseStorage storage = MakeStorage(12, ScrambledSixths(12));
    ExpectClusters(RunButina(storage, 2.0 / 6.0, false),
                   {{11, 1, 2, 7, 8}, {10, 3, 9}, {5}, {4}, {6, 0}});
    ExpectClusters(RunButina(storage, 0.5, false),
                   {{11, 1, 2, 3, 7, 8, 9}, {5}, {10, 4}, {6, 0}});
}

TEST(ButinaGoldenTest, ScrambledSixthsWithReordering) {
    const DenseStorage storage = MakeStorage(12, ScrambledSixths(12));
    ExpectClusters(RunButina(storage, 2.0 / 6.0, true),
                   {{11, 1, 2, 7, 8}, {9, 3, 4, 10}, {6, 0}, {5}});
    ExpectClusters(RunButina(storage, 0.5, true),
                   {{11, 1, 2, 3, 7, 8, 9}, {10, 4}, {6, 0}, {5}});
}

TEST(ButinaGoldenTest, HashedWithoutReordering) {
    const DenseStorage storage = MakeStorage(20, Hashed(20));
    ExpectClusters(RunButina(storage, 0.2, false),
                   {{18, 4, 5, 14, 15, 19}, {13, 0, 10}, {8, 9, 16}, {3},
                    {12, 6, 7, 17}, {2}, {11}, {1}});
    ExpectClusters(RunButina(storage, 0.35, false),
                   {{18, 4, 5, 6, 7, 14, 15, 16, 17, 19}, {13, 0, 2, 10, 12},
                    {8, 9}, {3}, {11, 1}});
}

TEST(ButinaGoldenTest, HashedWithReordering) {
    const DenseStorage storage = MakeStorage(20, Hashed(20));
    ExpectClusters(RunButina(storage, 0.2, true),
                   {{18, 4, 5, 14, 15, 19}, {9, 3, 8, 11, 16}, {17, 2, 7, 12},
                    {10, 0, 1, 13}, {6}});
    ExpectClusters(RunButina(storage, 0.35, true),
                   {{18, 4, 5, 6, 7, 14, 15, 16, 17, 19},
                    {10, 0, 1, 2, 3, 13}, {9, 8, 11, 12}});
}

TEST(ButinaGoldenTest, MatchesTheLegacyLoopOnRandomMatrices) {
    for (const size_t n : {size_t{2}, size_t{3}, size_t{5}, size_t{17},
                           size_t{40}}) {
        for (unsigned seed = 1; seed <= 8; ++seed) {
            const std::vector<double> condensed = Quantized(n, seed, 5);
            const DenseStorage storage = MakeStorage(n, condensed);
            for (const double threshold : {0.0, 0.2, 0.5, 1.0}) {
                for (const bool reordering : {false, true}) {
                    const Clusters expected =
                        LegacyButina(n, condensed, threshold, reordering);
                    for (const size_t threads : {size_t{1}, size_t{4}}) {
                        SCOPED_TRACE(::testing::Message()
                                     << "n " << n << ", seed " << seed
                                     << ", threshold " << threshold
                                     << ", reordering " << reordering
                                     << ", threads " << threads);
                        ExpectClusters(RunButina(storage, threshold,
                                                 reordering, threads, 7),
                                       expected);
                    }
                }
            }
        }
    }
}

TEST(SphereEngineTest, NeighborsFirstReportsButinaClustersAndTheirCenters) {
    const DenseStorage storage = MakeStorage(20, Hashed(20));
    for (const bool reordering : {false, true}) {
        ThresholdGraphOptions graph_options;
        graph_options.threshold = 0.2;
        const ThresholdNeighborGraph graph =
            BuildThresholdNeighborGraph(storage, graph_options);
        const detail::SphereEngineResult result =
            detail::sphere_neighbors_first(graph, reordering);
        const ButinaResult butina = RunButina(storage, 0.2, reordering);

        EXPECT_EQ(result.clusters, butina.Members());
        EXPECT_EQ(result.labels, butina.Labels());
        ASSERT_EQ(result.centers.size(), result.clusters.size());
        for (size_t c = 0; c < result.clusters.size(); ++c) {
            EXPECT_EQ(result.centers[c], result.clusters[c][0]);
        }
    }
}
