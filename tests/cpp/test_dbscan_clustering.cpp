#include <gtest/gtest.h>

#include <cmath>
#include <limits>
#include <vector>

#include "oecluster/StorageBackend.h"
#include "oecluster/clustering/ClusterTypes.h"
#include "oecluster/clustering/DBSCAN.h"

#include "diversity_test_support.h"

using namespace OECluster;
using namespace diversity_test;

TEST(ClusteringResultTest, ConvertsLabelsToClustersInLabelOrder) {
    const std::vector<ClusterLabel> labels{1, -1, 0, 1, 0, -1};

    const Clusters clusters = labels_to_clusters(labels);

    ASSERT_EQ(clusters.size(), 2);
    EXPECT_EQ(clusters[0], Cluster({2, 4}));
    EXPECT_EQ(clusters[1], Cluster({0, 3}));
}

TEST(ClusteringResultTest, ClusteringResultStoresLabelsAndClusters) {
    const std::vector<ClusterLabel> labels{0, -1, 0};
    const ClusteringResult result(labels, labels_to_clusters(labels));

    EXPECT_EQ(result.Labels(), std::vector<ClusterLabel>({0, -1, 0}));
    ASSERT_EQ(result.Members().size(), 1);
    EXPECT_EQ(result.Members()[0], Cluster({0, 2}));
    EXPECT_EQ(result.NumClusters(), 1u);
    EXPECT_EQ(result.NumSamples(), 3u);
}

TEST(DBSCANClusteringTest, MatchesScikitLearnToyCoreSamples) {
    DenseStorage storage(7);
    const std::vector<double> x{0, 2, 3, 4, 6, 8, 10};
    for (size_t i = 0; i < x.size(); ++i) {
        for (size_t j = i + 1; j < x.size(); ++j) {
            storage.Set(i, j, std::abs(x[i] - x[j]));
        }
    }

    DBSCANOptions options;
    options.eps = 1.0;
    options.min_samples = 3;

    const DBSCANResult result = dbscan_cluster(storage, options);

    EXPECT_EQ(result.CoreSampleIndices(), Cluster({2}));
    EXPECT_EQ(result.Labels(), std::vector<ClusterLabel>({-1, 0, 0, 0, -1, -1, -1}));
    ASSERT_EQ(result.Members().size(), 1);
    EXPECT_EQ(result.Members()[0], Cluster({1, 2, 3}));
}

TEST(DBSCANClusteringTest, AllNoiseWhenNoCoreSamples) {
    DenseStorage storage(3);
    storage.Set(0, 1, 10.0);
    storage.Set(0, 2, 10.0);
    storage.Set(1, 2, 10.0);

    DBSCANOptions options;
    options.eps = 1.0;
    options.min_samples = 2;

    const DBSCANResult result = dbscan_cluster(storage, options);

    EXPECT_TRUE(result.CoreSampleIndices().empty());
    EXPECT_EQ(result.Labels(), std::vector<ClusterLabel>({-1, -1, -1}));
    EXPECT_TRUE(result.Members().empty());
}

TEST(DBSCANClusteringTest, IgnoresAHugeChunkSize) {
    const DenseStorage storage = MakeStorage(5, Positions({0, 1, 3, 7, 8}));
    DBSCANOptions options;
    options.eps = 1.5;
    options.min_samples = 2;

    options.chunk_size = 1;
    options.num_threads = 1;
    const DBSCANResult expected = dbscan_cluster(storage, options);
    ASSERT_EQ(expected.Labels(), (std::vector<ClusterLabel>{0, 0, -1, 1, 1}));

    for (const size_t chunk : {size_t{1}, size_t{4096},
                               std::numeric_limits<size_t>::max()}) {
        SCOPED_TRACE(::testing::Message() << "chunk " << chunk);
        options.chunk_size = chunk;
        for (const size_t threads : {size_t{1}, size_t{4}}) {
            SCOPED_TRACE(::testing::Message() << "threads " << threads);
            options.num_threads = threads;
            const DBSCANResult result = dbscan_cluster(storage, options);
            EXPECT_EQ(result.Labels(), expected.Labels());
        }
    }
}
