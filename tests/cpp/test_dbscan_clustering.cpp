#include <gtest/gtest.h>

#include <cmath>
#include <vector>

#include "oecluster/StorageBackend.h"
#include "oecluster/clustering/ClusterTypes.h"
#include "oecluster/clustering/DBSCAN.h"

using namespace OECluster;

TEST(ClusteringResultTest, ConvertsLabelsToClustersInLabelOrder) {
    const std::vector<ClusterLabel> labels{1, -1, 0, 1, 0, -1};

    const Clusters clusters = labels_to_clusters(labels);

    ASSERT_EQ(clusters.size(), 2);
    EXPECT_EQ(clusters[0], Cluster({2, 4}));
    EXPECT_EQ(clusters[1], Cluster({0, 3}));
}

TEST(ClusteringResultTest, ClusteringResultStoresLabelsAndClusters) {
    ClusteringResult result;
    result.labels = {0, -1, 0};
    result.clusters = labels_to_clusters(result.labels);

    EXPECT_EQ(result.labels, std::vector<ClusterLabel>({0, -1, 0}));
    ASSERT_EQ(result.clusters.size(), 1);
    EXPECT_EQ(result.clusters[0], Cluster({0, 2}));
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

    EXPECT_EQ(result.core_sample_indices, Cluster({2}));
    EXPECT_EQ(result.labels, std::vector<ClusterLabel>({-1, 0, 0, 0, -1, -1, -1}));
    ASSERT_EQ(result.clusters.size(), 1);
    EXPECT_EQ(result.clusters[0], Cluster({1, 2, 3}));
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

    EXPECT_TRUE(result.core_sample_indices.empty());
    EXPECT_EQ(result.labels, std::vector<ClusterLabel>({-1, -1, -1}));
    EXPECT_TRUE(result.clusters.empty());
}
