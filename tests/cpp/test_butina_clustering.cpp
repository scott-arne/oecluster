#include <gtest/gtest.h>

#include <stdexcept>
#include <tuple>
#include <vector>

#include "oecluster/StorageBackend.h"
#include "oecluster/clustering/Butina.h"
#include "../../src/clustering/ThresholdGraph.h"

using namespace OECluster;

TEST(ButinaClusteringTest, SingletonInputReturnsSingletonCluster) {
    DenseStorage storage(1);
    ButinaOptions options;
    options.distance_threshold = 0.5;

    const Clusters clusters = butina_cluster(storage, options);

    ASSERT_EQ(clusters.size(), 1);
    ASSERT_EQ(clusters[0], Cluster({0}));
}

TEST(ThresholdGraphTest, DenseStorageIncludesSelfAndThresholdNeighbors) {
    DenseStorage storage(4);
    storage.Set(0, 1, 0.10);
    storage.Set(0, 2, 0.40);
    storage.Set(0, 3, 0.90);
    storage.Set(1, 2, 0.20);
    storage.Set(1, 3, 0.80);
    storage.Set(2, 3, 0.30);

    ThresholdGraphOptions options;
    options.threshold = 0.35;
    options.num_threads = 1;
    const ThresholdNeighborGraph graph = BuildThresholdNeighborGraph(storage, options);

    ASSERT_EQ(graph.Size(), 4);
    EXPECT_EQ(graph.Neighbors(0), std::vector<size_t>({0, 1}));
    EXPECT_EQ(graph.Neighbors(1), std::vector<size_t>({0, 1, 2}));
    EXPECT_EQ(graph.Neighbors(2), std::vector<size_t>({1, 2, 3}));
    EXPECT_EQ(graph.Neighbors(3), std::vector<size_t>({2, 3}));
}

TEST(ThresholdGraphTest, EmptyAndSingletonDenseStorageDoNotRequireData) {
    ThresholdGraphOptions options;
    options.threshold = 0.35;

    DenseStorage empty(0);
    const ThresholdNeighborGraph empty_graph = BuildThresholdNeighborGraph(empty, options);
    EXPECT_EQ(empty_graph.Size(), 0);

    DenseStorage singleton(1);
    const ThresholdNeighborGraph singleton_graph =
        BuildThresholdNeighborGraph(singleton, options);
    ASSERT_EQ(singleton_graph.Size(), 1);
    EXPECT_EQ(singleton_graph.Neighbors(0), std::vector<size_t>({0}));
}

TEST(ThresholdGraphTest, SparseStorageMatchesDenseStorage) {
    DenseStorage dense(4);
    SparseStorage sparse(4, 0.35);
    for (const auto& entry : std::vector<std::tuple<size_t, size_t, double>>{
             {0, 1, 0.10}, {0, 2, 0.40}, {0, 3, 0.90},
             {1, 2, 0.20}, {1, 3, 0.80}, {2, 3, 0.30}}) {
        const auto [i, j, value] = entry;
        dense.Set(i, j, value);
        sparse.Set(i, j, value);
    }
    sparse.Finalize();

    ThresholdGraphOptions options;
    options.threshold = 0.35;
    const ThresholdNeighborGraph dense_graph = BuildThresholdNeighborGraph(dense, options);
    const ThresholdNeighborGraph sparse_graph = BuildThresholdNeighborGraph(sparse, options);

    ASSERT_EQ(sparse_graph.Size(), dense_graph.Size());
    for (size_t i = 0; i < dense_graph.Size(); ++i) {
        EXPECT_EQ(sparse_graph.Neighbors(i), dense_graph.Neighbors(i));
    }
}

TEST(ThresholdGraphTest, SparseStorageRejectsThresholdAboveSparseCutoff) {
    SparseStorage sparse(3, 0.20);
    sparse.Set(0, 1, 0.10);
    sparse.Set(0, 2, 0.30);
    sparse.Set(1, 2, 0.30);
    sparse.Finalize();

    ThresholdGraphOptions options;
    options.threshold = 0.30;

    EXPECT_THROW(BuildThresholdNeighborGraph(sparse, options), std::invalid_argument);
}

TEST(ButinaClusteringTest, ChoosesLargestIndexOnNeighborCountTie) {
    DenseStorage storage(4);
    storage.Set(0, 1, 0.10);
    storage.Set(0, 2, 0.90);
    storage.Set(0, 3, 0.90);
    storage.Set(1, 2, 0.90);
    storage.Set(1, 3, 0.90);
    storage.Set(2, 3, 0.10);

    ButinaOptions options;
    options.distance_threshold = 0.2;

    const Clusters clusters = butina_cluster(storage, options);

    ASSERT_EQ(clusters.size(), 2);
    EXPECT_EQ(clusters[0], Cluster({3, 2}));
    EXPECT_EQ(clusters[1], Cluster({1, 0}));
}

TEST(ButinaClusteringTest, ClustersUnseenNeighborsInAscendingNeighborOrder) {
    DenseStorage storage(5);
    storage.Set(0, 1, 0.10);
    storage.Set(0, 2, 0.10);
    storage.Set(0, 3, 0.10);
    storage.Set(0, 4, 0.90);
    storage.Set(1, 2, 0.90);
    storage.Set(1, 3, 0.90);
    storage.Set(1, 4, 0.90);
    storage.Set(2, 3, 0.90);
    storage.Set(2, 4, 0.90);
    storage.Set(3, 4, 0.90);

    ButinaOptions options;
    options.distance_threshold = 0.2;

    const Clusters clusters = butina_cluster(storage, options);

    ASSERT_EQ(clusters.size(), 2);
    EXPECT_EQ(clusters[0], Cluster({0, 1, 2, 3}));
    EXPECT_EQ(clusters[1], Cluster({4}));
}

TEST(ButinaClusteringTest, ReorderingCanChangeLaterCenters) {
    DenseStorage storage(6);
    storage.Set(0, 1, 0.10);
    storage.Set(0, 2, 0.10);
    storage.Set(0, 3, 0.90);
    storage.Set(0, 4, 0.90);
    storage.Set(0, 5, 0.90);
    storage.Set(1, 2, 0.90);
    storage.Set(1, 3, 0.10);
    storage.Set(1, 4, 0.90);
    storage.Set(1, 5, 0.90);
    storage.Set(2, 3, 0.90);
    storage.Set(2, 4, 0.10);
    storage.Set(2, 5, 0.90);
    storage.Set(3, 4, 0.90);
    storage.Set(3, 5, 0.10);
    storage.Set(4, 5, 0.10);

    ButinaOptions no_reorder;
    no_reorder.distance_threshold = 0.2;
    no_reorder.reordering = false;

    ButinaOptions reorder = no_reorder;
    reorder.reordering = true;

    const Clusters clusters_no_reorder = butina_cluster(storage, no_reorder);
    const Clusters clusters_reorder = butina_cluster(storage, reorder);

    EXPECT_NE(clusters_no_reorder, clusters_reorder);
}

TEST(ButinaClusteringTest, NegativeThresholdThrows) {
    DenseStorage storage(2);
    ButinaOptions options;
    options.distance_threshold = -0.1;

    EXPECT_THROW(butina_cluster(storage, options), std::invalid_argument);
}
