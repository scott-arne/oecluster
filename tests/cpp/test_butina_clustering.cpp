#include <gtest/gtest.h>

#include <stdexcept>
#include <tuple>
#include <vector>

#include "oecluster/StorageBackend.h"
#include "oecluster/clustering/Butina.h"
#include "oecluster/clustering/Representative.h"
#include "../../src/clustering/ThresholdGraph.h"

using namespace OECluster;

namespace {

std::vector<size_t> neighbor_vector(const ThresholdNeighborGraph& graph, size_t index) {
    const auto neighbors = graph.Neighbors(index);
    return std::vector<size_t>(neighbors.begin(), neighbors.end());
}

}  // namespace

TEST(ButinaClusteringTest, SingletonInputReturnsSingletonCluster) {
    DenseStorage storage(1);
    ButinaOptions options;
    options.distance_threshold = 0.5;

    const ClusteringResult result = butina_cluster(storage, options);

    ASSERT_EQ(result.clusters.size(), 1);
    ASSERT_EQ(result.clusters[0], Cluster({0}));
    ASSERT_EQ(result.labels, std::vector<ClusterLabel>({0}));
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
    EXPECT_EQ(neighbor_vector(graph, 0), std::vector<size_t>({0, 1}));
    EXPECT_EQ(neighbor_vector(graph, 1), std::vector<size_t>({0, 1, 2}));
    EXPECT_EQ(neighbor_vector(graph, 2), std::vector<size_t>({1, 2, 3}));
    EXPECT_EQ(neighbor_vector(graph, 3), std::vector<size_t>({2, 3}));
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
    EXPECT_EQ(neighbor_vector(singleton_graph, 0), std::vector<size_t>({0}));
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
        EXPECT_EQ(neighbor_vector(sparse_graph, i), neighbor_vector(dense_graph, i));
    }
}

TEST(ThresholdGraphTest, DenseParallelBuildIsDeterministic) {
    DenseStorage storage(5);
    storage.Set(0, 1, 0.10);
    storage.Set(0, 2, 0.10);
    storage.Set(0, 3, 0.50);
    storage.Set(0, 4, 0.90);
    storage.Set(1, 2, 0.20);
    storage.Set(1, 3, 0.10);
    storage.Set(1, 4, 0.90);
    storage.Set(2, 3, 0.10);
    storage.Set(2, 4, 0.90);
    storage.Set(3, 4, 0.10);

    ThresholdGraphOptions serial_options;
    serial_options.threshold = 0.2;
    serial_options.num_threads = 1;

    ThresholdGraphOptions parallel_options = serial_options;
    parallel_options.num_threads = 4;

    const auto serial_graph = BuildThresholdNeighborGraph(storage, serial_options);
    const auto parallel_graph = BuildThresholdNeighborGraph(storage, parallel_options);

    ASSERT_EQ(serial_graph.Size(), parallel_graph.Size());
    for (size_t i = 0; i < serial_graph.Size(); ++i) {
        EXPECT_EQ(neighbor_vector(serial_graph, i), neighbor_vector(parallel_graph, i));
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

    const ClusteringResult result = butina_cluster(storage, options);

    ASSERT_EQ(result.clusters.size(), 2);
    EXPECT_EQ(result.clusters[0], Cluster({3, 2}));
    EXPECT_EQ(result.clusters[1], Cluster({1, 0}));

    // Labels view is consistent with cluster order: label == cluster position.
    ASSERT_EQ(result.labels.size(), storage.NumItems());
    for (size_t c = 0; c < result.clusters.size(); ++c) {
        for (const size_t member : result.clusters[c]) {
            EXPECT_EQ(result.labels[member], static_cast<ClusterLabel>(c));
        }
    }
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

    const ClusteringResult result = butina_cluster(storage, options);

    ASSERT_EQ(result.clusters.size(), 2);
    EXPECT_EQ(result.clusters[0], Cluster({0, 1, 2, 3}));
    EXPECT_EQ(result.clusters[1], Cluster({4}));

    // Every item is assigned a non-noise label equal to its cluster position.
    ASSERT_EQ(result.labels.size(), storage.NumItems());
    EXPECT_EQ(result.labels, std::vector<ClusterLabel>({0, 0, 0, 0, 1}));
}

TEST(ButinaClusteringTest, ReorderingCanChangeLaterRepresentatives) {
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

    const ClusteringResult result_no_reorder = butina_cluster(storage, no_reorder);
    const ClusteringResult result_reorder = butina_cluster(storage, reorder);

    EXPECT_NE(result_no_reorder.clusters, result_reorder.clusters);
}

TEST(ButinaClusteringTest, NegativeThresholdThrows) {
    DenseStorage storage(2);
    ButinaOptions options;
    options.distance_threshold = -0.1;

    EXPECT_THROW(butina_cluster(storage, options), std::invalid_argument);
}

TEST(ClusterRepresentativeTest, MedoidMinimizesMeanDistanceWithClusterOrderTieBreak) {
    DenseStorage storage(3);
    storage.Set(0, 1, 1.0);
    storage.Set(0, 2, 1.0);
    storage.Set(1, 2, 0.1);
    const Cluster cluster{0, 1, 2};

    EXPECT_EQ(cluster_representative(cluster, storage, RepresentativeMethod::Medoid), 1);
}

TEST(ClusterRepresentativeTest, MinimaxMinimizesWorstNeighborDistance) {
    DenseStorage storage(4);
    storage.Set(0, 1, 0.1);
    storage.Set(0, 2, 0.1);
    storage.Set(0, 3, 1.0);
    storage.Set(1, 2, 0.6);
    storage.Set(1, 3, 0.6);
    storage.Set(2, 3, 0.6);
    const Cluster cluster{0, 1, 2, 3};

    EXPECT_EQ(cluster_representative(cluster, storage, RepresentativeMethod::Medoid), 0);
    EXPECT_EQ(cluster_representative(cluster, storage, RepresentativeMethod::Minimax), 1);
}

TEST(ClusterRepresentativeTest, HighestNeighborhoodUsesThresholdNeighborCounts) {
    DenseStorage storage(4);
    storage.Set(0, 1, 0.1);
    storage.Set(0, 2, 0.1);
    storage.Set(0, 3, 1.0);
    storage.Set(1, 2, 0.6);
    storage.Set(1, 3, 0.6);
    storage.Set(2, 3, 0.6);
    const Cluster cluster{0, 1, 2, 3};

    RepresentativeOptions options;
    options.method = RepresentativeMethod::HighestNeighborhood;
    options.neighbor_threshold = 0.2;

    EXPECT_EQ(cluster_representative(cluster, storage, options), 0);
}

TEST(ClusterRepresentativeTest, RankRepresentativesReportsQualityMetrics) {
    DenseStorage storage(5);
    storage.Set(0, 1, 0.1);
    storage.Set(0, 2, 0.1);
    storage.Set(0, 3, 1.0);
    storage.Set(0, 4, 0.9);
    storage.Set(1, 2, 0.6);
    storage.Set(1, 3, 0.6);
    storage.Set(1, 4, 0.7);
    storage.Set(2, 3, 0.6);
    storage.Set(2, 4, 0.8);
    storage.Set(3, 4, 0.2);
    const Cluster cluster{0, 1, 2, 3};

    RepresentativeOptions options;
    options.method = RepresentativeMethod::Medoid;
    options.neighbor_threshold = 0.2;
    options.scaffold_labels = {"A", "A", "B", "B", "external"};

    const std::vector<ClusterRepresentative> ranked =
        rank_representatives(cluster, storage, options);

    ASSERT_EQ(ranked.size(), 4u);
    EXPECT_EQ(ranked[0].member, 0u);
    EXPECT_DOUBLE_EQ(ranked[0].score, 0.4);
    EXPECT_DOUBLE_EQ(ranked[0].metrics.mean_distance_to_cluster, 0.4);
    EXPECT_DOUBLE_EQ(ranked[0].metrics.max_distance_to_cluster, 1.0);
    EXPECT_DOUBLE_EQ(ranked[0].metrics.median_distance_to_cluster, 0.1);
    EXPECT_DOUBLE_EQ(ranked[0].metrics.neighbor_fraction_at_threshold, 0.75);
    EXPECT_DOUBLE_EQ(ranked[0].metrics.nearest_external_distance, 0.9);
    EXPECT_DOUBLE_EQ(ranked[0].metrics.cluster_radius, 1.0);
    EXPECT_DOUBLE_EQ(ranked[0].metrics.cluster_diameter, 1.0);
    EXPECT_NEAR(ranked[0].metrics.silhouette_like_score, (0.9 - 0.4) / 0.9, 1e-12);
    EXPECT_DOUBLE_EQ(ranked[0].metrics.scaffold_purity, 0.5);
    EXPECT_EQ(ranked[0].metrics.representative_rank, 1u);
}

TEST(ClusterRepresentativeTest, WeightedMedoidUsesPenaltyAndPriorityVectors) {
    DenseStorage storage(4);
    storage.Set(0, 1, 0.1);
    storage.Set(0, 2, 0.1);
    storage.Set(0, 3, 1.0);
    storage.Set(1, 2, 0.6);
    storage.Set(1, 3, 0.6);
    storage.Set(2, 3, 0.6);
    const Cluster cluster{0, 1, 2, 3};

    RepresentativeOptions options;
    options.method = RepresentativeMethod::WeightedMedoid;
    options.weights.alpha = 1.0;
    options.weights.beta = 1.0;
    options.weights.gamma = 1.0;
    options.liability_penalties = {0.0, 0.0, 0.0, 0.0};
    options.priority_scores = {0.0, 0.4, 0.0, 0.0};

    const std::vector<ClusterRepresentative> ranked =
        rank_representatives(cluster, storage, options);

    ASSERT_EQ(ranked.size(), 4u);
    EXPECT_EQ(ranked[0].member, 1u);
    EXPECT_NEAR(ranked[0].score, ((0.1 + 0.6 + 0.6) / 3.0) - 0.4, 1e-12);
}

TEST(ClusterRepresentativeTest, SelectRepresentativesSupportsScoreAndDiversityModes) {
    DenseStorage storage(4);
    storage.Set(0, 1, 0.1);
    storage.Set(0, 2, 0.8);
    storage.Set(0, 3, 0.8);
    storage.Set(1, 2, 0.8);
    storage.Set(1, 3, 0.8);
    storage.Set(2, 3, 0.2);
    const Cluster cluster{0, 1, 2, 3};

    RepresentativeOptions score_options;
    score_options.method = RepresentativeMethod::Medoid;
    score_options.selection = RepresentativeSelection::Score;

    RepresentativeOptions diversity_options = score_options;
    diversity_options.selection = RepresentativeSelection::Diversity;

    const std::vector<ClusterRepresentative> score_selected =
        select_representatives(cluster, storage, 2, score_options);
    const std::vector<ClusterRepresentative> diversity_selected =
        select_representatives(cluster, storage, 2, diversity_options);

    ASSERT_EQ(score_selected.size(), 2u);
    ASSERT_EQ(diversity_selected.size(), 2u);
    EXPECT_EQ(score_selected[0].member, 0u);
    EXPECT_EQ(score_selected[1].member, 1u);
    EXPECT_EQ(diversity_selected[0].member, 0u);
    EXPECT_EQ(diversity_selected[1].member, 2u);
    EXPECT_EQ(diversity_selected[1].metrics.representative_rank, 3u);
}

TEST(ClusterRepresentativeTest, EmptyClusterThrows) {
    DenseStorage storage(1);

    EXPECT_THROW(
        cluster_representative(Cluster{}, storage, RepresentativeMethod::Medoid),
        std::invalid_argument);
}

TEST(ClusterRepresentativeTest, InvalidMemberThrows) {
    DenseStorage storage(2);

    EXPECT_THROW(
        cluster_representative(Cluster({0, 2}), storage, RepresentativeMethod::Medoid),
        std::out_of_range);
}

TEST(ClusterRepresentativeTest, SparseStorageRejectsDistanceBasedMethods) {
    SparseStorage storage(3, 0.2);
    storage.Set(0, 1, 0.1);
    storage.Set(0, 2, 0.3);
    storage.Set(1, 2, 0.3);
    storage.Finalize();
    const Cluster cluster{1, 0, 2};

    EXPECT_EQ(cluster_representative(Cluster({2}), storage, RepresentativeMethod::Medoid), 2);
    EXPECT_THROW(
        cluster_representative(cluster, storage, RepresentativeMethod::Medoid),
        std::invalid_argument);
}
