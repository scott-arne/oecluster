#include <gtest/gtest.h>

#include <stdexcept>
#include <vector>

#include "oecluster/StorageBackend.h"
#include "oecluster/clustering/Agglomerative.h"
#include "oecluster/clustering/ClusterTypes.h"

using namespace OECluster;

namespace {

DenseStorage MakeTwoPairStorage() {
    DenseStorage storage(4);
    storage.Set(0, 1, 0.1);
    storage.Set(0, 2, 5.0);
    storage.Set(0, 3, 5.2);
    storage.Set(1, 2, 4.8);
    storage.Set(1, 3, 5.0);
    storage.Set(2, 3, 0.2);
    return storage;
}

DenseStorage MakeUnevenStorage() {
    DenseStorage storage(5);
    storage.Set(0, 1, 0.1);
    storage.Set(0, 2, 0.2);
    storage.Set(0, 3, 10.0);
    storage.Set(0, 4, 11.0);
    storage.Set(1, 2, 0.3);
    storage.Set(1, 3, 12.0);
    storage.Set(1, 4, 13.0);
    storage.Set(2, 3, 20.0);
    storage.Set(2, 4, 21.0);
    storage.Set(3, 4, 0.4);
    return storage;
}

}  // namespace

TEST(AgglomerativeClusteringTest, BuildsAverageLinkageTreeAndClusters) {
    AgglomerativeOptions options;
    options.n_clusters = 2;
    options.linkage = AgglomerativeLinkageMethod::Average;

    const AgglomerativeResult result = agglomerative_cluster(MakeTwoPairStorage(), options);

    EXPECT_EQ(result.Labels(), std::vector<ClusterLabel>({0, 0, 1, 1}));
    ASSERT_EQ(result.Members().size(), 2);
    EXPECT_EQ(result.Members()[0], Cluster({0, 1}));
    EXPECT_EQ(result.Members()[1], Cluster({2, 3}));

    EXPECT_EQ(result.ChildrenLeft(), std::vector<size_t>({0, 2, 4}));
    EXPECT_EQ(result.ChildrenRight(), std::vector<size_t>({1, 3, 5}));
    EXPECT_EQ(result.ClusterSizes(), std::vector<size_t>({2, 2, 4}));
    ASSERT_EQ(result.Distances().size(), 3);
    EXPECT_DOUBLE_EQ(result.Distances()[0], 0.1);
    EXPECT_DOUBLE_EQ(result.Distances()[1], 0.2);
    EXPECT_DOUBLE_EQ(result.Distances()[2], 5.0);
}

TEST(AgglomerativeClusteringTest, HonorsSingleAndCompleteLinkageDistances) {
    AgglomerativeOptions options;
    options.n_clusters = 2;

    options.linkage = AgglomerativeLinkageMethod::Single;
    const AgglomerativeResult single =
        agglomerative_cluster(MakeTwoPairStorage(), options);

    options.linkage = AgglomerativeLinkageMethod::Complete;
    const AgglomerativeResult complete =
        agglomerative_cluster(MakeTwoPairStorage(), options);

    ASSERT_EQ(single.Distances().size(), 3);
    ASSERT_EQ(complete.Distances().size(), 3);
    EXPECT_DOUBLE_EQ(single.Distances().back(), 4.8);
    EXPECT_DOUBLE_EQ(complete.Distances().back(), 5.2);
    EXPECT_EQ(single.Labels(), std::vector<ClusterLabel>({0, 0, 1, 1}));
    EXPECT_EQ(complete.Labels(), std::vector<ClusterLabel>({0, 0, 1, 1}));
}

TEST(AgglomerativeClusteringTest, HonorsWeightedLinkageForUnevenMerges) {
    AgglomerativeOptions options;
    options.n_clusters = 2;

    options.linkage = AgglomerativeLinkageMethod::Average;
    const AgglomerativeResult average =
        agglomerative_cluster(MakeUnevenStorage(), options);

    options.linkage = AgglomerativeLinkageMethod::Weighted;
    const AgglomerativeResult weighted =
        agglomerative_cluster(MakeUnevenStorage(), options);

    EXPECT_EQ(average.Labels(), std::vector<ClusterLabel>({0, 0, 0, 1, 1}));
    EXPECT_EQ(weighted.Labels(), std::vector<ClusterLabel>({0, 0, 0, 1, 1}));
    ASSERT_EQ(average.Distances().size(), 4);
    ASSERT_EQ(weighted.Distances().size(), 4);
    EXPECT_DOUBLE_EQ(average.Distances().back(), 14.5);
    EXPECT_DOUBLE_EQ(weighted.Distances().back(), 16.0);
}

TEST(AgglomerativeClusteringTest, CutsTreeByDistanceThreshold) {
    AgglomerativeOptions options;
    options.linkage = AgglomerativeLinkageMethod::Average;
    options.distance_threshold = 0.15;

    const AgglomerativeResult result = agglomerative_cluster(MakeTwoPairStorage(), options);

    EXPECT_EQ(result.Labels(), std::vector<ClusterLabel>({0, 0, 1, 2}));
    ASSERT_EQ(result.Members().size(), 3);
    EXPECT_EQ(result.Members()[0], Cluster({0, 1}));
    EXPECT_EQ(result.Members()[1], Cluster({2}));
    EXPECT_EQ(result.Members()[2], Cluster({3}));
}

TEST(AgglomerativeClusteringTest, CanStopBeforeFullTreeForClusterCountCut) {
    AgglomerativeOptions options;
    options.n_clusters = 2;
    options.linkage = AgglomerativeLinkageMethod::Average;
    options.compute_full_tree = false;

    const AgglomerativeResult result = agglomerative_cluster(MakeTwoPairStorage(), options);

    EXPECT_EQ(result.Labels(), std::vector<ClusterLabel>({0, 0, 1, 1}));
    EXPECT_EQ(result.ChildrenLeft(), std::vector<size_t>({0, 2}));
    EXPECT_EQ(result.ChildrenRight(), std::vector<size_t>({1, 3}));
    EXPECT_EQ(result.ClusterSizes(), std::vector<size_t>({2, 2}));
    ASSERT_EQ(result.Distances().size(), 2);
    EXPECT_DOUBLE_EQ(result.Distances()[0], 0.1);
    EXPECT_DOUBLE_EQ(result.Distances()[1], 0.2);
}

TEST(AgglomerativeClusteringTest, RejectsInvalidOptionsAndSparseStorage) {
    AgglomerativeOptions options;
    options.n_clusters = 0;
    EXPECT_THROW(agglomerative_cluster(MakeTwoPairStorage(), options), std::invalid_argument);

    options.n_clusters = 5;
    EXPECT_THROW(agglomerative_cluster(MakeTwoPairStorage(), options), std::invalid_argument);

    SparseStorage sparse(3, 0.5);
    sparse.Set(0, 1, 0.1);
    sparse.Finalize();
    options.n_clusters = 2;
    EXPECT_THROW(agglomerative_cluster(sparse, options), std::invalid_argument);
}
