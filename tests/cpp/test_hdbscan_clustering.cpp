#include <gtest/gtest.h>

#include <cmath>
#include <vector>

#include "oecluster/StorageBackend.h"
#include "oecluster/clustering/HDBSCAN.h"
#include "../../src/clustering/HDBSCANLinkage.h"

using namespace OECluster;

TEST(HDBSCANInfrastructureTest, ComputesCoreDistancesIncludingSelf) {
    DenseStorage storage(4);
    storage.Set(0, 1, 1.0);
    storage.Set(0, 2, 2.0);
    storage.Set(0, 3, 3.0);
    storage.Set(1, 2, 1.5);
    storage.Set(1, 3, 2.5);
    storage.Set(2, 3, 0.5);

    const std::vector<double> core = detail::compute_core_distances(storage, 2, 1);

    EXPECT_DOUBLE_EQ(core[0], 1.0);
    EXPECT_DOUBLE_EQ(core[1], 1.0);
    EXPECT_DOUBLE_EQ(core[2], 0.5);
    EXPECT_DOUBLE_EQ(core[3], 0.5);
}

TEST(HDBSCANInfrastructureTest, BuildsSingleLinkageTreeWithExpectedSize) {
    DenseStorage storage(4);
    storage.Set(0, 1, 0.1);
    storage.Set(0, 2, 2.0);
    storage.Set(0, 3, 2.0);
    storage.Set(1, 2, 2.0);
    storage.Set(1, 3, 2.0);
    storage.Set(2, 3, 0.1);

    const std::vector<double> core = detail::compute_core_distances(storage, 2, 1);
    const std::vector<detail::HDBSCANMSTEdge> mst =
        detail::hdbscan_mutual_reachability_mst(storage, core, 1.0);
    const std::vector<detail::HDBSCANLinkageNode> linkage =
        detail::make_hdbscan_single_linkage(mst, storage.NumSamples());

    ASSERT_EQ(mst.size(), 3);
    ASSERT_EQ(linkage.size(), 3);
    EXPECT_EQ(linkage.back().cluster_size, 4);
}

TEST(HDBSCANClusteringTest, FindsTwoDenseGroups) {
    DenseStorage storage(6);
    const std::vector<double> x{0.0, 0.1, 0.2, 5.0, 5.1, 5.2};
    for (size_t i = 0; i < x.size(); ++i) {
        for (size_t j = i + 1; j < x.size(); ++j) {
            storage.Set(i, j, std::abs(x[i] - x[j]));
        }
    }

    HDBSCANOptions options;
    options.min_cluster_size = 2;
    options.min_samples = 2;
    options.allow_single_cluster = false;
    options.cluster_selection_method = HDBSCANClusterSelectionMethod::Leaf;

    const HDBSCANResult result = hdbscan_cluster(storage, options);

    ASSERT_EQ(result.Labels().size(), 6);
    ASSERT_EQ(result.Members().size(), 2);
    EXPECT_EQ(result.Members()[0], Cluster({0, 1, 2}));
    EXPECT_EQ(result.Members()[1], Cluster({3, 4, 5}));
}

TEST(HDBSCANClusteringTest, RejectsSparseStorageForExactPrecomputedHDBSCAN) {
    SparseStorage storage(3, 0.5);
    storage.Set(0, 1, 0.1);
    storage.Finalize();

    HDBSCANOptions options;
    options.min_cluster_size = 2;

    EXPECT_THROW(hdbscan_cluster(storage, options), std::invalid_argument);
}
