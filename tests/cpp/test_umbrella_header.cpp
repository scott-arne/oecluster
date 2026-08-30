/**
 * @file test_umbrella_header.cpp
 * @brief The umbrella header must reach every public surface it advertises.
 */

#include <gtest/gtest.h>
#include <oecluster/oecluster.h>

// No other OECluster include: naming these types is the assertion. A header
// missing from the umbrella fails this at compile time, not at run time.
TEST(UmbrellaHeaderTest, ReachesDescriptorStatistics) {
    const OECluster::DescriptorStatisticsOptions options;
    EXPECT_TRUE(options.sources.empty());
    EXPECT_FALSE(options.inverse_covariance);

    const OECluster::DescriptorStatisticsResult result;
    EXPECT_EQ(result.num_rows, 0u);
    EXPECT_TRUE(result.columns.empty());
}

TEST(UmbrellaHeaderTest, ReachesDescriptorComparison) {
    const OECluster::DescriptorOptions options;
    EXPECT_TRUE(options.sources.empty());
    EXPECT_EQ(options.metric, "standardized_euclidean");
    EXPECT_EQ(options.missing, "complete_case");

    // The type name is the assertion; cannot construct without molecules.
    const OECluster::DescriptorComparison* ptr = nullptr;
    EXPECT_EQ(ptr, nullptr);
}

TEST(UmbrellaHeaderTest, ReachesRMSDComparison) {
    const OECluster::RMSDOptions options;
    EXPECT_FALSE(options.overlay);
    EXPECT_TRUE(options.automorph);
    EXPECT_TRUE(options.heavy_only);

    // The type name is the assertion; cannot construct without molecules.
    const OECluster::RMSDComparison* ptr = nullptr;
    EXPECT_EQ(ptr, nullptr);
}
