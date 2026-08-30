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
