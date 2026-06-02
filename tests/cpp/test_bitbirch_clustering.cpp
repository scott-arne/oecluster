#include <initializer_list>
#include <stdexcept>
#include <vector>

#include <gtest/gtest.h>

#include "oefp/batch.h"
#include "oefp/fingerprint.h"
#include "oecluster/clustering/BitBirch.h"

namespace {

OEFP::OEFP make_fp(const size_t size_bits, std::initializer_list<size_t> on_bits) {
    OEFP::FingerprintSpec spec;
    spec.size_bits = size_bits;
    spec.value_type = OEFP::FingerprintValueType::Binary;
    spec.source_name = "test";
    OEFP::OEFP fp(spec);
    for (const size_t bit : on_bits) {
        fp.SetBit(bit);
    }
    return fp;
}

OEFP::OEFPBatch make_batch(std::initializer_list<OEFP::OEFP> fps) {
    return OEFP::OEFPBatch::FromFingerprints(std::vector<OEFP::OEFP>(fps));
}

}  // namespace

TEST(BitBirchClusteringTest, ClustersDuplicateBlocksWithDiameterCriterion) {
    const auto batch = make_batch({
        make_fp(4, {0, 1}),
        make_fp(4, {0, 1}),
        make_fp(4, {2, 3}),
        make_fp(4, {2, 3}),
    });
    OECluster::BitBirchOptions options;
    options.threshold = 0.75;
    options.branching_factor = 2;
    options.merge_criterion = OECluster::BitBirchMergeCriterion::Diameter;

    const auto result = OECluster::bitbirch_cluster(batch, options);

    EXPECT_EQ(result.Labels(), (std::vector<OECluster::ClusterLabel>{0, 0, 1, 1}));
    ASSERT_EQ(result.Members().size(), 2u);
    EXPECT_EQ(result.Members()[0], (OECluster::Cluster{0, 1}));
    EXPECT_EQ(result.Members()[1], (OECluster::Cluster{2, 3}));
    EXPECT_EQ(result.ClusterSizes(), (std::vector<size_t>{2, 2}));
    ASSERT_EQ(result.Centroids().Size(), 2u);
    EXPECT_EQ(result.Centroids().PopCount(0), 2u);
    EXPECT_EQ(result.Centroids().PopCount(1), 2u);
}

TEST(BitBirchClusteringTest, PreservesReferenceSplitTieOrdering) {
    const auto batch = make_batch({
        make_fp(4, {0, 1}),
        make_fp(4, {0, 1}),
        make_fp(4, {2, 3}),
        make_fp(4, {2, 3}),
        make_fp(4, {0, 2}),
        make_fp(4, {0, 2}),
    });
    OECluster::BitBirchOptions options;
    options.threshold = 0.75;
    options.branching_factor = 2;
    options.merge_criterion = OECluster::BitBirchMergeCriterion::Diameter;

    const auto result = OECluster::bitbirch_cluster(batch, options);

    ASSERT_EQ(result.Members().size(), 4u);
    EXPECT_EQ(result.Members()[0], (OECluster::Cluster{0, 1}));
    EXPECT_EQ(result.Members()[1], (OECluster::Cluster{2, 3}));
    EXPECT_EQ(result.Members()[2], (OECluster::Cluster{5}));
    EXPECT_EQ(result.Members()[3], (OECluster::Cluster{4}));
    EXPECT_EQ(result.Labels(), (std::vector<OECluster::ClusterLabel>{0, 0, 1, 1, 3, 2}));
}

TEST(BitBirchClusteringTest, RejectsInvalidOptions) {
    const auto batch = make_batch({make_fp(4, {0})});
    OECluster::BitBirchOptions options;

    options.threshold = -0.1;
    EXPECT_THROW(OECluster::bitbirch_cluster(batch, options), std::invalid_argument);

    options.threshold = 0.65;
    options.branching_factor = 0;
    EXPECT_THROW(OECluster::bitbirch_cluster(batch, options), std::invalid_argument);

    options.branching_factor = 1;
    options.tolerance = -0.1;
    EXPECT_THROW(OECluster::bitbirch_cluster(batch, options), std::invalid_argument);
}

TEST(BitBirchClusteringTest, RefinePruneRequiresParentPointers) {
    const auto batch = make_batch({
        make_fp(4, {0, 1}),
        make_fp(4, {0, 1}),
        make_fp(4, {2, 3}),
    });
    OECluster::BitBirchRefinementOptions options;
    options.fit_options.singly = true;
    options.redistribute_largest_cluster = true;

    EXPECT_THROW(OECluster::bitbirch_refine(batch, options), std::invalid_argument);
}

TEST(BitBirchClusteringTest, RefinePruneRedistributesLargestCluster) {
    const auto batch = make_batch({
        make_fp(8, {2, 4}),
        make_fp(8, {1, 4, 6}),
        make_fp(8, {0, 2, 3, 5}),
        make_fp(8, {3, 4, 7}),
        make_fp(8, {4, 7}),
        make_fp(8, {3}),
        make_fp(8, {0, 4, 6, 7}),
        make_fp(8, {4, 5}),
        make_fp(8, {1, 2, 6}),
        make_fp(8, {3}),
        make_fp(8, {1, 3, 5, 7}),
        make_fp(8, {1, 3, 5, 6}),
    });
    OECluster::BitBirchRefinementOptions options;
    options.fit_options.threshold = 0.65;
    options.fit_options.branching_factor = 2;
    options.fit_options.merge_criterion = OECluster::BitBirchMergeCriterion::Diameter;
    options.fit_options.singly = false;
    options.redistribute_largest_cluster = true;

    const auto result = OECluster::bitbirch_refine(batch, options);

    EXPECT_EQ(result.Labels().size(), 12u);
    for (const OECluster::ClusterLabel label : result.Labels()) {
        EXPECT_GE(label, 0);
    }
    size_t assigned = 0;
    for (const auto& cluster : result.Members()) {
        assigned += cluster.size();
    }
    EXPECT_EQ(assigned, 12u);
}

TEST(BitBirchClusteringTest, RefinePruneMatchesZeroSampleReferenceCentroid) {
    const auto batch = make_batch({
        make_fp(3, {0, 1, 2}),
        make_fp(3, {1, 2}),
        make_fp(3, {0}),
        make_fp(3, {0, 1}),
    });
    OECluster::BitBirchRefinementOptions options;
    options.fit_options.threshold = 0.8;
    options.fit_options.branching_factor = 2;
    options.fit_options.merge_criterion = OECluster::BitBirchMergeCriterion::Diameter;
    options.fit_options.singly = false;
    options.redistribute_largest_cluster = true;

    const auto result = OECluster::bitbirch_refine(batch, options);

    EXPECT_EQ(result.Labels(), (std::vector<OECluster::ClusterLabel>{2, 1, 0, 3}));
    ASSERT_EQ(result.Members().size(), 5u);
    EXPECT_EQ(result.Members()[0], (OECluster::Cluster{2}));
    EXPECT_EQ(result.Members()[1], (OECluster::Cluster{1}));
    EXPECT_EQ(result.Members()[2], (OECluster::Cluster{0}));
    EXPECT_EQ(result.Members()[3], (OECluster::Cluster{3}));
    EXPECT_TRUE(result.Members()[4].empty());
}
