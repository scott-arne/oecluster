#include <cmath>
#include <cstdint>
#include <initializer_list>
#include <vector>

#include <gtest/gtest.h>

#include "../../src/clustering/BitBirchKernels.h"

namespace {

using OECluster::detail::BinaryCentroid;
using OECluster::detail::BitBirchLinearSum;
using OECluster::detail::AcceptBitBirchMerge;
using OECluster::detail::JaccardTanimotoISim;
using OECluster::detail::TanimotoWords;
using OECluster::detail::UpdateLinearSumFromWords;

std::vector<uint64_t> words(std::initializer_list<size_t> bits) {
    std::vector<uint64_t> out(1, 0);
    for (const size_t bit : bits) {
        out[0] |= (uint64_t{1} << bit);
    }
    return out;
}

}  // namespace

TEST(BitBirchKernelsTest, UpdatesDenseLinearSumFromWords) {
    BitBirchLinearSum sum(4, 0);
    const auto row = words({0, 2});

    UpdateLinearSumFromWords(sum, row.data(), 4);

    EXPECT_EQ(sum, (BitBirchLinearSum{1, 0, 1, 0}));
}

TEST(BitBirchKernelsTest, IgnoresPaddingBitsWhenUpdatingLinearSum) {
    BitBirchLinearSum sum(65, 0);
    const std::vector<uint64_t> row{
        (uint64_t{1} << 0u) | (uint64_t{1} << 63u),
        (uint64_t{1} << 0u) | (uint64_t{1} << 2u),
    };

    UpdateLinearSumFromWords(sum, row.data(), 65);

    EXPECT_EQ(sum[0], 1u);
    EXPECT_EQ(sum[63], 1u);
    EXPECT_EQ(sum[64], 1u);
    uint32_t total = 0;
    for (const uint32_t count : sum) {
        total += count;
    }
    EXPECT_EQ(total, 3u);
}

TEST(BitBirchKernelsTest, ComputesCentroidWithReferenceHalfTie) {
    const BitBirchLinearSum sum{2, 1, 0, 1};

    const auto centroid = BinaryCentroid(sum, 2);

    EXPECT_EQ(centroid, words({0, 1, 3}));
}

TEST(BitBirchKernelsTest, ComputesZeroSampleCentroidWithReferenceParity) {
    const BitBirchLinearSum sum{0, 0, 0, 0, 0};

    const auto centroid = BinaryCentroid(sum, 0);

    EXPECT_EQ(centroid, words({0, 1, 2, 3, 4}));
}

TEST(BitBirchKernelsTest, ComputesTanimotoFromPackedWords) {
    const auto first = words({0, 1});
    const auto second = words({1, 2});

    const double similarity = TanimotoWords(
        first.data(),
        2,
        second.data(),
        2,
        1);

    EXPECT_DOUBLE_EQ(similarity, 1.0 / 3.0);
}

TEST(BitBirchKernelsTest, MatchesReferenceJaccardTanimotoISim) {
    EXPECT_DOUBLE_EQ(JaccardTanimotoISim(BitBirchLinearSum{2, 0}, 2), 1.0);
    EXPECT_DOUBLE_EQ(JaccardTanimotoISim(BitBirchLinearSum{1, 1}, 2), 0.0);
    EXPECT_TRUE(std::isnan(JaccardTanimotoISim(BitBirchLinearSum{0, 0}, 2)));
}

TEST(BitBirchKernelsTest, RadiusCriterionAcceptsIdenticalSingletons) {
    const BitBirchLinearSum old_sum{1, 1, 0, 0};
    const BitBirchLinearSum nominee_sum{1, 1, 0, 0};
    const BitBirchLinearSum merged_sum{2, 2, 0, 0};
    const std::vector<uint64_t> merged_centroid = BinaryCentroid(merged_sum, 2);

    EXPECT_TRUE(AcceptBitBirchMerge(
        OECluster::BitBirchMergeCriterion::Radius,
        0.75,
        0.05,
        merged_sum,
        merged_centroid,
        2,
        old_sum,
        nominee_sum,
        1,
        1));
}
