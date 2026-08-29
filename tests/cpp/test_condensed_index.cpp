#include <gtest/gtest.h>
#include "oecluster/CondensedIndex.h"
#include "oecluster/Error.h"

using namespace OECluster;

TEST(CondensedIndexTest, RoundTripsEveryPair) {
    const size_t n = 7;
    size_t expected_index = 0;
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i + 1; j < n; ++j) {
            EXPECT_EQ(pair_to_condensed(i, j, n), expected_index);
            size_t decoded_i = 0;
            size_t decoded_j = 0;
            condensed_to_pair(expected_index, n, decoded_i, decoded_j);
            EXPECT_EQ(decoded_i, i);
            EXPECT_EQ(decoded_j, j);
            ++expected_index;
        }
    }
    EXPECT_EQ(expected_index, n * (n - 1) / 2);
}

TEST(CondensedIndexTest, PairOrderDoesNotMatter) {
    EXPECT_EQ(pair_to_condensed(2, 5, 8), pair_to_condensed(5, 2, 8));
}

TEST(CondensedIndexTest, RejectsOutOfRangeIndex) {
    size_t i = 0;
    size_t j = 0;
    EXPECT_THROW(condensed_to_pair(10, 5, i, j), ComparisonError);
}
