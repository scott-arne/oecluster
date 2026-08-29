#include <gtest/gtest.h>
#include "oecluster/CondensedIndex.h"
#include "oecluster/Error.h"

using namespace OECluster;

TEST(CondensedIndexTest, RoundTripsEveryPair) {
    for (size_t n : {2, 3, 5, 7, 12}) {
        size_t expected_index = 0;
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = i + 1; j < n; ++j) {
                EXPECT_EQ(pair_to_condensed(i, j, n), expected_index) << "n=" << n;
                size_t decoded_i = 0;
                size_t decoded_j = 0;
                condensed_to_pair(expected_index, n, decoded_i, decoded_j);
                EXPECT_EQ(decoded_i, i) << "n=" << n;
                EXPECT_EQ(decoded_j, j) << "n=" << n;
                ++expected_index;
            }
        }
        EXPECT_EQ(expected_index, n * (n - 1) / 2) << "n=" << n;
    }
}

TEST(CondensedIndexTest, PairOrderDoesNotMatter) {
    EXPECT_EQ(pair_to_condensed(2, 5, 8), pair_to_condensed(5, 2, 8));
}

TEST(CondensedIndexTest, RejectsOutOfRangeIndex) {
    size_t i = 0;
    size_t j = 0;
    EXPECT_THROW(condensed_to_pair(10, 5, i, j), ComparisonError);
}

TEST(CondensedIndexTest, RejectsIdenticalIndices) {
    EXPECT_THROW(pair_to_condensed(3, 3, 5), ComparisonError);
    EXPECT_THROW(pair_to_condensed(0, 0, 5), ComparisonError);
}

// Regression: before validation, (2, 5) at n=5 returned 9 — the real offset of
// pair (3, 4) — so an out-of-range pair silently aliased a valid slot.
TEST(CondensedIndexTest, RejectsIndexNotBelowN) {
    EXPECT_THROW(pair_to_condensed(2, 5, 5), ComparisonError);
    EXPECT_THROW(pair_to_condensed(5, 2, 5), ComparisonError);
    EXPECT_THROW(pair_to_condensed(0, 7, 5), ComparisonError);
}

TEST(CondensedIndexTest, RejectsEveryIndexAtDegenerateSizes) {
    size_t i = 0;
    size_t j = 0;
    EXPECT_THROW(condensed_to_pair(0, 0, i, j), ComparisonError);
    EXPECT_THROW(condensed_to_pair(0, 1, i, j), ComparisonError);
    EXPECT_THROW(pair_to_condensed(0, 1, 1), ComparisonError);
}
