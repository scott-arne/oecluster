#include <gtest/gtest.h>

#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <vector>

#include "isim_test_support.h"
#include "oecluster/clustering/ISimReport.h"
#include "../../src/clustering/BitBirchKernels.h"
#include "../../src/clustering/ISimKernels.h"

using namespace OECluster;
using namespace OECluster::detail;

namespace {
constexpr uint64_t MAX64 = std::numeric_limits<uint64_t>::max();
}  // namespace

TEST(ISimKernelsTest, AddCarriesAcrossTheLowWord) {
    UInt128 value = UInt128::FromU64(MAX64);
    value += UInt128::FromU64(1);
    EXPECT_EQ(value.high, 1u);
    EXPECT_EQ(value.low, 0u);
}

TEST(ISimKernelsTest, SubtractBorrowsAcrossTheLowWord) {
    UInt128 value{1u, 0u};
    value -= UInt128::FromU64(1);
    EXPECT_EQ(value.high, 0u);
    EXPECT_EQ(value.low, MAX64);
}

TEST(ISimKernelsTest, MultiplyAtTheExtremes) {
    const UInt128 product = multiply_u64(MAX64, MAX64);
    // (2^64 - 1)^2 = 2^128 - 2^65 + 1.
    EXPECT_EQ(product.high, MAX64 - 1u);
    EXPECT_EQ(product.low, 1u);
    const UInt128 shifted = multiply_u64(uint64_t{1} << 32, uint64_t{1} << 32);
    EXPECT_EQ(shifted.high, 1u);
    EXPECT_EQ(shifted.low, 0u);
    EXPECT_EQ(multiply_u64(0u, MAX64), UInt128{});
}

TEST(ISimKernelsTest, OrdersAcrossTheWordBoundary) {
    EXPECT_TRUE(UInt128::FromU64(MAX64) < (UInt128{1u, 0u}));
    EXPECT_FALSE((UInt128{1u, 0u}) < UInt128::FromU64(MAX64));
    EXPECT_TRUE((UInt128{1u, 5u}) < (UInt128{2u, 0u}));
    EXPECT_FALSE((UInt128{2u, 0u}) < (UInt128{2u, 0u}));
}

TEST(ISimKernelsTest, HalvesAcrossTheWordBoundary) {
    const UInt128 half = UInt128{1u, 0u}.Half();
    EXPECT_EQ(half.high, 0u);
    EXPECT_EQ(half.low, uint64_t{1} << 63);
}

TEST(ISimKernelsTest, ConvertsAbove2To64) {
    EXPECT_EQ((UInt128{1u, 0u}).ToDouble(), 18446744073709551616.0);
    EXPECT_EQ((UInt128{3u, 0u}).ToDouble(), 3.0 * 18446744073709551616.0);
    // 2^64 + 2^12 is exactly representable; the low word must not be lost.
    EXPECT_EQ((UInt128{1u, 4096u}).ToDouble(), 18446744073709551616.0 + 4096.0);
}

TEST(ISimKernelsTest, ToDoubleAvoidsDoubleRounding) {
    // high*2^64 is exact and ldexp is exact, but casting low to double on its
    // own would round low first whenever low >= 2^53; adding the (already
    // rounded) high contribution afterwards is a second rounding. The
    // correctly-rounded value of 2^64 + 18446744073709545473 is
    // 0x1.fffffffffffffp+64; the naive high*2^64 + (double)low formula
    // instead lands one ULP low, at 0x1.ffffffffffffep+64.
    EXPECT_EQ((UInt128{1u, 18446744073709545473ull}).ToDouble(), 0x1.fffffffffffffp+64);
    // An exact tie between 2^64 and the next double up, with every bit below
    // the rounding boundary clear: ties-to-even must land on 2^64, whose
    // trailing mantissa bit is 0.
    EXPECT_EQ((UInt128{1u, 2048ull}).ToDouble(), 0x1.0000000000000p+64);
    // The same tie, but with one more low bit set below the boundary: that
    // bit is a sticky bit the rounding must see, breaking the tie upward
    // instead of to even.
    EXPECT_EQ((UInt128{1u, 2049ull}).ToDouble(), 0x1.0000000000001p+64);
}

TEST(ISimKernelsTest, BatchSizeLimitIsTwoToThe32) {
    EXPECT_NO_THROW(check_isim_batch_size(static_cast<size_t>(4294967295ull)));
    EXPECT_THROW(check_isim_batch_size(static_cast<size_t>(4294967296ull)),
                 std::invalid_argument);
}

TEST(ISimKernelsTest, WidthLimitIsTwoToThe31) {
    EXPECT_NO_THROW(check_isim_width(2147483647ull));
    EXPECT_THROW(check_isim_width(2147483648ull), std::invalid_argument);
}

TEST(ISimKernelsTest, RatioTreatsAZeroUnionAsSimilarityOne) {
    EXPECT_EQ(isim_ratio(UInt128{}, UInt128{}), 1.0);
    EXPECT_EQ(isim_ratio(UInt128::FromU64(1), UInt128::FromU64(4)), 0.25);
}

TEST(ISimKernelsTest, ScoreComparisonIsExact) {
    // 1/3 against 333333333333/1000000000000: doubles of the two differ only
    // past the 12th digit, the cross products decide it exactly.
    const ISimScore third{1u, 3u};
    const ISimScore close{333333333333u, 1000000000000u};
    EXPECT_TRUE(score_greater(third, close));
    EXPECT_FALSE(score_greater(close, third));
    // A zero denominator scores 1.
    EXPECT_FALSE(score_greater(ISimScore{1u, 1u}, ISimScore{0u, 0u}));
    EXPECT_FALSE(score_greater(ISimScore{0u, 0u}, ISimScore{1u, 1u}));
    EXPECT_EQ(score_ratio(ISimScore{0u, 0u}), 1.0);
    EXPECT_EQ(score_ratio(ISimScore{1u, 4u}), 0.25);
}

TEST(ISimKernelsTest, RowHelpersMaskTheTailWord) {
    // 70 bits: one full word plus a 6-bit tail.
    const auto batch = isim_test::make_batch({
        isim_test::make_fp(70, {0, 63, 64, 69}),
        isim_test::make_fp(70, {0, 64, 65}),
    });
    EXPECT_EQ(row_intersection(batch.RowWords(0), batch.RowWords(1), 70), 2u);
    BitCounts counts(70, 0u);
    add_row_to_counts(counts, batch.RowWords(0), 70);
    add_row_to_counts(counts, batch.RowWords(1), 70);
    EXPECT_EQ(counts[0], 2u);
    EXPECT_EQ(counts[64], 2u);
    EXPECT_EQ(counts[69], 1u);
    EXPECT_EQ(row_dot_counts(batch.RowWords(0), 70, counts), 2u + 1u + 2u + 1u);
    EXPECT_DOUBLE_EQ(row_tanimoto_distance(batch.RowWords(0), 4, batch.RowWords(1), 3, 70),
                     1.0 - 2.0 / 5.0);
}

TEST(ISimKernelsTest, TwoAllZeroRowsAreAtDistanceZero) {
    const auto batch = isim_test::make_batch({
        isim_test::make_fp(16, {}), isim_test::make_fp(16, {})});
    EXPECT_EQ(row_tanimoto_distance(batch.RowWords(0), 0, batch.RowWords(1), 0, 16), 0.0);
}

TEST(ISimKernelsTest, SumsMatchBruteForcePairSums) {
    const auto batch = isim_test::make_random_batch(25, 130, 12345u);
    BitCounts counts(130, 0u);
    for (size_t i = 0; i < batch.Size(); ++i) add_row_to_counts(counts, batch.RowWords(i), 130);
    const ISimSums sums = isim_sums(batch.Size(), count_moments(counts));
    uint64_t inter = 0, uni = 0;
    for (size_t i = 0; i < batch.Size(); ++i)
        for (size_t j = i + 1; j < batch.Size(); ++j) {
            inter += isim_test::intersection(batch, i, j);
            uni += isim_test::union_count(batch, i, j);
        }
    EXPECT_EQ(sums.intersections, UInt128::FromU64(inter));
    EXPECT_EQ(sums.unions, UInt128::FromU64(uni));
}

TEST(ISimTest, MatchesTheOracleAndTheLegacyKernel) {
    for (const uint64_t seed : {1u, 7u, 99u}) {
        const auto batch = isim_test::make_random_batch(30, 200, seed);
        uint64_t inter = 0, uni = 0;
        for (size_t i = 0; i < batch.Size(); ++i)
            for (size_t j = i + 1; j < batch.Size(); ++j) {
                inter += isim_test::intersection(batch, i, j);
                uni += isim_test::union_count(batch, i, j);
            }
        ASSERT_GT(uni, 0u);
        const double value = isim(batch);
        isim_test::expect_rel(value, isim_test::ratio(inter, uni), "isim");

        BitBirchLinearSum linear(200, 0u);
        for (size_t i = 0; i < batch.Size(); ++i) {
            UpdateLinearSumFromWords(linear, batch.RowWords(i), 200);
        }
        EXPECT_NEAR(value, JaccardTanimotoISim(linear, batch.Size()), 1e-12);
    }
}

TEST(ISimTest, NaNBelowTwoFingerprints) {
    EXPECT_TRUE(std::isnan(isim(isim_test::make_batch({isim_test::make_fp(8, {1})}))));
}

TEST(ISimTest, AllZeroBatchIsSimilarityOne) {
    const auto batch = isim_test::make_batch({
        isim_test::make_fp(8, {}), isim_test::make_fp(8, {}), isim_test::make_fp(8, {})});
    EXPECT_EQ(isim(batch), 1.0);
}

TEST(ISimTest, RefusesAnyMetricButTanimoto) {
    const auto batch = isim_test::make_batch({
        isim_test::make_fp(8, {1}), isim_test::make_fp(8, {2})});
    ISimOptions options;
    options.metric = "dice";
    EXPECT_THROW(isim(batch, options), std::invalid_argument);
}

TEST(ISimKernelsTest, RowDotColumnsSumsEveryClusterColumn) {
    // Two clusters' counts over 70 bits, stored column-major.
    std::vector<uint32_t> columns(70 * 2, 0u);
    columns[0 * 2 + 0] = 3u;
    columns[0 * 2 + 1] = 1u;
    columns[69 * 2 + 1] = 4u;
    const auto batch = isim_test::make_batch({isim_test::make_fp(70, {0, 69})});
    std::vector<uint64_t> dots;
    row_dot_columns(batch.RowWords(0), 70, columns, 2, dots);
    EXPECT_EQ(dots, (std::vector<uint64_t>{3u, 5u}));
}

TEST(ISimKernelsTest, MultiplyU128CarriesAcrossEveryLimb) {
    constexpr uint64_t MAX = std::numeric_limits<uint64_t>::max();
    // (2^128 - 1)^2 = 2^256 - 2^129 + 1.
    const UInt256 square = multiply_u128(UInt128{MAX, MAX}, UInt128{MAX, MAX});
    EXPECT_EQ(square.limb[0], 1u);
    EXPECT_EQ(square.limb[1], 0u);
    EXPECT_EQ(square.limb[2], MAX - 1u);
    EXPECT_EQ(square.limb[3], MAX);
    // 2^64 * 2^64 = 2^128 lands exactly in limb 2.
    const UInt256 power = multiply_u128(UInt128{1u, 0u}, UInt128{1u, 0u});
    EXPECT_EQ(power.limb[2], 1u);
    EXPECT_EQ(power.limb[0] | power.limb[1] | power.limb[3], 0u);
    EXPECT_TRUE(power < square);
    EXPECT_FALSE(square < power);
}

TEST(ISimKernelsTest, RatioComparisonIsExactWhereDoublesDisagree) {
    // Both ratios are exactly 1/3, but converting each side to double first
    // gives adjacent doubles, which would let the later ordinal win a tie.
    const UInt128 a_num = UInt128::FromU64(8100000180000000ull);
    const UInt128 a_den = UInt128::FromU64(3ull * 8100000180000000ull);
    const UInt128 b_num = UInt128::FromU64(8100000450000006ull);
    const UInt128 b_den = UInt128::FromU64(3ull * 8100000450000006ull);
    ASSERT_NE(a_num.ToDouble() / a_den.ToDouble(), b_num.ToDouble() / b_den.ToDouble());
    EXPECT_FALSE(isim_ratio_greater(a_num, a_den, b_num, b_den));
    EXPECT_FALSE(isim_ratio_greater(b_num, b_den, a_num, a_den));
    // Strict order, and operands past 2^64 on both sides.
    const UInt128 big{1u, 0u};
    EXPECT_TRUE(isim_ratio_greater(big, big + big, big, big + big + big));
    EXPECT_FALSE(isim_ratio_greater(big, big + big + big, big, big + big));
    // A zero union reads as similarity 1.
    EXPECT_TRUE(isim_ratio_greater(UInt128{}, UInt128{}, big, big + big));
    EXPECT_FALSE(isim_ratio_greater(big, big, UInt128{}, UInt128{}));
    EXPECT_FALSE(isim_ratio_greater(UInt128{}, UInt128{}, big, big));
}
