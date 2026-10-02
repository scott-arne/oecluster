/**
 * @file ISimKernels.cpp
 * @brief Integer iSIM kernels.
 */

#include "ISimKernels.h"

#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

#if defined(_MSC_VER)
#include <intrin.h>
#endif

namespace OECluster::detail {
namespace {

// Copies of the helpers in BitBirchKernels.cpp, kept separate because this
// feature leaves BitBirch untouched. MSVC has no __builtin_* bit intrinsics.
uint32_t popcount64(const uint64_t word) {
#if defined(_MSC_VER)
    return static_cast<uint32_t>(__popcnt64(word));
#elif defined(__clang__) || defined(__GNUC__)
    return static_cast<uint32_t>(__builtin_popcountll(static_cast<unsigned long long>(word)));
#else
    uint32_t count = 0;
    uint64_t remaining = word;
    while (remaining != 0u) {
        remaining &= remaining - 1u;
        ++count;
    }
    return count;
#endif
}

// Undefined for word == 0; every caller tests the word first.
uint32_t count_trailing_zeros64(const uint64_t word) {
#if defined(_MSC_VER)
    unsigned long index = 0;
    _BitScanForward64(&index, word);
    return static_cast<uint32_t>(index);
#elif defined(__clang__) || defined(__GNUC__)
    return static_cast<uint32_t>(__builtin_ctzll(static_cast<unsigned long long>(word)));
#else
    uint32_t count = 0;
    uint64_t remaining = word;
    while ((remaining & 1u) == 0u) {
        remaining >>= 1;
        ++count;
    }
    return count;
#endif
}

uint64_t tail_mask(const size_t size_bits) {
    return (uint64_t{1} << (size_bits % 64u)) - 1u;
}

// Visits only set bits, which is what makes the core O(sum of popcounts)
// rather than O(N * bits) on sparse fingerprints.
template <class Visit>
void for_each_set_bit(const uint64_t* words, const size_t size_bits, Visit&& visit) {
    const size_t full_words = size_bits / 64u;
    for (size_t w = 0; w < full_words; ++w) {
        uint64_t word = words[w];
        while (word != 0u) {
            visit(w * 64u + static_cast<size_t>(count_trailing_zeros64(word)));
            word &= word - 1u;
        }
    }
    if (size_bits % 64u != 0u) {
        uint64_t word = words[full_words] & tail_mask(size_bits);
        while (word != 0u) {
            visit(full_words * 64u + static_cast<size_t>(count_trailing_zeros64(word)));
            word &= word - 1u;
        }
    }
}

}  // namespace

UInt128 UInt128::FromU64(const uint64_t value) {
    return UInt128{0u, value};
}

UInt128& UInt128::operator+=(const UInt128& other) {
    const uint64_t previous = low;
    low += other.low;
    high += other.high + (low < previous ? 1u : 0u);
    return *this;
}

UInt128& UInt128::operator-=(const UInt128& other) {
    const uint64_t borrow = low < other.low ? 1u : 0u;
    low -= other.low;
    high -= other.high + borrow;
    return *this;
}

UInt128 UInt128::Half() const {
    return UInt128{high >> 1, (low >> 1) | (high << 63)};
}

double UInt128::ToDouble() const {
    // high < 2^53 converts exactly and ldexp is exact, so the sum is the one
    // rounding.
    return std::ldexp(static_cast<double>(high), 64) + static_cast<double>(low);
}

UInt128 operator+(UInt128 lhs, const UInt128& rhs) { return lhs += rhs; }
UInt128 operator-(UInt128 lhs, const UInt128& rhs) { return lhs -= rhs; }
bool operator==(const UInt128& lhs, const UInt128& rhs) {
    return lhs.high == rhs.high && lhs.low == rhs.low;
}
bool operator!=(const UInt128& lhs, const UInt128& rhs) { return !(lhs == rhs); }
bool operator<(const UInt128& lhs, const UInt128& rhs) {
    return lhs.high != rhs.high ? lhs.high < rhs.high : lhs.low < rhs.low;
}

UInt128 multiply_u64(const uint64_t lhs, const uint64_t rhs) {
    constexpr uint64_t MASK32 = 0xffffffffu;
    const uint64_t lhs_lo = lhs & MASK32;
    const uint64_t lhs_hi = lhs >> 32;
    const uint64_t rhs_lo = rhs & MASK32;
    const uint64_t rhs_hi = rhs >> 32;
    const uint64_t lo_lo = lhs_lo * rhs_lo;
    const uint64_t hi_lo = lhs_hi * rhs_lo;
    const uint64_t lo_hi = lhs_lo * rhs_hi;
    const uint64_t hi_hi = lhs_hi * rhs_hi;
    // At most (2^32 - 1) * 2 + (2^32 - 1)^2 = 2^64 - 1, so this cannot wrap.
    const uint64_t cross = (lo_lo >> 32) + (hi_lo & MASK32) + lo_hi;
    return UInt128{hi_hi + (hi_lo >> 32) + (cross >> 32), (cross << 32) | (lo_lo & MASK32)};
}

CountMoments count_moments(const BitCounts& counts) {
    CountMoments moments;
    for (const uint32_t count : counts) {
        moments.sum += count;
        moments.square_sum += multiply_u64(count, count);
    }
    return moments;
}

ISimSums isim_sums(const uint64_t set_size, const CountMoments& moments) {
    const UInt128 sum = UInt128::FromU64(moments.sum);
    ISimSums sums;
    sums.intersections = (moments.square_sum - sum).Half();
    sums.unions = multiply_u64(set_size, moments.sum) - (moments.square_sum + sum).Half();
    return sums;
}

double isim_ratio(const UInt128& intersections, const UInt128& unions) {
    if (unions == UInt128{}) {
        return 1.0;
    }
    return intersections.ToDouble() / unions.ToDouble();
}

void check_isim_batch_size(const size_t num_fingerprints) {
    if (static_cast<uint64_t>(num_fingerprints) > std::numeric_limits<uint32_t>::max()) {
        throw std::invalid_argument(
            "iSIM supports at most 4294967295 fingerprints, got " +
            std::to_string(num_fingerprints));
    }
}

void check_isim_width(const uint64_t size_bits) {
    if (size_bits >= (uint64_t{1} << 31)) {
        throw std::invalid_argument(
            "iSIM supports fingerprints narrower than 2147483648 bits, got " +
            std::to_string(size_bits));
    }
}

void add_row_to_counts(BitCounts& counts, const uint64_t* words, const size_t size_bits) {
    for_each_set_bit(words, size_bits, [&counts](const size_t bit) { ++counts[bit]; });
}

uint64_t row_dot_counts(const uint64_t* words, const size_t size_bits, const BitCounts& counts) {
    uint64_t dot = 0;
    for_each_set_bit(words, size_bits, [&](const size_t bit) { dot += counts[bit]; });
    return dot;
}

uint32_t row_intersection(const uint64_t* lhs, const uint64_t* rhs, const size_t size_bits) {
    const size_t full_words = size_bits / 64u;
    uint32_t count = 0;
    for (size_t w = 0; w < full_words; ++w) {
        count += popcount64(lhs[w] & rhs[w]);
    }
    if (size_bits % 64u != 0u) {
        count += popcount64(lhs[full_words] & rhs[full_words] & tail_mask(size_bits));
    }
    return count;
}

double row_tanimoto_distance(const uint64_t* lhs, const uint32_t lhs_popcount,
                             const uint64_t* rhs, const uint32_t rhs_popcount,
                             const size_t size_bits) {
    const uint32_t shared = row_intersection(lhs, rhs, size_bits);
    const uint64_t unions = uint64_t{lhs_popcount} + rhs_popcount - shared;
    if (unions == 0u) {
        return 0.0;
    }
    return 1.0 - static_cast<double>(shared) / static_cast<double>(unions);
}

ISimScore isim_member_score(const uint64_t dot, const uint32_t popcount,
                            const uint64_t set_size, const uint64_t set_sum) {
    const uint64_t intersections = dot - popcount;
    const uint64_t unions =
        (set_size - 1u) * popcount + (set_sum - popcount) - intersections;
    return ISimScore{intersections, unions};
}

bool score_greater(const ISimScore& lhs, const ISimScore& rhs) {
    const ISimScore a = lhs.denominator == 0u ? ISimScore{1u, 1u} : lhs;
    const ISimScore b = rhs.denominator == 0u ? ISimScore{1u, 1u} : rhs;
    return multiply_u64(b.numerator, a.denominator) < multiply_u64(a.numerator, b.denominator);
}

double score_ratio(const ISimScore& score) {
    if (score.denominator == 0u) {
        return 1.0;
    }
    return static_cast<double>(score.numerator) / static_cast<double>(score.denominator);
}

}  // namespace OECluster::detail
