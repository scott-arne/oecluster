/**
 * @file ISimKernels.h
 * @brief Integer iSIM kernels shared by isim() and isim_report().
 *
 * Every count and sum is an integer, so a result never depends on the thread
 * count or summation order; only the final conversion to double and the
 * division round. Bounds hold for widths below 2^31 bits (check_isim_width) and batches
 * below 2^32 rows (check_isim_batch_size): per-bit counts then fit uint32_t,
 * set sums and medoid-score terms fit uint64_t, and squares, dot products of
 * count vectors and union sums need the 128-bit type.
 */

#ifndef OECLUSTER_CLUSTERING_ISIM_KERNELS_H
#define OECLUSTER_CLUSTERING_ISIM_KERNELS_H

#include <cstddef>
#include <cstdint>
#include <vector>

namespace OECluster::detail {

/**
 * @brief Portable unsigned 128-bit integer.
 *
 * __int128 is unavailable under MSVC, which builds the Windows wheel, and the
 * iSIM square sums pass 2^64 at about 67 million 4096-bit fingerprints.
 */
struct UInt128 {
    uint64_t high = 0;
    uint64_t low = 0;

    static UInt128 FromU64(uint64_t value);
    UInt128& operator+=(const UInt128& other);
    /// Requires *this >= other. Every iSIM difference is non-negative by
    /// construction, so a wrap would be a logic error rather than an input.
    UInt128& operator-=(const UInt128& other);
    /// floor(value / 2); every halved iSIM quantity is even.
    UInt128 Half() const;
    /// Correctly rounded while high < 2^53, i.e. for values below 2^117,
    /// which every iSIM sum is.
    double ToDouble() const;
};

UInt128 operator+(UInt128 lhs, const UInt128& rhs);
UInt128 operator-(UInt128 lhs, const UInt128& rhs);
bool operator==(const UInt128& lhs, const UInt128& rhs);
bool operator!=(const UInt128& lhs, const UInt128& rhs);
bool operator<(const UInt128& lhs, const UInt128& rhs);

UInt128 multiply_u64(uint64_t lhs, uint64_t rhs);

/// Unsigned 256-bit product holder, least significant limb first. Exists only
/// to compare cluster-to-cluster iSIM ratios exactly: their numerators and
/// unions reach about 2^94, so a cross product needs up to 2^188.
struct UInt256 {
    uint64_t limb[4] = {0u, 0u, 0u, 0u};
};

bool operator==(const UInt256& lhs, const UInt256& rhs);
bool operator<(const UInt256& lhs, const UInt256& rhs);

UInt256 multiply_u128(const UInt128& lhs, const UInt128& rhs);

/// Per-bit on-counts over a set of fingerprints, one entry per declared bit.
using BitCounts = std::vector<uint32_t>;

/// S = sum c and Q = c . c of a count vector.
struct CountMoments {
    uint64_t sum = 0;
    UInt128 square_sum;
};

CountMoments count_moments(const BitCounts& counts);

/// Summed pairwise intersections and unions of a set: the two sides of its
/// iSIM Tanimoto ratio.
struct ISimSums {
    UInt128 intersections;
    UInt128 unions;
};

/// A = (Q - S) / 2 and U = n S - (Q + S) / 2.
ISimSums isim_sums(uint64_t set_size, const CountMoments& moments);

/// intersections / unions, or 1.0 when unions == 0: a zero union sum means
/// every pair is two all-zero fingerprints, whose Tanimoto distance the
/// library defines as 0.
double isim_ratio(const UInt128& intersections, const UInt128& unions);

/// lhs_intersections / lhs_unions > rhs_intersections / rhs_unions, decided by
/// exact cross-multiplication; a zero union reads as 1, as in isim_ratio.
/// Converting each side to double first can round two equal ratios apart.
bool isim_ratio_greater(const UInt128& lhs_intersections, const UInt128& lhs_unions,
                        const UInt128& rhs_intersections, const UInt128& rhs_unions);

/// Refuses 2^32 or more fingerprints, where a uint32_t per-bit count could
/// wrap. A count rather than a batch so the boundary is testable without
/// materialising billions of rows.
void check_isim_batch_size(size_t num_fingerprints);

/// Refuses widths of 2^31 or more bits. A medoid-score union is below
/// 2 * set_size * bits, which fits uint64_t only under this bound and the
/// batch-size one. A width rather than a batch, for the same testability.
void check_isim_width(uint64_t size_bits);

/// Adds one row's set bits to counts; padding bits past size_bits are masked.
void add_row_to_counts(BitCounts& counts, const uint64_t* words, size_t size_bits);

/// x . c for one row against a count vector.
uint64_t row_dot_counts(const uint64_t* words, size_t size_bits, const BitCounts& counts);

/// x . c_l for every cluster l at once. columns holds the K count vectors
/// column-major (columns[bit * cluster_count + l] = c_l[bit]), so each set bit
/// reads one contiguous K-run. dots is resized to cluster_count.
void row_dot_columns(const uint64_t* words, size_t size_bits,
                     const std::vector<uint32_t>& columns, size_t cluster_count,
                     std::vector<uint64_t>& dots);

/// |x AND y| over the declared width.
uint32_t row_intersection(const uint64_t* lhs, const uint64_t* rhs, size_t size_bits);

/// Exact Tanimoto distance; 0.0 for two all-zero rows rather than the 0/0
/// that TanimotoWords in BitBirchKernels.cpp would divide.
double row_tanimoto_distance(const uint64_t* lhs, uint32_t lhs_popcount,
                             const uint64_t* rhs, uint32_t rhs_popcount,
                             size_t size_bits);

/// An iSIM score held as an exact fraction. A zero denominator scores 1.
struct ISimScore {
    uint64_t numerator = 0;
    uint64_t denominator = 0;
};

/**
 * @brief Score of one member against the rest of its set.
 *
 * dot is x . c over the set INCLUDING the member, so the member's own
 * contribution is removed here: I = dot - popcount is its summed intersection
 * with the other members and (set_size - 1) * popcount + (set_sum - popcount)
 * - I its summed union with them.
 */
ISimScore isim_member_score(uint64_t dot, uint32_t popcount, uint64_t set_size,
                            uint64_t set_sum);

/// lhs > rhs, decided by exact cross-multiplication.
bool score_greater(const ISimScore& lhs, const ISimScore& rhs);

/// The score as a double, 1.0 for a zero denominator.
double score_ratio(const ISimScore& score);

}  // namespace OECluster::detail

#endif  // OECLUSTER_CLUSTERING_ISIM_KERNELS_H
