/**
 * @file ExactMedian.h
 * @brief Exact median of a value stream in fixed-size state.
 *
 * Used by cluster_report for the two medians taken over pairs. Up to the
 * budget the values are stored and sorted, as before; above it the median is
 * selected by four 16-bit radix passes over a re-walkable stream, so memory
 * stays fixed however many pairs there are.
 */

#ifndef OECLUSTER_SRC_CLUSTERING_EXACTMEDIAN_H
#define OECLUSTER_SRC_CLUSTERING_EXACTMEDIAN_H

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace OECluster::detail {

constexpr size_t MEDIAN_DIRECT_BUDGET = size_t{1} << 20;

namespace median_bits {
constexpr uint64_t SIGN = uint64_t{1} << 63;
constexpr size_t DIGITS = size_t{1} << 16;
constexpr int NUM_DIGITS = 4;
}  // namespace median_bits

/// The signed total-order transform: keys compare as unsigned integers
/// exactly as the doubles compare, negatives included.
inline uint64_t median_key(double value) {
    uint64_t bits = 0;
    std::memcpy(&bits, &value, sizeof(bits));
    return (bits & median_bits::SIGN) != 0 ? ~bits : bits ^ median_bits::SIGN;
}

inline double median_value(uint64_t key) {
    const uint64_t bits =
        (key & median_bits::SIGN) != 0 ? key ^ median_bits::SIGN : ~key;
    double value = 0.0;
    std::memcpy(&value, &bits, sizeof(value));
    return value;
}

/**
 * @brief Multi-pass exact median over a stream of `count` values.
 *
 * Each pass visits every value once, in any order, then calls EndPass().
 * Its result equals detail::median_distance (sort; middle value, or the mean
 * of the two middle values) except that -0.0 is normalized to +0.0 first, so
 * a zero median is always +0.0. The values must be identical on every pass.
 */
class ExactMedian {
public:
    explicit ExactMedian(size_t count, size_t budget = MEDIAN_DIRECT_BUDGET)
        : count_(count),
          direct_(count <= budget),
          done_(count == 0),
          rank_{count == 0 ? 0 : (count - 1) / 2, count / 2} {
        if (done_) {
            return;
        }
        if (direct_) {
            values_.reserve(count);
        } else {
            hist_.assign(2 * median_bits::DIGITS, 0);
        }
    }

    void Visit(double value) {
        if (done_) {
            throw std::logic_error("ExactMedian: Visit after the median is known");
        }
        if (visited_ == count_) {
            throw std::logic_error("ExactMedian: more values than the declared count");
        }
        ++visited_;
        if (value == 0.0) {
            value = 0.0;
        }
        if (direct_) {
            values_.push_back(value);
            return;
        }
        const uint64_t key = median_key(value);
        const int shift = 48 - 16 * digit_;
        const size_t bucket = static_cast<size_t>((key >> shift) & 0xFFFF);
        // A key counts toward rank t's histogram only while it still shares
        // the digits already fixed for that rank.
        const int tracks = same_ ? 1 : 2;
        for (int t = 0; t < tracks; ++t) {
            if (digit_ == 0 || (key >> (shift + 16)) == prefix_[t]) {
                ++hist_[static_cast<size_t>(t) * median_bits::DIGITS + bucket];
            }
        }
    }

    void EndPass() {
        if (done_) {
            throw std::logic_error("ExactMedian: EndPass after the median is known");
        }
        if (visited_ != count_) {
            throw std::logic_error(
                "ExactMedian pass visited " + std::to_string(visited_) +
                " values, expected " + std::to_string(count_));
        }
        visited_ = 0;
        if (direct_) {
            std::sort(values_.begin(), values_.end());
            low_ = values_[rank_[0]];
            high_ = values_[rank_[1]];
            values_.clear();
            values_.shrink_to_fit();
            done_ = true;
            return;
        }
        uint64_t next[2] = {0, 0};
        for (int t = 0; t < 2; ++t) {
            const size_t* h = &hist_[static_cast<size_t>(same_ ? 0 : t) * median_bits::DIGITS];
            size_t r = rank_[t];
            size_t d = 0;
            while (d < median_bits::DIGITS && r >= h[d]) {
                r -= h[d];
                ++d;
            }
            if (d == median_bits::DIGITS) {
                throw std::logic_error("ExactMedian: the values changed between passes");
            }
            rank_[t] = r;
            next[t] = (prefix_[t] << 16) | static_cast<uint64_t>(d);
        }
        prefix_[0] = next[0];
        prefix_[1] = next[1];
        same_ = prefix_[0] == prefix_[1];
        std::fill(hist_.begin(), hist_.end(), size_t{0});
        ++digit_;
        if (digit_ == median_bits::NUM_DIGITS) {
            low_ = median_value(prefix_[0]);
            high_ = median_value(prefix_[1]);
            hist_.clear();
            hist_.shrink_to_fit();
            done_ = true;
        }
    }

    bool Done() const { return done_; }

    double Result() const {
        if (count_ == 0) {
            return std::numeric_limits<double>::quiet_NaN();
        }
        if (!done_) {
            throw std::logic_error("ExactMedian: Result before the last pass");
        }
        // Normalized after averaging as well as before keying: the mean of
        // -denorm_min and +0.0 underflows to -0.0.
        const double median = count_ % 2 == 1 ? low_ : (low_ + high_) / 2.0;
        return median == 0.0 ? 0.0 : median;
    }

    /// The direct route over a caller-owned buffer, with no budget: sorts it.
    /// Only the result's zero sign is normalized, never the buffer's entries:
    /// under pair-rank the buffer goes on to pair_rank_indices. Which zero
    /// sits at a rank cannot change a nonzero median, so this agrees with
    /// normalizing every value first.
    static double OfInPlace(std::vector<double>& values) {
        if (values.empty()) {
            return std::numeric_limits<double>::quiet_NaN();
        }
        std::sort(values.begin(), values.end());
        const size_t middle = values.size() / 2;
        double median = values.size() % 2 == 1
            ? values[middle]
            : (values[middle - 1] + values[middle]) / 2.0;
        if (median == 0.0) {
            median = 0.0;
        }
        return median;
    }

private:
    size_t count_;
    bool direct_;
    bool done_;
    size_t visited_ = 0;
    // The two target ranks, low and high; equal for an odd count. Within a
    // radix pass each is the rank remaining inside its current prefix.
    size_t rank_[2];
    uint64_t prefix_[2] = {0, 0};
    bool same_ = true;
    int digit_ = 0;
    double low_ = 0.0;
    double high_ = 0.0;
    std::vector<double> values_;
    std::vector<size_t> hist_;
};

}  // namespace OECluster::detail

#endif  // OECLUSTER_SRC_CLUSTERING_EXACTMEDIAN_H
