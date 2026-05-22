/**
 * @file BitBirchKernels.cpp
 * @brief Internal binary fingerprint kernels for BitBirch.
 */

#include "BitBirchKernels.h"

#include <limits>
#include <stdexcept>

namespace OECluster::detail {
namespace {

size_t dense_word_count(const size_t size_bits) {
    return (size_bits + 63u) / 64u;
}

uint32_t popcount64(const uint64_t word) {
    return static_cast<uint32_t>(__builtin_popcountll(word));
}

void IncrementLinearSumFromSetBits(
    BitBirchLinearSum& linear_sum,
    const uint64_t* words,
    const size_t size_bits) {
    // Fingerprints are sparse in the intended workloads, so touch only set bits
    // while still masking padding bits beyond the declared fingerprint width.
    const size_t full_words = size_bits / 64u;
    for (size_t word_index = 0; word_index < full_words; ++word_index) {
        uint64_t word = words[word_index];
        while (word != 0u) {
            const size_t bit = word_index * 64u +
                               static_cast<size_t>(__builtin_ctzll(word));
            ++linear_sum[bit];
            word &= word - 1u;
        }
    }

    const size_t tail_bits = size_bits % 64u;
    if (tail_bits == 0u) {
        return;
    }
    uint64_t word = words[full_words] & ((uint64_t{1} << tail_bits) - 1u);
    while (word != 0u) {
        const size_t bit = full_words * 64u +
                           static_cast<size_t>(__builtin_ctzll(word));
        ++linear_sum[bit];
        word &= word - 1u;
    }
}

void AddCentroidWordsToLinearSum(
    BitBirchLinearSum& linear_sum,
    const std::vector<uint64_t>& centroid_words) {
    IncrementLinearSumFromSetBits(
        linear_sum,
        centroid_words.data(),
        linear_sum.size());
}

}  // namespace

void UpdateLinearSumFromWords(
    BitBirchLinearSum& linear_sum,
    const uint64_t* words,
    const size_t size_bits) {
    if (linear_sum.size() != size_bits) {
        throw std::invalid_argument("BitBirch linear sum size does not match fingerprint size");
    }

    IncrementLinearSumFromSetBits(linear_sum, words, size_bits);
}

std::vector<uint64_t> BinaryCentroid(
    const BitBirchLinearSum& linear_sum,
    const size_t n_samples) {
    std::vector<uint64_t> centroid(dense_word_count(linear_sum.size()), 0);

    // Match the reference BitBirch rule exactly: linear_sum >= n_samples * 0.5.
    // Empty summaries therefore produce all-one centroids over the valid width.
    const double threshold = static_cast<double>(n_samples) * 0.5;
    for (size_t bit = 0; bit < linear_sum.size(); ++bit) {
        if (static_cast<double>(linear_sum[bit]) >= threshold) {
            centroid[bit / 64u] |= uint64_t{1} << (bit % 64u);
        }
    }
    return centroid;
}

double TanimotoWords(
    const uint64_t* lhs_words,
    const uint32_t lhs_popcount,
    const uint64_t* rhs_words,
    const uint32_t rhs_popcount,
    const size_t words_per_fingerprint) {
    uint32_t intersection = 0;
    for (size_t word = 0; word < words_per_fingerprint; ++word) {
        intersection += popcount64(lhs_words[word] & rhs_words[word]);
    }
    const uint32_t denominator = lhs_popcount + rhs_popcount - intersection;
    return static_cast<double>(intersection) / static_cast<double>(denominator);
}

double TanimotoPackedVectors(
    const std::vector<uint64_t>& lhs_words,
    const uint32_t lhs_popcount,
    const std::vector<uint64_t>& rhs_words,
    const uint32_t rhs_popcount) {
    if (lhs_words.size() != rhs_words.size()) {
        throw std::invalid_argument("BitBirch packed vectors must have the same width");
    }
    return TanimotoWords(
        lhs_words.data(),
        lhs_popcount,
        rhs_words.data(),
        rhs_popcount,
        lhs_words.size());
}

uint32_t PopCountWords(const std::vector<uint64_t>& words) {
    uint32_t count = 0;
    for (const uint64_t word : words) {
        count += popcount64(word);
    }
    return count;
}

double JaccardTanimotoISim(
    const BitBirchLinearSum& linear_sum,
    const size_t n_objects) {
    double sum_kq = 0.0;
    double sum_kqsq = 0.0;
    for (const uint32_t count : linear_sum) {
        const double value = static_cast<double>(count);
        sum_kq += value;
        sum_kqsq += value * value;
    }
    const double a = (sum_kqsq - sum_kq) / 2.0;
    return a / (a + static_cast<double>(n_objects) * sum_kq - sum_kqsq);
}

bool AcceptBitBirchMerge(
    const BitBirchMergeCriterion criterion,
    const double threshold,
    const double tolerance,
    const BitBirchLinearSum& new_linear_sum,
    const std::vector<uint64_t>& new_centroid,
    const size_t new_n,
    const BitBirchLinearSum& old_linear_sum,
    const BitBirchLinearSum& nominee_linear_sum,
    const size_t old_n,
    const size_t nominee_n) {
    switch (criterion) {
        case BitBirchMergeCriterion::Radius: {
            BitBirchLinearSum radius_linear_sum = new_linear_sum;
            AddCentroidWordsToLinearSum(radius_linear_sum, new_centroid);
            return (JaccardTanimotoISim(radius_linear_sum, new_n + 1) *
                        static_cast<double>(new_n + 1) -
                    JaccardTanimotoISim(new_linear_sum, new_n) *
                        static_cast<double>(new_n - 1)) >= threshold * 2.0;
        }
        case BitBirchMergeCriterion::Diameter:
            return JaccardTanimotoISim(new_linear_sum, new_n) >= threshold;
        case BitBirchMergeCriterion::Tolerance:
            if (JaccardTanimotoISim(new_linear_sum, new_n) < threshold) {
                return false;
            }
            if (old_n == 1 && nominee_n == 1) {
                return true;
            }
            if (nominee_n == 1) {
                return ((JaccardTanimotoISim(new_linear_sum, old_n + 1) *
                             static_cast<double>(old_n + 1) -
                         JaccardTanimotoISim(old_linear_sum, old_n) *
                             static_cast<double>(old_n - 1)) /
                        2.0) >=
                       JaccardTanimotoISim(old_linear_sum, old_n) - tolerance;
            }
            return true;
        case BitBirchMergeCriterion::ToleranceTough:
            if (JaccardTanimotoISim(new_linear_sum, new_n) < threshold) {
                return false;
            }
            if (old_n == 1 && nominee_n == 1) {
                return true;
            }
            if (nominee_n == 1) {
                return ((JaccardTanimotoISim(new_linear_sum, old_n + 1) *
                             static_cast<double>(old_n + 1) -
                         JaccardTanimotoISim(old_linear_sum, old_n) *
                             static_cast<double>(old_n - 1)) /
                        2.0) >=
                       JaccardTanimotoISim(old_linear_sum, old_n) - tolerance;
            }
            return ((JaccardTanimotoISim(new_linear_sum, old_n + nominee_n) *
                         static_cast<double>(old_n + nominee_n) *
                         static_cast<double>(old_n + nominee_n - 1) -
                     JaccardTanimotoISim(old_linear_sum, old_n) *
                         static_cast<double>(old_n) *
                         static_cast<double>(old_n - 1) -
                     JaccardTanimotoISim(nominee_linear_sum, nominee_n) *
                         static_cast<double>(nominee_n) *
                         static_cast<double>(nominee_n - 1)) /
                    (2.0 * static_cast<double>(old_n) *
                     static_cast<double>(nominee_n))) >=
                   JaccardTanimotoISim(old_linear_sum, old_n) - tolerance;
    }
    return false;
}

}  // namespace OECluster::detail
