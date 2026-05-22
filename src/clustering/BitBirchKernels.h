/**
 * @file BitBirchKernels.h
 * @brief Internal binary fingerprint kernels for BitBirch.
 */

#ifndef OECLUSTER_CLUSTERING_BITBIRCH_KERNELS_H
#define OECLUSTER_CLUSTERING_BITBIRCH_KERNELS_H

#include <cstddef>
#include <cstdint>
#include <vector>

#include "oecluster/clustering/BitBirch.h"

namespace OECluster::detail {

using BitBirchLinearSum = std::vector<uint32_t>;

void UpdateLinearSumFromWords(
    BitBirchLinearSum& linear_sum,
    const uint64_t* words,
    size_t size_bits);

std::vector<uint64_t> BinaryCentroid(
    const BitBirchLinearSum& linear_sum,
    size_t n_samples);

double TanimotoWords(
    const uint64_t* lhs_words,
    uint32_t lhs_popcount,
    const uint64_t* rhs_words,
    uint32_t rhs_popcount,
    size_t words_per_fingerprint);

double TanimotoPackedVectors(
    const std::vector<uint64_t>& lhs_words,
    uint32_t lhs_popcount,
    const std::vector<uint64_t>& rhs_words,
    uint32_t rhs_popcount);

uint32_t PopCountWords(const std::vector<uint64_t>& words);

double JaccardTanimotoISim(const BitBirchLinearSum& linear_sum, size_t n_objects);

bool AcceptBitBirchMerge(
    BitBirchMergeCriterion criterion,
    double threshold,
    double tolerance,
    const BitBirchLinearSum& new_linear_sum,
    const std::vector<uint64_t>& new_centroid,
    size_t new_n,
    const BitBirchLinearSum& old_linear_sum,
    const BitBirchLinearSum& nominee_linear_sum,
    size_t old_n,
    size_t nominee_n);

}  // namespace OECluster::detail

#endif  // OECLUSTER_CLUSTERING_BITBIRCH_KERNELS_H
