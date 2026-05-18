/**
 * @file CDist.h
 * @brief Cross-distance computation engine.
 */

#ifndef OECLUSTER_CDIST_H
#define OECLUSTER_CDIST_H

#include <cstddef>
#include <functional>

namespace OECluster {

class PairwiseComparison;

/**
 * @brief Options for controlling cross-distance computation.
 */
struct CDistOptions {
    size_t num_threads = 0;   ///< Number of threads (0 = auto-detect)
    size_t chunk_size = 256;  ///< Number of pairs per work unit
    double cutoff = 0.0;      ///< Distance cutoff (0.0 = no cutoff)

    /// Progress callback: (completed_pairs, total_pairs)
    std::function<void(size_t completed, size_t total)> progress;
};

/**
 * @brief Compute cross-distances between two item sets packed in one comparison.
 *
 * The comparison must have been constructed with items [set_a..., set_b...].
 * n_a items from set A, comparison.Size() - n_a items from set B.
 * Output is a flat array of size n_a * n_b in row-major order.
 *
 * :param comparison: Pairwise comparison constructed with concatenated items.
 * :param n_a: Number of items in set A (first n_a items in comparison).
 * :param output: Pre-allocated array of size n_a * (comparison.Size() - n_a).
 * :param options: Threading and progress options.
 */
void cdist(PairwiseComparison& comparison, size_t n_a, double* output,
           const CDistOptions& options = {});

}  // namespace OECluster

#endif  // OECLUSTER_CDIST_H
