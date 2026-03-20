/**
 * @file PDist.h
 * @brief Pairwise distance computation engine.
 */

#ifndef OECLUSTER_PDIST_H
#define OECLUSTER_PDIST_H

#include <cstddef>
#include <functional>

namespace OECluster {

class DistanceMetric;
class StorageBackend;

/**
 * @brief Options for controlling pairwise distance computation.
 */
struct PDistOptions {
    size_t num_threads = 0;   ///< Number of threads (0 = auto-detect)
    size_t chunk_size = 256;  ///< Number of pairs per work unit
    double cutoff = 0.0;      ///< Distance cutoff (0.0 = no cutoff)

    /// Progress callback: (completed_pairs, total_pairs)
    std::function<void(size_t completed, size_t total)> progress;
};

/**
 * @brief Compute all pairwise distances.
 *
 * Distributes work across threads using a ThreadPool. Each thread
 * gets its own Clone() of the metric for thread safety. Results
 * are written into the storage backend.
 *
 * :param metric: Distance metric (will be Clone()'d per thread).
 * :param storage: Storage backend to write results into.
 * :param options: Threading, chunking, cutoff, and progress options.
 */
void pdist(DistanceMetric& metric, StorageBackend& storage,
           const PDistOptions& options = {});

}  // namespace OECluster

#endif  // OECLUSTER_PDIST_H
