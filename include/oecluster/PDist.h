/**
 * @file PDist.h
 * @brief Pairwise distance computation engine.
 */

#ifndef OECLUSTER_PDIST_H
#define OECLUSTER_PDIST_H

#include <cstddef>
#include <functional>

namespace OECluster {

class PairwiseComparison;
class StorageBackend;

/**
 * @brief Options for controlling pairwise distance computation.
 */
struct PDistOptions {
    size_t num_threads = 0;   ///< Number of threads (0 = auto-detect)
    size_t chunk_size = 256;  ///< Number of pairs per work unit
    /// Not read by pdist or by the bundled TryPDist overrides: filtering is
    /// the storage backend's job, so pass a SparseStorage to drop distances
    /// above a cutoff. Kept for callers that record the cutoff alongside the
    /// options.
    double cutoff = 0.0;

    /// Progress callback: (completed_pairs, total_pairs)
    std::function<void(size_t completed, size_t total)> progress;
};

/**
 * @brief Compute all pairwise distances.
 *
 * Distributes work across threads using a ThreadPool. Each thread
 * gets its own Clone() of the comparison for thread safety. Results
 * are written into the storage backend.
 *
 * A SparseStorage keeps only values at or below its own cutoff, which for a
 * similarity would discard the highest values, so a comparison whose
 * ``Facts().is_distance`` is ``Capability::No`` is refused with any
 * SparseStorage, whatever its cutoff. ``Capability::Unknown`` is accepted.
 *
 * :param comparison: Pairwise comparison (will be Clone()'d per thread).
 * :param storage: Storage backend to write results into.
 * :param options: Threading, chunking, and progress options.
 * :raises ComparisonError: If ``storage.NumSamples()`` differs from
 *     ``comparison.Size()``, or if the comparison reports similarities and
 *     ``storage`` is a SparseStorage.
 */
void pdist(PairwiseComparison& comparison, StorageBackend& storage,
           const PDistOptions& options = {});

}  // namespace OECluster

#endif  // OECLUSTER_PDIST_H
