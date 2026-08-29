/**
 * @file PDist.cpp
 * @brief Pairwise distance computation engine implementation.
 */

#include "oecluster/PDist.h"

#include <atomic>
#include <mutex>
#include <vector>

#include "oecluster/CondensedIndex.h"
#include "oecluster/PairwiseComparison.h"
#include "oecluster/StorageBackend.h"
#include "oecluster/ThreadPool.h"

namespace OECluster {

void pdist(PairwiseComparison& comparison, StorageBackend& storage,
           const PDistOptions& options) {
    const size_t n = comparison.Size();
    const size_t total_pairs = n * (n - 1) / 2;

    if (total_pairs == 0) {
        storage.Finalize();
        return;
    }

    if (comparison.TryPDist(storage, options)) {
        storage.Finalize();
        return;
    }

    ThreadPool pool(options.num_threads);
    const size_t num_threads = pool.NumThreads();
    const size_t chunk_size = options.chunk_size > 0 ? options.chunk_size : 256;

    // Create one clone per thread
    std::vector<std::unique_ptr<PairwiseComparison>> clones;
    clones.reserve(num_threads);
    for (size_t t = 0; t < num_threads; ++t) {
        clones.push_back(comparison.Clone());
    }

    // Atomic counter for thread-local index assignment
    std::atomic<size_t> thread_ordinal{0};

    // Progress tracking
    std::atomic<size_t> completed_pairs{0};
    std::mutex progress_mutex;

    pool.ParallelFor(0, total_pairs, chunk_size,
        [&](size_t chunk_begin, size_t chunk_end) {
            // Assign each thread a unique ordinal on first entry
            thread_local size_t my_ordinal = thread_ordinal.fetch_add(1, std::memory_order_relaxed);
            auto& local_comparison = clones[my_ordinal];

            size_t i, j;
            for (size_t k = chunk_begin; k < chunk_end; ++k) {
                condensed_to_pair(k, n, i, j);
                double distance = local_comparison->Compare(i, j);
                storage.Set(i, j, distance);
            }

            if (options.progress) {
                // Keep the counter update and callback together so callers
                // observe monotonically increasing progress from worker chunks.
                std::lock_guard<std::mutex> lock(progress_mutex);
                size_t done = completed_pairs.fetch_add(
                    chunk_end - chunk_begin, std::memory_order_relaxed)
                    + (chunk_end - chunk_begin);
                options.progress(done, total_pairs);
            }
        });

    storage.Finalize();
}

}  // namespace OECluster
