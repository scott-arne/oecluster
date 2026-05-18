/**
 * @file CDist.cpp
 * @brief Cross-distance computation engine implementation.
 */

#include "oecluster/CDist.h"

#include <atomic>
#include <mutex>
#include <vector>

#include "oecluster/PairwiseComparison.h"
#include "oecluster/ThreadPool.h"

namespace OECluster {

void cdist(PairwiseComparison& comparison, size_t n_a, double* output,
           const CDistOptions& options) {
    const size_t n_total = comparison.Size();
    const size_t n_b = n_total - n_a;
    const size_t total_pairs = n_a * n_b;

    if (total_pairs == 0) return;

    if (comparison.TryCDist(n_a, output, options)) {
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
            thread_local size_t my_ordinal =
                thread_ordinal.fetch_add(1, std::memory_order_relaxed);
            auto& local_comparison = clones[my_ordinal];

            for (size_t k = chunk_begin; k < chunk_end; ++k) {
                size_t i = k / n_b;
                size_t j = n_a + (k % n_b);
                double distance = local_comparison->Compare(i, j);
                if (options.cutoff > 0.0 && distance > options.cutoff) {
                    output[k] = 0.0;
                } else {
                    output[k] = distance;
                }
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
}

}  // namespace OECluster
