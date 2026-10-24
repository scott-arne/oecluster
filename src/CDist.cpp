/**
 * @file CDist.cpp
 * @brief Cross-distance computation engine implementation.
 */

#include "oecluster/CDist.h"

#include <algorithm>
#include <atomic>
#include <mutex>
#include <vector>

#include "oecluster/Error.h"
#include "oecluster/PairwiseComparison.h"
#include "oecluster/PoolMode.h"
#include "oecluster/ThreadPool.h"

namespace OECluster {

void cdist(PairwiseComparison& comparison, size_t n_a, double* output,
           const CDistOptions& options) {
    const size_t n_total = comparison.Size();
    const size_t n_b = n_total - n_a;
    const size_t total_pairs = n_a * n_b;

    // A positive cutoff zeroes every value above it, which for a similarity
    // discards exactly the most similar pairs. Only a reported similarity is
    // refused: an Unknown orientation is accepted so comparisons without facts
    // keep working. Placed ahead of the empty-domain return and the bulk
    // TryCDist path, which applies the cutoff itself, so the contract does not
    // depend on how many pairs there are or which path scores them.
    if (options.cutoff > 0.0 &&
        comparison.Facts().is_distance == Capability::No) {
        throw ComparisonError(
            "cdist cannot apply a cutoff to a comparison that reports "
            "similarities: the cutoff zeroes values above the threshold, which "
            "would discard the highest similarities. Set cutoff to 0.0, or "
            "build the comparison to report distances");
    }

    if (total_pairs == 0) return;

    // The bulk hook sees the same capped count as the generic loop below, so
    // an override that runs OpenEye code is serialised too.
    CDistOptions effective = options;
    effective.num_threads = safe_worker_count(comparison, options.num_threads);

    if (comparison.TryCDist(n_a, output, effective)) {
        return;
    }

    // A progress callback is called from whichever worker finished a chunk.
    // When the caller can service callbacks only on its own thread
    // (callbacks_on_calling_thread(), which a CallerPoolModeScope sets), one
    // worker, which ThreadPool runs on the calling thread, keeps every call
    // there.
    const size_t requested = options.progress && callbacks_on_calling_thread()
                                 ? 1
                                 : effective.num_threads;
    ThreadPool pool(requested);
    const size_t chunk_size = options.chunk_size > 0 ? options.chunk_size : 256;
    const size_t total_chunks =
        total_pairs / chunk_size + (total_pairs % chunk_size != 0 ? 1 : 0);

    // One clone per participant. A one-thread pool runs every chunk on the
    // calling thread. Otherwise detail::launch_workers starts exactly
    // min(NumThreads(), total_chunks) workers, each a fresh thread that takes
    // one ordinal below, so every index stays in range.
    const size_t participants =
        pool.NumThreads() == 1 ? 1 : std::min(pool.NumThreads(), total_chunks);
    std::vector<std::unique_ptr<PairwiseComparison>> clones;
    clones.reserve(participants);
    for (size_t t = 0; t < participants; ++t) {
        clones.push_back(comparison.Clone());
    }

    // Atomic counter for worker ordinal assignment
    std::atomic<size_t> thread_ordinal{0};

    // Progress tracking
    std::atomic<size_t> completed_pairs{0};
    std::mutex progress_mutex;

    pool.ParallelFor(0, total_pairs, chunk_size,
        [&](size_t chunk_begin, size_t chunk_end) {
            PairwiseComparison* local_comparison = clones.front().get();
            if (participants > 1) {
                // Each worker takes a unique ordinal on its first chunk.
                thread_local size_t my_ordinal =
                    thread_ordinal.fetch_add(1, std::memory_order_relaxed);
                local_comparison = clones[my_ordinal].get();
            }

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
