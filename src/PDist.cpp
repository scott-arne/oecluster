/**
 * @file PDist.cpp
 * @brief Pairwise distance computation engine implementation.
 */

#include "oecluster/PDist.h"

#include <atomic>
#include <cmath>
#include <mutex>
#include <vector>

#include "oecluster/DistanceMetric.h"
#include "oecluster/StorageBackend.h"
#include "oecluster/ThreadPool.h"

namespace OECluster {

namespace {

/**
 * @brief Convert condensed index k back to pair (i, j).
 *
 * Given condensed index k for N items, recovers the row i and column j
 * such that i < j.
 *
 * :param k: Condensed index in [0, N*(N-1)/2).
 * :param n: Total number of items.
 * :param out_i: Output row index.
 * :param out_j: Output column index.
 */
void condensed_to_pair(size_t k, size_t n, size_t& out_i, size_t& out_j) {
    auto nd = static_cast<double>(n);
    auto kd = static_cast<double>(k);
    auto i_d = nd - 2.0 - std::floor(std::sqrt(-8.0 * kd + 4.0 * nd * (nd - 1.0) - 7.0) / 2.0 - 0.5);
    auto i = static_cast<size_t>(i_d);
    auto j = k + i + 1 - n * (n - 1) / 2 + (n - i) * ((n - i) - 1) / 2;
    out_i = i;
    out_j = j;
}

}  // namespace

void pdist(DistanceMetric& metric, StorageBackend& storage,
           const PDistOptions& options) {
    const size_t n = metric.Size();
    const size_t total_pairs = n * (n - 1) / 2;

    if (total_pairs == 0) {
        storage.Finalize();
        return;
    }

    ThreadPool pool(options.num_threads);
    const size_t num_threads = pool.NumThreads();
    const size_t chunk_size = options.chunk_size > 0 ? options.chunk_size : 256;

    // Create one clone per thread
    std::vector<std::unique_ptr<DistanceMetric>> clones;
    clones.reserve(num_threads);
    for (size_t t = 0; t < num_threads; ++t) {
        clones.push_back(metric.Clone());
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
            auto& local_metric = clones[my_ordinal];

            size_t i, j;
            for (size_t k = chunk_begin; k < chunk_end; ++k) {
                condensed_to_pair(k, n, i, j);
                double distance = local_metric->Distance(i, j);
                storage.Set(i, j, distance);
            }

            // Update progress
            if (options.progress) {
                size_t done = completed_pairs.fetch_add(
                    chunk_end - chunk_begin, std::memory_order_relaxed)
                    + (chunk_end - chunk_begin);
                std::lock_guard<std::mutex> lock(progress_mutex);
                options.progress(done, total_pairs);
            }
        });

    storage.Finalize();
}

}  // namespace OECluster
