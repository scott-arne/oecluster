/**
 * @file PDist.cpp
 * @brief Pairwise distance computation engine implementation.
 */

#include "oecluster/PDist.h"

#include <algorithm>
#include <atomic>
#include <mutex>
#include <vector>

#include "CompareUnordered.h"
#include "oecluster/CondensedIndex.h"
#include "oecluster/Error.h"
#include "oecluster/PairwiseComparison.h"
#include "oecluster/PoolMode.h"
#include "oecluster/StorageBackend.h"
#include "oecluster/ThreadPool.h"

namespace OECluster {

void pdist(PairwiseComparison& comparison, StorageBackend& storage,
           const PDistOptions& options) {
    const size_t n = comparison.Size();

    // Each backend's Set() re-maps the pair through its own sample count, so a
    // storage larger than the comparison silently scatters the results across
    // the wrong slots -- Set's own range check cannot see that, because every
    // index the loop produces is in range for the oversized storage. Checked
    // before the empty-domain shortcut so the contract does not depend on how
    // many pairs there happen to be.
    if (storage.NumSamples() != n) {
        throw ComparisonError("pdist storage size does not match the comparison size");
    }

    // A SparseStorage keeps only values at or below its cutoff, which for a
    // similarity discards exactly the most similar pairs. Every SparseStorage
    // filters, whatever the sign of its cutoff, so the refusal keys on the
    // backend type rather than on the value; PDistOptions::cutoff is not read
    // by the driver at all. Only a reported similarity is refused: an Unknown
    // orientation is accepted so comparisons without facts keep working. Placed
    // ahead of the empty-domain shortcut and the bulk TryPDist path, which
    // writes into the same storage, for the same reason as the size check.
    if (dynamic_cast<SparseStorage*>(&storage) != nullptr &&
        comparison.Facts().is_distance == Capability::No) {
        throw ComparisonError(
            "pdist cannot write a comparison that reports similarities into "
            "SparseStorage: its cutoff keeps only values at or below the "
            "threshold, which would discard the highest similarities. Use "
            "DenseStorage or MMapStorage, or build the comparison to report "
            "distances");
    }

    const size_t total_pairs = n * (n - 1) / 2;

    if (total_pairs == 0) {
        storage.Finalize();
        return;
    }

    // The bulk hook sees the same capped count as the generic loop below, so
    // an override that runs OpenEye code is serialised too.
    PDistOptions effective = options;
    effective.num_threads = safe_worker_count(comparison, options.num_threads);

    if (comparison.TryPDist(storage, effective)) {
        storage.Finalize();
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

            size_t i, j;
            for (size_t k = chunk_begin; k < chunk_end; ++k) {
                condensed_to_pair(k, n, i, j);
                // The unordered seam, not Compare: only one triangle is stored
                // and every reader presents it as symmetric, so a directional
                // comparison would otherwise make the matrix depend on the
                // caller's input order. The base implementation forwards to
                // Compare, so this changes nothing for a comparison that does
                // not override it.
                double distance =
                    detail::compare_unordered(*local_comparison, i, j);
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
