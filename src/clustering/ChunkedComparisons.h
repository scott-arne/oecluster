/**
 * @file ChunkedComparisons.h
 * @brief Chunked Compare() work with one comparison clone per running chunk.
 */

#ifndef OECLUSTER_SRC_CLUSTERING_CHUNKEDCOMPARISONS_H
#define OECLUSTER_SRC_CLUSTERING_CHUNKEDCOMPARISONS_H

#include <algorithm>
#include <cstddef>
#include <memory>
#include <mutex>
#include <utility>
#include <vector>

#include "DiversityValidation.h"
#include "oecluster/PairwiseComparison.h"
#include "oecluster/ThreadPool.h"

namespace OECluster::detail {

// Runs chunked Compare() work with one clone per concurrently running chunk.
// A thread_local ordinal (the pdist pattern) is safe only for a single
// ParallelFor call; these entry points make one per row or candidate, so
// clones circulate through a free list instead and their number is bounded by
// the concurrency actually reached.
class ChunkedComparisons {
public:
    ChunkedComparisons(const PairwiseComparison& prototype, size_t n,
                       size_t num_threads, size_t chunk_size)
        : prototype_(prototype),
          pool_(capped_threads(num_threads, n)),
          chunk_size_(chunk_size) {
        // At most one lease is live per pool thread, so returning a clone
        // never allocates in normal operation.
        free_.reserve(pool_.NumThreads());
    }

    // body(PairwiseComparison& clone, size_t begin, size_t end) over [0, range).
    template <typename Body>
    void Run(size_t range, Body&& body) {
        if (range == 0) {
            return;
        }
        // Capped at the range so a near-SIZE_MAX chunk_size cannot overflow
        // ParallelFor's ceiling arithmetic.
        const size_t chunk = std::min(chunk_size_, range);
        if (chunk == range) {
            // One chunk: ParallelFor would start every worker to run it on
            // one of them, once per row or candidate.
            Lease lease(*this);
            body(*lease.clone, 0, range);
            return;
        }
        pool_.ParallelFor(0, range, chunk, [&](size_t begin, size_t end) {
            Lease lease(*this);
            body(*lease.clone, begin, end);
        });
    }

private:
    struct Lease {
        explicit Lease(ChunkedComparisons& owner)
            : owner(owner), clone(owner.Acquire()) {}
        ~Lease() { owner.Release(std::move(clone)); }
        Lease(const Lease&) = delete;
        Lease& operator=(const Lease&) = delete;

        ChunkedComparisons& owner;
        std::unique_ptr<PairwiseComparison> clone;
    };

    // Clone() runs under the lock: it is const, but nothing promises that a
    // const call is safe to make concurrently on the prototype.
    std::unique_ptr<PairwiseComparison> Acquire() {
        std::lock_guard<std::mutex> lock(mutex_);
        if (free_.empty()) {
            return prototype_.Clone();
        }
        std::unique_ptr<PairwiseComparison> clone = std::move(free_.back());
        free_.pop_back();
        return clone;
    }

    // Called from ~Lease, which is noexcept, so a throw here would terminate
    // the process. If the clone cannot be returned to the free list, the
    // by-value parameter still owns it (push_back has the strong guarantee)
    // and destroys it; a later Acquire clones afresh.
    void Release(std::unique_ptr<PairwiseComparison> clone) noexcept {
        try {
            std::lock_guard<std::mutex> lock(mutex_);
            free_.push_back(std::move(clone));
        } catch (...) {
        }
    }

    const PairwiseComparison& prototype_;
    ThreadPool pool_;
    size_t chunk_size_;
    std::mutex mutex_;
    std::vector<std::unique_ptr<PairwiseComparison>> free_;
};

}  // namespace OECluster::detail

#endif  // OECLUSTER_SRC_CLUSTERING_CHUNKEDCOMPARISONS_H
