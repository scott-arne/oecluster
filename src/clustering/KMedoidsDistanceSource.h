/**
 * @file KMedoidsDistanceSource.h
 * @brief Where k-medoids reads a distance, so one kernel serves two inputs.
 */

#ifndef OECLUSTER_SRC_CLUSTERING_KMEDOIDSDISTANCESOURCE_H
#define OECLUSTER_SRC_CLUSTERING_KMEDOIDSDISTANCESOURCE_H

#include <algorithm>
#include <cstddef>

#include "DistanceAccess.h"
#include "oecluster/ThreadPool.h"

namespace OECluster::detail {

/**
 * @brief Clamp a chunk size to the item count.
 *
 * ThreadPool derives its chunk count as ``(range + chunk_size - 1) /
 * chunk_size``, which wraps to zero for a chunk size near SIZE_MAX and then
 * runs no chunk at all: every scan in the swap kernel would return its
 * default-initialized output and the call would report a zero-cost single
 * cluster as converged. Clamping to the item count is observationally free --
 * a chunk at least as wide as the range is one chunk either way -- and keeps
 * both the kernel's own ceiling arithmetic and the pool's below the overflow.
 * Validation guarantees ``chunk_size >= 1``, and every call site has
 * ``n >= n_clusters >= 1``, so the result is never zero.
 *
 * :param n: Number of items.
 * :param chunk_size: Caller-requested chunk size; at least one.
 * :returns: The chunk size to hand ThreadPool, never zero and never above n.
 */
inline size_t effective_chunk_size(size_t n, size_t chunk_size) {
    return std::min(chunk_size, n);
}

/**
 * @brief Reads a distance out of a materialized matrix.
 *
 * A thin wrapper over dense_distance() so the matrix path keeps compiling to
 * what it compiled to before the seam existed.
 */
class MatrixDistances {
public:
    MatrixDistances(const double* data, size_t n) : data_(data), n_(n) {}

    double operator()(size_t i, size_t j) const {
        return dense_distance(data_, n_, i, j);
    }

private:
    const double* data_;
    size_t n_;
};

/**
 * @brief Runs a k-medoids pass over a materialized matrix.
 *
 * A runner rather than a plain accessor because the streaming source added
 * later cannot hand one out: it leases a comparison clone per running chunk,
 * and a lease has to be scoped to the work that uses it. Owning the loop is
 * what lets both sources answer the same two questions -- how is this pass
 * partitioned, and what reads a distance inside it -- so the kernels never
 * learn which input they are reading.
 *
 * ``Run`` and ``WithAccessor`` are the only ways to obtain an accessor, and
 * both are non-const: running a pass is a mutating operation for the
 * streaming runner, which this one matches rather than tempting a caller to
 * hold a shared source const.
 */
class MatrixDistanceRunner {
public:
    /**
     * @brief Bind a matrix and the schedule every pass over it will use.
     *
     * :param data: Condensed distance array from ``StorageBackend::Data()``.
     * :param n: Number of items.
     * :param num_threads: Worker threads; 0 auto-detects hardware concurrency.
     *     Already capped by the caller, so no pass underneath can raise it.
     * :param chunk_size: Caller-requested chunk size; at least one.
     */
    MatrixDistanceRunner(const double* data, size_t n, size_t num_threads,
                         size_t chunk_size)
        : data_(data), n_(n), num_threads_(num_threads), chunk_size_(chunk_size) {}

    size_t Size() const { return n_; }

    /**
     * @brief The chunk width ``Run`` will partition a range into.
     *
     * Exposed because a kernel that stages one result per chunk has to size
     * and index that staging the same way the runner splits the work.
     *
     * :param range: Length of the range to be run.
     * :returns: The chunk width, never above the range.
     */
    size_t ChunkSize(size_t range) const {
        return effective_chunk_size(range, chunk_size_);
    }

    /**
     * @brief Run one pass over [0, range), chunk by chunk.
     *
     * The pool is built per pass because every kernel built one per call
     * before this seam existed, and ThreadPool's cancelled flag is sticky: a
     * pool reused across passes would silently do no work in every pass after
     * one body threw.
     *
     * :param range: Length of the range to cover.
     * :param body: Callable(accessor, chunk_begin, chunk_end).
     */
    template <typename Body>
    void Run(size_t range, Body&& body) {
        const MatrixDistances distance(data_, n_);
        ThreadPool pool(num_threads_);
        pool.ParallelFor(0, range, ChunkSize(range),
                         [&](size_t begin, size_t end) {
            body(distance, begin, end);
        });
    }

    /**
     * @brief Hand one accessor to serial work on the calling thread.
     *
     * :param body: Callable(accessor).
     */
    template <typename Body>
    void WithAccessor(Body&& body) {
        body(MatrixDistances(data_, n_));
    }

private:
    const double* data_;
    size_t n_;
    size_t num_threads_;
    size_t chunk_size_;
};

}  // namespace OECluster::detail

#endif  // OECLUSTER_SRC_CLUSTERING_KMEDOIDSDISTANCESOURCE_H
