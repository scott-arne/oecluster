/**
 * @file KMedoidsDistanceSource.h
 * @brief Where k-medoids reads a distance, so one kernel serves two inputs.
 */

#ifndef OECLUSTER_SRC_CLUSTERING_KMEDOIDSDISTANCESOURCE_H
#define OECLUSTER_SRC_CLUSTERING_KMEDOIDSDISTANCESOURCE_H

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <utility>

#include "ChunkedComparisons.h"
#include "DistanceAccess.h"
#include "ThresholdGraph.h"
#include "oecluster/PairwiseComparison.h"
#include "oecluster/ThreadPool.h"

namespace OECluster::detail {

/**
 * @brief Clamp a chunk size to the item count.
 *
 * The swap kernel sizes its per-chunk staging as ``(n + chunk - 1) / chunk``,
 * which wraps to zero for a chunk size near SIZE_MAX, leaving no slot for the
 * one chunk the pool runs. Clamping to the item count is observationally
 * free -- a chunk at least as wide as the range is one chunk either way --
 * and keeps that arithmetic below the overflow.
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

/**
 * @brief Reads a distance by scoring the pair, with no matrix behind it.
 *
 * Mirrors dense_distance() deliberately: zero on the diagonal *without*
 * scoring, and otherwise the pair read in ascending index order. pdist() fills
 * a matrix from condensed_to_pair(), which yields ``i < j``, so a differently
 * oriented read would see a different matrix for any comparison that is not
 * symmetric, and the streamed and matrix paths would stop agreeing.
 * PrimMST.cpp canonicalizes the same way for the same reason.
 *
 * Unlike the matrix path it checks finiteness as values arrive. The matrix
 * path may assume it, because the Python gate scans the stored matrix before
 * a kernel reads it; a streamed run has no stored matrix to scan.
 */
class StreamingDistances {
public:
    /**
     * @brief Bind one leased clone and the name diagnostics will use.
     *
     * :param comparison: Clone leased for the chunk being run; must outlive
     *     this accessor, which the lease guarantees.
     * :param caller: Entry point name for error messages.
     */
    StreamingDistances(PairwiseComparison& comparison, const std::string& caller)
        : comparison_(comparison), caller_(caller) {}

    double operator()(size_t i, size_t j) const {
        if (i == j) {
            return 0.0;
        }
        const size_t low = std::min(i, j);
        const size_t high = std::max(i, j);
        const double distance = comparison_.Compare(low, high);
        if (!std::isfinite(distance)) {
            throw non_finite_distance_error(caller_, low, high);
        }
        return distance;
    }

private:
    PairwiseComparison& comparison_;
    const std::string& caller_;
};

/**
 * @brief Runs a k-medoids pass by scoring pairs on demand.
 *
 * Interchangeable with MatrixDistanceRunner: same Size(), ChunkSize(), Run()
 * and WithAccessor(), so the kernels never learn which input they are reading.
 *
 * Clones come from ChunkedComparisons' free list rather than one per pass.
 * That is the whole reason the seam is a runner: a streamed run makes
 * ``k + iterations`` passes, ParallelFor on a pool of more than one thread
 * spawns fresh threads for each, so a per-pass ordinal scheme would mint a new
 * set of clones every pass. A lease bounds them by the pool size instead, and
 * a lease has to be scoped to the work that uses it.
 *
 * The partition matches the matrix runner's exactly for every chunk size
 * validation admits. Both hand ParallelFor the chunk size clamped to the
 * range, and ParallelFor derives every chunk boundary from that alone; the
 * worker count only decides who claims a chunk. ChunkedComparisons runs a
 * single-chunk pass inline rather than through ParallelFor, which is the same
 * one chunk either way.
 */
class StreamingDistanceRunner {
public:
    /**
     * @brief Bind a comparison and the schedule every pass over it will use.
     *
     * :param comparison: Prototype to clone from; must outlive this runner.
     * :param caller: Entry point name for error messages.
     * :param num_threads: Worker threads; 0 auto-detects hardware concurrency.
     *     Already capped by the caller, so no pass underneath can raise it.
     * :param chunk_size: Caller-requested chunk size; validation guarantees at
     *     least one, and a zero is raised to one rather than honoured. Zero is
     *     not a schedule: ParallelFor declines a zero chunk and runs nothing,
     *     so honouring it would make ``WithAccessor`` skip its body and return
     *     as though the work were done. Raising it keeps every pass running;
     *     the two runners then disagree only on an input validation rejects,
     *     and only by this one doing the work instead of silently skipping it.
     */
    StreamingDistanceRunner(const PairwiseComparison& comparison, std::string caller,
                            size_t num_threads, size_t chunk_size)
        : n_(comparison.Size()),
          caller_(std::move(caller)),
          chunk_size_(std::max<size_t>(chunk_size, 1)),
          work_(comparison, n_, num_threads, chunk_size_) {}

    size_t Size() const { return n_; }

    /**
     * @brief The chunk width ``Run`` will partition a range into.
     *
     * The same arithmetic the matrix runner uses, and the same ChunkedComparisons
     * applies internally, so a kernel staging one result per chunk sizes and
     * indexes that staging identically on either input.
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
     * :param range: Length of the range to cover.
     * :param body: Callable(accessor, chunk_begin, chunk_end).
     */
    template <typename Body>
    void Run(size_t range, Body&& body) {
        RunGuarded([&]() {
            work_.Run(range, [&](PairwiseComparison& local, size_t begin, size_t end) {
                body(StreamingDistances(local, caller_), begin, end);
            });
        });
    }

    /**
     * @brief Hand one accessor to serial work on the calling thread.
     *
     * A range of one leases a single clone and runs it inline, which is how
     * ChunkedComparisons already handles a pass it cannot split. The chunk
     * size cannot be zero, so this always reaches the body -- see the
     * constructor for why that is worth guaranteeing here.
     *
     * :param body: Callable(accessor).
     */
    template <typename Body>
    void WithAccessor(Body&& body) {
        RunGuarded([&]() {
            work_.Run(1, [&](PairwiseComparison& local, size_t, size_t) {
                body(StreamingDistances(local, caller_));
            });
        });
    }

private:
    /**
     * @brief Run one pass, and refuse every pass after one that threw.
     *
     * The matrix runner builds a pool per pass because ThreadPool's cancelled
     * flag is sticky. This runner cannot: the pool is what bounds the clone
     * count, so it is owned for the runner's lifetime, and a pass that threw
     * leaves it cancelled. A later pass would then do no work at all and the
     * kernel would read its default-initialized output as a real answer -- a
     * zero-cost single cluster, reported as converged. Refusing turns that
     * into a diagnosis. It also removes a schedule dependence: a single-chunk
     * pass runs inline and never reaches the pool, so without this the
     * behaviour after a throw would depend on the chunk width.
     *
     * The hazard is general to ChunkedComparisons rather than particular to
     * k-medoids -- every entry point that calls Run() more than once on one
     * instance has it -- so guarding here is the narrow fix, chosen over
     * changing when eight other callers throw.
     *
     * :param pass: Callable() running the pass.
     * :raises std::logic_error: When an earlier pass on this runner threw.
     */
    template <typename Pass>
    void RunGuarded(Pass&& pass) {
        if (failed_) {
            throw std::logic_error(
                caller_ + ": a pass failed, so this distance source is spent");
        }
        try {
            pass();
        } catch (...) {
            failed_ = true;
            throw;
        }
    }

    size_t n_;
    std::string caller_;
    size_t chunk_size_;
    bool failed_ = false;
    ChunkedComparisons work_;
};

}  // namespace OECluster::detail

#endif  // OECLUSTER_SRC_CLUSTERING_KMEDOIDSDISTANCESOURCE_H
