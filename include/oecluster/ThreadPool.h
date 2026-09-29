/**
 * @file ThreadPool.h
 * @brief Thread pool with dynamic chunk scheduling.
 */

#ifndef OECLUSTER_THREADPOOL_H
#define OECLUSTER_THREADPOOL_H

#include <atomic>
#include <cstddef>
#include <functional>
#include <memory>

namespace OECluster {

class ThreadPool {
public:
    /**
     * @brief Construct a thread pool with the specified worker count.
     *
     * :param num_threads: Number of worker threads; 0 auto-detects hardware concurrency.
     */
    explicit ThreadPool(size_t num_threads = 0);
    ~ThreadPool();

    ThreadPool(const ThreadPool&) = delete;
    ThreadPool& operator=(const ThreadPool&) = delete;

    /**
     * @brief Execute body(begin, end) for chunks spanning [begin, end).
     *
     * Chunks of size chunk_size are claimed by threads via atomic counter,
     * and at most one worker per chunk is started, so a call with fewer
     * chunks than NumThreads() starts only as many threads as it has chunks.
     * Exceptions from any worker are captured and re-thrown after all
     * workers complete.
     *
     * If a worker thread cannot be spawned, the workers already started are
     * joined and the std::system_error is propagated. A body exception or a
     * spawn failure leaves the pool cancelled, exactly as Cancel() does, so
     * later ParallelFor calls on the same pool do no work.
     *
     * :param begin: Start of range.
     * :param end: End of range (exclusive).
     * :param chunk_size: Items per work unit.
     * :param body: Callable(size_t chunk_begin, size_t chunk_end).
     * :raises std::system_error: When a worker thread cannot be spawned.
     */
    void ParallelFor(size_t begin, size_t end, size_t chunk_size,
                     std::function<void(size_t, size_t)> body);

    /**
     * @brief Get the configured worker count.
     *
     * :returns: Number of worker threads in the pool.
     */
    size_t NumThreads() const;

    /**
     * @brief Signal queued and running tasks to stop.
     */
    void Cancel();

    /**
     * @brief Check whether cancellation was requested.
     *
     * :returns: True if Cancel() has been called, false otherwise.
     */
    bool IsCancelled() const;

private:
    struct Impl;
    std::unique_ptr<Impl> pimpl_;
};

}  // namespace OECluster

#endif  // OECLUSTER_THREADPOOL_H
