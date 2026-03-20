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
    explicit ThreadPool(size_t num_threads = 0);
    ~ThreadPool();

    ThreadPool(const ThreadPool&) = delete;
    ThreadPool& operator=(const ThreadPool&) = delete;

    /**
     * @brief Execute body(begin, end) for chunks spanning [begin, end).
     *
     * Chunks of size chunk_size are claimed by threads via atomic counter.
     * Exceptions from any worker are captured and re-thrown after all
     * workers complete.
     *
     * :param begin: Start of range.
     * :param end: End of range (exclusive).
     * :param chunk_size: Items per work unit.
     * :param body: Callable(size_t chunk_begin, size_t chunk_end).
     */
    void ParallelFor(size_t begin, size_t end, size_t chunk_size,
                     std::function<void(size_t, size_t)> body);

    size_t NumThreads() const;
    void Cancel();
    bool IsCancelled() const;

private:
    struct Impl;
    std::unique_ptr<Impl> pimpl_;
};

}  // namespace OECluster

#endif  // OECLUSTER_THREADPOOL_H
