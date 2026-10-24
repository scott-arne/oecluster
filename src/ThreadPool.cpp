/**
 * @file ThreadPool.cpp
 * @brief Thread pool implementation with dynamic chunk scheduling.
 */

#include "oecluster/ThreadPool.h"

#include <algorithm>
#include <atomic>
#include <exception>
#include <mutex>
#include <thread>
#include <vector>

#include "ThreadPoolLaunch.h"

namespace OECluster {

struct ThreadPool::Impl {
    size_t num_threads_;
    std::atomic<bool> cancelled_{false};

    explicit Impl(size_t num_threads)
        : num_threads_(num_threads > 0
                           ? num_threads
                           : std::max<size_t>(1, std::thread::hardware_concurrency())) {}
};

ThreadPool::ThreadPool(size_t num_threads)
    : pimpl_(std::make_unique<Impl>(num_threads)) {}

ThreadPool::~ThreadPool() = default;

size_t ThreadPool::NumThreads() const {
    return pimpl_->num_threads_;
}

void ThreadPool::Cancel() {
    pimpl_->cancelled_.store(true, std::memory_order_relaxed);
}

bool ThreadPool::IsCancelled() const {
    return pimpl_->cancelled_.load(std::memory_order_relaxed);
}

void ThreadPool::ParallelFor(size_t begin, size_t end, size_t chunk_size,
                             std::function<void(size_t, size_t)> body) {
    if (begin >= end || chunk_size == 0) {
        return;
    }

    const size_t range = end - begin;
    const size_t total_chunks = (range + chunk_size - 1) / chunk_size;

    std::atomic<size_t> next_chunk{0};
    std::exception_ptr captured_exception;
    std::once_flag exception_flag;

    auto worker = [&]() {
        while (true) {
            if (pimpl_->cancelled_.load(std::memory_order_relaxed)) {
                break;
            }

            size_t chunk_idx = next_chunk.fetch_add(1, std::memory_order_relaxed);
            if (chunk_idx >= total_chunks) {
                break;
            }

            size_t chunk_begin = begin + chunk_idx * chunk_size;
            size_t chunk_end = std::min(chunk_begin + chunk_size, end);

            try {
                body(chunk_begin, chunk_end);
            } catch (...) {
                std::call_once(exception_flag, [&]() {
                    captured_exception = std::current_exception();
                });
                pimpl_->cancelled_.store(true, std::memory_order_relaxed);
                break;
            }
        }
    };

    if (pimpl_->num_threads_ == 1) {
        // A one-thread pool runs on the calling thread rather than on a worker
        // of its own. A caller that holds a lock a worker would need while it
        // waits -- the GIL, for a binding that keeps it -- would deadlock if a
        // worker had to take that lock, as a progress callback that calls into
        // Python does; on the calling thread, taking it again is re-entrant.
        // The worker loop is unchanged, so chunk order, cancellation and the
        // captured exception behave as before.
        worker();
    } else {
        detail::launch_workers(
            pimpl_->num_threads_, total_chunks,
            [&worker]() { return std::thread(worker); },
            [this]() { pimpl_->cancelled_.store(true, std::memory_order_relaxed); });
    }

    if (captured_exception) {
        std::rethrow_exception(captured_exception);
    }
}

}  // namespace OECluster
