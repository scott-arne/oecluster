/**
 * @file ThreadPoolLaunch.h
 * @brief Thread spawn-and-join launcher with failure recovery.
 */

#ifndef OECLUSTER_SRC_THREADPOOLLAUNCH_H
#define OECLUSTER_SRC_THREADPOOLLAUNCH_H

#include <algorithm>
#include <cstddef>
#include <functional>
#include <thread>
#include <utility>
#include <vector>

namespace OECluster {
namespace detail {

/**
 * @brief Spawn worker threads and join them all, with failure recovery.
 *
 * If spawning a thread throws, the failure callback runs, all threads already
 * started are joined, and the original exception is rethrown. Without this,
 * the threads vector would unwind with joinable threads in it, which calls
 * std::terminate.
 *
 * :param num_threads: Number of threads to spawn.
 * :param spawn: Callable that returns std::thread (typically constructs one).
 * :param on_failure: Callback invoked if spawn throws.
 */
template <typename Spawn, typename OnFailure>
void launch_threads(size_t num_threads, Spawn&& spawn, OnFailure&& on_failure) {
    std::vector<std::thread> threads;
    threads.reserve(num_threads);

    try {
        for (size_t i = 0; i < num_threads; ++i) {
            threads.emplace_back(spawn());
        }
    } catch (...) {
        on_failure();
        for (auto& t : threads) {
            if (t.joinable()) {
                t.join();
            }
        }
        throw;
    }

    for (auto& t : threads) {
        t.join();
    }
}

/**
 * @brief Launch at most one worker per chunk, then join them.
 *
 * Workers claim chunks from a shared counter, so any worker beyond the chunk
 * count would find nothing to do. Skipping those workers avoids the spawn
 * cost, which dominates when the per-chunk work is cheap. Spawn failures are
 * handled as in launch_threads.
 *
 * :param num_threads: Pool size, the upper bound on workers.
 * :param total_chunks: Number of chunks to be claimed; 0 starts no threads.
 * :param spawn: Callable that returns std::thread (typically constructs one).
 * :param on_failure: Callback invoked if spawn throws.
 */
template <typename Spawn, typename OnFailure>
void launch_workers(size_t num_threads, size_t total_chunks, Spawn&& spawn,
                    OnFailure&& on_failure) {
    launch_threads(std::min(num_threads, total_chunks), std::forward<Spawn>(spawn),
                   std::forward<OnFailure>(on_failure));
}

}  // namespace detail
}  // namespace OECluster

#endif  // OECLUSTER_SRC_THREADPOOLLAUNCH_H
