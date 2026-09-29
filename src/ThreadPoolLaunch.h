/**
 * @file ThreadPoolLaunch.h
 * @brief Thread spawn-and-join launcher with failure recovery.
 */

#ifndef OECLUSTER_SRC_THREADPOOLLAUNCH_H
#define OECLUSTER_SRC_THREADPOOLLAUNCH_H

#include <cstddef>
#include <functional>
#include <thread>
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

}  // namespace detail
}  // namespace OECluster

#endif  // OECLUSTER_SRC_THREADPOOLLAUNCH_H
