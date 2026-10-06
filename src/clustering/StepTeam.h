/**
 * @file StepTeam.h
 * @brief A persistent team of threads that runs one barrier-separated step at a time.
 */

#ifndef OECLUSTER_SRC_CLUSTERING_STEPTEAM_H
#define OECLUSTER_SRC_CLUSTERING_STEPTEAM_H

#include <algorithm>
#include <atomic>
#include <condition_variable>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <functional>
#include <limits>
#include <mutex>
#include <thread>
#include <utility>
#include <vector>

namespace OECluster::detail {

/**
 * @brief The number of threads, the caller included, that a run of n items uses.
 *
 * :param num_threads: Requested count; 0 asks for the hardware count.
 * :param n: Item count.
 * :param hardware: What std::thread::hardware_concurrency() reported, which the
 *     standard allows to be 0.
 * :returns: 0 for no items, otherwise min(n, max(1, requested or hardware)).
 */
inline size_t resolve_participants(size_t num_threads, size_t n, size_t hardware) {
    if (n == 0) {
        return 0;
    }
    const size_t requested = num_threads > 0 ? num_threads : hardware;
    return std::min(n, std::max<size_t>(1, requested));
}

inline size_t resolve_participants(size_t num_threads, size_t n) {
    return resolve_participants(num_threads, n, std::thread::hardware_concurrency());
}

/**
 * @brief Work-unit size that spreads range over about four units per participant.
 *
 * A unit as large as the range would leave one thread doing all of it, so the
 * caller's ceiling only ever shrinks the unit.
 *
 * :param range: Items to divide.
 * :param participants: Threads sharing them, at least 1.
 * :param ceiling: Largest unit allowed.
 * :returns: max(1, min(ceiling, ceil(range / (4 x participants)))).
 */
inline size_t work_unit(size_t range, size_t participants,
                        size_t ceiling = std::numeric_limits<size_t>::max()) {
    const size_t workers = std::max<size_t>(1, participants);
    const size_t units = workers > std::numeric_limits<size_t>::max() / 4
                             ? std::numeric_limits<size_t>::max()
                             : 4 * workers;
    const size_t spread = range == 0 ? 0 : 1 + (range - 1) / units;
    return std::max<size_t>(1, std::min(ceiling, spread));
}

/**
 * @brief Threads started once, then released step by step.
 *
 * ThreadPool::ParallelFor starts fresh threads on every call, about 139 us per
 * call for 14 threads, which a loop of one call per item cannot afford. This
 * team starts participants - 1 threads once; Run() releases them for one step
 * and the calling thread takes part as participant 0. Idle threads spin, then
 * yield, then sleep, so a long run of serial steps does not burn a core each.
 */
class StepTeam {
public:
    /** @brief Starts a thread from a body; injectable so tests can fail it. */
    using Spawn = std::function<std::thread(std::function<void()>)>;

    explicit StepTeam(size_t participants)
        : StepTeam(participants, [](std::function<void()> body) {
              return std::thread(std::move(body));
          }) {}

    /**
     * :raises std::system_error: If a thread fails to start; the threads already
     *     started, including any waiting for a step, are stopped and joined first.
     */
    StepTeam(size_t participants, const Spawn& spawn)
        : participants_(std::max<size_t>(1, participants)) {
        threads_.reserve(participants_ - 1);
        try {
            for (size_t index = 1; index < participants_; ++index) {
                threads_.push_back(spawn([this, index] { WorkerLoop(index); }));
            }
        } catch (...) {
            Stop();
            throw;
        }
    }

    ~StepTeam() { Stop(); }

    StepTeam(const StepTeam&) = delete;
    StepTeam& operator=(const StepTeam&) = delete;

    size_t Participants() const { return participants_; }

    /**
     * @brief Run body(participant) on every participant and wait for all of them.
     *
     * :raises: The first exception any participant threw, once every participant
     *     has finished the step.
     */
    void Run(const std::function<void(size_t)>& body) {
        body_ = &body;
        arrived_.store(0, std::memory_order_relaxed);
        Publish();
        Execute(0);
        // Workers finish quickly, so the caller spins and yields rather than sleeps.
        size_t spins = 0;
        while (arrived_.load(std::memory_order_acquire) + 1 < participants_) {
            if (++spins > SPINS) {
                std::this_thread::yield();
            }
        }
        body_ = nullptr;
        if (error_) {
            std::exception_ptr error = std::move(error_);
            error_ = nullptr;
            failed_.store(false, std::memory_order_relaxed);
            std::rethrow_exception(error);
        }
    }

    /** @brief True once any participant has thrown in the current step. */
    bool Failed() const { return failed_.load(std::memory_order_relaxed); }

private:
    static constexpr size_t SPINS = 2000;
    static constexpr size_t YIELDS = 2000;

    void Execute(size_t index) {
        try {
            (*body_)(index);
        } catch (...) {
            std::lock_guard<std::mutex> lock(error_mutex_);
            if (!error_) {
                error_ = std::current_exception();
            }
            failed_.store(true, std::memory_order_relaxed);
        }
    }

    void WorkerLoop(size_t index) {
        uint64_t seen = 0;
        while (true) {
            seen = WaitForGeneration(seen);
            if (stop_.load(std::memory_order_acquire)) {
                return;
            }
            Execute(index);
            arrived_.fetch_add(1, std::memory_order_acq_rel);
        }
    }

    uint64_t WaitForGeneration(uint64_t seen) {
        for (size_t i = 0; i < SPINS; ++i) {
            const uint64_t now = generation_.load(std::memory_order_acquire);
            if (now != seen) {
                return now;
            }
        }
        for (size_t i = 0; i < YIELDS; ++i) {
            const uint64_t now = generation_.load(std::memory_order_acquire);
            if (now != seen) {
                return now;
            }
            std::this_thread::yield();
        }
        std::unique_lock<std::mutex> lock(mutex_);
        wake_.wait(lock, [&] {
            return generation_.load(std::memory_order_acquire) != seen;
        });
        return generation_.load(std::memory_order_acquire);
    }

    // The increment happens under the mutex so a worker that has checked the
    // predicate but not yet slept cannot miss the notification.
    void Publish() {
        {
            std::lock_guard<std::mutex> lock(mutex_);
            generation_.fetch_add(1, std::memory_order_acq_rel);
        }
        wake_.notify_all();
    }

    void Stop() {
        stop_.store(true, std::memory_order_release);
        Publish();
        for (std::thread& thread : threads_) {
            if (thread.joinable()) {
                thread.join();
            }
        }
        threads_.clear();
    }

    const size_t participants_;
    std::vector<std::thread> threads_;
    // Each on its own cache line: idle threads poll generation_, finishing
    // threads write arrived_, and every chunk claim reads failed_.
    alignas(64) std::atomic<uint64_t> generation_{0};
    alignas(64) std::atomic<size_t> arrived_{0};
    alignas(64) std::atomic<bool> failed_{false};
    alignas(64) std::atomic<bool> stop_{false};
    std::mutex mutex_;
    std::condition_variable wake_;
    std::mutex error_mutex_;
    std::exception_ptr error_;
    const std::function<void(size_t)>* body_ = nullptr;
};

}  // namespace OECluster::detail

#endif  // OECLUSTER_SRC_CLUSTERING_STEPTEAM_H
