/**
 * @file PoolMode.h
 * @brief The OpenEye memory-pool mode a call runs under, and the worker counts it allows.
 *
 * OpenEye's molecule memory pool is safe to drive from several threads only
 * when its mode has the Mutexed or ThreadLocal bit. Under any other mode, work
 * that runs OpenEye code must stay on one thread. This header decides how many
 * workers such work may use.
 *
 * The mode that matters is the one the application set. A process normally
 * has one OpenEye library, so that is OESystem::OEGetMemPoolMode(). An
 * embedder whose OpenEye copy is not this library's, as the Python extension
 * is when it links OpenEye statically, records the application's mode for the
 * duration of a call with CallerPoolModeScope.
 *
 * What this guarantees for direct C++ callers is serial work within one call.
 * It does not serialise two calls made at once from different threads:
 * OpenEye requires single-threaded use under these modes, so a C++ application
 * that sets one must serialise its own calls into oecluster and OpenEye. Under
 * an unsafe mode the Python bindings keep the GIL for the call, except while a
 * progress callback runs Python code, so two calls made through them never run
 * their native work at once. That does not cover native code another thread
 * is already running without the GIL, an OpenEye call included: a Python
 * application that sets such a mode must still keep its OpenEye work off
 * other threads meanwhile.
 */

#ifndef OECLUSTER_POOLMODE_H
#define OECLUSTER_POOLMODE_H

#include <cstddef>
#include <optional>

namespace OECluster {

class PairwiseComparison;

/**
 * @brief Whether a memory-pool mode permits driving OpenEye from threads.
 *
 * :param mode: An ``OESystem::OEMemPoolMode`` word.
 * :returns: ``true`` when the Mutexed or the ThreadLocal bit is set.
 */
bool pool_is_thread_safe(unsigned int mode);

/**
 * @brief The mode the current call runs under.
 *
 * :returns: The mode recorded by the innermost live CallerPoolModeScope on
 *     the calling thread that recorded one, or ``OESystem::OEGetMemPoolMode()``
 *     when none did.
 */
unsigned int effective_pool_mode();

/**
 * @brief The mode a call would run under, given a mode its caller reports.
 *
 * :param caller_mode: The mode the caller reports, if it reports one.
 * :returns: ``caller_mode`` when present, otherwise effective_pool_mode().
 */
unsigned int resolve_pool_mode(std::optional<unsigned int> caller_mode);

/**
 * @brief Whether the caller can service callbacks only on its own thread.
 *
 * :returns: The flag of the innermost live CallerPoolModeScope on the calling
 *     thread, or ``false`` when there is none.
 */
bool callbacks_on_calling_thread();

/**
 * @brief The worker count a comparison may use under the current mode.
 *
 * :param comparison: The comparison whose ``Compare``, protected
 *     ``CompareUnordered`` seam, ``TryPDist`` or ``TryCDist`` the workers
 *     will run.
 * :param requested: The caller's request; 0 means hardware concurrency.
 * :returns: 1 when the comparison runs OpenEye code and effective_pool_mode()
 *     is not thread-safe, otherwise ``requested`` unchanged.
 */
size_t safe_worker_count(const PairwiseComparison& comparison, size_t requested);

/**
 * @brief Record, for the calling thread and one call, the mode the
 * application set.
 *
 * Scopes nest. An empty ``mode`` records nothing and leaves an outer scope's
 * mode in force; the callback flag always describes the innermost scope. The
 * destructor restores what the constructor replaced. Copy and move are
 * deleted, because a copy's destructor would restore the outer state while
 * the original scope is still live.
 *
 * Create it on the stack of the thread whose call it describes, so that it
 * is destroyed on that thread in reverse order of creation. Its ``operator
 * new`` and ``operator new[]`` are deleted, which rejects ``new`` and
 * ``std::make_unique``. A route that bypasses a class ``operator new`` --
 * ``::new``, ``std::make_shared``, or a heap-allocated object that holds a
 * scope, such as a ``std::optional`` -- still compiles, and is just as wrong.
 */
class CallerPoolModeScope {
public:
    /**
     * :param mode: The application's memory-pool mode, if known.
     * :param callbacks_on_calling_thread: Whether the caller can service a
     *     callback only on its own thread, as the Python bindings can while
     *     they hold the GIL.
     */
    explicit CallerPoolModeScope(std::optional<unsigned int> mode,
                                 bool callbacks_on_calling_thread = false);
    ~CallerPoolModeScope();

    CallerPoolModeScope(const CallerPoolModeScope&) = delete;
    CallerPoolModeScope& operator=(const CallerPoolModeScope&) = delete;
    CallerPoolModeScope(CallerPoolModeScope&&) = delete;
    CallerPoolModeScope& operator=(CallerPoolModeScope&&) = delete;

    // A scope must be destroyed on the thread that created it, in reverse
    // order of creation, which heap allocation cannot promise.
    static void* operator new(std::size_t) = delete;
    static void* operator new[](std::size_t) = delete;

private:
    std::optional<unsigned int> previous_mode_;
    bool previous_callbacks_;
};

}  // namespace OECluster

#endif  // OECLUSTER_POOLMODE_H
