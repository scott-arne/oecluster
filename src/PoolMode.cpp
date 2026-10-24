/**
 * @file PoolMode.cpp
 * @brief The memory-pool mode a call runs under, and the worker counts it allows.
 */

#include "oecluster/PoolMode.h"

#include <oesystem.h>

#include "oecluster/PairwiseComparison.h"

namespace OECluster {

namespace {

// Per thread, because a call's scope must not leak to another thread's call.
thread_local std::optional<unsigned int> caller_mode;
thread_local bool caller_callbacks_on_calling_thread = false;

}  // namespace

bool pool_is_thread_safe(unsigned int mode) {
    return (mode & (OESystem::OEMemPoolMode::ThreadLocal |
                    OESystem::OEMemPoolMode::Mutexed)) != 0u;
}

unsigned int effective_pool_mode() {
    return caller_mode ? *caller_mode : OESystem::OEGetMemPoolMode();
}

unsigned int resolve_pool_mode(std::optional<unsigned int> mode) {
    return mode ? *mode : effective_pool_mode();
}

bool callbacks_on_calling_thread() {
    return caller_callbacks_on_calling_thread;
}

size_t safe_worker_count(const PairwiseComparison& comparison, size_t requested) {
    if (comparison.RunsOpenEyeCode() && !pool_is_thread_safe(effective_pool_mode())) {
        return 1;
    }
    return requested;
}

CallerPoolModeScope::CallerPoolModeScope(std::optional<unsigned int> mode,
                                         bool callbacks)
    : previous_mode_(caller_mode),
      previous_callbacks_(caller_callbacks_on_calling_thread) {
    if (mode) {
        caller_mode = mode;
    }
    caller_callbacks_on_calling_thread = callbacks;
}

CallerPoolModeScope::~CallerPoolModeScope() {
    caller_mode = previous_mode_;
    caller_callbacks_on_calling_thread = previous_callbacks_;
}

}  // namespace OECluster
