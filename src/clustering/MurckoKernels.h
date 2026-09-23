/**
 * @file MurckoKernels.h
 * @brief Internal extraction kernel for Murcko scaffold assignment.
 *
 * src-private, following BitBirchKernels.h and KMedoidsSwapKernel.h: not
 * installed and not in the umbrella header, but included directly by the C++
 * tests so extraction correctness is unit-tested without expanding the public
 * surface.
 */

#ifndef OECLUSTER_CLUSTERING_MURCKO_KERNELS_H
#define OECLUSTER_CLUSTERING_MURCKO_KERNELS_H

#include <algorithm>
#include <cstddef>
#include <optional>
#include <string>
#include <thread>
#include <utility>
#include <vector>

#include <oesystem.h>

#include "oecluster/Error.h"
#include "oecluster/ThreadPool.h"
#include "oecluster/clustering/MurckoScaffold.h"

namespace OEChem {
class OEMolBase;
}

namespace OECluster::detail {

/**
 * @brief Extract one molecule's Bemis-Murcko scaffold.
 *
 * Never throws: the three outcomes are distinct in
 * the return type, which is what lets the threaded caller report the *first*
 * failing index in input order rather than whichever worker lost the race.
 *
 * :param mol: Molecule to read; it is copied and never modified.
 * :param type: Extraction level.
 * :returns: The scaffold SMILES, ``""`` for a molecule with no ring bonds, or
 *     ``std::nullopt`` when a ring-containing molecule yielded no framework or
 *     any SDK transformation step reported failure.
 */
std::optional<std::string> scaffold_of(const OEChem::OEMolBase& mol,
                                       ScaffoldType type) noexcept;

/**
 * @brief Run an index-aware extraction callable over ``count`` items.
 *
 * ``extract`` is a parameter rather than a hard-wired call to ``scaffold_of``
 * so that tests can force failures at chosen indices. No real molecule is known
 * to make ``OEGetBemisMurcko`` return nothing for a ring-containing input, so
 * without this seam the first-failing-index guarantee would ship with no
 * regression coverage at all.
 *
 * Each result is written to its own slot of a pre-sized buffer. This function
 * makes no decision about failures.
 *
 * :param count: Number of items.
 * :param num_threads: Worker threads; 0 auto-detects. Already clamped.
 * :param chunk: Items per work unit; 0 and over-large values are corrected.
 * :param extract: Callable taking an index and returning the slot's value. It
 *     is copied once and then invoked concurrently by every worker, so it must
 *     be safe to call from several threads at once -- a functor carrying
 *     mutable state would race.
 * :returns: One entry per item, in input order.
 */
template <typename Extract>
std::vector<std::optional<std::string>> extract_all(size_t count,
                                                    size_t num_threads,
                                                    size_t chunk,
                                                    Extract extract) {
    std::vector<std::optional<std::string>> raw(count);
    if (count == 0u) {
        return raw;
    }
    // ParallelFor early-returns on a zero chunk size, and an over-large value
    // overflows its chunk-count arithmetic into zero chunks. Either way no work
    // runs and every slot stays empty, so the failure scan then blames index 0
    // -- a confidently wrong answer rather than a crash.
    const size_t safe_chunk = std::min(chunk == 0u ? count : chunk, count);
    ThreadPool pool(num_threads);
    pool.ParallelFor(0u, count, safe_chunk, [&](size_t begin, size_t end) {
        for (size_t i = begin; i < end; ++i) {
            raw[i] = extract(i);
        }
    });
    return raw;
}

/** @brief Index of the first failed extraction in input order, else nullopt. */
inline std::optional<size_t> first_failure(
        const std::vector<std::optional<std::string>>& raw) {
    for (size_t i = 0; i < raw.size(); ++i) {
        if (!raw[i].has_value()) {
            return i;
        }
    }
    return std::nullopt;
}

/**
 * @brief Convert raw extraction slots into scaffolds, or report the first failure.
 *
 * Separated from murcko_scaffolds so that the failure conversion -- the
 * exception type, the index it names, and the message text -- is directly
 * testable. No real molecule is known to make the kernel fail, so without this
 * seam the branch would ship with no coverage at all.
 *
 * :param raw: Extraction slots in input order. Consumed: engaged slots are
 *     moved from.
 * :param caller: Name of the public entry point, used in the diagnostic.
 * :returns: One scaffold per slot, in input order.
 * :raises ComparisonError: If any slot is disengaged, naming the first such
 *     index in input order rather than whichever worker failed first.
 */
inline std::vector<std::string> finish_extraction(
        std::vector<std::optional<std::string>>&& raw, const char* caller) {
    if (const std::optional<size_t> bad = first_failure(raw)) {
        throw ComparisonError(std::string(caller) +
                              " could not extract a scaffold for molecule at index " +
                              std::to_string(*bad));
    }
    std::vector<std::string> scaffolds;
    scaffolds.reserve(raw.size());
    for (std::optional<std::string>& value : raw) {
        scaffolds.push_back(std::move(*value));
    }
    return scaffolds;
}

/**
 * @brief Whether a memory-pool mode permits driving the toolkit from threads.
 *
 * A pure function over the mode word rather than an inline expression, so the
 * gate is testable against every bit combination without mutating the
 * process-global setting.
 */
inline bool pool_is_thread_safe(unsigned int mode) {
    return (mode & (OESystem::OEMemPoolMode::ThreadLocal |
                    OESystem::OEMemPoolMode::Mutexed)) != 0u;
}

/**
 * @brief Clamp a requested worker count to what the work and machine support.
 *
 * ParallelFor reserves a thread vector of the requested size and then fills it,
 * so a failed thread creation partway unwinds a vector still holding joinable
 * threads and calls std::terminate. ``num_threads`` is an unbounded size_t on a
 * public options struct, so that path is reachable from Python; clamping by
 * item count alone does not close it, because that bound only binds when the
 * input is smaller than the request.
 */
inline size_t effective_thread_count(size_t requested, size_t n) {
    if (requested == 0u) {
        return 0u;  // 0 means "auto-detect"; ThreadPool resolves it.
    }
    const size_t hw = std::max<size_t>(1u, std::thread::hardware_concurrency());
    // Extraction is CPU-bound, so oversubscribing far past the core count buys
    // nothing. The multiple rather than hw exactly leaves room for a caller who
    // knows their machine.
    //
    // The max guards the n == 0 case, where the mins would otherwise return 0
    // and collide with the auto-detect sentinel above. For any n >= 1 the mins
    // already yield at least 1, so nothing reachable changes.
    return std::max<size_t>(1u, std::min(requested, std::min(n, 4u * hw)));
}

/**
 * @brief Resolve the worker count for one extraction run.
 *
 * Joins the two independent decisions -- whether the toolkit's memory pool is
 * safe to drive from threads, and how many workers this request justifies --
 * so that the join is a pure function a test can exercise at any mode.
 * Mutating the process-global mode is not an option: OESetMemPoolMode is fatal
 * on a second call, so a test that changed it could never restore it.
 *
 * :param mode: The value OEGetMemPoolMode() reported.
 * :param requested: The caller's MurckoOptions::num_threads.
 * :param n: Number of molecules.
 * :returns: 1 when the pool is not thread-safe, otherwise the clamped count.
 */
inline size_t dispatch_thread_count(unsigned int mode, size_t requested, size_t n) {
    if (!pool_is_thread_safe(mode)) {
        return 1u;
    }
    return effective_thread_count(requested, n);
}

}  // namespace OECluster::detail

#endif  // OECLUSTER_CLUSTERING_MURCKO_KERNELS_H
