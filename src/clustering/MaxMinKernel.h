/**
 * @file MaxMinKernel.h
 * @brief Deterministic farthest-first (MaxMin) selection over a row provider.
 */

#ifndef OECLUSTER_SRC_CLUSTERING_MAXMINKERNEL_H
#define OECLUSTER_SRC_CLUSTERING_MAXMINKERNEL_H

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "DistanceAccess.h"
// For the distance source shape. k-medoids initialization reaches these
// helpers through a runner, so they have to read the same way its kernels do.
#include "KMedoidsDistanceSource.h"

namespace OECluster::detail {

/** @brief Why a MaxMin run stopped; the public MaxMinStop mirrors it. */
enum class MaxMinKernelStop { Count, Threshold, Exhausted };

/**
 * @brief What one MaxMin run is asked to do.
 *
 * The kernel trusts its caller: ``initial`` is non-empty, distinct and in
 * range, ``count`` is 0 or in [initial.size(), n], and at least one of
 * ``count`` and ``threshold`` is set. Every public entry point validates
 * those before it builds a request.
 */
struct MaxMinKernelRequest {
    size_t count = 0;  ///< Total selection size, initial included; 0 = no limit.
    double threshold = std::numeric_limits<double>::quiet_NaN();  ///< NaN = unset.
    std::vector<size_t> initial;  ///< Selected before the first pick, in order.
};

/** @brief The selection one MaxMin run produced. */
struct MaxMinKernelResult {
    std::vector<size_t> indices;         ///< Selection order.
    std::vector<double> pick_distances;  ///< NaN for every initial entry.
    MaxMinKernelStop stop = MaxMinKernelStop::Exhausted;
};

/**
 * @brief The error every refusing provider throws for a NaN or infinite read.
 *
 * :param a: First item of the pair that was read.
 * :param b: Second item of the pair that was read.
 * :returns: An exception naming the pair.
 */
inline std::invalid_argument non_finite_distance_error(size_t a, size_t b) {
    return std::invalid_argument(
        "Diversity selection read a non-finite distance between items " +
        std::to_string(a) + " and " + std::to_string(b));
}

/**
 * @brief The fold every row provider performs, over whatever reads a distance.
 *
 * ``first`` assigns rather than improves. That is what keeps a NaN from the
 * first row in ``nearest`` -- no later comparison against it is true -- exactly
 * as the pre-refactor kernel behaved; initializing to +infinity and taking a
 * minimum would silently change that.
 *
 * :param distance: Callable taking two item indices and returning a distance.
 * :param n: Number of items.
 * :param p: Item whose row is folded; already marked selected.
 * :param selected: Items excluded from the fold.
 * :param nearest: Per-item distance to the selection, updated in place.
 * :param first: Assign instead of applying the strict-improvement rule.
 * :param refuse_non_finite: Refuse rather than fold a NaN or infinite read.
 * :raises std::invalid_argument: If refuse_non_finite is set and a read
 *     distance is NaN or infinite.
 */
template <typename Accessor>
inline void fold_row_into(const Accessor& distance, size_t n, size_t p,
                          const std::vector<bool>& selected,
                          std::vector<double>& nearest, bool first,
                          bool refuse_non_finite) {
    for (size_t j = 0; j < n; ++j) {
        if (selected[j]) {
            continue;
        }
        const double value = distance(p, j);
        if (refuse_non_finite && !std::isfinite(value)) {
            throw non_finite_distance_error(p, j);
        }
        if (first || value < nearest[j]) {
            nearest[j] = value;
        }
    }
}

/**
 * @brief Row provider over a condensed dense distance array.
 *
 * Serial, like the kernel it replaced. ``refuse_non_finite`` is off only for
 * k-medoids: its FarthestFirst initialization never validated finiteness, and
 * moving it onto this kernel must not change what a native caller observes.
 */
class MatrixRowProvider {
public:
    MatrixRowProvider(const double* data, size_t n, bool refuse_non_finite)
        : data_(data), n_(n), refuse_non_finite_(refuse_non_finite) {}

    /** @brief Fold the distances from item p into nearest; see fold_row_into. */
    void FoldRow(size_t p, const std::vector<bool>& selected,
                 std::vector<double>& nearest, bool first) const {
        fold_row_into(MatrixDistances(data_, n_), n_, p, selected, nearest,
                      first, refuse_non_finite_);
    }

private:
    const double* data_;
    size_t n_;
    bool refuse_non_finite_;
};

/**
 * @brief Row provider over a distance source runner.
 *
 * Serial, so it takes one accessor per row rather than holding one: a
 * streaming source leases a comparison clone for the duration, and holding a
 * lease across the whole selection would keep it out of circulation while the
 * kernel is between rows.
 */
template <typename Runner>
class RunnerRowProvider {
public:
    RunnerRowProvider(Runner& distances, bool refuse_non_finite)
        : distances_(distances), refuse_non_finite_(refuse_non_finite) {}

    /** @brief Fold the distances from item p into nearest; see fold_row_into. */
    void FoldRow(size_t p, const std::vector<bool>& selected,
                 std::vector<double>& nearest, bool first) const {
        const size_t n = distances_.Size();
        distances_.WithAccessor([&](const auto& distance) {
            fold_row_into(distance, n, p, selected, nearest, first,
                          refuse_non_finite_);
        });
    }

private:
    Runner& distances_;
    bool refuse_non_finite_;
};

/**
 * @brief Deterministic farthest-first (MaxMin) selection over a row provider.
 *
 * Ties resolve to the smaller item index at every step, so the result depends
 * only on the distances and the request, never on thread count or iteration
 * order.
 *
 * Already-selected items are excluded from every subsequent pick. That is not
 * an optimization: on a matrix where every distance is equal, an unmasked scan
 * would tie every candidate at zero and the smaller-index rule would return
 * the seed over and over.
 *
 * Every initial entry is marked selected before any row is folded, and each
 * pick is marked before its own row is folded, so no provider is ever asked
 * for a diagonal distance or for a pair between two initial entries. The row
 * of the pick that completes ``count`` is never folded: nothing reads it.
 *
 * Stop checks run after the initial set and after every pick, in the order
 * count, exhausted, threshold. A candidate at or within the threshold is not
 * added, so with no initial set the selection is a strict packing.
 *
 * :param rows: Row provider; see MatrixRowProvider::FoldRow.
 * :param n: Number of items.
 * :param request: A request that satisfies MaxMinKernelRequest's contract.
 * :returns: The selection, its pick distances and the stop reason.
 */
template <typename RowProvider>
MaxMinKernelResult maxmin_run(RowProvider& rows, size_t n,
                              const MaxMinKernelRequest& request) {
    const double unset = std::numeric_limits<double>::quiet_NaN();
    const bool has_count = request.count != 0;
    const bool has_threshold = !std::isnan(request.threshold);

    MaxMinKernelResult result;
    result.indices.reserve(has_count ? request.count : n);
    result.pick_distances.reserve(has_count ? request.count : n);

    std::vector<bool> selected(n, false);
    for (const size_t index : request.initial) {
        selected[index] = true;
        result.indices.push_back(index);
        result.pick_distances.push_back(unset);
    }
    if (has_count && result.indices.size() == request.count) {
        result.stop = MaxMinKernelStop::Count;
        return result;
    }

    std::vector<double> nearest(n, 0.0);
    for (size_t position = 0; position < request.initial.size(); ++position) {
        rows.FoldRow(request.initial[position], selected, nearest, position == 0);
    }

    while (true) {
        // Ascending scan with a strict improvement is the tie rule: the
        // smallest index wins any tie because a later equal never displaces it.
        size_t best = n;
        for (size_t j = 0; j < n; ++j) {
            if (selected[j]) {
                continue;
            }
            if (best == n || nearest[j] > nearest[best]) {
                best = j;
            }
        }
        if (best == n) {
            result.stop = MaxMinKernelStop::Exhausted;
            return result;
        }
        if (has_threshold && nearest[best] <= request.threshold) {
            result.stop = MaxMinKernelStop::Threshold;
            return result;
        }

        selected[best] = true;
        result.indices.push_back(best);
        result.pick_distances.push_back(nearest[best]);
        if (has_count && result.indices.size() == request.count) {
            result.stop = MaxMinKernelStop::Count;
            return result;
        }
        rows.FoldRow(best, selected, nearest, false);
    }
}

/**
 * @brief Deterministic farthest-first selection over any row provider.
 *
 * Does not validate finiteness; the provider decides that.
 *
 * :param rows: Row provider; see MatrixRowProvider::FoldRow.
 * :param n: Number of items.
 * :param count: Number of items to select; must be in [1, n].
 * :param seed: Index of the first selection; must be less than n.
 * :returns: ``count`` item indices in selection order.
 * :raises std::invalid_argument: If count is outside [1, n] or seed is not a
 *     valid item index.
 */
template <typename RowProvider>
inline std::vector<size_t> maxmin_select_rows(RowProvider& rows, size_t n,
                                              size_t count, size_t seed) {
    if (count == 0 || count > n) {
        throw std::invalid_argument(
            "MaxMin selection count must be in [1, item count]");
    }
    if (seed >= n) {
        throw std::invalid_argument(
            "MaxMin selection seed is outside the item range");
    }

    MaxMinKernelRequest request;
    request.count = count;
    request.initial = {seed};
    return maxmin_run(rows, n, request).indices;
}

/**
 * @brief Deterministic farthest-first selection over a distance source.
 *
 * k-medoids FarthestFirst initialization. Serial, so it never runs a pass on
 * the runner; it only borrows an accessor per row.
 *
 * :param distances: Distance source; see MatrixDistanceRunner.
 * :param count: Number of items to select; must be in [1, n].
 * :param seed: Index of the first selection; must be less than n.
 * :returns: ``count`` item indices in selection order.
 */
template <typename Runner>
inline std::vector<size_t> maxmin_select_from(Runner& distances, size_t count,
                                              size_t seed) {
    RunnerRowProvider<Runner> rows(distances, /*refuse_non_finite=*/false);
    return maxmin_select_rows(rows, distances.Size(), count, seed);
}

/**
 * @brief Deterministic farthest-first selection over a condensed dense array.
 *
 * The pre-refactor entry point, kept as an adapter so its white-box tests keep
 * pinning the kernel. It does not validate finiteness.
 *
 * :param data: Condensed distance array from ``StorageBackend::Data()``.
 * :param n: Number of items.
 * :param count: Number of items to select; must be in [1, n].
 * :param seed: Index of the first selection; must be less than n.
 * :returns: ``count`` item indices in selection order.
 */
inline std::vector<size_t> maxmin_select_from(const double* data, size_t n,
                                              size_t count, size_t seed) {
    MatrixRowProvider rows(data, n, /*refuse_non_finite=*/false);
    return maxmin_select_rows(rows, n, count, seed);
}

/**
 * @brief Each item's sum of distances to every item, the diagonal's 0 included.
 *
 * Neither finiteness nor overflow is checked; a caller that needs either
 * checks the returned sums or the input itself.
 *
 * :param distances: Distance source; see MatrixDistanceRunner.
 * :returns: One sum per item.
 */
template <typename Runner>
inline std::vector<double> medoid_row_sums(Runner& distances) {
    const size_t n = distances.Size();
    std::vector<double> sums(n, 0.0);

    distances.Run(n, [&](const auto& local, size_t begin, size_t end) {
        for (size_t candidate = begin; candidate < end; ++candidate) {
            double sum = 0.0;
            for (size_t j = 0; j < n; ++j) {
                sum += local(candidate, j);
            }
            sums[candidate] = sum;
        }
    });
    return sums;
}

/**
 * @brief Index of the smallest sum, ties to the smaller index.
 *
 * :param sums: Per-item sums from medoid_row_sums; non-empty.
 * :returns: The global medoid.
 */
inline size_t argmin_row_sum(const std::vector<double>& sums) {
    size_t best = 0;
    for (size_t candidate = 1; candidate < sums.size(); ++candidate) {
        if (sums[candidate] < sums[best]) {
            best = candidate;
        }
    }
    return best;
}

/**
 * @brief The item with the smallest distance sum, ties to the smaller index.
 *
 * :param distances: Distance source; see MatrixDistanceRunner.
 * :returns: The global medoid.
 */
template <typename Runner>
inline size_t global_medoid(Runner& distances) {
    return argmin_row_sum(medoid_row_sums(distances));
}

/**
 * @brief The global medoid of a condensed dense distance array.
 *
 * :param data: Condensed distance array from ``StorageBackend::Data()``.
 * :param n: Number of items; at least one.
 * :param num_threads: Worker threads; 0 auto-detects hardware concurrency.
 * :param chunk_size: Items per work unit; at least one.
 * :returns: The global medoid.
 */
inline size_t global_medoid(const double* data, size_t n, size_t num_threads,
                            size_t chunk_size) {
    MatrixDistanceRunner distances(data, n, num_threads, chunk_size);
    return global_medoid(distances);
}

}  // namespace OECluster::detail

#endif  // OECLUSTER_SRC_CLUSTERING_MAXMINKERNEL_H
