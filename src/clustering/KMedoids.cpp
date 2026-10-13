/**
 * @file KMedoids.cpp
 * @brief k-medoids (PAM) clustering implementation.
 */

#include "oecluster/clustering/KMedoids.h"

#include <algorithm>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "DistanceAccess.h"
#include "DiversityValidation.h"
#include "KMedoidsBudget.h"
#include "KMedoidsDistanceSource.h"
#include "KMedoidsSwapKernel.h"
#include "MaxMinKernel.h"

namespace OECluster {

namespace {

constexpr const char* KMEDOIDS_NAME = "k_medoids";

// The matrix path pays its passes in memory reads rather than in comparisons,
// so it runs with no swap budget at all. A sentinel rather than a merely huge
// allowance: sweep_comparisons() saturates at SIZE_MAX for a vast item count,
// and a plain SIZE_MAX allowance would then be spent by its first subtraction.
constexpr size_t UNBOUNDED_SWAP_COMPARISONS = static_cast<size_t>(-1);

// Everything validate_options() checks that does not need the storage itself,
// so both entry points refuse the same options the same way. The out-of-range
// message still names the storage: it is the one wording callers already match
// on, and k-medoids has only ever had one meaning for it.
void validate_options_for_count(size_t n, const KMedoidsOptions& options) {
    if (options.chunk_size == 0) {
        throw std::invalid_argument("K-medoids chunk_size must be at least one");
    }
    if (options.max_iterations == 0) {
        throw std::invalid_argument("K-medoids max_iterations must be at least one");
    }
    if (options.n_clusters == 0) {
        throw std::invalid_argument("K-medoids n_clusters must be at least one");
    }
    if (options.n_clusters > n) {
        throw std::invalid_argument(
            "K-medoids n_clusters must be at most the item count");
    }

    if (options.init != KMedoidsInit::Explicit) {
        if (!options.initial_medoids.empty()) {
            throw std::invalid_argument(
                "K-medoids initial_medoids requires an explicit initialization");
        }
    } else {
        if (options.initial_medoids.size() != options.n_clusters) {
            throw std::invalid_argument(
                "K-medoids initial_medoids must hold exactly n_clusters indices");
        }
        // Range before uniqueness, in two passes rather than one: a list that
        // is both out of range and duplicated must report the range first,
        // whichever position each problem occupies.
        for (const size_t index : options.initial_medoids) {
            if (index >= n) {
                throw std::out_of_range(
                    "K-medoids initial_medoids index is outside the storage range");
            }
        }
        std::vector<bool> seen(n, false);
        for (const size_t index : options.initial_medoids) {
            if (seen[index]) {
                throw std::invalid_argument(
                    "K-medoids initial_medoids must be unique");
            }
            seen[index] = true;
        }
    }

    switch (options.init) {
        case KMedoidsInit::Build:
        case KMedoidsInit::FarthestFirst:
        case KMedoidsInit::Explicit:
            return;
    }

    throw std::invalid_argument("Unknown k-medoids initialization method");
}

void validate_options(const StorageBackend& storage,
                      const KMedoidsOptions& options) {
    detail::validate_complete_distance_storage(storage, "K-medoids clustering");
    validate_options_for_count(storage.NumSamples(), options);
}

template <typename Runner>
std::vector<size_t> build_initialize(Runner& distances, size_t k) {
    const size_t n = distances.Size();

    std::vector<size_t> medoids;
    medoids.reserve(k);

    const size_t first = detail::global_medoid(distances);
    medoids.push_back(first);

    std::vector<bool> selected(n, false);
    selected[first] = true;

    std::vector<double> nearest(n);
    distances.WithAccessor([&](const auto& distance) {
        for (size_t j = 0; j < n; ++j) {
            nearest[j] = distance(j, first);
        }
    });

    std::vector<double> gains(n, 0.0);
    while (medoids.size() < k) {
        distances.Run(n, [&](const auto& local, size_t begin, size_t end) {
            for (size_t candidate = begin; candidate < end; ++candidate) {
                double gain = 0.0;
                for (size_t j = 0; j < n; ++j) {
                    const double distance = local(j, candidate);
                    if (distance < nearest[j]) {
                        gain += nearest[j] - distance;
                    }
                }
                gains[candidate] = gain;
            }
        });

        // Ranking only unselected items is load-bearing, not an optimization:
        // an already-selected item gains nothing, and so does every item on an
        // all-zero matrix, so an unmasked scan would re-pick medoid one for
        // every remaining slot.
        size_t best = n;
        for (size_t candidate = 0; candidate < n; ++candidate) {
            if (selected[candidate]) {
                continue;
            }
            if (best == n || gains[candidate] > gains[best]) {
                best = candidate;
            }
        }

        selected[best] = true;
        medoids.push_back(best);
        distances.WithAccessor([&](const auto& distance) {
            for (size_t j = 0; j < n; ++j) {
                const double to_best = distance(j, best);
                if (to_best < nearest[j]) {
                    nearest[j] = to_best;
                }
            }
        });
    }

    return medoids;
}

template <typename Runner>
std::vector<size_t> farthest_first_initialize(Runner& distances, size_t k) {
    // Seeding at the global medoid rather than item 0 costs one O(n^2) pass and
    // buys permutation invariance: an item-0 seed makes the whole clustering
    // depend on input row order.
    const size_t seed = detail::global_medoid(distances);
    return detail::maxmin_select_from(distances, k, seed);
}

template <typename Runner>
KMedoidsResult assemble(Runner& distances, std::vector<size_t> medoids,
                        std::vector<double> columns, size_t iterations,
                        bool converged) {
    const size_t n = distances.Size();

    // Sorting makes the labeling canonical: label i always belongs to the i-th
    // smallest medoid index, so two results compare without a permutation in
    // between. The cached columns are slot-indexed, so sorting carries them.
    detail::sort_medoids_with_columns(medoids, columns);

    const std::vector<detail::Assignment> assignments =
        detail::build_assignments(distances, medoids, columns);

    std::vector<ClusterLabel> labels(n);
    for (size_t j = 0; j < n; ++j) {
        labels[j] = static_cast<ClusterLabel>(assignments[j].nearest_slot);
    }
    Clusters members = labels_to_clusters(labels);

    // Recomputed from the returned assignment rather than carried forward as a
    // running sum of accepted deltas: accumulated deltas drift, and a reported
    // cost that disagrees with the returned labels is the worst failure this
    // API could have.
    const double cost = detail::total_cost(assignments);

    return KMedoidsResult(std::move(labels), std::move(members),
                          std::move(medoids), cost, iterations, converged);
}

KMedoidsResult identity_partition(size_t n) {
    std::vector<size_t> medoids(n);
    std::vector<ClusterLabel> labels(n);
    for (size_t j = 0; j < n; ++j) {
        medoids[j] = j;
        labels[j] = static_cast<ClusterLabel>(j);
    }
    Clusters members = labels_to_clusters(labels);
    return KMedoidsResult(std::move(labels), std::move(members),
                          std::move(medoids), 0.0, 0, true);
}

// Reads the worker count from the runner rather than from options.num_threads,
// so the cap applied in k_medoids_cluster cannot be bypassed by the
// initializers underneath: the runner is built with the capped count and is
// the only thing here that starts a thread.
template <typename Runner>
std::vector<size_t> initialize_medoids(Runner& distances,
                                       const KMedoidsOptions& options) {
    switch (options.init) {
        case KMedoidsInit::Build:
            return build_initialize(distances, options.n_clusters);
        case KMedoidsInit::FarthestFirst:
            return farthest_first_initialize(distances, options.n_clusters);
        case KMedoidsInit::Explicit:
            return options.initial_medoids;
    }

    throw std::invalid_argument("Unknown k-medoids initialization method");
}

// The medoid, the slot map and the cached column all describe one placement.
// Changing any of them alone silently mis-assigns every item, so they move
// together and no call site is left to remember the other two.
template <typename Runner>
void place_medoid(Runner& distances, std::vector<size_t>& medoids,
                  std::vector<size_t>& slot_of, std::vector<double>& columns,
                  size_t slot, size_t item) {
    medoids[slot] = item;
    detail::refresh_slot_map(slot_of, medoids);
    detail::refresh_medoid_column(distances, columns, slot, item);
}

/**
 * @brief Charge one pass to the swap allowance, or report that it does not fit.
 *
 * Called between passes and never inside one: a stop part-way through a pass
 * would depend on where the workers happened to be, which would make the
 * result vary with the thread count.
 *
 * :param swap_allowance: Remaining allowance, decremented when the pass fits.
 * :param cost: Comparisons the pass about to run will make.
 * :returns: True when the pass may run.
 */
bool afford_pass(size_t& swap_allowance, size_t cost) {
    if (swap_allowance == UNBOUNDED_SWAP_COMPARISONS) {
        return true;
    }
    if (swap_allowance < cost) {
        return false;
    }
    swap_allowance -= cost;
    return true;
}

/**
 * @brief The whole algorithm over whichever distance source it is handed.
 *
 * :param distances: Distance source; see MatrixDistanceRunner.
 * :param options: Already-validated k-medoids options.
 * :param swap_allowance: Comparisons the swap loop may spend, or
 *     ``UNBOUNDED_SWAP_COMPARISONS`` for a run that is not counting.
 * :returns: Labels, clusters, medoids, cost, iteration count and convergence.
 */
template <typename Runner>
KMedoidsResult cluster_with(Runner& distances, const KMedoidsOptions& options,
                            size_t swap_allowance) {
    const size_t n = distances.Size();

    std::vector<size_t> medoids = initialize_medoids(distances, options);

    std::vector<size_t> slot_of(n, n);
    detail::refresh_slot_map(slot_of, medoids);
    std::vector<double> columns =
        detail::build_medoid_columns(distances, medoids);
    std::vector<detail::Assignment> assignments =
        detail::build_assignments(distances, medoids, columns);
    double cost = detail::total_cost(assignments);

    // The two passes the loop can run are priced apart because they cost
    // differently: the predicted-delta scan reads each pair once, while the
    // verification pass rescans every item for each of the k slots.
    const size_t sweep = detail::sweep_comparisons(n);
    const size_t verification =
        detail::verification_comparisons(n, medoids.size());
    size_t iterations = 0;
    bool converged = false;

    while (iterations < options.max_iterations) {
        // A budget stop leaves converged false, which is already what the loop
        // does when it exits without a verified optimum.
        if (!afford_pass(swap_allowance, sweep)) {
            break;
        }

        const detail::SwapCandidate predicted = detail::best_predicted_swap(
            distances, medoids, assignments, slot_of);

        bool advanced = false;
        if (predicted.valid && predicted.score < 0.0) {
            // One slot changes, so one value is the whole undo state.
            const size_t displaced = medoids[predicted.leaving_slot];
            place_medoid(distances, medoids, slot_of, columns,
                         predicted.leaving_slot, predicted.entering_item);
            std::vector<detail::Assignment> trial =
                detail::build_assignments(distances, medoids, columns);
            const double trial_cost = detail::total_cost(trial);

            if (trial_cost < cost) {
                // Every accepted configuration has a strictly smaller
                // recomputed cost than its predecessor, computed by the
                // identical expression in the identical order, so no
                // configuration can repeat and the loop terminates.
                assignments = std::move(trial);
                cost = trial_cost;
                ++iterations;
                advanced = true;
            } else {
                // Prediction and recomputation disagree at rounding scale, so
                // neither alone can say whether the loop is finished. Undo and
                // let the verification pass decide.
                place_medoid(distances, medoids, slot_of, columns,
                             predicted.leaving_slot, displaced);
            }
        }

        if (advanced) {
            continue;
        }

        // The second stop, for the dearer pass. Stopping here rather than
        // partway through it is what keeps an exhausted budget giving the
        // same answer at every thread count.
        if (!afford_pass(swap_allowance, verification)) {
            break;
        }

        const detail::SwapCandidate verified = detail::verification_pass(
            distances, medoids, assignments, slot_of, cost);
        if (!verified.valid) {
            converged = true;
            break;
        }

        place_medoid(distances, medoids, slot_of, columns,
                     verified.leaving_slot, verified.entering_item);
        assignments = detail::build_assignments(distances, medoids, columns);
        cost = detail::total_cost(assignments);
        ++iterations;
    }

    // Reaching the cap skips the verification pass deliberately: the pass
    // exists to substantiate a local-optimality claim, and a capped run makes
    // no such claim.
    return assemble(distances, std::move(medoids), std::move(columns),
                    iterations, converged);
}

}  // namespace

KMedoidsResult k_medoids_cluster(const StorageBackend& storage,
                                 const KMedoidsOptions& options) {
    validate_options(storage, options);

    const size_t n = storage.NumSamples();

    // Capped at the item count once here and then passed to every helper
    // below, so no path from this entry point can hand ThreadPool the raw
    // field. num_threads is a size_t on a public options struct, so an
    // oversized request is attempted rather than refused. The pool itself now
    // absorbs the two ways that used to go badly: detail::launch_workers
    // starts at most one worker per chunk, and detail::launch_threads joins
    // every thread it managed to start before rethrowing a spawn failure. What
    // the cap still buys is that an extreme value never reaches the pool at
    // all, which keeps the resource cost bounded by the problem rather than by
    // what the caller typed. A num_threads of 0 means "use the hardware
    // concurrency" and has to keep meaning that; validation guarantees
    // n >= n_clusters >= 1, so the minimum passes that 0 through untouched and
    // can never manufacture one from a nonzero request.
    const size_t num_threads = std::min(options.num_threads, n);

    if (options.n_clusters == n) {
        // Every item is already a medoid, so the candidate set is empty and
        // there is no swap for the verification pass to examine. Without the
        // short-circuit BUILD would spend O(n^3) reaching the same answer.
        return identity_partition(n);
    }

    // Built once with the capped worker count, so every pass below is
    // partitioned and scheduled the same way whichever kernel asks for it.
    detail::MatrixDistanceRunner distances(storage.Data(), n, num_threads,
                                           options.chunk_size);
    return cluster_with(distances, options, UNBOUNDED_SWAP_COMPARISONS);
}

KMedoidsResult k_medoids_cluster(PairwiseComparison& comparison,
                                 const KMedoidsOptions& options) {
    const size_t n = comparison.Size();
    // k-medoids sums distances across pairs to form its objective, so a
    // subset-scored comparison is refused here where cluster_report accepts
    // one: its per-pair scores are not mutually comparable, let alone
    // summable. Neither validator asks about the triangle inequality, so this
    // leaves non-metric input admissible, as the matrix overload has it.
    detail::validate_comparison_facts(comparison, KMEDOIDS_NAME);
    detail::refuse_unrepeatable(
        comparison, KMEDOIDS_NAME, "cluster",
        "the same pair would score differently on each of the many passes "
        "k-medoids makes over it");
    validate_options_for_count(n, options);

    if (options.n_clusters == n) {
        // Matches the storage overload: every item is already a medoid, so
        // there is nothing to score. Checked before the budget so a trivial
        // call is never refused for work it will not do.
        return identity_partition(n);
    }

    const size_t budget = options.max_comparisons == 0
                              ? detail::KMEDOIDS_DEFAULT_MAX_COMPARISONS
                              : options.max_comparisons;
    const size_t unavoidable =
        detail::unavoidable_comparisons(n, options.n_clusters, options.init);
    if (unavoidable > budget) {
        // The lever is priced for the n_clusters actually asked for rather
        // than described by a rule of thumb, because no rule of thumb is true
        // across the range: Build's floor is 2 * n_clusters passes over every
        // pair against farthest-first's n_clusters + 1, a saving of nothing
        // at one medoid, a quarter at two, and under half however large
        // n_clusters grows -- the verification pass underneath is what
        // neither can shed. Naming the other floor also answers the question
        // the caller actually has, which is whether switching would clear
        // this limit. Suppressed when it would not help at all: at
        // n_clusters == 1 the two floors are equal, and the advised switch
        // would return this very message unchanged.
        std::string lever;
        if (options.init == KMedoidsInit::Build) {
            const size_t farthest_first = detail::unavoidable_comparisons(
                n, options.n_clusters, KMedoidsInit::FarthestFirst);
            if (farthest_first < unavoidable) {
                lever = " init=farthest_first would need " +
                        std::to_string(farthest_first) + " instead, " +
                        (farthest_first <= budget ? "which fits within this limit."
                                                  : "still over it.");
            }
        }
        throw std::invalid_argument(
            std::string(KMEDOIDS_NAME) + " needs at least " +
            std::to_string(unavoidable) +
            " pair comparisons before any swap, over the max_comparisons "
            "limit of " + std::to_string(budget) + "." + lever +
            " Raise max_comparisons, or cluster a matrix from pdist() "
            "instead");
    }

    // Capped exactly as the storage overload caps it, and for the same
    // reasons; see the comment there.
    detail::StreamingDistanceRunner distances(comparison, KMEDOIDS_NAME,
                                              std::min(options.num_threads, n),
                                              options.chunk_size);
    return cluster_with(distances, options, budget - unavoidable);
}

}  // namespace OECluster
