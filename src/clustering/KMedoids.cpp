/**
 * @file KMedoids.cpp
 * @brief k-medoids (PAM) clustering implementation.
 */

#include "oecluster/clustering/KMedoids.h"

#include <algorithm>
#include <stdexcept>
#include <vector>

#include "DistanceAccess.h"
#include "KMedoidsDistanceSource.h"
#include "KMedoidsSwapKernel.h"
#include "MaxMinKernel.h"

namespace OECluster {

namespace {

void validate_options(const StorageBackend& storage,
                      const KMedoidsOptions& options) {
    detail::validate_complete_distance_storage(storage, "K-medoids clustering");

    if (options.chunk_size == 0) {
        throw std::invalid_argument("K-medoids chunk_size must be at least one");
    }
    if (options.max_iterations == 0) {
        throw std::invalid_argument("K-medoids max_iterations must be at least one");
    }
    if (options.n_clusters == 0) {
        throw std::invalid_argument("K-medoids n_clusters must be at least one");
    }
    if (options.n_clusters > storage.NumSamples()) {
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
            if (index >= storage.NumSamples()) {
                throw std::out_of_range(
                    "K-medoids initial_medoids index is outside the storage range");
            }
        }
        std::vector<bool> seen(storage.NumSamples(), false);
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
                        size_t iterations, bool converged) {
    const size_t n = distances.Size();

    // Sorting makes the labeling canonical: label i always belongs to the i-th
    // smallest medoid index, so two results compare without a permutation in
    // between.
    std::sort(medoids.begin(), medoids.end());

    const std::vector<detail::Assignment> assignments =
        detail::build_assignments(distances, medoids);

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

}  // namespace

KMedoidsResult k_medoids_cluster(const StorageBackend& storage,
                                 const KMedoidsOptions& options) {
    validate_options(storage, options);

    const size_t n = storage.NumSamples();

    // Capped at the item count once here and then passed to every helper
    // below, so no path from this entry point can hand ThreadPool the raw
    // field. num_threads is a size_t on a public options struct, and uncapped
    // an oversized request is attempted rather than refused: ParallelFor
    // reserves a vector of that many threads and then fills it, to run what may
    // be a single chunk. Which way that goes wrong depends on the magnitude.
    // Past the vector's max_size the reserve throws length_error before any
    // thread exists, which is the benign case because the exception escapes
    // normally. Below that the request is genuinely attempted, so a merely
    // enormous count may consume extreme resources or fail to return promptly;
    // measured here, 100000 threads for a six-item problem did complete, while
    // 2^20 had not returned after ten minutes. Worst,
    // a creation that fails partway unwinds a vector still holding joinable
    // threads, and destroying a joinable thread calls std::terminate -- leaving
    // no exception for a caller, or the Python layer above it, to report. Only
    // the first of the three is what the regression tests observe; the cap
    // exists for the other two. A num_threads of 0 means "use the hardware
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

    std::vector<size_t> medoids = initialize_medoids(distances, options);

    std::vector<size_t> slot_of(n, n);
    detail::refresh_slot_map(slot_of, medoids);
    std::vector<detail::Assignment> assignments =
        detail::build_assignments(distances, medoids);
    double cost = detail::total_cost(assignments);

    size_t iterations = 0;
    bool converged = false;

    while (iterations < options.max_iterations) {
        const detail::SwapCandidate predicted = detail::best_predicted_swap(
            distances, medoids, assignments, slot_of);

        bool advanced = false;
        if (predicted.valid && predicted.score < 0.0) {
            // One slot changes, so one value is the whole undo state.
            const size_t displaced = medoids[predicted.leaving_slot];
            medoids[predicted.leaving_slot] = predicted.entering_item;
            detail::refresh_slot_map(slot_of, medoids);
            std::vector<detail::Assignment> trial =
                detail::build_assignments(distances, medoids);
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
                medoids[predicted.leaving_slot] = displaced;
                detail::refresh_slot_map(slot_of, medoids);
            }
        }

        if (advanced) {
            continue;
        }

        const detail::SwapCandidate verified = detail::verification_pass(
            distances, medoids, assignments, slot_of, cost);
        if (!verified.valid) {
            converged = true;
            break;
        }

        medoids[verified.leaving_slot] = verified.entering_item;
        detail::refresh_slot_map(slot_of, medoids);
        assignments = detail::build_assignments(distances, medoids);
        cost = detail::total_cost(assignments);
        ++iterations;
    }

    // Reaching the cap skips the verification pass deliberately: the pass
    // exists to substantiate a local-optimality claim, and a capped run makes
    // no such claim.
    return assemble(distances, medoids, iterations, converged);
}

}  // namespace OECluster
