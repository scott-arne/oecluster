/**
 * @file KMedoids.cpp
 * @brief k-medoids (PAM) clustering implementation.
 */

#include "oecluster/clustering/KMedoids.h"

#include <algorithm>
#include <stdexcept>
#include <vector>

#include "DistanceAccess.h"
#include "KMedoidsSwapKernel.h"
#include "MaxMinKernel.h"
#include "oecluster/ThreadPool.h"

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

size_t global_medoid(const double* data, size_t n, size_t num_threads,
                     size_t chunk_size) {
    std::vector<double> sums(n, 0.0);

    ThreadPool pool(num_threads);
    pool.ParallelFor(0, n, detail::effective_chunk_size(n, chunk_size),
                     [&](size_t begin, size_t end) {
        for (size_t candidate = begin; candidate < end; ++candidate) {
            double sum = 0.0;
            for (size_t j = 0; j < n; ++j) {
                sum += detail::dense_distance(data, n, candidate, j);
            }
            sums[candidate] = sum;
        }
    });

    size_t best = 0;
    for (size_t candidate = 1; candidate < n; ++candidate) {
        if (sums[candidate] < sums[best]) {
            best = candidate;
        }
    }
    return best;
}

std::vector<size_t> build_initialize(const double* data, size_t n, size_t k,
                                     size_t num_threads, size_t chunk_size) {
    std::vector<size_t> medoids;
    medoids.reserve(k);

    const size_t first = global_medoid(data, n, num_threads, chunk_size);
    medoids.push_back(first);

    std::vector<bool> selected(n, false);
    selected[first] = true;

    std::vector<double> nearest(n);
    for (size_t j = 0; j < n; ++j) {
        nearest[j] = detail::dense_distance(data, n, j, first);
    }

    std::vector<double> gains(n, 0.0);
    while (medoids.size() < k) {
        ThreadPool pool(num_threads);
        pool.ParallelFor(0, n, detail::effective_chunk_size(n, chunk_size),
                         [&](size_t begin, size_t end) {
            for (size_t candidate = begin; candidate < end; ++candidate) {
                double gain = 0.0;
                for (size_t j = 0; j < n; ++j) {
                    const double distance =
                        detail::dense_distance(data, n, j, candidate);
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
        for (size_t j = 0; j < n; ++j) {
            const double distance = detail::dense_distance(data, n, j, best);
            if (distance < nearest[j]) {
                nearest[j] = distance;
            }
        }
    }

    return medoids;
}

std::vector<size_t> farthest_first_initialize(const double* data, size_t n,
                                              size_t k, size_t num_threads,
                                              size_t chunk_size) {
    // Seeding at the global medoid rather than item 0 costs one O(n^2) pass and
    // buys permutation invariance: an item-0 seed makes the whole clustering
    // depend on input row order.
    const size_t seed = global_medoid(data, n, num_threads, chunk_size);
    return detail::maxmin_select_from(data, n, k, seed);
}

KMedoidsResult assemble(const double* data, size_t n,
                        std::vector<size_t> medoids, size_t iterations,
                        bool converged, size_t num_threads, size_t chunk_size) {
    // Sorting makes the labeling canonical: label i always belongs to the i-th
    // smallest medoid index, so two results compare without a permutation in
    // between.
    std::sort(medoids.begin(), medoids.end());

    const std::vector<detail::Assignment> assignments =
        detail::build_assignments(data, n, medoids, num_threads, chunk_size);

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

// Takes the worker count as an argument rather than reading
// options.num_threads, so the cap applied in k_medoids_cluster cannot be
// bypassed by the initializers underneath.
std::vector<size_t> initialize_medoids(const double* data, size_t n,
                                       const KMedoidsOptions& options,
                                       size_t num_threads) {
    switch (options.init) {
        case KMedoidsInit::Build:
            return build_initialize(data, n, options.n_clusters, num_threads,
                                    options.chunk_size);
        case KMedoidsInit::FarthestFirst:
            return farthest_first_initialize(data, n, options.n_clusters,
                                             num_threads, options.chunk_size);
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
    const double* data = storage.Data();

    // Capped at the item count once here and then passed to every helper
    // below, so no path from this entry point can hand ThreadPool the raw
    // field. num_threads is a size_t on a public options struct, and uncapped
    // an oversized request is attempted rather than refused:
    // ThreadPool::ParallelFor asks the OS for exactly that many threads to run
    // what may be a single chunk. The mild outcome is a call that never
    // returns, and the severe one kills the process, because the vector being
    // filled already holds joinable threads when a creation failure throws and
    // destroying a joinable thread calls std::terminate -- leaving no exception
    // for a caller, or the Python layer above it, to report. A num_threads of 0
    // means "use the hardware concurrency" and has to keep meaning that;
    // validation guarantees n >= n_clusters >= 1, so the minimum passes that 0
    // through untouched and can never manufacture one from a nonzero request.
    const size_t num_threads = std::min(options.num_threads, n);

    if (options.n_clusters == n) {
        // Every item is already a medoid, so the candidate set is empty and
        // there is no swap for the verification pass to examine. Without the
        // short-circuit BUILD would spend O(n^3) reaching the same answer.
        return identity_partition(n);
    }

    std::vector<size_t> medoids =
        initialize_medoids(data, n, options, num_threads);

    std::vector<size_t> slot_of(n, n);
    detail::refresh_slot_map(slot_of, medoids);
    std::vector<detail::Assignment> assignments = detail::build_assignments(
        data, n, medoids, num_threads, options.chunk_size);
    double cost = detail::total_cost(assignments);

    size_t iterations = 0;
    bool converged = false;

    while (iterations < options.max_iterations) {
        const detail::SwapCandidate predicted = detail::best_predicted_swap(
            data, n, medoids, assignments, slot_of, num_threads,
            options.chunk_size);

        bool advanced = false;
        if (predicted.valid && predicted.score < 0.0) {
            // One slot changes, so one value is the whole undo state.
            const size_t displaced = medoids[predicted.leaving_slot];
            medoids[predicted.leaving_slot] = predicted.entering_item;
            detail::refresh_slot_map(slot_of, medoids);
            std::vector<detail::Assignment> trial = detail::build_assignments(
                data, n, medoids, num_threads, options.chunk_size);
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
            data, n, medoids, assignments, slot_of, cost, num_threads,
            options.chunk_size);
        if (!verified.valid) {
            converged = true;
            break;
        }

        medoids[verified.leaving_slot] = verified.entering_item;
        detail::refresh_slot_map(slot_of, medoids);
        assignments = detail::build_assignments(data, n, medoids, num_threads,
                                                options.chunk_size);
        cost = detail::total_cost(assignments);
        ++iterations;
    }

    // Reaching the cap skips the verification pass deliberately: the pass
    // exists to substantiate a local-optimality claim, and a capped run makes
    // no such claim.
    return assemble(data, n, medoids, iterations, converged, num_threads,
                    options.chunk_size);
}

}  // namespace OECluster
