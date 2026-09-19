/**
 * @file KMedoids.cpp
 * @brief k-medoids (PAM) clustering implementation.
 */

#include "oecluster/clustering/KMedoids.h"

#include <algorithm>
#include <limits>
#include <stdexcept>
#include <vector>

#include "DistanceAccess.h"
#include "MaxMinKernel.h"
#include "oecluster/ThreadPool.h"

namespace OECluster {

namespace {

/**
 * @brief Cached nearest and second-nearest medoid facts for one item.
 */
struct Assignment {
    size_t nearest_slot = 0;          ///< Slot index, not an item index.
    double nearest_distance = 0.0;
    double second_nearest_distance =  ///< Where the item lands if its own medoid leaves.
        std::numeric_limits<double>::infinity();
};

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

void refresh_slot_map(std::vector<size_t>& slot_of,
                      const std::vector<size_t>& medoids) {
    std::fill(slot_of.begin(), slot_of.end(), slot_of.size());
    for (size_t slot = 0; slot < medoids.size(); ++slot) {
        slot_of[medoids[slot]] = slot;
    }
}

// ThreadPool derives its chunk count as ``(range + chunk_size - 1) /
// chunk_size``, which wraps to zero for a chunk size near SIZE_MAX and then
// runs no chunk at all: every scan below would return its default-initialized
// output and the call would report a zero-cost single cluster as converged.
// Clamping to the item count is observationally free -- a chunk at least as
// wide as the range is one chunk either way -- and keeps both this file's own
// ceiling arithmetic and the pool's below the overflow. Validation guarantees
// ``chunk_size >= 1``, and every call site has ``n >= n_clusters >= 1``, so the
// result is never zero.
size_t effective_chunk_size(size_t n, size_t chunk_size) {
    return std::min(chunk_size, n);
}

std::vector<Assignment> build_assignments(const double* data, size_t n,
                                          const std::vector<size_t>& medoids,
                                          size_t num_threads,
                                          size_t chunk_size) {
    const size_t k = medoids.size();
    const double infinity = std::numeric_limits<double>::infinity();

    std::vector<size_t> slot_of(n, n);
    refresh_slot_map(slot_of, medoids);

    std::vector<Assignment> assignments(n);
    ThreadPool pool(num_threads);
    pool.ParallelFor(0, n, effective_chunk_size(n, chunk_size),
                     [&](size_t begin, size_t end) {
        for (size_t j = begin; j < end; ++j) {
            Assignment entry;

            if (slot_of[j] != n) {
                // The self-assignment rule: a medoid always belongs to its own
                // slot. Ranking it against the others would hand two duplicate
                // medoids at distance 0 to the same slot and leave the other
                // cluster empty, breaking the exactly-k guarantee.
                entry.nearest_slot = slot_of[j];
                entry.nearest_distance = 0.0;
                entry.second_nearest_distance = infinity;
                for (size_t slot = 0; slot < k; ++slot) {
                    if (slot == entry.nearest_slot) {
                        continue;
                    }
                    const double distance =
                        detail::dense_distance(data, n, j, medoids[slot]);
                    if (distance < entry.second_nearest_distance) {
                        entry.second_nearest_distance = distance;
                    }
                }
            } else {
                entry.nearest_slot = 0;
                entry.nearest_distance = infinity;
                entry.second_nearest_distance = infinity;
                for (size_t slot = 0; slot < k; ++slot) {
                    const double distance =
                        detail::dense_distance(data, n, j, medoids[slot]);
                    // Slots are visited in slot order, but the tie rule ranks
                    // on the medoid's item index, so an equal distance only
                    // displaces the incumbent when its item index is smaller.
                    const bool wins =
                        distance < entry.nearest_distance ||
                        (distance == entry.nearest_distance &&
                         medoids[slot] < medoids[entry.nearest_slot]);
                    if (wins) {
                        entry.second_nearest_distance = entry.nearest_distance;
                        entry.nearest_distance = distance;
                        entry.nearest_slot = slot;
                    } else if (distance < entry.second_nearest_distance) {
                        entry.second_nearest_distance = distance;
                    }
                }
            }

            assignments[j] = entry;
        }
    });

    return assignments;
}

double total_cost(const std::vector<Assignment>& assignments) {
    // Ascending item order, so the sum is bit-reproducible whatever the
    // chunking that filled the cache.
    double cost = 0.0;
    for (const Assignment& entry : assignments) {
        cost += entry.nearest_distance;
    }
    return cost;
}

size_t global_medoid(const double* data, size_t n, size_t num_threads,
                     size_t chunk_size) {
    std::vector<double> sums(n, 0.0);

    ThreadPool pool(num_threads);
    pool.ParallelFor(0, n, effective_chunk_size(n, chunk_size),
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
        pool.ParallelFor(0, n, effective_chunk_size(n, chunk_size),
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

    const std::vector<Assignment> assignments =
        build_assignments(data, n, medoids, num_threads, chunk_size);

    std::vector<ClusterLabel> labels(n);
    for (size_t j = 0; j < n; ++j) {
        labels[j] = static_cast<ClusterLabel>(assignments[j].nearest_slot);
    }
    Clusters members = labels_to_clusters(labels);

    // Recomputed from the returned assignment rather than carried forward as a
    // running sum of accepted deltas: accumulated deltas drift, and a reported
    // cost that disagrees with the returned labels is the worst failure this
    // API could have.
    const double cost = total_cost(assignments);

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

std::vector<size_t> initialize_medoids(const double* data, size_t n,
                                       const KMedoidsOptions& options) {
    switch (options.init) {
        case KMedoidsInit::Build:
            return build_initialize(data, n, options.n_clusters,
                                    options.num_threads, options.chunk_size);
        case KMedoidsInit::FarthestFirst:
            return farthest_first_initialize(data, n, options.n_clusters,
                                             options.num_threads,
                                             options.chunk_size);
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

    if (options.n_clusters == n) {
        // Every item is already a medoid, so the candidate set is empty and
        // there is no swap for the verification pass to examine. Without the
        // short-circuit BUILD would spend O(n^3) reaching the same answer.
        return identity_partition(n);
    }

    const std::vector<size_t> medoids = initialize_medoids(data, n, options);

    return assemble(data, n, medoids, 0, false, options.num_threads,
                    options.chunk_size);
}

}  // namespace OECluster
