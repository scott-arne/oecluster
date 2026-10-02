#include "oecluster/clustering/SARCoherence.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <string>
#include <thread>
#include <utility>
#include <vector>

#include "ActivityMetrics.h"
#include "ChunkedComparisons.h"
#include "ContingencyTable.h"
#include "DistanceAccess.h"
#include "DiversityValidation.h"
#include "oecluster/ThreadPool.h"

namespace OECluster {

namespace {

/// Shared by the two continuous-activity entry points. owner names whatever
/// fixes the expected length, so the message reads the way the caller thinks
/// about the mismatch.
void validate_double_activity(const std::vector<double>& activity,
                              std::size_t expected, const std::string& caller,
                              const std::string& owner) {
    if (activity.empty()) {
        throw std::invalid_argument(caller +
                                    " requires a non-empty activity vector");
    }
    if (activity.size() != expected) {
        throw std::invalid_argument(
            caller + ": activity has " + std::to_string(activity.size()) +
            " entries but " + owner + " has " + std::to_string(expected) +
            " samples");
    }
}

SARCoherence sar_coherence_impl(const std::vector<ClusterLabel>& labels,
                                const std::vector<double>& activity,
                                const SARCoherenceOptions& options) {
    validate_double_activity(activity, labels.size(), "sar_coherence",
                             "the clustering");

    SARCoherence coherence;
    coherence.num_samples = activity.size();

    // Every exclusion reason is ORed into one mask before anything is
    // gathered. Gathering after only one of them leaves the value vector and
    // the group ids at different lengths.
    std::vector<bool> drop(activity.size(), false);
    detail::mark_excluded(labels, options.noise_handling, drop);
    detail::mark_missing_activity(activity, "sar_coherence", drop);

    const detail::ScoredActivity scored = detail::gather_scored(activity, drop);
    std::uint32_t num_ids = 0;
    const std::vector<std::uint32_t> ids =
        detail::intern_side(labels, options.noise_handling, drop, num_ids);

    coherence.num_scored = scored.values.size();
    coherence.num_clusters = num_ids;

    // The decomposition runs before the per-cluster table so that an
    // unformable mean is reported as a throw rather than written into a row.
    const detail::SumsOfSquares ss =
        detail::sums_of_squares(ids, scored.values, num_ids, "sar_coherence");
    coherence.eta_squared = detail::eta_squared(ss);
    coherence.omega_squared = detail::omega_squared(ss);

    std::vector<ClusterLabel> id_label(num_ids, 0);
    std::vector<bool> id_seen(num_ids, false);
    std::vector<std::vector<double>> group_values(num_ids);
    for (std::size_t k = 0; k < scored.values.size(); ++k) {
        const std::uint32_t id = ids[k];
        if (!id_seen[id]) {
            id_seen[id] = true;
            const ClusterLabel label = labels[scored.indices[k]];
            id_label[id] = (options.noise_handling == NoiseHandling::Grouped &&
                            detail::is_noise(label))
                               ? NOISE_LABEL
                               : label;
        }
        group_values[id].push_back(scored.values[k]);
    }

    coherence.clusters.reserve(num_ids);
    for (std::uint32_t id = 0; id < num_ids; ++id) {
        ClusterActivity row;
        row.label = id_label[id];
        row.num_scored = group_values[id].size();
        double sum = 0.0;
        for (const double value : group_values[id]) {
            sum += value;
        }
        const double n = static_cast<double>(row.num_scored);
        const double mean_0 = sum / n;

        // The same correction pass sums_of_squares applies to its group means.
        // group_values[id] holds this group's values in the order sums_of_squares
        // accumulates them, so the two land on the same value; without it the
        // published mean can sit an ulp off the mean the decomposition actually used.
        double correction = 0.0;
        for (const double value : group_values[id]) {
            correction += value - mean_0;
        }
        row.mean_activity = mean_0 + correction / n;
        row.stddev_activity = detail::population_stddev(group_values[id]);
        coherence.clusters.push_back(std::move(row));
    }

    return coherence;
}

void validate_landscape_options(const ActivityLandscapeOptions& options) {
    const std::pair<const char*, double> checks[] = {
        {"distance_threshold", options.distance_threshold},
        {"activity_threshold", options.activity_threshold},
        {"rmodi_delta", options.rmodi_delta},
    };
    for (const auto& check : checks) {
        if (!std::isfinite(check.second)) {
            throw std::invalid_argument(std::string("activity_landscape: ") +
                                        check.first + " must be finite");
        }
        if (check.second < 0.0) {
            throw std::invalid_argument(std::string("activity_landscape: ") +
                                        check.first + " must be non-negative");
        }
    }
}

void validate_class_annotation(const std::vector<std::string>& activity_classes,
                               std::size_t expected, const std::string& owner) {
    if (activity_classes.empty()) {
        throw std::invalid_argument(
            "modelability requires a non-empty class annotation");
    }
    if (activity_classes.size() != expected) {
        throw std::invalid_argument(
            "modelability: activity_classes has " +
            std::to_string(activity_classes.size()) + " entries but " + owner +
            " has " + std::to_string(expected) + " samples");
    }
}

// The MODI sweep reads each pair from both rows, so without a canonical order
// the message names "0 and 2" or "2 and 0" depending on which worker threw
// first. Sorting the pair makes the diagnostic reproducible; activity_landscape
// needs no equivalent because its sweep visits each pair once, from the lower
// row.
std::invalid_argument bad_distance(std::size_t left, std::size_t right) {
    return std::invalid_argument(
        "modelability: the distance between samples " +
        std::to_string(std::min(left, right)) + " and " +
        std::to_string(std::max(left, right)) +
        " must be finite and non-negative");
}

// Each activity_landscape chunk allocates and merges two length-n buffers, so
// a small chunk makes that O(n) bookkeeping dominate the O(n) row it was meant
// to serve. Both of its sweeps apply this floor.
constexpr std::size_t LANDSCAPE_MIN_ROWS = 64;

// The storage overloads' traversal, as it always was: rows split across a
// pool, each read straight from Data().
class MatrixSweep {
public:
    MatrixSweep(const StorageBackend& storage, std::size_t num_threads,
                std::size_t min_rows)
        : data_(storage.Data()),
          num_samples_(storage.NumSamples()),
          num_threads_(num_threads),
          min_rows_(min_rows) {}

    // body(read, begin, end) over rows [0, n), where read(a, b) is the
    // distance between samples a and b. Callers guarantee n >= 2.
    template <typename Body>
    void Run(std::size_t n, Body&& body) const {
        // Capped at the row count before the pool is built. num_threads is a
        // size_t on a public options struct, so this cap is the only thing
        // standing between a caller and ThreadPool trying to spawn 2^61 OS
        // threads. It is not what makes 8 * threads below safe to form: the
        // next statement bounds threads by n whatever the pool reports, so the
        // product is at most 8n. A num_threads of 0 bypasses the cap by
        // design, because it means "use the hardware concurrency" -- so a
        // small n still spawns that many workers for what is a single chunk,
        // matching how HDBSCAN and Agglomerative already size their pools.
        // Because n >= 2, no value other than that deliberate 0 can reach zero
        // here.
        ThreadPool pool(std::min<std::size_t>(num_threads_, n));
        const std::size_t threads =
            std::max<std::size_t>(1, std::min<std::size_t>(pool.NumThreads(), n));
        const std::size_t chunk_size =
            std::max<std::size_t>(min_rows_, n / (8 * threads));
        const auto read = [this](std::size_t a, std::size_t b) {
            return detail::dense_distance(data_, num_samples_, a, b);
        };
        pool.ParallelFor(0, n, chunk_size, [&](std::size_t begin, std::size_t end) {
            body(read, begin, end);
        });
    }

    // body(read, 0, n) on the calling thread.
    template <typename Body>
    void Serial(std::size_t n, Body&& body) const {
        const auto read = [this](std::size_t a, std::size_t b) {
            return detail::dense_distance(data_, num_samples_, a, b);
        };
        body(read, 0, n);
    }

private:
    const double* data_;
    std::size_t num_samples_;
    std::size_t num_threads_;
    std::size_t min_rows_;
};

// The comparison overloads' traversal: units of whole rows, about chunk_size
// distances each, each unit on its own clone. Every pair is read as
// Compare(min, max), so modelability's two visits of a pair see one value.
class ComparisonSweep {
public:
    ComparisonSweep(const PairwiseComparison& comparison, std::size_t num_threads,
                    std::size_t chunk_size, std::size_t min_rows)
        : comparison_(comparison),
          num_threads_(num_threads),
          chunk_size_(chunk_size),
          min_rows_(min_rows) {}

    // As MatrixSweep::Run. A row holds at most n - 1 distances.
    template <typename Body>
    void Run(std::size_t n, Body&& body) const {
        const std::size_t rows_per_unit =
            std::max<std::size_t>(min_rows_, chunk_size_ / (n - 1));
        // Zero is resolved here rather than in the pool so the item cap
        // applies to it too: an unresolved zero would let the hardware
        // concurrency, not n, bound the workers and clones.
        const std::size_t workers =
            num_threads_ > 0
                ? num_threads_
                : std::max<std::size_t>(1, std::thread::hardware_concurrency());
        detail::ChunkedComparisons units(comparison_, n, workers, rows_per_unit);
        units.Run(n, [&](PairwiseComparison& clone, std::size_t begin, std::size_t end) {
            const auto read = [&clone](std::size_t a, std::size_t b) {
                return clone.Compare(std::min(a, b), std::max(a, b));
            };
            body(read, begin, end);
        });
    }

    // body(read, 0, n) on the calling thread, over one clone: the caller's
    // comparison is a prototype and is never compared on directly.
    template <typename Body>
    void Serial(std::size_t n, Body&& body) const {
        const std::unique_ptr<PairwiseComparison> clone = comparison_.Clone();
        const auto read = [&clone](std::size_t a, std::size_t b) {
            return clone->Compare(std::min(a, b), std::max(a, b));
        };
        body(read, 0, n);
    }

private:
    const PairwiseComparison& comparison_;
    std::size_t num_threads_;
    std::size_t chunk_size_;
    std::size_t min_rows_;
};

template <typename Sweep>
ActivityLandscape landscape_engine(const Sweep& sweep,
                                   const std::vector<double>& activity,
                                   const ActivityLandscapeOptions& options) {
    ActivityLandscape landscape;
    landscape.num_samples = activity.size();

    std::vector<bool> drop(activity.size(), false);
    detail::mark_missing_activity(activity, "activity_landscape", drop);
    const detail::ScoredActivity scored = detail::gather_scored(activity, drop);

    const std::size_t n = scored.values.size();
    landscape.num_scored = n;
    landscape.num_pairs_scored = n < 2 ? 0 : n * (n - 1) / 2;

    // The band is an input to the sweep, so an unusable spread has to be
    // caught before any distance is read.
    landscape.activity_stddev = detail::population_stddev(scored.values);
    if (std::isinf(landscape.activity_stddev)) {
        throw std::invalid_argument(
            "activity_landscape: activity_stddev overflows to infinity; the "
            "supported range is values whose squared deviations sum finitely "
            "in double precision");
    }
    if (n < 2) {
        return landscape;
    }

    const double band = options.rmodi_delta * landscape.activity_stddev;
    constexpr double INFTY = std::numeric_limits<double>::infinity();

    std::vector<std::size_t> row_cliffs(n, 0);
    std::vector<std::size_t> row_zero_pairs(n, 0);
    std::vector<std::size_t> row_sali_count(n, 0);
    std::vector<double> row_sali_sum(n, 0.0);
    std::vector<double> row_max(n, -INFTY);
    std::vector<double> same_min(n, INFTY);
    std::vector<double> diff_min(n, INFTY);
    std::mutex merge_mutex;

    sweep.Run(n, [&](const auto& read, std::size_t begin, std::size_t end) {
        std::vector<double> local_same(n, INFTY);
        std::vector<double> local_diff(n, INFTY);
        for (std::size_t p = begin; p < end; ++p) {
            for (std::size_t q = p + 1; q < n; ++q) {
                const double distance = read(scored.indices[p], scored.indices[q]);
                if (!std::isfinite(distance) || distance < 0.0) {
                    throw std::invalid_argument(
                        "activity_landscape: the distance between samples " +
                        std::to_string(scored.indices[p]) + " and " +
                        std::to_string(scored.indices[q]) +
                        " must be finite and non-negative");
                }
                // Not checked for infinity, and that is a claim rather than an
                // oversight. |v_p - v_q| <= 2 * max_i |v_i - mean|, so a delta
                // that overflows forces a scale above DBL_MAX/2, whose square
                // is unrepresentable -- and population_stddev reports infinity
                // for exactly that, above, before any distance is read. A check
                // here would be unreachable. RejectsAnOverflowingSpread pins
                // the ordering with the case that would otherwise invert RMODI.
                // The unreachability is borrowed from population_stddev's
                // declared domain rather than from arithmetic alone -- see
                // ActivityMetrics.h, where refusing {-1e300, 1e300} is a
                // deliberate agreement with sums_of_squares even though the
                // scaling could represent that spread -- so widening that
                // domain makes an overflowing delta reachable here, and this
                // sweep would then need its own check.
                const double delta =
                    std::fabs(scored.values[p] - scored.values[q]);

                if (distance <= options.distance_threshold &&
                    delta >= options.activity_threshold) {
                    ++row_cliffs[p];
                }

                if (distance == 0.0) {
                    ++row_zero_pairs[p];
                } else {
                    const double sali = delta / distance;
                    if (sali > row_max[p]) {
                        row_max[p] = sali;
                    }
                    row_sali_sum[p] += sali;
                    ++row_sali_count[p];
                }

                if (delta <= band) {
                    if (distance < local_same[p]) {
                        local_same[p] = distance;
                    }
                    if (distance < local_same[q]) {
                        local_same[q] = distance;
                    }
                } else {
                    if (distance < local_diff[p]) {
                        local_diff[p] = distance;
                    }
                    if (distance < local_diff[q]) {
                        local_diff[q] = distance;
                    }
                }
            }
        }
        const std::lock_guard<std::mutex> lock(merge_mutex);
        for (std::size_t i = 0; i < n; ++i) {
            if (local_same[i] < same_min[i]) {
                same_min[i] = local_same[i];
            }
            if (local_diff[i] < diff_min[i]) {
                diff_min[i] = local_diff[i];
            }
        }
    });

    // Combined in ascending row order, after the join. This is what makes the
    // floating-point sums independent of the thread count.
    std::size_t num_cliffs = 0;
    std::size_t num_zero_pairs = 0;
    std::size_t sali_count = 0;
    double sali_sum = 0.0;
    double max_sali = -INFTY;
    for (std::size_t p = 0; p < n; ++p) {
        num_cliffs += row_cliffs[p];
        num_zero_pairs += row_zero_pairs[p];
        sali_count += row_sali_count[p];
        sali_sum += row_sali_sum[p];
        if (row_max[p] > max_sali) {
            max_sali = row_max[p];
        }
    }

    if (std::isinf(sali_sum)) {
        throw std::invalid_argument(
            "activity_landscape: mean_sali overflows to infinity; a pair's "
            "activity difference divided by its distance is not representable");
    }
    // max_sali starts at negative infinity, so only a positive infinity here
    // is an overflow rather than an empty set.
    if (std::isinf(max_sali) && max_sali > 0.0) {
        throw std::invalid_argument(
            "activity_landscape: max_sali overflows to infinity; a pair's "
            "activity difference divided by its distance is not representable");
    }

    landscape.num_cliffs = num_cliffs;
    landscape.num_zero_distance_pairs = num_zero_pairs;
    landscape.cliff_density = static_cast<double>(num_cliffs) /
                              static_cast<double>(landscape.num_pairs_scored);
    if (sali_count > 0) {
        landscape.max_sali = max_sali;
        landscape.mean_sali = sali_sum / static_cast<double>(sali_count);
    }

    std::size_t concordant = 0;
    for (std::size_t p = 0; p < n; ++p) {
        if (same_min[p] < diff_min[p]) {
            ++concordant;
        }
    }
    landscape.rmodi = static_cast<double>(concordant) / static_cast<double>(n);

    return landscape;
}

template <typename Sweep>
Modelability modelability_engine(const Sweep& sweep,
                                 const std::vector<std::string>& activity_classes) {
    Modelability model;
    model.num_samples = activity_classes.size();

    // An empty class string is this input's "missing", and Excluded is the
    // only reading of it that makes sense: a sample with no class cannot be a
    // class of its own, and cannot be pooled with the others either.
    std::vector<bool> drop(activity_classes.size(), false);
    detail::mark_excluded(activity_classes, NoiseHandling::Excluded, drop);

    const std::vector<std::size_t> indices = detail::gather_indices(drop);
    std::uint32_t num_ids = 0;
    const std::vector<std::uint32_t> class_ids = detail::intern_side(
        activity_classes, NoiseHandling::Excluded, drop, num_ids);

    const std::size_t n = indices.size();
    model.num_scored = n;
    model.num_classes = num_ids;

    std::vector<std::string> id_label(num_ids);
    std::vector<std::size_t> id_count(num_ids, 0);
    std::vector<bool> id_seen(num_ids, false);
    for (std::size_t p = 0; p < n; ++p) {
        const std::uint32_t id = class_ids[p];
        if (!id_seen[id]) {
            id_seen[id] = true;
            id_label[id] = activity_classes[indices[p]];
        }
        ++id_count[id];
    }

    model.classes.reserve(num_ids);
    if (num_ids < 2) {
        // Nothing to compute: with one class no neighbour can differ, so the
        // sweep would only confirm that every molecule matches itself. The
        // distances still have to be validated, though, or a caller with one
        // class and a corrupt matrix gets NaNs where every other input shape
        // gets a refusal. This scan is the only O(n^2) work on a path that
        // would otherwise read no distance at all, and it is serial, so its
        // diagnostic needs no canonicalisation.
        sweep.Serial(n, [&](const auto& read, std::size_t begin, std::size_t end) {
            for (std::size_t p = begin; p < end; ++p) {
                for (std::size_t q = p + 1; q < n; ++q) {
                    const double distance = read(indices[p], indices[q]);
                    if (!std::isfinite(distance) || distance < 0.0) {
                        throw bad_distance(indices[p], indices[q]);
                    }
                }
            }
        });
        for (std::uint32_t id = 0; id < num_ids; ++id) {
            ClassConcordance row;
            row.label = id_label[id];
            row.num_members = id_count[id];
            model.classes.push_back(std::move(row));
        }
        return model;
    }

    const std::size_t no_neighbour = n;
    // char, not bool: distinct elements of a vector<bool> share a word, so
    // writing different indices from different threads is a data race.
    std::vector<char> concordant(n, 0);

    sweep.Run(n, [&](const auto& read, std::size_t begin, std::size_t end) {
        for (std::size_t p = begin; p < end; ++p) {
            double best_distance = std::numeric_limits<double>::infinity();
            std::size_t best_q = no_neighbour;
            for (std::size_t q = 0; q < n; ++q) {
                if (q == p) {
                    continue;
                }
                const double distance = read(indices[p], indices[q]);
                if (!std::isfinite(distance) || distance < 0.0) {
                    throw bad_distance(indices[p], indices[q]);
                }
                // Strict, over an ascending scan: ties resolve to the lowest
                // scored index, which is what makes the result independent of
                // the thread count.
                if (distance < best_distance) {
                    best_distance = distance;
                    best_q = q;
                }
            }
            concordant[p] = (best_q != no_neighbour &&
                             class_ids[best_q] == class_ids[p])
                                ? 1
                                : 0;
        }
    });

    std::vector<std::size_t> id_concordant(num_ids, 0);
    for (std::size_t p = 0; p < n; ++p) {
        if (concordant[p] != 0) {
            ++id_concordant[class_ids[p]];
        }
    }

    double fraction_sum = 0.0;
    for (std::uint32_t id = 0; id < num_ids; ++id) {
        ClassConcordance row;
        row.label = id_label[id];
        row.num_members = id_count[id];
        row.fraction_same_class = static_cast<double>(id_concordant[id]) /
                                  static_cast<double>(id_count[id]);
        fraction_sum += row.fraction_same_class;
        model.classes.push_back(std::move(row));
    }
    model.modi = fraction_sum / static_cast<double>(num_ids);

    return model;
}

}  // namespace

SARCoherence sar_coherence(const ClusteringResult& result,
                           const std::vector<double>& activity,
                           const SARCoherenceOptions& options) {
    return sar_coherence_impl(result.Labels(), activity, options);
}

SARCoherence sar_coherence(const std::vector<ClusterLabel>& labels,
                           const std::vector<double>& activity,
                           const SARCoherenceOptions& options) {
    return sar_coherence_impl(labels, activity, options);
}

ActivityLandscape activity_landscape(const StorageBackend& storage,
                                     const std::vector<double>& activity,
                                     const ActivityLandscapeOptions& options) {
    detail::validate_complete_distance_storage(storage, "activity_landscape");
    validate_double_activity(activity, storage.NumSamples(),
                             "activity_landscape", "the storage");
    validate_landscape_options(options);
    const MatrixSweep sweep(storage, options.num_threads, LANDSCAPE_MIN_ROWS);
    return landscape_engine(sweep, activity, options);
}

ActivityLandscape activity_landscape(PairwiseComparison& comparison,
                                     const std::vector<double>& activity,
                                     const ActivityLandscapeOptions& options) {
    detail::validate_chunk_size(options.chunk_size, "activity_landscape");
    detail::validate_comparison_facts(comparison, "activity_landscape");
    validate_double_activity(activity, comparison.Size(),
                             "activity_landscape", "the comparison");
    validate_landscape_options(options);
    const ComparisonSweep sweep(comparison, options.num_threads,
                                options.chunk_size, LANDSCAPE_MIN_ROWS);
    return landscape_engine(sweep, activity, options);
}

Modelability modelability(const StorageBackend& storage,
                          const std::vector<std::string>& activity_classes,
                          const ModelabilityOptions& options) {
    detail::validate_complete_distance_storage(storage, "modelability");
    validate_class_annotation(activity_classes, storage.NumSamples(),
                              "the storage");
    const MatrixSweep sweep(storage, options.num_threads, 1);
    return modelability_engine(sweep, activity_classes);
}

Modelability modelability(PairwiseComparison& comparison,
                          const std::vector<std::string>& activity_classes,
                          const ModelabilityOptions& options) {
    detail::validate_chunk_size(options.chunk_size, "modelability");
    detail::validate_comparison_facts(comparison, "modelability");
    validate_class_annotation(activity_classes, comparison.Size(),
                              "the comparison");
    const ComparisonSweep sweep(comparison, options.num_threads,
                                options.chunk_size, 1);
    return modelability_engine(sweep, activity_classes);
}

}  // namespace OECluster
