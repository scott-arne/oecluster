#include "oecluster/clustering/SARCoherence.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <mutex>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "ActivityMetrics.h"
#include "ContingencyTable.h"
#include "DistanceAccess.h"
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

    const std::size_t num_samples = storage.NumSamples();
    const double* data = storage.Data();
    const double band = options.rmodi_delta * landscape.activity_stddev;
    constexpr double INFINITE = std::numeric_limits<double>::infinity();

    std::vector<std::size_t> row_cliffs(n, 0);
    std::vector<std::size_t> row_zero_pairs(n, 0);
    std::vector<std::size_t> row_sali_count(n, 0);
    std::vector<double> row_sali_sum(n, 0.0);
    std::vector<double> row_max(n, -INFINITE);
    std::vector<double> same_min(n, INFINITE);
    std::vector<double> diff_min(n, INFINITE);
    std::mutex merge_mutex;

    // Capped at the row count before the pool is built. num_threads is a
    // size_t on a public options struct, so the only thing standing between a
    // caller and ThreadPool trying to spawn 2^61 OS threads is this line; and
    // a worker with no row to take is pure overhead even at sane values. The
    // cap is also what makes 8 * threads below safe to form: a matrix with n
    // rows needs n^2/2 doubles, so n is nowhere near the value at which the
    // multiplication could wrap. A num_threads of 0 means "use the hardware
    // concurrency" and passes through the cap unchanged; the early return
    // above guarantees n >= 2, so no other value can reach zero here.
    ThreadPool pool(std::min<std::size_t>(options.num_threads, n));
    const std::size_t threads =
        std::max<std::size_t>(1, std::min<std::size_t>(pool.NumThreads(), n));
    // Each chunk allocates and merges two length-n buffers, so a small chunk
    // makes that O(n) bookkeeping dominate the O(n) row it was meant to serve.
    const std::size_t chunk_size = std::max<std::size_t>(64, n / (8 * threads));

    pool.ParallelFor(0, n, chunk_size, [&](std::size_t begin, std::size_t end) {
        std::vector<double> local_same(n, INFINITE);
        std::vector<double> local_diff(n, INFINITE);
        for (std::size_t p = begin; p < end; ++p) {
            for (std::size_t q = p + 1; q < n; ++q) {
                const double distance = detail::dense_distance(
                    data, num_samples, scored.indices[p], scored.indices[q]);
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
    double max_sali = -INFINITE;
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

}  // namespace OECluster
