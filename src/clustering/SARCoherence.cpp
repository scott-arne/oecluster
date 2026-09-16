#include "oecluster/clustering/SARCoherence.h"

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "ActivityMetrics.h"
#include "ContingencyTable.h"

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
        row.mean_activity = sum / static_cast<double>(row.num_scored);
        row.stddev_activity = detail::population_stddev(group_values[id]);
        coherence.clusters.push_back(std::move(row));
    }

    return coherence;
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

}  // namespace OECluster
