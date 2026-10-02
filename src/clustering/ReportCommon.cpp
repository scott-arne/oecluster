/**
 * @file ReportCommon.cpp
 * @brief Profile, partition validation and helpers shared by the cluster reports.
 */

#include "ReportCommon.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "DistanceAccess.h"

namespace OECluster::detail {

namespace {

// Gini coefficient via weighted sum (Lorenz curve area) not pairwise differences
// for O(n log n) not O(n²).
double size_gini(const std::vector<double>& sizes) {
    const size_t n = sizes.size();
    if (n <= 1) {
        return 0.0;
    }
    std::vector<double> sorted = sizes;
    std::sort(sorted.begin(), sorted.end());
    double total = 0.0;
    double weighted = 0.0;
    for (size_t i = 0; i < n; ++i) {
        total += sorted[i];
        weighted += static_cast<double>(i + 1) * sorted[i];
    }
    if (total == 0.0) {
        return 0.0;
    }
    return (2.0 * weighted) / (static_cast<double>(n) * total) -
           (static_cast<double>(n) + 1.0) / static_cast<double>(n);
}

double size_entropy(const std::vector<double>& sizes) {
    if (sizes.size() <= 1) {
        return 0.0;
    }
    double total = 0.0;
    for (const double s : sizes) {
        total += s;
    }
    if (total == 0.0) {
        return 0.0;
    }
    double entropy = 0.0;
    for (const double s : sizes) {
        if (s > 0.0) {
            const double p = s / total;
            entropy -= p * std::log2(p);
        }
    }
    return entropy;
}

}  // namespace

std::vector<double> preset_coverage_thresholds(const ClusterThreshold preset) {
    switch (preset) {
        case ClusterThreshold::Tight:
            return {0.20, 0.30, 0.40};
        case ClusterThreshold::Diversity:
            return {0.40, 0.50, 0.60};
        case ClusterThreshold::Default:
        default:
            return {0.25, 0.35, 0.45};
    }
}

// Fractional-rank percentile matches NumPy/Pandas default interpolation (not nearest-rank).
double percentile(std::vector<double> values, const double q) {
    if (values.empty()) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    std::sort(values.begin(), values.end());
    if (values.size() == 1) {
        return values.front();
    }
    const double rank = q * static_cast<double>(values.size() - 1);
    const size_t lo = static_cast<size_t>(std::floor(rank));
    const size_t hi = static_cast<size_t>(std::ceil(rank));
    const double frac = rank - static_cast<double>(lo);
    return values[lo] + frac * (values[hi] - values[lo]);
}

// One definition, called from the scorecard's mean over all points and from the
// per-cluster records. The two reductions differ -- the scalar averages over
// every clustered point, a record over its own members -- but the term they
// average must not.
//
// Rousseeuw defines s(i) = 0 for a point that is the only member of its
// cluster, and the size test is what implements that. A singleton has no own
// pair to average, so its a term arrives here as 0.0 for want of a distance
// rather than because its neighbours are coincident, and the general formula
// would read that as the perfect score 1.0. Awarding the maximum to a cluster
// of one inflates the mean on exactly the fragmented clusterings this
// scorecard exists to discriminate.
double silhouette_term(const double a_term, const double b_term, const size_t cluster_size) {
    if (cluster_size < 2) {
        return 0.0;
    }
    const double denominator = std::max(a_term, b_term);
    return denominator == 0.0 ? 0.0 : (b_term - a_term) / denominator;
}

ReportProfile report_profile(const ClusteringResult& result,
                             const bool treat_noise_as_singletons) {
    const std::vector<ClusterLabel>& labels = result.Labels();
    const Clusters& members = result.Members();

    ReportProfile profile;
    profile.num_samples = labels.size();
    profile.num_clusters = members.size();

    size_t num_noise = 0;
    for (const ClusterLabel label : labels) {
        if (label < 0) {
            ++num_noise;
        }
    }
    profile.num_noise = num_noise;

    std::vector<double> sizes;
    sizes.reserve(members.size());
    size_t largest = 0;
    size_t num_singletons = 0;
    for (const Cluster& cluster : members) {
        sizes.push_back(static_cast<double>(cluster.size()));
        largest = std::max(largest, cluster.size());
        if (cluster.size() == 1) {
            ++num_singletons;
        }
    }
    profile.num_singletons = num_singletons;

    const double n = static_cast<double>(profile.num_samples);
    profile.noise_fraction = n > 0.0 ? static_cast<double>(num_noise) / n : 0.0;
    profile.largest_cluster_fraction = n > 0.0 ? static_cast<double>(largest) / n : 0.0;

    // Under the flag, each noise point counts as its own singleton cluster, so
    // it is folded into both the numerator and the denominator.
    if (treat_noise_as_singletons) {
        const double denom = static_cast<double>(profile.num_clusters + num_noise);
        profile.singleton_fraction =
            denom > 0.0 ? static_cast<double>(num_singletons + num_noise) / denom : 0.0;
    } else {
        const double denom = static_cast<double>(profile.num_clusters);
        profile.singleton_fraction =
            denom > 0.0 ? static_cast<double>(num_singletons) / denom : 0.0;
    }

    if (sizes.empty()) {
        profile.cluster_size_median = std::numeric_limits<double>::quiet_NaN();
        profile.cluster_size_p90 = std::numeric_limits<double>::quiet_NaN();
    } else {
        profile.cluster_size_median = percentile(sizes, 0.5);
        profile.cluster_size_p90 = percentile(sizes, 0.9);
    }
    profile.size_gini = size_gini(sizes);
    profile.size_entropy = size_entropy(sizes);
    return profile;
}

std::vector<size_t> validate_report_partition(const ClusteringResult& result,
                                              const size_t num_items,
                                              const std::string& caller,
                                              const std::string& noun) {
    const std::vector<ClusterLabel>& labels = result.Labels();
    const Clusters& members = result.Members();

    // The label vector and the storage must describe the same sample set before
    // any per-sample check runs. Left to the pre-pass below, a surplus
    // *clustered* sample is reported as "appears in no cluster", which sends the
    // caller to fix a cluster list when the real error is that they paired a
    // result with the wrong storage. The members-non-empty conjunct preserves
    // the documented acceptance of a long all-noise label vector.
    if (!members.empty() && labels.size() > num_items) {
        throw std::out_of_range(
            caller + ": label count " + std::to_string(labels.size()) +
            " exceeds the " + noun + " sample count " +
            std::to_string(num_items));
    }

    // Every cluster is validated before the first distance read. Otherwise a
    // malformed cluster k would be named only after some cluster j < k had
    // already asked the backend for a pair it cannot answer, and the caller
    // would be told about a storage class instead of the bad cluster member
    // that is the error they have to fix. Checking up front makes that
    // ordering hold whichever cluster is malformed.
    //
    // This runs unconditionally rather than under the members-non-empty guard
    // it used to sit behind: labels = {0} with members = {} is exactly the
    // disagreement being checked for, and the guard would skip it.
    //
    // Checking members against labels alone is not enough. Labels() and
    // Members() must be two spellings of one partition; a one-directional check
    // leaves a clustered sample that no cluster lists -- counted in num_samples
    // and coverage_at, absent from every pair statistic. The owner array is
    // also the sample-to-ordinal lookup the silhouette and the record table
    // need later, so the check pays for itself.
    std::vector<size_t> owner(labels.size(), NO_OWNER);

    // Three staged passes rather than one interleaved loop, so that the layer a
    // refusal comes from does not depend on which cluster holds which error: a
    // malformed cluster is always reported ahead of a partition that
    // double-counts a sample, whatever order the two arrive in.
    //
    // Inside this first pass nothing is canonicalised. The shared validator
    // short-circuits on the first fault it meets, so which of "empty", "outside
    // the storage range" and "not unique" is named depends both on the order of
    // the clusters and on the order of the members within one, and the two
    // answers can differ in exception type as well as in message. That is
    // deliberate: all three tell the caller the same thing -- this cluster list
    // is malformed -- and ordering them here would mean reimplementing checks
    // that belong in DistanceAccess.h. The guarantee is between layers, not
    // inside one.
    for (const Cluster& cluster : members) {
        validate_cluster_members(cluster, num_items);
    }

    // A member can be inside the storage range and past the end of a shorter
    // label vector; without this, reading Labels()[member] below is undefined.
    // It runs to completion before ownership so that a bad index -- another
    // symptom of a result paired with the wrong labels -- is never masked by a
    // duplicate found in an earlier cluster.
    for (const Cluster& cluster : members) {
        for (const size_t member : cluster) {
            if (member >= labels.size()) {
                throw std::out_of_range(
                    caller + ": cluster member " + std::to_string(member) +
                    " is at or beyond the label count " +
                    std::to_string(labels.size()));
            }
        }
    }

    for (size_t k = 0; k < members.size(); ++k) {
        for (const size_t member : members[k]) {
            if (owner[member] != NO_OWNER) {
                throw std::invalid_argument(
                    caller + ": sample " + std::to_string(member) +
                    " appears in clusters " + std::to_string(owner[member]) +
                    " and " + std::to_string(k));
            }
            owner[member] = k;
        }
    }

    for (size_t i = 0; i < labels.size(); ++i) {
        if (labels[i] >= 0) {
            if (owner[i] != static_cast<size_t>(labels[i])) {
                throw std::invalid_argument(
                    owner[i] == NO_OWNER
                        ? caller + ": sample " + std::to_string(i) +
                              " has label " + std::to_string(labels[i]) +
                              " but appears in no cluster"
                        : caller + ": sample " + std::to_string(i) +
                              " has label " + std::to_string(labels[i]) +
                              " but appears in cluster " +
                              std::to_string(owner[i]));
            }
        } else if (owner[i] != NO_OWNER) {
            throw std::invalid_argument(
                caller + ": sample " + std::to_string(i) +
                " is labelled noise but appears in cluster " +
                std::to_string(owner[i]));
        }
    }
    return owner;
}

}  // namespace OECluster::detail
