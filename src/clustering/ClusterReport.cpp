/**
 * @file ClusterReport.cpp
 * @brief Clustering-quality report implementation.
 */

#include "oecluster/clustering/ClusterReport.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>

#include "ClusterMetrics.h"
#include "DistanceAccess.h"
#include "InternalIndices.h"

namespace OECluster {

ClusterReportOptions::ClusterReportOptions(ClusterThreshold preset) {
    switch (preset) {
        case ClusterThreshold::Tight:
            coverage_thresholds = {0.20, 0.30, 0.40};
            boundary_threshold = 0.25;
            break;
        case ClusterThreshold::Diversity:
            coverage_thresholds = {0.40, 0.50, 0.60};
            boundary_threshold = 0.40;
            break;
        case ClusterThreshold::Default:
        default:
            coverage_thresholds = {0.25, 0.35, 0.45};
            boundary_threshold = 0.30;
            break;
    }
}

namespace {

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

double nan_value() {
    return std::numeric_limits<double>::quiet_NaN();
}

// The number of unordered pairs among n items, refusing rather than wrapping.
// The three pair-scaled arrays below are reserved from this, and a wrapped
// count under-reserves in silence -- which would reinstate exactly the
// push_back growth the reservation exists to remove. Exactly one of n and
// n - 1 is even, so halving that one first is exact and the guard fires only
// when the pair count itself does not fit, not when an intermediate product
// would have overflowed on its way to a count that does.
//
// std::length_error rather than a bespoke type, following
// detail::add_couples in InternalIndices.h: a caller at this scale is already
// in memory trouble, and the SWIG layer maps length_error to Python
// MemoryError.
size_t pair_count(const size_t n) {
    if (n < 2) {
        return 0;
    }
    const bool even = n % 2 == 0;
    const size_t half = even ? n / 2 : (n - 1) / 2;
    const size_t other = even ? n - 1 : n;
    if (half > std::numeric_limits<size_t>::max() / other) {
        throw std::length_error(
            "cluster_report: the pair count of " + std::to_string(n) +
            " items exceeds the range of size_t");
    }
    return half * other;
}

// Accumulates pair counts under the same refusal, for the same reason.
void add_pair_count(size_t& total, const size_t increment) {
    if (total > std::numeric_limits<size_t>::max() - increment) {
        throw std::length_error(
            "cluster_report: the total pair count exceeds the range of size_t (" +
            std::to_string(total) + " + " + std::to_string(increment) + ")");
    }
    total += increment;
}

// One definition, called from the scorecard's mean over all points and from the
// per-cluster records. The two reductions differ -- the scalar averages over
// every clustered point, a record over its own members -- but the term they
// average must not.
double silhouette_term(const double a_term, const double b_term) {
    const double denominator = std::max(a_term, b_term);
    return denominator == 0.0 ? 0.0 : (b_term - a_term) / denominator;
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

ClusterReport cluster_report(
    const ClusteringResult& result,
    const StorageBackend& storage,
    const ClusterReportOptions& options) {
    detail::validate_complete_distance_storage(storage, "cluster_report");

    ClusterReport report;
    report.coverage_thresholds = options.coverage_thresholds;

    // Set before any branch, so an empty records table on a K == 0 result still
    // says whether anybody asked for one.
    report.requested.pair_rank_indices = options.compute_pair_rank_indices;
    report.requested.per_cluster_records = options.compute_per_cluster_records;

    const std::vector<ClusterLabel>& labels = result.Labels();
    const Clusters& members = result.Members();

    report.num_samples = labels.size();
    report.num_clusters = members.size();

    size_t num_noise = 0;
    for (const ClusterLabel label : labels) {
        if (label < 0) {
            ++num_noise;
        }
    }
    report.num_noise = num_noise;

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
    report.num_singletons = num_singletons;

    const double n = static_cast<double>(report.num_samples);
    report.noise_fraction = n > 0.0 ? static_cast<double>(num_noise) / n : 0.0;
    report.largest_cluster_fraction = n > 0.0 ? static_cast<double>(largest) / n : 0.0;

    // Under the flag, each noise point counts as its own singleton cluster, so
    // it is folded into both the numerator and the denominator.
    if (options.treat_noise_as_singletons) {
        const double denom = static_cast<double>(report.num_clusters + num_noise);
        report.singleton_fraction =
            denom > 0.0 ? static_cast<double>(num_singletons + num_noise) / denom : 0.0;
    } else {
        const double denom = static_cast<double>(report.num_clusters);
        report.singleton_fraction =
            denom > 0.0 ? static_cast<double>(num_singletons) / denom : 0.0;
    }

    if (sizes.empty()) {
        report.cluster_size_median = std::numeric_limits<double>::quiet_NaN();
        report.cluster_size_p90 = std::numeric_limits<double>::quiet_NaN();
    } else {
        report.cluster_size_median = percentile(sizes, 0.5);
        report.cluster_size_p90 = percentile(sizes, 0.9);
    }
    report.size_gini = size_gini(sizes);
    report.size_entropy = size_entropy(sizes);

    // The label vector and the storage must describe the same sample set before
    // any per-sample check runs. Left to the pre-pass below, a surplus
    // *clustered* sample is reported as "appears in no cluster", which sends the
    // caller to fix a cluster list when the real error is that they paired a
    // result with the wrong storage. The members-non-empty conjunct preserves
    // the documented acceptance of a long all-noise label vector.
    if (!members.empty() && labels.size() > storage.NumSamples()) {
        throw std::out_of_range(
            "cluster_report: label count " + std::to_string(labels.size()) +
            " exceeds the storage sample count " +
            std::to_string(storage.NumSamples()));
    }

    // Hoisted above the first storage read. cluster_representative validates
    // the cluster it is handed, but only when the intra pass reaches that
    // cluster -- so a malformed cluster k would be named only after some
    // cluster j < k had already asked the backend for a pair it cannot answer,
    // and the caller was told about a storage class instead of the bad cluster
    // member that is the error they have to fix. Checking every cluster before
    // any of them is processed makes that ordering hold whichever one is
    // malformed.
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
    constexpr size_t NO_OWNER = std::numeric_limits<size_t>::max();
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
        detail::validate_cluster_members(cluster, storage.NumSamples());
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
                    "cluster_report: cluster member " + std::to_string(member) +
                    " is at or beyond the label count " +
                    std::to_string(labels.size()));
            }
        }
    }

    for (size_t k = 0; k < members.size(); ++k) {
        for (const size_t member : members[k]) {
            if (owner[member] != NO_OWNER) {
                throw std::invalid_argument(
                    "cluster_report: sample " + std::to_string(member) +
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
                        ? "cluster_report: sample " + std::to_string(i) +
                              " has label " + std::to_string(labels[i]) +
                              " but appears in no cluster"
                        : "cluster_report: sample " + std::to_string(i) +
                              " has label " + std::to_string(labels[i]) +
                              " but appears in cluster " +
                              std::to_string(owner[i]));
            }
        } else if (owner[i] != NO_OWNER) {
            throw std::invalid_argument(
                "cluster_report: sample " + std::to_string(i) +
                " is labelled noise but appears in cluster " +
                std::to_string(owner[i]));
        }
    }

    if (!members.empty()) {
        const size_t cluster_count = members.size();

        // Refused before a single distance is read, so a method the caller
        // cannot configure outranks the finiteness refusal below without the
        // selection itself having to run early -- running it early would feed
        // an unchecked NaN to the selector's sorts. INVARIANT 1. Naming the
        // unsupported method here duplicates knowledge validate_options in
        // Representative.cpp also holds, and the build enables neither -Wall
        // nor -Wswitch, so a fifth RepresentativeMethod would be flagged at
        // neither site.
        if (options.representative_method == RepresentativeMethod::HighestNeighborhood) {
            throw std::invalid_argument(
                "cluster_report: representative_method HighestNeighborhood is "
                "unsupported because ClusterReportOptions carries no neighbor "
                "threshold to configure it");
        }

        // Ranked behind the method refusal above and ahead of every threshold
        // read. A NaN threshold does not propagate: `distance <= threshold` and
        // `nearest <= threshold` are both false for it, so the caller is handed
        // zero boundary violations and zero coverage -- a plausible number, not
        // an error, and indistinguishable from a well-separated clustering.
        // Passing an explicit NaN must not read as passing nothing. INVARIANT 2.
        //
        // std::isnan, not !std::isfinite. An infinite boundary_threshold asks
        // for every cross pair to count and an infinite coverage threshold for
        // every sample to be covered; both are answered exactly today, so
        // refusing them would be over-refusal. INVARIANT 3.
        //
        // These sit inside the members-non-empty guard, which means a partition
        // with no clusters is still accepted with a NaN threshold. That report
        // reads no distance and evaluates no threshold -- coverage_at stays
        // empty and boundary_violations is zero for want of a pair to count,
        // not for want of a comparison that held -- so there is no wrong number
        // for the NaN to hide behind, and refusing it would be over-refusal of
        // the same kind. INVARIANT 3 again.
        if (std::isnan(options.boundary_threshold)) {
            throw std::invalid_argument(
                "cluster_report: boundary_threshold must not be NaN");
        }
        for (size_t t = 0; t < options.coverage_thresholds.size(); ++t) {
            if (std::isnan(options.coverage_thresholds[t])) {
                throw std::invalid_argument(
                    "cluster_report: coverage threshold " + std::to_string(t) +
                    " must not be NaN");
            }
        }

        // Every pair-scaled reservation below is derived here, from the same
        // member lists the fill loops walk, so a reservation cannot disagree
        // with what is pushed into it. Computed once: the three counts differ
        // in reduction -- a sum, a maximum, and a complement -- but not in
        // input.
        const size_t clustered_count = std::accumulate(
            members.begin(),
            members.end(),
            size_t{0},
            [](const size_t total, const Cluster& cluster) {
                return total + cluster.size();
            });
        size_t intra_pair_count = 0;
        size_t largest_cluster_pair_count = 0;
        for (const Cluster& cluster : members) {
            const size_t pairs = pair_count(cluster.size());
            add_pair_count(intra_pair_count, pairs);
            largest_cluster_pair_count = std::max(largest_cluster_pair_count, pairs);
        }

        // ---- Intra pass: once per cluster. ----
        std::vector<double> intra_pairs;
        std::vector<double> radii;
        std::vector<double> diameters;
        std::vector<double> medoid_member_means;
        // Pair-scaled, not cluster-scaled like the three below it: this one
        // takes every within-cluster distance, sum_k C(n_k, 2) of them, and
        // takes them unconditionally because median_intra_distance reads it.
        intra_pairs.reserve(intra_pair_count);
        radii.reserve(cluster_count);
        diameters.reserve(cluster_count);
        medoid_member_means.reserve(cluster_count);

        std::vector<size_t> representatives;
        std::vector<size_t> true_medoids;
        representatives.reserve(cluster_count);
        true_medoids.reserve(cluster_count);

        std::vector<double> cluster_intra_sums(cluster_count, 0.0);
        std::vector<size_t> cluster_intra_counts(cluster_count, 0);
        std::vector<double> medoid_scatter(cluster_count, 0.0);
        std::vector<double> medoid_square_sums(cluster_count, 0.0);

        std::vector<double> cluster_medians(cluster_count, nan_value());
        // Cleared per cluster rather than grown across them, so its high-water
        // mark is the largest single cluster's pair count and not the sum.
        // clear() keeps capacity, so the one reservation here serves every
        // iteration of the loop below.
        std::vector<double> cluster_distances;
        if (options.compute_per_cluster_records) {
            cluster_distances.reserve(largest_cluster_pair_count);
        }

        std::vector<double> own_mean(labels.size(), 0.0);
        std::vector<double> point_total(labels.size(), 0.0);
        detail::DistanceMoments within_moments;

        double max_diameter = 0.0;
        for (size_t k = 0; k < cluster_count; ++k) {
            const Cluster& cluster = members[k];
            cluster_distances.clear();

            double diameter = 0.0;
            for (size_t i = 0; i < cluster.size(); ++i) {
                for (size_t j = i + 1; j < cluster.size(); ++j) {
                    const double distance =
                        detail::checked_distance(storage, cluster[i], cluster[j]);
                    intra_pairs.push_back(distance);
                    if (options.compute_per_cluster_records) {
                        cluster_distances.push_back(distance);
                    }
                    within_moments.Add(distance);
                    diameter = std::max(diameter, distance);
                    cluster_intra_sums[k] += distance;
                    ++cluster_intra_counts[k];
                    // Accumulated onto both endpoints in ascending partner
                    // order, which is the order mean_to_cluster used to sum in,
                    // so the silhouette's a term is bit-identical to 5.0.0's.
                    own_mean[cluster[i]] += distance;
                    own_mean[cluster[j]] += distance;
                }
            }
            if (options.compute_per_cluster_records && !cluster_distances.empty()) {
                cluster_medians[k] = detail::median_distance(cluster_distances);
            }
            diameters.push_back(diameter);
            max_diameter = std::max(max_diameter, diameter);

            // Selected only after every intra distance of cluster k has been
            // through checked_distance. cluster_representative sorts each
            // candidate's distance vector in median_distance and then
            // stable_sorts the candidate scores; a NaN in either range makes
            // operator< a non-strict-weak ordering, which is the undefined
            // behaviour the finiteness precondition exists to replace with a
            // named error.
            const size_t representative =
                cluster_representative(cluster, storage, options.representative_method);
            representatives.push_back(representative);

            double radius = 0.0;
            double representative_total = 0.0;
            size_t representative_count = 0;
            for (const size_t member : cluster) {
                const double distance =
                    detail::checked_distance(storage, representative, member);
                radius = std::max(radius, distance);
                if (member != representative) {
                    representative_total += distance;
                    ++representative_count;
                }
            }
            radii.push_back(radius);
            medoid_member_means.push_back(
                representative_count == 0
                    ? 0.0
                    : representative_total / static_cast<double>(representative_count));

            // The medoid-named indices always use the true medoid, whatever
            // representative_method selected, because a field called
            // calinski_harabasz_medoid must not be a minimax number. Under the
            // default method the two coincide and the second selection is
            // skipped.
            const size_t medoid =
                options.representative_method == RepresentativeMethod::Medoid
                    ? representative
                    : cluster_representative(cluster, storage, RepresentativeMethod::Medoid);
            true_medoids.push_back(medoid);

            double scatter_total = 0.0;
            for (const size_t member : cluster) {
                const double distance = detail::checked_distance(storage, medoid, member);
                scatter_total += distance;
                medoid_square_sums[k] += distance * distance;
            }
            // The n_k denominator, counting the medoid's own zero distance. This
            // is deliberately not medoid_member_means, which divides by
            // n_k - 1: reusing that field would inflate both Davies-Bouldin and
            // the medoid Dunn variant by n_k / (n_k - 1), a factor of two at
            // n_k == 2.
            medoid_scatter[k] =
                cluster.empty() ? 0.0 : scatter_total / static_cast<double>(cluster.size());

            const double own_denominator = static_cast<double>(cluster.size()) - 1.0;
            for (const size_t member : cluster) {
                point_total[member] = own_mean[member];
                own_mean[member] =
                    cluster.size() > 1 ? own_mean[member] / own_denominator : 0.0;
            }
        }

        report.mean_intra_distance =
            intra_pairs.empty() ? nan_value() : detail::mean_distance(intra_pairs);
        report.median_intra_distance =
            intra_pairs.empty() ? nan_value() : detail::median_distance(intra_pairs);
        report.median_radius = detail::median_distance(radii);
        report.p95_diameter = percentile(diameters, 0.95);
        report.median_medoid_member_distance = detail::median_distance(medoid_member_means);

        // ---- Cross pass: once per unordered cluster pair. ----
        // Three separate walks -- boundary violations, the silhouette b term and
        // the Dunn separation -- become one. Every quantity below comes off the
        // same distance read.
        std::vector<double> best_other_mean(
            labels.size(), std::numeric_limits<double>::infinity());
        std::vector<double> pair_sum(labels.size(), 0.0);
        // K-length, not K x K. A K x K matrix would be O(N^2) memory whenever
        // the clusters are near-singletons -- exactly the HDBSCAN and Butina
        // shapes this report is run on -- and no consumer reads a non-minimal
        // entry: the record table wants only each cluster's nearest neighbour,
        // and the Dunn variant wants only the global minimum mean.
        std::vector<double> nearest_cluster_distance(
            cluster_count, std::numeric_limits<double>::infinity());
        // An entry is meaningful only where nearest_cluster_distance[k] is
        // finite. size_t has no natural sentinel here, so a cluster with no
        // neighbour keeps the 0 it was built with, which is indistinguishable
        // from a genuine answer of "ordinal 0" -- at K == 1 the loop below never
        // runs and cluster 0 would name itself. The consumer maps the infinite
        // distance to NO_NEAREST_CLUSTER rather than reading this vector alone.
        std::vector<size_t> nearest_cluster(cluster_count, 0);
        std::vector<size_t> cluster_violations(cluster_count, 0);
        double min_mean_separation = std::numeric_limits<double>::infinity();
        detail::DistanceMoments between_moments;

        // Only materialised under the flag. This is the one allocation here
        // that changes the allocation class: it completes the pair-array
        // footprint to Nc(Nc-1)/2 doubles, since intra_pairs above already
        // holds the within-cluster half unconditionally -- together roughly
        // 400 MB at Nc = 10,000.
        //
        // Its own share is only the between-cluster half, which is what it is
        // reserved to. That share is not a fixed fraction of the total: it is
        // zero when the clustering is a single cluster, and nearly all of it
        // when the clusters are small. Reserving the full C(Nc, 2) here would
        // over-allocate by the within-cluster half on every call.
        std::vector<double> between_distances;
        if (options.compute_pair_rank_indices) {
            between_distances.reserve(pair_count(clustered_count) - intra_pair_count);
        }

        size_t violations = 0;
        double min_inter = std::numeric_limits<double>::infinity();
        for (size_t a = 0; a < cluster_count; ++a) {
            for (size_t b = a + 1; b < cluster_count; ++b) {
                for (const size_t i : members[a]) {
                    pair_sum[i] = 0.0;
                }
                for (const size_t j : members[b]) {
                    pair_sum[j] = 0.0;
                }

                double pair_min = std::numeric_limits<double>::infinity();
                double pair_total = 0.0;
                size_t pair_violations = 0;
                for (const size_t i : members[a]) {
                    for (const size_t j : members[b]) {
                        const double distance = detail::checked_distance(storage, i, j);
                        between_moments.Add(distance);
                        if (options.compute_pair_rank_indices) {
                            between_distances.push_back(distance);
                        }
                        pair_min = std::min(pair_min, distance);
                        pair_total += distance;
                        if (distance <= options.boundary_threshold) {
                            ++pair_violations;
                        }
                        pair_sum[i] += distance;
                        pair_sum[j] += distance;
                    }
                }

                violations += pair_violations;
                min_inter = std::min(min_inter, pair_min);

                const double size_a = static_cast<double>(members[a].size());
                const double size_b = static_cast<double>(members[b].size());
                const double cross_pair_count = size_a * size_b;
                const double mean_cross =
                    cross_pair_count > 0.0 ? pair_total / cross_pair_count : 0.0;
                min_mean_separation = std::min(min_mean_separation, mean_cross);

                // This loop visits cluster k's partners in strictly ascending
                // ordinal order -- 0..k-1 while k is the inner b, then
                // k+1..K-1 while k is the outer a -- so a strict < keeps the
                // lowest-ordinal winner on a tie, matching the documented rule
                // without a second scan.
                if (pair_min < nearest_cluster_distance[a]) {
                    nearest_cluster_distance[a] = pair_min;
                    nearest_cluster[a] = b;
                }
                if (pair_min < nearest_cluster_distance[b]) {
                    nearest_cluster_distance[b] = pair_min;
                    nearest_cluster[b] = a;
                }
                cluster_violations[a] += pair_violations;
                cluster_violations[b] += pair_violations;

                for (const size_t i : members[a]) {
                    best_other_mean[i] = std::min(best_other_mean[i], pair_sum[i] / size_b);
                    point_total[i] += pair_sum[i];
                }
                for (const size_t j : members[b]) {
                    best_other_mean[j] = std::min(best_other_mean[j], pair_sum[j] / size_a);
                    point_total[j] += pair_sum[j];
                }
            }
        }
        report.boundary_violations = violations;

        if (cluster_count >= 2) {
            double silhouette_sum = 0.0;
            size_t silhouette_count = 0;
            for (const Cluster& cluster : members) {
                for (const size_t point : cluster) {
                    const double a_term = own_mean[point];
                    const double b_term = best_other_mean[point];
                    silhouette_sum += silhouette_term(a_term, b_term);
                    ++silhouette_count;
                }
            }
            report.silhouette = silhouette_count == 0
                ? nan_value()
                : silhouette_sum / static_cast<double>(silhouette_count);
            report.dunn_index =
                max_diameter == 0.0 ? nan_value() : min_inter / max_diameter;
        } else {
            report.silhouette = nan_value();
            report.dunn_index = nan_value();
        }

        // ---- Finalize. ----
        // Under the default representative_method the configured representative
        // and the true medoid are the same point, and the Davies-Bouldin loop
        // in the internal-indices block below walks every pair of them anyway.
        // Redundancy is the median of that walk's row minima, so it is taken
        // from there and this loop is skipped -- one K^2 matrix, not two.
        // A different method means genuinely different points, and then this
        // loop is the only thing that reads them.
        if (representatives.size() < 2) {
            report.representative_redundancy = nan_value();
        } else if (options.representative_method != RepresentativeMethod::Medoid) {
            std::vector<double> nearest_representative_distance;
            nearest_representative_distance.reserve(representatives.size());
            for (size_t i = 0; i < representatives.size(); ++i) {
                double smallest = std::numeric_limits<double>::infinity();
                for (size_t j = 0; j < representatives.size(); ++j) {
                    if (i != j) {
                        smallest = std::min(
                            smallest,
                            detail::checked_distance(
                                storage, representatives[i], representatives[j]));
                    }
                }
                nearest_representative_distance.push_back(smallest);
            }
            report.representative_redundancy =
                detail::median_distance(nearest_representative_distance);
        }

        // Each sample's distance to its nearest configured representative,
        // computed once. Every coverage threshold is then a scan of this vector,
        // O(T*N) rather than the O(T*N*K) the per-threshold recomputation cost.
        //
        // The threshold list must be non-empty for the scan to run at all. With
        // no thresholds every distance it reads feeds coverage_at and
        // noise_coverage_at that stay empty, so the only thing the scan can
        // still produce is a refusal -- and refusing a report whose every field
        // is already determined is the over-refusal of INVARIANT 3, not a precondition.
        std::vector<double> nearest_representative(
            report.num_samples, std::numeric_limits<double>::infinity());
        if (report.num_samples > 0 && !representatives.empty() &&
            !options.coverage_thresholds.empty()) {
            for (size_t point = 0; point < report.num_samples; ++point) {
                for (const size_t representative : representatives) {
                    nearest_representative[point] = std::min(
                        nearest_representative[point],
                        detail::checked_distance(storage, point, representative));
                }
            }

            report.coverage_at.assign(options.coverage_thresholds.size(), 0.0);
            report.noise_coverage_at.assign(
                options.coverage_thresholds.size(), nan_value());
            for (size_t t = 0; t < options.coverage_thresholds.size(); ++t) {
                const double threshold = options.coverage_thresholds[t];
                size_t covered = 0;
                size_t covered_noise = 0;
                for (size_t point = 0; point < report.num_samples; ++point) {
                    if (nearest_representative[point] <= threshold) {
                        ++covered;
                        if (owner[point] == NO_OWNER) {
                            ++covered_noise;
                        }
                    }
                }
                report.coverage_at[t] =
                    static_cast<double>(covered) / static_cast<double>(report.num_samples);
                // Left NaN when the clustering has no noise: 0.0 would read as
                // "no noise point is covered" rather than "no noise point exists".
                if (num_noise > 0) {
                    report.noise_coverage_at[t] =
                        static_cast<double>(covered_noise) /
                        static_cast<double>(num_noise);
                }
            }
        }

        // ---- Internal indices (section 5.3). ----
        // clustered_count is the one the pair-array reservations were derived
        // from, hoisted above them rather than recomputed here, so the between
        // -cluster reservation and this denominator cannot drift apart.
        if (cluster_count >= 2) {
            // The global medoid M: the clustered point with the smallest total
            // distance to all clustered points. The lowest-index tiebreak is
            // deliberately not m_k's "earliest member" rule -- Butina emits
            // members in representative-first order, so the two can name
            // different points.
            //
            // That tiebreak governs ties in the COMPUTED totals, which is a
            // weaker guarantee than it reads as. point_total is seeded from
            // the intra pass and then accumulated once per foreign cluster, so
            // two points whose exact rational totals are equal can still land
            // a bit apart and resolve either way, and a measured instance put
            // two such totals one ULP apart -- enough to move M.
            //
            // The consequence is not small. Calinski-Harabasz forms a
            // cluster-size-weighted sum of squared distances to M, and equal
            // unweighted totals put no constraint on that quantity: two
            // candidates can tie exactly on the first and sit far apart on the
            // second. On a measured five-point case whose totals are equal in
            // exact arithmetic, calinski_harabasz_medoid reads 11.54 for one
            // choice of M and 16.72 for the other -- a factor of 1.45 between
            // two selections the tiebreak regards as equally correct. An index
            // of this family that shifts after a change which only reordered
            // accumulation should be traced back to here first.
            //
            // Exact or compensated summation is not the remedy. It would make
            // the computed ties coincide with the exact ones, but the tiebreak
            // still has to choose at an exact tie, and an input perturbed by
            // one ULP would flip M with the same swing. The instability is
            // intrinsic to taking an argmin across a near-tie, not to the
            // summation scheme, so paying for exact summation on this hot
            // O(N^2) reduction would buy nothing. Hence documented, not fixed.
            size_t global_medoid = 0;
            double smallest_total = std::numeric_limits<double>::infinity();
            for (size_t i = 0; i < labels.size(); ++i) {
                if (owner[i] != NO_OWNER && point_total[i] < smallest_total) {
                    smallest_total = point_total[i];
                    global_medoid = i;
                }
            }

            double between_scatter = 0.0;
            double within_scatter = 0.0;
            for (size_t k = 0; k < cluster_count; ++k) {
                const double to_global =
                    detail::checked_distance(storage, true_medoids[k], global_medoid);
                between_scatter +=
                    static_cast<double>(members[k].size()) * to_global * to_global;
                within_scatter += medoid_square_sums[k];
            }
            if (clustered_count > cluster_count) {
                const double numerator =
                    between_scatter / static_cast<double>(cluster_count - 1);
                const double denominator =
                    within_scatter /
                    static_cast<double>(clustered_count - cluster_count);
                report.calinski_harabasz_medoid =
                    denominator == 0.0 ? nan_value() : numerator / denominator;
            } else {
                report.calinski_harabasz_medoid = nan_value();
            }

            // This is the only K^2 walk over medoid pairs, and three consumers
            // read it: Davies-Bouldin, the medoid Dunn separation, and -- under
            // the default representative_method, where the configured
            // representative and the true medoid are the same point --
            // representative_redundancy, whose own loop above is switched off in
            // that case. Section 5.1 promises exactly one such matrix by
            // default and a second one only when the two identities differ.
            double davies_bouldin_total = 0.0;
            std::vector<double> nearest_medoid_distance;
            nearest_medoid_distance.reserve(cluster_count);
            for (size_t a = 0; a < cluster_count; ++a) {
                double worst_ratio = 0.0;
                double nearest = std::numeric_limits<double>::infinity();
                for (size_t b = 0; b < cluster_count; ++b) {
                    if (a == b) {
                        continue;
                    }
                    const double separation =
                        detail::checked_distance(storage, true_medoids[a], true_medoids[b]);
                    nearest = std::min(nearest, separation);
                    // Coincident medoids are reported as inf, and the zero
                    // separation is branched on rather than divided by. IEEE
                    // gives 0.0/0.0 == NaN, and std::max(0.0, NaN) returns 0.0
                    // because 0.0 < NaN is false -- so two coincident
                    // zero-scatter clusters would silently score a perfect
                    // Davies-Bouldin of 0 through the division. INVARIANT 3.
                    const double ratio =
                        separation == 0.0
                            ? std::numeric_limits<double>::infinity()
                            : (medoid_scatter[a] + medoid_scatter[b]) / separation;
                    worst_ratio = std::max(worst_ratio, ratio);
                }
                nearest_medoid_distance.push_back(nearest);
                davies_bouldin_total += worst_ratio;
            }
            report.davies_bouldin_medoid =
                davies_bouldin_total / static_cast<double>(cluster_count);

            // The smallest medoid-to-medoid distance is the smallest row
            // minimum: the matrix is symmetric, so scanning ordered pairs and
            // scanning unordered ones give the same answer, and reusing the
            // row minima avoids a second condition inside the hot loop.
            const double min_medoid_separation = *std::min_element(
                nearest_medoid_distance.begin(), nearest_medoid_distance.end());

            // Same reads, different reduction: the median of the row minima is
            // exactly what the representative_redundancy loop computes, so
            // under the default method it is answered here for free. Under any
            // other method that loop has already answered it against the
            // configured representatives, which are different points.
            if (options.representative_method == RepresentativeMethod::Medoid) {
                report.representative_redundancy =
                    detail::median_distance(nearest_medoid_distance);
            }

            double max_mean_within = 0.0;
            double max_medoid_spread = 0.0;
            for (size_t k = 0; k < cluster_count; ++k) {
                // A singleton contributes 0 to both denominators, which cannot
                // raise a maximum and so cannot move either index.
                const double mean_within = cluster_intra_counts[k] == 0
                    ? 0.0
                    : cluster_intra_sums[k] /
                          static_cast<double>(cluster_intra_counts[k]);
                max_mean_within = std::max(max_mean_within, mean_within);
                max_medoid_spread = std::max(max_medoid_spread, 2.0 * medoid_scatter[k]);
            }
            // min_mean_separation was accumulated by the cross pass above.
            report.dunn_mean_separation_mean_diameter =
                max_mean_within == 0.0 ? nan_value()
                                       : min_mean_separation / max_mean_within;
            report.dunn_medoid_separation_medoid_spread =
                max_medoid_spread == 0.0 ? nan_value()
                                         : min_medoid_separation / max_medoid_spread;
        } else {
            report.calinski_harabasz_medoid = nan_value();
            report.davies_bouldin_medoid = nan_value();
            report.dunn_mean_separation_mean_diameter = nan_value();
            report.dunn_medoid_separation_medoid_spread = nan_value();
        }

        // Point-biserial. Positive means between-cluster distances exceed
        // within-cluster ones; the convention is pinned here because published
        // sources differ on it.
        const detail::DistanceMoments all_moments =
            detail::merge_moments(within_moments, between_moments);
        const double spread = detail::population_stddev(all_moments);
        if (within_moments.count == 0 || between_moments.count == 0 || spread == 0.0) {
            report.point_biserial = nan_value();
        } else {
            const double within_pairs = static_cast<double>(within_moments.count);
            const double between_pairs = static_cast<double>(between_moments.count);
            const double total_pairs = static_cast<double>(all_moments.count);
            report.point_biserial =
                ((between_moments.mean - within_moments.mean) / spread) *
                std::sqrt(within_pairs * between_pairs) / total_pairs;
        }

        // ---- Per-cluster records (section 4.4). ----
        if (options.compute_per_cluster_records) {
            report.records.reserve(cluster_count);
            for (size_t k = 0; k < cluster_count; ++k) {
                ClusterRecord record;
                record.label = static_cast<ClusterLabel>(k);
                record.size = members[k].size();
                record.representative = representatives[k];
                record.mean_intra_distance = cluster_intra_counts[k] == 0
                    ? nan_value()
                    : cluster_intra_sums[k] /
                          static_cast<double>(cluster_intra_counts[k]);
                record.median_intra_distance = cluster_medians[k];
                record.radius = radii[k];
                record.diameter = diameters[k];
                record.mean_representative_distance = medoid_member_means[k];

                if (cluster_count >= 2) {
                    // At K >= 2 the cross loop visits every cluster in at least
                    // one pair, every cluster is non-empty (the validator
                    // rejects empty clusters before this point), and every
                    // cross distance is finite (detail::checked_distance
                    // rejects a non-finite read at the cross pass -- the
                    // storage precondition checks only the backend type, not
                    // the values). So nearest_cluster_distance[k] is finite for
                    // every k, and the count guard coincides with the
                    // producer's finiteness guard. nearest_cluster[k] is
                    // assigned in the same if bodies as the distance, so a
                    // finite distance also means the ordinal is a real answer
                    // rather than the default 0.
                    record.nearest_cluster =
                        static_cast<ClusterLabel>(nearest_cluster[k]);
                    record.nearest_cluster_distance = nearest_cluster_distance[k];

                    double silhouette_total = 0.0;
                    for (const size_t point : members[k]) {
                        const double a_term = own_mean[point];
                        const double b_term = best_other_mean[point];
                        silhouette_total += silhouette_term(a_term, b_term);
                    }
                    record.silhouette = members[k].empty()
                        ? nan_value()
                        : silhouette_total / static_cast<double>(members[k].size());
                } else {
                    record.nearest_cluster = NO_NEAREST_CLUSTER;
                    // NaN rather than 0.0: a zero would read as another cluster
                    // sitting at zero distance.
                    record.nearest_cluster_distance = nan_value();
                    record.silhouette = nan_value();
                }

                record.boundary_violations = cluster_violations[k];

                report.records.push_back(record);
            }
        }

        // ---- Pair-rank indices (section 5.4), opt-in. ----
        if (options.compute_pair_rank_indices) {
            // Both arrays are moved, not copied. pair_rank_indices takes them
            // by value and sorts in place, so passing lvalues would hold a
            // second pair-sized copy alive alongside the originals. Each was
            // also reserved to its exact final count where it was declared, so
            // the fill does not peak above the final size on the way here
            // either. Nothing reads either array after this point; the
            // moved-from state is never observed.
            const detail::PairRankIndices pair_rank = detail::pair_rank_indices(
                std::move(intra_pairs), std::move(between_distances));
            report.c_index = pair_rank.c_index;
            report.baker_hubert_gamma = pair_rank.baker_hubert_gamma;
        } else {
            report.c_index = nan_value();
            report.baker_hubert_gamma = nan_value();
        }

        // Any stage added below this point must not read intra_pairs or
        // between_distances. Both are moved-from above, so a read returns an
        // empty array rather than failing, and the metric it feeds would be
        // a plausible wrong number. A stage that needs either array belongs
        // above this block.
    } else {
        report.mean_intra_distance = nan_value();
        report.median_intra_distance = nan_value();
        report.median_radius = nan_value();
        report.p95_diameter = nan_value();
        report.silhouette = nan_value();
        report.dunn_index = nan_value();
        report.median_medoid_member_distance = nan_value();
        report.representative_redundancy = nan_value();
        report.calinski_harabasz_medoid = nan_value();
        report.davies_bouldin_medoid = nan_value();
        report.dunn_mean_separation_mean_diameter = nan_value();
        report.dunn_medoid_separation_medoid_spread = nan_value();
        report.point_biserial = nan_value();
        report.c_index = nan_value();
        report.baker_hubert_gamma = nan_value();
        // coverage_at and noise_coverage_at stay empty: with no clusters there
        // are no representatives, so no coverage question has an answer.
    }

    return report;
}

ClusterReportComparison compare_reports(const ClusterReport& a, const ClusterReport& b) {
    ClusterReportComparison comparison;
    comparison.a = a;
    comparison.b = b;
    return comparison;
}

}  // namespace OECluster
