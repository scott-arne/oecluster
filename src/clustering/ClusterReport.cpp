/**
 * @file ClusterReport.cpp
 * @brief Clustering-quality report implementation.
 */

#include "oecluster/clustering/ClusterReport.h"

#include <algorithm>
#include <cmath>
#include <limits>

#include "ClusterMetrics.h"
#include "DistanceAccess.h"

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

// Mean distance from `point` to every member of `cluster` except itself.
double mean_to_cluster(
    const size_t point,
    const Cluster& cluster,
    const StorageBackend& storage) {
    double total = 0.0;
    size_t count = 0;
    for (const size_t member : cluster) {
        if (member != point) {
            total += storage.Get(point, member);
            ++count;
        }
    }
    return count == 0 ? 0.0 : total / static_cast<double>(count);
}

double min_inter_cluster_distance(
    const Cluster& a,
    const Cluster& b,
    const StorageBackend& storage) {
    double smallest = std::numeric_limits<double>::infinity();
    for (const size_t i : a) {
        for (const size_t j : b) {
            smallest = std::min(smallest, storage.Get(i, j));
        }
    }
    return smallest;
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

    if (!members.empty()) {
        // Intra-cluster pair distances, radii, diameters.
        std::vector<double> intra_pairs;
        std::vector<double> radii;
        std::vector<double> diameters;
        double max_diameter = 0.0;
        for (const Cluster& cluster : members) {
            for (size_t i = 0; i < cluster.size(); ++i) {
                for (size_t j = i + 1; j < cluster.size(); ++j) {
                    intra_pairs.push_back(storage.Get(cluster[i], cluster[j]));
                }
            }
            const size_t medoid =
                cluster_representative(cluster, storage, options.representative_method);
            double radius = 0.0;
            for (const size_t member : cluster) {
                radius = std::max(radius, storage.Get(medoid, member));
            }
            radii.push_back(radius);
            const double diameter = detail::cluster_diameter(cluster, storage);
            diameters.push_back(diameter);
            max_diameter = std::max(max_diameter, diameter);
        }

        report.mean_intra_distance =
            intra_pairs.empty() ? nan_value() : detail::mean_distance(intra_pairs);
        report.median_intra_distance =
            intra_pairs.empty() ? nan_value() : detail::median_distance(intra_pairs);
        report.median_radius = detail::median_distance(radii);
        report.p95_diameter = percentile(diameters, 0.95);

        // Boundary violations: cross-cluster member pairs within boundary_threshold.
        size_t violations = 0;
        for (size_t a = 0; a < members.size(); ++a) {
            for (size_t b = a + 1; b < members.size(); ++b) {
                for (const size_t i : members[a]) {
                    for (const size_t j : members[b]) {
                        if (storage.Get(i, j) <= options.boundary_threshold) {
                            ++violations;
                        }
                    }
                }
            }
        }
        report.boundary_violations = violations;

        if (members.size() >= 2) {
            // True mean per-point silhouette over clustered points.
            double silhouette_sum = 0.0;
            size_t silhouette_count = 0;
            for (size_t k = 0; k < members.size(); ++k) {
                for (const size_t point : members[k]) {
                    const double a = mean_to_cluster(point, members[k], storage);
                    double b = std::numeric_limits<double>::infinity();
                    for (size_t other = 0; other < members.size(); ++other) {
                        if (other != k) {
                            b = std::min(b, mean_to_cluster(point, members[other], storage));
                        }
                    }
                    const double denom = std::max(a, b);
                    silhouette_sum += denom == 0.0 ? 0.0 : (b - a) / denom;
                    ++silhouette_count;
                }
            }
            report.silhouette = silhouette_count == 0
                ? nan_value()
                : silhouette_sum / static_cast<double>(silhouette_count);

            // Dunn index: min inter-cluster distance / max intra diameter.
            double min_inter = std::numeric_limits<double>::infinity();
            for (size_t a = 0; a < members.size(); ++a) {
                for (size_t b = a + 1; b < members.size(); ++b) {
                    min_inter = std::min(
                        min_inter,
                        min_inter_cluster_distance(members[a], members[b], storage));
                }
            }
            report.dunn_index =
                max_diameter == 0.0 ? nan_value() : min_inter / max_diameter;
        } else {
            report.silhouette = nan_value();
            report.dunn_index = nan_value();
        }

        // Representative / coverage.
        std::vector<size_t> medoids;
        medoids.reserve(members.size());
        std::vector<double> medoid_member_means;
        medoid_member_means.reserve(members.size());
        for (const Cluster& cluster : members) {
            const size_t medoid =
                cluster_representative(cluster, storage, options.representative_method);
            medoids.push_back(medoid);
            medoid_member_means.push_back(mean_to_cluster(medoid, cluster, storage));
        }
        report.median_medoid_member_distance = detail::median_distance(medoid_member_means);

        if (medoids.size() >= 2) {
            std::vector<double> nearest_medoid;
            nearest_medoid.reserve(medoids.size());
            for (size_t i = 0; i < medoids.size(); ++i) {
                double smallest = std::numeric_limits<double>::infinity();
                for (size_t j = 0; j < medoids.size(); ++j) {
                    if (i != j) {
                        smallest = std::min(smallest, storage.Get(medoids[i], medoids[j]));
                    }
                }
                nearest_medoid.push_back(smallest);
            }
            report.representative_redundancy = detail::median_distance(nearest_medoid);
        } else {
            report.representative_redundancy = nan_value();
        }

        // Coverage: fraction of ALL samples within threshold of some medoid.
        if (report.num_samples > 0 && !medoids.empty()) {
            report.coverage_at.assign(options.coverage_thresholds.size(), 0.0);
            for (size_t t = 0; t < options.coverage_thresholds.size(); ++t) {
                const double threshold = options.coverage_thresholds[t];
                size_t covered = 0;
                for (size_t point = 0; point < report.num_samples; ++point) {
                    double nearest = std::numeric_limits<double>::infinity();
                    for (const size_t medoid : medoids) {
                        nearest = std::min(nearest, storage.Get(point, medoid));
                    }
                    if (nearest <= threshold) {
                        ++covered;
                    }
                }
                report.coverage_at[t] =
                    static_cast<double>(covered) / static_cast<double>(report.num_samples);
            }
        }
    } else {
        report.mean_intra_distance = nan_value();
        report.median_intra_distance = nan_value();
        report.median_radius = nan_value();
        report.p95_diameter = nan_value();
        report.silhouette = nan_value();
        report.dunn_index = nan_value();
        report.median_medoid_member_distance = nan_value();
        report.representative_redundancy = nan_value();
        // coverage_at stays empty.
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
