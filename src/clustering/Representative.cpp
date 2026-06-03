/**
 * @file Representative.cpp
 * @brief Representative member selection and ranking for clusters.
 */

#include "oecluster/clustering/Representative.h"

#include "ClusterMetrics.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <unordered_set>

namespace OECluster {

namespace {

double missing_metric() {
    return std::numeric_limits<double>::quiet_NaN();
}

bool has_threshold(const double threshold) {
    return threshold >= 0.0;
}

void validate_cluster_members(const Cluster& cluster, const size_t num_samples) {
    if (cluster.empty()) {
        throw std::invalid_argument("Cluster representative requires at least one member");
    }

    std::unordered_set<size_t> seen;
    seen.reserve(cluster.size());
    for (const size_t member : cluster) {
        if (member >= num_samples) {
            throw std::out_of_range("Cluster member index is outside the storage range");
        }
        if (!seen.insert(member).second) {
            throw std::invalid_argument("Cluster members must be unique");
        }
    }
}

void validate_complete_distance_storage(
    const Cluster& cluster,
    const StorageBackend& storage) {
    if (cluster.size() <= 1) {
        return;
    }
    if (dynamic_cast<const SparseStorage*>(&storage) != nullptr) {
        throw std::invalid_argument(
            "SparseStorage cannot provide complete distances for representative methods");
    }
}

void validate_vector_size(
    const size_t actual,
    const size_t expected,
    const std::string& name) {
    if (actual != 0 && actual != expected) {
        throw std::invalid_argument(
            name + " must be empty or the same length as the distance matrix");
    }
}

void validate_options(const RepresentativeOptions& options) {
    if (options.method == RepresentativeMethod::HighestNeighborhood &&
        !has_threshold(options.neighbor_threshold)) {
        throw std::invalid_argument(
            "Highest-neighborhood representative requires a non-negative neighbor threshold");
    }
}

std::vector<double> in_cluster_distances(
    const size_t candidate,
    const Cluster& cluster,
    const StorageBackend& storage) {
    std::vector<double> distances;
    distances.reserve(cluster.size() > 0 ? cluster.size() - 1 : 0);
    for (const size_t member : cluster) {
        if (member != candidate) {
            distances.push_back(storage.Get(candidate, member));
        }
    }
    return distances;
}

std::vector<bool> cluster_membership(const Cluster& cluster, const size_t n_items) {
    std::vector<bool> in_cluster(n_items, false);
    for (const size_t member : cluster) {
        in_cluster[member] = true;
    }
    return in_cluster;
}

double nearest_external_distance(
    const size_t candidate,
    const StorageBackend& storage,
    const std::vector<bool>& in_cluster) {
    double nearest = std::numeric_limits<double>::infinity();
    for (size_t item = 0; item < storage.NumSamples(); ++item) {
        if (!in_cluster[item]) {
            nearest = std::min(nearest, storage.Get(candidate, item));
        }
    }
    return nearest;
}

double neighbor_fraction_at_threshold(
    const size_t candidate,
    const Cluster& cluster,
    const StorageBackend& storage,
    const double threshold) {
    if (!has_threshold(threshold)) {
        return missing_metric();
    }

    size_t neighbors = 0;
    for (const size_t member : cluster) {
        if (storage.Get(candidate, member) <= threshold) {
            ++neighbors;
        }
    }
    return static_cast<double>(neighbors) / static_cast<double>(cluster.size());
}

double scaffold_purity(
    const size_t candidate,
    const Cluster& cluster,
    const std::vector<std::string>& scaffold_labels) {
    if (scaffold_labels.empty()) {
        return missing_metric();
    }

    const std::string& representative_scaffold = scaffold_labels[candidate];
    size_t matching = 0;
    for (const size_t member : cluster) {
        if (scaffold_labels[member] == representative_scaffold) {
            ++matching;
        }
    }
    return static_cast<double>(matching) / static_cast<double>(cluster.size());
}

double silhouette_like_score(
    const double mean_in_cluster_distance,
    const double nearest_external) {
    if (std::isinf(nearest_external)) {
        return 1.0;
    }

    const double denominator = std::max(mean_in_cluster_distance, nearest_external);
    if (denominator == 0.0) {
        return 0.0;
    }
    return (nearest_external - mean_in_cluster_distance) / denominator;
}

RepresentativeMetrics make_metrics(
    const size_t candidate,
    const Cluster& cluster,
    const StorageBackend& storage,
    const RepresentativeOptions& options,
    const double diameter,
    const std::vector<bool>& in_cluster) {
    const std::vector<double> distances = in_cluster_distances(candidate, cluster, storage);

    RepresentativeMetrics metrics;
    metrics.mean_distance_to_cluster = detail::mean_distance(distances);
    metrics.max_distance_to_cluster = detail::max_distance(distances);
    metrics.median_distance_to_cluster = detail::median_distance(distances);
    metrics.neighbor_fraction_at_threshold = neighbor_fraction_at_threshold(
        candidate,
        cluster,
        storage,
        options.neighbor_threshold);
    metrics.nearest_external_distance =
        nearest_external_distance(candidate, storage, in_cluster);
    metrics.cluster_radius = metrics.max_distance_to_cluster;
    metrics.cluster_diameter = diameter;
    metrics.silhouette_like_score = silhouette_like_score(
        metrics.mean_distance_to_cluster,
        metrics.nearest_external_distance);
    metrics.scaffold_purity = scaffold_purity(candidate, cluster, options.scaffold_labels);
    return metrics;
}

double optional_value(const std::vector<double>& values, const size_t member) {
    return values.empty() ? 0.0 : values[member];
}

double representative_score(
    const size_t member,
    const RepresentativeMetrics& metrics,
    const RepresentativeOptions& options) {
    switch (options.method) {
        case RepresentativeMethod::Medoid:
            return metrics.mean_distance_to_cluster;
        case RepresentativeMethod::Minimax:
            return metrics.max_distance_to_cluster;
        case RepresentativeMethod::HighestNeighborhood:
            return -metrics.neighbor_fraction_at_threshold;
        case RepresentativeMethod::WeightedMedoid:
            return (
                options.weights.alpha * metrics.mean_distance_to_cluster +
                options.weights.beta * optional_value(options.liability_penalties, member) -
                options.weights.gamma * optional_value(options.priority_scores, member));
    }

    throw std::invalid_argument("Unknown representative method");
}

}  // namespace

std::vector<ClusterRepresentative> rank_representatives(
    const Cluster& cluster,
    const StorageBackend& storage,
    const RepresentativeOptions& options) {
    validate_cluster_members(cluster, storage.NumSamples());
    validate_complete_distance_storage(cluster, storage);
    validate_options(options);
    validate_vector_size(
        options.liability_penalties.size(),
        storage.NumSamples(),
        "liability_penalties");
    validate_vector_size(
        options.priority_scores.size(),
        storage.NumSamples(),
        "priority_scores");
    validate_vector_size(
        options.scaffold_labels.size(),
        storage.NumSamples(),
        "scaffold_labels");

    const double diameter = detail::cluster_diameter(cluster, storage);
    const std::vector<bool> in_cluster = cluster_membership(cluster, storage.NumSamples());

    std::vector<ClusterRepresentative> ranked;
    ranked.reserve(cluster.size());
    for (const size_t member : cluster) {
        ClusterRepresentative representative;
        representative.member = member;
        representative.metrics = make_metrics(
            member,
            cluster,
            storage,
            options,
            diameter,
            in_cluster);
        representative.score =
            representative_score(member, representative.metrics, options);
        ranked.push_back(representative);
    }

    std::stable_sort(
        ranked.begin(),
        ranked.end(),
        [](const ClusterRepresentative& lhs, const ClusterRepresentative& rhs) {
            return lhs.score < rhs.score;
        });

    for (size_t rank = 0; rank < ranked.size(); ++rank) {
        ranked[rank].metrics.representative_rank = rank + 1;
    }
    return ranked;
}

std::vector<ClusterRepresentative> select_representatives(
    const Cluster& cluster,
    const StorageBackend& storage,
    const size_t k,
    const RepresentativeOptions& options) {
    if (k == 0) {
        return {};
    }

    const std::vector<ClusterRepresentative> ranked =
        rank_representatives(cluster, storage, options);
    const size_t limit = std::min(k, ranked.size());
    if (options.selection == RepresentativeSelection::Score || limit <= 1) {
        return std::vector<ClusterRepresentative>(ranked.begin(), ranked.begin() + limit);
    }

    if (options.selection != RepresentativeSelection::Diversity) {
        throw std::invalid_argument("Unknown representative selection strategy");
    }

    std::vector<ClusterRepresentative> selected;
    selected.reserve(limit);
    std::vector<bool> used(ranked.size(), false);
    selected.push_back(ranked.front());
    used.front() = true;

    while (selected.size() < limit) {
        size_t best_index = ranked.size();
        double best_nearest_selected_distance = -1.0;
        for (size_t candidate_index = 0; candidate_index < ranked.size(); ++candidate_index) {
            if (used[candidate_index]) {
                continue;
            }

            double nearest_selected_distance =
                std::numeric_limits<double>::infinity();
            for (const ClusterRepresentative& representative : selected) {
                nearest_selected_distance = std::min(
                    nearest_selected_distance,
                    storage.Get(ranked[candidate_index].member, representative.member));
            }

            // Iterating score-ranked candidates preserves score order for diversity ties.
            if (nearest_selected_distance > best_nearest_selected_distance) {
                best_index = candidate_index;
                best_nearest_selected_distance = nearest_selected_distance;
            }
        }

        selected.push_back(ranked.at(best_index));
        used[best_index] = true;
    }

    return selected;
}

size_t cluster_representative(
    const Cluster& cluster,
    const StorageBackend& storage,
    const RepresentativeOptions& options) {
    const std::vector<ClusterRepresentative> ranked =
        select_representatives(cluster, storage, 1, options);
    return ranked.front().member;
}

size_t cluster_representative(
    const Cluster& cluster,
    const StorageBackend& storage,
    const RepresentativeMethod method) {
    RepresentativeOptions options;
    options.method = method;
    return cluster_representative(cluster, storage, options);
}

}  // namespace OECluster
