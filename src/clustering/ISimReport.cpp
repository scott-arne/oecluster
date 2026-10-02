/**
 * @file ISimReport.cpp
 * @brief iSIM set similarity and the approximate fingerprint-native cluster report.
 */

#include "oecluster/clustering/ISimReport.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#include "ClusterMetrics.h"
#include "ISimKernels.h"
#include "ReportCommon.h"
#include "oecluster/ThreadPool.h"

namespace OECluster {

ISimReportOptions::ISimReportOptions(const ClusterThreshold preset)
    : coverage_thresholds(detail::preset_coverage_thresholds(preset)) {}

namespace {

double nan_value() {
    return std::numeric_limits<double>::quiet_NaN();
}

// Shared by both entry points so they refuse the same inputs in the same
// order: metric, then zero width (as bitbirch refuses it), then width, then
// batch size.
void validate_isim_input(const std::string& metric, const OEFP::OEFPBatch& fingerprints,
                         const char* caller) {
    if (metric != "tanimoto") {
        throw std::invalid_argument(std::string(caller) +
                                    ": metric must be \"tanimoto\", got \"" + metric + "\"");
    }
    if (fingerprints.Size() > 0 && fingerprints.SizeBits() == 0) {
        throw std::invalid_argument(std::string(caller) +
                                    " requires non-zero-width fingerprints");
    }
    detail::check_isim_width(fingerprints.SizeBits());
    detail::check_isim_batch_size(fingerprints.Size());
}

// Per-cluster products of the core pass. Each worker writes only its own
// clusters' entries, and every reduction over them runs serially in ordinal
// order afterwards, so no value depends on the thread count.
struct ClusterCore {
    uint64_t size = 0;
    detail::CountMoments moments;
    detail::UInt128 dot_global;  ///< c_k . g, summed member by member.
    size_t medoid = 0;
    size_t global_candidate = 0;  ///< This cluster's best global-medoid candidate.
    detail::ISimScore global_score;
    double radius = 0.0;
    double mean_medoid_distance = 0.0;
    double medoid_scatter = 0.0;     ///< n_k denominator, for Davies-Bouldin and medoid Dunn.
    double medoid_square_sum = 0.0;  ///< For Calinski-Harabasz.
};

// A strictly better score, or an equal one at a lower sample index. Exact
// comparison, so the choice never rides on rounding.
bool better_candidate(const detail::ISimScore& score, const size_t index,
                      const detail::ISimScore& best_score, const size_t best_index) {
    if (detail::score_greater(score, best_score)) {
        return true;
    }
    return !detail::score_greater(best_score, score) && index < best_index;
}

double row_distance(const OEFP::OEFPBatch& fingerprints, const size_t lhs, const size_t rhs) {
    return detail::row_tanimoto_distance(fingerprints.RowWords(lhs), fingerprints.PopCount(lhs),
                                         fingerprints.RowWords(rhs), fingerprints.PopCount(rhs),
                                         fingerprints.SizeBits());
}

// Centroid-stage products per cluster, written only by the task that owns
// cluster k and reduced serially afterwards.
struct ClusterStage {
    double silhouette_sum = 0.0;
    ClusterLabel nearest = NO_NEAREST_CLUSTER;
    double nearest_similarity = std::numeric_limits<double>::quiet_NaN();
};

}  // namespace

double isim(const OEFP::OEFPBatch& fingerprints, const ISimOptions& options) {
    validate_isim_input(options.metric, fingerprints, "isim");
    const size_t num_fingerprints = fingerprints.Size();
    if (num_fingerprints < 2) {
        return nan_value();
    }
    const size_t size_bits = fingerprints.SizeBits();
    detail::BitCounts counts(size_bits, 0u);
    for (size_t row = 0; row < num_fingerprints; ++row) {
        detail::add_row_to_counts(counts, fingerprints.RowWords(row), size_bits);
    }
    const detail::ISimSums sums =
        detail::isim_sums(num_fingerprints, detail::count_moments(counts));
    return detail::isim_ratio(sums.intersections, sums.unions);
}

ISimReport isim_report(const ClusteringResult& result, const OEFP::OEFPBatch& fingerprints,
                       const ISimReportOptions& options) {
    validate_isim_input(options.metric, fingerprints, "isim_report");
    // Unlike cluster_report, refused at K == 0 as well: the native order puts
    // it ahead of every partition check, so one rule covers every call.
    for (size_t t = 0; t < options.coverage_thresholds.size(); ++t) {
        if (std::isnan(options.coverage_thresholds[t])) {
            throw std::invalid_argument("isim_report: coverage threshold " + std::to_string(t) +
                                        " must not be NaN");
        }
    }

    ISimReport report;
    report.coverage_thresholds = options.coverage_thresholds;
    report.requested.centroid_indices = options.compute_centroid_indices;
    report.requested.per_cluster_records = options.compute_per_cluster_records;
    detail::assign_profile(report,
                           detail::report_profile(result, options.treat_noise_as_singletons));

    const std::vector<size_t> owner = detail::validate_report_partition(
        result, fingerprints.Size(), "isim_report", "fingerprint batch");

    const Clusters& members = result.Members();
    const size_t cluster_count = members.size();
    if (cluster_count == 0) {
        return report;
    }
    const size_t size_bits = fingerprints.SizeBits();

    // g and S over clustered points only; noise never enters a quality field.
    // Serial: O(sum of popcounts), and integer, so order is immaterial.
    detail::BitCounts global_counts(size_bits, 0u);
    uint64_t clustered_count = 0;
    for (const Cluster& cluster : members) {
        for (const size_t member : cluster) {
            detail::add_row_to_counts(global_counts, fingerprints.RowWords(member), size_bits);
        }
        clustered_count += cluster.size();
    }
    const detail::CountMoments global_moments = detail::count_moments(global_counts);

    // K * bits * 4 bytes, allocated only when the stage runs.
    std::vector<uint32_t> columns;
    if (options.compute_centroid_indices) {
        columns.assign(size_bits * cluster_count, 0u);
    }

    std::vector<ClusterCore> core(cluster_count);
    // 0 is resolved to hardware concurrency before the cap, not left for
    // ThreadPool to expand: capping 0 first would leave it uncapped, and the
    // coverage scan's chunk count tracks N rather than K.
    const size_t requested_threads =
        options.num_threads > 0 ? options.num_threads
                                : std::max<size_t>(1, std::thread::hardware_concurrency());
    ThreadPool pool(std::min(requested_threads, cluster_count));
    pool.ParallelFor(0, cluster_count, 1, [&](const size_t begin, const size_t end) {
        detail::BitCounts counts(size_bits, 0u);
        for (size_t k = begin; k < end; ++k) {
            const Cluster& cluster = members[k];
            ClusterCore& out = core[k];
            std::fill(counts.begin(), counts.end(), 0u);
            for (const size_t member : cluster) {
                detail::add_row_to_counts(counts, fingerprints.RowWords(member), size_bits);
            }
            out.size = cluster.size();
            out.moments = detail::count_moments(counts);
            if (!columns.empty()) {
                // Each k belongs to exactly one task, so the strided writes
                // never collide.
                for (size_t bit = 0; bit < size_bits; ++bit) {
                    columns[bit * cluster_count + k] = counts[bit];
                }
            }

            detail::ISimScore best_score;
            bool first = true;
            for (const size_t member : cluster) {
                const uint64_t* row = fingerprints.RowWords(member);
                const uint32_t popcount = fingerprints.PopCount(member);
                const detail::ISimScore score = detail::isim_member_score(
                    detail::row_dot_counts(row, size_bits, counts), popcount, out.size,
                    out.moments.sum);
                const uint64_t dot_global = detail::row_dot_counts(row, size_bits, global_counts);
                out.dot_global += detail::UInt128::FromU64(dot_global);
                const detail::ISimScore global_score = detail::isim_member_score(
                    dot_global, popcount, clustered_count, global_moments.sum);
                if (first || better_candidate(score, member, best_score, out.medoid)) {
                    best_score = score;
                    out.medoid = member;
                }
                if (first || better_candidate(global_score, member, out.global_score,
                                              out.global_candidate)) {
                    out.global_score = global_score;
                    out.global_candidate = member;
                }
                first = false;
            }

            double radius = 0.0;
            double medoid_total = 0.0;
            double scatter_total = 0.0;
            double square_total = 0.0;
            for (const size_t member : cluster) {
                const double distance = row_distance(fingerprints, out.medoid, member);
                radius = std::max(radius, distance);
                if (member != out.medoid) {
                    medoid_total += distance;
                }
                scatter_total += distance;
                square_total += distance * distance;
            }
            out.radius = radius;
            out.mean_medoid_distance =
                out.size < 2 ? 0.0 : medoid_total / static_cast<double>(out.size - 1);
            // n_k, counting the medoid's own zero, as cluster_report's scatter does.
            out.medoid_scatter = scatter_total / static_cast<double>(out.size);
            out.medoid_square_sum = square_total;
        }
    });

    detail::UInt128 intra_intersections;
    detail::UInt128 intra_unions;
    detail::UInt128 square_sum_total;  // sum_k Q_k
    detail::UInt128 size_sum_total;    // sum_k n_k S_k
    bool has_intra_pair = false;
    std::vector<double> radii;
    std::vector<double> medoid_means;
    radii.reserve(cluster_count);
    medoid_means.reserve(cluster_count);
    for (const ClusterCore& cluster : core) {
        const detail::ISimSums sums = detail::isim_sums(cluster.size, cluster.moments);
        intra_intersections += sums.intersections;
        intra_unions += sums.unions;
        has_intra_pair = has_intra_pair || cluster.size >= 2;
        square_sum_total += cluster.moments.square_sum;
        size_sum_total += detail::multiply_u64(cluster.size, cluster.moments.sum);
        radii.push_back(cluster.radius);
        medoid_means.push_back(cluster.mean_medoid_distance);
    }
    report.isim_intra_distance =
        has_intra_pair ? 1.0 - detail::isim_ratio(intra_intersections, intra_unions) : nan_value();
    report.median_radius = detail::median_distance(radii);
    report.median_medoid_member_distance = detail::median_distance(medoid_means);

    std::vector<double> separation(cluster_count, nan_value());
    if (cluster_count >= 2) {
        detail::UInt128 global_square;  // g . g
        for (const uint32_t count : global_counts) {
            global_square += detail::multiply_u64(count, count);
        }
        // X = sum_{k<l} c_k . c_l and V = Nc S - sum_k n_k S_k - X.
        const detail::UInt128 cross_intersections = (global_square - square_sum_total).Half();
        const detail::UInt128 cross_unions =
            detail::multiply_u64(clustered_count, global_moments.sum) - size_sum_total -
            cross_intersections;
        report.isim_inter_distance = 1.0 - detail::isim_ratio(cross_intersections, cross_unions);

        for (size_t k = 0; k < cluster_count; ++k) {
            const ClusterCore& cluster = core[k];
            const detail::UInt128 outside_intersections =
                cluster.dot_global - cluster.moments.square_sum;
            const detail::UInt128 outside_unions =
                detail::multiply_u64(cluster.size, global_moments.sum - cluster.moments.sum) +
                detail::multiply_u64(clustered_count - cluster.size, cluster.moments.sum) -
                outside_intersections;
            separation[k] = 1.0 - detail::isim_ratio(outside_intersections, outside_unions);
        }

        size_t global_medoid = core[0].global_candidate;
        detail::ISimScore global_best = core[0].global_score;
        for (size_t k = 1; k < cluster_count; ++k) {
            if (better_candidate(core[k].global_score, core[k].global_candidate, global_best,
                                 global_medoid)) {
                global_best = core[k].global_score;
                global_medoid = core[k].global_candidate;
            }
        }

        double between_scatter = 0.0;
        double within_scatter = 0.0;
        for (const ClusterCore& cluster : core) {
            const double to_global = row_distance(fingerprints, cluster.medoid, global_medoid);
            between_scatter += static_cast<double>(cluster.size) * to_global * to_global;
            within_scatter += cluster.medoid_square_sum;
        }
        if (clustered_count > cluster_count) {
            const double numerator = between_scatter / static_cast<double>(cluster_count - 1);
            const double denominator =
                within_scatter / static_cast<double>(clustered_count - cluster_count);
            report.calinski_harabasz_medoid =
                denominator == 0.0 ? nan_value() : numerator / denominator;
        }
    }

    std::vector<ClusterStage> stage(cluster_count);
    if (options.compute_centroid_indices) {
        if (cluster_count >= 2) {
            pool.ParallelFor(0, cluster_count, 1, [&](const size_t begin, const size_t end) {
                std::vector<uint64_t> dots;
                std::vector<detail::UInt128> cross(cluster_count);
                for (size_t k = begin; k < end; ++k) {
                    const ClusterCore& own = core[k];
                    std::fill(cross.begin(), cross.end(), detail::UInt128{});
                    double silhouette_sum = 0.0;
                    for (const size_t member : members[k]) {
                        const uint32_t popcount = fingerprints.PopCount(member);
                        detail::row_dot_columns(fingerprints.RowWords(member), size_bits,
                                                columns, cluster_count, dots);
                        double a_term = 0.0;
                        if (own.size >= 2) {
                            a_term = 1.0 - detail::score_ratio(detail::isim_member_score(
                                               dots[k], popcount, own.size, own.moments.sum));
                        }
                        double b_term = std::numeric_limits<double>::infinity();
                        for (size_t l = 0; l < cluster_count; ++l) {
                            if (l == k) {
                                continue;
                            }
                            const detail::UInt128 shared = detail::UInt128::FromU64(dots[l]);
                            const detail::UInt128 unions =
                                detail::multiply_u64(core[l].size, popcount) +
                                detail::UInt128::FromU64(core[l].moments.sum) - shared;
                            b_term = std::min(b_term, 1.0 - detail::isim_ratio(shared, unions));
                        }
                        silhouette_sum += detail::silhouette_term(a_term, b_term, own.size);
                        for (size_t l = 0; l < cluster_count; ++l) {
                            cross[l] += detail::UInt128::FromU64(dots[l]);
                        }
                    }
                    stage[k].silhouette_sum = silhouette_sum;

                    // cross[l] is now c_k . c_l. Ratios are compared exactly, and
                    // only the winner is converted, so equal ratios tie and the
                    // strict comparison keeps the lowest ordinal.
                    detail::UInt128 best_shared;
                    detail::UInt128 best_unions;
                    for (size_t l = 0; l < cluster_count; ++l) {
                        if (l == k) {
                            continue;
                        }
                        const detail::UInt128 unions =
                            detail::multiply_u64(own.size, core[l].moments.sum) +
                            detail::multiply_u64(core[l].size, own.moments.sum) - cross[l];
                        if (stage[k].nearest == NO_NEAREST_CLUSTER ||
                            detail::isim_ratio_greater(cross[l], unions, best_shared,
                                                       best_unions)) {
                            best_shared = cross[l];
                            best_unions = unions;
                            stage[k].nearest = static_cast<ClusterLabel>(l);
                        }
                    }
                    stage[k].nearest_similarity = detail::isim_ratio(best_shared, best_unions);
                }
            });

            double silhouette_total = 0.0;
            for (const ClusterStage& entry : stage) {
                silhouette_total += entry.silhouette_sum;
            }
            report.isim_silhouette = silhouette_total / static_cast<double>(clustered_count);

            // K^2 exact medoid distances, one row per worker task, reduced
            // serially in ordinal order.
            std::vector<double> worst_ratio(cluster_count, 0.0);
            std::vector<double> nearest_medoid(cluster_count,
                                               std::numeric_limits<double>::infinity());
            pool.ParallelFor(0, cluster_count, 1, [&](const size_t begin, const size_t end) {
                for (size_t a = begin; a < end; ++a) {
                    for (size_t b = 0; b < cluster_count; ++b) {
                        if (a == b) {
                            continue;
                        }
                        const double separation_ab =
                            row_distance(fingerprints, core[a].medoid, core[b].medoid);
                        nearest_medoid[a] = std::min(nearest_medoid[a], separation_ab);
                        // Coincident medoids score inf, branched on rather than
                        // divided, as cluster_report does: 0/0 would be NaN and
                        // std::max(0.0, NaN) silently returns 0.
                        const double ratio =
                            separation_ab == 0.0
                                ? std::numeric_limits<double>::infinity()
                                : (core[a].medoid_scatter + core[b].medoid_scatter) / separation_ab;
                        worst_ratio[a] = std::max(worst_ratio[a], ratio);
                    }
                }
            });
            double davies_bouldin_total = 0.0;
            double max_medoid_spread = 0.0;
            for (size_t k = 0; k < cluster_count; ++k) {
                davies_bouldin_total += worst_ratio[k];
                max_medoid_spread = std::max(max_medoid_spread, 2.0 * core[k].medoid_scatter);
            }
            report.davies_bouldin_medoid =
                davies_bouldin_total / static_cast<double>(cluster_count);
            const double min_medoid_separation =
                *std::min_element(nearest_medoid.begin(), nearest_medoid.end());
            report.dunn_medoid_separation_medoid_spread =
                max_medoid_spread == 0.0 ? nan_value()
                                         : min_medoid_separation / max_medoid_spread;
        }

        // Every sample, noise included, against its nearest medoid: O(N K) exact reads.
        if (report.num_samples > 0 && !options.coverage_thresholds.empty()) {
            std::vector<double> nearest(report.num_samples,
                                        std::numeric_limits<double>::infinity());
            pool.ParallelFor(0, report.num_samples, 256,
                             [&](const size_t begin, const size_t end) {
                for (size_t point = begin; point < end; ++point) {
                    for (const ClusterCore& cluster : core) {
                        nearest[point] = std::min(nearest[point],
                                                  row_distance(fingerprints, point, cluster.medoid));
                    }
                }
            });
            const size_t threshold_count = options.coverage_thresholds.size();
            report.coverage_at.assign(threshold_count, 0.0);
            report.noise_coverage_at.assign(threshold_count, nan_value());
            for (size_t t = 0; t < threshold_count; ++t) {
                const double threshold = options.coverage_thresholds[t];
                size_t covered = 0;
                size_t covered_noise = 0;
                for (size_t point = 0; point < report.num_samples; ++point) {
                    if (nearest[point] <= threshold) {
                        ++covered;
                        if (owner[point] == detail::NO_OWNER) {
                            ++covered_noise;
                        }
                    }
                }
                report.coverage_at[t] =
                    static_cast<double>(covered) / static_cast<double>(report.num_samples);
                // NaN, not 0.0, when there is no noise: 0.0 would read as "no
                // noise point is covered".
                if (report.num_noise > 0) {
                    report.noise_coverage_at[t] =
                        static_cast<double>(covered_noise) / static_cast<double>(report.num_noise);
                }
            }
        }
    }

    if (options.compute_per_cluster_records) {
        report.records.reserve(cluster_count);
        for (size_t k = 0; k < cluster_count; ++k) {
            const ClusterCore& cluster = core[k];
            ISimClusterRecord record;
            record.label = static_cast<ClusterLabel>(k);
            record.size = static_cast<size_t>(cluster.size);
            record.medoid = cluster.medoid;
            if (cluster.size >= 2) {
                const detail::ISimSums sums = detail::isim_sums(cluster.size, cluster.moments);
                record.isim_intra_distance =
                    1.0 - detail::isim_ratio(sums.intersections, sums.unions);
            }
            record.isim_separation = separation[k];
            record.radius = cluster.radius;
            record.mean_medoid_distance = cluster.mean_medoid_distance;
            if (options.compute_centroid_indices && cluster_count >= 2) {
                record.isim_silhouette =
                    stage[k].silhouette_sum / static_cast<double>(cluster.size);
                record.nearest_cluster = stage[k].nearest;
                record.nearest_cluster_similarity = stage[k].nearest_similarity;
            }
            report.records.push_back(record);
        }
    }
    return report;
}

}  // namespace OECluster
