/**
 * @file report_test_support.h
 * @brief Bit-exact ClusterReport comparison shared by the report engine tests.
 */

#ifndef OECLUSTER_TESTS_CPP_REPORT_TEST_SUPPORT_H
#define OECLUSTER_TESTS_CPP_REPORT_TEST_SUPPORT_H

#include <gtest/gtest.h>

#include <cmath>
#include <cstddef>
#include <ios>
#include <sstream>
#include <string>
#include <vector>

#include "oecluster/clustering/ClusterReport.h"

namespace report_test {

// NaN equals NaN, and +0.0 differs from -0.0: the lazy path promises bit
// parity, and operator== would hide both a NaN and a sign flip.
inline bool SameDouble(double a, double b) {
    if (std::isnan(a) || std::isnan(b)) {
        return std::isnan(a) && std::isnan(b);
    }
    return a == b && std::signbit(a) == std::signbit(b);
}

inline std::string Hex(double value) {
    std::ostringstream out;
    out << std::hexfloat << value;
    return out.str();
}

inline ::testing::AssertionResult SameReport(const OECluster::ClusterReport& a,
                                             const OECluster::ClusterReport& b) {
    std::string diff;
    auto real = [&](const std::string& name, double x, double y) {
        if (!SameDouble(x, y)) {
            diff += name + ": " + Hex(x) + " vs " + Hex(y) + "\n";
        }
    };
    auto whole = [&](const std::string& name, long long x, long long y) {
        if (x != y) {
            diff += name + ": " + std::to_string(x) + " vs " + std::to_string(y) + "\n";
        }
    };
    auto reals = [&](const std::string& name, const std::vector<double>& x,
                     const std::vector<double>& y) {
        if (x.size() != y.size()) {
            diff += name + ": size " + std::to_string(x.size()) + " vs " +
                    std::to_string(y.size()) + "\n";
            return;
        }
        for (size_t t = 0; t < x.size(); ++t) {
            real(name + "[" + std::to_string(t) + "]", x[t], y[t]);
        }
    };

    whole("num_samples", a.num_samples, b.num_samples);
    whole("num_clusters", a.num_clusters, b.num_clusters);
    whole("num_noise", a.num_noise, b.num_noise);
    whole("num_singletons", a.num_singletons, b.num_singletons);
    real("noise_fraction", a.noise_fraction, b.noise_fraction);
    real("singleton_fraction", a.singleton_fraction, b.singleton_fraction);
    real("largest_cluster_fraction", a.largest_cluster_fraction, b.largest_cluster_fraction);
    real("cluster_size_median", a.cluster_size_median, b.cluster_size_median);
    real("cluster_size_p90", a.cluster_size_p90, b.cluster_size_p90);
    real("size_gini", a.size_gini, b.size_gini);
    real("size_entropy", a.size_entropy, b.size_entropy);
    real("mean_intra_distance", a.mean_intra_distance, b.mean_intra_distance);
    real("median_intra_distance", a.median_intra_distance, b.median_intra_distance);
    real("median_radius", a.median_radius, b.median_radius);
    real("p95_diameter", a.p95_diameter, b.p95_diameter);
    real("silhouette", a.silhouette, b.silhouette);
    real("dunn_index", a.dunn_index, b.dunn_index);
    whole("boundary_violations", a.boundary_violations, b.boundary_violations);
    real("median_medoid_member_distance", a.median_medoid_member_distance,
         b.median_medoid_member_distance);
    real("representative_redundancy", a.representative_redundancy,
         b.representative_redundancy);
    reals("coverage_thresholds", a.coverage_thresholds, b.coverage_thresholds);
    reals("coverage_at", a.coverage_at, b.coverage_at);
    reals("noise_coverage_at", a.noise_coverage_at, b.noise_coverage_at);
    real("calinski_harabasz_medoid", a.calinski_harabasz_medoid, b.calinski_harabasz_medoid);
    real("davies_bouldin_medoid", a.davies_bouldin_medoid, b.davies_bouldin_medoid);
    real("dunn_mean_separation_mean_diameter", a.dunn_mean_separation_mean_diameter,
         b.dunn_mean_separation_mean_diameter);
    real("dunn_medoid_separation_medoid_spread", a.dunn_medoid_separation_medoid_spread,
         b.dunn_medoid_separation_medoid_spread);
    real("point_biserial", a.point_biserial, b.point_biserial);
    real("c_index", a.c_index, b.c_index);
    real("baker_hubert_gamma", a.baker_hubert_gamma, b.baker_hubert_gamma);
    whole("requested.pair_rank_indices", a.requested.pair_rank_indices,
          b.requested.pair_rank_indices);
    whole("requested.per_cluster_records", a.requested.per_cluster_records,
          b.requested.per_cluster_records);

    whole("records.size", a.records.size(), b.records.size());
    for (size_t k = 0; k < a.records.size() && k < b.records.size(); ++k) {
        const OECluster::ClusterRecord& x = a.records[k];
        const OECluster::ClusterRecord& y = b.records[k];
        const std::string at = "records[" + std::to_string(k) + "].";
        whole(at + "label", x.label, y.label);
        whole(at + "size", x.size, y.size);
        whole(at + "representative", x.representative, y.representative);
        real(at + "mean_intra_distance", x.mean_intra_distance, y.mean_intra_distance);
        real(at + "median_intra_distance", x.median_intra_distance, y.median_intra_distance);
        real(at + "radius", x.radius, y.radius);
        real(at + "diameter", x.diameter, y.diameter);
        real(at + "mean_representative_distance", x.mean_representative_distance,
             y.mean_representative_distance);
        whole(at + "nearest_cluster", x.nearest_cluster, y.nearest_cluster);
        real(at + "nearest_cluster_distance", x.nearest_cluster_distance,
             y.nearest_cluster_distance);
        real(at + "silhouette", x.silhouette, y.silhouette);
        whole(at + "boundary_violations", x.boundary_violations, y.boundary_violations);
    }

    if (diff.empty()) {
        return ::testing::AssertionSuccess();
    }
    return ::testing::AssertionFailure() << diff;
}

}  // namespace report_test

#endif  // OECLUSTER_TESTS_CPP_REPORT_TEST_SUPPORT_H
