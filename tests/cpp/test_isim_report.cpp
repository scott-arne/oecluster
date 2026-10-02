#include <gtest/gtest.h>

#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "isim_test_support.h"
#include "oecluster/StorageBackend.h"
#include "oecluster/clustering/ClusterReport.h"
#include "oecluster/clustering/ISimReport.h"

using namespace OECluster;

namespace {

ClusteringResult make_result(std::vector<ClusterLabel> labels) {
    Clusters members = labels_to_clusters(labels);
    return ClusteringResult(std::move(labels), std::move(members));
}

// 40 rows over 70 bits (a tail word), six clusters plus noise and a singleton.
std::vector<int> mixed_labels() {
    std::vector<int> labels;
    for (int i = 0; i < 40; ++i) labels.push_back(i % 7 - 1);
    labels[39] = 6;
    return labels;
}

ISimReportOptions records_options(bool centroid) {
    ISimReportOptions options;
    options.compute_per_cluster_records = true;
    options.compute_centroid_indices = centroid;
    return options;
}

OECluster::DenseStorage tanimoto_storage(const OEFP::OEFPBatch& batch) {
    OECluster::DenseStorage storage(batch.Size());
    for (size_t i = 0; i < batch.Size(); ++i)
        for (size_t j = i + 1; j < batch.Size(); ++j)
            storage.Set(i, j, isim_test::distance(batch, i, j));
    return storage;
}

// Compares every field except singleton_fraction, the one field
// treat_noise_as_singletons is allowed to change.
void expect_same_report(const ISimReport& a, const ISimReport& b) {
    using isim_test::same_double;
    EXPECT_EQ(a.num_samples, b.num_samples);
    EXPECT_EQ(a.num_clusters, b.num_clusters);
    EXPECT_EQ(a.num_noise, b.num_noise);
    EXPECT_EQ(a.num_singletons, b.num_singletons);
    EXPECT_TRUE(same_double(a.noise_fraction, b.noise_fraction));
    EXPECT_TRUE(same_double(a.largest_cluster_fraction, b.largest_cluster_fraction));
    EXPECT_TRUE(same_double(a.cluster_size_median, b.cluster_size_median));
    EXPECT_TRUE(same_double(a.cluster_size_p90, b.cluster_size_p90));
    EXPECT_TRUE(same_double(a.size_gini, b.size_gini));
    EXPECT_TRUE(same_double(a.size_entropy, b.size_entropy));
    EXPECT_EQ(a.coverage_thresholds, b.coverage_thresholds);
    EXPECT_TRUE(same_double(a.isim_intra_distance, b.isim_intra_distance));
    EXPECT_TRUE(same_double(a.isim_inter_distance, b.isim_inter_distance));
    EXPECT_TRUE(same_double(a.median_radius, b.median_radius));
    EXPECT_TRUE(same_double(a.median_medoid_member_distance, b.median_medoid_member_distance));
    EXPECT_TRUE(same_double(a.calinski_harabasz_medoid, b.calinski_harabasz_medoid));
    EXPECT_TRUE(same_double(a.isim_silhouette, b.isim_silhouette));
    EXPECT_TRUE(same_double(a.davies_bouldin_medoid, b.davies_bouldin_medoid));
    EXPECT_TRUE(same_double(a.dunn_medoid_separation_medoid_spread,
                            b.dunn_medoid_separation_medoid_spread));
    ASSERT_EQ(a.coverage_at.size(), b.coverage_at.size());
    for (size_t t = 0; t < a.coverage_at.size(); ++t) {
        EXPECT_TRUE(same_double(a.coverage_at[t], b.coverage_at[t]));
        EXPECT_TRUE(same_double(a.noise_coverage_at[t], b.noise_coverage_at[t]));
    }
    ASSERT_EQ(a.records.size(), b.records.size());
    for (size_t k = 0; k < a.records.size(); ++k) {
        EXPECT_EQ(a.records[k].label, b.records[k].label);
        EXPECT_EQ(a.records[k].size, b.records[k].size);
        EXPECT_EQ(a.records[k].medoid, b.records[k].medoid);
        EXPECT_TRUE(same_double(a.records[k].isim_intra_distance, b.records[k].isim_intra_distance));
        EXPECT_TRUE(same_double(a.records[k].isim_separation, b.records[k].isim_separation));
        EXPECT_TRUE(same_double(a.records[k].radius, b.records[k].radius));
        EXPECT_TRUE(same_double(a.records[k].mean_medoid_distance, b.records[k].mean_medoid_distance));
        EXPECT_TRUE(same_double(a.records[k].isim_silhouette, b.records[k].isim_silhouette));
        EXPECT_EQ(a.records[k].nearest_cluster, b.records[k].nearest_cluster);
        EXPECT_TRUE(same_double(a.records[k].nearest_cluster_similarity,
                                b.records[k].nearest_cluster_similarity));
    }
}

}  // namespace

TEST(ISimReportOptionsTest, PresetsSeedTheReportThresholds) {
    EXPECT_EQ(ISimReportOptions().coverage_thresholds, (std::vector<double>{0.25, 0.35, 0.45}));
    EXPECT_EQ(ISimReportOptions(ClusterThreshold::Tight).coverage_thresholds,
              (std::vector<double>{0.20, 0.30, 0.40}));
    EXPECT_EQ(ISimReportOptions(ClusterThreshold::Diversity).coverage_thresholds,
              (std::vector<double>{0.40, 0.50, 0.60}));
}

TEST(ISimReportTest, CoreMatchesTheOracleWithNoise) {
    const auto labels = mixed_labels();
    const auto batch = isim_test::make_random_batch(40, 70, 2024u);
    const auto oracle = isim_test::oracle_report(batch, labels, {});
    const auto report = isim_report(make_result(labels), batch, records_options(false));

    isim_test::expect_rel(report.isim_intra_distance, oracle.isim_intra_distance, "intra");
    isim_test::expect_rel(report.isim_inter_distance, oracle.isim_inter_distance, "inter");
    isim_test::expect_rel(report.median_radius, oracle.median_radius, "median_radius");
    isim_test::expect_rel(report.median_medoid_member_distance,
                          oracle.median_medoid_member_distance, "median_medoid_member");
    isim_test::expect_rel(report.calinski_harabasz_medoid, oracle.calinski_harabasz_medoid, "ch");
    ASSERT_EQ(report.records.size(), oracle.medoids.size());
    for (size_t k = 0; k < report.records.size(); ++k) {
        const auto& record = report.records[k];
        const std::string at = " k=" + std::to_string(k);
        EXPECT_EQ(record.label, static_cast<ClusterLabel>(k));
        EXPECT_EQ(record.medoid, oracle.medoids[k]) << at;
        isim_test::expect_rel(record.isim_intra_distance, oracle.record_intra[k], "record intra" + at);
        isim_test::expect_rel(record.isim_separation, oracle.separation[k], "separation" + at);
        isim_test::expect_rel(record.radius, oracle.radius[k], "radius" + at);
        isim_test::expect_rel(record.mean_medoid_distance, oracle.mean_medoid_distance[k], "mean" + at);
    }
    EXPECT_EQ(report.num_noise, 6u);
}

TEST(ISimReportTest, ProfileMatchesClusterReport) {
    const auto labels = mixed_labels();
    const auto batch = isim_test::make_random_batch(40, 70, 2024u);
    const auto exact = cluster_report(make_result(labels), tanimoto_storage(batch));
    const auto approx = isim_report(make_result(labels), batch);
    EXPECT_EQ(approx.num_samples, exact.num_samples);
    EXPECT_EQ(approx.num_clusters, exact.num_clusters);
    EXPECT_EQ(approx.num_noise, exact.num_noise);
    EXPECT_EQ(approx.num_singletons, exact.num_singletons);
    EXPECT_EQ(approx.noise_fraction, exact.noise_fraction);
    EXPECT_EQ(approx.singleton_fraction, exact.singleton_fraction);
    EXPECT_EQ(approx.largest_cluster_fraction, exact.largest_cluster_fraction);
    EXPECT_EQ(approx.cluster_size_median, exact.cluster_size_median);
    EXPECT_EQ(approx.cluster_size_p90, exact.cluster_size_p90);
    EXPECT_EQ(approx.size_gini, exact.size_gini);
    EXPECT_EQ(approx.size_entropy, exact.size_entropy);
}

TEST(ISimReportTest, MeanMedoidDistanceIsAtLeastTheTrueMedoids) {
    const auto labels = mixed_labels();
    const auto batch = isim_test::make_random_batch(40, 70, 77u);
    ClusterReportOptions exact_options;
    exact_options.compute_per_cluster_records = true;
    exact_options.representative_method = RepresentativeMethod::Medoid;
    const auto exact = cluster_report(make_result(labels), tanimoto_storage(batch), exact_options);
    const auto approx = isim_report(make_result(labels), batch, records_options(false));
    ASSERT_EQ(exact.records.size(), approx.records.size());
    for (size_t k = 0; k < approx.records.size(); ++k) {
        EXPECT_GE(approx.records[k].mean_medoid_distance,
                  exact.records[k].mean_representative_distance - 1e-12) << "k=" << k;
    }
}

TEST(ISimReportTest, TreatNoiseAsSingletonsChangesOnlySingletonFraction) {
    const auto labels = mixed_labels();
    const auto batch = isim_test::make_random_batch(40, 70, 5u);
    ISimReportOptions on = records_options(true);
    ISimReportOptions off = records_options(true);
    off.treat_noise_as_singletons = false;
    const auto a = isim_report(make_result(labels), batch, on);
    const auto b = isim_report(make_result(labels), batch, off);
    EXPECT_NE(a.singleton_fraction, b.singleton_fraction);
    expect_same_report(a, b);
}

TEST(ISimReportTest, ThreadCountDoesNotChangeTheCore) {
    const auto labels = mixed_labels();
    const auto batch = isim_test::make_random_batch(40, 70, 31u);
    ISimReportOptions options = records_options(false);
    options.num_threads = 1;
    const auto one = isim_report(make_result(labels), batch, options);
    for (const size_t threads : {2u, 8u}) {
        options.num_threads = threads;
        expect_same_report(one, isim_report(make_result(labels), batch, options));
    }
}

TEST(ISimReportTest, WellSeparatedBeatsShuffledLabels) {
    // Three blocks of near-duplicates on disjoint bit ranges.
    std::vector<OEFP::OEFP> fps;
    std::vector<int> labels;
    for (int block = 0; block < 3; ++block)
        for (int i = 0; i < 6; ++i) {
            const size_t base = static_cast<size_t>(block) * 40u;
            fps.push_back(isim_test::make_fp(128, {base, base + 1, base + 2, base + 3,
                                                   base + 4 + static_cast<size_t>(i)}));
            labels.push_back(block);
        }
    const auto batch = isim_test::make_batch(fps);
    std::vector<int> shuffled;
    for (size_t i = 0; i < labels.size(); ++i) shuffled.push_back(static_cast<int>(i % 3));
    ISimReportOptions options;
    options.compute_centroid_indices = true;
    const auto good = isim_report(make_result(labels), batch, options);
    const auto bad = isim_report(make_result(shuffled), batch, options);
    EXPECT_LT(good.isim_intra_distance, bad.isim_intra_distance);
    EXPECT_GT(good.isim_inter_distance, bad.isim_inter_distance);
}

TEST(ISimReportTest, NoClustersLeavesEveryQualityFieldNaN) {
    const auto batch = isim_test::make_random_batch(4, 16, 3u);
    const auto report = isim_report(make_result({-1, -1, -1, -1}), batch, records_options(true));
    EXPECT_EQ(report.num_clusters, 0u);
    EXPECT_EQ(report.num_noise, 4u);
    EXPECT_TRUE(std::isnan(report.isim_intra_distance));
    EXPECT_TRUE(std::isnan(report.isim_inter_distance));
    EXPECT_TRUE(std::isnan(report.median_radius));
    EXPECT_TRUE(std::isnan(report.median_medoid_member_distance));
    EXPECT_TRUE(std::isnan(report.calinski_harabasz_medoid));
    EXPECT_TRUE(std::isnan(report.isim_silhouette));
    EXPECT_TRUE(std::isnan(report.davies_bouldin_medoid));
    EXPECT_TRUE(std::isnan(report.dunn_medoid_separation_medoid_spread));
    EXPECT_TRUE(report.coverage_at.empty());
    EXPECT_TRUE(report.noise_coverage_at.empty());
    EXPECT_TRUE(report.records.empty());
    EXPECT_EQ(report.coverage_thresholds, (std::vector<double>{0.25, 0.35, 0.45}));
    EXPECT_TRUE(report.requested.centroid_indices);
    EXPECT_TRUE(report.requested.per_cluster_records);
}

TEST(ISimReportTest, OneClusterCoreDegenerates) {
    const auto batch = isim_test::make_random_batch(5, 32, 9u);
    const auto report = isim_report(make_result({0, 0, 0, 0, -1}), batch, records_options(false));
    EXPECT_FALSE(std::isnan(report.isim_intra_distance));
    EXPECT_TRUE(std::isnan(report.isim_inter_distance));
    EXPECT_TRUE(std::isnan(report.calinski_harabasz_medoid));
    ASSERT_EQ(report.records.size(), 1u);
    EXPECT_TRUE(std::isnan(report.records[0].isim_separation));
    EXPECT_EQ(report.records[0].nearest_cluster, NO_NEAREST_CLUSTER);
}

TEST(ISimReportTest, AllSingletonsCoreDegenerates) {
    const auto batch = isim_test::make_random_batch(4, 32, 11u);
    const auto report = isim_report(make_result({0, 1, 2, 3}), batch, records_options(false));
    EXPECT_TRUE(std::isnan(report.isim_intra_distance));
    EXPECT_FALSE(std::isnan(report.isim_inter_distance));
    EXPECT_TRUE(std::isnan(report.calinski_harabasz_medoid));
    for (const auto& record : report.records) {
        EXPECT_TRUE(std::isnan(record.isim_intra_distance));
        EXPECT_EQ(record.radius, 0.0);
        EXPECT_EQ(record.mean_medoid_distance, 0.0);
    }
}

TEST(ISimReportTest, AllZeroClusterIsAtDistanceZero) {
    const auto batch = isim_test::make_batch({
        isim_test::make_fp(16, {}), isim_test::make_fp(16, {}), isim_test::make_fp(16, {}),
        isim_test::make_fp(16, {1, 2}), isim_test::make_fp(16, {1, 3})});
    const auto report = isim_report(make_result({0, 0, 0, 1, 1}), batch, records_options(false));
    ASSERT_EQ(report.records.size(), 2u);
    EXPECT_EQ(report.records[0].isim_intra_distance, 0.0);
    EXPECT_EQ(report.records[0].medoid, 0u);  // tie to the lowest sample index
    EXPECT_EQ(report.records[0].radius, 0.0);
    EXPECT_EQ(report.records[0].mean_medoid_distance, 0.0);
}

TEST(ISimReportTest, StageFieldsStayEmptyWithoutTheStage) {
    const auto batch = isim_test::make_random_batch(40, 70, 2024u);
    const auto report = isim_report(make_result(mixed_labels()), batch, records_options(false));
    EXPECT_TRUE(std::isnan(report.isim_silhouette));
    EXPECT_TRUE(std::isnan(report.davies_bouldin_medoid));
    EXPECT_TRUE(std::isnan(report.dunn_medoid_separation_medoid_spread));
    EXPECT_TRUE(report.coverage_at.empty());
    EXPECT_TRUE(report.noise_coverage_at.empty());
    EXPECT_EQ(report.coverage_thresholds.size(), 3u);
    EXPECT_FALSE(report.requested.centroid_indices);
    EXPECT_TRUE(report.requested.per_cluster_records);
    for (const auto& record : report.records) {
        EXPECT_TRUE(std::isnan(record.isim_silhouette));
        EXPECT_EQ(record.nearest_cluster, NO_NEAREST_CLUSTER);
        EXPECT_TRUE(std::isnan(record.nearest_cluster_similarity));
    }
}

TEST(ISimReportValidationTest, MetricIsCheckedFirst) {
    const auto batch = isim_test::make_random_batch(3, 16, 1u);
    ISimReportOptions options;
    options.metric = "dice";
    options.coverage_thresholds = {std::numeric_limits<double>::quiet_NaN()};
    // Labels longer than the batch would also fail; the metric is named first.
    try {
        isim_report(make_result({0, 0, 0, 0}), batch, options);
        FAIL() << "expected invalid_argument";
    } catch (const std::invalid_argument& error) {
        EXPECT_NE(std::string(error.what()).find("metric"), std::string::npos);
    }
}

TEST(ISimReportValidationTest, NaNThresholdPrecedesThePartitionChecks) {
    const auto batch = isim_test::make_random_batch(3, 16, 1u);
    ISimReportOptions options;
    options.coverage_thresholds = {0.3, std::numeric_limits<double>::quiet_NaN()};
    try {
        isim_report(make_result({0, 0, 0, 0}), batch, options);
        FAIL() << "expected invalid_argument";
    } catch (const std::invalid_argument& error) {
        EXPECT_EQ(std::string(error.what()), "isim_report: coverage threshold 1 must not be NaN");
    }
    // Refused at K == 0 too: the spec's native order puts it before any partition check.
    EXPECT_THROW(isim_report(make_result({-1, -1, -1}), batch, options), std::invalid_argument);
}

TEST(ISimReportValidationTest, LabelCountBeyondTheBatchIsOutOfRange) {
    const auto batch = isim_test::make_random_batch(3, 16, 1u);
    try {
        isim_report(make_result({0, 0, 0, 0}), batch);
        FAIL() << "expected out_of_range";
    } catch (const std::out_of_range& error) {
        EXPECT_EQ(std::string(error.what()),
                  "isim_report: label count 4 exceeds the fingerprint batch sample count 3");
    }
}

TEST(ISimReportValidationTest, MalformedPartitionIsRefused) {
    const auto batch = isim_test::make_random_batch(2, 16, 1u);
    const ClusteringResult duplicated({0, 0}, Clusters{{0, 1}, {1}});
    EXPECT_THROW(isim_report(duplicated, batch), std::invalid_argument);
}

TEST(ISimReportValidationTest, NegativeThresholdIsAcceptedNatively) {
    const auto batch = isim_test::make_random_batch(4, 16, 1u);
    ISimReportOptions options;
    options.coverage_thresholds = {-0.1};
    EXPECT_NO_THROW(isim_report(make_result({0, 0, 1, 1}), batch, options));
}
