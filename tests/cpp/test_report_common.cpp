#include <gtest/gtest.h>

#include <stdexcept>
#include <string>
#include <vector>

#include "oecluster/clustering/ClusterTypes.h"
#include "../../src/clustering/ReportCommon.h"

using namespace OECluster;

namespace {
ClusteringResult make_result(std::vector<ClusterLabel> labels) {
    Clusters members = labels_to_clusters(labels);
    return ClusteringResult(std::move(labels), std::move(members));
}
}  // namespace

TEST(ReportCommonTest, PresetThresholdsMatchTheReportTable) {
    EXPECT_EQ(detail::preset_coverage_thresholds(ClusterThreshold::Default),
              (std::vector<double>{0.25, 0.35, 0.45}));
    EXPECT_EQ(detail::preset_coverage_thresholds(ClusterThreshold::Tight),
              (std::vector<double>{0.20, 0.30, 0.40}));
    EXPECT_EQ(detail::preset_coverage_thresholds(ClusterThreshold::Diversity),
              (std::vector<double>{0.40, 0.50, 0.60}));
}

TEST(ReportCommonTest, ProfileCountsNoiseAndSingletons) {
    const auto profile = detail::report_profile(make_result({0, 0, 1, -1}), true);
    EXPECT_EQ(profile.num_samples, 4u);
    EXPECT_EQ(profile.num_clusters, 2u);
    EXPECT_EQ(profile.num_noise, 1u);
    EXPECT_EQ(profile.num_singletons, 1u);
    EXPECT_DOUBLE_EQ(profile.noise_fraction, 0.25);
    EXPECT_DOUBLE_EQ(profile.singleton_fraction, 2.0 / 3.0);
    const auto without = detail::report_profile(make_result({0, 0, 1, -1}), false);
    EXPECT_DOUBLE_EQ(without.singleton_fraction, 0.5);
}

TEST(ReportCommonTest, OwnerArrayMapsSamplesToOrdinals) {
    const auto owner = detail::validate_report_partition(
        make_result({1, 0, -1, 1}), 4, "isim_report", "fingerprint batch");
    EXPECT_EQ(owner, (std::vector<size_t>{1, 0, detail::NO_OWNER, 1}));
}

TEST(ReportCommonTest, MessagesCarryTheCallerAndNoun) {
    try {
        detail::validate_report_partition(make_result({0, 0, 0}), 2, "isim_report",
                                          "fingerprint batch");
        FAIL() << "expected out_of_range";
    } catch (const std::out_of_range& error) {
        EXPECT_EQ(std::string(error.what()),
                  "isim_report: label count 3 exceeds the fingerprint batch sample count 2");
    }
    const ClusteringResult duplicated({0, 0}, Clusters{{0, 1}, {1}});
    try {
        detail::validate_report_partition(duplicated, 2, "isim_report", "fingerprint batch");
        FAIL() << "expected invalid_argument";
    } catch (const std::invalid_argument& error) {
        EXPECT_EQ(std::string(error.what()),
                  "isim_report: sample 1 appears in clusters 0 and 1");
    }
}
