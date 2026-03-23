/**
 * @file test_bio_metrics.cpp
 * @brief Tests for SuperposeMetric using oespruce OESuperpose.
 */

#include <cmath>
#include <gtest/gtest.h>
#include "oecluster/oecluster.h"
#include "oecluster/metrics/SuperposeMetric.h"
#include <oechem.h>
#include <oebio.h>

using namespace OECluster;

static std::string asset_path(const std::string& name) {
    return std::string(TEST_ASSETS_DIR) + "/" + name;
}

class SuperposeMetricDUTest : public ::testing::Test {
protected:
    void SetUp() override {
        std::vector<std::string> files = {
            asset_path("spruce_5FQD_1_5FQD_1-ALIGNED_BC__DU__LVY_B-1438.oedu"),
            asset_path("spruce_8G66_1_8G66_1-ALIGNED_BC__DU__YOT_B-502.oedu"),
        };
        for (const auto& f : files) {
            auto du = std::make_shared<OEBio::OEDesignUnit>();
            if (!OEBio::OEReadDesignUnit(f, *du)) {
                GTEST_SKIP() << "Cannot read test asset: " << f;
            }
            dus_.push_back(du);
        }
    }

    std::vector<std::shared_ptr<OEBio::OEDesignUnit>> dus_;
};

// -- Method tests --

TEST_F(SuperposeMetricDUTest, GlobalCarbonAlpha_RMSD) {
    SuperposeOptions opts;
    opts.method = SuperposeMethod::GlobalCarbonAlpha;
    SuperposeMetric metric(dus_, opts);
    EXPECT_EQ(metric.Size(), 2u);
    EXPECT_EQ(metric.Name(), "superpose:global_carbon_alpha");
    double d = metric.Distance(0, 1);
    EXPECT_GT(d, 0.0);
    EXPECT_TRUE(std::isfinite(d));
}

TEST_F(SuperposeMetricDUTest, Global_RMSD) {
    SuperposeOptions opts;
    opts.method = SuperposeMethod::Global;
    SuperposeMetric metric(dus_, opts);
    EXPECT_EQ(metric.Name(), "superpose:global");
    double d = metric.Distance(0, 1);
    EXPECT_GT(d, 0.0);
}

TEST_F(SuperposeMetricDUTest, DDM_RMSD) {
    SuperposeOptions opts;
    opts.method = SuperposeMethod::DDM;
    SuperposeMetric metric(dus_, opts);
    EXPECT_EQ(metric.Name(), "superpose:ddm");
    double d = metric.Distance(0, 1);
    EXPECT_GT(d, 0.0);
}

TEST_F(SuperposeMetricDUTest, Weighted_RMSD) {
    SuperposeOptions opts;
    opts.method = SuperposeMethod::Weighted;
    SuperposeMetric metric(dus_, opts);
    EXPECT_EQ(metric.Name(), "superpose:weighted");
    double d = metric.Distance(0, 1);
    EXPECT_GT(d, 0.0);
}

TEST_F(SuperposeMetricDUTest, SSE_Tanimoto) {
    SuperposeOptions opts;
    opts.method = SuperposeMethod::SSE;
    SuperposeMetric metric(dus_, opts);
    EXPECT_EQ(metric.Name(), "superpose:sse");
    double d = metric.Distance(0, 1);
    // Distance = 1.0 - tanimoto, should be in [0, 1]
    EXPECT_GE(d, 0.0);
    EXPECT_LE(d, 1.0);
}

TEST_F(SuperposeMetricDUTest, SiteHopper_PatchScore) {
    SuperposeOptions opts;
    opts.method = SuperposeMethod::SiteHopper;
    SuperposeMetric metric(dus_, opts);
    EXPECT_EQ(metric.Name(), "superpose:sitehopper");
    double d = metric.Distance(0, 1);
    // Distance = 4.0 - patch_score, should be in [0, 4]
    EXPECT_GE(d, 0.0);
    EXPECT_LE(d, 4.0);
}

// -- Distance / Similarity mode tests --

TEST_F(SuperposeMetricDUTest, SSE_Similarity) {
    SuperposeOptions opts;
    opts.method = SuperposeMethod::SSE;
    opts.similarity = true;
    SuperposeMetric metric(dus_, opts);
    double sim = metric.Distance(0, 1);
    EXPECT_GE(sim, 0.0);
    EXPECT_LE(sim, 1.0);
}

TEST_F(SuperposeMetricDUTest, SiteHopper_Similarity) {
    SuperposeOptions opts;
    opts.method = SuperposeMethod::SiteHopper;
    opts.similarity = true;
    SuperposeMetric metric(dus_, opts);
    double sim = metric.Distance(0, 1);
    EXPECT_GE(sim, 0.0);
    EXPECT_LE(sim, 4.0);
}

TEST_F(SuperposeMetricDUTest, RMSD_SimilarityIsNoOp) {
    SuperposeOptions opts;
    opts.method = SuperposeMethod::GlobalCarbonAlpha;

    opts.similarity = false;
    SuperposeMetric dist_metric(dus_, opts);
    double d = dist_metric.Distance(0, 1);

    opts.similarity = true;
    SuperposeMetric sim_metric(dus_, opts);
    double s = sim_metric.Distance(0, 1);

    EXPECT_DOUBLE_EQ(d, s);
}

// -- Score type validation --

TEST_F(SuperposeMetricDUTest, IncompatibleScoreTypeThrows) {
    SuperposeOptions opts;
    opts.method = SuperposeMethod::SSE;
    opts.score_type = SuperposeScoreType::RMSD;
    EXPECT_THROW(SuperposeMetric(dus_, opts), MetricError);
}

TEST_F(SuperposeMetricDUTest, AutoScoreTypeResolvesCorrectly) {
    SuperposeOptions opts;
    opts.method = SuperposeMethod::GlobalCarbonAlpha;
    opts.score_type = SuperposeScoreType::Auto;
    SuperposeMetric metric(dus_, opts);
    double d = metric.Distance(0, 1);
    // Auto resolves to RMSD for GlobalCarbonAlpha
    EXPECT_GT(d, 0.0);
}

// -- Predicate tests --

TEST_F(SuperposeMetricDUTest, ValidPredicateConstructionSucceeds) {
    // Verify that a valid predicate is parsed and construction succeeds.
    // The actual superposition may fail depending on atom coverage, so
    // we only test construction here.
    SuperposeOptions opts;
    opts.method = SuperposeMethod::GlobalCarbonAlpha;
    opts.predicate = "backbone";
    EXPECT_NO_THROW({ SuperposeMetric m(dus_, opts); });
}

TEST_F(SuperposeMetricDUTest, InvalidPredicateThrows) {
    SuperposeOptions opts;
    opts.predicate = "!!!INVALID_PREDICATE!!!";
    EXPECT_THROW(SuperposeMetric(dus_, opts), MetricError);
}

// -- Clone and pdist integration --

TEST_F(SuperposeMetricDUTest, CloneProducesValidResults) {
    SuperposeMetric metric(dus_);
    auto clone = metric.Clone();
    EXPECT_EQ(clone->Size(), 2u);
    double d = clone->Distance(0, 1);
    EXPECT_TRUE(std::isfinite(d));
}

TEST_F(SuperposeMetricDUTest, IntegrationWithPDist) {
    SuperposeMetric metric(dus_);
    DenseStorage storage(2);
    pdist(metric, storage);
    EXPECT_EQ(storage.NumPairs(), 1u);
    double d = storage.Get(0, 1);
    EXPECT_TRUE(std::isfinite(d));
}

// -- Error cases --

TEST_F(SuperposeMetricDUTest, EmptyStructureListThrows) {
    std::vector<std::shared_ptr<OEBio::OEDesignUnit>> empty_dus;
    EXPECT_THROW({ SuperposeMetric m(empty_dus); }, MetricError);
}
