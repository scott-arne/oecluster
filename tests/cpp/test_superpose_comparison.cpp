/**
 * @file test_superpose_comparison.cpp
 * @brief Tests for SuperposeComparison using oespruce OESuperpose.
 */

#include <cmath>
#include <gtest/gtest.h>
#include "oecluster/oecluster.h"
#include "oecluster/comparisons/SuperposeComparison.h"
#include <oechem.h>
#include <oebio.h>

using namespace OECluster;

static std::string asset_path(const std::string& name) {
    return std::string(TEST_ASSETS_DIR) + "/" + name;
}

class SuperposeComparisonDUTest : public ::testing::Test {
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

TEST_F(SuperposeComparisonDUTest, GlobalCarbonAlpha_RMSD) {
    SuperposeOptions opts;
    opts.method = SuperposeMethod::GlobalCarbonAlpha;
    SuperposeComparison comparison(dus_, opts);
    EXPECT_EQ(comparison.Size(), 2u);
    EXPECT_EQ(comparison.ComparisonName(), "superpose:global_carbon_alpha");
    double d = comparison.Compare(0, 1);
    EXPECT_GT(d, 0.0);
    EXPECT_TRUE(std::isfinite(d));
}

TEST_F(SuperposeComparisonDUTest, Global_RMSD) {
    SuperposeOptions opts;
    opts.method = SuperposeMethod::Global;
    SuperposeComparison comparison(dus_, opts);
    EXPECT_EQ(comparison.ComparisonName(), "superpose:global");
    double d = comparison.Compare(0, 1);
    EXPECT_GT(d, 0.0);
}

TEST_F(SuperposeComparisonDUTest, DDM_RMSD) {
    SuperposeOptions opts;
    opts.method = SuperposeMethod::DDM;
    SuperposeComparison comparison(dus_, opts);
    EXPECT_EQ(comparison.ComparisonName(), "superpose:ddm");
    double d = comparison.Compare(0, 1);
    EXPECT_GT(d, 0.0);
}

TEST_F(SuperposeComparisonDUTest, Weighted_RMSD) {
    SuperposeOptions opts;
    opts.method = SuperposeMethod::Weighted;
    SuperposeComparison comparison(dus_, opts);
    EXPECT_EQ(comparison.ComparisonName(), "superpose:weighted");
    double d = comparison.Compare(0, 1);
    EXPECT_GT(d, 0.0);
}

TEST_F(SuperposeComparisonDUTest, SSE_Tanimoto) {
    SuperposeOptions opts;
    opts.method = SuperposeMethod::SSE;
    SuperposeComparison comparison(dus_, opts);
    EXPECT_EQ(comparison.ComparisonName(), "superpose:sse");
    double d = comparison.Compare(0, 1);
    // Distance = 1.0 - tanimoto, should be in [0, 1]
    EXPECT_GE(d, 0.0);
    EXPECT_LE(d, 1.0);
}

TEST_F(SuperposeComparisonDUTest, SiteHopper_PatchScore) {
    SuperposeOptions opts;
    opts.method = SuperposeMethod::SiteHopper;
    SuperposeComparison comparison(dus_, opts);
    EXPECT_EQ(comparison.ComparisonName(), "superpose:sitehopper");
    double d = comparison.Compare(0, 1);
    // Distance = 4.0 - patch_score, should be in [0, 4]
    EXPECT_GE(d, 0.0);
    EXPECT_LE(d, 4.0);
}

// -- Distance / Similarity mode tests --

TEST_F(SuperposeComparisonDUTest, SSE_Similarity) {
    SuperposeOptions opts;
    opts.method = SuperposeMethod::SSE;
    opts.similarity = true;
    SuperposeComparison comparison(dus_, opts);
    double sim = comparison.Compare(0, 1);
    EXPECT_GE(sim, 0.0);
    EXPECT_LE(sim, 1.0);
}

TEST_F(SuperposeComparisonDUTest, SiteHopper_Similarity) {
    SuperposeOptions opts;
    opts.method = SuperposeMethod::SiteHopper;
    opts.similarity = true;
    SuperposeComparison comparison(dus_, opts);
    double sim = comparison.Compare(0, 1);
    EXPECT_GE(sim, 0.0);
    EXPECT_LE(sim, 4.0);
}

TEST_F(SuperposeComparisonDUTest, RMSD_SimilarityIsNoOp) {
    SuperposeOptions opts;
    opts.method = SuperposeMethod::GlobalCarbonAlpha;

    opts.similarity = false;
    SuperposeComparison distance_comparison(dus_, opts);
    double d = distance_comparison.Compare(0, 1);

    opts.similarity = true;
    SuperposeComparison similarity_comparison(dus_, opts);
    double s = similarity_comparison.Compare(0, 1);

    EXPECT_DOUBLE_EQ(d, s);
}

// -- Score type validation --

TEST_F(SuperposeComparisonDUTest, IncompatibleScoreTypeThrows) {
    SuperposeOptions opts;
    opts.method = SuperposeMethod::SSE;
    opts.score_type = SuperposeScoreType::RMSD;
    EXPECT_THROW(SuperposeComparison(dus_, opts), ComparisonError);
}

TEST_F(SuperposeComparisonDUTest, AutoScoreTypeResolvesCorrectly) {
    SuperposeOptions opts;
    opts.method = SuperposeMethod::GlobalCarbonAlpha;
    opts.score_type = SuperposeScoreType::Auto;
    SuperposeComparison comparison(dus_, opts);
    double d = comparison.Compare(0, 1);
    // Auto resolves to RMSD for GlobalCarbonAlpha
    EXPECT_GT(d, 0.0);
}

// -- Predicate tests --

TEST_F(SuperposeComparisonDUTest, ValidPredicateConstructionSucceeds) {
    // Verify that a valid predicate is parsed and construction succeeds.
    // The actual superposition may fail depending on atom coverage, so
    // we only test construction here.
    SuperposeOptions opts;
    opts.method = SuperposeMethod::GlobalCarbonAlpha;
    opts.predicate = "backbone";
    EXPECT_NO_THROW({ SuperposeComparison m(dus_, opts); });
}

TEST_F(SuperposeComparisonDUTest, InvalidPredicateThrows) {
    SuperposeOptions opts;
    opts.predicate = "!!!INVALID_PREDICATE!!!";
    EXPECT_THROW(SuperposeComparison(dus_, opts), ComparisonError);
}

// -- Clone and pdist integration --

TEST_F(SuperposeComparisonDUTest, CloneProducesValidResults) {
    SuperposeComparison comparison(dus_);
    auto clone = comparison.Clone();
    EXPECT_EQ(clone->Size(), 2u);
    double d = clone->Compare(0, 1);
    EXPECT_TRUE(std::isfinite(d));
}

TEST_F(SuperposeComparisonDUTest, IntegrationWithPDist) {
    SuperposeComparison comparison(dus_);
    DenseStorage storage(2);
    pdist(comparison, storage);
    EXPECT_EQ(storage.NumPairs(), 1u);
    double d = storage.Get(0, 1);
    EXPECT_TRUE(std::isfinite(d));
}

// -- Error cases --

TEST_F(SuperposeComparisonDUTest, EmptyStructureListThrows) {
    std::vector<std::shared_ptr<OEBio::OEDesignUnit>> empty_dus;
    EXPECT_THROW({ SuperposeComparison m(empty_dus); }, ComparisonError);
}
