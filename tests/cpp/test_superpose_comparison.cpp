/**
 * @file test_superpose_comparison.cpp
 * @brief Tests for SuperposeComparison using oespruce OESuperpose.
 */

#include <cmath>
#include <string>
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

TEST_F(SuperposeComparisonDUTest, CompareRefusesAnIndexPastTheEnd) {
    // Bounded on the design-unit vector here, and on the molecule vector when
    // the comparison was built from molecules; only one of the two is ever
    // populated.
    SuperposeComparison comparison(dus_);
    try {
        comparison.Compare(0, 1000000);
        FAIL() << "expected ComparisonError";
    } catch (const ComparisonError& error) {
        const std::string message(error.what());
        EXPECT_NE(message.find("1000000"), std::string::npos) << message;
        EXPECT_NE(message.find("2 items"), std::string::npos) << message;
    }
    EXPECT_THROW(comparison.Compare(1000000, 0), ComparisonError);
    EXPECT_THROW(comparison.Compare(dus_.size(), dus_.size()), ComparisonError);
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

// Compile-time tripwire for the hand-written method list below. An exhaustive
// switch with no default case makes the compiler diagnose a SuperposeMethod
// that has gained an enumerator, which an assertion on the list's own size
// cannot do.
static Capability expected_similarity_zero_self(SuperposeMethod method) {
    switch (method) {
        // These four resolve to RMSD, where the similarity flag is a
        // documented no-op, so both directions are Yes.
        case SuperposeMethod::GlobalCarbonAlpha:
        case SuperposeMethod::Global:
        case SuperposeMethod::DDM:
        case SuperposeMethod::Weighted:
            return Capability::Yes;
        case SuperposeMethod::SSE:
        case SuperposeMethod::SiteHopper:
            return Capability::No;
    }
    return Capability::Unknown;  // unreachable; silences a missing-return warning
}

// Second tripwire, same reasoning as the one above: an exhaustive switch with
// no default case makes a new SuperposeMethod a compile error rather than a
// silently wrong orientation.
static Capability expected_similarity_is_distance(SuperposeMethod method) {
    switch (method) {
        // RMSD ignores the similarity flag, so these stay distances.
        case SuperposeMethod::GlobalCarbonAlpha:
        case SuperposeMethod::Global:
        case SuperposeMethod::DDM:
        case SuperposeMethod::Weighted:
            return Capability::Yes;
        case SuperposeMethod::SSE:
        case SuperposeMethod::SiteHopper:
            return Capability::No;
    }
    return Capability::Unknown;  // unreachable; silences a missing-return warning
}

TEST_F(SuperposeComparisonDUTest, FactsFollowTheResolvedScoreType) {
    const std::vector<SuperposeMethod> methods{
        SuperposeMethod::GlobalCarbonAlpha,
        SuperposeMethod::Global,
        SuperposeMethod::DDM,
        SuperposeMethod::Weighted,
        SuperposeMethod::SSE,
        SuperposeMethod::SiteHopper,
    };

    for (SuperposeMethod method : methods) {
        SuperposeOptions distance_opts;
        distance_opts.method = method;
        distance_opts.similarity = false;
        SuperposeComparison distance_comparison(dus_, distance_opts);
        const GateFacts distance_facts = distance_comparison.Facts();
        EXPECT_EQ(distance_facts.is_distance, Capability::Yes)
            << static_cast<int>(method);
        EXPECT_EQ(distance_facts.zero_self, Capability::Yes)
            << static_cast<int>(method);
        EXPECT_EQ(distance_facts.triangle, Capability::Unknown)
            << static_cast<int>(method);
        EXPECT_EQ(distance_facts.data_integrity, DataIntegrity::Complete)
            << static_cast<int>(method);

        SuperposeOptions similarity_opts = distance_opts;
        similarity_opts.similarity = true;
        SuperposeComparison similarity_comparison(dus_, similarity_opts);
        const GateFacts similarity_facts = similarity_comparison.Facts();
        EXPECT_EQ(similarity_facts.is_distance,
                  expected_similarity_is_distance(method))
            << static_cast<int>(method);
        EXPECT_EQ(similarity_facts.zero_self,
                  expected_similarity_zero_self(method))
            << static_cast<int>(method);
        EXPECT_EQ(similarity_facts.data_integrity, DataIntegrity::Complete)
            << static_cast<int>(method);
    }
}

TEST_F(SuperposeComparisonDUTest, AsymmetricSelectionsDoNotBreakTheDiagonal) {
    // Facts() stamps zero_self = Yes without reading the selection strings.
    // That is only sound because superposing a structure onto itself yields
    // the identity transform whatever atoms each side selects: a nested pair
    // still scores 0.0, and a pair too disjoint to align raises rather than
    // returning a nonzero diagonal.
    SuperposeOptions nested_opts;
    nested_opts.method = SuperposeMethod::GlobalCarbonAlpha;
    nested_opts.ref_predicate = "protein";
    nested_opts.fit_predicate = "backbone";
    SuperposeComparison nested(dus_, nested_opts);
    EXPECT_EQ(nested.Facts().is_distance, Capability::Yes);
    EXPECT_EQ(nested.Facts().zero_self, Capability::Yes);
    EXPECT_DOUBLE_EQ(nested.Compare(0, 0), 0.0);

    SuperposeOptions disjoint_opts;
    disjoint_opts.method = SuperposeMethod::GlobalCarbonAlpha;
    disjoint_opts.ref_predicate = "protein";
    disjoint_opts.fit_predicate = "ligand";
    SuperposeComparison disjoint(dus_, disjoint_opts);
    EXPECT_THROW(disjoint.Compare(0, 0), ComparisonError);
}
