#include <gtest/gtest.h>
#include "oecluster/oecluster.h"
#include "oecluster/comparisons/ROCSComparison.h"
#include <oechem.h>
#include <oeomega2.h>
#include <oeshape.h>
#include <string>

using namespace OECluster;

class ROCSComparisonTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Omega, not OEGenerate2DCoordinates: a planar molecule has zero shape
        // volume, so every ROCS Tanimoto degenerates to 0/0 and the whole suite
        // measures nothing.
        OEConfGen::OEOmega omega;
        omega.SetMaxConfs(1);
        omega.SetStrictStereo(false);
        auto make_mol = [&omega](const char* smi) -> std::shared_ptr<OEChem::OEMol> {
            auto mol = std::make_shared<OEChem::OEMol>();
            OEChem::OESmilesToMol(*mol, smi);
            // EXPECT_TRUE rather than ASSERT_TRUE: ASSERT_* expands to a bare
            // return, which cannot compile in a value-returning lambda. A failed
            // build leaves the molecule empty, which surfaces loudly downstream.
            EXPECT_TRUE(omega(*mol)) << smi;
            return mol;
        };
        mols_.push_back(make_mol("c1ccccc1"));      // benzene
        mols_.push_back(make_mol("c1ccc(O)cc1"));    // phenol
        mols_.push_back(make_mol("CCCCCCCC"));        // octane
    }

    std::vector<std::shared_ptr<OEChem::OEMol>> mols_;
};

TEST_F(ROCSComparisonTest, ConstructAndSize) {
    ROCSComparison comparison(mols_);
    EXPECT_EQ(comparison.Size(), 3);
    EXPECT_EQ(comparison.ComparisonName(), "rocs");
}

TEST_F(ROCSComparisonTest, DistanceRange) {
    ROCSComparison comparison(mols_);
    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = i + 1; j < 3; ++j) {
            double d = comparison.Compare(i, j);
            // Default is ComboNorm: distance ∈ [0,1]
            EXPECT_GE(d, 0.0);
            EXPECT_LE(d, 1.0);
        }
    }
}

TEST_F(ROCSComparisonTest, ShapeOnlyDistanceRange) {
    ROCSOptions opts;
    opts.score_type = ROCSScoreType::Shape;
    ROCSComparison comparison(mols_, opts);
    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = i + 1; j < 3; ++j) {
            double d = comparison.Compare(i, j);
            // Shape Tanimoto is in [0,1], so distance is in [0, 1]
            EXPECT_GE(d, 0.0);
            EXPECT_LE(d, 1.0);
        }
    }
}

TEST_F(ROCSComparisonTest, CloneCreatesIndependentScorer) {
    ROCSComparison comparison(mols_);
    auto clone = comparison.Clone();
    EXPECT_EQ(clone->Size(), 3);
    double d = clone->Compare(0, 1);
    // Default is ComboNorm: distance ∈ [0,1]
    EXPECT_GE(d, 0.0);
    EXPECT_LE(d, 1.0);
}

TEST_F(ROCSComparisonTest, IntegrationWithPDist) {
    ROCSComparison comparison(mols_);
    DenseStorage storage(3);
    pdist(comparison, storage);
    EXPECT_EQ(storage.NumPairs(), 3);
}

TEST_F(ROCSComparisonTest, ComboNormScoreType) {
    ROCSOptions opts;
    opts.score_type = ROCSScoreType::ComboNorm;
    ROCSComparison comparison(mols_, opts);
    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = i + 1; j < 3; ++j) {
            double d = comparison.Compare(i, j);
            // ComboNorm: distance = 1.0 - combo/2.0 ∈ [0,1]
            EXPECT_GE(d, 0.0);
            EXPECT_LE(d, 1.0);
        }
    }
}

TEST_F(ROCSComparisonTest, ComboScoreType) {
    ROCSOptions opts;
    opts.score_type = ROCSScoreType::Combo;
    ROCSComparison comparison(mols_, opts);
    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = i + 1; j < 3; ++j) {
            double d = comparison.Compare(i, j);
            // Combo: distance = 2.0 - combo ∈ [0,2]
            EXPECT_GE(d, 0.0);
            EXPECT_LE(d, 2.0);
        }
    }
}

TEST_F(ROCSComparisonTest, ColorScoreType) {
    ROCSOptions opts;
    opts.score_type = ROCSScoreType::Color;
    ROCSComparison comparison(mols_, opts);
    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = i + 1; j < 3; ++j) {
            double d = comparison.Compare(i, j);
            // Color: distance = 1.0 - color ∈ [0,1]
            EXPECT_GE(d, 0.0);
            EXPECT_LE(d, 1.0);
        }
    }
}

TEST_F(ROCSComparisonTest, ColorForceFieldConfiguration) {
    ROCSOptions opts;
    opts.score_type = ROCSScoreType::Color;
    opts.color_ff_type = 2;  // ExplicitMillsDean
    ROCSComparison comparison(mols_, opts);
    EXPECT_EQ(comparison.Size(), 3);
    // Just verify it constructs and runs without error
    double d = comparison.Compare(0, 1);
    EXPECT_GE(d, 0.0);
    EXPECT_LE(d, 1.0);
}

// Index 1 is the phenol: its hydroxyl carries donor and acceptor color atoms,
// where benzene has almost none, so the color term is genuinely exercised.
// The tolerance is 1e-3 rather than something tighter because BestOverlay
// optimizes from several starting poses rather than being handed the identity,
// so the diagonal is a converged result and not an exact one by construction.

TEST_F(ROCSComparisonTest, ColorSelfSimilarityIsOne) {
    ROCSOptions opts;
    opts.score_type = ROCSScoreType::Color;
    opts.similarity = true;
    ROCSComparison comparison(mols_, opts);
    EXPECT_NEAR(comparison.Compare(1, 1), 1.0, 1e-3);
}

TEST_F(ROCSComparisonTest, ComboSelfSimilarityIsTwo) {
    ROCSOptions opts;
    opts.score_type = ROCSScoreType::Combo;
    opts.similarity = true;
    ROCSComparison comparison(mols_, opts);
    EXPECT_NEAR(comparison.Compare(1, 1), 2.0, 1e-3);
}

TEST_F(ROCSComparisonTest, ComboSelfScoreIsZeroDistance) {
    ROCSOptions opts;
    opts.score_type = ROCSScoreType::Combo;
    ROCSComparison comparison(mols_, opts);
    EXPECT_NEAR(comparison.Compare(1, 1), 0.0, 1e-3);
}

TEST_F(ROCSComparisonTest, InputMoleculesAreNotMutated) {
    const unsigned int before = mols_[1]->NumAtoms();
    ROCSComparison comparison(mols_, ROCSOptions());
    comparison.Compare(0, 1);
    EXPECT_EQ(mols_[1]->NumAtoms(), before);
}

TEST_F(ROCSComparisonTest, FactsCoverEveryScoreTypeAndDirection) {
    // Self-comparison values on valid 3D conformers: shape and color Tanimoto
    // both reach 1.0 now that color atoms are prepared, so combo reaches 2.0.
    // Every distance form subtracts its own saturation value and vanishes on the
    // diagonal; every similarity form returns that value instead. The Color
    // similarity cell used to read Yes, as an artifact of a color term that was
    // identically zero, and correctly reads No now that it is not.
    struct Row {
        ROCSScoreType score_type;
        bool similarity;
        Capability is_distance;
        Capability zero_self;
    };
    const std::vector<Row> rows{
        {ROCSScoreType::ComboNorm, false, Capability::Yes, Capability::Yes},
        {ROCSScoreType::Combo, false, Capability::Yes, Capability::Yes},
        {ROCSScoreType::Shape, false, Capability::Yes, Capability::Yes},
        {ROCSScoreType::Color, false, Capability::Yes, Capability::Yes},
        {ROCSScoreType::ComboNorm, true, Capability::No, Capability::No},
        {ROCSScoreType::Combo, true, Capability::No, Capability::No},
        {ROCSScoreType::Shape, true, Capability::No, Capability::No},
        {ROCSScoreType::Color, true, Capability::No, Capability::No},
    };

    for (const Row& row : rows) {
        ROCSOptions opts;
        opts.score_type = row.score_type;
        opts.similarity = row.similarity;
        ROCSComparison comparison(mols_, opts);
        const GateFacts facts = comparison.Facts();
        const std::string label =
            std::to_string(static_cast<int>(row.score_type)) +
            (row.similarity ? " similarity" : " distance");
        EXPECT_EQ(facts.is_distance, row.is_distance) << label;
        EXPECT_EQ(facts.zero_self, row.zero_self) << label;
        EXPECT_EQ(facts.triangle, Capability::Unknown) << label;
        EXPECT_EQ(facts.data_integrity, DataIntegrity::Complete) << label;
    }
}
