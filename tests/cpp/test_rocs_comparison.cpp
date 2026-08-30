#include <gtest/gtest.h>
#include "oecluster/oecluster.h"
#include "oecluster/comparisons/ROCSComparison.h"
#include <oechem.h>
#include <oeshape.h>
#include <string>

using namespace OECluster;

class ROCSComparisonTest : public ::testing::Test {
protected:
    void SetUp() override {
        auto make_mol = [](const char* smi) -> std::shared_ptr<OEChem::OEMol> {
            auto mol = std::make_shared<OEChem::OEMol>();
            OEChem::OESmilesToMol(*mol, smi);
            OEChem::OEAddExplicitHydrogens(*mol);
            OEChem::OEGenerate2DCoordinates(*mol);
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

TEST_F(ROCSComparisonTest, FactsCoverEveryScoreTypeAndDirection) {
    // Self-comparison values on valid 3D conformers: shape Tanimoto is 1.0 and
    // color Tanimoto is 0.0, because nothing prepares color atoms. Exactly two
    // of the eight configurations therefore vanish on the diagonal -- the Shape
    // distance and, as an artifact of that same broken color term, the Color
    // similarity. is_distance follows the flag alone; zero_self does not.
    struct Row {
        ROCSScoreType score_type;
        bool similarity;
        Capability is_distance;
        Capability zero_self;
    };
    const std::vector<Row> rows{
        {ROCSScoreType::ComboNorm, false, Capability::Yes, Capability::No},
        {ROCSScoreType::Combo, false, Capability::Yes, Capability::No},
        {ROCSScoreType::Shape, false, Capability::Yes, Capability::Yes},
        {ROCSScoreType::Color, false, Capability::Yes, Capability::No},
        {ROCSScoreType::ComboNorm, true, Capability::No, Capability::No},
        {ROCSScoreType::Combo, true, Capability::No, Capability::No},
        {ROCSScoreType::Shape, true, Capability::No, Capability::No},
        {ROCSScoreType::Color, true, Capability::No, Capability::Yes},
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
