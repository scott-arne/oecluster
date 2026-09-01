#include <gtest/gtest.h>
#include "oecluster/oecluster.h"
#include "oecluster/comparisons/ROCSComparison.h"
#include <oechem.h>
#include <oeomega2.h>
#include <oeshape.h>
#include <cmath>
#include <limits>
#include <string>

using namespace OECluster;

class ROCSComparisonTest : public ::testing::Test {
protected:
    /// A single Omega conformer for ``smi``. Shared by SetUp and by the tests
    /// that need molecules outside the standard fixture, so that every molecule
    /// in this file is built under one Omega configuration.
    ///
    /// Omega, not OEGenerate2DCoordinates: the overlay finds no coordinates it
    /// can use on planar input, warns, and returns a Tanimoto of exactly 0.0 for
    /// every pair, so the whole suite would measure saturation values rather
    /// than shape.
    static std::shared_ptr<OEChem::OEMol> MakeConformer(const char* smi) {
        OEConfGen::OEOmega omega;
        omega.SetMaxConfs(1);
        omega.SetStrictStereo(false);
        auto mol = std::make_shared<OEChem::OEMol>();
        OEChem::OESmilesToMol(*mol, smi);
        // EXPECT_TRUE rather than ASSERT_TRUE: ASSERT_* expands to a bare
        // return, which cannot compile in a value-returning function. A failed
        // build leaves the molecule empty, which surfaces loudly downstream.
        EXPECT_TRUE(omega(*mol)) << smi;
        return mol;
    }

    void SetUp() override {
        mols_.push_back(MakeConformer("c1ccccc1"));     // benzene
        mols_.push_back(MakeConformer("c1ccc(O)cc1"));  // phenol
        mols_.push_back(MakeConformer("CCCCCCCC"));     // octane
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

TEST_F(ROCSComparisonTest, ZeroSelfIsNoWhenAMoleculeHasNoColorAtoms) {
    // ImplicitMillsDean gives methane no color atom at all, so its color
    // self-Tanimoto is 0.0 however well the prep runs, and ComboNorm's diagonal
    // sits at 0.5. Asserting the distance as well as the stamp keeps this test
    // pinned to the reason rather than to the verdict.
    auto methane = MakeConformer("C");
    ROCSComparison comparison({methane}, ROCSOptions());
    EXPECT_NEAR(comparison.Compare(0, 0), 0.5, 1e-3);
    EXPECT_EQ(comparison.Facts().zero_self, Capability::No);
}

TEST_F(ROCSComparisonTest, ZeroSelfFollowsTheMeasurementNotTheSimilarityFlag) {
    // The stamp is a measurement, and this is the case that proves it: methane
    // has no color atoms, so its color *similarity* to itself is genuinely 0.0
    // and Yes is the honest stamp -- the opposite of what every other
    // similarity configuration gets. A rule keyed on opts_.similarity would
    // return No here and would be wrong.
    auto methane = MakeConformer("C");
    ROCSOptions opts;
    opts.score_type = ROCSScoreType::Color;
    opts.similarity = true;
    ROCSComparison comparison({methane}, opts);
    // 1e-6, matching SELF_SCORE_TOLERANCE: at a looser tolerance a value in
    // (1e-6, 1e-3] would satisfy this assertion while the stamp below read No,
    // and the two assertions would contradict each other. Color Tanimoto with no
    // color atoms is structurally exact, so the tighter bound costs nothing.
    EXPECT_NEAR(comparison.Compare(0, 0), 0.0, 1e-6);
    EXPECT_EQ(comparison.Facts().zero_self, Capability::Yes);
}

TEST_F(ROCSComparisonTest, TwoDimensionalInputIsRefused) {
    auto flat = std::make_shared<OEChem::OEMol>();
    OEChem::OESmilesToMol(*flat, "c1ccc(O)cc1");
    OEChem::OEAddExplicitHydrogens(*flat);
    OEChem::OEGenerate2DCoordinates(*flat);
    ASSERT_EQ(flat->GetDimension(), 2u);
    EXPECT_THROW(ROCSComparison({flat}, ROCSOptions()), ComparisonError);
}

TEST_F(ROCSComparisonTest, CloneInheritsTheMeasuredStamp) {
    // The measurement is cached in SharedData precisely so clones do not repeat
    // it. This catches a clone that inherits no stamp and reports Unknown. It
    // cannot catch a clone that re-measures, which would score methane again and
    // reach the same No by the expensive route.
    auto methane = MakeConformer("C");
    ROCSComparison comparison({methane}, ROCSOptions());
    auto clone = comparison.Clone();
    EXPECT_EQ(clone->Facts().zero_self, Capability::No);
}

TEST_F(ROCSComparisonTest, ANonFiniteDiagonalStampsNoAndFlagsNaN) {
    // std::abs(NaN) > tol is false, so a naive threshold reads a NaN diagonal as
    // "vanished" and stamps a tier-1 Yes on a matrix that is not a number.
    auto mol = MakeConformer("c1ccc(O)cc1");
    OESystem::OEIter<OEChem::OEAtomBase> atom = mol->GetAtoms();
    ASSERT_TRUE(atom);
    float coords[3];
    ASSERT_TRUE(mol->GetCoords(&*atom, coords));
    coords[0] = std::numeric_limits<float>::quiet_NaN();
    ASSERT_TRUE(mol->SetCoords(&*atom, coords));
    ROCSComparison comparison({mol}, ROCSOptions());
    ASSERT_FALSE(std::isfinite(comparison.Compare(0, 0)));
    EXPECT_EQ(comparison.Facts().zero_self, Capability::No);
    EXPECT_EQ(comparison.Facts().data_integrity, DataIntegrity::NaNPresent);
}

TEST_F(ROCSComparisonTest, RealCoordinatesWithAStaleDimensionAreAdmitted) {
    // SetCoords does not refresh the dimension attribute, so a molecule with a
    // perfectly good conformer can report 0. Refusing it would reject input that
    // scores exactly right, and would advise a remedy the caller does not need.
    auto mol = MakeConformer("c1ccc(O)cc1");
    ASSERT_TRUE(mol->SetDimension(0));
    ROCSComparison comparison({mol}, ROCSOptions());
    EXPECT_NEAR(comparison.Compare(0, 0), 0.0, 1e-6);
    EXPECT_EQ(comparison.Facts().zero_self, Capability::Yes);
}

TEST_F(ROCSComparisonTest, ZeroCoordinatesAreRefusedDespiteAThreeDimensionAttribute) {
    // An SDF written with all-zero coordinates round-trips to dimension 3 and
    // scores 1.0 for every pair including the diagonal. The attribute alone
    // cannot catch it; recomputing from the coordinates can.
    auto mol = MakeConformer("c1ccc(O)cc1");
    const float origin[3] = {0.0f, 0.0f, 0.0f};
    for (OESystem::OEIter<OEChem::OEAtomBase> atom = mol->GetAtoms(); atom; ++atom) {
        ASSERT_TRUE(mol->SetCoords(&*atom, origin));
    }
    ASSERT_TRUE(mol->SetDimension(3));
    EXPECT_THROW(ROCSComparison({mol}, ROCSOptions()), ComparisonError);
}

TEST_F(ROCSComparisonTest, FactsCoverEveryScoreTypeAndDirection) {
    // zero_self is measured on this fixture at construction rather than derived
    // from the direction flag, so these rows record what benzene, phenol and
    // octane actually do: every distance form vanishes on the diagonal and every
    // similarity form saturates instead. The stamps are a property of these three
    // molecules, not of the score types -- methane stamps differently, which the
    // tests above pin down. The Color similarity cell used to read Yes, as an
    // artifact of a color term that was identically zero, and correctly reads No
    // now that it is not.
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
