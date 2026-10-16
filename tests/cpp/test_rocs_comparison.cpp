#include <gtest/gtest.h>
#include "oecluster/oecluster.h"
#include "oecluster/comparisons/ROCSComparison.h"
// Private to the ROCS implementation and reached by relative path, as
// test_storage.cpp reaches StorageOverflow.h: the test target has no ``src``
// include directory, and the helper is deliberately given external linkage so
// that its non-finite rule can be tested on numbers rather than on an overlay
// fixture that may not be able to produce them.
#include "../../src/comparisons/ROCSScoreCombine.h"
#include <oechem.h>
#include <oeomega2.h>
#include <oeshape.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

using namespace OECluster;

class ROCSComparisonTest : public ::testing::Test {
protected:
    /// An Omega ensemble of up to ``max_confs`` conformers for ``smi``, one by
    /// default. Shared by SetUp and by the tests that need molecules outside the
    /// standard fixture, so that every molecule in this file is built under one
    /// Omega configuration.
    ///
    /// Omega, not OEGenerate2DCoordinates: the overlay finds no coordinates it
    /// can use on planar input, warns, and returns a Tanimoto of exactly 0.0 for
    /// every pair, so the whole suite would measure saturation values rather
    /// than shape.
    static std::shared_ptr<OEChem::OEMol> MakeConformer(const char* smi,
                                                        unsigned int max_confs = 1) {
        OEConfGen::OEOmega omega;
        omega.SetMaxConfs(max_confs);
        omega.SetStrictStereo(false);
        auto mol = std::make_shared<OEChem::OEMol>();
        OEChem::OESmilesToMol(*mol, smi);
        // EXPECT_TRUE rather than ASSERT_TRUE: ASSERT_* expands to a bare
        // return, which cannot compile in a value-returning function. A failed
        // build leaves the molecule empty, which surfaces loudly downstream.
        EXPECT_TRUE(omega(*mol)) << smi;
        return mol;
    }

    /// Collapse one conformer of ``mol`` onto z = 0. ``flatten_active`` selects
    /// between the conformer OEShape actually reads and the first one it does
    /// not, which is the distinction the dimension guard turns on.
    static void FlattenOneConformer(OEChem::OEMol& mol, bool flatten_active) {
        const OEChem::OEConfBase* const active = mol.GetActive();
        for (OESystem::OEIter<OEChem::OEConfBase> conf = mol.GetConfs(); conf; ++conf) {
            if ((&*conf == active) != flatten_active) {
                continue;
            }
            for (OESystem::OEIter<OEChem::OEAtomBase> atom = conf->GetAtoms(); atom; ++atom) {
                float coords[3];
                EXPECT_TRUE(conf->GetCoords(&*atom, coords));
                coords[2] = 0.0f;
                EXPECT_TRUE(conf->SetCoords(&*atom, coords));
            }
            // One non-active conformer is the whole point of that case; leave
            // the rest of the ensemble sound so the molecule still has a
            // conformer worth overlaying.
            if (!flatten_active) {
                return;
            }
        }
    }

    /// Move one atom of one conformer of ``mol`` to ``offset`` on ``axis``,
    /// stretching that conformer's bounding box to roughly ``offset``.
    /// ``stretch_active`` selects between the conformer OEShape reads through
    /// the OEMolBase view and the first one it does not, which is the
    /// distinction the extent guard has to be insensitive to.
    ///
    /// Axis 0 by default, and deliberately. Measured through oecluster.pdist
    /// before this guard existed, the same 1e10 displacement crashed on x (exit
    /// 139) but returned a saturated 1.0 on y and on z. Neither outcome is a
    /// score, but only the x case proves the guard is what stops the crash.
    static void StretchOneConformer(OEChem::OEMol& mol, bool stretch_active,
                                    double offset, unsigned int axis = 0) {
        const OEChem::OEConfBase* const active = mol.GetActive();
        for (OESystem::OEIter<OEChem::OEConfBase> conf = mol.GetConfs(); conf; ++conf) {
            if ((&*conf == active) != stretch_active) {
                continue;
            }
            OESystem::OEIter<OEChem::OEAtomBase> atom = conf->GetAtoms();
            EXPECT_TRUE(atom);
            double coords[3];
            EXPECT_TRUE(conf->GetCoords(&*atom, coords));
            coords[axis] = offset;
            EXPECT_TRUE(conf->SetCoords(&*atom, coords));
            return;
        }
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

// The value of this test is that it compiles. Adding the OEMolBase* overload
// made all three braced forms below ambiguous -- either vector type can be
// brace-initialized from them, so no candidate wins -- and a source break is
// invisible to a suite that never spells the call that way. If the ambiguity
// returns, this file stops building.
TEST_F(ROCSComparisonTest, BracedInitializersStayUnambiguous) {
    ROCSComparison empty({});
    EXPECT_EQ(empty.Size(), 0u);

    EXPECT_THROW(ROCSComparison comparison({nullptr}), ComparisonError);

    // A braced list of real molecules binds the initializer-list overload
    // rather than the vector one it used to reach, so the score is pinned
    // against the vector's. A multiconformer molecule is used deliberately:
    // it is what makes the two distinguishable, since the wrong delegation
    // would reach the OEMolBase view and leave BestOverlay one pose instead of
    // an ensemble. The last assertion measures that the input really does
    // separate the two paths, so the equality above cannot hold vacuously.
    const std::shared_ptr<OEChem::OEMol> flexible = MakeConformer("c1ccccc1CCCCc1ccccc1", 3);
    ASSERT_GT(flexible->NumConfs(), 1u);
    const std::vector<std::shared_ptr<OEChem::OEMol>> pair{mols_[0], flexible};

    ROCSComparison via_vector(pair);
    ROCSComparison braced({mols_[0], flexible});
    EXPECT_DOUBLE_EQ(braced.Compare(0, 1), via_vector.Compare(0, 1));

    std::vector<OEChem::OEMolBase*> as_base{&static_cast<OEChem::OEMolBase&>(*mols_[0]),
                                            &static_cast<OEChem::OEMolBase&>(*flexible)};
    ROCSComparison collapsed(as_base);
    EXPECT_NE(collapsed.Compare(0, 1), via_vector.Compare(0, 1));
}

// The ``OEMolBase*`` overload's own null check, in ``to_oemol_snapshots``.
// Nothing else in the suite reaches it: the braced ``{nullptr}`` case above
// binds the initializer-list constructor, which outranks both vector ones, so
// it enters the ``shared_ptr`` shim and stops at the *strict* constructor's
// check instead; and the bindings refuse ``None`` in the typemap before C++
// sees it. Without this case, deleting the check leaves the whole suite green
// while a legal C++ call shape dereferences null.
//
// The vector is named rather than braced for exactly that reason. Brace it and
// the call routes back to the initializer-list shim and asserts nothing new,
// so this is not a spelling to tidy away.
TEST_F(ROCSComparisonTest, NullRawPointerIsRejected) {
    std::vector<OEChem::OEMolBase*> as_base{&static_cast<OEChem::OEMolBase&>(*mols_[0]), nullptr};
    EXPECT_THROW(ROCSComparison comparison(as_base), ComparisonError);
}

TEST_F(ROCSComparisonTest, CompareRefusesAnIndexPastTheEnd) {
    // SetupRef dereferences ``*shared_->mols[i]`` before the overlay runs, so
    // an out-of-range index is an out-of-bounds read on the molecule vector.
    ROCSComparison comparison(mols_);
    try {
        comparison.Compare(0, 1000000);
        FAIL() << "expected ComparisonError";
    } catch (const ComparisonError& error) {
        const std::string message(error.what());
        EXPECT_NE(message.find("1000000"), std::string::npos) << message;
        EXPECT_NE(message.find("3 items"), std::string::npos) << message;
    }
    EXPECT_THROW(comparison.Compare(1000000, 0), ComparisonError);
    EXPECT_THROW(comparison.Compare(mols_.size(), mols_.size()), ComparisonError);
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

TEST_F(ROCSComparisonTest, ANonFiniteCoordinateIsRefused) {
    // Admitted until the coordinate guard existed, and it scored: every pair
    // including the diagonal came back NaN, stamped honestly as
    // DataIntegrity::NaNPresent but still handed to the caller as a matrix. The
    // guard turns that into a refusal at construction. MeasureDiagonal's own
    // non-finite branch stays as a backstop -- ``std::abs(NaN) > tol`` is false,
    // so a naive threshold there would stamp a tier-1 Yes on a matrix that is
    // not a number -- but no input reaching it through this constructor can
    // carry a NaN coordinate any more.
    auto mol = MakeConformer("c1ccc(O)cc1");
    OESystem::OEIter<OEChem::OEAtomBase> atom = mol->GetAtoms();
    ASSERT_TRUE(atom);
    float coords[3];
    ASSERT_TRUE(mol->GetCoords(&*atom, coords));
    coords[0] = std::numeric_limits<float>::quiet_NaN();
    ASSERT_TRUE(mol->SetCoords(&*atom, coords));
    EXPECT_THROW(ROCSComparison({mol}, ROCSOptions()), ComparisonError);
}

TEST_F(ROCSComparisonTest, ANonFiniteCoordinateRefusalNamesItsMoleculeIndex) {
    // Methane ahead of the broken phenol so the message has to identify which
    // molecule to fix rather than reporting the first index it looked at. Sound
    // input in front of the offender is also what distinguishes a guard that
    // scans the whole set from one that stops early.
    auto methane = MakeConformer("C");
    auto broken = MakeConformer("c1ccc(O)cc1");
    OESystem::OEIter<OEChem::OEAtomBase> atom = broken->GetAtoms();
    ASSERT_TRUE(atom);
    float coords[3];
    ASSERT_TRUE(broken->GetCoords(&*atom, coords));
    coords[0] = std::numeric_limits<float>::quiet_NaN();
    ASSERT_TRUE(broken->SetCoords(&*atom, coords));

    try {
        ROCSComparison comparison({methane, broken}, ROCSOptions());
        FAIL() << "expected ComparisonError for a non-finite coordinate";
    } catch (const ComparisonError& exc) {
        const std::string message(exc.what());
        EXPECT_NE(message.find("index 1"), std::string::npos) << message;
        EXPECT_NE(message.find("non-finite"), std::string::npos) << message;
    }
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
    // The refresh has to land on the comparison's snapshot. Moving
    // OESetDimensionFromCoords above the copy loop would still pass every other
    // test in this file: InputMoleculesAreNotMutated only compares NumAtoms().
    EXPECT_EQ(mol->GetDimension(), 0u);
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

TEST_F(ROCSComparisonTest, OnlyTheActiveConformerDecidesTheDimensionRefusal) {
    // OESetDimensionFromCoords takes an OEMolBase&, and the OEMolBase view of a
    // multiconformer OEMol is its active conformer, so the refresh never sees
    // the rest of the ensemble. Both halves matter and neither is obvious, so
    // they are pinned together: the guard is permissive about the conformers
    // OEShape will not read, and unconditional about the one it will.
    auto ensemble = MakeConformer("CCCCCCCC", 6);
    ASSERT_GE(ensemble->NumConfs(), 2u);

    // A degenerate non-active conformer is admitted, and scores as though it
    // were not there: BestOverlay picks a sound conformer and seats octane
    // exactly on itself.
    auto non_active_flat = std::make_shared<OEChem::OEMol>(*ensemble);
    FlattenOneConformer(*non_active_flat, false);
    ROCSComparison admitted({non_active_flat}, ROCSOptions());
    EXPECT_NEAR(admitted.Compare(0, 0), 0.0, 1e-6);
    EXPECT_EQ(admitted.Facts().zero_self, Capability::Yes);

    // A degenerate active conformer is refused, and the five sound conformers
    // behind it do not rescue the molecule: OEShape scores it 0.0 against
    // everything, as reference and as fit alike.
    auto active_flat = std::make_shared<OEChem::OEMol>(*ensemble);
    FlattenOneConformer(*active_flat, true);
    EXPECT_THROW(ROCSComparison({active_flat}, ROCSOptions()), ComparisonError);
}

TEST_F(ROCSComparisonTest, AnExtremeCoordinateExtentIsRefused) {
    // Before the guard this exact input took the process down: driven through
    // oecluster.pdist it exited 139 (SIGSEGV) rather than raising, so the
    // refusal is the difference between an error and a crash, not between two
    // error messages.
    auto mol = MakeConformer("c1ccc(O)cc1");
    StretchOneConformer(*mol, true, 1e10);
    ASSERT_EQ(mol->GetDimension(), 3u) << "the dimension guard must not be what refuses this";

    try {
        ROCSComparison comparison({mol}, ROCSOptions());
        FAIL() << "expected ComparisonError for an extreme coordinate extent";
    } catch (const ComparisonError& exc) {
        const std::string message(exc.what());
        EXPECT_NE(message.find("index 0"), std::string::npos) << message;
        EXPECT_NE(message.find("extent"), std::string::npos) << message;
    }
}

TEST_F(ROCSComparisonTest, AStretchedNonActiveConformerIsRefused) {
    // The counterpart of OnlyTheActiveConformerDecidesTheDimensionRefusal, and
    // deliberately the opposite verdict. A flat non-active conformer is harmless
    // because BestOverlay picks a sound one; a stretched non-active conformer is
    // not, because BestOverlay grids the whole fit ensemble before it chooses.
    // Driven through oecluster.pdist before the guard existed, this exact damage
    // exited 138 (SIGBUS) even though the rest of the ensemble was sound and the
    // dimension guard never looks past the active conformer. That is why the
    // extent guard scans the whole ensemble where the guard above it does not.
    auto ensemble = MakeConformer("CCCCCCCC", 6);
    ASSERT_GE(ensemble->NumConfs(), 2u);
    StretchOneConformer(*ensemble, false, 1e10);
    ASSERT_EQ(ensemble->GetDimension(), 3u);
    EXPECT_THROW(ROCSComparison({ensemble}, ROCSOptions()), ComparisonError);
}

TEST_F(ROCSComparisonTest, ALargeButPhysicallyRealExtentIsAdmitted) {
    // 200 angstroms is the span of a large protein: an order of magnitude under
    // MAX_COORDINATE_EXTENT and well inside what OEShape can grid. Refusing it
    // would be as much of a defect as admitting the 1e10 case above, so this
    // asserts more than "did not throw" -- the self-overlay lands exactly on
    // 0.0, which for ComboNorm means a combo of 2.0, so the molecule is
    // genuinely overlaid rather than handed the 1.0 an unusable grid returns.
    //
    // The cross pair is deliberately not asserted. Measured on this fixture it
    // is 0.546 on an overlay being used for the first time and exactly 1.0 on
    // one that has scored anything before it, a split that appears somewhere
    // between 100 and 200 angstroms. That belongs to OEOverlay reuse rather than
    // to this guard, so pinning a number here would pin the call order instead.
    auto stretched = MakeConformer("c1ccc(O)cc1");
    StretchOneConformer(*stretched, true, 200.0);
    auto benzene = MakeConformer("c1ccccc1");

    ROCSComparison comparison({stretched, benzene}, ROCSOptions());
    EXPECT_NEAR(comparison.Compare(0, 0), 0.0, 1e-6);
    EXPECT_EQ(comparison.Facts().zero_self, Capability::Yes);
}

TEST_F(ROCSComparisonTest, ARigidlyTranslatedFrameIsAdmitted) {
    // The case that separates extent from magnitude. Every coordinate here is
    // around 1e6, six orders of magnitude past anything a guard on
    // ``abs(coordinate)`` would tolerate, but the bounding box is untouched and
    // OEShape scores it as though it had never moved. Asserting agreement with
    // the untranslated pair rather than a hard-coded number keeps the test
    // pinned to that invariance.
    auto phenol = MakeConformer("c1ccc(O)cc1");
    auto benzene = MakeConformer("c1ccccc1");
    auto translated = std::make_shared<OEChem::OEMol>(*phenol);
    for (OESystem::OEIter<OEChem::OEConfBase> conf = translated->GetConfs(); conf; ++conf) {
        for (OESystem::OEIter<OEChem::OEAtomBase> atom = conf->GetAtoms(); atom; ++atom) {
            double coords[3];
            ASSERT_TRUE(conf->GetCoords(&*atom, coords));
            coords[0] += 1e6;
            ASSERT_TRUE(conf->SetCoords(&*atom, coords));
        }
    }

    ROCSComparison reference({phenol, benzene}, ROCSOptions());
    ROCSComparison moved({translated, benzene}, ROCSOptions());
    EXPECT_NEAR(moved.Compare(0, 1), reference.Compare(0, 1), 1e-2);
    // Not vacuous only if the shared value is a real score rather than a
    // saturation constant both sides could reach by failing identically.
    EXPECT_GT(reference.Compare(0, 1), 0.0);
    EXPECT_LT(reference.Compare(0, 1), 1.0);
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

/// Reaches the protected seam ``pdist`` calls. ``CompareUnordered`` is
/// protected so that it cannot become a second public way to compare two
/// items; deriving is the access the language offers a caller that genuinely
/// needs it, and the call it makes is the one ``pdist`` makes through its
/// friendship. That the ``using`` declaration below compiles is itself part of
/// the test: it would not if the override had been left public, because the
/// base declares it protected.
struct ROCSUnorderedProbe : ROCSComparison {
    using ROCSComparison::ROCSComparison;
    using ROCSComparison::CompareUnordered;
};

TEST_F(ROCSComparisonTest, UnorderedKeepsTheBetterOverlayForEveryScoreType) {
    // Eight cases rather than one. The two orientations are combined on the raw
    // Tanimoto, and the conversion from that to the reported value differs per
    // score type and reverses direction on the similarity flag, so a
    // single-type test would prove an eighth of the rule. An implementation
    // that took the maximum of the converted values instead of the raw ones
    // would pass all four similarity rows and fail all four distance rows.
    const std::vector<std::shared_ptr<OEChem::OEMol>> flexible{
        MakeConformer("CCCCCCCC", 20), MakeConformer("c1ccc(CCCCO)cc1", 20)};
    ASSERT_GT(flexible[0]->NumConfs(), 1u);
    ASSERT_GT(flexible[1]->NumConfs(), 1u);

    // Not exact agreement: the overlay keeps state between calls, so the pair
    // of searches run for {0, 1} is not bit-identical to the pair run for
    // {1, 0}, nor to the two standalone Compare calls below.
    constexpr double tolerance = 1e-6;
    double largest_orientation_gap = 0.0;

    for (ROCSScoreType score_type : {ROCSScoreType::ComboNorm, ROCSScoreType::Combo,
                                     ROCSScoreType::Shape, ROCSScoreType::Color}) {
        for (bool similarity : {false, true}) {
            ROCSOptions opts;
            opts.score_type = score_type;
            opts.similarity = similarity;
            ROCSUnorderedProbe probe(flexible, opts);

            const double unordered = probe.CompareUnordered(0, 1);
            const double swapped = probe.CompareUnordered(1, 0);
            const double forward = probe.Compare(0, 1);
            const double reverse = probe.Compare(1, 0);

            const std::string label =
                std::to_string(static_cast<int>(score_type)) +
                (similarity ? " similarity" : " distance");

            EXPECT_NEAR(unordered, swapped, tolerance) << label;
            if (similarity) {
                EXPECT_GE(unordered, std::max(forward, reverse) - tolerance) << label;
            } else {
                EXPECT_LE(unordered, std::min(forward, reverse) + tolerance) << label;
            }

            largest_orientation_gap =
                std::max(largest_orientation_gap, std::abs(forward - reverse));
        }
    }

    // The assertions above hold for free on a fixture whose two orientations
    // agree, so the fixture has to be shown to disagree. Measured 0.070 on
    // 2026-10-15, on the Shape rows.
    EXPECT_GT(largest_orientation_gap, 0.01);
}

TEST_F(ROCSComparisonTest, UnorderedAgreesWithCompareOnTheDiagonal) {
    // pdist stores no diagonal, so this is for a future caller. It also keeps
    // the self-score a direct caller sees here equal to the one MeasureDiagonal
    // stamps zero_self from, which reads Compare.
    ROCSUnorderedProbe probe(mols_, ROCSOptions());
    for (size_t i = 0; i < mols_.size(); ++i) {
        EXPECT_DOUBLE_EQ(probe.CompareUnordered(i, i), probe.Compare(i, i)) << i;
    }
}

TEST(ROCSScoreCombine, NonFiniteOperandsCollapseToOneCanonicalNaN) {
    // std::max compares false against a NaN and so returns whichever operand it
    // was handed first, which would leave the combine order-dependent on
    // exactly the inputs it exists to make order-independent. Every row is
    // asserted in both argument orders for that reason.
    constexpr double quiet_nan = std::numeric_limits<double>::quiet_NaN();
    constexpr double infinity = std::numeric_limits<double>::infinity();

    struct Row {
        double a;
        double b;
        const char* why;
    };
    const std::vector<Row> rows{
        {quiet_nan, 0.5, "the basic case"},
        {quiet_nan, -infinity, "two non-finite values, neither preferable"},
        {infinity, 0.5,
         "a rule testing isnan rather than !isfinite passes every NaN row "
         "above and still returns inf here"},
        {-infinity, 0.5, "the same hole on the other side"},
    };

    for (const Row& row : rows) {
        EXPECT_TRUE(std::isnan(better_raw_score(row.a, row.b))) << row.why;
        EXPECT_TRUE(std::isnan(better_raw_score(row.b, row.a))) << row.why;
    }
}

TEST(ROCSScoreCombine, FiniteOperandsTakeTheLargerInEitherOrder) {
    // Larger is better for every ROCS raw score: all four are Tanimotos that
    // rise with overlay quality, which is what lets one combine serve the
    // distance and the similarity forms alike.
    EXPECT_DOUBLE_EQ(better_raw_score(0.25, 0.75), 0.75);
    EXPECT_DOUBLE_EQ(better_raw_score(0.75, 0.25), 0.75);
    EXPECT_DOUBLE_EQ(better_raw_score(-1.0, -2.0), -1.0);
    EXPECT_DOUBLE_EQ(better_raw_score(-2.0, -1.0), -1.0);
    EXPECT_DOUBLE_EQ(better_raw_score(0.5, 0.5), 0.5);
}
