#include <cmath>
#include <gtest/gtest.h>
#include <oechem.h>
#include "oecluster/Error.h"
#include "oecluster/PDist.h"
#include "oecluster/StorageBackend.h"
#include "oecluster/comparisons/RMSDComparison.h"

using namespace OECluster;

namespace {

/// Build a 2D-embedded molecule translated along x, so the in-frame RMSD
/// against the untranslated copy is exactly the shift.
std::shared_ptr<OEChem::OEMol> make_shifted(const char* smiles, double shift) {
    auto mol = std::make_shared<OEChem::OEMol>();
    OEChem::OESmilesToMol(*mol, smiles);
    OEChem::OEGenerate2DCoordinates(*mol);
    for (OESystem::OEIter<OEChem::OEAtomBase> atom = mol->GetAtoms(); atom; ++atom) {
        float coords[3] = {0.0f, 0.0f, 0.0f};
        mol->GetCoords(atom, coords);
        coords[0] += static_cast<float>(shift);
        mol->SetCoords(atom, coords);
    }
    return mol;
}

/// The atomic numbers of a molecule's atoms in index order, so a test can prove
/// that a step meant to reorder atoms actually reordered them.
std::vector<unsigned int> element_sequence(const OEChem::OEMol& mol) {
    std::vector<unsigned int> elements;
    for (OESystem::OEIter<OEChem::OEAtomBase> atom = mol.GetAtoms(); atom; ++atom) {
        elements.push_back(atom->GetAtomicNum());
    }
    return elements;
}

std::vector<float> all_coords(const std::vector<std::shared_ptr<OEChem::OEMol>>& mols) {
    std::vector<float> flat;
    for (const auto& mol : mols) {
        for (OESystem::OEIter<OEChem::OEAtomBase> atom = mol->GetAtoms(); atom; ++atom) {
            float coords[3] = {0.0f, 0.0f, 0.0f};
            mol->GetCoords(atom, coords);
            flat.push_back(coords[0]);
            flat.push_back(coords[1]);
            flat.push_back(coords[2]);
        }
    }
    return flat;
}

std::vector<double> run_pdist(RMSDComparison& comparison, size_t num_threads) {
    DenseStorage storage(comparison.Size());
    PDistOptions options;
    options.num_threads = num_threads;
    options.chunk_size = 1;  // Force multiple chunks to exercise concurrency.
    pdist(comparison, storage, options);
    std::vector<double> values;
    for (size_t i = 0; i < comparison.Size(); ++i) {
        for (size_t j = i + 1; j < comparison.Size(); ++j) {
            values.push_back(storage.Get(i, j));
        }
    }
    return values;
}

}  // namespace

class RMSDComparisonTest : public ::testing::Test {
protected:
    void SetUp() override {
        for (int i = 0; i < 4; ++i) {
            mols_.push_back(make_shifted("c1ccccc1", static_cast<double>(i)));
        }
    }

    std::vector<std::shared_ptr<OEChem::OEMol>> mols_;
};

TEST_F(RMSDComparisonTest, SelfDistanceIsZero) {
    RMSDComparison comparison(mols_);
    for (size_t i = 0; i < mols_.size(); ++i) {
        EXPECT_NEAR(comparison.Compare(i, i), 0.0, 1e-6);
    }
}

TEST_F(RMSDComparisonTest, InFrameRMSDIsTheTranslation) {
    RMSDOptions opts;
    opts.automorph = false;  // a translated hexagon has cheaper symmetry matches
    RMSDComparison comparison(mols_, opts);
    for (size_t i = 0; i < mols_.size(); ++i) {
        for (size_t j = i + 1; j < mols_.size(); ++j) {
            const double expected = static_cast<double>(j - i);
            EXPECT_NEAR(comparison.Compare(i, j), expected, 1e-4);
        }
    }
}

TEST_F(RMSDComparisonTest, OverlayRemovesTheTranslation) {
    RMSDOptions opts;
    opts.overlay = true;
    RMSDComparison comparison(mols_, opts);
    EXPECT_NEAR(comparison.Compare(0, 3), 0.0, 1e-4);
}

TEST_F(RMSDComparisonTest, OverlayDiffersFromInFrame) {
    RMSDOptions in_frame;
    RMSDOptions overlaid;
    overlaid.overlay = true;
    RMSDComparison a(mols_, in_frame);
    RMSDComparison b(mols_, overlaid);
    EXPECT_GT(std::abs(a.Compare(0, 3) - b.Compare(0, 3)), 1e-3);
}

TEST_F(RMSDComparisonTest, PDistDoesNotModifyTheInput) {
    const std::vector<float> before = all_coords(mols_);
    RMSDOptions opts;
    opts.overlay = true;
    RMSDComparison comparison(mols_, opts);
    (void)run_pdist(comparison, 4);
    EXPECT_EQ(all_coords(mols_), before);
}

TEST_F(RMSDComparisonTest, RepeatedPDistIsIdentical) {
    RMSDOptions opts;
    opts.overlay = true;
    RMSDComparison comparison(mols_, opts);
    const std::vector<double> first = run_pdist(comparison, 4);
    const std::vector<double> second = run_pdist(comparison, 4);
    EXPECT_EQ(first, second);
}

TEST_F(RMSDComparisonTest, ThreadCountDoesNotChangeTheResult) {
    RMSDOptions opts;
    opts.overlay = true;
    RMSDComparison comparison(mols_, opts);
    EXPECT_EQ(run_pdist(comparison, 1), run_pdist(comparison, 8));
}

TEST_F(RMSDComparisonTest, TopologyMismatchIsRejected) {
    std::vector<std::shared_ptr<OEChem::OEMol>> mixed = mols_;
    mixed.push_back(make_shifted("CCCCCCCC", 0.0));
    try {
        RMSDComparison comparison(mixed);
        FAIL() << "expected a topology mismatch to be rejected";
    } catch (const ComparisonError& exc) {
        const std::string message = exc.what();
        EXPECT_NE(message.find("rocs"), std::string::npos) << message;
        EXPECT_NE(message.find('4'), std::string::npos) << message;
    }
}

TEST_F(RMSDComparisonTest, NullMoleculeIsRejected) {
    std::vector<std::shared_ptr<OEChem::OEMol>> with_null = mols_;
    with_null[2].reset();
    EXPECT_THROW((RMSDComparison{with_null}), ComparisonError);
}

TEST_F(RMSDComparisonTest, FactsAreZeroSelfAndComplete) {
    RMSDComparison comparison(mols_);
    const GateFacts facts = comparison.Facts();
    EXPECT_EQ(facts.is_distance, Capability::Yes);
    EXPECT_EQ(facts.zero_self, Capability::Yes);
    // Unknown, not Yes: see the comment on RMSDComparison::Facts.
    EXPECT_EQ(facts.triangle, Capability::Unknown);
    EXPECT_EQ(facts.data_integrity, DataIntegrity::Complete);
}

TEST_F(RMSDComparisonTest, CloneScoresIdentically) {
    RMSDComparison comparison(mols_);
    const std::unique_ptr<PairwiseComparison> clone = comparison.Clone();
    EXPECT_EQ(clone->Size(), comparison.Size());
    EXPECT_NEAR(clone->Compare(0, 2), comparison.Compare(0, 2), 1e-9);
    EXPECT_EQ(clone->ComparisonName(), "rmsd");
}

TEST_F(RMSDComparisonTest, NonFiniteCoordinateIsRejected) {
    auto mol_with_nan = make_shifted("c1ccccc1", 0.0);
    // Set a NaN coordinate on the first atom.
    OESystem::OEIter<OEChem::OEAtomBase> atom = mol_with_nan->GetAtoms();
    float coords[3] = {std::nanf(""), 0.0f, 0.0f};
    mol_with_nan->SetCoords(atom, coords);

    std::vector<std::shared_ptr<OEChem::OEMol>> with_nan = {mol_with_nan, mols_[0]};
    RMSDComparison comparison(with_nan);
    try {
        comparison.Compare(0, 1);
        FAIL() << "expected a non-finite coordinate to be rejected";
    } catch (const ComparisonError& exc) {
        const std::string message = exc.what();
        EXPECT_NE(message.find("non-finite"), std::string::npos) << message;
        EXPECT_NE(message.find("coordinate"), std::string::npos) << message;
    }
}

TEST_F(RMSDComparisonTest, HeavyOnlySkipsHydrogens) {
    // Build a methanol molecule with explicit hydrogens.
    auto mol_original = std::make_shared<OEChem::OEMol>();
    OEChem::OESmilesToMol(*mol_original, "CO");
    OEChem::OEGenerate2DCoordinates(*mol_original);
    OEChem::OEAddExplicitHydrogens(*mol_original);

    // Clone it and displace only the hydrogens.
    auto mol_displaced_h = std::make_shared<OEChem::OEMol>(*mol_original);
    for (OESystem::OEIter<OEChem::OEAtomBase> atom = mol_displaced_h->GetAtoms(); atom; ++atom) {
        if (atom->GetAtomicNum() == 1) {  // hydrogen
            float coords[3] = {0.0f, 0.0f, 0.0f};
            mol_displaced_h->GetCoords(atom, coords);
            coords[0] += 3.0f;
            mol_displaced_h->SetCoords(atom, coords);
        }
    }

    std::vector<std::shared_ptr<OEChem::OEMol>> mols = {mol_original, mol_displaced_h};

    // With heavy_only=true, hydrogens are skipped and RMSD is zero.
    RMSDOptions opts_heavy;
    opts_heavy.heavy_only = true;
    RMSDComparison comparison_heavy(mols, opts_heavy);
    EXPECT_NEAR(comparison_heavy.Compare(0, 1), 0.0, 1e-6);

    // With heavy_only=false, the displaced hydrogens contribute and RMSD is nonzero.
    RMSDOptions opts_all;
    opts_all.heavy_only = false;
    RMSDComparison comparison_all(mols, opts_all);
    const double rmsd_all = comparison_all.Compare(0, 1);
    EXPECT_GT(rmsd_all, 0.0);
}

TEST_F(RMSDComparisonTest, UnembeddedMoleculesAreRejected) {
    // Build molecules from SMILES without embedding.
    auto mol_a = std::make_shared<OEChem::OEMol>();
    OEChem::OESmilesToMol(*mol_a, "CC");
    auto mol_b = std::make_shared<OEChem::OEMol>();
    OEChem::OESmilesToMol(*mol_b, "CC");

    std::vector<std::shared_ptr<OEChem::OEMol>> unembedded = {mol_a, mol_b};
    try {
        RMSDComparison comparison(unembedded);
        FAIL() << "expected unembedded molecules to be rejected";
    } catch (const ComparisonError& exc) {
        const std::string message = exc.what();
        EXPECT_NE(message.find("no coordinates"), std::string::npos) << message;
        EXPECT_NE(message.find("GetDimension"), std::string::npos) << message;
    }
}

TEST_F(RMSDComparisonTest, EmbeddedMixedWithUnembeddedIsRejected) {
    auto embedded = make_shifted("c1ccccc1", 0.0);
    auto unembedded = std::make_shared<OEChem::OEMol>();
    OEChem::OESmilesToMol(*unembedded, "c1ccccc1");

    std::vector<std::shared_ptr<OEChem::OEMol>> mixed = {embedded, unembedded};
    try {
        RMSDComparison comparison(mixed);
        FAIL() << "expected mixed embedded/unembedded to be rejected";
    } catch (const ComparisonError& exc) {
        const std::string message = exc.what();
        EXPECT_NE(message.find("no coordinates"), std::string::npos) << message;
        // The unembedded molecule is item 1, so the message must name that index
        // rather than blaming the embedded item 0.
        EXPECT_NE(message.find("index 1"), std::string::npos) << message;
    }
}

TEST_F(RMSDComparisonTest, DifferentAtomOrderingRejectedWithAutomorphOff) {
    // Build two ethanols with identical geometry per element but different atom ordering.
    auto mol_a = std::make_shared<OEChem::OEMol>();
    auto mol_b = std::make_shared<OEChem::OEMol>();

    // mol_a: NewAtom order C, C, O -> [6, 6, 8]
    auto c1_a = mol_a->NewAtom(6);
    auto c2_a = mol_a->NewAtom(6);
    auto o_a = mol_a->NewAtom(8);
    mol_a->NewBond(c1_a, c2_a, 1);
    mol_a->NewBond(c2_a, o_a, 1);
    // Manually set coordinates: C at (0,0,0), C at (1,0,0), O at (2,0,0)
    float origin[3] = {0.0f, 0.0f, 0.0f};
    float middle[3] = {1.0f, 0.0f, 0.0f};
    float far_end[3] = {2.0f, 0.0f, 0.0f};
    mol_a->SetCoords(c1_a, origin);
    mol_a->SetCoords(c2_a, middle);
    mol_a->SetCoords(o_a, far_end);
    mol_a->SetDimension(2);

    // mol_b: NewAtom order O, C, C -> [8, 6, 6]
    auto o_b = mol_b->NewAtom(8);
    auto c1_b = mol_b->NewAtom(6);
    auto c2_b = mol_b->NewAtom(6);
    mol_b->NewBond(c1_b, c2_b, 1);
    mol_b->NewBond(c2_b, o_b, 1);
    // Set the same geometry: O at (2,0,0), C at (0,0,0), C at (1,0,0)
    mol_b->SetCoords(o_b, far_end);
    mol_b->SetCoords(c1_b, origin);
    mol_b->SetCoords(c2_b, middle);
    mol_b->SetDimension(2);

    std::vector<std::shared_ptr<OEChem::OEMol>> different_order = {mol_a, mol_b};

    // With automorph=false, the different atom ordering is rejected at construction.
    RMSDOptions opts_no_auto;
    opts_no_auto.automorph = false;
    try {
        RMSDComparison comparison(different_order, opts_no_auto);
        FAIL() << "expected different atom ordering to be rejected with automorph=false";
    } catch (const ComparisonError& exc) {
        const std::string message = exc.what();
        EXPECT_NE(message.find("automorph=false"), std::string::npos) << message;
        EXPECT_NE(message.find("index"), std::string::npos) << message;
    }

    // With automorph=true, the pair is accepted and scores correctly.
    RMSDOptions opts_auto;
    opts_auto.automorph = true;
    RMSDComparison comparison_auto(different_order, opts_auto);
    EXPECT_NEAR(comparison_auto.Compare(0, 1), 0.0, 1e-6);
}

TEST_F(RMSDComparisonTest, OverlayAutomorphCombinationUnderThreads) {
    // Exercise the most expensive and stateful option combination under concurrency.
    RMSDOptions opts;
    opts.overlay = true;
    opts.automorph = true;
    RMSDComparison comparison(mols_, opts);
    EXPECT_EQ(run_pdist(comparison, 1), run_pdist(comparison, 8));
}

TEST_F(RMSDComparisonTest, ExplicitVsSuppressedHydrogensRejectedWithAutomorphOff) {
    // Build methanol with hydrogens suppressed.
    auto mol_suppressed = std::make_shared<OEChem::OEMol>();
    OEChem::OESmilesToMol(*mol_suppressed, "CO");
    OEChem::OEGenerate2DCoordinates(*mol_suppressed);

    // Build methanol with explicit hydrogens.
    auto mol_explicit = std::make_shared<OEChem::OEMol>();
    OEChem::OESmilesToMol(*mol_explicit, "CO");
    OEChem::OEGenerate2DCoordinates(*mol_explicit);
    OEChem::OEAddExplicitHydrogens(*mol_explicit);

    std::vector<std::shared_ptr<OEChem::OEMol>> mixed_h = {mol_suppressed, mol_explicit};

    // With automorph=false, the different atom counts are rejected at construction.
    RMSDOptions opts;
    opts.automorph = false;
    try {
        RMSDComparison comparison(mixed_h, opts);
        FAIL() << "expected differing hydrogen treatment to be rejected with automorph=false";
    } catch (const ComparisonError& exc) {
        const std::string message = exc.what();
        EXPECT_NE(message.find("automorph=false"), std::string::npos) << message;
        EXPECT_NE(message.find("atoms"), std::string::npos) << message;
    }
}

TEST_F(RMSDComparisonTest, MixedHydrogenRepresentationRejectedWhenHeavyOnlyIsFalse) {
    // Build methanol with hydrogens suppressed.
    auto mol_suppressed = std::make_shared<OEChem::OEMol>();
    OEChem::OESmilesToMol(*mol_suppressed, "CO");
    OEChem::OEGenerate2DCoordinates(*mol_suppressed);

    // Build methanol with explicit hydrogens, displaced so a non-zero distance is expected.
    auto mol_explicit = std::make_shared<OEChem::OEMol>();
    OEChem::OESmilesToMol(*mol_explicit, "CO");
    OEChem::OEGenerate2DCoordinates(*mol_explicit);
    OEChem::OEAddExplicitHydrogens(*mol_explicit);
    // Displace only the hydrogens.
    for (OESystem::OEIter<OEChem::OEAtomBase> atom = mol_explicit->GetAtoms(); atom; ++atom) {
        if (atom->GetAtomicNum() == 1) {
            float coords[3] = {0.0f, 0.0f, 0.0f};
            mol_explicit->GetCoords(atom, coords);
            coords[0] += 3.0f;
            mol_explicit->SetCoords(atom, coords);
        }
    }

    std::vector<std::shared_ptr<OEChem::OEMol>> mixed_h = {mol_suppressed, mol_explicit};

    // With automorph=true and heavy_only=false, the different hydrogen representations
    // are rejected because hydrogens are counted.
    RMSDOptions opts_reject;
    opts_reject.automorph = true;
    opts_reject.heavy_only = false;
    try {
        RMSDComparison comparison(mixed_h, opts_reject);
        FAIL() << "expected mixed hydrogen representation to be rejected with heavy_only=false";
    } catch (const ComparisonError& exc) {
        const std::string message = exc.what();
        EXPECT_NE(message.find("heavy_only=false"), std::string::npos) << message;
        EXPECT_NE(message.find("hydrogen"), std::string::npos) << message;
    }

    // With automorph=true and heavy_only=true, the pair is ACCEPTED (hydrogens are
    // excluded by request) and scores 0.0 because the heavy atoms coincide.
    RMSDOptions opts_accept;
    opts_accept.automorph = true;
    opts_accept.heavy_only = true;
    RMSDComparison comparison_accept(mixed_h, opts_accept);
    EXPECT_NEAR(comparison_accept.Compare(0, 1), 0.0, 1e-6);
}

TEST_F(RMSDComparisonTest, TwoDVsThreeDIsRejected) {
    // Build a 2D-embedded molecule.
    auto mol_2d = std::make_shared<OEChem::OEMol>();
    OEChem::OESmilesToMol(*mol_2d, "CCO");
    OEChem::OEGenerate2DCoordinates(*mol_2d);

    // Build a 3D molecule manually (avoid pulling in Omega).
    auto mol_3d = std::make_shared<OEChem::OEMol>();
    auto c1 = mol_3d->NewAtom(6);
    auto c2 = mol_3d->NewAtom(6);
    auto o = mol_3d->NewAtom(8);
    mol_3d->NewBond(c1, c2, 1);
    mol_3d->NewBond(c2, o, 1);
    float coords_c1[3] = {0.0f, 0.0f, 0.0f};
    float coords_c2[3] = {1.0f, 0.5f, 0.2f};
    float coords_o[3] = {2.0f, 0.3f, 0.7f};
    mol_3d->SetCoords(c1, coords_c1);
    mol_3d->SetCoords(c2, coords_c2);
    mol_3d->SetCoords(o, coords_o);
    mol_3d->SetDimension(3);

    std::vector<std::shared_ptr<OEChem::OEMol>> mixed_dim = {mol_2d, mol_3d};

    try {
        RMSDComparison comparison(mixed_dim);
        FAIL() << "expected 2D vs 3D to be rejected";
    } catch (const ComparisonError& exc) {
        const std::string message = exc.what();
        EXPECT_NE(message.find("2D"), std::string::npos) << message;
        EXPECT_NE(message.find("3D"), std::string::npos) << message;
    }
}

TEST_F(RMSDComparisonTest, AtomCountChangeAfterConstructionDoesNotChangeTheScore) {
    // Two methanols with explicit hydrogens, displaced on one side only, so under
    // heavy_only=false the whole score is carried by the hydrogens.
    auto mol_stable = make_shifted("CO", 0.0);
    OEChem::OEAddExplicitHydrogens(*mol_stable);
    auto mol_mutated = std::make_shared<OEChem::OEMol>(*mol_stable);
    for (OESystem::OEIter<OEChem::OEAtomBase> atom = mol_mutated->GetAtoms(); atom; ++atom) {
        if (atom->GetAtomicNum() == 1) {
            float coords[3] = {0.0f, 0.0f, 0.0f};
            mol_mutated->GetCoords(atom, coords);
            coords[0] += 3.0f;
            mol_mutated->SetCoords(atom, coords);
        }
    }

    std::vector<std::shared_ptr<OEChem::OEMol>> mols = {mol_stable, mol_mutated};
    RMSDOptions opts;
    opts.automorph = true;
    opts.heavy_only = false;
    RMSDComparison comparison(mols, opts);
    const double before = comparison.Compare(0, 1);
    ASSERT_GT(before, 1.0);

    // The caller still owns these pointers. Suppressing the hydrogens reaches the
    // mixed-representation state the constructor rejects under heavy_only=false;
    // scoring the caller's molecules there drives OERMSD into its atom-matching
    // failure. The comparison holds its own copies, so the answer cannot move.
    ASSERT_TRUE(OEChem::OESuppressHydrogens(*mol_mutated));
    ASSERT_LT(mol_mutated->NumAtoms(), mol_stable->NumAtoms());

    EXPECT_DOUBLE_EQ(comparison.Compare(0, 1), before);
}

TEST_F(RMSDComparisonTest, DimensionChangeAfterConstructionDoesNotChangeTheScore) {
    auto mol_stable = make_shifted("CCO", 0.0);
    auto mol_mutated = make_shifted("CCO", 1.0);

    std::vector<std::shared_ptr<OEChem::OEMol>> mols = {mol_stable, mol_mutated};
    std::vector<std::shared_ptr<OEChem::OEMol>> swapped = {mol_mutated, mol_stable};
    RMSDComparison comparison(mols);
    // The same molecule at the other index, so both Compare arguments are covered.
    RMSDComparison swapped_comparison(swapped);
    const double before = comparison.Compare(0, 1);
    const double before_swapped = swapped_comparison.Compare(0, 1);
    ASSERT_NEAR(before, 1.0, 1e-4);

    // Re-embedding one item out of the plane reaches the 2D-versus-3D state the
    // constructor rejects, and would silently change the measured distance.
    for (OESystem::OEIter<OEChem::OEAtomBase> atom = mol_mutated->GetAtoms(); atom; ++atom) {
        float coords[3] = {0.0f, 0.0f, 0.0f};
        mol_mutated->GetCoords(atom, coords);
        coords[2] += 2.0f;
        mol_mutated->SetCoords(atom, coords);
    }
    mol_mutated->SetDimension(3);
    ASSERT_EQ(mol_mutated->GetDimension(), 3u);

    EXPECT_DOUBLE_EQ(comparison.Compare(0, 1), before);
    EXPECT_DOUBLE_EQ(swapped_comparison.Compare(0, 1), before_swapped);
}

TEST_F(RMSDComparisonTest, AtomReorderAfterConstructionDoesNotChangeTheScore) {
    // With automorph=false atoms are matched by index, so reordering the caller's
    // atoms is the mutation that most directly corrupts the score: it leaves the
    // atom count, the dimension and the canonical SMILES all untouched.
    auto mol_stable = make_shifted("C(=O)N", 0.0);
    auto mol_reordered = make_shifted("C(=O)N", 1.0);

    std::vector<std::shared_ptr<OEChem::OEMol>> mols = {mol_stable, mol_reordered};
    RMSDOptions opts;
    opts.automorph = false;
    RMSDComparison comparison(mols, opts);
    const double before = comparison.Compare(0, 1);
    ASSERT_NEAR(before, 1.0, 1e-4);

    const std::vector<unsigned int> elements_before = element_sequence(*mol_reordered);
    OEChem::OECanonicalOrderAtoms(*mol_reordered);
    // Guard the fixture: a reordering that did not reorder would prove nothing.
    ASSERT_NE(element_sequence(*mol_reordered), elements_before);
    ASSERT_EQ(mol_reordered->NumAtoms(), mol_stable->NumAtoms());

    EXPECT_DOUBLE_EQ(comparison.Compare(0, 1), before);
}

TEST_F(RMSDComparisonTest, AddingHydrogensAfterConstructionStillScoresOnTheDefaultPath) {
    // The default path deliberately admits mixed hydrogen representations, because
    // hydrogens are excluded from the score. A comparison built from this pair must
    // therefore keep scoring it after one side gains hydrogens, not refuse it.
    auto mol_stable = make_shifted("CCO", 0.0);
    auto mol_with_hydrogens = make_shifted("CCO", 0.0);

    std::vector<std::shared_ptr<OEChem::OEMol>> mols = {mol_stable, mol_with_hydrogens};
    RMSDComparison comparison(mols);
    ASSERT_NEAR(comparison.Compare(0, 1), 0.0, 1e-6);

    OEChem::OEAddExplicitHydrogens(*mol_with_hydrogens);
    ASSERT_GT(mol_with_hydrogens->NumAtoms(), mol_stable->NumAtoms());

    double value = -1.0;
    EXPECT_NO_THROW(value = comparison.Compare(0, 1));
    EXPECT_NEAR(value, 0.0, 1e-6);
}

TEST_F(RMSDComparisonTest, RepeatedScoringReturnsTheSameValues) {
    // Scoring one comparison twice must give identical numbers, including on the
    // heavy_only=false path where atom counts enter the score.
    RMSDOptions opts;
    opts.heavy_only = false;
    RMSDComparison comparison(mols_, opts);
    std::vector<double> passes[2];
    for (size_t pass = 0; pass < 2; ++pass) {
        for (size_t i = 0; i < mols_.size(); ++i) {
            for (size_t j = i; j < mols_.size(); ++j) {
                passes[pass].push_back(comparison.Compare(i, j));
            }
        }
    }
    EXPECT_EQ(passes[0], passes[1]);
}
