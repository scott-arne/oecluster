/**
 * @file test_mcs_comparison.cpp
 * @brief Maximum-common-substructure comparison.
 *
 * Every expected value is derived from the bond counts in the fixture table
 * below and nothing else. Bond counts are after
 * ``OESuppressHydrogens(mol, false, false, false)``.
 *
 * | benzene 6 | toluene 7 | cyclohexane 6 | testosterone 24 | morphine 25 |
 * | penicillin G 25 | sucrose 24 | macrolide fragment 29 | benzene-d1 6 |
 */

#include <algorithm>
#include <memory>
#include <string>
#include <vector>
#include <gtest/gtest.h>
#include <oechem.h>
#include "oecluster/Error.h"
#include "oecluster/PDist.h"
#include "oecluster/StorageBackend.h"
#include "oecluster/comparisons/MCSComparison.h"

using namespace OECluster;

namespace {

const char* const BENZENE = "c1ccccc1";
const char* const TOLUENE = "Cc1ccccc1";
const char* const CYCLOHEXANE = "C1CCCCC1";
const char* const TESTOSTERONE =
    "C[C@]12CC[C@H]3[C@@H](CCC4=CC(=O)CC[C@]34C)[C@@H]1CC[C@@H]2O";
const char* const MORPHINE = "CN1CC[C@]23c4c5ccc(O)c4O[C@H]2[C@@H](O)C=C[C@H]3[C@H]1C5";
const char* const PENICILLIN_G = "CC1(C)S[C@@H]2[C@H](NC(=O)Cc3ccccc3)C(=O)N2[C@H]1C(=O)O";
const char* const SUCROSE =
    "OC[C@H]1O[C@@](CO)(O[C@H]2O[C@H](CO)[C@@H](O)[C@H](O)[C@H]2O)[C@@H](O)[C@@H]1O";
const char* const MACROLIDE =
    "CC[C@H]1OC(=O)[C@H](C)[C@@H](O)[C@H](C)[C@@H](O)[C@](C)(O)C[C@@H](C)C(=O)"
    "[C@H](C)[C@@H](O)[C@]1(C)O";
const char* const BENZENE_D1 = "[2H]c1ccccc1";
const char* const METHANE = "C";

/// Parse a SMILES into a titled molecule. MCS is topological, so no embedding
/// step is needed and none is done.
std::shared_ptr<OEChem::OEMol> from_smiles(const char* smiles, const char* title = "mol") {
    auto mol = std::make_shared<OEChem::OEMol>();
    OEChem::OESmilesToMol(*mol, smiles);
    mol->SetTitle(title);
    return mol;
}

std::vector<std::shared_ptr<OEChem::OEMol>> pair_of(const char* first, const char* second) {
    return {from_smiles(first, "first"), from_smiles(second, "second")};
}

/// The upper triangle of a pdist run, in row-major order.
std::vector<double> run_pdist(MCSComparison& comparison, size_t num_threads) {
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

/// One directed approximate MCS search, written out here rather than reused
/// from the implementation. The symmetry test's first assertion has to be
/// independent evidence that the two directions genuinely disagree; taking that
/// evidence from the code under test would make it circular.
unsigned int directed_bonds(const char* pattern_smiles, const char* target_smiles) {
    OEChem::OEMol pattern;
    OEChem::OESmilesToMol(pattern, pattern_smiles);
    OEChem::OESuppressHydrogens(pattern, false, false, false);
    OEChem::OEMol target;
    OEChem::OESmilesToMol(target, target_smiles);
    OEChem::OESuppressHydrogens(target, false, false, false);

    OEChem::OEMCSSearch search(pattern, OEChem::OEExprOpts::DefaultAtoms,
                               OEChem::OEExprOpts::DefaultBonds, OEChem::OEMCSType::Approximate);
    search.SetMCSFunc(OEChem::OEMCSMaxBondsCompleteCycles(1.0));
    search.SetMaxMatches(1024);
    unsigned int best = 0;
    for (OESystem::OEIter<const OEChem::OEMatchBase> match = search.Match(target, /*umatch=*/true);
         match; ++match) {
        best = std::max<unsigned int>(best, match->NumBonds());
    }
    return best;
}

}  // namespace

class MCSComparisonTest : public ::testing::Test {
protected:
    void SetUp() override {
        mols_.push_back(from_smiles(BENZENE, "benzene"));
        mols_.push_back(from_smiles(TOLUENE, "toluene"));
        mols_.push_back(from_smiles(CYCLOHEXANE, "cyclohexane"));
        mols_.push_back(from_smiles(MORPHINE, "morphine"));
    }

    std::vector<std::shared_ptr<OEChem::OEMol>> mols_;
};

// --- Score correctness ------------------------------------------------------

TEST_F(MCSComparisonTest, SelfDistanceIsZero) {
    MCSComparison comparison(mols_);
    for (size_t i = 0; i < comparison.Size(); ++i) {
        EXPECT_DOUBLE_EQ(comparison.Compare(i, i), 0.0);
    }
}

TEST_F(MCSComparisonTest, IdenticalMoleculesScoreZero) {
    // Two separately parsed copies, so this goes through a real search rather
    // than the diagonal short-circuit.
    MCSComparison comparison(pair_of(MORPHINE, MORPHINE));
    EXPECT_NEAR(comparison.Compare(0, 1), 0.0, 1e-9);
}

TEST_F(MCSComparisonTest, DisjointMoleculesScoreOne) {
    // Benzene against cyclohexane at the default match level: aromatic bonds do
    // not match single bonds, so nothing matches at all.
    MCSComparison comparison(pair_of(BENZENE, CYCLOHEXANE));
    EXPECT_NEAR(comparison.Compare(0, 1), 1.0, 1e-9);
}

TEST_F(MCSComparisonTest, KnownPairScoresTheExpectedTanimoto) {
    // Benzene (6 bonds) against toluene (7): the ring matches, c = 6, so the
    // similarity is 6 / (6 + 7 - 6) = 6/7.
    MCSComparison comparison(pair_of(BENZENE, TOLUENE));
    EXPECT_NEAR(comparison.Compare(0, 1), 0.142857, 1e-6);
}

// --- Validation -------------------------------------------------------------

TEST_F(MCSComparisonTest, ZeroBondMoleculeIsRejected) {
    std::vector<std::shared_ptr<OEChem::OEMol>> mols = {from_smiles(BENZENE, "benzene"),
                                                        from_smiles(METHANE, "methane")};
    EXPECT_THROW(MCSComparison comparison(mols), ComparisonError);
}

TEST_F(MCSComparisonTest, MaxMatchesZeroIsRejected) {
    MCSOptions opts;
    opts.max_matches = 0;
    EXPECT_THROW(MCSComparison comparison(mols_, opts), ComparisonError);
}

TEST_F(MCSComparisonTest, NullMoleculeIsRejected) {
    std::vector<std::shared_ptr<OEChem::OEMol>> mols = {from_smiles(BENZENE, "benzene"), nullptr};
    EXPECT_THROW(MCSComparison comparison(mols), ComparisonError);
}

// The value of this test is that it compiles. Adding the OEMolBase* overload
// made all three braced forms below ambiguous -- either vector type can be
// brace-initialized from them, so no candidate wins -- and a source break is
// invisible to a suite that never spells the call that way. The runtime
// assertions are incidental; if the ambiguity returns, this file stops
// building.
TEST_F(MCSComparisonTest, BracedInitializersStayUnambiguous) {
    MCSComparison empty({});
    EXPECT_EQ(empty.Size(), 0u);

    EXPECT_THROW(MCSComparison comparison({nullptr}), ComparisonError);

    // A braced list of real molecules binds the initializer-list overload
    // rather than the vector one it used to reach, so the score is pinned to
    // confirm the delegation lands on the strict constructor: the same 6/7
    // KnownPairScoresTheExpectedTanimoto measures through the vector.
    MCSComparison braced({from_smiles(BENZENE, "first"), from_smiles(TOLUENE, "second")});
    EXPECT_EQ(braced.Size(), 2u);
    EXPECT_NEAR(braced.Compare(0, 1), 0.142857, 1e-6);
}

TEST_F(MCSComparisonTest, CompareRefusesAnIndexPastTheEnd) {
    MCSComparison comparison(mols_);
    EXPECT_THROW((void)comparison.Compare(0, mols_.size()), ComparisonError);
    EXPECT_THROW((void)comparison.Compare(mols_.size(), 0), ComparisonError);
}

TEST_F(MCSComparisonTest, InvalidSearchModeIsRejected) {
    MCSOptions opts;
    opts.search_mode = static_cast<MCSSearchMode>(7);
    EXPECT_THROW(MCSComparison comparison(mols_, opts), ComparisonError);
}

TEST_F(MCSComparisonTest, InvalidMatchLevelIsRejected) {
    MCSOptions opts;
    opts.match_level = static_cast<MCSMatchLevel>(9);
    EXPECT_THROW(MCSComparison comparison(mols_, opts), ComparisonError);
}

// --- Infrastructure ---------------------------------------------------------

TEST_F(MCSComparisonTest, FactsInDistanceMode) {
    MCSComparison comparison(mols_);
    const GateFacts facts = comparison.Facts();
    EXPECT_EQ(facts.is_distance, Capability::Yes);
    EXPECT_EQ(facts.zero_self, Capability::Yes);
    EXPECT_EQ(facts.triangle, Capability::Unknown);
    EXPECT_EQ(facts.data_integrity, DataIntegrity::Complete);
}

TEST_F(MCSComparisonTest, FactsInSimilarityMode) {
    // Paired with the distance-mode case above so neither flag can be
    // hardcoded: both orientations flip, and the other two do not.
    MCSOptions opts;
    opts.similarity = true;
    MCSComparison comparison(mols_, opts);
    const GateFacts facts = comparison.Facts();
    EXPECT_EQ(facts.is_distance, Capability::No);
    EXPECT_EQ(facts.zero_self, Capability::No);
    EXPECT_EQ(facts.triangle, Capability::Unknown);
    EXPECT_EQ(facts.data_integrity, DataIntegrity::Complete);
    // The similarity-mode diagonal is what zero_self = No reports.
    EXPECT_DOUBLE_EQ(comparison.Compare(0, 0), 1.0);
}

TEST_F(MCSComparisonTest, CloneScoresIdentically) {
    MCSComparison comparison(mols_);
    std::unique_ptr<PairwiseComparison> clone = comparison.Clone();
    ASSERT_EQ(clone->Size(), comparison.Size());
    EXPECT_EQ(clone->ComparisonName(), "mcs");
    for (size_t i = 0; i < comparison.Size(); ++i) {
        for (size_t j = i + 1; j < comparison.Size(); ++j) {
            EXPECT_DOUBLE_EQ(clone->Compare(i, j), comparison.Compare(i, j));
        }
    }
}

TEST_F(MCSComparisonTest, CloneOutlivesItsParent) {
    // A lifetime test, and worth being explicit about what it does not prove:
    // it passes under the aliasing model too, because a shared_ptr alias keeps
    // the parent's snapshots alive by itself. The copy-per-clone property is not
    // observable through the public surface at all -- the snapshots have no
    // accessor -- so it is enforced by review of Clone(), not by this.
    std::unique_ptr<PairwiseComparison> clone;
    double expected = 0.0;
    {
        MCSComparison comparison(pair_of(BENZENE, TOLUENE));
        expected = comparison.Compare(0, 1);
        clone = comparison.Clone();
    }
    EXPECT_NEAR(clone->Compare(0, 1), expected, 1e-9);
    EXPECT_NEAR(clone->Compare(0, 1), 0.142857, 1e-6);
}

// --- Symmetry ---------------------------------------------------------------

TEST_F(MCSComparisonTest, SymmetryTakesTheLargerDirection) {
    // The three pairs form a cycle in the winning direction: testosterone beats
    // morphine, morphine beats penicillin G, penicillin G beats testosterone.
    // The winner relation is therefore not transitive and no ordering of these
    // three molecules exists at all, so any implementation that picks its one
    // pattern by ranking the pair -- on canonical SMILES, bond count,
    // heavy-atom count, weight, address or anything else -- imposes a total
    // order and must get at least one row wrong in at least one list order.
    // Only running both searches passes. This is the only cycle among the 14
    // asymmetric pairs in the 190-pair scan the design ran, so the fixtures are
    // not substitutable.
    struct Case {
        const char* name;
        const char* first;
        const char* second;
        unsigned int first_as_pattern;
        unsigned int second_as_pattern;
        double expected;  ///< The max-derived distance.
        double wrong;     ///< What the losing direction alone would give.
    };
    const Case cases[] = {
        // 24 and 25 bonds: 8/(24+25-8) = 8/41 against 7/42.
        {"testosterone/morphine", TESTOSTERONE, MORPHINE, 8, 7, 0.804878, 0.833333},
        // 25 and 25 bonds: 11/(50-11) = 11/39 against 10/40.
        {"morphine/penicillinG", MORPHINE, PENICILLIN_G, 11, 10, 0.717949, 0.750000},
        // 25 and 24 bonds: 5/(49-5) = 5/44 against 4/45.
        {"penicillinG/testosterone", PENICILLIN_G, TESTOSTERONE, 5, 4, 0.886364, 0.911111},
    };

    for (const Case& test_case : cases) {
        SCOPED_TRACE(test_case.name);

        // 1. The two directions really do disagree. If a future toolkit makes
        //    any of these symmetric, the cycle is gone and the test says so
        //    loudly rather than quietly passing on nothing.
        EXPECT_EQ(directed_bonds(test_case.first, test_case.second), test_case.first_as_pattern);
        EXPECT_EQ(directed_bonds(test_case.second, test_case.first), test_case.second_as_pattern);
        EXPECT_NE(test_case.first_as_pattern, test_case.second_as_pattern);

        // 2. Compare returns the max-derived value, not the other direction's.
        MCSComparison forward(pair_of(test_case.first, test_case.second));
        EXPECT_NEAR(forward.Compare(0, 1), test_case.expected, 1e-6);
        EXPECT_GT(std::abs(test_case.expected - test_case.wrong), 1e-4);

        // 3. And it does so from the reversed list too, which rules out any
        //    canonicalization on the index pair.
        MCSComparison reversed(pair_of(test_case.second, test_case.first));
        EXPECT_NEAR(reversed.Compare(0, 1), test_case.expected, 1e-6);
    }
}

// --- Options observability --------------------------------------------------
//
// Every advertised control has a case that fails if the control is ignored.

TEST_F(MCSComparisonTest, LooseMatchesAcrossAromaticity) {
    // Benzene (6 bonds) against cyclohexane (6). Loose leaves bonds
    // unconstrained, so the whole ring matches: c = 6, similarity 6/6.
    MCSOptions loose;
    loose.match_level = MCSMatchLevel::Loose;
    MCSComparison loose_comparison(pair_of(BENZENE, CYCLOHEXANE), loose);
    EXPECT_NEAR(loose_comparison.Compare(0, 1), 0.000000, 1e-6);

    MCSComparison default_comparison(pair_of(BENZENE, CYCLOHEXANE));
    EXPECT_NEAR(default_comparison.Compare(0, 1), 1.000000, 1e-6);
}

TEST_F(MCSComparisonTest, ExactRequiresSubstitutionPattern) {
    // Benzene (6) against toluene (7). Exact adds hydrogen count and degree, so
    // toluene's substituted ring carbon no longer matches a benzene CH and the
    // match drops from the full ring to c = 4: 4/(6+7-4) = 4/9.
    MCSOptions exact;
    exact.match_level = MCSMatchLevel::Exact;
    MCSComparison exact_comparison(pair_of(BENZENE, TOLUENE), exact);
    EXPECT_NEAR(exact_comparison.Compare(0, 1), 0.555556, 1e-6);

    MCSComparison default_comparison(pair_of(BENZENE, TOLUENE));
    EXPECT_NEAR(default_comparison.Compare(0, 1), 0.142857, 1e-6);
}

TEST_F(MCSComparisonTest, ExhaustiveFindsALargerMatchThanApproximate) {
    // Sucrose (24 bonds) against a macrolide fragment (29). Exhaustive finds 17
    // matched bonds against approximate's 15: 17/36 against 15/38. This pair
    // was the cheapest of the eleven measured pairs that discriminate the two
    // modes, at 6.4 ms exhaustive; the same scan found pairs costing seconds.
    MCSComparison approximate(pair_of(SUCROSE, MACROLIDE));
    EXPECT_NEAR(approximate.Compare(0, 1), 0.605263, 1e-6);

    MCSOptions exhaustive;
    exhaustive.search_mode = MCSSearchMode::Exhaustive;
    MCSComparison exhaustive_comparison(pair_of(SUCROSE, MACROLIDE), exhaustive);
    EXPECT_NEAR(exhaustive_comparison.Compare(0, 1), 0.527778, 1e-6);
}

TEST_F(MCSComparisonTest, MaxMatchesBoundsTheSearch) {
    // Morphine against penicillin G, both 25 bonds. A budget of one match cuts
    // the count from 11 to 2: 2/48 against 11/39. Without this case only
    // max_matches == 0 is exercised, and a build that ignored every positive
    // value would pass the whole suite.
    MCSOptions budget;
    budget.max_matches = 1;
    MCSComparison bounded(pair_of(MORPHINE, PENICILLIN_G), budget);
    EXPECT_NEAR(bounded.Compare(0, 1), 0.958333, 1e-6);

    MCSComparison unbounded(pair_of(MORPHINE, PENICILLIN_G));
    EXPECT_NEAR(unbounded.Compare(0, 1), 0.717949, 1e-6);
}

TEST_F(MCSComparisonTest, HydrogenRepresentationDoesNotChangeTheScore) {
    // Explicit hydrogens are suppressed, so they cannot reach the denominator.
    auto explicit_benzene = from_smiles(BENZENE, "explicit");
    OEChem::OEAddExplicitHydrogens(*explicit_benzene);
    ASSERT_GT(explicit_benzene->NumBonds(), 6u);
    std::vector<std::shared_ptr<OEChem::OEMol>> with_hydrogens = {
        from_smiles(BENZENE, "plain"), explicit_benzene};
    MCSComparison hydrogen_comparison(with_hydrogens);
    EXPECT_NEAR(hydrogen_comparison.Compare(0, 1), 0.000000, 1e-6);

    // Isotopic hydrogens too. This one fails if retainIsotope is left at the
    // toolkit default of true: benzene-d1 would keep 7 bonds against benzene's
    // 6 and score 6/7 similarity, a distance of 0.142857.
    MCSComparison deuterium_comparison(pair_of(BENZENE_D1, BENZENE));
    EXPECT_NEAR(deuterium_comparison.Compare(0, 1), 0.000000, 1e-6);
}

// --- Reproducibility --------------------------------------------------------

TEST_F(MCSComparisonTest, ThreadCountDoesNotChangeTheResult) {
    // Reproducibility only. Equal results across thread counts cannot detect a
    // data race, and nothing here should be read as evidence of one's absence;
    // the isolation argument rests on Clone() giving each worker private
    // molecules, not on this test.
    MCSComparison comparison(mols_);
    EXPECT_EQ(run_pdist(comparison, 1), run_pdist(comparison, 8));
}

TEST_F(MCSComparisonTest, RepeatedPDistIsIdentical) {
    MCSComparison comparison(mols_);
    EXPECT_EQ(run_pdist(comparison, 4), run_pdist(comparison, 4));
}
