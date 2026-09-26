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
