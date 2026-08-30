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
