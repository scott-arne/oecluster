#include <gtest/gtest.h>
#include "oecluster/oecluster.h"
#include "oecluster/comparisons/FingerprintComparison.h"
#include <oechem.h>

using namespace OECluster;

class FingerprintComparisonTest : public ::testing::Test {
protected:
    void SetUp() override {
        auto add_mol = [this](const char* smi) {
            graph_mols_.emplace_back();
            OEChem::OESmilesToMol(graph_mols_.back(), smi);
        };
        add_mol("c1ccccc1");      // benzene
        add_mol("c1ccc(O)cc1");   // phenol
        add_mol("CCCCCCCC");      // octane

        for (auto& gm : graph_mols_) {
            mols_.push_back(&static_cast<OEChem::OEMolBase&>(gm));
        }
    }

    std::vector<OEChem::OEGraphMol> graph_mols_;
    std::vector<OEChem::OEMolBase*> mols_;
};

TEST_F(FingerprintComparisonTest, ConstructAndSize) {
    FingerprintOptions opts;
    EXPECT_EQ(opts.fp_type, "morgan");

    FingerprintComparison comparison(mols_);
    EXPECT_EQ(comparison.Size(), 3);
    EXPECT_EQ(comparison.ComparisonName(), "fingerprint");
}

TEST_F(FingerprintComparisonTest, SimilarMoleculesCloser) {
    FingerprintComparison comparison(mols_);
    double d_benzene_phenol = comparison.Compare(0, 1);
    double d_benzene_octane = comparison.Compare(0, 2);
    EXPECT_LT(d_benzene_phenol, d_benzene_octane);
}

TEST_F(FingerprintComparisonTest, DistanceRange) {
    FingerprintComparison comparison(mols_);
    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = i + 1; j < 3; ++j) {
            double d = comparison.Compare(i, j);
            EXPECT_GE(d, 0.0);
            EXPECT_LE(d, 1.0);
        }
    }
}

TEST_F(FingerprintComparisonTest, CloneSharesData) {
    FingerprintComparison comparison(mols_);
    auto clone = comparison.Clone();
    EXPECT_EQ(clone->Size(), 3);
    EXPECT_DOUBLE_EQ(clone->Compare(0, 1), comparison.Compare(0, 1));
}

TEST_F(FingerprintComparisonTest, IntegrationWithPDist) {
    FingerprintComparison comparison(mols_);
    DenseStorage storage(3);
    pdist(comparison, storage);
    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = i + 1; j < 3; ++j) {
            EXPECT_GT(storage.Get(i, j), 0.0);
        }
    }
}

TEST_F(FingerprintComparisonTest, MorganFingerprintType) {
    FingerprintOptions opts;
    opts.fp_type = "morgan";
    opts.numbits = 2048;
    opts.max_distance = 2;
    FingerprintComparison comparison(mols_, opts);
    EXPECT_EQ(comparison.Size(), 3);
    double d = comparison.Compare(0, 1);
    EXPECT_GE(d, 0.0);
    EXPECT_LE(d, 1.0);
}

TEST_F(FingerprintComparisonTest, AtomPairFingerprintType) {
    FingerprintOptions opts;
    opts.fp_type = "atom_pair";
    FingerprintComparison comparison(mols_, opts);
    double d = comparison.Compare(0, 1);
    EXPECT_GE(d, 0.0);
    EXPECT_LE(d, 1.0);
}

TEST_F(FingerprintComparisonTest, TanimotoSimilarityComplementsJaccardDistance) {
    FingerprintOptions opts;
    opts.similarity = false;
    FingerprintComparison distance_comparison(mols_, opts);

    opts.similarity = true;
    FingerprintComparison similarity_comparison(mols_, opts);

    EXPECT_NEAR(
        distance_comparison.Compare(0, 1) + similarity_comparison.Compare(0, 1),
        1.0,
        1.0e-12);
}

TEST_F(FingerprintComparisonTest, DiceMetric) {
    FingerprintOptions opts;
    opts.metric = "dice";
    FingerprintComparison comparison(mols_, opts);
    double d = comparison.Compare(0, 1);
    EXPECT_GE(d, 0.0);
    EXPECT_LE(d, 1.0);
}

TEST_F(FingerprintComparisonTest, RemovedOpenEyeFingerprintTypesThrow) {
    FingerprintOptions opts;

    for (const auto* fp_type : {"circular", "tree", "path", "maccs", "lingo"}) {
        opts.fp_type = fp_type;
        EXPECT_THROW(FingerprintComparison(mols_, opts), ComparisonError)
            << "fp_type=" << fp_type;
    }
}

TEST_F(FingerprintComparisonTest, EuclideanSimilarityFunctionThrows) {
    FingerprintOptions opts;
    opts.metric = "euclidean";
    EXPECT_THROW(FingerprintComparison(mols_, opts), ComparisonError);
}

TEST_F(FingerprintComparisonTest, InvalidFpTypeThrows) {
    FingerprintOptions opts;
    opts.fp_type = "invalid";
    EXPECT_THROW(FingerprintComparison(mols_, opts), ComparisonError);
}
