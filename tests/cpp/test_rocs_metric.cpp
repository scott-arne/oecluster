#include <gtest/gtest.h>
#include "oecluster/oecluster.h"
#include "oecluster/metrics/ROCSMetric.h"
#include <oechem.h>
#include <oeshape.h>

using namespace OECluster;

class ROCSMetricTest : public ::testing::Test {
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

TEST_F(ROCSMetricTest, ConstructAndSize) {
    ROCSMetric metric(mols_);
    EXPECT_EQ(metric.Size(), 3);
    EXPECT_EQ(metric.Name(), "rocs");
}

TEST_F(ROCSMetricTest, DistanceRange) {
    ROCSMetric metric(mols_);
    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = i + 1; j < 3; ++j) {
            double d = metric.Distance(i, j);
            // Default is ComboNorm: distance ∈ [0,1]
            EXPECT_GE(d, 0.0);
            EXPECT_LE(d, 1.0);
        }
    }
}

TEST_F(ROCSMetricTest, ShapeOnlyDistanceRange) {
    ROCSOptions opts;
    opts.score_type = ROCSScoreType::Shape;
    ROCSMetric metric(mols_, opts);
    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = i + 1; j < 3; ++j) {
            double d = metric.Distance(i, j);
            // Shape Tanimoto is in [0,1], so distance is in [0, 1]
            EXPECT_GE(d, 0.0);
            EXPECT_LE(d, 1.0);
        }
    }
}

TEST_F(ROCSMetricTest, CloneCreatesIndependentScorer) {
    ROCSMetric metric(mols_);
    auto clone = metric.Clone();
    EXPECT_EQ(clone->Size(), 3);
    double d = clone->Distance(0, 1);
    // Default is ComboNorm: distance ∈ [0,1]
    EXPECT_GE(d, 0.0);
    EXPECT_LE(d, 1.0);
}

TEST_F(ROCSMetricTest, IntegrationWithPDist) {
    ROCSMetric metric(mols_);
    DenseStorage storage(3);
    pdist(metric, storage);
    EXPECT_EQ(storage.NumPairs(), 3);
}

TEST_F(ROCSMetricTest, ComboNormScoreType) {
    ROCSOptions opts;
    opts.score_type = ROCSScoreType::ComboNorm;
    ROCSMetric metric(mols_, opts);
    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = i + 1; j < 3; ++j) {
            double d = metric.Distance(i, j);
            // ComboNorm: distance = 1.0 - combo/2.0 ∈ [0,1]
            EXPECT_GE(d, 0.0);
            EXPECT_LE(d, 1.0);
        }
    }
}

TEST_F(ROCSMetricTest, ComboScoreType) {
    ROCSOptions opts;
    opts.score_type = ROCSScoreType::Combo;
    ROCSMetric metric(mols_, opts);
    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = i + 1; j < 3; ++j) {
            double d = metric.Distance(i, j);
            // Combo: distance = 2.0 - combo ∈ [0,2]
            EXPECT_GE(d, 0.0);
            EXPECT_LE(d, 2.0);
        }
    }
}

TEST_F(ROCSMetricTest, ColorScoreType) {
    ROCSOptions opts;
    opts.score_type = ROCSScoreType::Color;
    ROCSMetric metric(mols_, opts);
    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = i + 1; j < 3; ++j) {
            double d = metric.Distance(i, j);
            // Color: distance = 1.0 - color ∈ [0,1]
            EXPECT_GE(d, 0.0);
            EXPECT_LE(d, 1.0);
        }
    }
}

TEST_F(ROCSMetricTest, ColorForceFieldConfiguration) {
    ROCSOptions opts;
    opts.score_type = ROCSScoreType::Color;
    opts.color_ff_type = 2;  // ExplicitMillsDean
    ROCSMetric metric(mols_, opts);
    EXPECT_EQ(metric.Size(), 3);
    // Just verify it constructs and runs without error
    double d = metric.Distance(0, 1);
    EXPECT_GE(d, 0.0);
    EXPECT_LE(d, 1.0);
}
