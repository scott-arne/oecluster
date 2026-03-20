#include <gtest/gtest.h>
#include "oecluster/oecluster.h"

#ifdef OECLUSTER_HAS_GRAPHSIM
#include "oecluster/metrics/FingerprintMetric.h"
#include <oechem.h>
#include <oegraphsim.h>

using namespace OECluster;

class FingerprintMetricTest : public ::testing::Test {
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

TEST_F(FingerprintMetricTest, ConstructAndSize) {
    FingerprintMetric metric(mols_);
    EXPECT_EQ(metric.Size(), 3);
    EXPECT_EQ(metric.Name(), "fingerprint");
}

TEST_F(FingerprintMetricTest, SimilarMoleculesCloser) {
    FingerprintMetric metric(mols_);
    double d_benzene_phenol = metric.Distance(0, 1);
    double d_benzene_octane = metric.Distance(0, 2);
    EXPECT_LT(d_benzene_phenol, d_benzene_octane);
}

TEST_F(FingerprintMetricTest, DistanceRange) {
    FingerprintMetric metric(mols_);
    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = i + 1; j < 3; ++j) {
            double d = metric.Distance(i, j);
            EXPECT_GE(d, 0.0);
            EXPECT_LE(d, 1.0);
        }
    }
}

TEST_F(FingerprintMetricTest, CloneSharesData) {
    FingerprintMetric metric(mols_);
    auto clone = metric.Clone();
    EXPECT_EQ(clone->Size(), 3);
    EXPECT_DOUBLE_EQ(clone->Distance(0, 1), metric.Distance(0, 1));
}

TEST_F(FingerprintMetricTest, IntegrationWithPDist) {
    FingerprintMetric metric(mols_);
    DenseStorage storage(3);
    pdist(metric, storage);
    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = i + 1; j < 3; ++j) {
            EXPECT_GT(storage.Get(i, j), 0.0);
        }
    }
}

#endif  // OECLUSTER_HAS_GRAPHSIM
