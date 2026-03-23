/**
 * @file test_cdist.cpp
 * @brief Tests for cross-distance computation.
 */

#include <gtest/gtest.h>
#include "oecluster/oecluster.h"
#include "oecluster/metrics/FingerprintMetric.h"
#include <oechem.h>
#include <vector>

using namespace OECluster;

class CDistTest : public ::testing::Test {
protected:
    void SetUp() override {
        auto add_mol = [this](const char* smi) {
            graph_mols_.emplace_back();
            OEChem::OESmilesToMol(graph_mols_.back(), smi);
        };
        // Set A: 2 molecules
        add_mol("c1ccccc1");      // benzene
        add_mol("c1ccc(O)cc1");   // phenol
        // Set B: 3 molecules
        add_mol("CCCCCCCC");      // octane
        add_mol("c1ccncc1");      // pyridine
        add_mol("CC(=O)O");       // acetic acid

        for (auto& gm : graph_mols_) {
            mols_.push_back(&static_cast<OEChem::OEMolBase&>(gm));
        }
    }

    std::vector<OEChem::OEGraphMol> graph_mols_;
    std::vector<OEChem::OEMolBase*> mols_;
};

TEST_F(CDistTest, OutputShape) {
    FingerprintMetric metric(mols_);
    size_t n_a = 2;
    size_t n_b = 3;
    std::vector<double> output(n_a * n_b, -1.0);
    cdist(metric, n_a, output.data());
    // All values should be filled (not -1.0)
    for (size_t k = 0; k < n_a * n_b; ++k) {
        EXPECT_GE(output[k], 0.0);
        EXPECT_LE(output[k], 1.0);
    }
}

TEST_F(CDistTest, MatchesDirectDistance) {
    FingerprintMetric metric(mols_);
    size_t n_a = 2;
    size_t n_b = 3;
    std::vector<double> output(n_a * n_b, 0.0);
    cdist(metric, n_a, output.data());

    // Verify output[i * n_b + j] == metric.Distance(i, n_a + j)
    for (size_t i = 0; i < n_a; ++i) {
        for (size_t j = 0; j < n_b; ++j) {
            double expected = metric.Distance(i, n_a + j);
            EXPECT_DOUBLE_EQ(output[i * n_b + j], expected);
        }
    }
}

TEST_F(CDistTest, WithCutoff) {
    FingerprintMetric metric(mols_);
    size_t n_a = 2;
    size_t n_b = 3;
    std::vector<double> output(n_a * n_b, -1.0);
    CDistOptions opts;
    opts.cutoff = 0.3;
    cdist(metric, n_a, output.data(), opts);
    for (size_t k = 0; k < n_a * n_b; ++k) {
        EXPECT_TRUE(output[k] == 0.0 || output[k] <= 0.3);
    }
}

TEST_F(CDistTest, ProgressCallback) {
    FingerprintMetric metric(mols_);
    size_t n_a = 2;
    size_t n_b = 3;
    std::vector<double> output(n_a * n_b, 0.0);
    CDistOptions opts;
    size_t last_total = 0;
    opts.progress = [&](size_t completed, size_t total) {
        last_total = total;
    };
    cdist(metric, n_a, output.data(), opts);
    EXPECT_EQ(last_total, n_a * n_b);
}
