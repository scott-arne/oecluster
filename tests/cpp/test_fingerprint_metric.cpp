#include <gtest/gtest.h>
#include "oecluster/oecluster.h"
#include "oecluster/metrics/FingerprintMetric.h"
#include <oechem.h>

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
    FingerprintOptions opts;
    EXPECT_EQ(opts.fp_type, "morgan");

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

TEST_F(FingerprintMetricTest, MorganFingerprintType) {
    FingerprintOptions opts;
    opts.fp_type = "morgan";
    opts.numbits = 2048;
    opts.max_distance = 2;
    FingerprintMetric metric(mols_, opts);
    EXPECT_EQ(metric.Size(), 3);
    double d = metric.Distance(0, 1);
    EXPECT_GE(d, 0.0);
    EXPECT_LE(d, 1.0);
}

TEST_F(FingerprintMetricTest, AtomPairFingerprintType) {
    FingerprintOptions opts;
    opts.fp_type = "atom_pair";
    FingerprintMetric metric(mols_, opts);
    double d = metric.Distance(0, 1);
    EXPECT_GE(d, 0.0);
    EXPECT_LE(d, 1.0);
}

TEST_F(FingerprintMetricTest, TanimotoSimilarityComplementsJaccardDistance) {
    FingerprintOptions opts;
    opts.similarity = false;
    FingerprintMetric distance_metric(mols_, opts);

    opts.similarity = true;
    FingerprintMetric similarity_metric(mols_, opts);

    EXPECT_NEAR(
        distance_metric.Distance(0, 1) + similarity_metric.Distance(0, 1),
        1.0,
        1.0e-12);
}

TEST_F(FingerprintMetricTest, DiceSimilarity) {
    FingerprintOptions opts;
    opts.similarity_func = "dice";
    FingerprintMetric metric(mols_, opts);
    double d = metric.Distance(0, 1);
    EXPECT_GE(d, 0.0);
    EXPECT_LE(d, 1.0);
}

TEST_F(FingerprintMetricTest, RemovedOpenEyeFingerprintTypesThrow) {
    FingerprintOptions opts;

    for (const auto* fp_type : {"circular", "tree", "path", "maccs", "lingo"}) {
        opts.fp_type = fp_type;
        EXPECT_THROW(FingerprintMetric(mols_, opts), MetricError)
            << "fp_type=" << fp_type;
    }
}

TEST_F(FingerprintMetricTest, EuclideanSimilarityFunctionThrows) {
    FingerprintOptions opts;
    opts.similarity_func = "euclidean";
    EXPECT_THROW(FingerprintMetric(mols_, opts), MetricError);
}

TEST_F(FingerprintMetricTest, InvalidFpTypeThrows) {
    FingerprintOptions opts;
    opts.fp_type = "invalid";
    EXPECT_THROW(FingerprintMetric(mols_, opts), MetricError);
}
