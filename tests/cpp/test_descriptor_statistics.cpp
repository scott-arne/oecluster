#include <algorithm>
#include <cmath>
#include <gtest/gtest.h>
#include <oechem.h>
#include "oecluster/DescriptorStatistics.h"
#include "oecluster/Error.h"

using namespace OECluster;

class DescriptorStatisticsTest : public ::testing::Test {
protected:
    void SetUp() override {
        auto add_mol = [this](const char* smi) {
            graph_mols_.emplace_back();
            OEChem::OESmilesToMol(graph_mols_.back(), smi);
        };
        add_mol("c1ccccc1");
        add_mol("c1ccc(O)cc1");
        add_mol("CCCCCCCC");
        add_mol("CC(=O)Oc1ccccc1C(=O)O");

        for (auto& gm : graph_mols_) {
            mols_.push_back(&static_cast<OEChem::OEMolBase&>(gm));
        }
    }

    std::vector<OEChem::OEGraphMol> graph_mols_;
    std::vector<OEChem::OEMolBase*> mols_;
};

TEST_F(DescriptorStatisticsTest, DefaultSourceIsOpenEye) {
    DescriptorStatisticsOptions options;
    const DescriptorStatisticsResult stats = descriptor_statistics(mols_, options);
    EXPECT_EQ(stats.num_rows, 4u);
    EXPECT_FALSE(stats.columns.empty());
    EXPECT_EQ(stats.mean.size(), stats.columns.size());
    EXPECT_EQ(stats.variance.size(), stats.columns.size());
    EXPECT_EQ(stats.minimum.size(), stats.columns.size());
    EXPECT_EQ(stats.maximum.size(), stats.columns.size());
    EXPECT_EQ(stats.present_count.size(), stats.columns.size());
}

TEST_F(DescriptorStatisticsTest, SurvivingColumnsAllHavePositiveVariance) {
    DescriptorStatisticsOptions options;
    const DescriptorStatisticsResult stats = descriptor_statistics(mols_, options);
    for (size_t i = 0; i < stats.columns.size(); ++i) {
        EXPECT_TRUE(std::isfinite(stats.variance[i])) << stats.columns[i];
        EXPECT_GT(stats.variance[i], 0.0) << stats.columns[i];
    }
}

TEST_F(DescriptorStatisticsTest, ZeroVarianceColumnsAreDroppedAndReported) {
    // The fixture above drops nothing: all eleven OpenEye columns have positive
    // variance over those four molecules, so a test written against mols_ would
    // assert 0 == 0 and pass against an implementation that never populates the
    // drop report at all. These four are each a single aromatic ring with no
    // rotatable bonds, which collapses exactly three columns.
    std::vector<OEChem::OEGraphMol> ring_mols(4);
    OEChem::OESmilesToMol(ring_mols[0], "c1ccccc1");
    OEChem::OESmilesToMol(ring_mols[1], "c1ccc(O)cc1");
    OEChem::OESmilesToMol(ring_mols[2], "Cc1ccccc1");
    OEChem::OESmilesToMol(ring_mols[3], "Nc1ccccc1");
    std::vector<OEChem::OEMolBase*> rings;
    for (auto& gm : ring_mols) {
        rings.push_back(&static_cast<OEChem::OEMolBase&>(gm));
    }

    DescriptorStatisticsOptions options;
    const DescriptorStatisticsResult stats = descriptor_statistics(rings, options);

    ASSERT_EQ(stats.dropped_columns.size(), stats.dropped_reasons.size());
    EXPECT_EQ(stats.dropped_columns.size(), 3u);
    EXPECT_EQ(stats.columns.size(), 8u);
    for (const char* name : {"HBA", "AromaticRingCount", "RotatableBondCount"}) {
        EXPECT_NE(std::find(stats.dropped_columns.begin(), stats.dropped_columns.end(), name),
                  stats.dropped_columns.end())
            << name;
        // A dropped column must be gone from the surviving vectors, not merely
        // named in the report.
        EXPECT_EQ(std::find(stats.columns.begin(), stats.columns.end(), name),
                  stats.columns.end())
            << name;
    }
    for (const std::string& reason : stats.dropped_reasons) {
        EXPECT_EQ(reason, "zero-variance") << reason;
    }
}

TEST_F(DescriptorStatisticsTest, ExplicitColumnsAreHonored) {
    DescriptorStatisticsOptions probe;
    const DescriptorStatisticsResult all = descriptor_statistics(mols_, probe);
    ASSERT_GE(all.columns.size(), 2u);

    DescriptorStatisticsOptions options;
    options.columns = {all.columns[0], all.columns[1]};
    const DescriptorStatisticsResult subset = descriptor_statistics(mols_, options);
    EXPECT_EQ(subset.columns.size(), 2u);
    EXPECT_EQ(subset.columns[0], all.columns[0]);
    EXPECT_EQ(subset.columns[1], all.columns[1]);
    EXPECT_NEAR(subset.variance[0], all.variance[0], 1e-9);
}

TEST_F(DescriptorStatisticsTest, InverseCovarianceIsSquareOverSurvivingColumns) {
    DescriptorStatisticsOptions options;
    options.inverse_covariance = true;
    const DescriptorStatisticsResult stats = descriptor_statistics(mols_, options);
    // Measured against OEFP 0.3.0: the OpenEye source has 11 columns, all
    // surviving on this fixture, and four rows give a rank-3 centred matrix.
    // Pinning 121 and 3 is what makes this bite -- k*k restates whatever k is,
    // and rank > 0 passes for 1, 2 or 3 alike.
    EXPECT_EQ(stats.columns.size(), 11u);
    EXPECT_EQ(stats.inverse_covariance.size(), 121u);
    EXPECT_EQ(stats.inverse_covariance_rank, 3u);
}

TEST_F(DescriptorStatisticsTest, InverseCovarianceIsSkippedByDefault) {
    DescriptorStatisticsOptions options;
    const DescriptorStatisticsResult stats = descriptor_statistics(mols_, options);
    EXPECT_TRUE(stats.inverse_covariance.empty());
    EXPECT_EQ(stats.inverse_covariance_rank, 0u);
}

TEST_F(DescriptorStatisticsTest, UnknownSourceIsRejected) {
    DescriptorStatisticsOptions options;
    options.sources = {"chemaxon"};
    EXPECT_THROW(descriptor_statistics(mols_, options), ComparisonError);
}

TEST_F(DescriptorStatisticsTest, UnknownColumnIsRejected) {
    DescriptorStatisticsOptions options;
    options.columns = {"NoSuchDescriptor"};
    EXPECT_THROW(descriptor_statistics(mols_, options), ComparisonError);
}

TEST_F(DescriptorStatisticsTest, UnknownGroupIsRejected) {
    DescriptorStatisticsOptions options;
    options.groups = {"rdkit:NoSuchGroup"};
    EXPECT_THROW(descriptor_statistics(mols_, options), ComparisonError);
}
