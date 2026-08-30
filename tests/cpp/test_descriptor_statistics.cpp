#include <algorithm>
#include <cmath>
#include <string>
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

TEST_F(DescriptorStatisticsTest, TooFewMoleculesIsRejectedByCount) {
    // The pre-guard behavior was to report every column as zero-variance, which
    // is a different and false claim. Match on the count so a regression back to
    // that message fails here rather than passing on the exception type alone.
    for (size_t n : {size_t{0}, size_t{1}}) {
        const std::vector<OEChem::OEMolBase*> few(mols_.begin(), mols_.begin() + n);
        DescriptorStatisticsOptions options;
        try {
            descriptor_statistics(few, options);
            ADD_FAILURE() << "expected a throw for " << n << " molecules";
        } catch (const ComparisonError& err) {
            EXPECT_NE(std::string(err.what()).find("at least two molecules"),
                      std::string::npos)
                << n << ": " << err.what();
            EXPECT_NE(std::string(err.what()).find("got " + std::to_string(n)),
                      std::string::npos)
                << n << ": " << err.what();
        }
    }
}

TEST_F(DescriptorStatisticsTest, TwoMoleculesAreEnoughToFit) {
    // The floor is exactly two, not "several". Without this, raising the guard
    // to three would go unnoticed.
    const std::vector<OEChem::OEMolBase*> pair(mols_.begin(), mols_.begin() + 2);
    DescriptorStatisticsOptions options;
    const DescriptorStatisticsResult stats = descriptor_statistics(pair, options);
    EXPECT_EQ(stats.num_rows, 2u);
    EXPECT_FALSE(stats.columns.empty());
}

TEST_F(DescriptorStatisticsTest, StatisticPayloadMatchesKnownValues) {
    // Shape assertions cannot catch a transposition among the six parallel
    // vectors. HeavyAtomCount's four statistics are mutually distinct
    // (6, 13, 8.5, 29/3), so any pair of them swapped fails here.
    DescriptorStatisticsOptions options;
    options.columns = {"HeavyAtomCount", "TotalAtomCount", "LipinskiHBD"};
    const DescriptorStatisticsResult stats = descriptor_statistics(mols_, options);

    ASSERT_EQ(stats.columns.size(), 3u);
    EXPECT_EQ(stats.columns[0], "HeavyAtomCount");
    EXPECT_EQ(stats.columns[1], "TotalAtomCount");
    EXPECT_EQ(stats.columns[2], "LipinskiHBD");

    EXPECT_NEAR(stats.mean[0], 8.5, 1e-9);
    EXPECT_NEAR(stats.variance[0], 29.0 / 3.0, 1e-9);
    EXPECT_NEAR(stats.minimum[0], 6.0, 1e-9);
    EXPECT_NEAR(stats.maximum[0], 13.0, 1e-9);

    EXPECT_NEAR(stats.mean[1], 18.0, 1e-9);
    EXPECT_NEAR(stats.variance[1], 134.0 / 3.0, 1e-9);
    EXPECT_NEAR(stats.minimum[1], 12.0, 1e-9);
    EXPECT_NEAR(stats.maximum[1], 26.0, 1e-9);

    EXPECT_NEAR(stats.mean[2], 0.5, 1e-9);
    EXPECT_NEAR(stats.variance[2], 1.0 / 3.0, 1e-9);
    EXPECT_NEAR(stats.minimum[2], 0.0, 1e-9);
    EXPECT_NEAR(stats.maximum[2], 1.0, 1e-9);

    for (size_t k = 0; k < stats.columns.size(); ++k) {
        EXPECT_EQ(stats.present_count[k], 4u) << stats.columns[k];
    }
}

TEST_F(DescriptorStatisticsTest, ResultsComeBackInSchemaOrderNotRequestOrder) {
    // TotalAtomCount is schema index 3 and HeavyAtomCount is index 2, so
    // requesting them in this order proves the result is re-sorted rather than
    // echoing the request, and that the statistics follow the returned names.
    DescriptorStatisticsOptions options;
    options.columns = {"TotalAtomCount", "HeavyAtomCount"};
    const DescriptorStatisticsResult stats = descriptor_statistics(mols_, options);

    ASSERT_EQ(stats.columns.size(), 2u);
    EXPECT_EQ(stats.columns[0], "HeavyAtomCount");
    EXPECT_EQ(stats.columns[1], "TotalAtomCount");
    EXPECT_NEAR(stats.mean[0], 8.5, 1e-9);
    EXPECT_NEAR(stats.mean[1], 18.0, 1e-9);
}

TEST_F(DescriptorStatisticsTest, InverseCovarianceReportsItsRowCount) {
    DescriptorStatisticsOptions options;
    options.inverse_covariance = true;
    const DescriptorStatisticsResult stats = descriptor_statistics(mols_, options);
    EXPECT_EQ(stats.inverse_covariance_rows, 4u);
    EXPECT_EQ(stats.num_rows, 4u);

    DescriptorStatisticsOptions without;
    const DescriptorStatisticsResult skipped = descriptor_statistics(mols_, without);
    EXPECT_EQ(skipped.inverse_covariance_rows, 0u);
}
