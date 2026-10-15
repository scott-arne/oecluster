#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <gtest/gtest.h>
#include <oechem.h>
#include "oecluster/DescriptorStatistics.h"
#include "oecluster/Error.h"
#include "../../src/DescriptorBuild.h"

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
    // Matched against the refusal the guard pre-empts, so removing the guard
    // fails here. That refusal used to report every column as zero-variance;
    // it now names the causes present, and under two molecules every column
    // trips present_count < 2, so the cause would be insufficient-data. Both
    // refusals throw ComparisonError, so match on the count, not the type.
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

TEST(VarianceDropReasonTest, DistinguishesTheThreeCauses) {
    // Ordered, because -inf is both non-finite and negative: the predicates
    // overlap unless present_count is asked about first and finiteness
    // second. An earlier draft of this rule had them overlapping.
    EXPECT_EQ(variance_drop_reason(std::numeric_limits<double>::quiet_NaN(), 0u),
              "insufficient-data");
    EXPECT_EQ(variance_drop_reason(std::numeric_limits<double>::quiet_NaN(), 1u),
              "insufficient-data");
    EXPECT_EQ(variance_drop_reason(std::numeric_limits<double>::quiet_NaN(), 2u),
              "non-finite-variance");
    EXPECT_EQ(variance_drop_reason(-std::numeric_limits<double>::infinity(), 2u),
              "non-finite-variance");
    EXPECT_EQ(variance_drop_reason(0.0, 2u), "zero-variance");
    EXPECT_EQ(variance_drop_reason(-1.0, 2u), "zero-variance");
    EXPECT_TRUE(variance_drop_reason(1.5, 2u).empty());
}

TEST(DescriptorStatisticsNonFiniteTest, NonFiniteVarianceColumnsAreDropped) {
    // Same present-NaN route as the comparison test: this pins the
    // non-finite-variance reason on the statistics surface, which until now
    // was only covered for variance == 0.
    std::vector<OEChem::OEGraphMol> graph_mols(2);
    OEChem::OESmilesToMol(graph_mols[0], "c1ccccc1");
    OEChem::OESmilesToMol(graph_mols[1], "[Na+]");
    std::vector<OEChem::OEMolBase*> mols;
    for (auto& gm : graph_mols) {
        mols.push_back(&static_cast<OEChem::OEMolBase&>(gm));
    }

    DescriptorStatisticsOptions opts;
    opts.sources = {"rdkit"};
    opts.columns = {"BCUT2D_MWHI", "BCUT2D_LOGPLOW", "MaxPartialCharge"};
    const DescriptorStatisticsResult result = descriptor_statistics(mols, opts);

    ASSERT_EQ(result.columns.size(), 1u);
    EXPECT_EQ(result.columns[0], "MaxPartialCharge");
    EXPECT_TRUE(std::isfinite(result.variance[0]));
    EXPECT_EQ(result.dropped_columns.size(), 2u);
    for (const std::string& reason : result.dropped_reasons) {
        // Both ions are present -- RDKit marks the value valid and writes a
        // NaN -- so present_count is 2 and the cause the layer can establish
        // is the variance, not a shortage of rows.
        EXPECT_EQ(reason, "non-finite-variance") << reason;
    }
}

TEST(DescriptorStatisticsEmptySelectionTest, TheAllDroppedRefusalNamesTheCauses) {
    // The refusal fires before any result is returned, so it is the whole of
    // what the caller is told -- the per-column reasons never reach them. It
    // used to report that every column was constant, and over these ions,
    // where XLogP is computable exactly once, constant is the one thing the
    // column provably is not.
    std::vector<OEChem::OEGraphMol> ion_mols(4);
    OEChem::OESmilesToMol(ion_mols[0], "[Na+]");
    OEChem::OESmilesToMol(ion_mols[1], "[K+]");
    OEChem::OESmilesToMol(ion_mols[2], "[Li+]");
    OEChem::OESmilesToMol(ion_mols[3], "CCO");
    std::vector<OEChem::OEMolBase*> ions;
    for (auto& gm : ion_mols) {
        ions.push_back(&static_cast<OEChem::OEMolBase&>(gm));
    }

    DescriptorStatisticsOptions logp_only;
    logp_only.columns = {"XLogP"};
    std::string message;
    try {
        descriptor_statistics(ions, logp_only);
        ADD_FAILURE() << "an empty surviving selection must be refused";
    } catch (const ComparisonError& err) {
        message = err.what();
    }
    EXPECT_NE(message.find("Every selected descriptor column was dropped over these "
                           "molecules (insufficient-data)"),
              std::string::npos)
        << message;
    EXPECT_EQ(message.find("zero"), std::string::npos) << message;

    // The control, because a reason that always appears carries no more
    // information than the wording it replaced: two identical molecules make
    // the same column genuinely constant, and that refusal must still say so.
    std::vector<OEChem::OEGraphMol> identical_mols(2);
    OEChem::OESmilesToMol(identical_mols[0], "CCO");
    OEChem::OESmilesToMol(identical_mols[1], "CCO");
    std::vector<OEChem::OEMolBase*> identical;
    for (auto& gm : identical_mols) {
        identical.push_back(&static_cast<OEChem::OEMolBase&>(gm));
    }

    std::string constant_message;
    try {
        descriptor_statistics(identical, logp_only);
        ADD_FAILURE() << "an empty surviving selection must be refused";
    } catch (const ComparisonError& err) {
        constant_message = err.what();
    }
    EXPECT_NE(constant_message.find("Every selected descriptor column was dropped over these "
                                    "molecules (zero-variance)"),
              std::string::npos)
        << constant_message;
}

TEST(JoinDropReasonsTest, DeduplicatesAndOrdersCanonically) {
    // The joined list is part of two public refusal messages, so it is
    // ordered here rather than by however the reasons were collected. The
    // ledgers append every non-numeric drop before any variance drop, and a
    // future collector could append in some third order; sorting at the join
    // means none of that reaches the message. Caller column order is not the
    // reason -- resolve_column_indices sorts the selection long before this.
    EXPECT_EQ(join_drop_reasons({}), "");
    EXPECT_EQ(join_drop_reasons({"zero-variance"}), "zero-variance");
    EXPECT_EQ(join_drop_reasons({"zero-variance", "zero-variance"}), "zero-variance");
    EXPECT_EQ(join_drop_reasons({"zero-variance", "insufficient-data"}),
              "insufficient-data, zero-variance");
    EXPECT_EQ(join_drop_reasons({"zero-variance", "non-finite-variance", "insufficient-data"}),
              "insufficient-data, non-finite-variance, zero-variance");
    // The shapes below are simulated, not collected. No source
    // make_descriptor_calculator can build emits a non-numeric column, so no
    // schema produces a selection carrying this reason and these inputs are
    // handed to the helper directly.
    //
    // What they simulate is the mixed selection a future source would make:
    // both refusals are built from the surface's whole drop ledger, which the
    // non-numeric filter writes to before the variance filter does, so the
    // helper would receive "non-numeric" alongside a variance reason whenever
    // some columns were dropped as non-numeric and the remaining numeric ones
    // all failed the variance test. The refusal claims to say why every
    // selected column was dropped, so the reason would have to reach the
    // message rather than be silently discarded for sitting outside the
    // canonical order. It is placed after the causes that are in it.
    //
    // Pinned here because the filter the reason comes from is kept against
    // schema drift, and a message that over-claims the moment it first fires
    // is the defect this helper exists to prevent.
    EXPECT_EQ(join_drop_reasons({"non-numeric"}), "non-numeric");
    EXPECT_EQ(join_drop_reasons({"non-numeric", "zero-variance"}),
              "zero-variance, non-numeric");
    EXPECT_EQ(join_drop_reasons({"non-numeric", "insufficient-data", "non-numeric"}),
              "insufficient-data, non-numeric");
}
