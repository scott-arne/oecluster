#include <algorithm>
#include <cmath>
#include <iterator>
#include <string>
#include <utility>
#include <vector>
#include <gtest/gtest.h>
#include <oechem.h>
#include "oecluster/CDist.h"
#include "oecluster/Error.h"
#include "oecluster/PDist.h"
#include "oecluster/StorageBackend.h"
#include "oecluster/comparisons/DescriptorComparison.h"

using namespace OECluster;

class DescriptorComparisonTest : public ::testing::Test {
protected:
    void SetUp() override {
        const char* smiles[] = {"c1ccccc1", "c1ccc(O)cc1", "CCCCCCCC",
                                "CC(=O)Oc1ccccc1C(=O)O", "CCN(CC)CCOC(=O)c1ccccc1"};
        for (const char* smi : smiles) {
            graph_mols_.emplace_back();
            OEChem::OESmilesToMol(graph_mols_.back(), smi);
        }
        for (auto& gm : graph_mols_) {
            mols_.push_back(&static_cast<OEChem::OEMolBase&>(gm));
        }
    }

    std::vector<OEChem::OEGraphMol> graph_mols_;
    std::vector<OEChem::OEMolBase*> mols_;
};

TEST_F(DescriptorComparisonTest, DefaultsToStandardizedEuclidean) {
    DescriptorComparison comparison(mols_);
    EXPECT_EQ(comparison.ComparisonName(), "descriptor");
    EXPECT_EQ(comparison.Size(), mols_.size());
    EXPECT_FALSE(comparison.Columns().empty());
    EXPECT_EQ(comparison.Variances().size(), comparison.Columns().size());
    EXPECT_TRUE(comparison.InverseCovariance().empty());
}

TEST_F(DescriptorComparisonTest, SelfDistanceIsZero) {
    DescriptorComparison comparison(mols_);
    for (size_t i = 0; i < mols_.size(); ++i) {
        EXPECT_NEAR(comparison.Compare(i, i), 0.0, 1e-12);
    }
}

TEST_F(DescriptorComparisonTest, PDistAgreesWithCompare) {
    DescriptorComparison comparison(mols_);
    DenseStorage storage(mols_.size());
    ASSERT_TRUE(comparison.TryPDist(storage, PDistOptions()));
    for (size_t i = 0; i < mols_.size(); ++i) {
        for (size_t j = i + 1; j < mols_.size(); ++j) {
            EXPECT_NEAR(storage.Get(i, j), comparison.Compare(i, j), 1e-9);
        }
    }
}

TEST_F(DescriptorComparisonTest, CDistAgreesWithCompare) {
    DescriptorComparison comparison(mols_);
    const size_t n_a = 2;
    const size_t n_b = mols_.size() - n_a;
    std::vector<double> output(n_a * n_b, -1.0);
    ASSERT_TRUE(comparison.TryCDist(n_a, output.data(), CDistOptions()));
    for (size_t i = 0; i < n_a; ++i) {
        for (size_t j = 0; j < n_b; ++j) {
            EXPECT_NEAR(output[i * n_b + j], comparison.Compare(i, n_a + j), 1e-9);
        }
    }
}

TEST_F(DescriptorComparisonTest, DefaultFactsAreCompleteAndMetric) {
    DescriptorComparison comparison(mols_);
    const GateFacts facts = comparison.Facts();
    EXPECT_EQ(facts.is_distance, Capability::Yes);
    EXPECT_EQ(facts.zero_self, Capability::Yes);
    EXPECT_EQ(facts.triangle, Capability::Yes);
    EXPECT_EQ(facts.data_integrity, DataIntegrity::Complete);
}

TEST_F(DescriptorComparisonTest, OverflowingDistancesDowngradeTheIntegrityStamp) {
    // The one case the input mask cannot predict: every descriptor is present
    // and finite, but the accumulator overflows while scoring. A large
    // Minkowski exponent reaches it deterministically, since raw descriptor
    // values differ by far more than 1. This is the reason Facts() must be
    // read after the computation rather than at construction.
    DescriptorOptions opts;
    opts.metric = "minkowski";
    opts.p = 400.0;
    DescriptorComparison comparison(mols_, opts);
    EXPECT_EQ(comparison.Facts().data_integrity, DataIntegrity::Complete);

    DenseStorage storage(mols_.size());
    ASSERT_TRUE(comparison.TryPDist(storage, PDistOptions()));
    EXPECT_FALSE(std::isfinite(storage.Get(0, 4)));
    EXPECT_EQ(comparison.Facts().data_integrity, DataIntegrity::NaNPresent);
}

TEST_F(DescriptorComparisonTest, PropagateStampsNaNPresent) {
    DescriptorOptions opts;
    opts.metric = "euclidean";
    opts.missing = "propagate";
    DescriptorComparison comparison(mols_, opts);
    EXPECT_EQ(comparison.Facts().data_integrity, DataIntegrity::NaNPresent);
}

TEST_F(DescriptorComparisonTest, IgnoreStampsSubsetScored) {
    DescriptorOptions opts;
    opts.metric = "euclidean";
    opts.missing = "ignore";
    DescriptorComparison comparison(mols_, opts);
    EXPECT_EQ(comparison.Facts().data_integrity, DataIntegrity::SubsetScored);
}

TEST_F(DescriptorComparisonTest, IgnoreIsRejectedForFittedMetrics) {
    DescriptorOptions opts;
    opts.missing = "ignore";  // metric defaults to standardized_euclidean
    EXPECT_THROW(DescriptorComparison(mols_, opts), ComparisonError);

    DescriptorOptions mahalanobis;
    mahalanobis.metric = "mahalanobis";
    mahalanobis.missing = "ignore";
    EXPECT_THROW(DescriptorComparison(mols_, mahalanobis), ComparisonError);
}

TEST_F(DescriptorComparisonTest, UnknownMissingPolicyIsRejected) {
    DescriptorOptions opts;
    opts.missing = "drop";
    EXPECT_THROW(DescriptorComparison(mols_, opts), ComparisonError);
}

TEST_F(DescriptorComparisonTest, BooleanMetricIsRejectedOnDescriptors) {
    DescriptorOptions opts;
    opts.metric = "tanimoto";
    EXPECT_THROW(DescriptorComparison(mols_, opts), ComparisonError);
}

TEST_F(DescriptorComparisonTest, BothOverridesTogetherAreRejected) {
    DescriptorOptions opts;
    opts.variances = {1.0, 1.0};
    opts.inverse_covariance = {1.0, 0.0, 0.0, 1.0};
    EXPECT_THROW(DescriptorComparison(mols_, opts), ComparisonError);
}

TEST_F(DescriptorComparisonTest, OverrideIsRejectedForTheWrongMetric) {
    DescriptorOptions with_variances;
    with_variances.metric = "mahalanobis";
    with_variances.variances = {1.0};
    EXPECT_THROW(DescriptorComparison(mols_, with_variances), ComparisonError);

    DescriptorOptions with_inverse;
    with_inverse.metric = "euclidean";
    with_inverse.inverse_covariance = {1.0};
    EXPECT_THROW(DescriptorComparison(mols_, with_inverse), ComparisonError);
}

TEST_F(DescriptorComparisonTest, OverrideLengthMismatchIsRejected) {
    DescriptorOptions opts;
    opts.variances = {1.0};
    EXPECT_THROW(DescriptorComparison(mols_, opts), ComparisonError);
}

TEST_F(DescriptorComparisonTest, OverridePathDropsNothingAndReproducesTheFit) {
    DescriptorComparison fitted(mols_);
    DescriptorOptions opts;
    opts.columns = fitted.Columns();
    opts.variances = fitted.Variances();
    DescriptorComparison overridden(mols_, opts);

    EXPECT_TRUE(overridden.DroppedColumns().empty());
    EXPECT_EQ(overridden.Columns(), fitted.Columns());
    EXPECT_NEAR(overridden.Compare(0, 1), fitted.Compare(0, 1), 1e-9);
}

TEST_F(DescriptorComparisonTest, DropReportIsParallel) {
    DescriptorComparison comparison(mols_);
    EXPECT_EQ(comparison.DroppedColumns().size(), comparison.DroppedReasons().size());
}

TEST_F(DescriptorComparisonTest, CloneScoresIdentically) {
    DescriptorComparison comparison(mols_);
    const std::unique_ptr<PairwiseComparison> clone = comparison.Clone();
    EXPECT_EQ(clone->Size(), comparison.Size());
    EXPECT_NEAR(clone->Compare(0, 2), comparison.Compare(0, 2), 1e-12);
}

TEST_F(DescriptorComparisonTest, ExcludedIndicesIsEmptyForCompleteInput) {
    const std::vector<size_t> excluded = descriptor_excluded_indices(mols_, DescriptorOptions());
    EXPECT_TRUE(excluded.empty());
}

TEST_F(DescriptorComparisonTest, NullMoleculeIsRejected) {
    std::vector<OEChem::OEMolBase*> with_null = mols_;
    with_null[1] = nullptr;
    EXPECT_THROW(DescriptorComparison(with_null, DescriptorOptions()), ComparisonError);
}

// Spec section 2.4's measured missingness set: 19 molecules, of which OpenEye
// cannot assign XLogP types to 6. The empty-mask test above proves the mask is
// not over-eager; only a set with real gaps proves it fires at all. Verified
// against OEFP 0.3.0 and OpenEye 2026.1.0: the 11-column OpenEye default has
// two gapped columns, XLogP and FractionCsp3, and both are gapped on exactly
// the same 6 rows, leaving 13 complete molecules.
class DescriptorMissingnessTest : public ::testing::Test {
protected:
    void SetUp() override {
        const std::pair<const char*, const char*> molecules[] = {
            {"benzene", "c1ccccc1"},
            {"phenol", "c1ccc(O)cc1"},
            {"octane", "CCCCCCCC"},
            {"aspirin", "CC(=O)Oc1ccccc1C(=O)O"},
            {"procaine", "CCN(CC)CCOC(=O)c1ccccc1"},
            {"pyridine", "c1ccncc1"},
            {"ethanol", "CCO"},
            {"toluene", "Cc1ccccc1"},
            {"aniline", "Nc1ccccc1"},
            {"caffeine", "Cn1cnc2c1c(=O)n(C)c(=O)n2C"},
            {"acetaminophen", "CC(=O)Nc1ccc(O)cc1"},
            {"ibuprofen", "CC(C)Cc1ccc(cc1)C(C)C(=O)O"},
            {"naphthalene", "c1ccc2ccccc2c1"},
            {"sodium", "[Na+]"},
            {"iron", "[Fe]"},
            {"helium", "[He]"},
            {"platinum", "[Pt]"},
            {"uranium", "[U]"},
            {"silicon", "[Si]"},
        };
        // Every molecule is constructed before any pointer is taken: growing
        // the vector afterwards would invalidate the pointers already stored.
        graph_mols_.reserve(std::size(molecules));
        for (const auto& entry : molecules) {
            graph_mols_.emplace_back();
            OEChem::OESmilesToMol(graph_mols_.back(), entry.second);
            graph_mols_.back().SetTitle(entry.first);
        }
        for (auto& gm : graph_mols_) {
            mols_.push_back(&static_cast<OEChem::OEMolBase&>(gm));
        }
    }

    std::vector<std::string> TitlesAt(const std::vector<size_t>& indices) const {
        std::vector<std::string> titles;
        for (const size_t index : indices) {
            titles.emplace_back(mols_[index]->GetTitle());
        }
        return titles;
    }

    std::vector<OEChem::OEMolBase*> Retained() const {
        const std::vector<size_t> excluded =
            descriptor_excluded_indices(mols_, DescriptorOptions());
        std::vector<OEChem::OEMolBase*> kept;
        for (size_t i = 0; i < mols_.size(); ++i) {
            if (std::find(excluded.begin(), excluded.end(), i) == excluded.end()) {
                kept.push_back(mols_[i]);
            }
        }
        return kept;
    }

    std::vector<OEChem::OEGraphMol> graph_mols_;
    std::vector<OEChem::OEMolBase*> mols_;
};

TEST_F(DescriptorMissingnessTest, ExcludedIndicesNamesTheSixUntypedSpecies) {
    const std::vector<size_t> excluded =
        descriptor_excluded_indices(mols_, DescriptorOptions());
    // Compared by title rather than by count. A count-only assertion still
    // passes when the mask slides onto the wrong rows, and tells the reader
    // nothing when it fails.
    EXPECT_EQ(TitlesAt(excluded),
              (std::vector<std::string>{"sodium", "iron", "helium", "platinum",
                                        "uranium", "silicon"}));
}

TEST_F(DescriptorMissingnessTest, TheRetainedSubsetScoresCompletely) {
    const std::vector<OEChem::OEMolBase*> kept = Retained();
    ASSERT_EQ(kept.size(), 13u);

    DescriptorComparison comparison(kept);
    DenseStorage storage(kept.size());
    ASSERT_TRUE(comparison.TryPDist(storage, PDistOptions()));
    for (size_t i = 0; i < kept.size(); ++i) {
        for (size_t j = i + 1; j < kept.size(); ++j) {
            EXPECT_TRUE(std::isfinite(storage.Get(i, j)));
        }
    }
    EXPECT_EQ(comparison.Facts().data_integrity, DataIntegrity::Complete);
}

TEST_F(DescriptorMissingnessTest, UnfilteredCompleteCaseInputIsRejected) {
    // C++ refuses rather than filtering (deviation 4), and the message names
    // the call that produces the mask.
    try {
        DescriptorComparison comparison(mols_, DescriptorOptions());
        FAIL() << "expected ComparisonError";
    } catch (const ComparisonError& exc) {
        EXPECT_NE(std::string(exc.what()).find("descriptor_excluded_indices"),
                  std::string::npos);
    }
}

TEST_F(DescriptorMissingnessTest, PropagateReproducesTheMeasuredNaNCount) {
    // Section 2.4 measured 93 of 171 pairs NaN, which is 171 - C(13, 2): every
    // pair touching one of the 6 untyped species.
    DescriptorOptions opts;
    opts.metric = "euclidean";
    opts.missing = "propagate";
    DescriptorComparison comparison(mols_, opts);
    DenseStorage storage(mols_.size());
    ASSERT_TRUE(comparison.TryPDist(storage, PDistOptions()));

    size_t nan_pairs = 0;
    for (size_t i = 0; i < mols_.size(); ++i) {
        for (size_t j = i + 1; j < mols_.size(); ++j) {
            if (std::isnan(storage.Get(i, j))) {
                ++nan_pairs;
            }
        }
    }
    EXPECT_EQ(nan_pairs, 93u);
    EXPECT_EQ(comparison.Facts().data_integrity, DataIntegrity::NaNPresent);
}

TEST_F(DescriptorMissingnessTest, IgnoreScoresEveryPairAndStampsSubsetScored) {
    // The trap: ignore returns a finite, plausible number for every pair while
    // scoring different pairs over different dimension subsets. Nothing in the
    // values reveals it, which is why the gate reads the integrity stamp.
    //
    // The selection is deliberate. OEFP treats missingness as a property of
    // the validity mask and never of the value, so ignore can only skip a cell
    // whose validity bit is clear. Over this fixture the six untyped species
    // have exactly two gaps each: XLogP is absent, and FractionCsp3 is present
    // and NaN. Including FractionCsp3 would make every pair touching an
    // untyped species NaN under ignore -- the same 93 pairs propagate produces
    // -- and would test nothing about subset scoring. XLogP is the column that
    // makes the policy observable.
    DescriptorOptions opts;
    opts.metric = "euclidean";
    opts.missing = "ignore";
    opts.columns = {"MolecularWeight", "TopologicalPSA", "XLogP"};
    DescriptorComparison comparison(mols_, opts);
    DenseStorage storage(mols_.size());
    ASSERT_TRUE(comparison.TryPDist(storage, PDistOptions()));

    for (size_t i = 0; i < mols_.size(); ++i) {
        for (size_t j = i + 1; j < mols_.size(); ++j) {
            EXPECT_TRUE(std::isfinite(storage.Get(i, j)));
        }
    }
    EXPECT_EQ(comparison.Facts().data_integrity, DataIntegrity::SubsetScored);
}

// A gap that is present-and-NaN rather than absent, on the direct constructor
// path that bypasses the Python normalizer. This is the case the validation
// ordering exists for, and no test above reaches it: every gap in the fixture
// above is absent-typed, so XLogP keeps a finite variance and the row check
// fires no matter when it runs.
//
// Measured against OEFP 0.3.0 and OpenEye 2026.1.0 over water plus the five
// organics: zero absent cells, exactly one present-and-NaN cell (FractionCsp3
// on water), and FractionCsp3 is the only one of the 11 columns whose variance
// is non-finite. Validate after the zero-variance drop and that column is gone,
// the remaining 10 are complete, and the constructor silently scores a molecule
// that descriptor_excluded_indices excludes.
class DescriptorPresentNaNTest : public ::testing::Test {
protected:
    void SetUp() override {
        const char* smiles[] = {"O", "c1ccccc1", "c1ccc(O)cc1", "CCCCCCCC",
                                "CC(=O)Oc1ccccc1C(=O)O", "CCN(CC)CCOC(=O)c1ccccc1"};
        for (const char* smi : smiles) {
            graph_mols_.emplace_back();
            OEChem::OESmilesToMol(graph_mols_.back(), smi);
        }
        for (auto& gm : graph_mols_) {
            mols_.push_back(&static_cast<OEChem::OEMolBase&>(gm));
        }
    }

    std::vector<OEChem::OEGraphMol> graph_mols_;
    std::vector<OEChem::OEMolBase*> mols_;
};

TEST_F(DescriptorPresentNaNTest, WaterIsExcludedForAPresentButNonFiniteValue) {
    EXPECT_EQ(descriptor_excluded_indices(mols_, DescriptorOptions()),
              (std::vector<size_t>{0}));
}

TEST_F(DescriptorPresentNaNTest, TheConstructorRefusesItBeforeFittingDropsTheColumn) {
    try {
        DescriptorComparison comparison(mols_, DescriptorOptions());
        FAIL() << "expected ComparisonError; the NaN column was dropped as "
                  "zero-variance before the row check could see it";
    } catch (const ComparisonError& exc) {
        const std::string message(exc.what());
        EXPECT_NE(message.find("FractionCsp3"), std::string::npos) << message;
        EXPECT_NE(message.find("descriptor_excluded_indices"), std::string::npos)
            << message;
    }
}

TEST_F(DescriptorPresentNaNTest, TheRetainedFiveKeepEveryColumn) {
    // The complement: with water gone nothing is missing, so no column drops and
    // the refusal above is about the row, not about the selection.
    const std::vector<OEChem::OEMolBase*> kept(mols_.begin() + 1, mols_.end());
    DescriptorComparison comparison(kept, DescriptorOptions());
    EXPECT_EQ(comparison.Columns().size(), 11u);
    EXPECT_TRUE(comparison.DroppedColumns().empty());
    EXPECT_EQ(comparison.Facts().data_integrity, DataIntegrity::Complete);
}

TEST(DescriptorComparisonNonFiniteTest, NonFiniteVarianceColumnsAreDropped) {
    // OEFP's RDKit source emits a *present* NaN for BCUT2D and partial-charge
    // descriptors on elements with no Gasteiger parameters, which makes those
    // columns' variances NaN while they still count as present. That is the
    // only route to the !isfinite half of the drop condition; a zero variance
    // alone would not distinguish the two branches.
    //
    // missing='propagate' is what makes the branch reachable at all. Under
    // the default complete_case the row check runs first and refuses this
    // input outright -- see
    // DescriptorPresentNaNTest.TheConstructorRefusesItBeforeFittingDropsTheColumn,
    // which pins that ordering.
    std::vector<OEChem::OEGraphMol> graph_mols(2);
    OEChem::OESmilesToMol(graph_mols[0], "c1ccccc1");
    OEChem::OESmilesToMol(graph_mols[1], "[Na+]");
    std::vector<OEChem::OEMolBase*> mols;
    for (auto& gm : graph_mols) {
        mols.push_back(&static_cast<OEChem::OEMolBase&>(gm));
    }

    DescriptorOptions opts;
    opts.sources = {"rdkit"};
    opts.columns = {"BCUT2D_MWHI", "BCUT2D_LOGPLOW", "MaxPartialCharge"};
    opts.missing = "propagate";
    DescriptorComparison comparison(mols, opts);

    ASSERT_EQ(comparison.Columns().size(), 1u);
    EXPECT_EQ(comparison.Columns()[0], "MaxPartialCharge");
    ASSERT_EQ(comparison.DroppedColumns().size(), 2u);
    EXPECT_NE(std::find(comparison.DroppedColumns().begin(),
                        comparison.DroppedColumns().end(), "BCUT2D_MWHI"),
              comparison.DroppedColumns().end());
    EXPECT_NE(std::find(comparison.DroppedColumns().begin(),
                        comparison.DroppedColumns().end(), "BCUT2D_LOGPLOW"),
              comparison.DroppedColumns().end());
    for (const std::string& reason : comparison.DroppedReasons()) {
        EXPECT_EQ(reason, "zero-variance") << reason;
    }
    EXPECT_TRUE(std::isfinite(comparison.Compare(0, 1)));
}
