#include <algorithm>
#include <cmath>
#include <iterator>
#include <limits>
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

namespace {

/// Return the ComparisonError message a callable throws, or "" if it throws none.
template <typename Callable>
std::string refusal_message(Callable&& call) {
    try {
        call();
    } catch (const ComparisonError& exc) {
        return std::string(exc.what());
    }
    return std::string();
}

/// Option mistakes ``validate_descriptor_options`` refuses, paired with options
/// that produce each. Some need the descriptor schema, which is resolved from
/// the source names alone, so they belong here too. Shared by the two halves of
/// the extraction test: the validator must catch each one, and the constructor
/// must still report each one with the same words.
///
/// Neither a complete list of what the validator refuses nor of what needs no
/// molecules. An unknown source, column or group name needs none either, and
/// what the validator does with one depends on the override: with none it
/// accepts, and the name is reported downstream; with one supplied it refuses
/// the name here. Both halves are pinned by
/// TheSchemaIsResolvedOnlyWhenAnOverrideNeedsIt below, not by this table.
///
/// What this table detects is a rule the validator has *lost*, and
/// EveryTabledMistakeIsRefusedWithoutMolecules is what detects it: the validator
/// falls silent and that test fails on its own, whatever the constructor does
/// next. The byte-for-byte comparison in
/// TheConstructorRepeatsEveryValidatorMessageVerbatim is a second net only where
/// construction still reaches an independent refusal, and for some cases it
/// reaches none -- with the inverse_covariance-under-euclidean rule deleted the
/// constructor stores the override, resolves a euclidean metric that ignores it,
/// and falls silent too, so the two messages still agree. Neither test detects
/// the reverse: a duplicate rule added to the constructor below its validate
/// call is unreachable for every case listed here, because the validator refuses
/// each of them first.
std::vector<std::pair<std::string, DescriptorOptions>> molecule_independent_mistakes() {
    std::vector<std::pair<std::string, DescriptorOptions>> cases;

    DescriptorOptions unknown_metric;
    unknown_metric.metric = "bogus";
    cases.emplace_back("unknown metric", unknown_metric);

    DescriptorOptions wrong_surface;
    wrong_surface.metric = "tanimoto";
    cases.emplace_back("bit-set metric on the descriptor surface", wrong_surface);

    DescriptorOptions bad_exponent;
    bad_exponent.metric = "minkowski";
    bad_exponent.p = -1.0;
    cases.emplace_back("non-positive Minkowski exponent", bad_exponent);

    DescriptorOptions unknown_policy;
    unknown_policy.missing = "drop";
    cases.emplace_back("unknown missing-value policy", unknown_policy);

    DescriptorOptions both_overrides;
    both_overrides.variances = {1.0, 1.0};
    both_overrides.inverse_covariance = {1.0, 0.0, 0.0, 1.0};
    cases.emplace_back("both overrides at once", both_overrides);

    DescriptorOptions variances_on_mahalanobis;
    variances_on_mahalanobis.metric = "mahalanobis";
    variances_on_mahalanobis.variances = {1.0};
    cases.emplace_back("variances under mahalanobis", variances_on_mahalanobis);

    DescriptorOptions inverse_on_euclidean;
    inverse_on_euclidean.metric = "euclidean";
    inverse_on_euclidean.inverse_covariance = {1.0};
    cases.emplace_back("inverse_covariance under euclidean", inverse_on_euclidean);

    DescriptorOptions ignore_on_fitted;
    ignore_on_fitted.missing = "ignore";  // metric defaults to standardized_euclidean
    cases.emplace_back("ignore under a fitted metric", ignore_on_fitted);

    // The seven cases below cover the six rules applied once the selection is
    // resolved: the caller's column order, the two length rules, the two
    // entry-value rules, and the semidefinite verdict. None needs a molecule.
    // All but the inverse_covariance finiteness rule read the schema; that one
    // is applied there only to keep it beside the length rule it follows. Where
    // the fault is in the override's values rather than its length, ``columns``
    // is named explicitly, which keeps the case independent of how many columns
    // the default selection happens to have.
    DescriptorOptions variances_length;
    variances_length.variances = {1.0, 2.0};
    cases.emplace_back("variances length against the selection", variances_length);

    DescriptorOptions inverse_length;
    inverse_length.metric = "mahalanobis";
    inverse_length.inverse_covariance = {1.0, 2.0, 3.0};
    cases.emplace_back("inverse_covariance is not square over the selection", inverse_length);

    DescriptorOptions non_positive_variance;
    non_positive_variance.columns = {"MolecularWeight", "XLogP"};
    non_positive_variance.variances = {0.0, 1.0};
    cases.emplace_back("a zero variance entry", non_positive_variance);

    DescriptorOptions non_finite_variance;
    non_finite_variance.columns = {"MolecularWeight", "XLogP"};
    non_finite_variance.variances = {std::numeric_limits<double>::quiet_NaN(), 1.0};
    cases.emplace_back("a non-finite variance entry", non_finite_variance);

    DescriptorOptions non_finite_inverse;
    non_finite_inverse.metric = "mahalanobis";
    non_finite_inverse.columns = {"MolecularWeight", "XLogP"};
    non_finite_inverse.inverse_covariance = {std::numeric_limits<double>::quiet_NaN(), 0.0, 0.0,
                                             1.0};
    cases.emplace_back("a non-finite inverse_covariance entry", non_finite_inverse);

    DescriptorOptions non_semidefinite_inverse;
    non_semidefinite_inverse.metric = "mahalanobis";
    non_semidefinite_inverse.columns = {"MolecularWeight", "XLogP"};
    non_semidefinite_inverse.inverse_covariance = {-1.0, 0.0, 0.0, -1.0};
    cases.emplace_back("a non-semidefinite inverse_covariance", non_semidefinite_inverse);

    DescriptorOptions descending_columns;
    descending_columns.columns = {"XLogP", "MolecularWeight"};
    descending_columns.variances = {1.5, 2.5};
    cases.emplace_back("columns out of ascending schema order", descending_columns);

    return cases;
}

}  // namespace

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

TEST_F(DescriptorComparisonTest, CompareRefusesAnIndexPastTheEnd) {
    // The row offset ``i * columns`` is read straight out of the descriptor
    // matrix, so an out-of-range index answers about a nonexistent molecule.
    DescriptorComparison comparison(mols_);
    try {
        comparison.Compare(0, 1000000);
        FAIL() << "expected ComparisonError";
    } catch (const ComparisonError& error) {
        const std::string message(error.what());
        EXPECT_NE(message.find("1000000"), std::string::npos) << message;
        EXPECT_NE(message.find("5 items"), std::string::npos) << message;
    }
    EXPECT_THROW(comparison.Compare(1000000, 0), ComparisonError);
    EXPECT_THROW(comparison.Compare(mols_.size(), mols_.size()), ComparisonError);
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

TEST_F(DescriptorComparisonTest, CompareFeedsTheIntegrityStamp) {
    // TryPDist and TryCDist both record a non-finite result; the per-pair path
    // must too, or a caller that only ever calls Compare reads a clean stamp
    // off a matrix that overflowed.
    DescriptorOptions opts;
    opts.metric = "minkowski";
    opts.p = 400.0;
    DescriptorComparison comparison(mols_, opts);
    ASSERT_EQ(comparison.Facts().data_integrity, DataIntegrity::Complete);

    bool saw_non_finite = false;
    for (size_t i = 0; i < comparison.Size() && !saw_non_finite; ++i) {
        for (size_t j = i + 1; j < comparison.Size(); ++j) {
            if (!std::isfinite(comparison.Compare(i, j))) {
                saw_non_finite = true;
                break;
            }
        }
    }
    ASSERT_TRUE(saw_non_finite) << "fixture no longer overflows; pick a larger p";
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

TEST_F(DescriptorComparisonTest, ZeroVarianceInOverrideIsRejected) {
    DescriptorComparison fitted(mols_);
    DescriptorOptions opts;
    opts.columns = fitted.Columns();
    opts.variances = fitted.Variances();
    opts.variances[0] = 0.0;
    try {
        DescriptorComparison comparison(mols_, opts);
        FAIL() << "expected ComparisonError for zero variance";
    } catch (const ComparisonError& exc) {
        const std::string message(exc.what());
        EXPECT_NE(message.find("variances[0]"), std::string::npos) << message;
        EXPECT_NE(message.find(fitted.Columns()[0]), std::string::npos) << message;
    }
}

TEST_F(DescriptorComparisonTest, NaNInInverseCovarianceIsRejected) {
    DescriptorComparison fitted(mols_);
    DescriptorOptions opts;
    opts.metric = "mahalanobis";
    opts.columns = fitted.Columns();
    opts.inverse_covariance = std::vector<double>(fitted.Columns().size() * fitted.Columns().size(), 0.0);
    opts.inverse_covariance[0] = std::numeric_limits<double>::quiet_NaN();
    try {
        DescriptorComparison comparison(mols_, opts);
        FAIL() << "expected ComparisonError for NaN in inverse_covariance";
    } catch (const ComparisonError& exc) {
        const std::string message(exc.what());
        EXPECT_NE(message.find("inverse_covariance[0]"), std::string::npos) << message;
    }
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

TEST_F(DescriptorComparisonTest, TheSuppliedInverseCovarianceIsTheMatrixTheScoreUses) {
    // Every other test that supplies a valid inverse_covariance either expects a
    // refusal or asserts only that a score came back finite, so nothing read the
    // matrix itself: dropping the assignment that carries it from the options
    // into the metric left the whole suite green.
    //
    // A diagonal matrix is what makes the value checkable by hand. Mahalanobis
    // over one reduces to a weighted Euclidean sum, so each column's gap can be
    // measured on its own -- a single-column euclidean comparison is exactly
    // that gap -- and recombined here. The two weights are distinct and neither
    // is 1, so a matrix read transposed, replaced by the identity, or refitted
    // from the molecules all move the answer.
    DescriptorOptions weight_only;
    weight_only.metric = "euclidean";
    weight_only.columns = {"MolecularWeight"};
    DescriptorComparison weight_gap(mols_, weight_only);

    DescriptorOptions logp_only;
    logp_only.metric = "euclidean";
    logp_only.columns = {"XLogP"};
    DescriptorComparison logp_gap(mols_, logp_only);

    DescriptorOptions opts;
    opts.metric = "mahalanobis";
    opts.columns = {"MolecularWeight", "XLogP"};
    opts.inverse_covariance = {4.0, 0.0, 0.0, 9.0};
    DescriptorComparison scored(mols_, opts);

    for (size_t i = 0; i < mols_.size(); ++i) {
        for (size_t j = i + 1; j < mols_.size(); ++j) {
            const double weight_delta = weight_gap.Compare(i, j);
            const double logp_delta = logp_gap.Compare(i, j);
            const double expected = std::sqrt(4.0 * weight_delta * weight_delta +
                                              9.0 * logp_delta * logp_delta);
            ASSERT_GT(expected, 0.0) << "pair (" << i << ", " << j << ") is degenerate";
            EXPECT_NEAR(scored.Compare(i, j), expected, 1e-9 * expected)
                << "pair (" << i << ", " << j << ")";
        }
    }

    // Reported last, and deliberately not as the gate on the loop above: the
    // accessor and the metric read the same member, so an accessor check on its
    // own would pass for a matrix that never reached the metric.
    EXPECT_EQ(scored.InverseCovariance(), opts.inverse_covariance);
}

TEST_F(DescriptorComparisonTest, DropReportIsParallel) {
    // Not written against mols_: all eleven columns have positive variance over
    // that fixture, so the two sizes would both be zero and the assertion would
    // hold against an implementation that never populates either vector. These
    // four are each a single aromatic ring with no rotatable bonds, which
    // collapses three columns to zero variance under the fitted default metric.
    std::vector<OEChem::OEGraphMol> ring_mols(4);
    OEChem::OESmilesToMol(ring_mols[0], "c1ccccc1");
    OEChem::OESmilesToMol(ring_mols[1], "c1ccc(O)cc1");
    OEChem::OESmilesToMol(ring_mols[2], "Cc1ccccc1");
    OEChem::OESmilesToMol(ring_mols[3], "Nc1ccccc1");
    std::vector<OEChem::OEMolBase*> rings;
    for (auto& gm : ring_mols) {
        rings.push_back(&static_cast<OEChem::OEMolBase&>(gm));
    }

    DescriptorComparison comparison(rings);
    ASSERT_FALSE(comparison.DroppedColumns().empty());
    EXPECT_EQ(comparison.DroppedColumns().size(), comparison.DroppedReasons().size());
}

TEST_F(DescriptorComparisonTest, DescendingColumnsWithAVarianceOverrideAreRefused) {
    // The two orderings have the same length, so only an order check catches
    // this. Scoring it would standardize XLogP by MolecularWeight's variance.
    DescriptorOptions opts;
    opts.metric = "standardized_euclidean";
    opts.columns = {"XLogP", "MolecularWeight"};
    opts.variances = {1.5, 2.5};
    try {
        DescriptorComparison comparison(mols_, opts);
        FAIL() << "expected ComparisonError for a descending column order";
    } catch (const ComparisonError& exc) {
        const std::string message(exc.what());
        EXPECT_NE(message.find("ascending schema order"), std::string::npos) << message;
        EXPECT_NE(message.find("MolecularWeight, XLogP"), std::string::npos) << message;
    }
}

TEST_F(DescriptorComparisonTest, AscendingColumnsWithAVarianceOverrideAreAccepted) {
    DescriptorOptions opts;
    opts.metric = "standardized_euclidean";
    opts.columns = {"MolecularWeight", "XLogP"};
    opts.variances = {2.5, 1.5};
    DescriptorComparison comparison(mols_, opts);
    EXPECT_EQ(comparison.Columns().size(), 2u);
    EXPECT_TRUE(std::isfinite(comparison.Compare(0, 1)));
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

TEST(DescriptorOptionValidationTest, DefaultOptionsAreAccepted) {
    EXPECT_NO_THROW(validate_descriptor_options(DescriptorOptions()));
}

TEST(DescriptorOptionValidationTest, EveryTabledMistakeIsRefusedWithoutMolecules) {
    // No molecules anywhere in this test: that is the point of the entry
    // point. The Python layer calls it before the complete-case filter, which
    // can empty the item list and answer for a filter the caller never asked
    // for instead of the option they have to change.
    for (const auto& entry : molecule_independent_mistakes()) {
        const std::string message =
            refusal_message([&] { validate_descriptor_options(entry.second); });
        EXPECT_FALSE(message.empty()) << "not refused: " << entry.first;
    }
}

TEST(DescriptorOptionValidationTest, TheUnknownMetricMessageNamesTheSurfacesMetrics) {
    DescriptorOptions opts;
    opts.metric = "bogus";
    const std::string message =
        refusal_message([&] { validate_descriptor_options(opts); });
    EXPECT_NE(message.find("Unknown metric 'bogus'"), std::string::npos) << message;
    EXPECT_NE(message.find("standardized_euclidean"), std::string::npos) << message;
    // The fingerprint-only names must not be advertised here.
    EXPECT_EQ(message.find("tanimoto"), std::string::npos) << message;
}

TEST(DescriptorOptionValidationTest, TheMetricNameIsFoldedBeforeItIsReported) {
    // The constructor lowercases before resolving, so the validator must too
    // or the same input would be refused with a differently-spelled message
    // depending on which check reached it first.
    DescriptorOptions opts;
    opts.metric = "BOGUS";
    const std::string message =
        refusal_message([&] { validate_descriptor_options(opts); });
    EXPECT_NE(message.find("Unknown metric 'bogus'"), std::string::npos) << message;
}

TEST(DescriptorOptionValidationTest, TheSchemaIsResolvedOnlyWhenAnOverrideNeedsIt) {
    // The boundary is "needs no molecules", so the schema is fair game -- but
    // no rule reads it unless an override is supplied, and building a
    // calculator to check nothing would tax every plain request. An unknown
    // name is therefore still reported downstream on the plain path, which is
    // where descriptor_excluded_indices and the constructor both resolve the
    // same selection and report it themselves.
    DescriptorOptions unknown_source;
    unknown_source.sources = {"nosuchsource"};
    EXPECT_NO_THROW(validate_descriptor_options(unknown_source));

    DescriptorOptions unknown_column;
    unknown_column.columns = {"NoSuchColumn"};
    EXPECT_NO_THROW(validate_descriptor_options(unknown_column));

    DescriptorOptions unknown_group;
    unknown_group.groups = {"nosuchgroup"};
    EXPECT_NO_THROW(validate_descriptor_options(unknown_group));

    // With an override the schema has to be resolved, so the same names are
    // reported here instead. Same message either way; only the messenger moves.
    unknown_source.variances = {1.0};
    EXPECT_NE(refusal_message([&] { validate_descriptor_options(unknown_source); })
                  .find("Unknown descriptor source: nosuchsource"),
              std::string::npos);

    unknown_column.variances = {1.0};
    EXPECT_NE(refusal_message([&] { validate_descriptor_options(unknown_column); })
                  .find("Unknown descriptor column: NoSuchColumn"),
              std::string::npos);

    unknown_group.variances = {1.0};
    EXPECT_NE(refusal_message([&] { validate_descriptor_options(unknown_group); })
                  .find("Unknown or empty descriptor group: nosuchgroup"),
              std::string::npos);
}

TEST(DescriptorOptionValidationTest, ANonSemidefiniteInverseCovarianceIsRefusedWithoutMolecules) {
    // The verdict is read off an eigendecomposition, which only OEFP performs,
    // so the validator obtains it by scoring a probe rather than by keeping a
    // copy of the rule. Left to the first real scoring call it would surface
    // long after the call that supplied the matrix, and through a wrapper that
    // blames the pdist.
    DescriptorOptions opts;
    opts.metric = "mahalanobis";
    opts.columns = {"MolecularWeight", "XLogP"};
    opts.inverse_covariance = {-1.0, 0.0, 0.0, -1.0};
    EXPECT_NE(refusal_message([&] { validate_descriptor_options(opts); })
                  .find("positive semidefinite"),
              std::string::npos);

    // Indefinite with a positive diagonal, so a diagonal-only screen would wave
    // it through; separating it takes a verdict on the matrix as a whole.
    opts.inverse_covariance = {1.0, 2.0, 2.0, 1.0};
    EXPECT_NE(refusal_message([&] { validate_descriptor_options(opts); })
                  .find("positive semidefinite"),
              std::string::npos);

    // An asymmetric matrix is read as (M + M^T) / 2 rather than refused, so it
    // is accepted whenever that symmetric part qualifies. Refusing here would
    // reject a matrix OEFP goes on to score.
    opts.inverse_covariance = {1.0, 0.5, -0.5, 1.0};
    EXPECT_NO_THROW(validate_descriptor_options(opts));

    // A cheaper rule still outranks it. The decomposition is the most expensive
    // check here, and a caller whose matrix is also the wrong shape has to fix
    // the shape first.
    opts.inverse_covariance = {-1.0, 0.0, 0.0};
    EXPECT_NE(refusal_message([&] { validate_descriptor_options(opts); }).find("square 2x2"),
              std::string::npos);
}

TEST(DescriptorOptionValidationTest, TheInputSizeRuleIsLeftToTheConstructor) {
    // A fitted metric over one molecule is a valid *option* set; only the
    // input makes it impossible. Refusing it here would refuse it for every
    // caller, including the ones passing enough molecules.
    EXPECT_NO_THROW(validate_descriptor_options(DescriptorOptions()));

    std::vector<OEChem::OEGraphMol> graph_mols(1);
    OEChem::OESmilesToMol(graph_mols[0], "c1ccccc1");
    std::vector<OEChem::OEMolBase*> mols{&static_cast<OEChem::OEMolBase&>(graph_mols[0])};
    const std::string message =
        refusal_message([&] { DescriptorComparison(mols, DescriptorOptions()); });
    EXPECT_NE(message.find("at least two molecules"), std::string::npos) << message;
}

TEST(DescriptorOptionValidationTest, TheAllColumnsConstantRuleIsLeftToTheConstructor) {
    // One of the rules the constructor keeps because it has to read the input
    // to decide, and the one DescriptorComparison.h's prose once omitted: a
    // fitted metric needs spread to fit against, and whether there is any is a
    // fact about the molecules. Deliberately no ordinal -- NullMoleculeIsRejected
    // above covers another input-reading refusal, and the set is open.
    // Identical inputs make every selected column constant, which no option
    // value can be blamed for and no validator could have foreseen.
    DescriptorOptions single_column;
    single_column.columns = {"MolecularWeight"};
    EXPECT_NO_THROW(validate_descriptor_options(single_column));

    std::vector<OEChem::OEGraphMol> graph_mols(2);
    OEChem::OESmilesToMol(graph_mols[0], "CCO");
    OEChem::OESmilesToMol(graph_mols[1], "CCO");
    std::vector<OEChem::OEMolBase*> mols{&static_cast<OEChem::OEMolBase&>(graph_mols[0]),
                                         &static_cast<OEChem::OEMolBase&>(graph_mols[1])};
    const std::string message =
        refusal_message([&] { DescriptorComparison(mols, single_column); });
    EXPECT_NE(message.find("Every selected descriptor column has zero variance"),
              std::string::npos)
        << message;

    // It is a property of the fit, not of the molecules alone: an unfitted
    // metric asks nothing of the spread and must still score the same input.
    DescriptorOptions unfitted = single_column;
    unfitted.metric = "euclidean";
    EXPECT_NO_THROW(DescriptorComparison(mols, unfitted));
}

TEST_F(DescriptorComparisonTest, TheConstructorRepeatsEveryValidatorMessageVerbatim) {
    // The extraction must not have changed what a caller is told, only when:
    // for every tabled refusal the constructor path reports the validator's own
    // wording. It does not also keep the constructor from growing a second copy
    // of a rule -- a duplicate below the validate call is unreachable for every
    // case in the table, and the suite passed with one inserted. The caption on
    // molecule_independent_mistakes() has the rest.
    for (const auto& entry : molecule_independent_mistakes()) {
        const std::string from_validator =
            refusal_message([&] { validate_descriptor_options(entry.second); });
        const std::string from_constructor =
            refusal_message([&] { DescriptorComparison(mols_, entry.second); });
        EXPECT_EQ(from_constructor, from_validator) << entry.first;
    }
}

TEST(DescriptorComparisonMinimumInputTest, AFittedMetricNeedsTwoMolecules) {
    std::vector<OEChem::OEGraphMol> graph_mols(1);
    OEChem::OESmilesToMol(graph_mols[0], "c1ccccc1");
    std::vector<OEChem::OEMolBase*> mols{&static_cast<OEChem::OEMolBase&>(graph_mols[0])};

    try {
        DescriptorComparison comparison(mols, DescriptorOptions());
        FAIL() << "expected ComparisonError naming the input size";
    } catch (const ComparisonError& exc) {
        const std::string message(exc.what());
        EXPECT_NE(message.find("at least two molecules"), std::string::npos) << message;
        EXPECT_EQ(message.find("zero variance"), std::string::npos)
            << "the error must not blame the descriptors: " << message;
    }
}

TEST(DescriptorComparisonMinimumInputTest, AnUnfittedMetricAcceptsOneMolecule) {
    std::vector<OEChem::OEGraphMol> graph_mols(1);
    OEChem::OESmilesToMol(graph_mols[0], "c1ccccc1");
    std::vector<OEChem::OEMolBase*> mols{&static_cast<OEChem::OEMolBase&>(graph_mols[0])};

    DescriptorOptions opts;
    opts.metric = "euclidean";
    DescriptorComparison comparison(mols, opts);
    EXPECT_EQ(comparison.Size(), 1u);
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

TEST_F(DescriptorMissingnessTest, IgnoreEscalatesToNaNPresentWhenNaNActuallyReaches) {
    // ignore is a tier-2 stamp only while it stays NaN-free. Over the full
    // default selection it does not: FractionCsp3 is present-and-NaN for the
    // untyped species, and OEFP never drops a present value, so NaN reaches the
    // matrix. Stamping SubsetScored here would let allow_nonmetric=True admit a
    // matrix with NaN in it.
    DescriptorOptions opts;
    opts.metric = "euclidean";
    opts.missing = "ignore";
    DescriptorComparison comparison(mols_, opts);

    // Before scoring, the declared stamp stands.
    EXPECT_EQ(comparison.Facts().data_integrity, DataIntegrity::SubsetScored);

    DenseStorage storage(mols_.size());
    ASSERT_TRUE(comparison.TryPDist(storage, PDistOptions()));
    size_t nan_count = 0;
    for (size_t i = 0; i < mols_.size(); ++i) {
        for (size_t j = i + 1; j < mols_.size(); ++j) {
            if (!std::isfinite(storage.Get(i, j))) {
                ++nan_count;
            }
        }
    }
    ASSERT_GT(nan_count, 0u) << "fixture no longer produces NaN under ignore; "
                                "the escalation is untested";

    EXPECT_EQ(comparison.Facts().data_integrity, DataIntegrity::NaNPresent);
}

// A gap that is present-and-NaN rather than absent, on the direct constructor
// path that bypasses the Python normalizer. This is the case the validation
// ordering exists for, and no test above reaches it -- though not for want of a
// present-and-NaN gap. The missingness fixture has one, FractionCsp3 on its six
// untyped species. It also has an absent gap on those same six rows, XLogP, so
// the row check fires whether it runs before or after the zero-variance drop
// removes FractionCsp3. Isolating the ordering needs a fixture whose only gap
// is the present-and-NaN one, which is what this one is.
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
