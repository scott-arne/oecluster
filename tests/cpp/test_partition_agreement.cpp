#include <gtest/gtest.h>

#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "oecluster/clustering/ClusterTypes.h"
#include "oecluster/clustering/PartitionAgreement.h"
#include "../../src/clustering/ContingencyTable.h"

using namespace OECluster;

namespace {

using Pair = std::tuple<uint64_t, uint64_t, uint64_t>;

// Record every callback for_each_marginal_pair makes, in order.
std::vector<Pair> RecordPairs(const std::vector<detail::SizeMultiplicity>& a,
                              const std::vector<detail::SizeMultiplicity>& b) {
    std::vector<Pair> recorded;
    detail::for_each_marginal_pair(
        a, b, [&recorded](uint64_t u, uint64_t w, uint64_t multiplicity) {
            recorded.emplace_back(u, w, multiplicity);
        });
    return recorded;
}

std::vector<std::pair<uint64_t, uint32_t>> Flatten(
    const std::vector<detail::SizeMultiplicity>& histogram) {
    std::vector<std::pair<uint64_t, uint32_t>> flat;
    for (const detail::SizeMultiplicity& entry : histogram) {
        flat.emplace_back(entry.size, entry.multiplicity);
    }
    return flat;
}

}  // namespace

TEST(ContingencyTableTest, InternsInFirstAppearanceOrder) {
    const std::vector<ClusterLabel> labels{7, 7, 3, 7, 3};
    const std::vector<bool> drop(labels.size(), false);
    uint32_t num_ids = 0;
    const std::vector<uint32_t> ids =
        detail::intern_side(labels, NoiseHandling::Singletons, drop, num_ids);

    EXPECT_EQ(num_ids, 2u);
    EXPECT_EQ(ids, (std::vector<uint32_t>{0, 0, 1, 0, 1}));
}

TEST(ContingencyTableTest, SingletonsGiveEachNoiseSampleItsOwnId) {
    const std::vector<ClusterLabel> labels{0, -1, 0, -2};
    const std::vector<bool> drop(labels.size(), false);
    uint32_t num_ids = 0;
    const std::vector<uint32_t> ids =
        detail::intern_side(labels, NoiseHandling::Singletons, drop, num_ids);

    EXPECT_EQ(num_ids, 3u);
    EXPECT_EQ(ids, (std::vector<uint32_t>{0, 1, 0, 2}));
}

TEST(ContingencyTableTest, GroupedCollapsesNoiseIntoOneId) {
    const std::vector<ClusterLabel> labels{0, -1, 0, -2};
    const std::vector<bool> drop(labels.size(), false);
    uint32_t num_ids = 0;
    const std::vector<uint32_t> ids =
        detail::intern_side(labels, NoiseHandling::Grouped, drop, num_ids);

    EXPECT_EQ(num_ids, 2u);
    EXPECT_EQ(ids, (std::vector<uint32_t>{0, 1, 0, 1}));
}

TEST(ContingencyTableTest, ExcludedDropMaskIsTheUnionOfBothSides) {
    // Side A marks sample 1, side B marks sample 3. Both must go.
    const std::vector<ClusterLabel> a{0, -1, 1, 1};
    const std::vector<std::string> b{"x", "y", "y", ""};

    std::vector<bool> drop(a.size(), false);
    detail::mark_excluded(a, NoiseHandling::Excluded, drop);
    detail::mark_excluded(b, NoiseHandling::Excluded, drop);

    EXPECT_EQ(drop, (std::vector<bool>{false, true, false, true}));
}

TEST(ContingencyTableTest, CellsAreSortedByRowThenColumn) {
    // Interning in first-appearance order yields unsorted (row, col) keys, so
    // the sort in step 4 has visible work to do here.
    const std::vector<ClusterLabel> a{2, 1, 0, 2, 1, 0};
    const std::vector<ClusterLabel> b{2, 1, 0, 0, 2, 1};
    const detail::ContingencyTable table =
        detail::build_contingency(a, b, NoiseHandling::Singletons);

    ASSERT_FALSE(table.cells.empty());
    for (size_t i = 1; i < table.cells.size(); ++i) {
        const detail::ContingencyTable::Cell& previous = table.cells[i - 1];
        const detail::ContingencyTable::Cell& current = table.cells[i];
        EXPECT_TRUE(previous.row < current.row ||
                    (previous.row == current.row &&
                     previous.col < current.col));
    }
}

TEST(ContingencyTableTest, MarginalsAndSampleCountFollowExclusion) {
    const std::vector<ClusterLabel> a{0, 0, 0, -1, 1, 1, 1, 1, 2, 2};
    const std::vector<ClusterLabel> b{0, 0, 1, 1, 1, 1, 2, -1, 2, 2};
    const detail::ContingencyTable table =
        detail::build_contingency(a, b, NoiseHandling::Excluded);

    EXPECT_EQ(table.num_samples, 8u);
    EXPECT_EQ(table.marginals_a.size(), 3u);
    EXPECT_EQ(table.marginals_b.size(), 3u);
    EXPECT_EQ(table.marginals_a, (std::vector<uint64_t>{3, 3, 2}));
    EXPECT_EQ(table.marginals_b, (std::vector<uint64_t>{2, 3, 3}));
}

TEST(ContingencyTableTest, EmptyScaffoldStringIsNoise) {
    const std::vector<ClusterLabel> labels{0, 0, 1, 1};
    const std::vector<std::string> scaffolds{"ar", "", "pi", ""};
    const detail::ContingencyTable grouped =
        detail::build_contingency(labels, scaffolds, NoiseHandling::Grouped);
    const detail::ContingencyTable singletons =
        detail::build_contingency(labels, scaffolds, NoiseHandling::Singletons);

    EXPECT_EQ(grouped.marginals_b.size(), 3u);
    EXPECT_EQ(singletons.marginals_b.size(), 4u);
}

TEST(ContingencyTableTest, MarginalHistogramIsAscendingAndDeduplicated) {
    using Entry = std::pair<uint64_t, uint32_t>;

    // The AMIFIX marginals.
    EXPECT_EQ(Flatten(detail::marginal_histogram({1, 2, 3, 3, 4, 5})),
              (std::vector<Entry>{{1, 1}, {2, 1}, {3, 2}, {4, 1}, {5, 1}}));
    EXPECT_EQ(Flatten(detail::marginal_histogram({1, 3, 4, 4, 6})),
              (std::vector<Entry>{{1, 1}, {3, 1}, {4, 2}, {6, 1}}));
    // Every size distinct, given out of order: every multiplicity is 1.
    EXPECT_EQ(Flatten(detail::marginal_histogram({4, 1, 3, 2})),
              (std::vector<Entry>{{1, 1}, {2, 1}, {3, 1}, {4, 1}}));
    // Every size equal: one entry.
    EXPECT_EQ(Flatten(detail::marginal_histogram({4, 4, 4})),
              (std::vector<Entry>{{4, 3}}));
}

TEST(ContingencyTableTest, MarginalPairTraversalIsCompleteAndOrdered) {
    const std::vector<detail::SizeMultiplicity> hist_a{
        {1, 1}, {2, 1}, {3, 2}, {4, 1}, {5, 1}};
    const std::vector<detail::SizeMultiplicity> hist_b{
        {1, 1}, {3, 1}, {4, 2}, {6, 1}};

    // The whole sequence: every Cartesian pair exactly once, each carrying
    // cnt_a[u] * cnt_b[w], outer-ascending u and inner-ascending w. An
    // aggregate -- a pair count, or a multiplicity total -- would survive a
    // per-pair error that cancels, and would not see the order at all.
    EXPECT_EQ(RecordPairs(hist_a, hist_b),
              (std::vector<Pair>{
                  {1, 1, 1}, {1, 3, 1}, {1, 4, 2}, {1, 6, 1},
                  {2, 1, 1}, {2, 3, 1}, {2, 4, 2}, {2, 6, 1},
                  {3, 1, 2}, {3, 3, 2}, {3, 4, 4}, {3, 6, 2},
                  {4, 1, 1}, {4, 3, 1}, {4, 4, 2}, {4, 6, 1},
                  {5, 1, 1}, {5, 3, 1}, {5, 4, 2}, {5, 6, 1}}));
}

TEST(ContingencyTableTest, MarginalPairMultiplicityDoesNotWrapAt32Bits) {
    // Two 70004-sample partitions of 70000 singleton clusters each reach this.
    // A 32-bit multiply gives 605032704.
    const std::vector<detail::SizeMultiplicity> hist{{1, 70000}};
    const std::vector<Pair> recorded = RecordPairs(hist, hist);

    ASSERT_EQ(recorded.size(), 1u);
    EXPECT_EQ(std::get<2>(recorded[0]), 4900000000ull);
}

namespace {

// A ClusteringResult with members derived from the labels, as the other
// clustering tests build one.
ClusteringResult MakeResult(std::vector<ClusterLabel> labels) {
    Clusters members = labels_to_clusters(labels);
    return ClusteringResult(std::move(labels), std::move(members));
}

PartitionAgreementOptions WithNoise(NoiseHandling noise_handling) {
    PartitionAgreementOptions options;
    options.noise_handling = noise_handling;
    return options;
}

PartitionAgreementOptions WithAmi() {
    PartitionAgreementOptions options;
    options.compute_adjusted_mutual_information = true;
    return options;
}

// 18 samples, cluster sizes {1,2,3,3,4,5} against {1,3,4,4,6}: 9 nonzero cells
// of 30, so 21 zero cells contribute to E[MI].
const std::vector<ClusterLabel> kAmiA{0, 1, 1, 2, 2, 2, 3, 3, 3,
                                      4, 4, 4, 4, 5, 5, 5, 5, 5};
const std::vector<ClusterLabel> kAmiB{0, 1, 1, 1, 2, 2, 2, 2, 3,
                                      3, 3, 3, 4, 4, 4, 4, 4, 4};

/// A deliberately naive K_a x K_b E[MI]: every marginal pair visited on its
/// own, with no grouping by marginal value. The per-pair term comes from
/// production so that what this reference contradicts is the grouping, not a
/// second transcription of Vinh et al. (2010) that could drift from the first.
double NaiveExpectedMutualInformation(const detail::ContingencyTable& table) {
    const std::vector<double> logfact =
        detail::log_factorials(table.num_samples);
    double expected = 0.0;
    for (uint64_t u : table.marginals_a) {
        for (uint64_t w : table.marginals_b) {
            expected += detail::expected_term(u, w, table.num_samples, logfact);
        }
    }
    return expected;
}

// NaN != NaN, so pairing isnan is the only way to compare a field that is
// undefined on both sides. That is what lets the tests below assert on every
// field from Task 2 onwards: a field this task leaves at its NaN default
// compares equal to itself now, and the same assertion widens by itself to the
// real value once Task 3 or Task 4 fills the field in. Defined values are
// compared bitwise rather than with EXPECT_DOUBLE_EQ, because the two calls
// under comparison run identical arithmetic and anything but an exact match is
// a bug.
void ExpectSameDouble(double actual, double expected) {
    if (std::isnan(expected)) {
        EXPECT_TRUE(std::isnan(actual));
    } else {
        EXPECT_EQ(actual, expected);
    }
}

void ExpectSameAgreement(const PartitionAgreement& actual,
                         const PartitionAgreement& expected) {
    EXPECT_EQ(actual.num_samples, expected.num_samples);
    EXPECT_EQ(actual.num_clusters_a, expected.num_clusters_a);
    EXPECT_EQ(actual.num_clusters_b, expected.num_clusters_b);
    EXPECT_EQ(actual.requested.adjusted_mutual_information,
              expected.requested.adjusted_mutual_information);
    ExpectSameDouble(actual.adjusted_rand_index, expected.adjusted_rand_index);
    ExpectSameDouble(actual.fowlkes_mallows, expected.fowlkes_mallows);
    ExpectSameDouble(actual.normalized_mutual_information,
                     expected.normalized_mutual_information);
    ExpectSameDouble(actual.homogeneity, expected.homogeneity);
    ExpectSameDouble(actual.completeness, expected.completeness);
    ExpectSameDouble(actual.v_measure, expected.v_measure);
    ExpectSameDouble(actual.adjusted_mutual_information,
                     expected.adjusted_mutual_information);
}

const std::vector<ClusterLabel> kMainA{0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2};
const std::vector<ClusterLabel> kMainB{0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3};
const std::vector<ClusterLabel> kNoiseA{0, 0, 0, -1, 1, 1, 1, -1, 2, 2, -1, 2};
const std::vector<ClusterLabel> kNoiseB{0, 0, 1, 1, 1, -1, 2, 2, 2, 3, 3, -1};

}  // namespace

// sklearn.metrics.adjusted_rand_score / fowlkes_mallows_score on kMainA/kMainB.
TEST(PartitionAgreementTest, MainFixturePairMetrics) {
    const PartitionAgreement agreement = partition_agreement(kMainA, kMainB);

    EXPECT_EQ(agreement.num_samples, 12u);
    EXPECT_EQ(agreement.num_clusters_a, 3u);
    EXPECT_EQ(agreement.num_clusters_b, 4u);
    EXPECT_NEAR(agreement.adjusted_rand_index, 0.40310077519379844, 1e-12);
    EXPECT_NEAR(agreement.fowlkes_mallows, 0.5443310539518174, 1e-12);
}

TEST(PartitionAgreementTest, SelfAgreementIsOneEverywhere) {
    const std::vector<std::vector<ClusterLabel>> fixtures{
        {0, 0, 1, 1, 2, 2},  // a normal partition
        {0, 1, 2, 3, 4, 5},  // all singletons
        {0, 0, 0, 0, 0, 0},  // a single cluster
    };
    for (const std::vector<ClusterLabel>& labels : fixtures) {
        const PartitionAgreement agreement = partition_agreement(labels, labels);
        EXPECT_DOUBLE_EQ(agreement.adjusted_rand_index, 1.0);
        EXPECT_DOUBLE_EQ(agreement.fowlkes_mallows, 1.0);
        EXPECT_DOUBLE_EQ(agreement.normalized_mutual_information, 1.0);
        EXPECT_DOUBLE_EQ(agreement.homogeneity, 1.0);
        EXPECT_DOUBLE_EQ(agreement.completeness, 1.0);
        EXPECT_DOUBLE_EQ(agreement.v_measure, 1.0);
    }
}

// Rule 2 is about the grouping, not the label values.
TEST(PartitionAgreementTest, RelabelledPartitionsAreIdentical) {
    const PartitionAgreement agreement =
        partition_agreement({0, 0, 1, 1}, {9, 9, 4, 4});
    EXPECT_DOUBLE_EQ(agreement.adjusted_rand_index, 1.0);
    EXPECT_DOUBLE_EQ(agreement.fowlkes_mallows, 1.0);
}

// Rule 2 outranks rule 3: scikit-learn reports FM = 0.0 here.
TEST(PartitionAgreementTest, BothSidesAllSingletonsAreIdentical) {
    const PartitionAgreement agreement =
        partition_agreement({0, 1, 2, 3}, {3, 2, 1, 0});
    EXPECT_DOUBLE_EQ(agreement.fowlkes_mallows, 1.0);
    EXPECT_DOUBLE_EQ(agreement.adjusted_rand_index, 1.0);
}

TEST(PartitionAgreementTest, PairMetricsAreSymmetric) {
    const PartitionAgreement forward = partition_agreement(kMainA, kMainB);
    const PartitionAgreement backward = partition_agreement(kMainB, kMainA);
    EXPECT_DOUBLE_EQ(forward.adjusted_rand_index, backward.adjusted_rand_index);
    EXPECT_DOUBLE_EQ(forward.fowlkes_mallows, backward.fowlkes_mallows);
    EXPECT_EQ(forward.num_clusters_a, backward.num_clusters_b);
    EXPECT_EQ(forward.num_clusters_b, backward.num_clusters_a);
}

TEST(PartitionAgreementTest, ThreeNoiseModesGiveDistinctResults) {
    const PartitionAgreement singletons =
        partition_agreement(kNoiseA, kNoiseB, WithNoise(NoiseHandling::Singletons));
    EXPECT_EQ(singletons.num_samples, 12u);
    EXPECT_EQ(singletons.num_clusters_a, 6u);
    EXPECT_EQ(singletons.num_clusters_b, 6u);
    EXPECT_NEAR(singletons.adjusted_rand_index, -0.012269938650306749, 1e-12);
    EXPECT_NEAR(singletons.fowlkes_mallows, 0.11785113019775792, 1e-12);

    const PartitionAgreement grouped =
        partition_agreement(kNoiseA, kNoiseB, WithNoise(NoiseHandling::Grouped));
    EXPECT_EQ(grouped.num_samples, 12u);
    EXPECT_EQ(grouped.num_clusters_a, 4u);
    EXPECT_EQ(grouped.num_clusters_b, 5u);
    EXPECT_NEAR(grouped.adjusted_rand_index, -0.07179487179487179, 1e-12);
    EXPECT_NEAR(grouped.fowlkes_mallows, 0.09622504486493762, 1e-12);

    const PartitionAgreement excluded =
        partition_agreement(kNoiseA, kNoiseB, WithNoise(NoiseHandling::Excluded));
    EXPECT_EQ(excluded.num_samples, 7u);
    EXPECT_EQ(excluded.num_clusters_a, 3u);
    EXPECT_EQ(excluded.num_clusters_b, 4u);
    EXPECT_NEAR(excluded.adjusted_rand_index, 0.08695652173913043, 1e-12);
    EXPECT_NEAR(excluded.fowlkes_mallows, 0.2581988897471611, 1e-12);
}

TEST(PartitionAgreementTest, ExcludedDropsSamplesNoisyOnEitherSide) {
    const std::vector<ClusterLabel> a{0, 0, 0, -1, 1, 1, 1, 1, 2, 2};
    const std::vector<ClusterLabel> b{0, 0, 1, 1, 1, 1, 2, -1, 2, 2};

    EXPECT_EQ(partition_agreement(a, b, WithNoise(NoiseHandling::Singletons))
                  .num_samples,
              10u);
    const PartitionAgreement excluded =
        partition_agreement(a, b, WithNoise(NoiseHandling::Excluded));
    EXPECT_EQ(excluded.num_samples, 8u);
    EXPECT_NEAR(excluded.adjusted_rand_index, 0.23809523809523808, 1e-12);
    EXPECT_NEAR(excluded.fowlkes_mallows, 0.4285714285714285, 1e-12);
}

TEST(PartitionAgreementTest, LabelsBelowMinusOneAreNoise) {
    const std::vector<ClusterLabel> minus_one{0, -1, 1, -1, 2, 2};
    const std::vector<ClusterLabel> assorted{0, -2, 1, -3, 2, 2};
    const std::vector<ClusterLabel> other{0, 0, 1, 1, 2, 2};

    // Every field, not just the pair metrics: the requirement is that the two
    // fixtures score identically, and a divergence confined to the entropy
    // metrics would be just as much a defect. AMI is requested so this covers
    // the last field too — at this task both sides are still NaN and the
    // NaN-aware compare passes; from Task 4 the same assertion compares real
    // values.
    for (NoiseHandling mode : {NoiseHandling::Singletons, NoiseHandling::Grouped,
                               NoiseHandling::Excluded}) {
        PartitionAgreementOptions options = WithNoise(mode);
        options.compute_adjusted_mutual_information = true;
        ExpectSameAgreement(partition_agreement(assorted, other, options),
                            partition_agreement(minus_one, other, options));
    }
}

TEST(PartitionAgreementTest, FowlkesMallowsIsNaNWhenEitherSideIsAllSingletons) {
    // sum_a == 0 and sum_b == 0 are separate guards; cover both.
    const PartitionAgreement side_a =
        partition_agreement({0, 1, 2, 3}, {0, 0, 1, 1});
    EXPECT_TRUE(std::isnan(side_a.fowlkes_mallows));
    EXPECT_FALSE(std::isnan(side_a.adjusted_rand_index));

    const PartitionAgreement side_b =
        partition_agreement({0, 0, 1, 1}, {0, 1, 2, 3});
    EXPECT_TRUE(std::isnan(side_b.fowlkes_mallows));
    EXPECT_FALSE(std::isnan(side_b.adjusted_rand_index));
}

TEST(PartitionAgreementTest, RuleOneReportsEveryMetricNaN) {
    // Directly: one sample in, one sample out.
    const PartitionAgreement direct = partition_agreement({0}, {0});
    // And through Excluded, which leaves one survivor of three.
    const PartitionAgreement reduced = partition_agreement(
        {0, -1, -1}, {0, 1, 1}, WithNoise(NoiseHandling::Excluded));

    for (const PartitionAgreement& agreement : {direct, reduced}) {
        EXPECT_EQ(agreement.num_samples, 1u);
        EXPECT_TRUE(std::isnan(agreement.adjusted_rand_index));
        EXPECT_TRUE(std::isnan(agreement.fowlkes_mallows));
        EXPECT_TRUE(std::isnan(agreement.normalized_mutual_information));
        EXPECT_TRUE(std::isnan(agreement.homogeneity));
        EXPECT_TRUE(std::isnan(agreement.completeness));
        EXPECT_TRUE(std::isnan(agreement.v_measure));
        EXPECT_TRUE(std::isnan(agreement.adjusted_mutual_information));
    }
}

TEST(PartitionAgreementTest, ScaffoldAgreementFollowsNoiseHandling) {
    const std::vector<ClusterLabel> labels{0, 0, 0, 1, 1, 1, 2, 2, 2};
    const std::vector<std::string> scaffolds{"ar", "ar", "ar", "pi", "",
                                             "al", "al", "al", ""};

    const PartitionAgreement singletons =
        scaffold_agreement(labels, scaffolds, WithNoise(NoiseHandling::Singletons));
    EXPECT_EQ(singletons.num_samples, 9u);
    EXPECT_EQ(singletons.num_clusters_b, 5u);
    EXPECT_NEAR(singletons.adjusted_rand_index, 0.4166666666666667, 1e-12);
    EXPECT_NEAR(singletons.fowlkes_mallows, 0.5443310539518174, 1e-12);

    const PartitionAgreement grouped =
        scaffold_agreement(labels, scaffolds, WithNoise(NoiseHandling::Grouped));
    EXPECT_EQ(grouped.num_clusters_b, 4u);
    EXPECT_NEAR(grouped.adjusted_rand_index, 0.36, 1e-12);
    EXPECT_NEAR(grouped.fowlkes_mallows, 0.5039526306789696, 1e-12);

    const PartitionAgreement excluded =
        scaffold_agreement(labels, scaffolds, WithNoise(NoiseHandling::Excluded));
    EXPECT_EQ(excluded.num_samples, 7u);
    EXPECT_EQ(excluded.num_clusters_b, 3u);
    EXPECT_NEAR(excluded.adjusted_rand_index, 0.631578947368421, 1e-12);
    EXPECT_NEAR(excluded.fowlkes_mallows, 0.7302967433402214, 1e-12);
}

// Field for field, so a forwarder that drops an argument is caught.
TEST(PartitionAgreementTest, ClusteringResultOverloadsForward) {
    ExpectSameAgreement(
        partition_agreement(MakeResult(kMainA), MakeResult(kMainB)),
        partition_agreement(kMainA, kMainB));

    const std::vector<ClusterLabel> labels{0, 0, 1, 1};
    const std::vector<std::string> scaffolds{"ar", "pi", "pi", "pi"};
    ExpectSameAgreement(scaffold_agreement(MakeResult(labels), scaffolds),
                        scaffold_agreement(labels, scaffolds));
}

TEST(PartitionAgreementTest, ValidationRejectsMismatchedAndEmptyInputs) {
    // A braced `{}` is ambiguous between the ClusteringResult and raw-vector
    // overloads -- both are default-constructible -- so name the type.
    const std::vector<ClusterLabel> no_labels;
    const std::vector<std::string> no_scaffolds;

    EXPECT_THROW(partition_agreement({0, 0, 1}, {0, 1}), std::invalid_argument);
    EXPECT_THROW(partition_agreement(no_labels, no_labels),
                 std::invalid_argument);
    EXPECT_THROW(partition_agreement(MakeResult({0, 0, 1}), MakeResult({0, 1})),
                 std::invalid_argument);
    EXPECT_THROW(partition_agreement(ClusteringResult(), ClusteringResult()),
                 std::invalid_argument);
    EXPECT_THROW(scaffold_agreement({0, 0, 1}, {"a", "b"}),
                 std::invalid_argument);
    EXPECT_THROW(scaffold_agreement({0, 0}, no_scaffolds),
                 std::invalid_argument);
    EXPECT_THROW(scaffold_agreement(MakeResult({0, 0, 1}), {"a", "b"}),
                 std::invalid_argument);
    EXPECT_THROW(scaffold_agreement(ClusteringResult(), no_scaffolds),
                 std::invalid_argument);
}

TEST(PartitionAgreementTest, MismatchMessageNamesBothSizes) {
    try {
        scaffold_agreement({0, 0, 1}, {"a", "b"});
        FAIL() << "expected std::invalid_argument";
    } catch (const std::invalid_argument& error) {
        const std::string message = error.what();
        EXPECT_NE(message.find("2"), std::string::npos);
        EXPECT_NE(message.find("3"), std::string::npos);
        EXPECT_NE(message.find("scaffold_labels"), std::string::npos);
    }
}

// sklearn.metrics.normalized_mutual_info_score and
// homogeneity_completeness_v_measure on kMainA/kMainB.
TEST(PartitionAgreementTest, MainFixtureEntropyMetrics) {
    const PartitionAgreement agreement = partition_agreement(kMainA, kMainB);
    EXPECT_NEAR(agreement.normalized_mutual_information, 0.6280760724651604,
                1e-12);
    EXPECT_NEAR(agreement.homogeneity, 0.7103099178571527, 1e-12);
    EXPECT_NEAR(agreement.completeness, 0.5629072918469558, 1e-12);
    EXPECT_NEAR(agreement.v_measure, 0.6280760724651604, 1e-12);
}

TEST(PartitionAgreementTest, NoiseModesEntropyMetrics) {
    const PartitionAgreement singletons =
        partition_agreement(kNoiseA, kNoiseB, WithNoise(NoiseHandling::Singletons));
    EXPECT_NEAR(singletons.normalized_mutual_information, 0.5919578611106004,
                1e-12);
    EXPECT_NEAR(singletons.homogeneity, 0.5997280461117795, 1e-12);
    EXPECT_NEAR(singletons.completeness, 0.5843864446715682, 1e-12);

    const PartitionAgreement grouped =
        partition_agreement(kNoiseA, kNoiseB, WithNoise(NoiseHandling::Grouped));
    EXPECT_NEAR(grouped.normalized_mutual_information, 0.4073100686148156, 1e-12);
    EXPECT_NEAR(grouped.homogeneity, 0.4370927081530442, 1e-12);
    EXPECT_NEAR(grouped.completeness, 0.3813271825750331, 1e-12);

    const PartitionAgreement excluded =
        partition_agreement(kNoiseA, kNoiseB, WithNoise(NoiseHandling::Excluded));
    EXPECT_NEAR(excluded.normalized_mutual_information, 0.561884803785396, 1e-12);
    EXPECT_NEAR(excluded.homogeneity, 0.6329129160661658, 1e-12);
    EXPECT_NEAR(excluded.completeness, 0.505190257900143, 1e-12);
}

TEST(PartitionAgreementTest, NmiAndVMeasureAreBitwiseEqual) {
    // `<utility>` is already in the include block from Task 1.
    const std::vector<std::pair<std::vector<ClusterLabel>,
                                std::vector<ClusterLabel>>> fixtures{
        {kMainA, kMainB},
        {kNoiseA, kNoiseB},
        {{0, 0, 1, 1}, {0, 1, 0, 1}},
        {{0, 0, 0, 0}, {0, 0, 1, 1}},
        {{0, 1, 2, 3}, {0, 0, 1, 1}},
    };
    for (const auto& fixture : fixtures) {
        const PartitionAgreement agreement =
            partition_agreement(fixture.first, fixture.second);
        EXPECT_EQ(agreement.normalized_mutual_information, agreement.v_measure);
    }
}

// The harmonic form 2hc/(h+c) is 0/0 here; the shared-value definition is 0.
TEST(PartitionAgreementTest, IndependentPartitionsGiveZeroNotNaN) {
    const PartitionAgreement agreement =
        partition_agreement({0, 0, 1, 1}, {0, 1, 0, 1});
    EXPECT_DOUBLE_EQ(agreement.homogeneity, 0.0);
    EXPECT_DOUBLE_EQ(agreement.completeness, 0.0);
    EXPECT_DOUBLE_EQ(agreement.v_measure, 0.0);
    EXPECT_DOUBLE_EQ(agreement.normalized_mutual_information, 0.0);
    EXPECT_NEAR(agreement.adjusted_rand_index, -0.5, 1e-12);
    EXPECT_DOUBLE_EQ(agreement.fowlkes_mallows, 0.0);
}

// The harmonic form inherits homogeneity's NaN here; the composite is defined.
TEST(PartitionAgreementTest, ZeroEntropySideLeavesTheCompositeDefined) {
    const PartitionAgreement agreement =
        partition_agreement({0, 0, 0, 0}, {0, 0, 1, 1});
    EXPECT_TRUE(std::isnan(agreement.homogeneity));
    EXPECT_DOUBLE_EQ(agreement.completeness, 0.0);
    EXPECT_DOUBLE_EQ(agreement.v_measure, 0.0);
    EXPECT_DOUBLE_EQ(agreement.normalized_mutual_information, 0.0);
}

TEST(PartitionAgreementTest, ZeroEntropySideBIsCompletenessNaN) {
    const PartitionAgreement agreement =
        partition_agreement({0, 0, 1, 1}, {0, 0, 0, 0});
    EXPECT_DOUBLE_EQ(agreement.homogeneity, 0.0);
    EXPECT_TRUE(std::isnan(agreement.completeness));
    EXPECT_DOUBLE_EQ(agreement.v_measure, 0.0);
}

TEST(PartitionAgreementTest, HomogeneityAndCompletenessSwapWithTheArguments) {
    const PartitionAgreement forward = partition_agreement(kMainA, kMainB);
    const PartitionAgreement backward = partition_agreement(kMainB, kMainA);
    EXPECT_DOUBLE_EQ(forward.homogeneity, backward.completeness);
    EXPECT_DOUBLE_EQ(forward.completeness, backward.homogeneity);
    EXPECT_DOUBLE_EQ(forward.normalized_mutual_information,
                     backward.normalized_mutual_information);
    EXPECT_DOUBLE_EQ(forward.v_measure, backward.v_measure);
}

// kMainA/kMainB happens to be symmetric enough that an order-dependent
// accumulation still swaps exactly. This contingency table -- [[2, 7],
// [5, 1]], marginals 9/6 against 7/8 -- is not, and the MI term ordering
// must be invariant under transposition or these values differ in the last
// bit. The assertion is bitwise because the difference is a single ULP, which
// EXPECT_DOUBLE_EQ's 4-ULP tolerance would hide.
TEST(PartitionAgreementTest, ArgumentSwapIsExactOnAsymmetricMarginals) {
    const std::vector<ClusterLabel> a{0, 0, 0, 0, 0, 0, 0, 0, 0,
                                      1, 1, 1, 1, 1, 1};
    const std::vector<ClusterLabel> b{0, 0, 1, 1, 1, 1, 1, 1, 1,
                                      0, 0, 0, 0, 0, 1};
    ASSERT_EQ(a.size(), b.size());
    ASSERT_EQ(a.size(), 15u);

    const PartitionAgreement forward = partition_agreement(a, b);
    const PartitionAgreement backward = partition_agreement(b, a);
    ExpectSameDouble(forward.homogeneity, backward.completeness);
    ExpectSameDouble(forward.completeness, backward.homogeneity);
    ExpectSameDouble(forward.normalized_mutual_information,
                     backward.normalized_mutual_information);
    ExpectSameDouble(forward.v_measure, backward.v_measure);
}

TEST(PartitionAgreementTest, ScaffoldEntropyMetrics) {
    const std::vector<ClusterLabel> labels{0, 0, 0, 1, 1, 1, 2, 2, 2};
    const std::vector<std::string> scaffolds{"ar", "ar", "ar", "pi", "",
                                             "al", "al", "al", ""};
    const PartitionAgreement singletons = scaffold_agreement(labels, scaffolds);
    EXPECT_NEAR(singletons.normalized_mutual_information, 0.6916056673469443,
                1e-12);
    EXPECT_NEAR(singletons.homogeneity, 0.8068732785714351, 1e-12);
    EXPECT_NEAR(singletons.completeness, 0.6051549589285762, 1e-12);
}

// The unordered_map traversal is reproducible inside one build, so scoring the
// same input twice cannot fail. Two orderings of the same table can: interned
// ids follow first appearance, so permuting the samples renumbers the clusters
// and reorders every id-driven traversal. This fixture was verified to differ
// in the last bits when the MI and entropy accumulations are left unsorted, so
// it exercises the canonical summation order requirement.
TEST(PartitionAgreementTest, PermutedInputsGiveBitwiseEqualResults) {
    std::vector<ClusterLabel> a;
    std::vector<ClusterLabel> b;
    int index = 0;
    for (int cluster = 0; cluster < 12; ++cluster) {
        // Unequal cluster sizes and unequal cell counts are what make the
        // accumulation order observable: a fixture whose rows are all the
        // same size sums to the same double in any order.
        const int size = 4 + cluster % 3;
        for (int member = 0; member < size; ++member) {
            a.push_back(cluster);
            b.push_back(index % 7);
            ++index;
        }
    }

    // A stride permutation: same contingency table, different sample order.
    std::vector<ClusterLabel> permuted_a;
    std::vector<ClusterLabel> permuted_b;
    for (size_t offset = 0; offset < 5; ++offset) {
        for (size_t i = offset; i < a.size(); i += 5) {
            permuted_a.push_back(a[i]);
            permuted_b.push_back(b[i]);
        }
    }
    ASSERT_EQ(permuted_a.size(), a.size());

    // Every field, AMI included: the spec asks for bitwise equality across the
    // whole struct. AMI is still NaN on both sides at this task, so
    // ExpectSameAgreement passes it now and compares the real values from
    // Task 4 on without this test needing to change.
    PartitionAgreementOptions options;
    options.compute_adjusted_mutual_information = true;
    ExpectSameAgreement(partition_agreement(permuted_a, permuted_b, options),
                        partition_agreement(a, b, options));
}

// sklearn.metrics.adjusted_mutual_info_score on kMainA/kMainB.
TEST(PartitionAgreementTest, MainFixtureAdjustedMutualInformation) {
    const PartitionAgreement agreement =
        partition_agreement(kMainA, kMainB, WithAmi());
    EXPECT_TRUE(agreement.requested.adjusted_mutual_information);
    EXPECT_NEAR(agreement.adjusted_mutual_information, 0.47492716044734345,
                1e-12);
}

// A sparse-only E[MI] passes everything else and fails this: 21 of the 30
// (i, j) pairs have a zero observed count and still contribute.
TEST(PartitionAgreementTest, AdjustedMutualInformationIteratesZeroCells) {
    const PartitionAgreement agreement =
        partition_agreement(kAmiA, kAmiB, WithAmi());
    EXPECT_EQ(agreement.num_samples, 18u);
    EXPECT_EQ(agreement.num_clusters_a, 6u);
    EXPECT_EQ(agreement.num_clusters_b, 5u);
    EXPECT_NEAR(agreement.adjusted_mutual_information, 0.5531261900574347,
                1e-12);
    // The always-computed metrics on the same fixture, for good measure.
    EXPECT_NEAR(agreement.adjusted_rand_index, 0.5225144895229603, 1e-12);
    EXPECT_NEAR(agreement.normalized_mutual_information, 0.7261677779961542,
                1e-12);
}

// The grouping of equal marginal values is exact in exact arithmetic but not
// bitwise identical to the naive loop: the grouped form multiplies once where
// the naive form adds repeatedly, and IEEE-754 addition is not associative.
TEST(PartitionAgreementTest, MarginalGroupingMatchesTheNaiveLoop) {
    const detail::ContingencyTable table =
        detail::build_contingency(kAmiA, kAmiB, NoiseHandling::Singletons);

    const std::vector<double> logfact =
        detail::log_factorials(table.num_samples);
    const double grouped =
        detail::expected_mutual_information(table, logfact);

    const double naive = NaiveExpectedMutualInformation(table);
    EXPECT_NEAR(grouped, naive, 1e-12 * std::abs(naive));
}

TEST(PartitionAgreementTest, AmiIsNaNWhenNobodyAsked) {
    const PartitionAgreement agreement = partition_agreement(kMainA, kMainB);
    EXPECT_FALSE(agreement.requested.adjusted_mutual_information);
    EXPECT_TRUE(std::isnan(agreement.adjusted_mutual_information));
}

TEST(PartitionAgreementTest, AmiRequestedOnDegenerateInputIsNaNAndRequested) {
    const PartitionAgreement agreement = partition_agreement({0}, {0}, WithAmi());
    EXPECT_TRUE(agreement.requested.adjusted_mutual_information);
    EXPECT_TRUE(std::isnan(agreement.adjusted_mutual_information));
}

TEST(PartitionAgreementTest, IdenticalPartitionsReportAmiOnlyWhenRequested) {
    // The three shapes SelfAgreementIsOneEverywhere uses, so the spec's
    // self-agreement requirement covers AMI on all of them. The degenerate two
    // matter most: an all-singletons or single-cluster partition against itself
    // is where a chance-corrected metric computed rather than short-circuited
    // would return 0/0, and rule 2 is what keeps them at 1.0.
    const std::vector<std::vector<ClusterLabel>> fixtures{
        {0, 0, 1, 1, 2, 2},
        {0, 1, 2, 3, 4, 5},
        {0, 0, 0, 0, 0, 0},
    };
    for (const std::vector<ClusterLabel>& labels : fixtures) {
        EXPECT_TRUE(std::isnan(partition_agreement(labels, labels)
                                   .adjusted_mutual_information));
        EXPECT_DOUBLE_EQ(partition_agreement(labels, labels, WithAmi())
                             .adjusted_mutual_information,
                         1.0);
    }
}

// A nearly independent table cancels to a true MI far below the roundoff of
// the terms that produced it, so without a clamp the sum lands slightly
// negative and carries the four metrics derived from it below zero. This 2x2
// table is the smallest reproduction found: the exact MI is +3.1e-18 and the
// unclamped sum is -2.7e-17. scikit-learn reports 0.0 on the same input.
TEST(PartitionAgreementTest, NearIndependentTableDoesNotGoNegative) {
    struct Cell {
        ClusterLabel a;
        ClusterLabel b;
        int count;
    };
    const std::vector<Cell> cells{
        {0, 0, 10000}, {0, 1, 9999}, {1, 0, 10001}, {1, 1, 10000}};
    std::vector<ClusterLabel> a;
    std::vector<ClusterLabel> b;
    for (const Cell& cell : cells) {
        a.insert(a.end(), cell.count, cell.a);
        b.insert(b.end(), cell.count, cell.b);
    }

    const PartitionAgreement agreement = partition_agreement(a, b, WithAmi());
    EXPECT_GE(agreement.normalized_mutual_information, 0.0);
    EXPECT_GE(agreement.homogeneity, 0.0);
    EXPECT_GE(agreement.completeness, 0.0);
    EXPECT_GE(agreement.v_measure, 0.0);
    EXPECT_NEAR(agreement.normalized_mutual_information, 0.0, 1e-15);
    EXPECT_TRUE(std::isfinite(agreement.adjusted_mutual_information));
}

// The smallest AMI denominator this fixture's eight samples can reach, 0.0866.
// The clamp itself is defensive: the denominator is bounded below by log(2)/N
// over every input rule 2 does not intercept, so eps would take about 3e15
// samples. The assertion is finite and near zero, as the spec states.
TEST(PartitionAgreementTest, AmiStaysFiniteWhenExpectedMiConsumesTheNormalizer) {
    const PartitionAgreement agreement = partition_agreement(
        {0, 1, 2, 3, 4, 5, 6, 7}, {0, 0, 1, 2, 3, 4, 5, 6}, WithAmi());
    EXPECT_TRUE(std::isfinite(agreement.adjusted_mutual_information));
    EXPECT_NEAR(agreement.adjusted_mutual_information, 0.0, 1e-9);
}

TEST(PartitionAgreementTest, AmiIsSymmetric) {
    EXPECT_NEAR(partition_agreement(kMainA, kMainB, WithAmi())
                    .adjusted_mutual_information,
                partition_agreement(kMainB, kMainA, WithAmi())
                    .adjusted_mutual_information,
                1e-12);
}

TEST(PartitionAgreementTest, AmiUnderTheThreeNoiseModes) {
    // WithAmi() already defaults noise_handling to Singletons.
    EXPECT_NEAR(partition_agreement(kNoiseA, kNoiseB, WithAmi())
                    .adjusted_mutual_information,
                -0.012595285524420317, 1e-12);

    PartitionAgreementOptions grouped = WithAmi();
    grouped.noise_handling = NoiseHandling::Grouped;
    EXPECT_NEAR(partition_agreement(kNoiseA, kNoiseB, grouped)
                    .adjusted_mutual_information,
                -0.08744260969263824, 1e-12);

    PartitionAgreementOptions excluded = WithAmi();
    excluded.noise_handling = NoiseHandling::Excluded;
    EXPECT_NEAR(partition_agreement(kNoiseA, kNoiseB, excluded)
                    .adjusted_mutual_information,
                0.09605662055753134, 1e-12);
}
