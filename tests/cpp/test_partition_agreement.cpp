#include <gtest/gtest.h>

#include <cstdint>
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
