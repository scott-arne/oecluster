#include <gtest/gtest.h>

#include <chrono>
#include <cmath>
#include <condition_variable>
#include <cstddef>
#include <limits>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <string>
#include <vector>

#include "oecluster/StorageBackend.h"
#include "oecluster/clustering/ClusterTypes.h"

#include "../../src/clustering/InternalIndices.h"
#include "../../src/clustering/ReportDistanceSource.h"
#include "diversity_test_support.h"

using namespace OECluster;
using namespace OECluster::detail;
using diversity_test::Condensed;
using diversity_test::Hashed;
using diversity_test::IsolationComparison;
using diversity_test::MakeStorage;
using diversity_test::TableComparison;

namespace {

// Drains a plan through a source, returning (anchor, target, value) triples.
struct Read {
    size_t anchor;
    size_t target;
    double value;
};

template <class Source, class Plan>
std::vector<Read> Drain(Source& source, Plan plan) {
    Plan walk = plan;
    auto reader = source.Open(plan);
    std::vector<Read> out;
    ReportRow row{};
    while (walk.Next(row)) {
        for (size_t k = 0; k < row.count; ++k) {
            out.push_back({row.anchor, row.targets[k], reader.Next(row.anchor, row.targets[k])});
        }
    }
    reader.Finish();
    return out;
}

const Clusters CLUSTERS = {{0, 3, 5}, {1, 2}, {4}, {6, 7, 8, 9}};

// Throws at chosen canonical pairs; otherwise a table.
class ThrowingComparison : public PairwiseComparison {
public:
    ThrowingComparison(size_t n, std::vector<double> condensed,
                       std::vector<std::pair<size_t, size_t>> bad)
        : table_(n, std::move(condensed)), bad_(std::move(bad)) {}
    double Compare(size_t i, size_t j) override {
        for (const auto& pair : bad_) {
            if (pair.first == i && pair.second == j) {
                throw std::runtime_error("bad pair " + std::to_string(i) + "," +
                                         std::to_string(j));
            }
        }
        return table_.Compare(i, j);
    }
    std::unique_ptr<PairwiseComparison> Clone() const override {
        return std::make_unique<ThrowingComparison>(*this);
    }
    size_t Size() const override { return table_.Size(); }
    std::string ComparisonName() const override { return "throwing"; }

private:
    TableComparison table_;
    std::vector<std::pair<size_t, size_t>> bad_;
};

// Forces the adverse completion order: the earlier canonical bad pair waits
// until the later one has thrown (or two seconds pass, so a serial run cannot
// deadlock), then throws itself. Clones share the gate.
class GatedComparison : public PairwiseComparison {
public:
    GatedComparison(size_t n, std::vector<double> condensed,
                    std::pair<size_t, size_t> early, std::pair<size_t, size_t> late)
        : table_(n, std::move(condensed)), early_(early), late_(late),
          gate_(std::make_shared<Gate>()) {}
    double Compare(size_t i, size_t j) override {
        if (i == late_.first && j == late_.second) {
            {
                const std::lock_guard<std::mutex> lock(gate_->mutex);
                gate_->late_thrown = true;
            }
            gate_->cv.notify_all();
            throw std::runtime_error("bad pair " + std::to_string(i) + "," +
                                     std::to_string(j));
        }
        if (i == early_.first && j == early_.second) {
            std::unique_lock<std::mutex> lock(gate_->mutex);
            gate_->late_first = gate_->cv.wait_for(lock, std::chrono::seconds(2),
                                                   [&] { return gate_->late_thrown; });
            throw std::runtime_error("bad pair " + std::to_string(i) + "," +
                                     std::to_string(j));
        }
        return table_.Compare(i, j);
    }
    std::unique_ptr<PairwiseComparison> Clone() const override {
        return std::make_unique<GatedComparison>(*this);
    }
    size_t Size() const override { return table_.Size(); }
    std::string ComparisonName() const override { return "gated"; }

    bool LateFirst() const {
        const std::lock_guard<std::mutex> lock(gate_->mutex);
        return gate_->late_first;
    }

private:
    struct Gate {
        std::mutex mutex;
        std::condition_variable cv;
        bool late_thrown = false;
        bool late_first = false;
    };
    TableComparison table_;
    std::pair<size_t, size_t> early_;
    std::pair<size_t, size_t> late_;
    std::shared_ptr<Gate> gate_;
};

}  // namespace

TEST(ReportDistanceSourceTest, FiniteReportDistanceKeepsTheMessage) {
    EXPECT_EQ(finite_report_distance(0.5, 1, 2), 0.5);
    EXPECT_EQ(finite_report_distance(-0.5, 1, 2), -0.5);
    try {
        finite_report_distance(std::numeric_limits<double>::infinity(), 3, 7);
        FAIL();
    } catch (const std::invalid_argument& e) {
        EXPECT_STREQ(e.what(), "cluster_report: distance between samples 3 and 7 is not finite");
    }
}

TEST(ReportDistanceSourceTest, PlansYieldTheDocumentedRows) {
    IntraRows intra(CLUSTERS, 0, CLUSTERS.size());
    ReportRow row{};
    std::vector<std::pair<size_t, size_t>> rows;
    while (intra.Next(row)) {
        rows.emplace_back(row.anchor, row.count);
    }
    const std::vector<std::pair<size_t, size_t>> expected = {
        {0, 2}, {3, 1}, {1, 1}, {6, 3}, {7, 2}, {8, 1}};
    EXPECT_EQ(rows, expected);

    CrossRows cross(CLUSTERS);
    size_t cross_pairs = 0;
    while (cross.Next(row)) {
        cross_pairs += row.count;
    }
    EXPECT_EQ(cross_pairs, 3u * 2 + 3 * 1 + 3 * 4 + 2 * 1 + 2 * 4 + 1 * 4);

    const std::vector<size_t> none;
    AllPointsRows empty(10, none);
    EXPECT_FALSE(empty.Next(row));
}

TEST(ReportDistanceSourceTest, ComparisonMatchesMatrixForEveryPlanThreadAndChunk) {
    const size_t n = 10;
    const std::vector<double> condensed = Hashed(n);
    const DenseStorage storage = MakeStorage(n, condensed);
    MatrixSource matrix(storage);
    const std::vector<size_t> reps = {3, 2, 4, 7};
    const std::vector<size_t> medoids = {0, 1, 4, 8};
    const size_t target = 5;

    for (const size_t threads : {size_t{1}, size_t{2}, size_t{4}, size_t{0}}) {
        for (const size_t chunk : {size_t{1}, size_t{3}, size_t{4096}, SIZE_MAX}) {
            for (const size_t block : {size_t{1}, size_t{5}, FILL_BLOCK_DISTANCES}) {
                TableComparison table(n, condensed);
                ComparisonSource lazy(table, threads, chunk, block);
                auto same = [&](auto plan) {
                    const std::vector<Read> a = Drain(matrix, plan);
                    const std::vector<Read> b = Drain(lazy, plan);
                    ASSERT_EQ(a.size(), b.size());
                    for (size_t k = 0; k < a.size(); ++k) {
                        EXPECT_EQ(a[k].value, b[k].value);
                    }
                };
                same(IntraRows(CLUSTERS, 0, CLUSTERS.size()));
                same(CrossRows(CLUSTERS));
                same(MemberRows(CLUSTERS, reps, &medoids));
                same(MemberRows(CLUSTERS, reps, nullptr));
                same(PeerRows(reps));
                same(FixedTargetRows(medoids, &target));
                same(AllPointsRows(n, reps));
            }
        }
    }
}

TEST(ReportDistanceSourceTest, BlocksRespectTheCap) {
    const size_t n = 10;
    TableComparison table(n, Hashed(n));
    for (const size_t cap : {size_t{1}, size_t{2}, size_t{5}, size_t{7}}) {
        size_t worst = 0;
        size_t blocks = 0;
        ComparisonSource lazy(table, 2, 1, cap, [&](size_t distances, size_t rows) {
            EXPECT_GE(rows, 1u);
            worst = std::max(worst, distances);
            ++blocks;
        });
        Drain(lazy, CrossRows(CLUSTERS));
        // The longest cross row is 4 (cluster 3), so a block may exceed a
        // smaller cap only by holding that one row.
        EXPECT_LE(worst, std::max<size_t>(cap, 4));
        EXPECT_GT(blocks, 1u);
    }
}

TEST(ReportDistanceSourceTest, EarliestPairWinsWhateverFilledFirst) {
    const size_t n = 10;
    // Two bad pairs in cross order: (0, 1) is the first cross pair, (5, 9) --
    // cluster 0's last member against cluster 3 -- much later and in a
    // different unit.
    for (const size_t threads : {size_t{1}, size_t{2}, size_t{4}, size_t{0}}) {
        for (const size_t chunk : {size_t{1}, size_t{2}}) {
            ThrowingComparison throwing(n, Hashed(n), {{5, 9}, {0, 1}});
            ComparisonSource lazy(throwing, threads, chunk, SIZE_MAX);
            try {
                Drain(lazy, CrossRows(CLUSTERS));
                FAIL();
            } catch (const std::runtime_error& e) {
                EXPECT_STREQ(e.what(), "bad pair 0,1");
            }
        }
    }
}

TEST(ReportDistanceSourceTest, EarliestPairWinsWhenTheLaterOneCompletesFirst) {
    const size_t n = 10;
    for (const size_t threads : {size_t{1}, size_t{2}, size_t{4}}) {
        GatedComparison gated(n, Hashed(n), {0, 1}, {5, 9});
        ComparisonSource lazy(gated, threads, 1, SIZE_MAX);
        try {
            Drain(lazy, CrossRows(CLUSTERS));
            FAIL();
        } catch (const std::runtime_error& e) {
            EXPECT_STREQ(e.what(), "bad pair 0,1") << "threads " << threads;
        }
        if (threads > 1) {
            // Proves the adverse order was actually staged, not just allowed.
            EXPECT_TRUE(gated.LateFirst()) << "threads " << threads;
        }
    }
}

TEST(ReportDistanceSourceTest, ClonesAreIsolatedAndBoundedByTheItems) {
    const size_t n = 10;
    const std::vector<double> condensed = Hashed(n);
    const DenseStorage storage = MakeStorage(n, condensed);
    MatrixSource matrix(storage);
    IsolationComparison isolation(n, condensed);
    ComparisonSource lazy(isolation, 4, 1, SIZE_MAX);
    const std::vector<Read> expected = Drain(matrix, CrossRows(CLUSTERS));
    const std::vector<Read> got = Drain(lazy, CrossRows(CLUSTERS));
    ASSERT_EQ(expected.size(), got.size());
    for (size_t k = 0; k < got.size(); ++k) {
        EXPECT_EQ(expected[k].value, got[k].value);
    }
    EXPECT_EQ(isolation.Violations(), 0u);
    EXPECT_TRUE(isolation.OverlapObserved())
        << "Overlap not observed; test may be flaky on this machine";

    // More workers requested than items, and the automatic count: the cross
    // plan has more units than items, so only the item cap bounds the clones.
    for (const size_t threads : {size_t{64}, size_t{0}}) {
        TableComparison table(n, condensed);
        ComparisonSource capped(table, threads, 1, SIZE_MAX);
        Drain(capped, CrossRows(CLUSTERS));
        EXPECT_GE(table.NumClones(), 1u);
        EXPECT_LE(table.NumClones(), n) << "threads " << threads;
    }
}

TEST(ReportDistanceSourceTest, NonFiniteSurfacesAtItsPairWithTheReportMessage) {
    const size_t n = 10;
    std::vector<double> condensed = Hashed(n);
    // Pair (2, 3) at condensed index n*i - i*(i+1)/2 + j - i - 1; the cross
    // pass reads it as anchor 3 against target 2.
    const size_t i = 2;
    const size_t j = 3;
    condensed[n * i - i * (i + 1) / 2 + j - i - 1] = std::numeric_limits<double>::quiet_NaN();
    const DenseStorage storage = MakeStorage(n, condensed);
    MatrixSource matrix(storage);
    TableComparison table(n, condensed);
    ComparisonSource lazy(table, 4, 1, SIZE_MAX);
    std::string matrix_message;
    std::string lazy_message;
    try {
        Drain(matrix, CrossRows(CLUSTERS));
    } catch (const std::invalid_argument& e) {
        matrix_message = e.what();
    }
    try {
        Drain(lazy, CrossRows(CLUSTERS));
    } catch (const std::invalid_argument& e) {
        lazy_message = e.what();
    }
    EXPECT_FALSE(matrix_message.empty());
    EXPECT_EQ(matrix_message, lazy_message);
}

TEST(ReportDistanceSourceTest, ReaderContractViolationsAreLogicErrors) {
    const size_t n = 10;
    TableComparison table(n, Hashed(n));
    ComparisonSource lazy(table, 1, 4096);
    {
        auto reader = lazy.Open(IntraRows(CLUSTERS, 0, 1));
        EXPECT_THROW(reader.Next(3, 5), std::logic_error);  // the first pair is (0, 3)
    }
    {
        auto reader = lazy.Open(IntraRows(CLUSTERS, 0, 1));
        reader.Next(0, 3);
        EXPECT_THROW(reader.Finish(), std::logic_error);
    }
    {
        auto reader = lazy.Open(IntraRows(CLUSTERS, 2, 3));  // a singleton: no rows
        EXPECT_THROW(reader.Next(4, 4), std::logic_error);
        reader.Finish();
    }
}

TEST(ReportDistanceSourceTest, NeverReadsReversedOrDiagonalPairs) {
    const size_t n = 10;
    TableComparison table(n, Hashed(n));  // throws logic_error on either
    ComparisonSource lazy(table, 2, 3);
    const std::vector<size_t> reps = {9, 2, 4, 6};
    EXPECT_NO_THROW(Drain(lazy, PeerRows(reps)));
    EXPECT_NO_THROW(Drain(lazy, MemberRows(CLUSTERS, reps, nullptr)));
    EXPECT_NO_THROW(Drain(lazy, AllPointsRows(n, reps)));
}
