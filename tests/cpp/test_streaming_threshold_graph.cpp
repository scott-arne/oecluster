/**
 * @file test_streaming_threshold_graph.cpp
 * @brief The threshold graph built from a comparison, and its memory guard.
 */

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include <oechem.h>

#include "oecluster/StorageBackend.h"
#include "oecluster/comparisons/DescriptorComparison.h"
#include "oecluster/comparisons/FingerprintComparison.h"

#include "diversity_test_support.h"
#include "streaming_test_support.h"

using namespace OECluster;
using namespace diversity_test;
using namespace streaming_test;

static_assert(sizeof(size_t) == 8,
              "the guard figures below assume a 64-bit size_t");

namespace {

const double NaN = std::numeric_limits<double>::quiet_NaN();
const double INF = std::numeric_limits<double>::infinity();
constexpr size_t ONE_GIB = size_t{1} << 30;

ThresholdGraphOptions GraphOptions(double threshold, size_t threads = 1,
                                   size_t chunk = 4096) {
    ThresholdGraphOptions options;
    options.threshold = threshold;
    options.num_threads = threads;
    options.chunk_size = chunk;
    return options;
}

std::vector<std::vector<size_t>> StorageRows(size_t n,
                                             const std::vector<double>& condensed,
                                             double threshold) {
    return Rows(BuildThresholdNeighborGraph(MakeStorage(n, condensed),
                                            GraphOptions(threshold)));
}

void ExpectMessage(const std::function<void()>& call, const std::string& kind,
                   const std::string& message) {
    try {
        call();
        FAIL() << "expected " << kind << ": " << message;
    } catch (const std::length_error& error) {
        EXPECT_EQ(kind, "length_error");
        EXPECT_EQ(std::string(error.what()), message);
    } catch (const std::logic_error& error) {
        EXPECT_EQ(kind, "logic_error");
        EXPECT_EQ(std::string(error.what()), message);
    } catch (const std::runtime_error& error) {
        EXPECT_EQ(kind, "runtime_error");
        EXPECT_EQ(std::string(error.what()), message);
    }
}

std::string LimitMessage(const std::string& caller, size_t n, size_t edges,
                         const std::string& kind, size_t limit) {
    return caller + " would build a threshold graph of " +
           std::to_string(detail::threshold_graph_bytes(n, edges)) +
           " bytes for " + std::to_string(n) + " items and " +
           std::to_string(edges) + " edges, above its " + kind + " of " +
           std::to_string(limit) +
           " bytes; use a tighter threshold, a larger max_graph_bytes, or a "
           "memory-mapped matrix from pdist(output=...)";
}

std::vector<OEChem::OEGraphMol> Molecules() {
    const char* smiles[] = {"CCO",       "CCCO",       "CCCCO",     "c1ccccc1",
                            "Cc1ccccc1", "CCc1ccccc1", "CC(=O)O",   "CC(=O)OC",
                            "CCN",       "CCCN",       "C1CCCCC1",  "c1ccncc1",
                            "CCCCCO",    "c1ccc2ccccc2c1", "CC(C)O", "OCCO"};
    std::vector<OEChem::OEGraphMol> mols;
    for (const char* smi : smiles) {
        mols.emplace_back();
        EXPECT_TRUE(OEChem::OESmilesToMol(mols.back(), smi)) << smi;
    }
    return mols;
}

std::vector<OEChem::OEMolBase*> Pointers(std::vector<OEChem::OEGraphMol>& mols) {
    std::vector<OEChem::OEMolBase*> pointers;
    for (auto& mol : mols) {
        pointers.push_back(&static_cast<OEChem::OEMolBase&>(mol));
    }
    return pointers;
}

// The median of a matrix's distances, so a threshold admits about half the
// pairs whatever the comparison's scale.
double MedianDistance(const DenseStorage& storage) {
    std::vector<double> values(storage.Data(),
                               storage.Data() + storage.NumPairs());
    std::nth_element(values.begin(), values.begin() + values.size() / 2,
                     values.end());
    return values[values.size() / 2];
}

}  // namespace

TEST(StreamingThresholdGraphTest, MatchesTheStorageBuilderOnTables) {
    struct Case {
        size_t n;
        std::vector<double> condensed;
    };
    const std::vector<Case> cases{{2, Scrambled(2)},
                                  {7, Scrambled(7)},
                                  {12, ScrambledSixths(12)},
                                  {25, Hashed(25)},
                                  {40, Quantized(40, 3, 5)}};
    for (const Case& fixture : cases) {
        const size_t pairs = fixture.n * (fixture.n - 1) / 2;
        for (const double threshold : {0.0, 0.2, 0.5, 2.5, 4.0}) {
            const auto expected = StorageRows(fixture.n, fixture.condensed,
                                              threshold);
            for (const size_t threads : {size_t{1}, size_t{2}, size_t{8}}) {
                for (const size_t chunk : {size_t{1}, size_t{7}, pairs + 1}) {
                    SCOPED_TRACE(::testing::Message()
                                 << "n " << fixture.n << ", threshold "
                                 << threshold << ", threads " << threads
                                 << ", chunk " << chunk);
                    TableComparison table(fixture.n, fixture.condensed);
                    EXPECT_EQ(Rows(BuildThresholdNeighborGraph(
                                  table, GraphOptions(threshold, threads, chunk))),
                              expected);
                }
            }
        }
    }
}

TEST(StreamingThresholdGraphTest, MatchesTheStorageBuilderOnRealComparisons) {
    std::vector<OEChem::OEGraphMol> mols = Molecules();
    const std::vector<OEChem::OEMolBase*> pointers = Pointers(mols);
    FingerprintComparison fingerprint(pointers);
    DescriptorComparison descriptor(pointers);
    for (PairwiseComparison* comparison :
         {static_cast<PairwiseComparison*>(&fingerprint),
          static_cast<PairwiseComparison*>(&descriptor)}) {
        const DenseStorage storage = CompareFilled(*comparison);
        const double threshold = MedianDistance(storage);
        const auto expected =
            Rows(BuildThresholdNeighborGraph(storage, GraphOptions(threshold)));
        const size_t pairs = storage.NumPairs();
        for (const size_t threads : {size_t{1}, size_t{2}, size_t{8}}) {
            for (const size_t chunk : {size_t{1}, size_t{7}, pairs + 1}) {
                SCOPED_TRACE(::testing::Message()
                             << comparison->ComparisonName() << ", threads "
                             << threads << ", chunk " << chunk);
                EXPECT_EQ(Rows(BuildThresholdNeighborGraph(
                              *comparison,
                              GraphOptions(threshold, threads, chunk))),
                          expected);
            }
        }
    }
}

TEST(StreamingThresholdGraphTest, ZeroOneAndTwoItems) {
    TableComparison empty(0, {});
    EXPECT_EQ(BuildThresholdNeighborGraph(empty, GraphOptions(1.0)).Size(), 0u);
    TableComparison single(1, {});
    EXPECT_EQ(Rows(BuildThresholdNeighborGraph(single, GraphOptions(1.0))),
              (std::vector<std::vector<size_t>>{{0}}));
    TableComparison pair(2, {0.5});
    EXPECT_EQ(Rows(BuildThresholdNeighborGraph(pair, GraphOptions(0.5))),
              (std::vector<std::vector<size_t>>{{0, 1}, {0, 1}}));
    EXPECT_EQ(Rows(BuildThresholdNeighborGraph(pair, GraphOptions(0.4))),
              (std::vector<std::vector<size_t>>{{0}, {1}}));
}

TEST(StreamingThresholdGraphTest, ThresholdsAdmittingNoPairsAndEveryPair) {
    const size_t n = 6;
    TableComparison table(n, Line(n));
    const auto none = Rows(BuildThresholdNeighborGraph(table, GraphOptions(0.5)));
    const auto every = Rows(BuildThresholdNeighborGraph(table, GraphOptions(10.0)));
    for (size_t i = 0; i < n; ++i) {
        EXPECT_EQ(none[i], std::vector<size_t>{i});
        EXPECT_EQ(every[i], (std::vector<size_t>{0, 1, 2, 3, 4, 5}));
    }
}

TEST(StreamingThresholdGraphTest, RefusesANegativeThreshold) {
    TableComparison table(3, Line(3));
    try {
        BuildThresholdNeighborGraph(table, GraphOptions(-0.1));
        FAIL() << "expected std::invalid_argument";
    } catch (const std::invalid_argument& error) {
        EXPECT_EQ(std::string(error.what()),
                  "ThresholdGraph threshold cannot be negative");
    }
}

TEST(StreamingThresholdGraphTest, ZeroChunkSizeUsesTheDefault) {
    const size_t n = 25;
    const std::vector<double> condensed = Hashed(n);
    TableComparison table(n, condensed);
    const auto expected = StorageRows(n, condensed, 0.5);
    EXPECT_EQ(Rows(BuildThresholdNeighborGraph(table, GraphOptions(0.5, 4, 0))),
              expected);
    // Not vacuous: the fixture's graph has rows beyond the item itself.
    EXPECT_GT(expected[0].size(), 1u);
}

// chunk_size 1 keeps the scoring on the ThreadPool path.
TEST(StreamingThresholdGraphTest, CapsAnAbsurdThreadCount) {
    const size_t n = 9;
    const std::vector<double> condensed =
        Positions({0, 1, 3, 7, 8, 12, 13, 20, 21});
    TableComparison table(n, condensed);
    EXPECT_EQ(Rows(BuildThresholdNeighborGraph(
                  table, GraphOptions(1.5, std::size_t{1} << 61, 1))),
              StorageRows(n, condensed, 1.5));
}

TEST(StreamingThresholdGraphTest, ProvesCloneIsolation) {
    const size_t n = 40;
    const std::vector<double> condensed = Scrambled(n);
    IsolationComparison comparison(n, condensed);
    const auto rows =
        Rows(BuildThresholdNeighborGraph(comparison, GraphOptions(2.0, 4, 1)));
    EXPECT_EQ(comparison.Violations(), 0u);
    EXPECT_TRUE(comparison.OverlapObserved())
        << "Overlap not observed; test may be flaky on this machine";
    EXPECT_EQ(rows, StorageRows(n, condensed, 2.0));
}

TEST(ThresholdGraphSizeTest, GraphBytesFollowTheFormula) {
    EXPECT_EQ(detail::threshold_graph_bytes(0, 0), sizeof(size_t));
    EXPECT_EQ(detail::threshold_graph_bytes(5, 4), sizeof(size_t) * 19);
    EXPECT_EQ(detail::threshold_graph_bytes(100, 4950), size_t{80808});
    EXPECT_EQ(detail::threshold_graph_bytes(12000, 71994000),
              size_t{1152096008});
}

TEST(ThresholdGraphSizeTest, DefaultLimitIsTheMatrixOrOneGibibyte) {
    EXPECT_EQ(detail::default_threshold_graph_limit(0), ONE_GIB);
    EXPECT_EQ(detail::default_threshold_graph_limit(1), ONE_GIB);
    EXPECT_EQ(detail::default_threshold_graph_limit(100), ONE_GIB);
    EXPECT_EQ(detail::default_threshold_graph_limit(12000), ONE_GIB);
    EXPECT_EQ(detail::default_threshold_graph_limit(16500),
              sizeof(double) * size_t{136116750});
}

TEST(ThresholdGraphSizeTest, AnExplicitBudgetReplacesTheDefault) {
    EXPECT_EQ(detail::threshold_graph_limit(100, 0),
              detail::default_threshold_graph_limit(100));
    // Below the default and above it: the budget is not a cap on the default.
    EXPECT_EQ(detail::threshold_graph_limit(100, 5), size_t{5});
    EXPECT_EQ(detail::threshold_graph_limit(100, ONE_GIB * 2), ONE_GIB * 2);
    // A budget never computes the default, so a count whose matrix size
    // cannot be represented is still served.
    EXPECT_EQ(detail::threshold_graph_limit(std::size_t{1} << 33, 7), size_t{7});
}

TEST(ThresholdGraphSizeTest, SizesThatDoNotFitThrowLengthError) {
    const size_t max = std::numeric_limits<size_t>::max();
    EXPECT_THROW(detail::threshold_graph_bytes(max / 2 + 1, 0), std::length_error);
    EXPECT_THROW(detail::threshold_graph_bytes(0, max / 4), std::length_error);
    // The pair count itself overflows.
    EXPECT_THROW(detail::default_threshold_graph_limit(std::size_t{1} << 33),
                 std::length_error);
    // The pair count fits; its size in bytes does not.
    EXPECT_THROW(detail::default_threshold_graph_limit(std::size_t{1} << 32),
                 std::length_error);
    EXPECT_THROW(detail::threshold_graph_limit(std::size_t{1} << 33, 0),
                 std::length_error);
}

// Pure-function tests alone would pass if the builder treated a zero budget
// as unlimited. Every pair of a constant-zero comparison is an edge, so the
// graph's size depends only on N, and a refusal comes after the degree pass,
// before any graph memory exists.
TEST(StreamingThresholdGraphTest, TheDefaultLimitFloorAdmitsASmallDenseGraph) {
    const size_t n = 100;
    ConstantComparison constant(n, 0.0);
    // Larger than the 39,600-byte matrix it replaces, so only the floor
    // admits it.
    ASSERT_GT(detail::threshold_graph_bytes(n, n * (n - 1) / 2),
              sizeof(double) * n * (n - 1) / 2);
    const ThresholdNeighborGraph graph =
        BuildThresholdNeighborGraph(constant, GraphOptions(0.0, 0));
    ASSERT_EQ(graph.Size(), n);
    for (size_t i = 0; i < n; ++i) {
        EXPECT_EQ(graph.Neighbors(i).size(), n);
    }
}

TEST(StreamingThresholdGraphTest, TheDefaultLimitRefusesAboveTheFloor) {
    const size_t n = 12000;
    ConstantComparison constant(n, 0.0);
    ASSERT_EQ(detail::default_threshold_graph_limit(n), ONE_GIB);
    ExpectMessage(
        [&] { BuildThresholdNeighborGraph(constant, GraphOptions(0.0, 0)); },
        "length_error",
        LimitMessage("threshold_graph", n, n * (n - 1) / 2, "default limit",
                     ONE_GIB));
}

TEST(StreamingThresholdGraphTest, TheDefaultLimitRefusesAboveTheMatrixSize) {
    const size_t n = 16500;
    ConstantComparison constant(n, 0.0);
    const size_t matrix = sizeof(double) * (n * (n - 1) / 2);
    ASSERT_GT(matrix, ONE_GIB);
    ASSERT_EQ(detail::default_threshold_graph_limit(n), matrix);
    ExpectMessage(
        [&] { BuildThresholdNeighborGraph(constant, GraphOptions(0.0, 0)); },
        "length_error",
        LimitMessage("threshold_graph", n, n * (n - 1) / 2, "default limit",
                     matrix));
}

TEST(StreamingThresholdGraphTest, AnExplicitBudgetIsExact) {
    const size_t n = 5;
    // Adjacent items only: four edges.
    TableComparison table(n, Line(n));
    const size_t exact = detail::threshold_graph_bytes(n, 4);
    ThresholdGraphOptions options = GraphOptions(1.0);
    options.max_graph_bytes = exact;
    EXPECT_EQ(BuildThresholdNeighborGraph(table, options).Size(), n);

    options.max_graph_bytes = exact - 1;
    options.caller = "butina";
    ExpectMessage([&] { BuildThresholdNeighborGraph(table, options); },
                  "length_error",
                  LimitMessage("butina", n, 4, "max_graph_bytes limit",
                               exact - 1));
}

// The budget is a contract on every input, so even the trivial graphs below
// two items are refused when it is smaller than they are.
TEST(StreamingThresholdGraphTest, AnExplicitBudgetHoldsBelowTwoItems) {
    for (const size_t n : {size_t{0}, size_t{1}}) {
        SCOPED_TRACE(::testing::Message() << "n " << n);
        TableComparison table(n, {});
        const size_t exact = detail::threshold_graph_bytes(n, 0);
        ThresholdGraphOptions options = GraphOptions(1.0);
        options.max_graph_bytes = exact;
        EXPECT_EQ(BuildThresholdNeighborGraph(table, options).Size(), n);
        options.max_graph_bytes = exact - 1;
        ExpectMessage([&] { BuildThresholdNeighborGraph(table, options); },
                      "length_error",
                      LimitMessage("threshold_graph", n, 0,
                                   "max_graph_bytes limit", exact - 1));
    }
}

// Item 1 and item 3 sit far apart in the table; the script decides whether
// the fill pass sees them as neighbors. Every other pair is far, so rows 1
// and 3 are the only ones the mismatch touches, and row 1 is reported.
TEST(StreamingThresholdGraphTest, RefusesAPairTheFillPassGains) {
    for (const size_t threads : {size_t{1}, size_t{4}}) {
        for (const size_t chunk : {size_t{1}, size_t{4096}}) {
            ScriptedPairComparison scripted(5, std::vector<double>(10, 9.0), 1, 3,
                                            {9.0, 0.0});
            ExpectMessage(
                [&] {
                    BuildThresholdNeighborGraph(
                        scripted, GraphOptions(1.0, threads, chunk));
                },
                "logic_error",
                "threshold_graph: item 1 gained a neighbor between the two "
                "threshold graph passes (counted 1, written at least 2); "
                "Compare must return the same value for a pair on every call");
        }
    }
}

TEST(StreamingThresholdGraphTest, RefusesAPairTheFillPassLoses) {
    for (const size_t threads : {size_t{1}, size_t{4}}) {
        for (const size_t chunk : {size_t{1}, size_t{4096}}) {
            ScriptedPairComparison scripted(5, std::vector<double>(10, 9.0), 1, 3,
                                            {0.0, 9.0});
            ExpectMessage(
                [&] {
                    BuildThresholdNeighborGraph(
                        scripted, GraphOptions(1.0, threads, chunk));
                },
                "logic_error",
                "threshold_graph: item 1 lost a neighbor between the two "
                "threshold graph passes (counted 2, written 1); Compare must "
                "return the same value for a pair on every call");
        }
    }
}

TEST(StreamingThresholdGraphTest, RefusesANonFiniteDistanceOnEitherPass) {
    for (const double bad : {NaN, INF, -INF}) {
        // {bad} reaches the degree pass; {0.5, bad} passes it and reaches the
        // fill pass.
        for (const std::vector<double>& script :
             {std::vector<double>{bad}, std::vector<double>{0.5, bad}}) {
            SCOPED_TRACE(::testing::Message()
                         << "value " << bad << ", pass " << script.size());
            ScriptedPairComparison scripted(5, std::vector<double>(10, 9.0), 1, 3,
                                            script);
            ThresholdGraphOptions options = GraphOptions(1.0);
            options.caller = "dbscan";
            ExpectMessage([&] { BuildThresholdNeighborGraph(scripted, options); },
                          "runtime_error",
                          "dbscan read a non-finite distance between items 1 "
                          "and 3");
            EXPECT_EQ(scripted.Evaluations(), script.size());
        }
    }
}
