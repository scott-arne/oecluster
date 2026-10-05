/**
 * @file test_streaming_clustering.cpp
 * @brief Butina, DBSCAN and neighbor-order sphere exclusion over a comparison.
 */

#include <gtest/gtest.h>

#include <cstddef>
#include <functional>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include <oechem.h>
#include <oeomega2.h>

#include "oecluster/Error.h"
#include "oecluster/GateFacts.h"
#include "oecluster/StorageBackend.h"
#include "oecluster/clustering/Butina.h"
#include "oecluster/clustering/DBSCAN.h"
#include "oecluster/clustering/SphereExclusion.h"
#include "oecluster/comparisons/ROCSComparison.h"

#include "diversity_test_support.h"
#include "streaming_test_support.h"

using namespace OECluster;
using namespace diversity_test;
using namespace streaming_test;

namespace {

const double NaN = std::numeric_limits<double>::quiet_NaN();
const double INF = std::numeric_limits<double>::infinity();

struct Fixture {
    size_t n;
    std::vector<double> condensed;
};

// Tie-heavy tables, so neighbor-count ties and Butina's larger-index rule
// are exercised, plus seeded random levels.
std::vector<Fixture> Fixtures() {
    return {{2, Scrambled(2)},
            {7, Scrambled(7)},
            {12, ScrambledSixths(12)},
            {25, Hashed(25)},
            {40, Quantized(40, 3, 5)}};
}

ButinaOptions Butina(double threshold, bool reordering = false) {
    ButinaOptions options;
    options.distance_threshold = threshold;
    options.reordering = reordering;
    return options;
}

DBSCANOptions Dbscan(double eps, size_t min_samples) {
    DBSCANOptions options;
    options.eps = eps;
    options.min_samples = min_samples;
    return options;
}

SphereExclusionOptions Neighbors(double threshold, bool reordering,
                                 SphereAssignment assignment) {
    SphereExclusionOptions options;
    options.distance_threshold = threshold;
    options.order = SphereOrder::Neighbors;
    options.reordering = reordering;
    options.assignment = assignment;
    return options;
}

template <typename Exception>
void ExpectThrowWithMessage(const std::function<void()>& call,
                            const std::string& message) {
    try {
        call();
        FAIL() << "expected an exception: " << message;
    } catch (const Exception& error) {
        EXPECT_EQ(std::string(error.what()), message);
    }
}

}  // namespace

TEST(StreamingButinaTest, MatchesTheMatrixAtEveryThreadCountAndChunkSize) {
    for (const Fixture& fixture : Fixtures()) {
        const DenseStorage storage = MakeStorage(fixture.n, fixture.condensed);
        for (const double threshold : {0.0, 0.2, 0.5, 2.5, 4.0}) {
            for (const bool reordering : {false, true}) {
                const ButinaResult expected =
                    butina_cluster(storage, Butina(threshold, reordering));
                for (const size_t threads : {size_t{1}, size_t{2}, size_t{8}}) {
                    for (const size_t chunk : {size_t{1}, size_t{7}, size_t{4096}}) {
                        SCOPED_TRACE(::testing::Message()
                                     << "n " << fixture.n << ", threshold "
                                     << threshold << ", reordering "
                                     << reordering << ", threads " << threads
                                     << ", chunk " << chunk);
                        ButinaOptions options = Butina(threshold, reordering);
                        options.num_threads = threads;
                        options.chunk_size = chunk;
                        TableComparison table(fixture.n, fixture.condensed);
                        const ButinaResult result = butina_cluster(table, options);
                        EXPECT_EQ(result.Labels(), expected.Labels());
                        EXPECT_EQ(result.Members(), expected.Members());
                    }
                }
            }
        }
    }
}

TEST(StreamingDBSCANTest, MatchesTheMatrixAtEveryThreadCountAndChunkSize) {
    for (const Fixture& fixture : Fixtures()) {
        const DenseStorage storage = MakeStorage(fixture.n, fixture.condensed);
        // n + 1 makes every item noise: no row can hold more than n entries.
        for (const size_t min_samples :
             {size_t{1}, size_t{2}, size_t{3}, size_t{5}, fixture.n + 1}) {
            for (const double eps : {0.0, 0.2, 0.5, 2.5}) {
                const DBSCANResult expected =
                    dbscan_cluster(storage, Dbscan(eps, min_samples));
                for (const size_t threads : {size_t{1}, size_t{2}, size_t{8}}) {
                    for (const size_t chunk : {size_t{1}, size_t{7}, size_t{4096}}) {
                        SCOPED_TRACE(::testing::Message()
                                     << "n " << fixture.n << ", eps " << eps
                                     << ", min_samples " << min_samples
                                     << ", threads " << threads << ", chunk "
                                     << chunk);
                        DBSCANOptions options = Dbscan(eps, min_samples);
                        options.num_threads = threads;
                        options.chunk_size = chunk;
                        TableComparison table(fixture.n, fixture.condensed);
                        const DBSCANResult result = dbscan_cluster(table, options);
                        EXPECT_EQ(result.Labels(), expected.Labels());
                        EXPECT_EQ(result.Members(), expected.Members());
                        EXPECT_EQ(result.CoreSampleIndices(),
                                  expected.CoreSampleIndices());
                    }
                }
            }
        }
    }
}

TEST(StreamingDBSCANTest, AllNoiseWhenNoItemIsCore) {
    const size_t n = 6;
    TableComparison table(n, Line(n));
    const DBSCANResult result = dbscan_cluster(table, Dbscan(0.5, 2));
    EXPECT_EQ(result.Labels(), std::vector<ClusterLabel>(n, NOISE_LABEL));
    EXPECT_TRUE(result.CoreSampleIndices().empty());
}

TEST(StreamingSphereExclusionTest, NeighborsOrderMatchesTheMatrix) {
    for (const Fixture& fixture : Fixtures()) {
        const DenseStorage storage = MakeStorage(fixture.n, fixture.condensed);
        for (const double threshold : {0.0, 0.2, 0.5, 2.5}) {
            for (const bool reordering : {false, true}) {
                for (const SphereAssignment assignment :
                     {SphereAssignment::First, SphereAssignment::Nearest}) {
                    const SphereExclusionResult expected = sphere_exclusion(
                        storage, Neighbors(threshold, reordering, assignment));
                    for (const size_t threads : {size_t{1}, size_t{2}, size_t{8}}) {
                        for (const size_t chunk : {size_t{1}, size_t{7}, size_t{4096}}) {
                            SCOPED_TRACE(::testing::Message()
                                         << "n " << fixture.n << ", threshold "
                                         << threshold << ", reordering "
                                         << reordering << ", assignment "
                                         << static_cast<int>(assignment)
                                         << ", threads " << threads
                                         << ", chunk " << chunk);
                            SphereExclusionOptions options =
                                Neighbors(threshold, reordering, assignment);
                            options.num_threads = threads;
                            options.chunk_size = chunk;
                            TableComparison table(fixture.n, fixture.condensed);
                            const SphereExclusionResult result =
                                sphere_exclusion(table, options);
                            EXPECT_EQ(result.Labels(), expected.Labels());
                            EXPECT_EQ(result.Members(), expected.Members());
                            EXPECT_EQ(result.Centers(), expected.Centers());
                        }
                    }
                }
            }
        }
    }
}

TEST(StreamingClusteringTest, ZeroAndOneItems) {
    TableComparison empty(0, {});
    TableComparison single(1, {});
    EXPECT_EQ(butina_cluster(empty, Butina(1.0)).NumSamples(), 0u);
    EXPECT_EQ(butina_cluster(single, Butina(1.0)).Members(), Clusters{Cluster{0}});
    EXPECT_EQ(dbscan_cluster(empty, Dbscan(1.0, 1)).NumSamples(), 0u);
    EXPECT_EQ(dbscan_cluster(single, Dbscan(1.0, 1)).Members(), Clusters{Cluster{0}});
    EXPECT_EQ(dbscan_cluster(single, Dbscan(1.0, 2)).Labels(),
              std::vector<ClusterLabel>{NOISE_LABEL});
    const SphereExclusionOptions neighbors =
        Neighbors(1.0, false, SphereAssignment::First);
    EXPECT_EQ(sphere_exclusion(empty, neighbors).NumSamples(), 0u);
    EXPECT_EQ(sphere_exclusion(single, neighbors).Centers(),
              std::vector<size_t>{0});
}

TEST(StreamingClusteringTest, ZeroChunkSizeUsesTheDefault) {
    const size_t n = 25;
    const std::vector<double> condensed = Hashed(n);
    ButinaOptions butina = Butina(0.5);
    butina.num_threads = 4;
    DBSCANOptions dbscan = Dbscan(0.5, 3);
    dbscan.num_threads = 4;
    TableComparison table(n, condensed);
    const ButinaResult butina_default = butina_cluster(table, butina);
    const DBSCANResult dbscan_default = dbscan_cluster(table, dbscan);
    butina.chunk_size = 0;
    dbscan.chunk_size = 0;
    EXPECT_EQ(butina_cluster(table, butina).Members(), butina_default.Members());
    EXPECT_EQ(dbscan_cluster(table, dbscan).Members(), dbscan_default.Members());
    // Not vacuous: a self-only graph would leave every item a singleton.
    EXPECT_LT(butina_default.NumClusters(), n);

    SphereExclusionOptions sphere = Neighbors(0.5, false, SphereAssignment::First);
    sphere.chunk_size = 0;
    ExpectThrowWithMessage<std::invalid_argument>(
        [&] { sphere_exclusion(table, sphere); },
        "sphere_exclusion chunk_size must be at least one");
}

TEST(StreamingClusteringTest, RefusesComparisonsItsFactsRuleOutBeforeScoring) {
    struct Case {
        GateFacts facts;
        std::string reason;
    };
    GateFacts similarity;
    similarity.is_distance = Capability::No;
    GateFacts nonzero_self;
    nonzero_self.zero_self = Capability::No;
    GateFacts nan_present;
    nan_present.data_integrity = DataIntegrity::NaNPresent;
    const std::vector<Case> cases{
        {similarity, " requires distances, but the comparison reports similarities"},
        {nonzero_self,
         " requires a zero self-distance, but the comparison reports that "
         "d(x, x) is not zero"},
        {nan_present,
         " cannot rank distances the comparison declares may be non-finite "
         "(missing='propagate')"},
    };
    for (const Case& c : cases) {
        SCOPED_TRACE(c.reason);
        CountingComparison butina_counter(6, c.facts);
        ExpectThrowWithMessage<ComparisonError>(
            [&] { butina_cluster(butina_counter, Butina(1.0)); },
            "butina" + c.reason);
        EXPECT_EQ(butina_counter.Count(), 0u);
        CountingComparison dbscan_counter(6, c.facts);
        ExpectThrowWithMessage<ComparisonError>(
            [&] { dbscan_cluster(dbscan_counter, Dbscan(1.0, 2)); },
            "dbscan" + c.reason);
        EXPECT_EQ(dbscan_counter.Count(), 0u);
    }
}

// The subset and triangle tier is Python's, where allow_nonmetric can lift
// it; the C++ overloads leave it there, as the matrix overloads do.
TEST(StreamingClusteringTest, LeavesTheOverridableTierToPython) {
    GateFacts subset;
    subset.data_integrity = DataIntegrity::SubsetScored;
    subset.triangle = Capability::No;
    CountingComparison counter(4, subset);
    EXPECT_EQ(butina_cluster(counter, Butina(1.0)).NumSamples(), 4u);
    EXPECT_EQ(dbscan_cluster(counter, Dbscan(1.0, 2)).NumSamples(), 4u);
    EXPECT_GT(counter.Count(), 0u);
}

TEST(StreamingClusteringTest, ANonFiniteDistanceNamesTheEntryPoint) {
    const auto scripted = [](const std::vector<double>& script) {
        return ScriptedPairComparison(5, std::vector<double>(10, 9.0), 1, 3,
                                      script);
    };
    for (const double bad : {NaN, INF, -INF}) {
        for (const std::vector<double>& script :
             {std::vector<double>{bad}, std::vector<double>{0.5, bad}}) {
            SCOPED_TRACE(::testing::Message()
                         << "value " << bad << ", pass " << script.size());
            ScriptedPairComparison for_butina = scripted(script);
            ExpectThrowWithMessage<std::runtime_error>(
                [&] { butina_cluster(for_butina, Butina(1.0)); },
                "butina read a non-finite distance between items 1 and 3");
            ScriptedPairComparison for_dbscan = scripted(script);
            ExpectThrowWithMessage<std::runtime_error>(
                [&] { dbscan_cluster(for_dbscan, Dbscan(1.0, 2)); },
                "dbscan read a non-finite distance between items 1 and 3");
            ScriptedPairComparison for_sphere = scripted(script);
            ExpectThrowWithMessage<std::runtime_error>(
                [&] {
                    sphere_exclusion(
                        for_sphere, Neighbors(1.0, false, SphereAssignment::First));
                },
                "sphere_exclusion read a non-finite distance between items 1 "
                "and 3");
        }
    }
}

TEST(StreamingClusteringTest, TheBudgetReachesTheGraph) {
    const size_t n = 5;
    // Adjacent items only: four edges.
    TableComparison table(n, Line(n));
    const size_t exact = detail::threshold_graph_bytes(n, 4);
    const auto message = [&](const std::string& caller) {
        return caller + " would build a threshold graph of " +
               std::to_string(exact) + " bytes for 5 items and 4 edges, above "
               "its max_graph_bytes limit of " + std::to_string(exact - 1) +
               " bytes; use a tighter threshold, a larger max_graph_bytes, or "
               "a memory-mapped matrix from pdist(output=...)";
    };

    ButinaOptions butina = Butina(1.0);
    butina.max_graph_bytes = exact;
    EXPECT_EQ(butina_cluster(table, butina).NumClusters(), 2u);
    butina.max_graph_bytes = exact - 1;
    ExpectThrowWithMessage<std::length_error>(
        [&] { butina_cluster(table, butina); }, message("butina"));

    DBSCANOptions dbscan = Dbscan(1.0, 2);
    dbscan.max_graph_bytes = exact;
    EXPECT_EQ(dbscan_cluster(table, dbscan).NumClusters(), 1u);
    dbscan.max_graph_bytes = exact - 1;
    ExpectThrowWithMessage<std::length_error>(
        [&] { dbscan_cluster(table, dbscan); }, message("dbscan"));

    SphereExclusionOptions sphere = Neighbors(1.0, false, SphereAssignment::First);
    sphere.max_graph_bytes = exact;
    EXPECT_EQ(sphere_exclusion(table, sphere).NumClusters(), 2u);
    sphere.max_graph_bytes = exact - 1;
    ExpectThrowWithMessage<std::length_error>(
        [&] { sphere_exclusion(table, sphere); }, message("sphere_exclusion"));
}

TEST(StreamingClusteringTest, TheBudgetHoldsBelowTwoItems) {
    TableComparison single(1, {});
    const size_t exact = detail::threshold_graph_bytes(1, 0);
    const auto message = [&](const std::string& caller) {
        return caller + " would build a threshold graph of " +
               std::to_string(exact) + " bytes for 1 items and 0 edges, above "
               "its max_graph_bytes limit of " + std::to_string(exact - 1) +
               " bytes; use a tighter threshold, a larger max_graph_bytes, or "
               "a memory-mapped matrix from pdist(output=...)";
    };

    ButinaOptions butina = Butina(1.0);
    butina.max_graph_bytes = exact;
    EXPECT_EQ(butina_cluster(single, butina).Members(), Clusters{Cluster{0}});
    butina.max_graph_bytes = exact - 1;
    ExpectThrowWithMessage<std::length_error>(
        [&] { butina_cluster(single, butina); }, message("butina"));

    DBSCANOptions dbscan = Dbscan(1.0, 1);
    dbscan.max_graph_bytes = exact;
    EXPECT_EQ(dbscan_cluster(single, dbscan).Members(), Clusters{Cluster{0}});
    dbscan.max_graph_bytes = exact - 1;
    ExpectThrowWithMessage<std::length_error>(
        [&] { dbscan_cluster(single, dbscan); }, message("dbscan"));

    SphereExclusionOptions sphere = Neighbors(1.0, false, SphereAssignment::First);
    sphere.max_graph_bytes = exact;
    EXPECT_EQ(sphere_exclusion(single, sphere).Centers(), std::vector<size_t>{0});
    sphere.max_graph_bytes = exact - 1;
    ExpectThrowWithMessage<std::length_error>(
        [&] { sphere_exclusion(single, sphere); }, message("sphere_exclusion"));
}

TEST(StreamingClusteringTest, TheStorageOverloadsRefuseABudget) {
    const DenseStorage storage = MakeStorage(4, Line(4));
    ButinaOptions butina = Butina(1.0);
    butina.max_graph_bytes = 1;
    ExpectThrowWithMessage<std::invalid_argument>(
        [&] { butina_cluster(storage, butina); },
        "butina max_graph_bytes applies only when clustering from a comparison");
    DBSCANOptions dbscan = Dbscan(1.0, 2);
    dbscan.max_graph_bytes = 1;
    ExpectThrowWithMessage<std::invalid_argument>(
        [&] { dbscan_cluster(storage, dbscan); },
        "dbscan max_graph_bytes applies only when clustering from a comparison");
    SphereExclusionOptions sphere = Neighbors(1.0, false, SphereAssignment::First);
    sphere.max_graph_bytes = 1;
    ExpectThrowWithMessage<std::invalid_argument>(
        [&] { sphere_exclusion(storage, sphere); },
        "sphere_exclusion max_graph_bytes applies only when clustering from a "
        "comparison");
}

TEST(StreamingClusteringTest, SphereRefusesABudgetWithoutTheNeighborsOrder) {
    TableComparison table(4, Line(4));
    SphereExclusionOptions input;
    input.distance_threshold = 1.0;
    input.max_graph_bytes = 1 << 20;
    SphereExclusionOptions permutation = input;
    permutation.order = SphereOrder::Permutation;
    permutation.permutation = {3, 2, 1, 0};
    for (const SphereExclusionOptions& options : {input, permutation}) {
        ExpectThrowWithMessage<std::invalid_argument>(
            [&] { sphere_exclusion(table, options); },
            "sphere_exclusion max_graph_bytes requires the Neighbors order");
    }
}

// ROCS fails the repeatability precondition
// (test_comparison_repeatability.cpp, ROCSDependsOnItsCloneHistory), so every
// entry point that builds a threshold graph from a comparison refuses it by
// name before scoring a pair. Sphere exclusion's other orders build no graph
// and keep accepting it.
TEST(StreamingClusteringTest, RefusesROCSWhereverItWouldBuildAGraph) {
    OEConfGen::OEOmega omega;
    omega.SetMaxConfs(1);
    omega.SetStrictStereo(false);
    std::vector<std::shared_ptr<OEChem::OEMol>> mols;
    for (const char* smi : {"c1ccccc1", "Cc1ccccc1", "c1ccc(O)cc1"}) {
        auto mol = std::make_shared<OEChem::OEMol>();
        ASSERT_TRUE(OEChem::OESmilesToMol(*mol, smi)) << smi;
        ASSERT_TRUE(omega(*mol)) << smi;
        mols.push_back(mol);
    }
    ROCSComparison rocs(mols);
    const auto message = [](const std::string& caller) {
        return caller +
               " cannot build a threshold graph from a ROCS comparison: a ROCS "
               "score depends on what its overlay scored before, so the graph's "
               "two passes can disagree; cluster a matrix from pdist() instead";
    };
    ExpectThrowWithMessage<ComparisonError>(
        [&] { butina_cluster(rocs, Butina(0.5)); }, message("butina"));
    ExpectThrowWithMessage<ComparisonError>(
        [&] { dbscan_cluster(rocs, Dbscan(0.5, 2)); }, message("dbscan"));
    ExpectThrowWithMessage<ComparisonError>(
        [&] {
            sphere_exclusion(rocs, Neighbors(0.5, false, SphereAssignment::First));
        },
        message("sphere_exclusion"));

    SphereExclusionOptions input;
    input.distance_threshold = 0.5;
    EXPECT_EQ(sphere_exclusion(rocs, input).NumSamples(), 3u);
}
