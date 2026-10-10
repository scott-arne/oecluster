/**
 * @file test_streaming_single_linkage.cpp
 * @brief Single-linkage agglomerative clustering on the spanning-tree kernel.
 */

#include <gtest/gtest.h>

#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include <oechem.h>

#include "oecluster/Error.h"
#include "oecluster/GateFacts.h"
#include "oecluster/StorageBackend.h"
#include "oecluster/clustering/Agglomerative.h"
#include "oecluster/comparisons/DescriptorComparison.h"
#include "oecluster/comparisons/FingerprintComparison.h"

#include "agglomerative_oracle.h"
#include "diversity_test_support.h"
#include "mst_test_support.h"
#include "streaming_test_support.h"

using namespace OECluster;
using namespace diversity_test;
using namespace mst_test;
using namespace streaming_test;

namespace {

const double NaN = std::numeric_limits<double>::quiet_NaN();

struct Fixture {
    size_t n;
    std::vector<double> condensed;
};

std::vector<Fixture> Fixtures() {
    return {{2, Scrambled(2)},
            {7, Scrambled(7)},
            {12, ScrambledSixths(12)},
            {25, Hashed(25)},
            {40, Quantized(40, 3, 5)},
            {90, Quantized(90, 8, 1000)}};
}

// Distinct pairwise distances: no two merges tie, so single linkage has one
// valid tree and the heap and the spanning tree must agree on all of it.
std::vector<double> Distinct(size_t n) {
    return Condensed(n, [n](size_t i, size_t j) {
        return 1.0 + static_cast<double>((i * 7919 + j * 104729) % 1000003) / 1000003.0 +
               static_cast<double>(i * n + j) * 1e-9;
    });
}


AgglomerativeOptions Single(size_t n_clusters, double threshold = -1.0,
                            size_t threads = 1) {
    AgglomerativeOptions options;
    options.linkage = AgglomerativeLinkageMethod::Single;
    options.n_clusters = n_clusters;
    options.distance_threshold = threshold;
    options.num_threads = threads;
    return options;
}


void ExpectSame(const AgglomerativeResult& a, const AgglomerativeResult& b) {
    EXPECT_EQ(a.Labels(), b.Labels());
    EXPECT_EQ(a.Members(), b.Members());
    EXPECT_EQ(a.ChildrenLeft(), b.ChildrenLeft());
    EXPECT_EQ(a.ChildrenRight(), b.ChildrenRight());
    EXPECT_EQ(a.Distances(), b.Distances());
    EXPECT_EQ(a.ClusterSizes(), b.ClusterSizes());
}


}  // namespace

TEST(StreamingSingleLinkageTest, MatchesTheHeapWhenNoMergesTie) {
    for (size_t n : {2, 3, 8, 30, 75}) {
        const DenseStorage storage = MakeStorage(n, Distinct(n));
        for (size_t n_clusters : {size_t{1}, std::min<size_t>(3, n), n}) {
            for (bool full : {true, false}) {
                AgglomerativeOptions options = Single(n_clusters);
                options.compute_full_tree = full;
                const AgglomerativeResult expected = agglomerative_oracle::heap_cluster(storage, options);
                for (size_t threads : {1, 4}) {
                    options.num_threads = threads;
                    ExpectSame(agglomerative_cluster(storage, options), expected);
                }
            }
        }
        const AgglomerativeOptions cut = Single(1, 1.4);
        ExpectSame(agglomerative_cluster(storage, cut), agglomerative_oracle::heap_cluster(storage, cut));
    }
}

TEST(StreamingSingleLinkageTest, KeepsTheHeapsHeightsAndThresholdCutsWhenMergesTie) {
    for (const Fixture& fixture : Fixtures()) {
        const DenseStorage storage = MakeStorage(fixture.n, fixture.condensed);
        const AgglomerativeResult heap = agglomerative_oracle::heap_cluster(storage, Single(1));
        const AgglomerativeResult tree = agglomerative_cluster(storage, Single(1));
        EXPECT_EQ(tree.Distances(), heap.Distances()) << "n=" << fixture.n;
        for (double threshold : {0.0, 0.2, 0.5, 2.0, 3.5}) {
            EXPECT_EQ(agglomerative_cluster(storage, Single(1, threshold)).Labels(),
                      agglomerative_oracle::heap_cluster(storage, Single(1, threshold)).Labels())
                << "n=" << fixture.n << " threshold=" << threshold;
        }
    }
}

// The behavior change at ties, pinned: d(0,1) = d(1,2) = d(2,3) = 1, others 2.
TEST(StreamingSingleLinkageTest, MergesAtATiedHeightComeInTheTreesOrder) {
    const std::vector<double> condensed{1.0, 2.0, 2.0, 1.0, 2.0, 1.0};
    const DenseStorage storage = MakeStorage(4, condensed);
    const AgglomerativeResult heap = agglomerative_oracle::heap_cluster(storage, Single(2));
    EXPECT_EQ(heap.Members(), Clusters({{0, 1}, {2, 3}}));
    const AgglomerativeResult tree = agglomerative_cluster(storage, Single(2));
    EXPECT_EQ(tree.ChildrenLeft(), std::vector<size_t>({0, 2, 3}));
    EXPECT_EQ(tree.ChildrenRight(), std::vector<size_t>({1, 4, 5}));
    EXPECT_EQ(tree.Distances(), std::vector<double>({1.0, 1.0, 1.0}));
    EXPECT_EQ(tree.Members(), Clusters({{0, 1, 2}, {3}}));
}

TEST(StreamingSingleLinkageTest, TheComparisonPathMatchesACompareFilledMatrix) {
    for (const Fixture& fixture : Fixtures()) {
        TableComparison comparison(fixture.n, fixture.condensed);
        const DenseStorage filled = CompareFilled(comparison);
        for (const AgglomerativeOptions& base :
             {Single(1), Single(std::min<size_t>(3, fixture.n)), Single(1, 0.5)}) {
            const AgglomerativeResult expected = agglomerative_cluster(filled, base);
            for (size_t threads : {1, 4, 8}) {
                AgglomerativeOptions options = base;
                options.num_threads = threads;
                ExpectSame(agglomerative_cluster(comparison, options), expected);
            }
        }
    }
}

TEST(StreamingSingleLinkageTest, FingerprintsMatchACompareFilledMatrix) {
    std::vector<OEChem::OEGraphMol> mols = FingerprintMolecules();
    FingerprintComparison comparison(Pointers(mols));
    const DenseStorage filled = CompareFilled(comparison);
    for (const AgglomerativeOptions& base : {Single(4), Single(1, 0.5)}) {
        const AgglomerativeResult expected = agglomerative_cluster(filled, base);
        for (size_t threads : {1, 4}) {
            AgglomerativeOptions options = base;
            options.num_threads = threads;
            ExpectSame(agglomerative_cluster(comparison, options), expected);
        }
    }
}

TEST(StreamingSingleLinkageTest, DescriptorsMatchACompareFilledMatrix) {
    std::vector<OEChem::OEGraphMol> mols = FingerprintMolecules();
    DescriptorComparison comparison(Pointers(mols));
    const DenseStorage filled = CompareFilled(comparison);
    for (const AgglomerativeOptions& base : {Single(4), Single(1, 1.0)}) {
        const AgglomerativeResult expected = agglomerative_cluster(filled, base);
        for (size_t threads : {1, 4}) {
            AgglomerativeOptions options = base;
            options.num_threads = threads;
            ExpectSame(agglomerative_cluster(comparison, options), expected);
        }
    }
}

TEST(StreamingSingleLinkageTest, SmallInputsMatchTheMatrixOverload) {
    for (size_t n : {0, 1}) {
        PairCountingComparison comparison(n, {});
        ExpectSame(agglomerative_cluster(comparison, Single(n == 0 ? 1 : 1, 0.5)),
                   agglomerative_cluster(DenseStorage(n), Single(1, 0.5)));
    }
}

TEST(StreamingSingleLinkageTest, ComparesEveryPairExactlyOnce) {
    const size_t n = 45;
    PairCountingComparison comparison(n, Quantized(n, 2, 9));
    agglomerative_cluster(comparison, Single(3, -1.0, 4));
    EXPECT_EQ(comparison.Counts(), std::vector<size_t>(n * (n - 1) / 2, 1));
}

TEST(StreamingSingleLinkageTest, RefusesANonFiniteDistanceOnBothForms) {
    std::vector<double> condensed = Line(6);
    condensed[7] = NaN;
    TableComparison comparison(6, condensed);
    EXPECT_THROW(agglomerative_cluster(comparison, Single(2)), std::runtime_error);
    EXPECT_THROW(agglomerative_cluster(MakeStorage(6, condensed), Single(2)),
                 std::runtime_error);
}

TEST(StreamingSingleLinkageTest, NegativeZeroMergesAtPositiveZero) {
    std::vector<double> condensed = Line(4);
    condensed[0] = -0.0;
    TableComparison comparison(4, condensed);
    for (const AgglomerativeResult& result :
         {agglomerative_cluster(MakeStorage(4, condensed), Single(1)),
          agglomerative_cluster(comparison, Single(1))}) {
        ASSERT_EQ(result.Distances().front(), 0.0);
        EXPECT_FALSE(std::signbit(result.Distances().front()));
    }
}

TEST(StreamingSingleLinkageTest, RefusesOtherLinkagesRocsAndBadFactsBeforeReading) {
    for (auto linkage : {AgglomerativeLinkageMethod::Complete,
                         AgglomerativeLinkageMethod::Average,
                         AgglomerativeLinkageMethod::Weighted}) {
        AgglomerativeOptions options = Single(2);
        options.linkage = linkage;
        NamedRocsComparison rocs(6, Line(6));
        try {
            agglomerative_cluster(rocs, options);
            FAIL() << "accepted a non-single linkage";
        } catch (const std::invalid_argument& error) {
            EXPECT_NE(std::string(error.what()).find("only single linkage"),
                      std::string::npos)
                << error.what();
        }
        EXPECT_EQ(rocs.Total(), 0u);
    }

    NamedRocsComparison rocs(6, Line(6));
    try {
        agglomerative_cluster(rocs, Single(2));
        FAIL() << "accepted ROCS";
    } catch (const ComparisonError& error) {
        EXPECT_NE(std::string(error.what()).find(
                      "agglomerative cannot cluster from a ROCS comparison"),
                  std::string::npos)
            << error.what();
    }

    GateFacts similarity;
    similarity.is_distance = Capability::No;
    PairCountingComparison scores(6, Line(6), similarity);
    EXPECT_THROW(agglomerative_cluster(scores, Single(50)), ComparisonError);

    PairCountingComparison small(6, Line(6));
    try {
        agglomerative_cluster(small, Single(50));
        FAIL() << "accepted n_clusters above the item count";
    } catch (const std::invalid_argument& error) {
        EXPECT_STREQ(error.what(), "Agglomerative n_clusters must be at most the item count");
    }
    // A threshold cut ignores n_clusters, so the bound does not apply.
    EXPECT_NO_THROW(agglomerative_cluster(small, Single(50, 2.5)));
    EXPECT_EQ(rocs.Total() + scores.Total(), 0u);
}

TEST(StreamingSingleLinkageTest, ArgumentChecksComeFirst) {
    AgglomerativeOptions options = Single(2);
    options.chunk_size = 0;
    options.linkage = AgglomerativeLinkageMethod::Average;
    NamedRocsComparison rocs(6, Line(6));
    try {
        agglomerative_cluster(rocs, options);
        FAIL() << "accepted chunk_size 0";
    } catch (const std::invalid_argument& error) {
        EXPECT_STREQ(error.what(), "Agglomerative chunk_size must be at least one");
    }
}
