/**
 * @file test_streaming_hdbscan.cpp
 * @brief HDBSCAN over a comparison and on the shared spanning-tree kernel.
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
#include "oecluster/clustering/HDBSCAN.h"
#include "oecluster/comparisons/DescriptorComparison.h"
#include "oecluster/comparisons/FingerprintComparison.h"

#include "../../src/clustering/HDBSCANLinkage.h"
#include "../../src/clustering/HDBSCANTree.h"
#include "diversity_test_support.h"
#include "mst_test_support.h"
#include "streaming_test_support.h"

using namespace OECluster;
using namespace diversity_test;
using namespace mst_test;
using namespace streaming_test;

namespace {

const double NaN = std::numeric_limits<double>::quiet_NaN();
const double INF = std::numeric_limits<double>::infinity();

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


HDBSCANOptions Hdbscan(size_t min_cluster_size, size_t min_samples, size_t threads = 1) {
    HDBSCANOptions options;
    options.min_cluster_size = min_cluster_size;
    options.min_samples = min_samples;
    options.num_threads = threads;
    return options;
}


void ExpectSame(const HDBSCANResult& a, const HDBSCANResult& b) {
    EXPECT_EQ(a.Labels(), b.Labels());
    EXPECT_EQ(a.Members(), b.Members());
    EXPECT_EQ(a.Probabilities(), b.Probabilities());
}


// HDBSCAN as 5.19.0 assembled it, on the legacy scan.
HDBSCANResult LegacyHdbscan(size_t n, const std::vector<double>& condensed,
                            const HDBSCANOptions& options) {
    const size_t min_samples =
        options.min_samples == 0 ? options.min_cluster_size : options.min_samples;
    const std::vector<double> core =
        detail::compute_core_distances(MakeStorage(n, condensed), min_samples, 1);
    const auto linkage = detail::make_hdbscan_single_linkage(
        LegacyPrim(n, condensed, core, options.alpha), n);
    const auto selection = detail::select_clusters(
        detail::condense_tree(linkage, options.min_cluster_size),
        options.cluster_selection_method, options.allow_single_cluster,
        options.cluster_selection_epsilon, options.max_cluster_size);
    std::vector<ClusterLabel> labels = selection.labels;
    std::vector<double> probabilities = selection.probabilities;
    if (labels.empty()) {
        labels.assign(n, NOISE_LABEL);
    }
    if (probabilities.empty()) {
        probabilities.assign(n, 0.0);
    }
    Clusters members = labels_to_clusters(labels);
    return HDBSCANResult(std::move(labels), std::move(members), std::move(probabilities));
}

}  // namespace

TEST(StreamingHDBSCANTest, TheMatrixPathMatchesTheLegacyPipeline) {
    for (const Fixture& fixture : Fixtures()) {
        if (fixture.n < 3) {
            continue;
        }
        const DenseStorage storage = MakeStorage(fixture.n, fixture.condensed);
        for (size_t min_samples : {1, 2, 3}) {
            for (auto method : {HDBSCANClusterSelectionMethod::EOM,
                                HDBSCANClusterSelectionMethod::Leaf}) {
                HDBSCANOptions options = Hdbscan(2, min_samples);
                options.cluster_selection_method = method;
                const HDBSCANResult expected =
                    LegacyHdbscan(fixture.n, fixture.condensed, options);
                for (size_t threads : {1, 4}) {
                    options.num_threads = threads;
                    ExpectSame(hdbscan_cluster(storage, options), expected);
                }
            }
        }
    }
}

TEST(StreamingHDBSCANTest, TheComparisonPathMatchesACompareFilledMatrix) {
    for (const Fixture& fixture : Fixtures()) {
        if (fixture.n < 3) {
            continue;
        }
        TableComparison comparison(fixture.n, fixture.condensed);
        const DenseStorage filled = CompareFilled(comparison);
        for (size_t min_samples : {1, 2, 3}) {
            const HDBSCANResult expected = hdbscan_cluster(filled, Hdbscan(2, min_samples));
            for (size_t threads : {1, 4, 8}) {
                ExpectSame(hdbscan_cluster(comparison, Hdbscan(2, min_samples, threads)),
                           expected);
            }
        }
    }
}

TEST(StreamingHDBSCANTest, FingerprintsMatchACompareFilledMatrix) {
    std::vector<OEChem::OEGraphMol> mols = FingerprintMolecules();
    FingerprintComparison comparison(Pointers(mols));
    const DenseStorage filled = CompareFilled(comparison);
    for (size_t min_samples : {1, 2, 4}) {
        const HDBSCANResult expected = hdbscan_cluster(filled, Hdbscan(2, min_samples));
        for (size_t threads : {1, 4}) {
            ExpectSame(hdbscan_cluster(comparison, Hdbscan(2, min_samples, threads)),
                       expected);
        }
    }
}

TEST(StreamingHDBSCANTest, DescriptorsMatchACompareFilledMatrix) {
    std::vector<OEChem::OEGraphMol> mols = FingerprintMolecules();
    DescriptorComparison comparison(Pointers(mols));
    const DenseStorage filled = CompareFilled(comparison);
    for (size_t min_samples : {1, 3}) {
        const HDBSCANResult expected = hdbscan_cluster(filled, Hdbscan(2, min_samples));
        for (size_t threads : {1, 4}) {
            ExpectSame(hdbscan_cluster(comparison, Hdbscan(2, min_samples, threads)),
                       expected);
        }
    }
}

TEST(StreamingHDBSCANTest, ComparesEveryPairOnceOrTwice) {
    const size_t n = 45;
    const std::vector<double> condensed = Quantized(n, 6, 9);
    // min_samples = 1 has no core pass, and its Prim pass is unpruned.
    PairCountingComparison single_pass(n, condensed);
    hdbscan_cluster(single_pass, Hdbscan(2, 1, 4));
    EXPECT_EQ(single_pass.Counts(), std::vector<size_t>(condensed.size(), 1));
    // Otherwise the core pass reads every pair once and the pruned Prim pass
    // at most once more.
    PairCountingComparison two_passes(n, condensed);
    hdbscan_cluster(two_passes, Hdbscan(2, 4, 4));
    for (const size_t count : two_passes.Counts()) {
        EXPECT_GE(count, 1u);
        EXPECT_LE(count, 2u);
    }
}

TEST(StreamingHDBSCANTest, RefusesValuesOutsideTheDomainOnBothForms) {
    for (double bad : {NaN, INF, -0.5}) {
        std::vector<double> condensed = Line(8);
        condensed[9] = bad;  // the pair (1, 4)
        for (size_t min_samples : {1, 3}) {
            TableComparison comparison(8, condensed);
            EXPECT_THROW(hdbscan_cluster(comparison, Hdbscan(2, min_samples)),
                         std::runtime_error);
            EXPECT_THROW(hdbscan_cluster(MakeStorage(8, condensed), Hdbscan(2, min_samples)),
                         std::runtime_error);
        }
    }
}

TEST(StreamingHDBSCANTest, AnOverflowingAlphaIsRefusedOnBothFormsAndRegimes) {
    for (size_t min_samples : {1, 3}) {
        HDBSCANOptions options = Hdbscan(2, min_samples);
        options.alpha = 1e-310;
        TableComparison comparison(8, Line(8));
        EXPECT_THROW(hdbscan_cluster(comparison, options), std::invalid_argument);
        EXPECT_THROW(hdbscan_cluster(MakeStorage(8, Line(8)), options),
                     std::invalid_argument);
    }
}

// Only the pre-check on the largest distance can refuse this input. Every core
// distance is zero, so at the second Prim step the one pair big enough to
// overflow, (1, 2), is pruned before it is read and the kernel's own quotient
// check never sees it.
TEST(StreamingHDBSCANTest, AnOverflowingAlphaIsRefusedOnAPairThePruningSkips) {
    const std::vector<double> condensed{0.0, 0.0, 1e9};  // (0, 1), (0, 2), (1, 2)
    HDBSCANOptions options = Hdbscan(2, 2);
    options.alpha = 1e-300;
    const char* expected =
        "hdbscan alpha=1e-300 is too small: a distance divided by alpha is not finite";
    TableComparison comparison(3, condensed);
    try {
        hdbscan_cluster(comparison, options);
        ADD_FAILURE() << "the comparison form accepted the alpha";
    } catch (const std::invalid_argument& error) {
        EXPECT_STREQ(error.what(), expected);
    }
    try {
        hdbscan_cluster(MakeStorage(3, condensed), options);
        ADD_FAILURE() << "the matrix form accepted the alpha";
    } catch (const std::invalid_argument& error) {
        EXPECT_STREQ(error.what(), expected);
    }
}

TEST(StreamingHDBSCANTest, ANanAlphaIsRefusedBeforeAnyPairIsRead) {
    for (size_t min_samples : {1, 3}) {
        HDBSCANOptions options = Hdbscan(2, min_samples);
        options.alpha = NaN;
        PairCountingComparison comparison(8, Line(8));
        try {
            hdbscan_cluster(comparison, options);
            FAIL() << "accepted a NaN alpha";
        } catch (const std::invalid_argument& error) {
            EXPECT_STREQ(error.what(), "HDBSCAN alpha must be positive");
        }
        EXPECT_EQ(comparison.Total(), 0u);
        EXPECT_THROW(hdbscan_cluster(MakeStorage(8, Line(8)), options),
                     std::invalid_argument);
    }
}

TEST(StreamingHDBSCANTest, RefusesRocsAndBadFactsBeforeAnyPairIsRead) {
    NamedRocsComparison rocs(6, Line(6));
    try {
        hdbscan_cluster(rocs, Hdbscan(2, 2));
        FAIL() << "accepted ROCS";
    } catch (const ComparisonError& error) {
        EXPECT_NE(std::string(error.what()).find(
                      "hdbscan cannot cluster from a ROCS comparison"),
                  std::string::npos)
            << error.what();
    }
    EXPECT_EQ(rocs.Total(), 0u);

    GateFacts similarity;
    similarity.is_distance = Capability::No;
    PairCountingComparison scores(6, Line(6), similarity);
    EXPECT_THROW(hdbscan_cluster(scores, Hdbscan(2, 2)), ComparisonError);
    EXPECT_EQ(scores.Total(), 0u);
}

TEST(StreamingHDBSCANTest, TheErrorsComeInTheirOrder) {
    // A bad alpha beats ROCS: local arguments come before the gate.
    HDBSCANOptions bad_alpha = Hdbscan(2, 2);
    bad_alpha.alpha = NaN;
    NamedRocsComparison rocs(6, Line(6));
    EXPECT_THROW(hdbscan_cluster(rocs, bad_alpha), std::invalid_argument);

    // ROCS beats bad facts, and both beat the min_samples bound.
    GateFacts similarity;
    similarity.is_distance = Capability::No;
    NamedRocsComparison rocs_scores(6, Line(6), similarity);
    try {
        hdbscan_cluster(rocs_scores, Hdbscan(2, 50));
        FAIL() << "accepted ROCS";
    } catch (const ComparisonError& error) {
        EXPECT_NE(std::string(error.what()).find("ROCS"), std::string::npos);
    }
    PairCountingComparison scores(6, Line(6), similarity);
    try {
        hdbscan_cluster(scores, Hdbscan(2, 50));
        FAIL() << "accepted a similarity";
    } catch (const ComparisonError& error) {
        EXPECT_NE(std::string(error.what()).find("similarities"), std::string::npos);
    }

    // The bound is checked before any pair is read.
    PairCountingComparison small(6, Line(6));
    try {
        hdbscan_cluster(small, Hdbscan(2, 50));
        FAIL() << "accepted min_samples above the item count";
    } catch (const std::invalid_argument& error) {
        EXPECT_STREQ(error.what(), "HDBSCAN min_samples must be at most the item count");
    }
    EXPECT_EQ(small.Total() + scores.Total() + rocs_scores.Total() + rocs.Total(), 0u);
}

