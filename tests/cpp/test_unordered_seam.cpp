/**
 * @file test_unordered_seam.cpp
 * @brief The entry points that read the unordered-pair seam.
 *
 * ROCS is the only comparison in the library that overrides
 * ``PairwiseComparison::CompareUnordered``, and the ROCS tests elsewhere
 * assert to 1e-6 rather than to the bit, because OEShape promises no
 * bit-reproducibility. The double below stands in for it: a directional
 * ``Compare`` over an asymmetric table, and an override that keeps the
 * smaller of the two orientations, exactly as ROCS keeps the better overlay.
 * Routing a reader through the seam is then visible as a different number
 * rather than as a different call, which is what these tests assert --
 * against the storage overload over each of the two matrices, to the bit.
 *
 * Every test here is non-vacuous only because the two matrices disagree, and
 * ``TheTwoMatricesDisagree`` is what establishes that for all of them.
 *
 * ``activity_landscape`` and ``logdet_diversity`` are not tested here because
 * each shares its one read with an entry point that is: ``modelability``'s
 * ``ComparisonSweep`` and ``vendi_score``'s ``ComparisonPairs`` respectively.
 */

#include <gtest/gtest.h>

#include <algorithm>
#include <cstddef>
#include <memory>
#include <string>
#include <vector>

#include "oecluster/GateFacts.h"
#include "oecluster/PairwiseComparison.h"
#include "oecluster/StorageBackend.h"
#include "oecluster/clustering/ClusterReport.h"
#include "oecluster/clustering/ClusterTypes.h"
#include "oecluster/clustering/DiversitySelection.h"
#include "oecluster/clustering/KNNGraph.h"
#include "oecluster/clustering/SARCoherence.h"
#include "oecluster/clustering/SetDiversity.h"
#include "oecluster/clustering/SphereExclusion.h"

using namespace OECluster;

namespace {

constexpr size_t N = 10;

/// An asymmetric table of distances in (0, 1). Both triangles are needed, so
/// this is the full n-by-n square; the diagonal is zero.
double Directed(size_t i, size_t j) {
    if (i == j) {
        return 0.0;
    }
    return static_cast<double>((i * 37 + j * 11 + i * j * 3) % 97 + 1) / 100.0;
}

/// What a caller reading Compare(min, max) sees: the upper triangle only.
double UpperTriangle(size_t i, size_t j) {
    return Directed(std::min(i, j), std::max(i, j));
}

/// What the seam reports: the better -- here the smaller -- of the two
/// orientations, as ROCS keeps the better overlay.
double Unordered(size_t i, size_t j) {
    return std::min(Directed(i, j), Directed(j, i));
}

DenseStorage MakeStorage(double (*value)(size_t, size_t)) {
    DenseStorage storage(N);
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = i + 1; j < N; ++j) {
            storage.Set(i, j, value(i, j));
        }
    }
    return storage;
}

/// Directional, and symmetrizing only through the seam. Named something other
/// than "rocs" so that none of the refusals in DiversityValidation.h fires.
class DirectionalComparison : public PairwiseComparison {
public:
    double Compare(size_t i, size_t j) override { return Directed(i, j); }

    GateFacts Facts() const override { return GateFacts(); }

    std::unique_ptr<PairwiseComparison> Clone() const override {
        return std::make_unique<DirectionalComparison>(*this);
    }

    size_t Size() const override { return N; }
    std::string ComparisonName() const override { return "directional"; }

protected:
    double CompareUnordered(size_t i, size_t j) override {
        return Unordered(i, j);
    }
};

ClusteringResult TwoClusters() {
    std::vector<ClusterLabel> labels(N, 0);
    for (size_t i = N / 2; i < N; ++i) {
        labels[i] = 1;
    }
    return ClusteringResult(labels, labels_to_clusters(labels));
}

std::vector<std::string> TwoClasses() {
    std::vector<std::string> classes(N, "a");
    for (size_t i = N / 2; i < N; ++i) {
        classes[i] = "b";
    }
    return classes;
}

TEST(UnorderedSeam, TheTwoMatricesDisagree) {
    // Without this the tests below would pass on a build that routed nothing.
    size_t differing = 0;
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = i + 1; j < N; ++j) {
            if (UpperTriangle(i, j) != Unordered(i, j)) {
                ++differing;
            }
        }
    }
    // 21 of the 45 pairs, counted on this table rather than designed into it.
    EXPECT_EQ(differing, 21u) << "the fixture's two matrices no longer differ "
                                 "on the pairs the tests below rely on";
}

TEST(UnorderedSeam, ClusterReportReadsTheSeam) {
    DirectionalComparison comparison;
    const ClusteringResult result = TwoClusters();
    const DenseStorage seam = MakeStorage(&Unordered);
    const DenseStorage direct = MakeStorage(&UpperTriangle);

    const double lazy = cluster_report(result, comparison).silhouette;
    EXPECT_EQ(lazy, cluster_report(result, seam).silhouette);
    EXPECT_NE(lazy, cluster_report(result, direct).silhouette);
}

TEST(UnorderedSeam, ModelabilityReadsTheSeam) {
    DirectionalComparison comparison;
    const std::vector<std::string> classes = TwoClasses();
    const DenseStorage seam = MakeStorage(&Unordered);
    const DenseStorage direct = MakeStorage(&UpperTriangle);

    const double lazy = modelability(comparison, classes).modi;
    EXPECT_EQ(lazy, modelability(seam, classes).modi);
    EXPECT_NE(lazy, modelability(direct, classes).modi);
}

TEST(UnorderedSeam, MaxMinSelectionReadsTheSeam) {
    DirectionalComparison comparison;
    MaxMinOptions options;
    options.count = 4;
    const DenseStorage seam = MakeStorage(&Unordered);
    const DenseStorage direct = MakeStorage(&UpperTriangle);

    const std::vector<size_t> lazy = maxmin_select(comparison, options).indices;
    EXPECT_EQ(lazy, maxmin_select(seam, options).indices);
    EXPECT_NE(lazy, maxmin_select(direct, options).indices);
}

TEST(UnorderedSeam, CirclesReadsTheSeamUnderBothMethods) {
    DirectionalComparison comparison;
    const DenseStorage seam = MakeStorage(&Unordered);
    const DenseStorage direct = MakeStorage(&UpperTriangle);
    // Both methods, because they read through two different call sites.
    for (CirclesMethod method : {CirclesMethod::MaxMin, CirclesMethod::Sequential}) {
        CirclesOptions options;
        options.method = method;
        const std::vector<size_t> lazy =
            circles(comparison, 0.5, options).members;
        EXPECT_EQ(lazy, circles(seam, 0.5, options).members)
            << static_cast<int>(method);
        EXPECT_NE(lazy, circles(direct, 0.5, options).members)
            << static_cast<int>(method);
    }
}

TEST(UnorderedSeam, KNNGraphReadsTheSeam) {
    DirectionalComparison comparison;
    KNNGraphOptions options;
    options.k = 3;
    const DenseStorage seam = MakeStorage(&Unordered);
    const DenseStorage direct = MakeStorage(&UpperTriangle);

    const std::vector<double> lazy = knn_graph(comparison, options).Distances();
    EXPECT_EQ(lazy, knn_graph(seam, options).Distances());
    EXPECT_NE(lazy, knn_graph(direct, options).Distances());
}

TEST(UnorderedSeam, SphereExclusionReadsTheSeamAwayFromTheNeighborOrder) {
    DirectionalComparison comparison;
    SphereExclusionOptions options;
    options.distance_threshold = 0.4;
    options.order = SphereOrder::Input;
    const DenseStorage seam = MakeStorage(&Unordered);
    const DenseStorage direct = MakeStorage(&UpperTriangle);

    const std::vector<ClusterLabel> lazy =
        sphere_exclusion(comparison, options).Labels();
    EXPECT_EQ(lazy, sphere_exclusion(seam, options).Labels());
    EXPECT_NE(lazy, sphere_exclusion(direct, options).Labels());
}

TEST(UnorderedSeam, VendiScoreReadsTheSeam) {
    DirectionalComparison comparison;
    VendiOptions options;
    const DenseStorage seam = MakeStorage(&Unordered);
    const DenseStorage direct = MakeStorage(&UpperTriangle);

    const double lazy = vendi_score(comparison, options).score;
    EXPECT_EQ(lazy, vendi_score(seam, options).score);
    EXPECT_NE(lazy, vendi_score(direct, options).score);
}

}  // namespace
