#include <gtest/gtest.h>

#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

#include "oecluster/StorageBackend.h"
#include "oecluster/clustering/ClusterReport.h"
#include "oecluster/clustering/ClusterTypes.h"
#include "oecluster/clustering/Representative.h"
#include "../../src/clustering/ClusterMetrics.h"
#include "../../src/clustering/ClusterReportTuned.h"
#include "diversity_test_support.h"
#include "report_test_support.h"

using namespace OECluster;
using diversity_test::Condensed;
using diversity_test::Hashed;
using diversity_test::MakeStorage;
using diversity_test::Quantized;
using report_test::SameDouble;
using report_test::SameReport;

namespace {

constexpr size_t UNBOUNDED = std::numeric_limits<size_t>::max();

ClusteringResult MakeResult(std::vector<ClusterLabel> labels) {
    Clusters members = labels_to_clusters(labels);
    return ClusteringResult(std::move(labels), std::move(members));
}

// 24 samples with members interleaved: cluster 0 has 12 members (66 pairs),
// cluster 1 has 5 (10 pairs), cluster 2 has 3 (3 pairs), cluster 3 is a
// singleton, and 6, 15 and 22 are noise. Budgets 0, 1 and 5 push clusters 0
// and 1 (and the 79-pair global median) onto the radix route while cluster 2
// stays direct at budget 5.
std::vector<ClusterLabel> MixedLabels() {
    return {0, 0, 1, 0, 2, 0, -1, 0, 0, 1, 0, 3, 0, 1, 0, -1, 0, 2, 0, 1, 0, 2, -1, 1};
}

// Hashed has no zeros; the shifted copy spans negative and positive values,
// which the radix key transform must order correctly.
std::vector<std::vector<double>> Fixtures(size_t n) {
    std::vector<double> shifted = Hashed(n);
    for (double& distance : shifted) {
        distance -= 0.5;
    }
    return {Hashed(n), shifted};
}

ClusterReport Tuned(const ClusteringResult& result, const StorageBackend& storage,
                    const ClusterReportOptions& options, size_t budget) {
    detail::ReportTuning tuning;
    tuning.median_direct_budget = budget;
    return detail::cluster_report_tuned(result, storage, options, tuning);
}

const RepresentativeMethod METHODS[] = {
    RepresentativeMethod::Medoid,
    RepresentativeMethod::Minimax,
    RepresentativeMethod::WeightedMedoid,
};

}  // namespace

TEST(ClusterReportEngineTest, SmallMedianBudgetsMatchTheDirectRoute) {
    const ClusteringResult result = MakeResult(MixedLabels());
    const size_t n = result.Labels().size();
    for (const std::vector<double>& condensed : Fixtures(n)) {
        const DenseStorage storage = MakeStorage(n, condensed);
        for (const RepresentativeMethod method : METHODS) {
            ClusterReportOptions options;
            options.representative_method = method;
            options.compute_per_cluster_records = true;
            const ClusterReport direct = Tuned(result, storage, options, UNBOUNDED);
            EXPECT_TRUE(SameReport(direct, cluster_report(result, storage, options)));
            for (const size_t budget : {size_t{0}, size_t{1}, size_t{5}}) {
                EXPECT_TRUE(SameReport(direct, Tuned(result, storage, options, budget)))
                    << "method " << static_cast<int>(method) << " budget " << budget;
            }
        }
    }
}

TEST(ClusterReportEngineTest, InlineSelectionMatchesClusterRepresentative) {
    const ClusteringResult result = MakeResult(MixedLabels());
    const size_t n = result.Labels().size();
    std::vector<std::vector<double>> fixtures = Fixtures(n);
    // Four distance levels give exact ties in both the means and the maxima.
    for (const unsigned seed : {1u, 2u, 3u}) {
        fixtures.push_back(Quantized(n, seed, 3));
    }
    for (const std::vector<double>& condensed : fixtures) {
        const DenseStorage storage = MakeStorage(n, condensed);
        for (const RepresentativeMethod method : METHODS) {
            ClusterReportOptions options;
            options.representative_method = method;
            options.compute_per_cluster_records = true;
            const ClusterReport report = cluster_report(result, storage, options);
            ASSERT_EQ(report.records.size(), result.Members().size());
            for (size_t k = 0; k < result.Members().size(); ++k) {
                EXPECT_EQ(report.records[k].representative,
                          cluster_representative(result.Members()[k], storage, method))
                    << "method " << static_cast<int>(method) << " cluster " << k;
            }
        }
    }
}

TEST(ClusterReportEngineTest, MedoidIndicesIgnoreTheConfiguredMethod) {
    const ClusteringResult result = MakeResult(MixedLabels());
    const size_t n = result.Labels().size();
    for (const std::vector<double>& condensed : Fixtures(n)) {
        const DenseStorage storage = MakeStorage(n, condensed);
        ClusterReportOptions medoid_options;
        ClusterReportOptions minimax_options;
        minimax_options.representative_method = RepresentativeMethod::Minimax;
        const ClusterReport medoid = cluster_report(result, storage, medoid_options);
        const ClusterReport minimax = cluster_report(result, storage, minimax_options);
        EXPECT_TRUE(SameDouble(medoid.calinski_harabasz_medoid, minimax.calinski_harabasz_medoid));
        EXPECT_TRUE(SameDouble(medoid.davies_bouldin_medoid, minimax.davies_bouldin_medoid));
        EXPECT_TRUE(SameDouble(medoid.dunn_medoid_separation_medoid_spread,
                               minimax.dunn_medoid_separation_medoid_spread));
    }
}

TEST(ClusterReportEngineTest, MediansEqualMedianDistanceOfTheIntraPairs) {
    const ClusteringResult result = MakeResult(MixedLabels());
    const size_t n = result.Labels().size();
    const Clusters& members = result.Members();
    for (const std::vector<double>& condensed : Fixtures(n)) {
        const DenseStorage storage = MakeStorage(n, condensed);
        ClusterReportOptions options;
        options.compute_per_cluster_records = true;
        for (const size_t budget : {size_t{0}, UNBOUNDED}) {
            const ClusterReport report = Tuned(result, storage, options, budget);
            std::vector<double> all;
            for (size_t k = 0; k < members.size(); ++k) {
                std::vector<double> mine;
                for (size_t i = 0; i < members[k].size(); ++i) {
                    for (size_t j = i + 1; j < members[k].size(); ++j) {
                        mine.push_back(storage.Get(members[k][i], members[k][j]));
                    }
                }
                all.insert(all.end(), mine.begin(), mine.end());
                if (mine.empty()) {
                    EXPECT_TRUE(std::isnan(report.records[k].median_intra_distance));
                } else {
                    EXPECT_EQ(report.records[k].median_intra_distance,
                              detail::median_distance(mine))
                        << "cluster " << k << " budget " << budget;
                }
            }
            EXPECT_EQ(report.median_intra_distance, detail::median_distance(all))
                << "budget " << budget;
        }
    }
}

TEST(ClusterReportEngineTest, ZeroMedianIsPositiveZero) {
    const ClusteringResult result = MakeResult({0, 0, 0, 1, 1});
    const DenseStorage storage = MakeStorage(5, Condensed(5, [](size_t i, size_t j) {
        return (i < 3) == (j < 3) ? -0.0 : 1.0;
    }));
    ClusterReportOptions options;
    options.compute_per_cluster_records = true;
    for (const size_t budget : {size_t{0}, UNBOUNDED}) {
        const ClusterReport report = Tuned(result, storage, options, budget);
        EXPECT_EQ(report.median_intra_distance, 0.0);
        EXPECT_FALSE(std::signbit(report.median_intra_distance)) << "budget " << budget;
        EXPECT_FALSE(std::signbit(report.records[0].median_intra_distance));
        EXPECT_FALSE(std::signbit(report.records[1].median_intra_distance));
    }
}

TEST(ClusterReportEngineTest, PairRankUsesEveryPairAtAnyBudget) {
    const ClusteringResult result = MakeResult(MixedLabels());
    const size_t n = result.Labels().size();
    const DenseStorage storage = MakeStorage(n, Hashed(n));
    ClusterReportOptions options;
    options.compute_pair_rank_indices = true;
    options.compute_per_cluster_records = true;
    const ClusterReport direct = Tuned(result, storage, options, UNBOUNDED);
    const ClusterReport tight = Tuned(result, storage, options, 0);
    EXPECT_TRUE(SameReport(direct, tight));
    EXPECT_TRUE(std::isfinite(tight.c_index));
}
