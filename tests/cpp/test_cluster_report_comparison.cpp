#include <gtest/gtest.h>

#include <oechem.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "oecluster/Error.h"
#include "oecluster/GateFacts.h"
#include "oecluster/PairwiseComparison.h"
#include "oecluster/StorageBackend.h"
#include "oecluster/clustering/ClusterReport.h"
#include "oecluster/clustering/ClusterTypes.h"
#include "oecluster/comparisons/FingerprintComparison.h"
#include "../../src/clustering/ClusterReportTuned.h"
#include "diversity_test_support.h"
#include "report_test_support.h"

using namespace OECluster;
using diversity_test::CountingComparison;
using diversity_test::Hashed;
using diversity_test::MakeStorage;
using diversity_test::TableComparison;
using report_test::SameReport;

namespace {

constexpr size_t UNBOUNDED = std::numeric_limits<size_t>::max();

ClusteringResult MakeResult(std::vector<ClusterLabel> labels) {
    Clusters members = labels_to_clusters(labels);
    return ClusteringResult(std::move(labels), std::move(members));
}

// The engine test's interleaved partition: clusters of 12, 5, 3 and 1, and
// noise at 6, 15 and 22.
std::vector<ClusterLabel> MixedLabels() {
    return {0, 0, 1, 0, 2, 0, -1, 0, 0, 1, 0, 3, 0, 1, 0, -1, 0, 2, 0, 1, 0, 2, -1, 1};
}

std::vector<double> Shifted(size_t n) {
    std::vector<double> shifted = Hashed(n);
    for (double& distance : shifted) {
        distance -= 0.5;
    }
    return shifted;
}

size_t CondensedIndex(size_t n, size_t i, size_t j) {
    return n * i - i * (i + 1) / 2 + j - i - 1;
}

// Fails on chosen pairs, either by throwing or by returning NaN, and reads
// every other pair from a table that enforces Compare(min, max).
class FaultyComparison : public PairwiseComparison {
public:
    FaultyComparison(size_t n, std::vector<double> condensed,
                     std::vector<std::pair<size_t, size_t>> bad, bool throws)
        : table_(n, std::move(condensed)), bad_(std::move(bad)), throws_(throws) {}

    double Compare(size_t i, size_t j) override {
        for (const auto& pair : bad_) {
            if (pair.first == i && pair.second == j) {
                if (throws_) {
                    throw std::runtime_error(
                        "bad pair " + std::to_string(i) + " " + std::to_string(j));
                }
                return std::numeric_limits<double>::quiet_NaN();
            }
        }
        return table_.Compare(i, j);
    }

    GateFacts Facts() const override { return GateFacts(); }
    std::unique_ptr<PairwiseComparison> Clone() const override {
        return std::make_unique<FaultyComparison>(*this);
    }
    size_t Size() const override { return table_.Size(); }
    std::string ComparisonName() const override { return "faulty"; }

private:
    TableComparison table_;
    std::vector<std::pair<size_t, size_t>> bad_;
    bool throws_;
};

const RepresentativeMethod METHODS[] = {
    RepresentativeMethod::Medoid,
    RepresentativeMethod::Minimax,
    RepresentativeMethod::WeightedMedoid,
};

std::string MessageOf(const std::function<void()>& call) {
    try {
        call();
    } catch (const std::exception& error) {
        return error.what();
    }
    return "<no exception>";
}

}  // namespace

TEST(ClusterReportComparisonTest, MatchesTheMatrixForEveryThreadAndChunk) {
    const ClusteringResult result = MakeResult(MixedLabels());
    const size_t n = result.Labels().size();
    for (const std::vector<double>& condensed : {Hashed(n), Shifted(n)}) {
        const DenseStorage storage = MakeStorage(n, condensed);
        for (const RepresentativeMethod method : METHODS) {
            ClusterReportOptions options;
            options.representative_method = method;
            options.compute_per_cluster_records = true;
            const ClusterReport expected = cluster_report(result, storage, options);
            for (const size_t threads : {size_t{1}, size_t{2}, size_t{4}, size_t{0}}) {
                for (const size_t chunk : {size_t{1}, size_t{7}, size_t{4096}, UNBOUNDED}) {
                    TableComparison comparison(n, condensed);
                    options.num_threads = threads;
                    options.chunk_size = chunk;
                    EXPECT_TRUE(SameReport(expected, cluster_report(result, comparison, options)))
                        << "method " << static_cast<int>(method) << " threads " << threads
                        << " chunk " << chunk;
                }
            }
        }
    }
}

TEST(ClusterReportComparisonTest, SmallBudgetsAndBlocksMatchTheMatrix) {
    const ClusteringResult result = MakeResult(MixedLabels());
    const size_t n = result.Labels().size();
    for (const std::vector<double>& condensed : {Hashed(n), Shifted(n)}) {
        const DenseStorage storage = MakeStorage(n, condensed);
        ClusterReportOptions options;
        options.compute_per_cluster_records = true;
        options.num_threads = 4;
        options.chunk_size = 3;
        const ClusterReport expected = cluster_report(result, storage, options);

        std::mutex guard;
        size_t blocks = 0;
        size_t largest_block = 0;
        detail::ReportTuning tuning;
        tuning.median_direct_budget = 0;
        tuning.fill_block_distances = 16;
        tuning.on_block = [&](size_t distances, size_t) {
            std::lock_guard<std::mutex> lock(guard);
            ++blocks;
            largest_block = std::max(largest_block, distances);
        };
        TableComparison comparison(n, condensed);
        EXPECT_TRUE(SameReport(
            expected, detail::cluster_report_tuned(result, comparison, options, tuning)));
        EXPECT_GT(blocks, 0u);
        // No row in this fixture is longer than n, so a block exceeds the cap
        // only by finishing the row that crossed it.
        EXPECT_LE(largest_block, std::max<size_t>(16, n));
    }
}

TEST(ClusterReportComparisonTest, MatchesAMatrixFilledFromAFingerprintComparison) {
    const std::vector<const char*> smiles = {
        "c1ccccc1", "c1ccc(O)cc1", "c1ccc(N)cc1", "c1ccc(Cl)cc1",
        "CCCCCCCC", "CCCCCCO", "CCCCCCN", "CC(C)CCO",
        "C1CCCCC1", "C1CCCCC1O", "c1ccncc1", "CC(=O)Oc1ccccc1C(=O)O",
    };
    std::vector<OEChem::OEGraphMol> graph_mols(smiles.size());
    std::vector<OEChem::OEMolBase*> mols;
    for (size_t i = 0; i < smiles.size(); ++i) {
        ASSERT_TRUE(OEChem::OESmilesToMol(graph_mols[i], smiles[i]));
        mols.push_back(&static_cast<OEChem::OEMolBase&>(graph_mols[i]));
    }
    FingerprintComparison comparison(mols);
    const size_t n = mols.size();
    DenseStorage storage(n);
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i + 1; j < n; ++j) {
            storage.Set(i, j, comparison.Compare(i, j));
        }
    }
    const ClusteringResult result = MakeResult({0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 0, -1});
    for (const RepresentativeMethod method : METHODS) {
        ClusterReportOptions options;
        options.representative_method = method;
        options.compute_per_cluster_records = true;
        const ClusterReport expected = cluster_report(result, storage, options);
        for (const size_t threads : {size_t{1}, size_t{4}}) {
            options.num_threads = threads;
            EXPECT_TRUE(SameReport(expected, cluster_report(result, comparison, options)))
                << "method " << static_cast<int>(method) << " threads " << threads;
        }
    }
}

TEST(ClusterReportComparisonTest, EarliestCanonicalFailureWins) {
    const ClusteringResult result = MakeResult(MixedLabels());
    const size_t n = result.Labels().size();
    // (0, 2) is a cross pair and (2, 9) an intra pair of cluster 1, both
    // index-wise ahead of (16, 18); but cluster 0's intra rows are read first,
    // so (16, 18) is the earliest pair in canonical order.
    const std::vector<std::pair<size_t, size_t>> bad = {{0, 2}, {2, 9}, {16, 18}};

    std::vector<double> poisoned = Hashed(n);
    for (const auto& pair : bad) {
        poisoned[CondensedIndex(n, pair.first, pair.second)] =
            std::numeric_limits<double>::quiet_NaN();
    }
    const DenseStorage storage = MakeStorage(n, poisoned);
    ClusterReportOptions options;
    const std::string matrix_message =
        MessageOf([&] { cluster_report(result, storage, options); });
    EXPECT_EQ(matrix_message,
              "cluster_report: distance between samples 16 and 18 is not finite");

    for (const size_t threads : {size_t{1}, size_t{2}, size_t{4}, size_t{8}}) {
        for (const size_t chunk : {size_t{1}, size_t{7}}) {
            options.num_threads = threads;
            options.chunk_size = chunk;
            FaultyComparison throwing(n, Hashed(n), bad, true);
            EXPECT_EQ(MessageOf([&] { cluster_report(result, throwing, options); }),
                      "bad pair 16 18")
                << "threads " << threads << " chunk " << chunk;
            FaultyComparison nan(n, Hashed(n), bad, false);
            EXPECT_EQ(MessageOf([&] { cluster_report(result, nan, options); }), matrix_message)
                << "threads " << threads << " chunk " << chunk;
        }
    }
}

TEST(ClusterReportComparisonTest, ANonFiniteValueOnAnUnreadPairIsNeverSeen) {
    const ClusteringResult result = MakeResult(MixedLabels());
    const size_t n = result.Labels().size();
    // 6 and 15 are both noise: no metric reads a noise-to-noise pair.
    FaultyComparison comparison(n, Hashed(n), {{6, 15}}, false);
    ClusterReportOptions options;
    options.compute_per_cluster_records = true;
    const ClusterReport expected = cluster_report(result, MakeStorage(n, Hashed(n)), options);
    EXPECT_TRUE(SameReport(expected, cluster_report(result, comparison, options)));
}

TEST(ClusterReportComparisonTest, RefusesFactsBeforeAnyComparison) {
    const ClusteringResult result = MakeResult(MixedLabels());
    GateFacts similarity;
    similarity.is_distance = Capability::No;
    GateFacts nonzero_self;
    nonzero_self.zero_self = Capability::No;
    GateFacts nan_present;
    nan_present.data_integrity = DataIntegrity::NaNPresent;
    for (const GateFacts& facts : {similarity, nonzero_self, nan_present}) {
        CountingComparison comparison(24, facts);
        ClusterReportOptions options;
        options.compute_pair_rank_indices = true;  // facts outrank the pair-rank refusal
        EXPECT_THROW(cluster_report(result, comparison, options), ComparisonError);
        options.chunk_size = 0;  // and the chunk refusal outranks facts
        EXPECT_EQ(MessageOf([&] { cluster_report(result, comparison, options); }),
                  "cluster_report chunk_size must be at least one");
        EXPECT_EQ(comparison.Count(), 0u);
    }
}

TEST(ClusterReportComparisonTest, AcceptsTheMetricTierFacts) {
    // The C++ overload is gate-free beyond tier 1, as the matrix overload is;
    // the triangle and subset refusals live in Python behind allow_nonmetric.
    const ClusteringResult result = MakeResult(MixedLabels());
    GateFacts nonmetric;
    nonmetric.triangle = Capability::No;
    GateFacts subset;
    subset.data_integrity = DataIntegrity::SubsetScored;
    for (const GateFacts& facts : {nonmetric, subset}) {
        CountingComparison comparison(24, facts);
        EXPECT_NO_THROW(cluster_report(result, comparison));
        EXPECT_GT(comparison.Count(), 0u);
    }
}

TEST(ClusterReportComparisonTest, RefusesChunkThenPairRankBeforeAnyComparison) {
    const ClusteringResult result = MakeResult(MixedLabels());
    CountingComparison comparison(24, GateFacts());
    ClusterReportOptions options;
    options.chunk_size = 0;
    options.compute_pair_rank_indices = true;
    EXPECT_EQ(MessageOf([&] { cluster_report(result, comparison, options); }),
              "cluster_report chunk_size must be at least one");
    options.chunk_size = 4096;
    EXPECT_EQ(MessageOf([&] { cluster_report(result, comparison, options); }),
              "cluster_report cannot compute pair-rank indices from a comparison; pass a "
              "SymmetricDistanceMatrix or set compute_pair_rank_indices=False");
    EXPECT_THROW(cluster_report(result, comparison, options), std::invalid_argument);
    EXPECT_EQ(comparison.Count(), 0u);
}

TEST(ClusterReportComparisonTest, NamesTheComparisonInTheLabelCountRefusal) {
    const ClusteringResult result = MakeResult(MixedLabels());
    TableComparison comparison(20, Hashed(20));
    EXPECT_EQ(MessageOf([&] { cluster_report(result, comparison); }),
              "cluster_report: label count 24 exceeds the comparison sample count 20");
}
