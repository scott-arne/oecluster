#include <cmath>
#include <initializer_list>
#include <stdexcept>
#include <vector>

#include <gtest/gtest.h>

#include "oefp/batch.h"
#include "oefp/fingerprint.h"
#include "oecluster/clustering/BitBirch.h"
#include "oecluster/clustering/ClusterReport.h"
#include "oecluster/StorageBackend.h"
#include "../../src/clustering/BitBirchTree.h"

namespace OECluster::detail {
size_t partition_count(size_t n);
}  // namespace OECluster::detail

namespace {

OEFP::OEFP make_fp(const size_t size_bits, std::initializer_list<size_t> on_bits) {
    OEFP::FingerprintSpec spec;
    spec.size_bits = size_bits;
    spec.value_type = OEFP::FingerprintValueType::Binary;
    spec.source_name = "test";
    OEFP::OEFP fp(spec);
    for (const size_t bit : on_bits) {
        fp.SetBit(bit);
    }
    return fp;
}

OEFP::OEFPBatch make_batch(std::initializer_list<OEFP::OEFP> fps) {
    return OEFP::OEFPBatch::FromFingerprints(std::vector<OEFP::OEFP>(fps));
}

OEFP::OEFPBatch make_random_batch(const size_t rows, const size_t bits) {
    std::vector<OEFP::OEFP> fps;
    fps.reserve(rows);
    uint64_t state = 88172645463325252ull;  // fixed seed, xorshift64
    auto next = [&state]() {
        state ^= state << 13; state ^= state >> 7; state ^= state << 17;
        return state;
    };
    for (size_t i = 0; i < rows; ++i) {
        OEFP::FingerprintSpec spec;
        spec.size_bits = bits;
        spec.value_type = OEFP::FingerprintValueType::Binary;
        spec.source_name = "test";
        OEFP::OEFP fp(spec);
        for (size_t b = 0; b < bits; ++b) {
            if ((next() & 7u) == 0u) fp.SetBit(b);  // ~12.5% density
        }
        fp.SetBit(i % bits);  // guarantee at least one bit
        fps.push_back(fp);
    }
    return OEFP::OEFPBatch::FromFingerprints(fps);
}

OECluster::DenseStorage tanimoto_storage(const OEFP::OEFPBatch& batch) {
    const size_t n = batch.Size();
    const size_t words = batch.WordsPerFingerprint();
    OECluster::DenseStorage storage(n);
    for (size_t i = 0; i < n; ++i) {
        const uint64_t* wi = batch.RowWords(i);
        const uint32_t pi = batch.PopCount(i);
        for (size_t j = i + 1; j < n; ++j) {
            const uint64_t* wj = batch.RowWords(j);
            uint32_t inter = 0;
            for (size_t w = 0; w < words; ++w) {
                inter += static_cast<uint32_t>(__builtin_popcountll(wi[w] & wj[w]));
            }
            const uint32_t uni = pi + batch.PopCount(j) - inter;
            const double sim = uni == 0u ? 1.0 : static_cast<double>(inter) /
                                                  static_cast<double>(uni);
            storage.Set(i, j, 1.0 - sim);
        }
    }
    return storage;
}

}  // namespace

TEST(BitBirchClusteringTest, ClustersDuplicateBlocksWithDiameterCriterion) {
    const auto batch = make_batch({
        make_fp(4, {0, 1}),
        make_fp(4, {0, 1}),
        make_fp(4, {2, 3}),
        make_fp(4, {2, 3}),
    });
    OECluster::BitBirchOptions options;
    options.threshold = 0.75;
    options.branching_factor = 2;
    options.merge_criterion = OECluster::BitBirchMergeCriterion::Diameter;

    const auto result = OECluster::bitbirch_cluster(batch, options);

    EXPECT_EQ(result.Labels(), (std::vector<OECluster::ClusterLabel>{0, 0, 1, 1}));
    ASSERT_EQ(result.Members().size(), 2u);
    EXPECT_EQ(result.Members()[0], (OECluster::Cluster{0, 1}));
    EXPECT_EQ(result.Members()[1], (OECluster::Cluster{2, 3}));
    EXPECT_EQ(result.ClusterSizes(), (std::vector<size_t>{2, 2}));
    ASSERT_EQ(result.Centroids().Size(), 2u);
    EXPECT_EQ(result.Centroids().PopCount(0), 2u);
    EXPECT_EQ(result.Centroids().PopCount(1), 2u);
}

TEST(BitBirchClusteringTest, PreservesReferenceSplitTieOrdering) {
    const auto batch = make_batch({
        make_fp(4, {0, 1}),
        make_fp(4, {0, 1}),
        make_fp(4, {2, 3}),
        make_fp(4, {2, 3}),
        make_fp(4, {0, 2}),
        make_fp(4, {0, 2}),
    });
    OECluster::BitBirchOptions options;
    options.threshold = 0.75;
    options.branching_factor = 2;
    options.merge_criterion = OECluster::BitBirchMergeCriterion::Diameter;

    const auto result = OECluster::bitbirch_cluster(batch, options);

    ASSERT_EQ(result.Members().size(), 4u);
    EXPECT_EQ(result.Members()[0], (OECluster::Cluster{0, 1}));
    EXPECT_EQ(result.Members()[1], (OECluster::Cluster{2, 3}));
    EXPECT_EQ(result.Members()[2], (OECluster::Cluster{5}));
    EXPECT_EQ(result.Members()[3], (OECluster::Cluster{4}));
    EXPECT_EQ(result.Labels(), (std::vector<OECluster::ClusterLabel>{0, 0, 1, 1, 3, 2}));
}

TEST(BitBirchClusteringTest, RejectsInvalidOptions) {
    const auto batch = make_batch({make_fp(4, {0})});
    OECluster::BitBirchOptions options;

    options.threshold = -0.1;
    EXPECT_THROW(OECluster::bitbirch_cluster(batch, options), std::invalid_argument);

    options.threshold = 0.65;
    options.branching_factor = 0;
    EXPECT_THROW(OECluster::bitbirch_cluster(batch, options), std::invalid_argument);

    options.branching_factor = 1;
    options.tolerance = -0.1;
    EXPECT_THROW(OECluster::bitbirch_cluster(batch, options), std::invalid_argument);
}

TEST(BitBirchClusteringTest, RefinePruneRequiresParentPointers) {
    const auto batch = make_batch({
        make_fp(4, {0, 1}),
        make_fp(4, {0, 1}),
        make_fp(4, {2, 3}),
    });
    OECluster::BitBirchRefinementOptions options;
    options.fit_options.singly = true;
    options.redistribute_largest_cluster = true;

    EXPECT_THROW(OECluster::bitbirch_refine(batch, options), std::invalid_argument);
}

TEST(BitBirchClusteringTest, RefinePruneRedistributesLargestCluster) {
    const auto batch = make_batch({
        make_fp(8, {2, 4}),
        make_fp(8, {1, 4, 6}),
        make_fp(8, {0, 2, 3, 5}),
        make_fp(8, {3, 4, 7}),
        make_fp(8, {4, 7}),
        make_fp(8, {3}),
        make_fp(8, {0, 4, 6, 7}),
        make_fp(8, {4, 5}),
        make_fp(8, {1, 2, 6}),
        make_fp(8, {3}),
        make_fp(8, {1, 3, 5, 7}),
        make_fp(8, {1, 3, 5, 6}),
    });
    OECluster::BitBirchRefinementOptions options;
    options.fit_options.threshold = 0.65;
    options.fit_options.branching_factor = 2;
    options.fit_options.merge_criterion = OECluster::BitBirchMergeCriterion::Diameter;
    options.fit_options.singly = false;
    options.redistribute_largest_cluster = true;

    const auto result = OECluster::bitbirch_refine(batch, options);

    EXPECT_EQ(result.Labels().size(), 12u);
    for (const OECluster::ClusterLabel label : result.Labels()) {
        EXPECT_GE(label, 0);
    }
    size_t assigned = 0;
    for (const auto& cluster : result.Members()) {
        assigned += cluster.size();
    }
    EXPECT_EQ(assigned, 12u);
}

TEST(BitBirchClusteringTest, RefinePruneMatchesZeroSampleReferenceCentroid) {
    const auto batch = make_batch({
        make_fp(3, {0, 1, 2}),
        make_fp(3, {1, 2}),
        make_fp(3, {0}),
        make_fp(3, {0, 1}),
    });
    OECluster::BitBirchRefinementOptions options;
    options.fit_options.threshold = 0.8;
    options.fit_options.branching_factor = 2;
    options.fit_options.merge_criterion = OECluster::BitBirchMergeCriterion::Diameter;
    options.fit_options.singly = false;
    options.redistribute_largest_cluster = true;

    const auto result = OECluster::bitbirch_refine(batch, options);

    EXPECT_EQ(result.Labels(), (std::vector<OECluster::ClusterLabel>{2, 1, 0, 3}));
    ASSERT_EQ(result.Members().size(), 5u);
    EXPECT_EQ(result.Members()[0], (OECluster::Cluster{2}));
    EXPECT_EQ(result.Members()[1], (OECluster::Cluster{1}));
    EXPECT_EQ(result.Members()[2], (OECluster::Cluster{0}));
    EXPECT_EQ(result.Members()[3], (OECluster::Cluster{3}));
    EXPECT_TRUE(result.Members()[4].empty());
}

TEST(BitBirchFastTest, PartitionCountIsDeterministicInNAlone) {
    using OECluster::detail::partition_count;
    EXPECT_EQ(partition_count(0u), 1u);
    EXPECT_EQ(partition_count(1u), 1u);
    EXPECT_EQ(partition_count(2048u), 1u);     // at target -> single partition
    EXPECT_EQ(partition_count(2049u), 2u);     // just over target -> two
    EXPECT_EQ(partition_count(4096u), 2u);
    EXPECT_EQ(partition_count(4097u), 3u);
    // Capped at BITBIRCH_FAST_MAX_CHUNKS (64).
    EXPECT_EQ(partition_count(64u * 2048u * 4u), 64u);
}

TEST(BitBirchFastTest, FitRangeOverFullRangeMatchesFitAll) {
    const auto batch = make_batch({
        make_fp(4, {0, 1}),
        make_fp(4, {0, 1}),
        make_fp(4, {2, 3}),
        make_fp(4, {2, 3}),
    });
    OECluster::BitBirchOptions options;
    options.threshold = 0.75;
    options.branching_factor = 2;
    options.merge_criterion = OECluster::BitBirchMergeCriterion::Diameter;

    OECluster::detail::BitBirchTree tree_all(options);
    tree_all.Fit(batch);
    const auto result_all = tree_all.Result(batch.Spec(), batch.Size());

    OECluster::detail::BitBirchTree tree_range(options);
    tree_range.Fit(batch, 0, batch.Size());
    const auto result_range = tree_range.Result(batch.Spec(), batch.Size());

    EXPECT_EQ(result_all.Labels(), result_range.Labels());
    EXPECT_EQ(result_all.Members(), result_range.Members());
}

TEST(BitBirchFastTest, BuildFastTreeSinglePartitionMatchesStrict) {
    const auto batch = make_batch({
        make_fp(4, {0, 1}),
        make_fp(4, {0, 1}),
        make_fp(4, {2, 3}),
        make_fp(4, {2, 3}),
        make_fp(4, {0, 2}),
        make_fp(4, {0, 2}),
    });
    OECluster::BitBirchOptions options;
    options.threshold = 0.75;
    options.branching_factor = 2;
    options.merge_criterion = OECluster::BitBirchMergeCriterion::Diameter;

    OECluster::detail::BitBirchTree strict_tree(options);
    strict_tree.Fit(batch);
    const auto strict = strict_tree.Result(batch.Spec(), batch.Size());

    OECluster::detail::BitBirchTree fast_tree(options);
    OECluster::detail::BitBirchTree::BuildFastTree(batch, options, fast_tree);
    const auto fast = fast_tree.Result(batch.Spec(), batch.Size());

    EXPECT_EQ(strict.Labels(), fast.Labels());
    EXPECT_EQ(strict.Members(), fast.Members());
    EXPECT_EQ(strict.ClusterSizes(), fast.ClusterSizes());
}

TEST(BitBirchFastTest, FastResultIsDeterministicAcrossThreadCounts) {
    const auto batch = make_random_batch(5000, 64);

    OECluster::BitBirchOptions base;
    base.threshold = 0.5;
    base.branching_factor = 50;
    base.merge_criterion = OECluster::BitBirchMergeCriterion::Diameter;

    auto run = [&](size_t threads) {
        OECluster::BitBirchOptions opts = base;
        opts.num_threads = threads;
        OECluster::detail::BitBirchTree tree(opts);
        OECluster::detail::BitBirchTree::BuildFastTree(batch, opts, tree);
        return tree.Result(batch.Spec(), batch.Size());
    };

    const auto r1 = run(1);
    for (const size_t threads : {size_t{2}, size_t{4}, size_t{0}}) {
        const auto r = run(threads);
        EXPECT_EQ(r1.Labels(), r.Labels()) << "threads=" << threads;
        EXPECT_EQ(r1.Members(), r.Members()) << "threads=" << threads;
    }
}

TEST(BitBirchFastTest, BitBirchClusterFastSmallNMatchesStrict) {
    const auto batch = make_batch({
        make_fp(4, {0, 1}),
        make_fp(4, {0, 1}),
        make_fp(4, {2, 3}),
        make_fp(4, {2, 3}),
    });
    OECluster::BitBirchOptions strict_opts;
    strict_opts.threshold = 0.75;
    strict_opts.branching_factor = 2;
    strict_opts.merge_criterion = OECluster::BitBirchMergeCriterion::Diameter;
    strict_opts.mode = OECluster::BitBirchMode::StrictParity;

    OECluster::BitBirchOptions fast_opts = strict_opts;
    fast_opts.mode = OECluster::BitBirchMode::Fast;

    const auto strict = OECluster::bitbirch_cluster(batch, strict_opts);
    const auto fast = OECluster::bitbirch_cluster(batch, fast_opts);

    EXPECT_EQ(strict.Labels(), fast.Labels());
    EXPECT_EQ(strict.Members(), fast.Members());
}

TEST(BitBirchFastTest, BitBirchClusterFastRoutesThroughEngine) {
    const auto batch = make_random_batch(5000, 64);
    OECluster::BitBirchOptions opts;
    opts.threshold = 0.5;
    opts.branching_factor = 50;
    opts.merge_criterion = OECluster::BitBirchMergeCriterion::Diameter;
    opts.mode = OECluster::BitBirchMode::Fast;

    const auto via_api = OECluster::bitbirch_cluster(batch, opts);

    OECluster::detail::BitBirchTree tree(opts);
    OECluster::detail::BitBirchTree::BuildFastTree(batch, opts, tree);
    const auto via_engine = tree.Result(batch.Spec(), batch.Size());

    EXPECT_EQ(via_api.Labels(), via_engine.Labels());
    EXPECT_EQ(via_api.Members(), via_engine.Members());
}

TEST(BitBirchFastTest, BitBirchClusterFastQualityEquivalentToStrict) {
    const auto batch = make_random_batch(3000, 64);  // > chunk target -> P >= 2
    const auto storage = tanimoto_storage(batch);

    OECluster::BitBirchOptions strict_opts;
    strict_opts.threshold = 0.5;
    strict_opts.branching_factor = 50;
    strict_opts.merge_criterion = OECluster::BitBirchMergeCriterion::Diameter;
    OECluster::BitBirchOptions fast_opts = strict_opts;
    fast_opts.mode = OECluster::BitBirchMode::Fast;

    const auto strict = OECluster::bitbirch_cluster(batch, strict_opts);
    const auto fast = OECluster::bitbirch_cluster(batch, fast_opts);

    const auto rs = OECluster::cluster_report(strict, storage, OECluster::ClusterReportOptions());
    const auto rf = OECluster::cluster_report(fast, storage, OECluster::ClusterReportOptions());

    // Non-degenerate in both modes.
    ASSERT_GE(rs.num_clusters, 2u);
    ASSERT_GE(rf.num_clusters, 2u);
    ASSERT_FALSE(std::isnan(rs.silhouette));
    ASSERT_FALSE(std::isnan(rf.silhouette));

    // Fast must be quality-equivalent: tighter-or-equal intra distance and
    // not-much-lower silhouette, cluster count within 15%. Tolerances mirror
    // the Python FAST_QUALITY_TOL constants (Task 6); keep the two in sync.
    EXPECT_LE(rf.mean_intra_distance, rs.mean_intra_distance + 0.05);
    EXPECT_GE(rf.silhouette, rs.silhouette - 0.05);
    EXPECT_LE(std::abs(static_cast<double>(rf.num_clusters) -
                       static_cast<double>(rs.num_clusters)),
              0.15 * static_cast<double>(rs.num_clusters));
}
