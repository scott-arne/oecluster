/**
 * @file test_consensus.cpp
 * @brief The three consensus kernels: values, extraction, strength, refusals.
 */
#include <gtest/gtest.h>

#include <cmath>
#include <cstdint>
#include <filesystem>
#include <functional>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "oecluster/Consensus.h"
#include "oecluster/StorageBackend.h"

using namespace OECluster;

namespace {

// Three members over four items. Member 0 splits {0,1} from {2,3}; member 1
// groups {0,1,2} against {3}; member 2 observes only {0,1,3} and splits
// {0,1} from {3}. Every co-association below is hand-computed from these.
struct Ensemble {
    std::vector<size_t> offsets{0, 4, 8, 11};
    std::vector<size_t> positions{0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 3};
    std::vector<int> labels{0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 1};
};

// A temp-file MMapStorage. The file is removed first, because MMapStorage
// reuses an existing file of the right size, and again on destruction so a
// failing test leaves nothing behind.
class TempMMap {
public:
    TempMMap(const std::string& name, size_t n)
        : path_(std::filesystem::temp_directory_path() / name) {
        std::filesystem::remove(path_);
        storage_ = std::make_unique<MMapStorage>(path_.string(), n);
    }
    ~TempMMap() {
        storage_.reset();
        std::filesystem::remove(path_);
    }
    TempMMap(const TempMMap&) = delete;
    TempMMap& operator=(const TempMMap&) = delete;

    MMapStorage& Storage() { return *storage_; }

private:
    std::filesystem::path path_;
    std::unique_ptr<MMapStorage> storage_;
};

/// Assert the refusal and the reason: a bare EXPECT_THROW would pass
/// whichever check won, which is exactly what these two cases pin down.
void ExpectRefusal(const std::function<void()>& call,
                   const std::string& needle) {
    try {
        call();
        ADD_FAILURE() << "expected std::invalid_argument mentioning "
                      << needle;
    } catch (const std::invalid_argument& error) {
        EXPECT_NE(std::string(error.what()).find(needle), std::string::npos)
            << "message was: " << error.what();
    }
}

}  // namespace

TEST(CoassociationTest, DistancesMatchTheHandComputedEnsemble) {
    const Ensemble ensemble;
    DenseStorage destination(4);
    const ConsensusMatrixSummary summary = coassociation_distances(
        4, ensemble.offsets, ensemble.positions, ensemble.labels, destination);

    EXPECT_EQ(summary.num_partitions, 3u);
    EXPECT_EQ(summary.unobserved_pairs, 0u);
    // co/obs: (0,1)=3/3 (0,2)=1/2 (0,3)=0/3 (1,2)=1/2 (1,3)=0/3 (2,3)=1/2
    EXPECT_DOUBLE_EQ(destination.Get(0, 1), 0.0);
    EXPECT_DOUBLE_EQ(destination.Get(0, 2), 0.5);
    EXPECT_DOUBLE_EQ(destination.Get(0, 3), 1.0);
    EXPECT_DOUBLE_EQ(destination.Get(1, 2), 0.5);
    EXPECT_DOUBLE_EQ(destination.Get(1, 3), 1.0);
    EXPECT_DOUBLE_EQ(destination.Get(2, 3), 0.5);
}

TEST(CoassociationTest, NoiseIsObservedButNeverCoClustered) {
    // One member of four items, two of them noise.
    const std::vector<size_t> offsets{0, 4};
    const std::vector<size_t> positions{0, 1, 2, 3};
    const std::vector<int> labels{0, 0, -1, -1};
    DenseStorage destination(4);
    const ConsensusMatrixSummary summary =
        coassociation_distances(4, offsets, positions, labels, destination);

    EXPECT_EQ(summary.unobserved_pairs, 0u);  // noise items were still observed
    EXPECT_DOUBLE_EQ(destination.Get(0, 1), 0.0);
    EXPECT_DOUBLE_EQ(destination.Get(2, 3), 1.0);
    EXPECT_DOUBLE_EQ(destination.Get(0, 2), 1.0);
}

TEST(CoassociationTest, PairsNoMemberObservedTogetherAreCounted) {
    // Two members covering disjoint halves: the four crossing pairs have no
    // evidence either way.
    const std::vector<size_t> offsets{0, 2, 4};
    const std::vector<size_t> positions{0, 1, 2, 3};
    const std::vector<int> labels{0, 0, 0, 0};
    DenseStorage destination(4);
    const ConsensusMatrixSummary summary =
        coassociation_distances(4, offsets, positions, labels, destination);

    EXPECT_EQ(summary.unobserved_pairs, 4u);
    EXPECT_DOUBLE_EQ(destination.Get(0, 1), 0.0);
    EXPECT_DOUBLE_EQ(destination.Get(2, 3), 0.0);
    EXPECT_DOUBLE_EQ(destination.Get(0, 2), 1.0);
}

TEST(CoassociationTest, MoreThanSixtyFourMembersUseEveryMaskWord) {
    // The membership bitmask is ceil(R / 64) words per item, so an ensemble
    // past 64 members is the first that must read and intersect a second
    // word. It is also the documented default: cluster_stability resamples
    // 100 times. Every member observes items 0, 1 and 2 and splits {0,1}
    // from {2}; only the last five, which live in the second word, observe
    // item 3 at all, two of them with {0,1} and three with {2}.
    std::vector<size_t> offsets{0};
    std::vector<size_t> positions;
    std::vector<int> labels;
    for (size_t member = 0; member < 70; ++member) {
        positions.insert(positions.end(), {0u, 1u, 2u});
        labels.insert(labels.end(), {0, 0, 1});
        if (member >= 65) {
            positions.push_back(3);
            labels.push_back(member < 67 ? 0 : 1);
        }
        offsets.push_back(positions.size());
    }

    DenseStorage destination(4);
    const ConsensusMatrixSummary summary =
        coassociation_distances(4, offsets, positions, labels, destination);

    EXPECT_EQ(summary.num_partitions, 70u);
    EXPECT_EQ(summary.unobserved_pairs, 0u);
    // co/obs: (0,1)=70/70 (0,2)=0/70 (1,2)=0/70, and item 3 is observed by
    // five members only, so its denominator is 5 and not 70.
    EXPECT_DOUBLE_EQ(destination.Get(0, 1), 0.0);
    EXPECT_DOUBLE_EQ(destination.Get(0, 2), 1.0);
    EXPECT_DOUBLE_EQ(destination.Get(1, 2), 1.0);
    EXPECT_DOUBLE_EQ(destination.Get(0, 3), 1.0 - 2.0 / 5.0);
    EXPECT_DOUBLE_EQ(destination.Get(1, 3), 1.0 - 2.0 / 5.0);
    EXPECT_DOUBLE_EQ(destination.Get(2, 3), 1.0 - 3.0 / 5.0);
}

TEST(CoassociationTest, ThreadCountDoesNotChangeTheResult) {
    // Sixteen items in one member of four clusters, so the accumulation
    // spreads over several chunks with more than one worker.
    std::vector<size_t> offsets{0, 16};
    std::vector<size_t> positions;
    std::vector<int> labels;
    for (size_t item = 0; item < 16; ++item) {
        positions.push_back(item);
        labels.push_back(static_cast<int>(item % 4));
    }
    DenseStorage one(16);
    DenseStorage two(16);
    DenseStorage automatic(16);
    ConsensusOptions single;
    single.num_threads = 1;
    single.chunk_size = 1;
    ConsensusOptions pair;
    pair.num_threads = 2;
    pair.chunk_size = 1;
    coassociation_distances(16, offsets, positions, labels, one, single);
    coassociation_distances(16, offsets, positions, labels, two, pair);
    coassociation_distances(16, offsets, positions, labels, automatic);
    for (size_t i = 0; i < 16; ++i) {
        for (size_t j = i + 1; j < 16; ++j) {
            EXPECT_DOUBLE_EQ(one.Get(i, j), two.Get(i, j));
            EXPECT_DOUBLE_EQ(one.Get(i, j), automatic.Get(i, j));
        }
    }
}

TEST(CoassociationTest, APrefilledDestinationIsZeroedFirst) {
    const Ensemble ensemble;
    DenseStorage reference(4);
    coassociation_distances(4, ensemble.offsets, ensemble.positions,
                            ensemble.labels, reference);

    // MMapStorage keeps an existing file of the right size without clearing
    // it, so arbitrary prior contents must not survive into the counts.
    TempMMap mapped("oecluster_test_consensus_reuse.bin", 4);
    for (size_t i = 0; i < 4; ++i) {
        for (size_t j = i + 1; j < 4; ++j) {
            mapped.Storage().Set(i, j, 7.5);
        }
    }
    coassociation_distances(4, ensemble.offsets, ensemble.positions,
                            ensemble.labels, mapped.Storage());
    for (size_t i = 0; i < 4; ++i) {
        for (size_t j = i + 1; j < 4; ++j) {
            EXPECT_DOUBLE_EQ(mapped.Storage().Get(i, j), reference.Get(i, j));
        }
    }

    // And running the kernel twice over one destination is the same.
    coassociation_distances(4, ensemble.offsets, ensemble.positions,
                            ensemble.labels, mapped.Storage());
    EXPECT_DOUBLE_EQ(mapped.Storage().Get(0, 1), reference.Get(0, 1));
}

TEST(CoassociationTest, AMemberOfTwoLargeClustersStillSpreadsOverThePool) {
    // The unit of work is a pair row, not a cluster: with cluster-sized
    // chunks this member's whole quadratic accumulation would land on one
    // worker, and the result would still be correct, so only the values are
    // asserted here while the shape of the work is the reason for the test.
    std::vector<size_t> offsets{0, 400};
    std::vector<size_t> positions;
    std::vector<int> labels;
    for (size_t item = 0; item < 400; ++item) {
        positions.push_back(item);
        labels.push_back(static_cast<int>(item % 2));
    }
    DenseStorage serial(400);
    DenseStorage parallel(400);
    ConsensusOptions one;
    one.num_threads = 1;
    ConsensusOptions many;
    many.num_threads = 8;
    coassociation_distances(400, offsets, positions, labels, serial, one);
    coassociation_distances(400, offsets, positions, labels, parallel, many);
    for (size_t i = 0; i < 400; ++i) {
        for (size_t j = i + 1; j < 400; ++j) {
            EXPECT_DOUBLE_EQ(serial.Get(i, j), parallel.Get(i, j));
        }
    }
    EXPECT_DOUBLE_EQ(serial.Get(0, 2), 0.0);   // same cluster
    EXPECT_DOUBLE_EQ(serial.Get(0, 1), 1.0);   // different clusters
}

TEST(CoassociationTest, ExtremeThreadAndChunkCountsStillDoTheWork) {
    // Both of these are accepted by the contract and both overflow a chunk
    // count if they reach the arithmetic unclamped: 4 * workers wraps to
    // zero, and ThreadPool's (range + chunk_size - 1) wraps to fewer chunks
    // than the range needs, which drops a whole pass in silence.
    const Ensemble ensemble;
    DenseStorage plain(4);
    DenseStorage extreme(4);
    ConsensusOptions wild;
    wild.num_threads = (SIZE_MAX / 4) + 1;  // 4 * workers == 0 if uncapped
    wild.chunk_size = SIZE_MAX;
    coassociation_distances(4, ensemble.offsets, ensemble.positions,
                            ensemble.labels, plain);
    coassociation_distances(4, ensemble.offsets, ensemble.positions,
                            ensemble.labels, extreme, wild);
    for (size_t i = 0; i < 4; ++i) {
        for (size_t j = i + 1; j < 4; ++j) {
            EXPECT_DOUBLE_EQ(plain.Get(i, j), extreme.Get(i, j));
        }
    }
    // Not merely equal to each other: a dropped finalization pass would
    // leave the raw co-association counts behind.
    EXPECT_DOUBLE_EQ(extreme.Get(0, 1), 0.0);
}

TEST(CoassociationTest, RefusesMalformedInput) {
    const Ensemble ensemble;
    DenseStorage destination(4);

    EXPECT_THROW(coassociation_distances(1, {0, 1}, {0}, {0}, destination),
                 std::invalid_argument);
    DenseStorage wrong_size(3);
    EXPECT_THROW(coassociation_distances(4, ensemble.offsets,
                                         ensemble.positions, ensemble.labels,
                                         wrong_size),
                 std::invalid_argument);
    SparseStorage sparse(4, 0.5);
    EXPECT_THROW(coassociation_distances(4, ensemble.offsets,
                                         ensemble.positions, ensemble.labels,
                                         sparse),
                 std::invalid_argument);
    ConsensusOptions zero_chunk;
    zero_chunk.chunk_size = 0;
    EXPECT_THROW(coassociation_distances(4, ensemble.offsets,
                                         ensemble.positions, ensemble.labels,
                                         destination, zero_chunk),
                 std::invalid_argument);
    // Refusal order: the destination's storage kind outranks the chunk size,
    // so a call that is wrong in both ways must name the storage.
    ExpectRefusal([&ensemble, &sparse, &zero_chunk] {
        coassociation_distances(4, ensemble.offsets, ensemble.positions,
                                ensemble.labels, sparse, zero_chunk);
    }, "contiguous");
    EXPECT_THROW(coassociation_distances(4, {}, {}, {}, destination),
                 std::invalid_argument);
    EXPECT_THROW(coassociation_distances(4, {1, 4}, {0, 1, 2, 3},
                                         {0, 0, 1, 1}, destination),
                 std::invalid_argument);  // offsets must start at 0
    EXPECT_THROW(coassociation_distances(4, {0, 3}, {0, 1, 2, 3},
                                         {0, 0, 1, 1}, destination),
                 std::invalid_argument);  // offsets must end at the length
    EXPECT_THROW(coassociation_distances(4, {0, 2, 1, 4}, {0, 1, 2, 3},
                                         {0, 0, 1, 1}, destination),
                 std::invalid_argument);  // offsets must not decrease
    // The whole offset structure is checked before any span is read: this
    // one ends correctly, so only the out-of-range middle entry refuses it,
    // and it must refuse before member 0's span reads positions[4]. The
    // message is asserted because an implementation that read out of bounds
    // first and threw afterwards would satisfy a bare EXPECT_THROW.
    ExpectRefusal([&destination] {
        coassociation_distances(4, {0, 5, 4}, {0, 1, 2, 3}, {0, 0, 1, 1},
                                destination);
    }, "exceeds");
    EXPECT_THROW(coassociation_distances(4, {0, 0, 4}, {0, 1, 2, 3},
                                         {0, 0, 1, 1}, destination),
                 std::invalid_argument);  // an empty member
    EXPECT_THROW(coassociation_distances(4, {0, 4}, {0, 1, 2, 3}, {0, 0, 1},
                                         destination),
                 std::invalid_argument);  // labels shorter than positions
    EXPECT_THROW(coassociation_distances(4, {0, 4}, {0, 1, 2, 9},
                                         {0, 0, 1, 1}, destination),
                 std::invalid_argument);  // position out of range
    EXPECT_THROW(coassociation_distances(4, {0, 4}, {0, 1, 1, 3},
                                         {0, 0, 1, 1}, destination),
                 std::invalid_argument);  // repeated position in one member
}

TEST(ConsensusComponentsTest, ThresholdDecidesTheComponents) {
    const Ensemble ensemble;
    DenseStorage matrix(4);
    coassociation_distances(4, ensemble.offsets, ensemble.positions,
                            ensemble.labels, matrix);

    // At 0.5 every stored pair but (0,3) and (1,3) is at or above the
    // threshold, and those two join through item 2.
    const std::vector<int> half = consensus_components(matrix, 0.5);
    EXPECT_EQ(half, (std::vector<int>{0, 0, 0, 0}));

    // At 1.0 only pairs co-clustered in every member that saw them merge.
    const std::vector<int> unanimous = consensus_components(matrix, 1.0);
    EXPECT_EQ(unanimous, (std::vector<int>{0, 0, 1, 2}));

    // At 0.0 everything merges, including pairs with no evidence.
    const std::vector<int> all = consensus_components(matrix, 0.0);
    EXPECT_EQ(all, (std::vector<int>{0, 0, 0, 0}));
}

TEST(ConsensusComponentsTest, LabelsFollowTheSmallestMemberPosition) {
    // Items 0 and 3 agree; 1 and 2 agree; the two groups never do.
    const std::vector<size_t> offsets{0, 4};
    const std::vector<size_t> positions{0, 1, 2, 3};
    const std::vector<int> labels{7, 9, 9, 7};
    DenseStorage matrix(4);
    coassociation_distances(4, offsets, positions, labels, matrix);
    EXPECT_EQ(consensus_components(matrix, 0.5), (std::vector<int>{0, 1, 1, 0}));
}

TEST(ConsensusComponentsTest, RefusesABadMatrixOrThreshold) {
    DenseStorage matrix(4);
    SparseStorage sparse(4, 0.5);
    EXPECT_THROW(consensus_components(sparse, 0.5), std::invalid_argument);
    DenseStorage single(1);
    EXPECT_THROW(consensus_components(single, 0.5), std::invalid_argument);
    EXPECT_THROW(consensus_components(matrix, -0.1), std::invalid_argument);
    EXPECT_THROW(consensus_components(matrix, 1.1), std::invalid_argument);
    EXPECT_THROW(consensus_components(matrix, std::nan("")),
                 std::invalid_argument);
}

TEST(ConsensusStrengthTest, MeansMatchTheHandComputedEnsemble) {
    const Ensemble ensemble;
    DenseStorage matrix(4);
    coassociation_distances(4, ensemble.offsets, ensemble.positions,
                            ensemble.labels, matrix);
    const std::vector<int> labels = consensus_components(matrix, 0.5);
    const ConsensusStrength strength = consensus_strength(matrix, labels);

    // Co-associations: (0,1)=1 (0,2)=.5 (0,3)=0 (1,2)=.5 (1,3)=0 (2,3)=.5
    ASSERT_EQ(strength.item_consensus.size(), 4u);
    EXPECT_DOUBLE_EQ(strength.item_consensus[0], 1.5 / 3.0);
    EXPECT_DOUBLE_EQ(strength.item_consensus[3], 0.5 / 3.0);
    ASSERT_EQ(strength.cluster_consensus.size(), 1u);
    EXPECT_DOUBLE_EQ(strength.cluster_consensus[0], 2.5 / 6.0);
}

TEST(ConsensusStrengthTest, SingletonsAndNoiseAreNotANumber) {
    DenseStorage matrix(3);
    matrix.Set(0, 1, 0.25);
    matrix.Set(0, 2, 1.0);
    matrix.Set(1, 2, 1.0);
    // Item 2 is noise; items 0 and 1 share a cluster; label 4 is a singleton.
    const ConsensusStrength strength =
        consensus_strength(matrix, std::vector<int>{0, 0, -1});
    EXPECT_DOUBLE_EQ(strength.item_consensus[0], 0.75);
    EXPECT_TRUE(std::isnan(strength.item_consensus[2]));
    ASSERT_EQ(strength.cluster_consensus.size(), 1u);
    EXPECT_DOUBLE_EQ(strength.cluster_consensus[0], 0.75);

    const ConsensusStrength singletons =
        consensus_strength(matrix, std::vector<int>{0, 1, 4});
    EXPECT_TRUE(std::isnan(singletons.item_consensus[0]));
    ASSERT_EQ(singletons.cluster_consensus.size(), 3u);
    EXPECT_TRUE(std::isnan(singletons.cluster_consensus[2]));
}

TEST(ConsensusStrengthTest, ExtremeThreadAndChunkCountsStillDoTheWork) {
    const Ensemble ensemble;
    DenseStorage matrix(4);
    coassociation_distances(4, ensemble.offsets, ensemble.positions,
                            ensemble.labels, matrix);
    const std::vector<int> labels{0, 0, 1, 1};
    ConsensusOptions wild;
    wild.num_threads = (SIZE_MAX / 4) + 1;
    wild.chunk_size = SIZE_MAX;
    const auto plain = consensus_strength(matrix, labels);
    const auto extreme = consensus_strength(matrix, labels, wild);
    for (size_t i = 0; i < 4; ++i) {
        EXPECT_DOUBLE_EQ(plain.item_consensus[i], extreme.item_consensus[i]);
    }
    EXPECT_DOUBLE_EQ(extreme.item_consensus[0], 1.0);
    EXPECT_DOUBLE_EQ(extreme.cluster_consensus[0], 1.0);
}

TEST(ConsensusStrengthTest, RefusesAMismatchedLabelCount) {
    DenseStorage matrix(4);
    SparseStorage sparse(4, 0.5);
    EXPECT_THROW(consensus_strength(matrix, std::vector<int>{0, 0, 0}),
                 std::invalid_argument);
    EXPECT_THROW(consensus_strength(sparse, std::vector<int>{0, 0, 0, 0}),
                 std::invalid_argument);
    ConsensusOptions zero_chunk;
    zero_chunk.chunk_size = 0;
    EXPECT_THROW(
        consensus_strength(matrix, std::vector<int>{0, 0, 0, 0}, zero_chunk),
        std::invalid_argument);
}
