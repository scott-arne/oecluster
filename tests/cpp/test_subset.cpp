/**
 * @file test_subset.cpp
 * @brief take_pairs and take_fingerprints: gathers, refusals and invariance.
 */
#include <gtest/gtest.h>

#include <algorithm>
#include <cstdint>
#include <filesystem>
#include <initializer_list>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "oecluster/StorageBackend.h"
#include "oecluster/Subset.h"
#include "oefp/batch.h"
#include "oefp/fingerprint.h"

using namespace OECluster;

namespace {

// d(i, j) = 100 i + j for i < j: every pair carries its own identity, so a
// gather that lands a distance in the wrong slot is caught by value.
double pair_value(size_t i, size_t j) {
    if (i > j) {
        std::swap(i, j);
    }
    return 100.0 * static_cast<double>(i) + static_cast<double>(j);
}

DenseStorage make_dense(size_t n) {
    DenseStorage storage(n);
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i + 1; j < n; ++j) {
            storage.Set(i, j, pair_value(i, j));
        }
    }
    return storage;
}

// SparseStorage owns a mutex and cannot be returned by value.
void fill_sparse(SparseStorage& storage) {
    storage.Set(0, 1, 0.2);
    storage.Set(2, 3, 0.2);
    storage.Set(4, 5, 0.2);
    storage.Set(1, 4, 0.4);
    storage.Finalize();
}

void expect_gathered(const StorageBackend& source,
                     const std::vector<size_t>& indices,
                     const StorageBackend& destination) {
    ASSERT_EQ(destination.NumSamples(), indices.size());
    for (size_t a = 0; a < indices.size(); ++a) {
        for (size_t b = a + 1; b < indices.size(); ++b) {
            EXPECT_DOUBLE_EQ(destination.Get(a, b),
                             source.Get(indices[a], indices[b]))
                << "pair " << a << "," << b;
        }
    }
}

// A temp-file MMapStorage holding pair_value distances. The file is removed
// first, because MMapStorage reuses an existing file of the right size, and
// removed again on destruction so a failing test leaves nothing behind.
class TempMMap {
public:
    TempMMap(const std::string& name, size_t n)
        : path_(std::filesystem::temp_directory_path() / name) {
        std::filesystem::remove(path_);
        storage_ = std::make_unique<MMapStorage>(path_.string(), n);
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = i + 1; j < n; ++j) {
                storage_->Set(i, j, pair_value(i, j));
            }
        }
    }
    ~TempMMap() {
        storage_.reset();
        std::filesystem::remove(path_);
    }
    TempMMap(const TempMMap&) = delete;
    TempMMap& operator=(const TempMMap&) = delete;

    const MMapStorage& Storage() const { return *storage_; }

private:
    std::filesystem::path path_;
    std::unique_ptr<MMapStorage> storage_;
};

OEFP::OEFP make_fp(size_t size_bits, std::initializer_list<size_t> on_bits) {
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

}  // namespace

TEST(TakePairsTest, DenseGatherFollowsTheIndexOrder) {
    DenseStorage source = make_dense(8);
    const std::vector<size_t> indices{6, 1, 4, 0};
    DenseStorage destination(4);
    take_pairs(source, indices, destination);
    expect_gathered(source, indices, destination);
    EXPECT_DOUBLE_EQ(destination.Get(0, 1), pair_value(6, 1));
}

TEST(TakePairsTest, MMapSourceGathersIntoDense) {
    TempMMap mmap("oecluster_test_subset_source.bin", 7);
    const std::vector<size_t> indices{5, 2, 6};
    DenseStorage destination(3);
    take_pairs(mmap.Storage(), indices, destination);
    expect_gathered(mmap.Storage(), indices, destination);
}

TEST(TakePairsTest, SparseGatherKeepsOnlyStoredPairsAndTheCutoff) {
    SparseStorage source(6, 0.5);
    fill_sparse(source);
    const std::vector<size_t> indices{4, 1, 0, 3};
    SparseStorage destination(4, 0.5);
    take_pairs(source, indices, destination);
    EXPECT_DOUBLE_EQ(destination.Get(1, 2), 0.2);  // source (1, 0)
    EXPECT_DOUBLE_EQ(destination.Get(0, 1), 0.4);  // source (4, 1)
    EXPECT_DOUBLE_EQ(destination.Get(0, 3), 0.0);  // source (4, 3): never stored
    EXPECT_DOUBLE_EQ(destination.Get(2, 3), 0.0);  // source (0, 3): never stored
    EXPECT_EQ(destination.Entries().size(), 2u);
    EXPECT_DOUBLE_EQ(destination.Cutoff(), 0.5);
}

TEST(TakePairsTest, DenseDestinationIsFullyOverwritten) {
    DenseStorage source = make_dense(5);
    const std::vector<size_t> indices{4, 2, 0};
    DenseStorage destination(3);
    for (size_t a = 0; a < 3; ++a) {
        for (size_t b = a + 1; b < 3; ++b) {
            destination.Set(a, b, 99.0);
        }
    }
    take_pairs(source, indices, destination);
    for (size_t a = 0; a < 3; ++a) {
        for (size_t b = a + 1; b < 3; ++b) {
            EXPECT_NE(destination.Get(a, b), 99.0);
        }
    }
    expect_gathered(source, indices, destination);
}

TEST(TakePairsTest, ThreadCountDoesNotChangeTheResult) {
    DenseStorage source = make_dense(100);
    // m = 64 descending positions: with two workers the chunk rule gives
    // ceil(64 / 8) = 8 rows per chunk, so several chunks share the work.
    std::vector<size_t> indices;
    for (size_t i = 0; i < 64; ++i) {
        indices.push_back(99 - i);
    }
    DenseStorage one(64);
    DenseStorage two(64);
    DenseStorage automatic(64);
    take_pairs(source, indices, one, 1);
    take_pairs(source, indices, two, 2);
    take_pairs(source, indices, automatic, 0);
    for (size_t a = 0; a < 64; ++a) {
        for (size_t b = a + 1; b < 64; ++b) {
            EXPECT_DOUBLE_EQ(one.Get(a, b), two.Get(a, b));
            EXPECT_DOUBLE_EQ(one.Get(a, b), automatic.Get(a, b));
        }
    }
    expect_gathered(source, indices, one);
}

TEST(TakePairsTest, RefusesBadIndicesAndShapes) {
    DenseStorage source = make_dense(5);
    DenseStorage destination(2);
    EXPECT_THROW(take_pairs(source, {}, destination), std::invalid_argument);
    EXPECT_THROW(take_pairs(source, {0, 5}, destination), std::invalid_argument);
    EXPECT_THROW(take_pairs(source, {1, 1}, destination), std::invalid_argument);
    EXPECT_THROW(take_pairs(source, {0, 1, 2}, destination), std::invalid_argument);
    EXPECT_THROW(take_pairs(source, {0, 1}, destination, 0, 0), std::invalid_argument);
    EXPECT_THROW(take_pairs(source, {0, 1, 2, 3, 4}, source), std::invalid_argument);
}

TEST(TakePairsTest, RefusesStorageKindAndCutoffMismatches) {
    DenseStorage dense_source = make_dense(5);
    SparseStorage sparse_source(6, 0.5);
    fill_sparse(sparse_source);

    SparseStorage sparse_for_dense(2, 0.5);
    EXPECT_THROW(take_pairs(dense_source, {0, 1}, sparse_for_dense),
                 std::invalid_argument);
    DenseStorage dense_for_sparse(2);
    EXPECT_THROW(take_pairs(sparse_source, {0, 1}, dense_for_sparse),
                 std::invalid_argument);

    SparseStorage lower(2, 0.4);
    EXPECT_THROW(take_pairs(sparse_source, {0, 1}, lower), std::invalid_argument);
    SparseStorage higher(2, 0.6);
    EXPECT_THROW(take_pairs(sparse_source, {0, 1}, higher), std::invalid_argument);

    SparseStorage stale(2, 0.5);
    stale.Set(0, 1, 0.1);
    stale.Finalize();
    EXPECT_THROW(take_pairs(sparse_source, {0, 1}, stale), std::invalid_argument);

    // A write still in a thread buffer is invisible to Entries() until
    // Finalize(); the gather must refuse it all the same.
    SparseStorage buffered(2, 0.5);
    buffered.Set(0, 1, 0.1);
    EXPECT_THROW(take_pairs(sparse_source, {0, 1}, buffered),
                 std::invalid_argument);

    SparseStorage fresh(2, 0.5);
    EXPECT_NO_THROW(take_pairs(sparse_source, {0, 1}, fresh));
    EXPECT_DOUBLE_EQ(fresh.Get(0, 1), 0.2);
}

TEST(TakePairsTest, RefusesAMemoryMappedDestination) {
    TempMMap source("oecluster_test_subset_alias_source.bin", 4);
    const std::filesystem::path path =
        std::filesystem::temp_directory_path()
        / "oecluster_test_subset_alias_destination.bin";
    std::filesystem::remove(path);
    {
        MMapStorage destination(path.string(), 2);
        EXPECT_THROW(take_pairs(source.Storage(), {1, 0}, destination),
                     std::invalid_argument);
    }
    std::filesystem::remove(path);
}

TEST(TakeFingerprintsTest, CopiesRowsInOrderAndKeepsTheSpec) {
    OEFP::OEFPBatch batch = make_batch({
        make_fp(128, {1, 2}), make_fp(128, {3, 64}), make_fp(128, {127}),
        make_fp(128, {5, 6, 7}), make_fp(128, {0})});
    const std::vector<size_t> indices{3, 0, 4};
    OEFP::OEFPBatch subset = take_fingerprints(batch, indices);
    ASSERT_EQ(subset.Size(), 3u);
    EXPECT_TRUE(subset.Spec() == batch.Spec());
    ASSERT_EQ(subset.WordsPerFingerprint(), batch.WordsPerFingerprint());
    for (size_t a = 0; a < indices.size(); ++a) {
        const std::uint64_t* got = subset.RowWords(a);
        const std::uint64_t* want = batch.RowWords(indices[a]);
        for (size_t w = 0; w < batch.WordsPerFingerprint(); ++w) {
            EXPECT_EQ(got[w], want[w]) << "row " << a << " word " << w;
        }
        EXPECT_EQ(subset.PopCount(a), batch.PopCount(indices[a]));
    }
}

TEST(TakeFingerprintsTest, RefusesBadIndices) {
    OEFP::OEFPBatch batch = make_batch({
        make_fp(64, {1}), make_fp(64, {2}), make_fp(64, {3})});
    EXPECT_THROW(take_fingerprints(batch, {}), std::invalid_argument);
    EXPECT_THROW(take_fingerprints(batch, {0, 3}), std::invalid_argument);
    EXPECT_THROW(take_fingerprints(batch, {2, 2}), std::invalid_argument);
}
