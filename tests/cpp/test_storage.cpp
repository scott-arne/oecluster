#include <gtest/gtest.h>
#include <filesystem>
#include <stdexcept>
#include "oecluster/StorageBackend.h"

using namespace OECluster;

TEST(DenseStorageTest, ConstructorSetsSize) {
    DenseStorage storage(5);
    EXPECT_EQ(storage.NumSamples(), 5);
    EXPECT_EQ(storage.NumPairs(), 10);  // 5*4/2
}

TEST(DenseStorageTest, SetAndGet) {
    DenseStorage storage(4);
    storage.Set(0, 1, 0.5);
    storage.Set(0, 2, 0.8);
    storage.Set(1, 3, 0.2);
    EXPECT_DOUBLE_EQ(storage.Get(0, 1), 0.5);
    EXPECT_DOUBLE_EQ(storage.Get(0, 2), 0.8);
    EXPECT_DOUBLE_EQ(storage.Get(1, 3), 0.2);
}

TEST(DenseStorageTest, GetSymmetric) {
    DenseStorage storage(4);
    storage.Set(0, 1, 0.5);
    EXPECT_DOUBLE_EQ(storage.Get(1, 0), 0.5);
}

TEST(DenseStorageTest, GetDiagonalIsZero) {
    DenseStorage storage(4);
    EXPECT_DOUBLE_EQ(storage.Get(0, 0), 0.0);
    EXPECT_DOUBLE_EQ(storage.Get(2, 2), 0.0);
}

TEST(DenseStorageTest, DataPointerNotNull) {
    DenseStorage storage(4);
    EXPECT_NE(storage.Data(), nullptr);
}

TEST(DenseStorageTest, DataPointerMatchesCondensedIndex) {
    DenseStorage storage(4);  // N=4, pairs=6
    storage.Set(0, 1, 1.0);
    storage.Set(0, 2, 2.0);
    storage.Set(0, 3, 3.0);
    storage.Set(1, 2, 4.0);
    storage.Set(1, 3, 5.0);
    storage.Set(2, 3, 6.0);
    const double* data = storage.Data();
    // Condensed order: (0,1), (0,2), (0,3), (1,2), (1,3), (2,3)
    EXPECT_DOUBLE_EQ(data[0], 1.0);
    EXPECT_DOUBLE_EQ(data[1], 2.0);
    EXPECT_DOUBLE_EQ(data[2], 3.0);
    EXPECT_DOUBLE_EQ(data[3], 4.0);
    EXPECT_DOUBLE_EQ(data[4], 5.0);
    EXPECT_DOUBLE_EQ(data[5], 6.0);
}

TEST(DenseStorageTest, InitializedToZero) {
    DenseStorage storage(4);
    for (size_t i = 0; i < storage.NumPairs(); ++i) {
        EXPECT_DOUBLE_EQ(storage.Data()[i], 0.0);
    }
}

TEST(DenseStorageTest, GetOutOfRangeThrows) {
    DenseStorage storage(2);
    EXPECT_THROW(storage.Get(0, 2), std::out_of_range);
    EXPECT_THROW(storage.Get(1, 1399), std::out_of_range);
    EXPECT_THROW(storage.Get(500, 1399), std::out_of_range);
}

TEST(DenseStorageTest, GetOutOfRangeDiagonalThrows) {
    // The diagonal shortcut used to answer 0.0 for an index naming no item.
    DenseStorage storage(2);
    EXPECT_THROW(storage.Get(1000, 1000), std::out_of_range);
}

TEST(DenseStorageTest, GetOutOfRangeMessageNamesIndexAndCount) {
    DenseStorage storage(2);
    try {
        storage.Get(1, 1399);
        FAIL() << "Expected std::out_of_range";
    } catch (const std::out_of_range& e) {
        EXPECT_STREQ(e.what(),
                     "DenseStorage index 1399 is outside the storage range of 2 samples");
    }
}

TEST(DenseStorageTest, SetOutOfRangeThrowsRatherThanAliasingAnotherPair) {
    // CondensedIndex(5, 2, 5) == 9 == CondensedIndex(5, 3, 4): just past the
    // end, the computed index lands on a real pair, so the write used to
    // succeed and replace a distance the caller never named.
    DenseStorage storage(5);
    storage.Set(3, 4, 0.25);
    EXPECT_THROW(storage.Set(2, 5, 0.99), std::out_of_range);
    EXPECT_DOUBLE_EQ(storage.Get(3, 4), 0.25);
}

TEST(DenseStorageTest, SetFarOutOfRangeThrowsRatherThanWritingPastTheBuffer) {
    // Far past the end there is no pair to collide with: the index leaves the
    // allocation and the write corrupts whatever follows it, undiagnosed.
    DenseStorage storage(5);
    EXPECT_THROW(storage.Set(0, 5000, 1.0), std::out_of_range);
}

TEST(DenseStorageTest, SetOutOfRangeMessageNamesIndexAndCount) {
    DenseStorage storage(2);
    try {
        storage.Set(1, 1399, 0.5);
        FAIL() << "Expected std::out_of_range";
    } catch (const std::out_of_range& e) {
        EXPECT_STREQ(e.what(),
                     "DenseStorage index 1399 is outside the storage range of 2 samples");
    }
}

TEST(DenseStorageTest, SetDiagonalThrowsRatherThanAliasingAnotherPair) {
    // CondensedIndex(5, 2, 2) == 6 == CondensedIndex(5, 1, 4). Get answers the
    // diagonal 0.0 without consulting storage, so the diagonal owns no slot to
    // write; the condensed formula hands the write a neighbour's slot instead.
    DenseStorage storage(5);
    storage.Set(1, 4, 0.75);
    EXPECT_THROW(storage.Set(2, 2, 0.99), std::invalid_argument);
    EXPECT_DOUBLE_EQ(storage.Get(1, 4), 0.75);
}

TEST(DenseStorageTest, SetDiagonalMessageNamesTheIndex) {
    DenseStorage storage(5);
    try {
        storage.Set(3, 3, 0.5);
        FAIL() << "Expected std::invalid_argument";
    } catch (const std::invalid_argument& e) {
        EXPECT_STREQ(e.what(),
                     "DenseStorage cannot store the diagonal pair (3, 3); "
                     "only distances between distinct items are stored");
    }
}

TEST(MMapStorageTest, CreateAndWrite) {
    auto path = std::filesystem::temp_directory_path() / "test_mmap.bin";
    {
        MMapStorage storage(path.string(), 4);
        storage.Set(0, 1, 0.5);
        storage.Set(2, 3, 0.8);
        EXPECT_DOUBLE_EQ(storage.Get(0, 1), 0.5);
        EXPECT_DOUBLE_EQ(storage.Get(2, 3), 0.8);
        EXPECT_NE(storage.Data(), nullptr);
    }
    EXPECT_TRUE(std::filesystem::exists(path));
    std::filesystem::remove(path);
}

TEST(MMapStorageTest, PersistsAfterDestruction) {
    auto path = std::filesystem::temp_directory_path() / "test_mmap_persist.bin";
    {
        MMapStorage storage(path.string(), 4);
        storage.Set(0, 1, 0.5);
        storage.Set(1, 2, 0.9);
    }
    {
        MMapStorage storage(path.string(), 4);
        EXPECT_DOUBLE_EQ(storage.Get(0, 1), 0.5);
        EXPECT_DOUBLE_EQ(storage.Get(1, 2), 0.9);
    }
    std::filesystem::remove(path);
}

TEST(MMapStorageTest, GetOutOfRangeThrows) {
    auto path = std::filesystem::temp_directory_path() / "test_mmap_range.bin";
    {
        MMapStorage storage(path.string(), 2);
        EXPECT_THROW(storage.Get(0, 2), std::out_of_range);
        EXPECT_THROW(storage.Get(1000, 1000), std::out_of_range);
    }
    std::filesystem::remove(path);
}

TEST(MMapStorageTest, SetOutOfRangeThrows) {
    // The mapping is sized for exactly NumPairs doubles, so an unchecked write
    // past the end is a write past the mapping, not merely past a vector.
    auto path = std::filesystem::temp_directory_path() / "test_mmap_set_range.bin";
    {
        MMapStorage storage(path.string(), 2);
        EXPECT_THROW(storage.Set(0, 2, 1.0), std::out_of_range);
        // Out of range on the diagonal is reported as a range error: the
        // indices name no item at all, which is the more basic complaint.
        EXPECT_THROW(storage.Set(1000, 1000, 1.0), std::out_of_range);
        EXPECT_THROW(storage.Set(1, 1, 1.0), std::invalid_argument);
    }
    std::filesystem::remove(path);
}

TEST(SparseStorageTest, StoresOnlyBelowCutoff) {
    SparseStorage storage(4, 0.5);  // cutoff = 0.5
    storage.Set(0, 1, 0.3);  // below cutoff, stored
    storage.Set(0, 2, 0.8);  // above cutoff, ignored
    storage.Set(1, 2, 0.5);  // at cutoff, stored
    storage.Finalize();
    auto& entries = storage.Entries();
    EXPECT_EQ(entries.size(), 2);
}

TEST(SparseStorageTest, DataReturnsNullptr) {
    SparseStorage storage(4, 0.5);
    EXPECT_EQ(storage.Data(), nullptr);
}

TEST(SparseStorageTest, GetReturnsStoredValues) {
    SparseStorage storage(4, 0.5);
    storage.Set(0, 1, 0.3);
    storage.Finalize();
    EXPECT_DOUBLE_EQ(storage.Get(0, 1), 0.3);
}

TEST(SparseStorageTest, GetReturnsZeroForUnstored) {
    SparseStorage storage(4, 0.5);
    storage.Finalize();
    EXPECT_DOUBLE_EQ(storage.Get(0, 1), 0.0);
}

TEST(SparseStorageTest, GetOutOfRangeThrowsRatherThanMissing) {
    // "Not stored" answers 0.0; "not an item" must not borrow that answer.
    SparseStorage storage(4, 0.5);
    storage.Finalize();
    EXPECT_THROW(storage.Get(0, 4), std::out_of_range);
    EXPECT_THROW(storage.Get(1000, 1000), std::out_of_range);
}

TEST(SparseStorageTest, SetOutOfRangeThrowsWhicheverSideOfTheCutoff) {
    // The cutoff shortcut returns before anything is stored. Checking the
    // indices behind it would make the diagnosis depend on the value the bad
    // call happened to carry: below the cutoff a refusal, above it silence.
    SparseStorage storage(4, 0.5);
    EXPECT_THROW(storage.Set(0, 4, 0.3), std::out_of_range);
    EXPECT_THROW(storage.Set(0, 4, 0.9), std::out_of_range);
    EXPECT_THROW(storage.Set(2, 2, 0.3), std::invalid_argument);
    EXPECT_THROW(storage.Set(2, 2, 0.9), std::invalid_argument);
}

TEST(SparseStorageTest, SetOutOfRangeStoresNothing) {
    // A refused Set must not leave a half-written entry behind for Finalize to
    // pick up -- an entry whose condensed index collides with a real pair.
    SparseStorage storage(4, 0.5);
    EXPECT_THROW(storage.Set(0, 4, 0.3), std::out_of_range);
    storage.Finalize();
    EXPECT_EQ(storage.Entries().size(), 0u);
}
