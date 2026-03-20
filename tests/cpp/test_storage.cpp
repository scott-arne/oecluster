#include <gtest/gtest.h>
#include <filesystem>
#include "oecluster/StorageBackend.h"

using namespace OECluster;

TEST(DenseStorageTest, ConstructorSetsSize) {
    DenseStorage storage(5);
    EXPECT_EQ(storage.NumItems(), 5);
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
