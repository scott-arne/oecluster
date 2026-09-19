#include <gtest/gtest.h>

#include <stdexcept>
#include <vector>

#include "oecluster/StorageBackend.h"

#include "../../src/clustering/MaxMinKernel.h"

using namespace OECluster;

namespace {

// Four items on a line at 0, 0.5, 10, 11: farthest-first from item 0 must
// reach the far pair before the near one.
//
// Item 1 sits at 0.5 rather than 1.0 so the third selection is decided by
// distance instead of by the tie rule. At 1.0, items 1 and 2 would both be 1.0
// from the selected set {0, 3} and the smaller-index rule would pick item 1 --
// a correct answer, but one that tests the tie rule rather than farthest-first
// coverage. At 0.5 item 2 wins outright at 1.0 against item 1's 0.5.
// Every position is dyadic, so every distance below is exact.
DenseStorage MakeLineStorage() {
    DenseStorage storage(4);
    storage.Set(0, 1, 0.5);
    storage.Set(0, 2, 10.0);
    storage.Set(0, 3, 11.0);
    storage.Set(1, 2, 9.5);
    storage.Set(1, 3, 10.5);
    storage.Set(2, 3, 1.0);
    return storage;
}

DenseStorage MakeAllZeroStorage(size_t n) {
    DenseStorage storage(n);
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i + 1; j < n; ++j) {
            storage.Set(i, j, 0.0);
        }
    }
    return storage;
}

// Four items at 0, 4, 9, 10: the pair (2, 3) is clustered together (distance 1.0)
// while item 1 is midway between them. This fixture separates the MaxMin recurrence
// from "farthest from the seed" — the nearest-distance update after selecting item 3
// lowers item 2's distance from 9.0 to 1.0, which causes item 1 (at 4.0) to be picked
// next instead of item 2. Without that update item 2 would be chosen. Every position
// is an integer, so every distance is exact.
DenseStorage MakeClusteredLineStorage() {
    DenseStorage storage(4);
    storage.Set(0, 1, 4.0);
    storage.Set(0, 2, 9.0);
    storage.Set(0, 3, 10.0);
    storage.Set(1, 2, 5.0);
    storage.Set(1, 3, 6.0);
    storage.Set(2, 3, 1.0);
    return storage;
}

}  // namespace

TEST(MaxMinKernelTest, SelectsFarthestFirstFromTheSeed) {
    const DenseStorage storage = MakeLineStorage();

    const std::vector<size_t> selection =
        detail::maxmin_select_from(storage.Data(), storage.NumSamples(), 3, 0);

    EXPECT_EQ(selection, std::vector<size_t>({0, 3, 2}));
}

TEST(MaxMinKernelTest, LowersTheNearestDistanceAfterEachSelection) {
    const DenseStorage storage = MakeClusteredLineStorage();

    const std::vector<size_t> selection =
        detail::maxmin_select_from(storage.Data(), storage.NumSamples(), 3, 0);

    EXPECT_EQ(selection, std::vector<size_t>({0, 3, 1}));
}

TEST(MaxMinKernelTest, HonorsTheSeed) {
    const DenseStorage storage = MakeLineStorage();

    const std::vector<size_t> selection =
        detail::maxmin_select_from(storage.Data(), storage.NumSamples(), 2, 2);

    EXPECT_EQ(selection, std::vector<size_t>({2, 0}));
}

// On an all-zero matrix every candidate ties at every step. Without the
// selected mask the smaller-index rule would return the seed `count` times.
TEST(MaxMinKernelTest, ExcludesSelectedItemsWhenEveryDistanceTies) {
    const DenseStorage storage = MakeAllZeroStorage(5);

    const std::vector<size_t> selection =
        detail::maxmin_select_from(storage.Data(), storage.NumSamples(), 3, 2);

    EXPECT_EQ(selection, std::vector<size_t>({2, 0, 1}));
}

TEST(MaxMinKernelTest, SelectsEveryItemWhenCountEqualsTheItemCount) {
    const DenseStorage storage = MakeAllZeroStorage(3);

    const std::vector<size_t> selection =
        detail::maxmin_select_from(storage.Data(), storage.NumSamples(), 3, 1);

    EXPECT_EQ(selection, std::vector<size_t>({1, 0, 2}));
}

TEST(MaxMinKernelTest, RejectsAnOutOfRangeCount) {
    const DenseStorage storage = MakeLineStorage();

    EXPECT_THROW(
        detail::maxmin_select_from(storage.Data(), storage.NumSamples(), 0, 0),
        std::invalid_argument);
    EXPECT_THROW(
        detail::maxmin_select_from(storage.Data(), storage.NumSamples(), 5, 0),
        std::invalid_argument);
}

TEST(MaxMinKernelTest, RejectsAnOutOfRangeSeed) {
    const DenseStorage storage = MakeLineStorage();

    EXPECT_THROW(
        detail::maxmin_select_from(storage.Data(), storage.NumSamples(), 2, 4),
        std::invalid_argument);
}
