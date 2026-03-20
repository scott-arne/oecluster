#include <gtest/gtest.h>
#include <cmath>
#include "oecluster/PDist.h"
#include "oecluster/StorageBackend.h"
#include "oecluster/DistanceMetric.h"

using namespace OECluster;

/// Mock metric: distance = |i - j| / n (always in [0, 1])
class MockMetric : public DistanceMetric {
public:
    explicit MockMetric(size_t n) : n_(n) {}

    double Distance(size_t i, size_t j) override {
        return static_cast<double>(std::abs(static_cast<int>(i) - static_cast<int>(j)))
               / static_cast<double>(n_);
    }

    std::unique_ptr<DistanceMetric> Clone() const override {
        return std::make_unique<MockMetric>(n_);
    }

    size_t Size() const override { return n_; }
    std::string Name() const override { return "mock"; }

private:
    size_t n_;
};

TEST(PDistTest, BasicComputation) {
    MockMetric metric(4);
    DenseStorage storage(4);
    pdist(metric, storage);

    // Verify all pairs
    EXPECT_DOUBLE_EQ(storage.Get(0, 1), 0.25);
    EXPECT_DOUBLE_EQ(storage.Get(0, 2), 0.50);
    EXPECT_DOUBLE_EQ(storage.Get(0, 3), 0.75);
    EXPECT_DOUBLE_EQ(storage.Get(1, 2), 0.25);
    EXPECT_DOUBLE_EQ(storage.Get(1, 3), 0.50);
    EXPECT_DOUBLE_EQ(storage.Get(2, 3), 0.25);
}

TEST(PDistTest, MultiThreaded) {
    MockMetric metric(100);
    DenseStorage storage(100);
    PDistOptions opts;
    opts.num_threads = 4;
    opts.chunk_size = 64;
    pdist(metric, storage, opts);

    // Spot check
    EXPECT_DOUBLE_EQ(storage.Get(0, 50), 0.50);
    EXPECT_DOUBLE_EQ(storage.Get(99, 0), 0.99);
}

TEST(PDistTest, WithCutoff) {
    MockMetric metric(4);
    SparseStorage storage(4, 0.3);
    PDistOptions opts;
    opts.cutoff = 0.3;
    pdist(metric, storage, opts);

    // Only pairs with distance <= 0.3 should be stored
    // (0,1)=0.25, (1,2)=0.25, (2,3)=0.25 are <= 0.3
    // (0,2)=0.50, (0,3)=0.75, (1,3)=0.50 are > 0.3
    auto& entries = storage.Entries();
    EXPECT_EQ(entries.size(), 3);
}

TEST(PDistTest, ProgressCallback) {
    MockMetric metric(10);
    DenseStorage storage(10);
    PDistOptions opts;
    opts.chunk_size = 5;
    size_t last_completed = 0;
    size_t callback_count = 0;
    opts.progress = [&](size_t completed, size_t total) {
        EXPECT_GE(completed, last_completed);
        EXPECT_EQ(total, 45);  // 10*9/2
        last_completed = completed;
        callback_count++;
    };
    pdist(metric, storage, opts);
    EXPECT_GT(callback_count, 0);
}
