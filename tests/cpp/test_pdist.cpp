#include <gtest/gtest.h>
#include <cmath>
#include "oecluster/PDist.h"
#include "oecluster/StorageBackend.h"
#include "oecluster/PairwiseComparison.h"

using namespace OECluster;

/// Mock comparison: distance = |i - j| / n (always in [0, 1])
class MockComparison : public PairwiseComparison {
public:
    explicit MockComparison(size_t n) : n_(n) {}

    double Compare(size_t i, size_t j) override {
        return static_cast<double>(std::abs(static_cast<int>(i) - static_cast<int>(j)))
               / static_cast<double>(n_);
    }

    std::unique_ptr<PairwiseComparison> Clone() const override {
        return std::make_unique<MockComparison>(n_);
    }

    size_t Size() const override { return n_; }
    std::string ComparisonName() const override { return "mock"; }

private:
    size_t n_;
};

class BulkPDistComparison : public MockComparison {
public:
    BulkPDistComparison() : MockComparison(2) {}

    bool TryPDist(StorageBackend& storage, const PDistOptions& options) override {
        storage.Set(0, 1, 0.125);
        if (options.progress) {
            options.progress(1, 1);
        }
        used_bulk = true;
        return true;
    }

    double Compare(size_t, size_t) override {
        used_compare = true;
        return 1.0;
    }

    bool used_bulk = false;
    bool used_compare = false;
};

TEST(PDistTest, BasicComputation) {
    MockComparison comparison(4);
    DenseStorage storage(4);
    pdist(comparison, storage);

    // Verify all pairs
    EXPECT_DOUBLE_EQ(storage.Get(0, 1), 0.25);
    EXPECT_DOUBLE_EQ(storage.Get(0, 2), 0.50);
    EXPECT_DOUBLE_EQ(storage.Get(0, 3), 0.75);
    EXPECT_DOUBLE_EQ(storage.Get(1, 2), 0.25);
    EXPECT_DOUBLE_EQ(storage.Get(1, 3), 0.50);
    EXPECT_DOUBLE_EQ(storage.Get(2, 3), 0.25);
}

TEST(PDistTest, MultiThreaded) {
    MockComparison comparison(100);
    DenseStorage storage(100);
    PDistOptions opts;
    opts.num_threads = 4;
    opts.chunk_size = 64;
    pdist(comparison, storage, opts);

    // Spot check
    EXPECT_DOUBLE_EQ(storage.Get(0, 50), 0.50);
    EXPECT_DOUBLE_EQ(storage.Get(99, 0), 0.99);
}

TEST(PDistTest, WithCutoff) {
    MockComparison comparison(4);
    SparseStorage storage(4, 0.3);
    PDistOptions opts;
    opts.cutoff = 0.3;
    pdist(comparison, storage, opts);

    // Only pairs with distance <= 0.3 should be stored
    // (0,1)=0.25, (1,2)=0.25, (2,3)=0.25 are <= 0.3
    // (0,2)=0.50, (0,3)=0.75, (1,3)=0.50 are > 0.3
    auto& entries = storage.Entries();
    EXPECT_EQ(entries.size(), 3);
}

TEST(PDistTest, ProgressCallback) {
    MockComparison comparison(10);
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
    pdist(comparison, storage, opts);
    EXPECT_GT(callback_count, 0);
}

TEST(PDistTest, UsesBulkComparisonWhenAvailable) {
    BulkPDistComparison comparison;
    DenseStorage storage(2);
    size_t callback_count = 0;
    PDistOptions opts;
    opts.progress = [&](size_t completed, size_t total) {
        EXPECT_EQ(completed, 1);
        EXPECT_EQ(total, 1);
        callback_count++;
    };

    pdist(comparison, storage, opts);

    EXPECT_TRUE(comparison.used_bulk);
    EXPECT_FALSE(comparison.used_compare);
    EXPECT_EQ(callback_count, 1);
    EXPECT_DOUBLE_EQ(storage.Get(0, 1), 0.125);
}
