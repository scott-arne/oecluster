#include <gtest/gtest.h>
#include <cmath>
#include "oecluster/PDist.h"
#include "oecluster/Error.h"
#include "oecluster/StorageBackend.h"
#include "oecluster/PairwiseComparison.h"
#include "oecluster/comparisons/FingerprintComparison.h"
#include <oechem.h>
#include <vector>

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

// A backend sized differently from the comparison maps (i, j) through its own
// sample count: at n=4 into DenseStorage(3), pair (0, 3) overwrites pair
// (1, 2)'s slot and pair (2, 3) writes index 3 of a 3-element vector. Only a
// debug-only assert stood between that and a release-build heap write.
TEST(PDistTest, RejectsStorageSizeMismatch) {
    MockComparison comparison(4);

    DenseStorage too_small(3);
    EXPECT_THROW(pdist(comparison, too_small), ComparisonError);

    DenseStorage too_large(5);
    EXPECT_THROW(pdist(comparison, too_large), ComparisonError);
}

// A sparse backend keeps only values at or below its cutoff. For a similarity
// that discards exactly the closest pairs, so the driver refuses the pairing
// outright. Fingerprint comparisons also take the bulk TryPDist path, which is
// why the refusal has to land ahead of it.
class PDistCutoffOrientationTest : public ::testing::Test {
protected:
    void SetUp() override {
        for (const char* smiles : {"c1ccccc1", "Cc1ccccc1", "c1ccncc1"}) {
            graph_mols_.emplace_back();
            ASSERT_TRUE(OEChem::OESmilesToMol(graph_mols_.back(), smiles)) << smiles;
        }
        for (auto& mol : graph_mols_) {
            mols_.push_back(&static_cast<OEChem::OEMolBase&>(mol));
        }
    }

    FingerprintComparison MakeComparison(bool similarity) {
        FingerprintOptions opts;
        opts.similarity = similarity;
        return FingerprintComparison(mols_, opts);
    }

    std::vector<OEChem::OEGraphMol> graph_mols_;
    std::vector<OEChem::OEMolBase*> mols_;
};

TEST_F(PDistCutoffOrientationTest, SimilarityIntoSparseStorageThrows) {
    FingerprintComparison comparison = MakeComparison(true);
    ASSERT_EQ(comparison.Facts().is_distance, Capability::No);
    SparseStorage storage(comparison.Size(), 0.5);
    size_t progress_calls = 0;
    PDistOptions opts;
    opts.progress = [&](size_t, size_t) { ++progress_calls; };

    EXPECT_THROW(pdist(comparison, storage, opts), ComparisonError);
    EXPECT_EQ(progress_calls, 0u);
}

TEST_F(PDistCutoffOrientationTest, SimilarityIntoZeroCutoffSparseStorageThrows) {
    // A zero cutoff still filters: SparseStorage drops every value above it.
    FingerprintComparison comparison = MakeComparison(true);
    SparseStorage storage(comparison.Size(), 0.0);
    EXPECT_THROW(pdist(comparison, storage), ComparisonError);
}

TEST_F(PDistCutoffOrientationTest, DistanceIntoSparseStorageStillWorks) {
    FingerprintComparison comparison = MakeComparison(false);
    ASSERT_EQ(comparison.Facts().is_distance, Capability::Yes);
    const double cutoff = 0.9;
    SparseStorage storage(comparison.Size(), cutoff);

    pdist(comparison, storage);

    size_t expected = 0;
    for (size_t i = 0; i < comparison.Size(); ++i) {
        for (size_t j = i + 1; j < comparison.Size(); ++j) {
            const double distance = comparison.Compare(i, j);
            if (distance <= cutoff) {
                ++expected;
                EXPECT_NEAR(storage.Get(i, j), distance, 1e-12);
            }
        }
    }
    EXPECT_EQ(storage.Entries().size(), expected);
}

TEST_F(PDistCutoffOrientationTest, SimilarityIntoDenseStorageStillWorks) {
    FingerprintComparison comparison = MakeComparison(true);
    DenseStorage storage(comparison.Size());

    pdist(comparison, storage);

    for (size_t i = 0; i < comparison.Size(); ++i) {
        for (size_t j = i + 1; j < comparison.Size(); ++j) {
            EXPECT_NEAR(storage.Get(i, j), comparison.Compare(i, j), 1e-12);
        }
    }
}
