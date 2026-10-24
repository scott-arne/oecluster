/**
 * @file test_cdist.cpp
 * @brief Tests for cross-distance computation.
 */

#include <gtest/gtest.h>
#include <atomic>
#include <memory>
#include <string>
#include "oecluster/oecluster.h"
#include "oecluster/comparisons/FingerprintComparison.h"
#include <oechem.h>
#include <vector>

using namespace OECluster;

class BulkCDistComparison : public PairwiseComparison {
public:
    double Compare(size_t, size_t) override {
        used_compare = true;
        return 1.0;
    }

    bool TryCDist(size_t n_a, double* output, const CDistOptions& options) override {
        if (n_a != 1) {
            return false;
        }
        output[0] = 0.25;
        if (options.progress) {
            options.progress(1, 1);
        }
        used_bulk = true;
        return true;
    }

    std::unique_ptr<PairwiseComparison> Clone() const override {
        return std::make_unique<BulkCDistComparison>();
    }

    size_t Size() const override { return 2; }
    std::string ComparisonName() const override { return "bulk"; }

    bool used_bulk = false;
    bool used_compare = false;
};

class CDistTest : public ::testing::Test {
protected:
    void SetUp() override {
        auto add_mol = [this](const char* smi) {
            graph_mols_.emplace_back();
            OEChem::OESmilesToMol(graph_mols_.back(), smi);
        };
        // Set A: 2 molecules
        add_mol("c1ccccc1");      // benzene
        add_mol("c1ccc(O)cc1");   // phenol
        // Set B: 3 molecules
        add_mol("CCCCCCCC");      // octane
        add_mol("c1ccncc1");      // pyridine
        add_mol("CC(=O)O");       // acetic acid

        for (auto& gm : graph_mols_) {
            mols_.push_back(&static_cast<OEChem::OEMolBase&>(gm));
        }
    }

    std::vector<OEChem::OEGraphMol> graph_mols_;
    std::vector<OEChem::OEMolBase*> mols_;
};

TEST_F(CDistTest, OutputShape) {
    FingerprintComparison comparison(mols_);
    size_t n_a = 2;
    size_t n_b = 3;
    std::vector<double> output(n_a * n_b, -1.0);
    cdist(comparison, n_a, output.data());
    // All values should be filled (not -1.0)
    for (size_t k = 0; k < n_a * n_b; ++k) {
        EXPECT_GE(output[k], 0.0);
        EXPECT_LE(output[k], 1.0);
    }
}

TEST_F(CDistTest, MatchesDirectDistance) {
    FingerprintComparison comparison(mols_);
    size_t n_a = 2;
    size_t n_b = 3;
    std::vector<double> output(n_a * n_b, 0.0);
    cdist(comparison, n_a, output.data());

    // Verify output[i * n_b + j] == comparison.Compare(i, n_a + j)
    for (size_t i = 0; i < n_a; ++i) {
        for (size_t j = 0; j < n_b; ++j) {
            double expected = comparison.Compare(i, n_a + j);
            EXPECT_DOUBLE_EQ(output[i * n_b + j], expected);
        }
    }
}

TEST_F(CDistTest, WithCutoff) {
    FingerprintComparison comparison(mols_);
    size_t n_a = 2;
    size_t n_b = 3;
    std::vector<double> output(n_a * n_b, -1.0);
    CDistOptions opts;
    opts.cutoff = 0.3;
    cdist(comparison, n_a, output.data(), opts);
    for (size_t k = 0; k < n_a * n_b; ++k) {
        EXPECT_TRUE(output[k] == 0.0 || output[k] <= 0.3);
    }
}

TEST_F(CDistTest, ProgressCallback) {
    FingerprintComparison comparison(mols_);
    size_t n_a = 2;
    size_t n_b = 3;
    std::vector<double> output(n_a * n_b, 0.0);
    CDistOptions opts;
    opts.num_threads = 2;
    opts.chunk_size = 1;
    size_t last_completed = 0;
    size_t last_total = 0;
    opts.progress = [&](size_t completed, size_t total) {
        EXPECT_GE(completed, last_completed);
        last_completed = completed;
        last_total = total;
    };
    cdist(comparison, n_a, output.data(), opts);
    EXPECT_EQ(last_total, n_a * n_b);
}

TEST(CDistBulkTest, UsesBulkComparisonWhenAvailable) {
    BulkCDistComparison comparison;
    std::vector<double> output(1, 0.0);
    size_t callback_count = 0;
    CDistOptions opts;
    opts.progress = [&](size_t completed, size_t total) {
        EXPECT_EQ(completed, 1);
        EXPECT_EQ(total, 1);
        callback_count++;
    };

    cdist(comparison, 1, output.data(), opts);

    EXPECT_TRUE(comparison.used_bulk);
    EXPECT_FALSE(comparison.used_compare);
    EXPECT_EQ(callback_count, 1);
    EXPECT_DOUBLE_EQ(output[0], 0.25);
}

// A positive cutoff zeroes every value above it, which for a similarity is
// exactly the closest pairs, so the driver refuses the pairing. Fingerprint
// comparisons take the bulk TryCDist path, so the refusal must precede it.
TEST_F(CDistTest, SimilarityWithCutoffThrows) {
    FingerprintOptions fp_opts;
    fp_opts.similarity = true;
    FingerprintComparison comparison(mols_, fp_opts);
    ASSERT_EQ(comparison.Facts().is_distance, Capability::No);
    std::vector<double> output(2 * 3, -1.0);
    size_t progress_calls = 0;
    CDistOptions opts;
    opts.cutoff = 0.5;
    opts.progress = [&](size_t, size_t) { ++progress_calls; };

    EXPECT_THROW(cdist(comparison, 2, output.data(), opts), ComparisonError);
    EXPECT_EQ(progress_calls, 0u);
    for (double value : output) {
        EXPECT_EQ(value, -1.0);
    }
}

TEST_F(CDistTest, SimilarityWithCutoffThrowsEvenWithNoPairs) {
    // The contract must not depend on how many pairs there are.
    FingerprintOptions fp_opts;
    fp_opts.similarity = true;
    FingerprintComparison comparison(mols_, fp_opts);
    CDistOptions opts;
    opts.cutoff = 0.5;
    EXPECT_THROW(cdist(comparison, mols_.size(), nullptr, opts), ComparisonError);
}

TEST_F(CDistTest, DistanceWithCutoffStillWorks) {
    FingerprintComparison comparison(mols_);
    ASSERT_EQ(comparison.Facts().is_distance, Capability::Yes);
    const size_t n_a = 2;
    const size_t n_b = 3;
    const double cutoff = 0.9;
    std::vector<double> output(n_a * n_b, -1.0);
    CDistOptions opts;
    opts.cutoff = cutoff;

    cdist(comparison, n_a, output.data(), opts);

    for (size_t i = 0; i < n_a; ++i) {
        for (size_t j = 0; j < n_b; ++j) {
            const double distance = comparison.Compare(i, n_a + j);
            EXPECT_NEAR(output[i * n_b + j], distance > cutoff ? 0.0 : distance, 1e-12);
        }
    }
}

TEST_F(CDistTest, SimilarityWithoutCutoffStillWorks) {
    FingerprintOptions fp_opts;
    fp_opts.similarity = true;
    FingerprintComparison comparison(mols_, fp_opts);
    const size_t n_a = 2;
    const size_t n_b = 3;
    std::vector<double> output(n_a * n_b, -1.0);

    cdist(comparison, n_a, output.data());

    for (size_t i = 0; i < n_a; ++i) {
        for (size_t j = 0; j < n_b; ++j) {
            EXPECT_NEAR(output[i * n_b + j], comparison.Compare(i, n_a + j), 1e-12);
        }
    }
}

/// Distances |i - j| / n over four items, whose clones share one Clone()
/// counter.
class CloneCountingCDistComparison : public PairwiseComparison {
public:
    CloneCountingCDistComparison() : clones_(std::make_shared<std::atomic<size_t>>(0)) {}

    double Compare(size_t i, size_t j) override {
        return static_cast<double>(i > j ? i - j : j - i) / static_cast<double>(Size());
    }

    std::unique_ptr<PairwiseComparison> Clone() const override {
        clones_->fetch_add(1);
        return std::make_unique<CloneCountingCDistComparison>(*this);
    }

    size_t Size() const override { return 4; }
    std::string ComparisonName() const override { return "clone-counting"; }

    size_t Clones() const { return clones_->load(); }

private:
    std::shared_ptr<std::atomic<size_t>> clones_;
};

// TEST_F because CDistTest names the fixture above, and a suite cannot mix
// TEST and TEST_F; the fixture's molecules are unused here.
TEST_F(CDistTest, AOneChunkProblemBuildsOneClone) {
    CloneCountingCDistComparison comparison;  // 2x2 = 4 pairs, one chunk of 256
    std::vector<double> out(2 * 2, -1.0);
    CDistOptions opts;
    opts.num_threads = 8;
    cdist(comparison, 2, out.data(), opts);
    EXPECT_EQ(comparison.Clones(), 1u);
    EXPECT_DOUBLE_EQ(out[1], 0.75);  // pair (0, 3)
}

TEST_F(CDistTest, AnAbsurdThreadCountBuildsAtMostOneClonePerChunk) {
    CloneCountingCDistComparison comparison;
    std::vector<double> out(2 * 2, -1.0);
    CDistOptions opts;
    opts.num_threads = size_t{1} << 61;
    opts.chunk_size = 2;  // 4 pairs, 2 chunks
    EXPECT_NO_THROW(cdist(comparison, 2, out.data(), opts));
    EXPECT_LE(comparison.Clones(), 2u);
    EXPECT_DOUBLE_EQ(out[1], 0.75);  // pair (0, 3)
}
