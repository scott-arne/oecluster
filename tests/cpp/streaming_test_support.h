/**
 * @file streaming_test_support.h
 * @brief Fixtures shared by the comparison-built threshold graph tests.
 */

#ifndef OECLUSTER_TESTS_CPP_STREAMING_TEST_SUPPORT_H
#define OECLUSTER_TESTS_CPP_STREAMING_TEST_SUPPORT_H

#include <algorithm>
#include <atomic>
#include <cstddef>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "oecluster/PairwiseComparison.h"
#include "oecluster/StorageBackend.h"

#include "../../src/clustering/ThresholdGraph.h"

namespace streaming_test {

// Every row of a graph, so two graphs compare with one EXPECT_EQ.
inline std::vector<std::vector<size_t>> Rows(
    const OECluster::ThresholdNeighborGraph& graph) {
    std::vector<std::vector<size_t>> rows;
    for (size_t i = 0; i < graph.Size(); ++i) {
        const OECluster::NeighborRange range = graph.Neighbors(i);
        rows.emplace_back(range.begin(), range.end());
    }
    return rows;
}

// The matrix the comparison path must reproduce: every pair read once
// through Compare(i, j) on a clone, as the spec defines exactness.
inline OECluster::DenseStorage CompareFilled(
    const OECluster::PairwiseComparison& comparison) {
    std::unique_ptr<OECluster::PairwiseComparison> local = comparison.Clone();
    const size_t n = comparison.Size();
    OECluster::DenseStorage storage(n);
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i + 1; j < n; ++j) {
            storage.Set(i, j, local->Compare(i, j));
        }
    }
    return storage;
}

/**
 * @brief The same distance for every pair, and no bookkeeping.
 *
 * The default-limit tests score up to 1.4e8 pairs; a shared counter such as
 * CountingComparison's would cost more than the scoring.
 */
class ConstantComparison : public OECluster::PairwiseComparison {
public:
    ConstantComparison(size_t n, double value) : n_(n), value_(value) {}

    double Compare(size_t, size_t) override { return value_; }
    std::unique_ptr<OECluster::PairwiseComparison> Clone() const override {
        return std::make_unique<ConstantComparison>(*this);
    }
    size_t Size() const override { return n_; }
    std::string ComparisonName() const override { return "constant"; }

private:
    size_t n_;
    double value_;
};

/**
 * @brief A condensed table, except that one pair answers from a script.
 *
 * Evaluation k of the scripted pair (counting from 1, across every clone)
 * returns script[k - 1], and every later evaluation repeats the last entry.
 * The builder scores each pair once per pass, so whatever the thread
 * schedule, entry 1 is what the degree pass sees and entry 2 what the fill
 * pass sees. That makes a deliberately non-repeatable comparison
 * deterministic enough to test.
 */
class ScriptedPairComparison : public OECluster::PairwiseComparison {
public:
    ScriptedPairComparison(size_t n, std::vector<double> condensed, size_t i,
                           size_t j, std::vector<double> script)
        : n_(n),
          condensed_(std::make_shared<const std::vector<double>>(
              std::move(condensed))),
          i_(i),
          j_(j),
          script_(std::make_shared<const std::vector<double>>(std::move(script))),
          evaluations_(std::make_shared<std::atomic<size_t>>(0)) {}

    double Compare(size_t i, size_t j) override {
        if (i == i_ && j == j_) {
            const size_t k = evaluations_->fetch_add(1);
            return (*script_)[std::min(k, script_->size() - 1)];
        }
        return (*condensed_)[n_ * i - i * (i + 1) / 2 + j - i - 1];
    }
    std::unique_ptr<OECluster::PairwiseComparison> Clone() const override {
        return std::make_unique<ScriptedPairComparison>(*this);
    }
    size_t Size() const override { return n_; }
    std::string ComparisonName() const override { return "scripted"; }

    size_t Evaluations() const { return evaluations_->load(); }

private:
    size_t n_;
    std::shared_ptr<const std::vector<double>> condensed_;
    size_t i_;
    size_t j_;
    std::shared_ptr<const std::vector<double>> script_;
    std::shared_ptr<std::atomic<size_t>> evaluations_;
};

}  // namespace streaming_test

#endif  // OECLUSTER_TESTS_CPP_STREAMING_TEST_SUPPORT_H
