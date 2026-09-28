/**
 * @file diversity_test_support.h
 * @brief Fixtures shared by the maxmin_select and circles tests.
 */

#ifndef OECLUSTER_TESTS_CPP_DIVERSITY_TEST_SUPPORT_H
#define OECLUSTER_TESTS_CPP_DIVERSITY_TEST_SUPPORT_H

#include <atomic>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "oecluster/GateFacts.h"
#include "oecluster/PairwiseComparison.h"
#include "oecluster/StorageBackend.h"

namespace diversity_test {

// The condensed upper triangle of distance(i, j) for i < j, in the order
// StorageBackend::Data() lays it out.
inline std::vector<double> Condensed(
    size_t n, const std::function<double(size_t, size_t)>& distance) {
    std::vector<double> condensed;
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i + 1; j < n; ++j) {
            condensed.push_back(distance(i, j));
        }
    }
    return condensed;
}

inline OECluster::DenseStorage MakeStorage(size_t n,
                                           const std::vector<double>& condensed) {
    OECluster::DenseStorage storage(n);
    size_t k = 0;
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i + 1; j < n; ++j) {
            storage.Set(i, j, condensed[k++]);
        }
    }
    return storage;
}

// Items on a line at their own index, d(i, j) = j - i. Integer distances keep
// every comparison, sum and pick distance exact.
inline std::vector<double> Line(size_t n) {
    return Condensed(n, [](size_t i, size_t j) {
        return static_cast<double>(j - i);
    });
}

// Items at arbitrary positions on a line, d(i, j) = |x_i - x_j|.
inline std::vector<double> Positions(const std::vector<double>& x) {
    return Condensed(x.size(), [&](size_t i, size_t j) {
        return std::abs(x[i] - x[j]);
    });
}

// A deterministic non-metric fixture with plenty of ties: integer distances in
// [1, 6], so brute-force oracles exercise the smallest-index rule often.
inline std::vector<double> Scrambled(size_t n) {
    return Condensed(n, [](size_t i, size_t j) {
        return static_cast<double>((i * 7 + j * 13) % 6 + 1);
    });
}

/**
 * @brief A comparison over a fixed condensed table that enforces the lazy
 * path's reading contract.
 *
 * A reversed pair (i > j) always throws, which pins the Compare(min, max)
 * rule. A self-pair throws unless a self-distance is configured, in which case
 * it returns that value; a large one would change the Farthest seed if the
 * diagonal were ever read. Clones share a counter so a test can bound them.
 */
class TableComparison : public OECluster::PairwiseComparison {
public:
    TableComparison(size_t n, std::vector<double> condensed,
                    OECluster::GateFacts facts = OECluster::GateFacts(),
                    double self_distance = std::numeric_limits<double>::quiet_NaN())
        : n_(n),
          condensed_(std::make_shared<const std::vector<double>>(std::move(condensed))),
          facts_(facts),
          self_distance_(self_distance),
          clones_(std::make_shared<std::atomic<size_t>>(0)) {}

    double Compare(size_t i, size_t j) override {
        if (i == j) {
            if (std::isnan(self_distance_)) {
                throw std::logic_error("TableComparison read the diagonal at " +
                                       std::to_string(i));
            }
            return self_distance_;
        }
        if (i > j) {
            throw std::logic_error("TableComparison read the reversed pair (" +
                                   std::to_string(i) + ", " +
                                   std::to_string(j) + ")");
        }
        return (*condensed_)[n_ * i - i * (i + 1) / 2 + j - i - 1];
    }

    OECluster::GateFacts Facts() const override { return facts_; }

    std::unique_ptr<OECluster::PairwiseComparison> Clone() const override {
        clones_->fetch_add(1);
        return std::make_unique<TableComparison>(*this);
    }

    size_t Size() const override { return n_; }
    std::string ComparisonName() const override { return "table"; }

    size_t NumClones() const { return clones_->load(); }

private:
    size_t n_;
    std::shared_ptr<const std::vector<double>> condensed_;
    OECluster::GateFacts facts_;
    double self_distance_;
    std::shared_ptr<std::atomic<size_t>> clones_;
};

// Storage whose Data() is null while pairs exist, which every diversity entry
// point must refuse rather than dereference.
class NullDataStorage : public OECluster::StorageBackend {
public:
    explicit NullDataStorage(size_t n) : n_(n) {}

    void Set(size_t, size_t, double) override {}
    double Get(size_t, size_t) const override { return 0.0; }
    size_t NumSamples() const override { return n_; }
    size_t NumPairs() const override { return n_ * (n_ - 1) / 2; }
    double* Data() override { return nullptr; }
    const double* Data() const override { return nullptr; }

private:
    size_t n_;
};

}  // namespace diversity_test

#endif  // OECLUSTER_TESTS_CPP_DIVERSITY_TEST_SUPPORT_H
