/**
 * @file knn_test_support.h
 * @brief Fixtures shared by the knn_graph and jarvis_patrick tests.
 */

#ifndef OECLUSTER_TESTS_CPP_KNN_TEST_SUPPORT_H
#define OECLUSTER_TESTS_CPP_KNN_TEST_SUPPORT_H

#include <gtest/gtest.h>

#include <algorithm>
#include <cstddef>
#include <filesystem>
#include <functional>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "oecluster/StorageBackend.h"
#include "oecluster/clustering/KNNGraph.h"

namespace knn_test {

const double NaN = std::numeric_limits<double>::quiet_NaN();
const double INF = std::numeric_limits<double>::infinity();

inline double At(size_t n, const std::vector<double>& condensed, size_t a,
                 size_t b) {
    if (a == b) {
        return 0.0;
    }
    const size_t i = std::min(a, b);
    const size_t j = std::max(a, b);
    return condensed[n * i - i * (i + 1) / 2 + j - i - 1];
}

// The brute-force graph: every other item sorted by (distance, index) and the
// first k kept. It shares no code with the builder's bounded heap.
inline OECluster::KNNGraph Oracle(size_t n, const std::vector<double>& condensed,
                                  size_t k) {
    std::vector<size_t> indices;
    std::vector<double> distances;
    for (size_t i = 0; i < n; ++i) {
        std::vector<std::pair<double, size_t>> row;
        for (size_t j = 0; j < n; ++j) {
            if (j != i) {
                row.emplace_back(At(n, condensed, i, j), j);
            }
        }
        std::sort(row.begin(), row.end());
        for (size_t m = 0; m < k; ++m) {
            distances.push_back(row[m].first);
            indices.push_back(row[m].second);
        }
    }
    return OECluster::KNNGraph(n, k, std::move(indices), std::move(distances));
}

inline void ExpectSameGraph(const OECluster::KNNGraph& actual,
                            const OECluster::KNNGraph& expected) {
    EXPECT_EQ(actual.NumItems(), expected.NumItems());
    EXPECT_EQ(actual.K(), expected.K());
    EXPECT_EQ(actual.Indices(), expected.Indices());
    EXPECT_EQ(actual.Distances(), expected.Distances());
}

// Finalized sparse storage holding every pair at or within the cutoff, as
// pdist() with a cutoff writes it. SparseStorage holds a mutex and cannot be
// moved, hence the pointer.
inline std::unique_ptr<OECluster::SparseStorage> MakeSparse(
    size_t n, const std::vector<double>& condensed, double cutoff) {
    auto storage = std::make_unique<OECluster::SparseStorage>(n, cutoff);
    size_t k = 0;
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i + 1; j < n; ++j) {
            storage->Set(i, j, condensed[k++]);
        }
    }
    storage->Finalize();
    return storage;
}

// A temp-file MMapStorage holding the given condensed distances. The file is
// removed first, because MMapStorage reuses an existing file of the right
// size, and removed again on destruction so a failing test leaves nothing.
class TempMMap {
public:
    TempMMap(const std::string& name, size_t n,
             const std::vector<double>& condensed)
        : path_(std::filesystem::temp_directory_path() / name) {
        std::filesystem::remove(path_);
        storage_ = std::make_unique<OECluster::MMapStorage>(path_.string(), n);
        std::copy(condensed.begin(), condensed.end(), storage_->Data());
    }
    ~TempMMap() {
        storage_.reset();
        std::filesystem::remove(path_);
    }
    TempMMap(const TempMMap&) = delete;
    TempMMap& operator=(const TempMMap&) = delete;

    const OECluster::MMapStorage& Storage() const { return *storage_; }

private:
    std::filesystem::path path_;
    std::unique_ptr<OECluster::MMapStorage> storage_;
};

inline void ExpectInvalidArgument(const std::function<void()>& call,
                                  const std::string& message) {
    try {
        call();
        FAIL() << "expected std::invalid_argument: " << message;
    } catch (const std::invalid_argument& error) {
        EXPECT_EQ(std::string(error.what()), message);
    }
}

inline void ExpectRuntimeError(const std::function<void()>& call,
                               const std::string& message) {
    try {
        call();
        FAIL() << "expected std::runtime_error: " << message;
    } catch (const std::runtime_error& error) {
        EXPECT_EQ(std::string(error.what()), message);
    }
}

}  // namespace knn_test

#endif  // OECLUSTER_TESTS_CPP_KNN_TEST_SUPPORT_H
