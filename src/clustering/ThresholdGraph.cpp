/**
 * @file ThresholdGraph.cpp
 * @brief Internal threshold-neighbor graph utilities for clustering algorithms.
 */

#include "ThresholdGraph.h"

#include <algorithm>
#include <functional>
#include <mutex>
#include <stdexcept>
#include <tuple>
#include <utility>

#include "oecluster/ThreadPool.h"

namespace OECluster {

namespace {

void for_each_condensed_pair(size_t begin, size_t end, size_t n,
                             const std::function<void(size_t, size_t, size_t)>& body) {
    size_t row_start = 0;
    size_t i = 0;
    while (i + 1 < n && row_start + (n - i - 1) <= begin) {
        row_start += n - i - 1;
        ++i;
    }

    for (size_t k = begin; k < end; ++k) {
        while (i + 1 < n && k >= row_start + (n - i - 1)) {
            row_start += n - i - 1;
            ++i;
        }
        const size_t j = i + 1 + (k - row_start);
        body(k, i, j);
    }
}

void sort_unique_neighbors(std::vector<size_t>& neighbors) {
    std::sort(neighbors.begin(), neighbors.end());
    neighbors.erase(std::unique(neighbors.begin(), neighbors.end()), neighbors.end());
}

}  // namespace

ThresholdNeighborGraph::ThresholdNeighborGraph(std::vector<std::vector<size_t>> neighbors)
    : neighbors_(std::move(neighbors)) {}

size_t ThresholdNeighborGraph::Size() const {
    return neighbors_.size();
}

const std::vector<size_t>& ThresholdNeighborGraph::Neighbors(size_t index) const {
    return neighbors_.at(index);
}

std::vector<size_t>& ThresholdNeighborGraph::MutableNeighbors(size_t index) {
    return neighbors_.at(index);
}

ThresholdNeighborGraph BuildThresholdNeighborGraph(
    const StorageBackend& storage,
    const ThresholdGraphOptions& options) {
    if (options.threshold < 0.0) {
        throw std::invalid_argument("ThresholdGraph threshold cannot be negative");
    }

    const size_t n = storage.NumItems();
    std::vector<std::vector<size_t>> neighbors(n);
    for (size_t i = 0; i < n; ++i) {
        neighbors[i].push_back(i);
    }
    if (n < 2) {
        return ThresholdNeighborGraph(std::move(neighbors));
    }

    if (const auto* sparse = dynamic_cast<const SparseStorage*>(&storage)) {
        if (options.threshold > sparse->Cutoff()) {
            throw std::invalid_argument(
                "SparseStorage cutoff is lower than the Butina distance threshold");
        }

        for (const auto& [i, j, value] : sparse->Entries()) {
            if (value <= options.threshold) {
                neighbors[i].push_back(j);
                neighbors[j].push_back(i);
            }
        }
    } else {
        const double* data = storage.Data();
        if (data == nullptr) {
            throw std::invalid_argument(
                "Threshold graph requires contiguous or sparse storage");
        }

        std::vector<std::mutex> row_mutexes(n);
        ThreadPool pool(options.num_threads);
        const size_t chunk_size = options.chunk_size == 0 ? 4096 : options.chunk_size;
        pool.ParallelFor(0, storage.NumPairs(), chunk_size,
                         [&](size_t begin, size_t end) {
                             for_each_condensed_pair(
                                 begin, end, n,
                                 [&](size_t k, size_t i, size_t j) {
                                     if (data[k] <= options.threshold) {
                                         std::scoped_lock lock(row_mutexes[i],
                                                               row_mutexes[j]);
                                         neighbors[i].push_back(j);
                                         neighbors[j].push_back(i);
                                     }
                                 });
                         });
    }

    for (auto& row : neighbors) {
        sort_unique_neighbors(row);
    }

    return ThresholdNeighborGraph(std::move(neighbors));
}

}  // namespace OECluster
