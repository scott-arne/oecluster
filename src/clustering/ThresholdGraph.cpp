/**
 * @file ThresholdGraph.cpp
 * @brief Internal threshold-neighbor graph utilities for clustering algorithms.
 */

#include "ThresholdGraph.h"

#include <algorithm>
#include <atomic>
#include <functional>
#include <numeric>
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

ThresholdNeighborGraph make_compact_graph(std::vector<std::vector<size_t>> neighbors) {
    std::vector<size_t> offsets(neighbors.size() + 1, 0);
    size_t total = 0;
    for (size_t i = 0; i < neighbors.size(); ++i) {
        sort_unique_neighbors(neighbors[i]);
        offsets[i] = total;
        total += neighbors[i].size();
    }
    offsets[neighbors.size()] = total;

    std::vector<size_t> indices;
    indices.reserve(total);
    for (const auto& row : neighbors) {
        indices.insert(indices.end(), row.begin(), row.end());
    }

    return ThresholdNeighborGraph(std::move(offsets), std::move(indices));
}

}  // namespace

NeighborRange::NeighborRange(const size_t* begin, const size_t* end)
    : begin_(begin), end_(end) {}

const size_t* NeighborRange::begin() const {
    return begin_;
}

const size_t* NeighborRange::end() const {
    return end_;
}

size_t NeighborRange::size() const {
    return static_cast<size_t>(end_ - begin_);
}

bool NeighborRange::empty() const {
    return begin_ == end_;
}

ThresholdNeighborGraph::ThresholdNeighborGraph(std::vector<std::vector<size_t>> neighbors) {
    ThresholdNeighborGraph compact = make_compact_graph(std::move(neighbors));
    offsets_ = std::move(compact.offsets_);
    indices_ = std::move(compact.indices_);
}

ThresholdNeighborGraph::ThresholdNeighborGraph(
    std::vector<size_t> offsets,
    std::vector<size_t> indices)
    : offsets_(std::move(offsets)), indices_(std::move(indices)) {}

size_t ThresholdNeighborGraph::Size() const {
    return offsets_.empty() ? 0 : offsets_.size() - 1;
}

NeighborRange ThresholdNeighborGraph::Neighbors(size_t index) const {
    if (index >= Size()) {
        throw std::out_of_range("ThresholdNeighborGraph index is outside the graph");
    }
    const size_t begin = offsets_[index];
    const size_t end = offsets_[index + 1];
    return NeighborRange(indices_.data() + begin, indices_.data() + end);
}

ThresholdNeighborGraph BuildThresholdNeighborGraph(
    const StorageBackend& storage,
    const ThresholdGraphOptions& options) {
    if (options.threshold < 0.0) {
        throw std::invalid_argument("ThresholdGraph threshold cannot be negative");
    }

    const size_t n = storage.NumItems();
    std::vector<size_t> offsets(n + 1, 0);
    if (n < 2) {
        std::vector<size_t> indices;
        if (n == 1) {
            offsets[1] = 1;
            indices.push_back(0);
        }
        return ThresholdNeighborGraph(std::move(offsets), std::move(indices));
    }

    if (const auto* sparse = dynamic_cast<const SparseStorage*>(&storage)) {
        if (options.threshold > sparse->Cutoff()) {
            throw std::invalid_argument(
                "SparseStorage cutoff is lower than the Butina distance threshold");
        }

        std::vector<std::vector<size_t>> neighbors(n);
        for (size_t i = 0; i < n; ++i) {
            neighbors[i].push_back(i);
        }
        for (const auto& [i, j, value] : sparse->Entries()) {
            if (value <= options.threshold) {
                neighbors[i].push_back(j);
                neighbors[j].push_back(i);
            }
        }
        return make_compact_graph(std::move(neighbors));
    } else {
        const double* data = storage.Data();
        if (data == nullptr) {
            throw std::invalid_argument(
                "Threshold graph requires contiguous or sparse storage");
        }

        std::vector<std::atomic<size_t>> counts(n);
        for (size_t i = 0; i < n; ++i) {
            counts[i].store(1, std::memory_order_relaxed);
        }

        ThreadPool pool(options.num_threads);
        const size_t chunk_size = options.chunk_size == 0 ? 4096 : options.chunk_size;
        pool.ParallelFor(0, storage.NumPairs(), chunk_size,
                         [&](size_t begin, size_t end) {
                             for_each_condensed_pair(
                                 begin, end, n,
                                 [&](size_t k, size_t i, size_t j) {
                                     if (data[k] <= options.threshold) {
                                         counts[i].fetch_add(1, std::memory_order_relaxed);
                                         counts[j].fetch_add(1, std::memory_order_relaxed);
                                     }
                                 });
                         });

        for (size_t i = 0; i < n; ++i) {
            offsets[i + 1] = offsets[i] + counts[i].load(std::memory_order_relaxed);
        }

        std::vector<size_t> indices(offsets.back());
        std::vector<std::atomic<size_t>> positions(n);
        for (size_t i = 0; i < n; ++i) {
            indices[offsets[i]] = i;
            positions[i].store(offsets[i] + 1, std::memory_order_relaxed);
        }

        pool.ParallelFor(0, storage.NumPairs(), chunk_size,
                         [&](size_t begin, size_t end) {
                             for_each_condensed_pair(
                                 begin, end, n,
                                 [&](size_t k, size_t i, size_t j) {
                                     if (data[k] <= options.threshold) {
                                         const size_t i_pos =
                                             positions[i].fetch_add(1, std::memory_order_relaxed);
                                         const size_t j_pos =
                                             positions[j].fetch_add(1, std::memory_order_relaxed);
                                         indices[i_pos] = j;
                                         indices[j_pos] = i;
                                     }
                                 });
                         });

        for (size_t i = 0; i < n; ++i) {
            std::sort(indices.begin() + static_cast<std::ptrdiff_t>(offsets[i]),
                      indices.begin() + static_cast<std::ptrdiff_t>(offsets[i + 1]));
        }

        return ThresholdNeighborGraph(std::move(offsets), std::move(indices));
    }
}

}  // namespace OECluster
