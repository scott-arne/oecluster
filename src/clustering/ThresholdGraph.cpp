/**
 * @file ThresholdGraph.cpp
 * @brief Internal threshold-neighbor graph utilities for clustering algorithms.
 */

#include "ThresholdGraph.h"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <functional>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>

#include "ChunkedComparisons.h"
#include "oecluster/CondensedIndex.h"
#include "oecluster/Error.h"
#include "oecluster/ThreadPool.h"

namespace OECluster {

namespace {

// Maps condensed linear index k back to (i,j) pair by computing row start offsets.
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

constexpr size_t DEFAULT_CHUNK_SIZE = 4096;
constexpr size_t ONE_GIB = size_t{1} << 30;

// Visits the pairs at condensed offsets [begin, end), finding the first by
// condensed_to_pair's binary search. for_each_condensed_pair walks rows up
// from item 0 to find it, O(N) per chunk, which at the item counts the
// comparison build exists for would outweigh the comparisons themselves.
template <typename Body>
void for_each_pair_from(size_t begin, size_t end, size_t n, Body&& body) {
    size_t i = 0;
    size_t j = 0;
    condensed_to_pair(begin, n, i, j);
    for (size_t k = begin; k < end; ++k) {
        body(i, j);
        if (++j == n) {
            ++i;
            j = i + 1;
        }
    }
}

// Graph sizes refuse rather than wrap, as ClusterReport.cpp's pair_count
// does: a wrapped size would read as a small graph and pass the guard.
size_t checked_add(size_t a, size_t b) {
    if (a > std::numeric_limits<size_t>::max() - b) {
        throw std::length_error(
            "threshold graph size exceeds the range of size_t");
    }
    return a + b;
}

size_t checked_multiply(size_t a, size_t b) {
    if (b != 0 && a > std::numeric_limits<size_t>::max() / b) {
        throw std::length_error(
            "threshold graph size exceeds the range of size_t");
    }
    return a * b;
}

// Exactly one of n and n - 1 is even, so halving that one first keeps the
// product exact; the guard fires only when the count itself does not fit.
size_t checked_pair_count(size_t n) {
    if (n < 2) {
        return 0;
    }
    return n % 2 == 0 ? checked_multiply(n / 2, n - 1)
                      : checked_multiply(n, (n - 1) / 2);
}

std::string limit_message(const ThresholdGraphOptions& options, size_t n,
                          size_t edges, size_t bytes, size_t limit) {
    return std::string(options.caller) + " would build a threshold graph of " +
           std::to_string(bytes) + " bytes for " + std::to_string(n) +
           " items and " + std::to_string(edges) + " edges, above its " +
           (options.max_graph_bytes == 0 ? "default limit"
                                         : "max_graph_bytes limit") +
           " of " + std::to_string(limit) +
           " bytes; use a tighter threshold, a larger max_graph_bytes, or a "
           "memory-mapped matrix from pdist(output=...)";
}

// Refuses before anything is allocated: on Linux an oversized allocation
// tends to succeed, and the process is killed later, while filling it.
void enforce_graph_limit(const ThresholdGraphOptions& options, size_t n,
                         size_t edges) {
    const size_t bytes = detail::threshold_graph_bytes(n, edges);
    const size_t limit = detail::threshold_graph_limit(n, options.max_graph_bytes);
    if (bytes > limit) {
        throw std::length_error(limit_message(options, n, edges, bytes, limit));
    }
}

// ROCS fails the repeatability precondition: its overlay keeps state between
// calls, so a clone's score for a pair depends on what that clone scored
// before (test_comparison_repeatability.cpp, ROCSDependsOnItsCloneHistory),
// and the two passes could disagree. Refused by name, the remedy the design
// gives for a family that fails.
void refuse_unrepeatable(const PairwiseComparison& comparison,
                         const char* caller) {
    if (comparison.ComparisonName() == "rocs") {
        throw ComparisonError(
            std::string(caller) +
            " cannot build a threshold graph from a ROCS comparison: a ROCS "
            "score depends on what its overlay scored before, so the graph's "
            "two passes can disagree; cluster a matrix from pdist() instead");
    }
}

std::logic_error pass_mismatch_error(const char* caller, size_t row,
                                     const char* change, size_t counted,
                                     const std::string& written) {
    return std::logic_error(
        std::string(caller) + ": item " + std::to_string(row) + " " + change +
        " a neighbor between the two threshold graph passes (counted " +
        std::to_string(counted) + ", written " + written +
        "); Compare must return the same value for a pair on every call");
}

}  // namespace

namespace detail {

size_t threshold_graph_bytes(size_t n, size_t edges) {
    const size_t entries = checked_add(
        checked_add(checked_multiply(2, n), checked_multiply(2, edges)), 1);
    return checked_multiply(sizeof(size_t), entries);
}

size_t default_threshold_graph_limit(size_t n) {
    return std::max(checked_multiply(sizeof(double), checked_pair_count(n)),
                    ONE_GIB);
}

size_t threshold_graph_limit(size_t n, size_t max_graph_bytes) {
    return max_graph_bytes != 0 ? max_graph_bytes
                                : default_threshold_graph_limit(n);
}

}  // namespace detail

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

    const size_t n = storage.NumSamples();
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

        // Two-pass build: count neighbors in parallel, then atomically write indices
        // at pre-computed offsets to avoid races.
        std::vector<std::atomic<size_t>> counts(n);
        for (size_t i = 0; i < n; ++i) {
            counts[i].store(1, std::memory_order_relaxed);
        }

        ThreadPool pool(options.num_threads);
        // Clamped to the pair count, as KMedoidsSwapKernel.h and
        // ChunkedComparisons clamp theirs: ThreadPool's ceiling division wraps
        // to zero chunks for a chunk size near SIZE_MAX and would skip every
        // pair. A chunk at least as wide as the range is one chunk either way.
        const size_t chunk_size = std::min(
            options.chunk_size == 0 ? size_t{4096} : options.chunk_size,
            storage.NumPairs());
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

ThresholdNeighborGraph BuildThresholdNeighborGraph(
    PairwiseComparison& comparison,
    const ThresholdGraphOptions& options) {
    if (options.threshold < 0.0) {
        throw std::invalid_argument("ThresholdGraph threshold cannot be negative");
    }
    refuse_unrepeatable(comparison, options.caller);

    const size_t n = comparison.Size();
    if (n < 2) {
        // A trivial graph still answers to an explicit budget, so the
        // refusal holds on every input.
        enforce_graph_limit(options, n, 0);
        std::vector<size_t> offsets(n + 1, 0);
        std::vector<size_t> indices;
        if (n == 1) {
            offsets[1] = 1;
            indices.push_back(0);
        }
        return ThresholdNeighborGraph(std::move(offsets), std::move(indices));
    }

    const size_t pairs = checked_pair_count(n);
    // Normalized as the storage overload normalizes it: ChunkedComparisons
    // would hand a zero chunk to ParallelFor, which does no work, and every
    // item would come back alone.
    const size_t chunk_size = std::min(
        options.chunk_size == 0 ? DEFAULT_CHUNK_SIZE : options.chunk_size, pairs);
    detail::ChunkedComparisons work(comparison, n, options.num_threads,
                                    chunk_size);
    const std::string caller = options.caller;
    const double threshold = options.threshold;
    // Checked on both passes: a NaN would read as "not a neighbor" and -inf
    // as a neighbor, silently, where the matrix path refuses both.
    const auto within = [&caller, threshold](PairwiseComparison& local,
                                             size_t i, size_t j) {
        const double distance = local.Compare(i, j);
        if (!std::isfinite(distance)) {
            throw detail::non_finite_distance_error(caller, i, j);
        }
        return distance <= threshold;
    };

    std::vector<std::atomic<size_t>> counts(n);
    for (size_t i = 0; i < n; ++i) {
        counts[i].store(1, std::memory_order_relaxed);
    }
    work.Run(pairs, [&](PairwiseComparison& local, size_t begin, size_t end) {
        for_each_pair_from(begin, end, n, [&](size_t i, size_t j) {
            if (within(local, i, j)) {
                counts[i].fetch_add(1, std::memory_order_relaxed);
                counts[j].fetch_add(1, std::memory_order_relaxed);
            }
        });
    });

    size_t total = 0;
    for (size_t i = 0; i < n; ++i) {
        total = checked_add(total, counts[i].load(std::memory_order_relaxed));
    }
    // The counts are exact, so the graph's size is known before any of it,
    // offsets included, is allocated.
    enforce_graph_limit(options, n, (total - n) / 2);

    std::vector<size_t> offsets(n + 1, 0);
    for (size_t i = 0; i < n; ++i) {
        offsets[i + 1] = offsets[i] + counts[i].load(std::memory_order_relaxed);
    }
    std::vector<size_t> indices(total);
    std::vector<std::atomic<size_t>> positions(n);
    for (size_t i = 0; i < n; ++i) {
        indices[offsets[i]] = i;
        positions[i].store(offsets[i] + 1, std::memory_order_relaxed);
    }
    // A slot beyond the row's counted end belongs to the next row, so the
    // claim is refused before anything is written there.
    const auto claim = [&](size_t row) {
        const size_t position =
            positions[row].fetch_add(1, std::memory_order_relaxed);
        if (position >= offsets[row + 1]) {
            throw pass_mismatch_error(
                options.caller, row, "gained",
                offsets[row + 1] - offsets[row],
                "at least " + std::to_string(offsets[row + 1] - offsets[row] + 1));
        }
        return position;
    };
    work.Run(pairs, [&](PairwiseComparison& local, size_t begin, size_t end) {
        for_each_pair_from(begin, end, n, [&](size_t i, size_t j) {
            if (within(local, i, j)) {
                const size_t i_pos = claim(i);
                const size_t j_pos = claim(j);
                indices[i_pos] = j;
                indices[j_pos] = i;
            }
        });
    });

    // The index array starts zeroed, so a slot the fill pass left unwritten
    // would sort into place as a false neighbor "item 0".
    for (size_t i = 0; i < n; ++i) {
        const size_t counted = offsets[i + 1] - offsets[i];
        const size_t written =
            positions[i].load(std::memory_order_relaxed) - offsets[i];
        if (written != counted) {
            throw pass_mismatch_error(options.caller, i, "lost", counted,
                                      std::to_string(written));
        }
    }

    for (size_t i = 0; i < n; ++i) {
        std::sort(indices.begin() + static_cast<std::ptrdiff_t>(offsets[i]),
                  indices.begin() + static_cast<std::ptrdiff_t>(offsets[i + 1]));
    }
    return ThresholdNeighborGraph(std::move(offsets), std::move(indices));
}

}  // namespace OECluster
