/**
 * @file StreamingCoreDistances.cpp
 * @brief One pass over every pair, scheduled so that no two tiles share an item.
 */

#include "StreamingCoreDistances.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <mutex>
#include <stdexcept>

#include "ChunkedComparisons.h"
#include "DistanceAccess.h"
#include "StepTeam.h"
#include "ThresholdGraph.h"
#include "oecluster/ThreadPool.h"

namespace OECluster::detail {

namespace {

constexpr size_t DEFAULT_CORE_CHUNK = 64;

// Raised before any pair is read, with the matrix path's wording.
void validate_min_samples(size_t min_samples, size_t n) {
    if (min_samples == 0) {
        throw std::invalid_argument("HDBSCAN min_samples must be at least one");
    }
    if (min_samples > n) {
        throw std::invalid_argument("HDBSCAN min_samples must be at most the item count");
    }
}

double checked_core_distance(double distance, const std::string& caller,
                             size_t a, size_t b) {
    if (!std::isfinite(distance)) {
        throw non_finite_distance_error(caller, a, b);
    }
    if (distance < 0.0) {
        throw negative_distance_error(caller, a, b);
    }
    return distance == 0.0 ? 0.0 : distance;
}

// Item i's q smallest distances as a max-heap in buffer[i * q, i * q + count).
class CoreHeaps {
public:
    // The size check runs before either member allocates.
    CoreHeaps(size_t n, size_t q)
        : q_(q), buffer_(checked_size(n, q)), counts_(n, 0) {}

    // Equal values are interchangeable, so declining a value equal to the
    // current maximum leaves the same multiset of q smallest values.
    void Offer(size_t item, double distance) {
        double* heap = buffer_.data() + item * q_;
        size_t& count = counts_[item];
        if (count < q_) {
            heap[count++] = distance;
            std::push_heap(heap, heap + count);
        } else if (distance < heap[0]) {
            std::pop_heap(heap, heap + q_);
            heap[q_ - 1] = distance;
            std::push_heap(heap, heap + q_);
        }
    }

    // Every item has met n - 1 >= q others by the end of the pass.
    double Maximum(size_t item) const { return buffer_[item * q_]; }

private:
    static size_t checked_size(size_t n, size_t q) {
        if (n > std::numeric_limits<size_t>::max() / sizeof(double) / q) {
            throw std::length_error(
                "HDBSCAN core distance heaps for " + std::to_string(n) +
                " items and min_samples " + std::to_string(q + 1) +
                " exceed the addressable size");
        }
        return n * q;
    }

    size_t q_;
    std::vector<double> buffer_;
    std::vector<size_t> counts_;
};

}  // namespace

CoreBlockSchedule core_block_schedule(size_t n, size_t participants) {
    CoreBlockSchedule schedule;
    if (n == 0) {
        return schedule;
    }
    const size_t workers = std::max<size_t>(1, participants);
    const size_t blocks = std::min(n, workers > n / 4 ? n : 4 * workers);
    schedule.bounds.resize(blocks + 1);
    for (size_t b = 0; b <= blocks; ++b) {
        // b * n cannot overflow in practice; blocks <= n and the item count of
        // anything that fits in memory is far below sqrt(SIZE_MAX).
        schedule.bounds[b] = b * n / blocks;
    }

    std::vector<std::pair<size_t, size_t>> diagonal;
    for (size_t b = 0; b < blocks; ++b) {
        diagonal.emplace_back(b, b);
    }
    schedule.rounds.push_back(std::move(diagonal));

    // The circle method over an even count, the last position fixed; a padding
    // block past the real ones stands in for an odd count and is skipped.
    const size_t even = blocks + (blocks % 2);
    const size_t rotating = even - 1;
    for (size_t round = 0; round < rotating; ++round) {
        std::vector<std::pair<size_t, size_t>> tiles;
        auto add = [&](size_t a, size_t b) {
            if (a < blocks && b < blocks) {
                tiles.emplace_back(std::min(a, b), std::max(a, b));
            }
        };
        add(round, rotating);
        for (size_t k = 1; k < even / 2; ++k) {
            add((round + k) % rotating, (round + rotating - k) % rotating);
        }
        if (!tiles.empty()) {
            schedule.rounds.push_back(std::move(tiles));
        }
    }
    return schedule;
}

CoreDistances streaming_core_distances(PairwiseComparison& comparison,
                                       size_t min_samples, size_t num_threads,
                                       const std::string& caller) {
    const size_t n = comparison.Size();
    validate_min_samples(min_samples, n);
    CoreDistances result;
    const size_t q = min_samples - 1;
    if (q == 0) {
        result.values.assign(n, 0.0);
        return result;
    }

    CoreHeaps heaps(n, q);
    result.values.assign(n, 0.0);
    const size_t participants = resolve_participants(num_threads, n);
    const CoreBlockSchedule schedule = core_block_schedule(n, participants);
    ChunkedComparisons work(comparison, n, participants, 1);
    std::mutex max_mutex;
    double max_distance = 0.0;

    for (const auto& tiles : schedule.rounds) {
        work.Run(tiles.size(), [&](PairwiseComparison& local, size_t begin, size_t end) {
            double tile_max = 0.0;
            for (size_t t = begin; t < end; ++t) {
                const auto [a, b] = tiles[t];
                for (size_t i = schedule.bounds[a]; i < schedule.bounds[a + 1]; ++i) {
                    const size_t first = (a == b) ? i + 1 : schedule.bounds[b];
                    for (size_t j = first; j < schedule.bounds[b + 1]; ++j) {
                        const double distance =
                            checked_core_distance(local.Compare(i, j), caller, i, j);
                        tile_max = std::max(tile_max, distance);
                        heaps.Offer(i, distance);
                        heaps.Offer(j, distance);
                    }
                }
            }
            std::lock_guard<std::mutex> lock(max_mutex);
            max_distance = std::max(max_distance, tile_max);
        });
    }

    for (size_t i = 0; i < n; ++i) {
        result.values[i] = heaps.Maximum(i);
    }
    result.max_distance = max_distance;
    return result;
}

CoreDistances matrix_core_distances(const StorageBackend& storage,
                                    size_t min_samples, size_t num_threads,
                                    size_t chunk_size, const std::string& caller) {
    validate_complete_distance_storage(storage, "HDBSCAN");
    const size_t n = storage.NumSamples();
    validate_min_samples(min_samples, n);

    CoreDistances result;
    result.values.assign(n, 0.0);
    if (min_samples == 1) {
        return result;
    }

    const double* data = storage.Data();
    const size_t participants = resolve_participants(num_threads, n);
    const size_t unit = work_unit(
        n, participants, chunk_size == 0 ? DEFAULT_CORE_CHUNK : chunk_size);
    std::mutex max_mutex;
    double max_distance = 0.0;
    ThreadPool pool(participants);
    pool.ParallelFor(0, n, unit, [&](size_t begin, size_t end) {
        std::vector<double> row(n);
        double chunk_max = 0.0;
        for (size_t i = begin; i < end; ++i) {
            row[0] = 0.0;
            size_t write_index = 1;
            for (size_t j = 0; j < n; ++j) {
                if (i == j) {
                    continue;
                }
                const double distance = checked_core_distance(
                    dense_distance(data, n, i, j), caller, std::min(i, j), std::max(i, j));
                chunk_max = std::max(chunk_max, distance);
                row[write_index++] = distance;
            }
            const auto nth = row.begin() + static_cast<std::ptrdiff_t>(min_samples - 1);
            std::nth_element(row.begin(), nth, row.end());
            result.values[i] = *nth;
        }
        std::lock_guard<std::mutex> lock(max_mutex);
        max_distance = std::max(max_distance, chunk_max);
    });
    result.max_distance = max_distance;
    return result;
}

}  // namespace OECluster::detail
