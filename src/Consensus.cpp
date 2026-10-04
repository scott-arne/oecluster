/**
 * @file Consensus.cpp
 * @brief Co-association consensus over an ensemble of partitions.
 */
#include "oecluster/Consensus.h"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

#if defined(_MSC_VER)
#include <intrin.h>
#endif

#include "oecluster/CondensedIndex.h"
#include "oecluster/ThreadPool.h"

namespace OECluster {

namespace {

// MSVC provides no __builtin_* bit intrinsics, so dispatch per compiler the
// way BitBirchKernels.cpp does; C++20's <bit> would replace this, but the
// target is pinned to C++17.
uint32_t popcount64(const uint64_t word) {
#if defined(_MSC_VER)
    return static_cast<uint32_t>(__popcnt64(word));
#elif defined(__clang__) || defined(__GNUC__)
    return static_cast<uint32_t>(
        __builtin_popcountll(static_cast<unsigned long long>(word)));
#else
    uint32_t count = 0;
    uint64_t remaining = word;
    while (remaining != 0u) {
        remaining &= remaining - 1u;
        ++count;
    }
    return count;
#endif
}

const double* contiguous_data(const StorageBackend& matrix,
                              const char* function_name) {
    const double* data = matrix.Data();
    if (data == nullptr) {
        throw std::invalid_argument(
            std::string(function_name) + " requires contiguous storage; a "
            "SparseStorage matrix holds only the pairs below its cutoff");
    }
    return data;
}

/// Offsets must describe a partition of `positions`, and every member must be
/// a non-empty set of in-range, distinct positions. Checked before anything
/// is written, so a malformed call cannot read out of bounds.
void check_members(size_t num_items,
                   const std::vector<size_t>& offsets,
                   const std::vector<size_t>& positions,
                   const std::vector<int>& labels) {
    if (offsets.size() < 2) {
        throw std::invalid_argument(
            "coassociation_distances requires at least one member");
    }
    if (offsets.front() != 0) {
        throw std::invalid_argument(
            "coassociation_distances: offsets must start at 0");
    }
    if (positions.size() != labels.size()) {
        throw std::invalid_argument(
            "coassociation_distances: positions and labels differ in length ("
            + std::to_string(positions.size()) + " and "
            + std::to_string(labels.size()) + ")");
    }
    // The whole offset structure is checked before any span is indexed: a
    // later decreasing entry must not be discovered only after an earlier
    // span has already read past the end of `positions`.
    for (size_t member = 0; member + 1 < offsets.size(); ++member) {
        if (offsets[member + 1] < offsets[member]) {
            throw std::invalid_argument(
                "coassociation_distances: offsets must not decrease");
        }
        if (offsets[member + 1] > positions.size()) {
            throw std::invalid_argument(
                "coassociation_distances: offset "
                + std::to_string(offsets[member + 1]) + " exceeds the "
                + std::to_string(positions.size()) + " positions given");
        }
    }
    if (offsets.back() != positions.size()) {
        throw std::invalid_argument(
            "coassociation_distances: offsets end at "
            + std::to_string(offsets.back()) + " but "
            + std::to_string(positions.size()) + " positions were given");
    }
    std::unordered_set<size_t> seen;
    for (size_t member = 0; member + 1 < offsets.size(); ++member) {
        const size_t begin = offsets[member];
        const size_t end = offsets[member + 1];
        if (end == begin) {
            throw std::invalid_argument(
                "coassociation_distances: member " + std::to_string(member)
                + " observes no item");
        }
        seen.clear();
        seen.reserve(end - begin);
        for (size_t slot = begin; slot < end; ++slot) {
            const size_t position = positions[slot];
            if (position >= num_items) {
                throw std::invalid_argument(
                    "coassociation_distances: member " + std::to_string(member)
                    + " names position " + std::to_string(position)
                    + ", out of range for " + std::to_string(num_items)
                    + " items");
            }
            if (!seen.insert(position).second) {
                throw std::invalid_argument(
                    "coassociation_distances: member " + std::to_string(member)
                    + " names position " + std::to_string(position)
                    + " more than once");
            }
        }
    }
}

}  // namespace

ConsensusMatrixSummary coassociation_distances(
    size_t num_items,
    const std::vector<size_t>& offsets,
    const std::vector<size_t>& positions,
    const std::vector<int>& labels,
    StorageBackend& destination,
    const ConsensusOptions& options) {
    if (num_items < 2) {
        throw std::invalid_argument(
            "coassociation_distances requires at least 2 items");
    }
    if (destination.NumSamples() != num_items) {
        throw std::invalid_argument(
            "coassociation_distances: destination holds "
            + std::to_string(destination.NumSamples()) + " items but "
            + std::to_string(num_items) + " were given");
    }
    double* data = destination.Data();
    if (data == nullptr) {
        throw std::invalid_argument(
            "coassociation_distances requires a contiguous destination; a "
            "SparseStorage destination would drop pairs above its cutoff");
    }
    check_members(num_items, offsets, positions, labels);
    if (options.chunk_size == 0) {
        throw std::invalid_argument(
            "coassociation_distances: chunk_size must be positive");
    }

    const size_t num_members = offsets.size() - 1;
    const size_t num_pairs = num_items * (num_items - 1) / 2;
    // MMapStorage reuses a file of the right size without clearing it, so the
    // counts must start from a known zero rather than from whatever a prior
    // run left on disk.
    std::fill(data, data + num_pairs, 0.0);

    const size_t words = (num_members + 63) / 64;
    std::vector<uint64_t> masks(num_items * words, 0u);

    ThreadPool pool(options.num_threads);
    std::vector<std::pair<int, size_t>> grouped;
    std::vector<size_t> cluster_end;

    for (size_t member = 0; member < num_members; ++member) {
        const size_t begin = offsets[member];
        const size_t end = offsets[member + 1];
        const size_t word = member / 64;
        const uint64_t bit = uint64_t{1} << (member % 64);

        grouped.clear();
        grouped.reserve(end - begin);
        for (size_t slot = begin; slot < end; ++slot) {
            masks[positions[slot] * words + word] |= bit;
            if (labels[slot] >= 0) {
                grouped.emplace_back(labels[slot], positions[slot]);
            }
        }
        if (grouped.size() < 2) {
            continue;  // nothing in this member can share a cluster
        }
        std::sort(grouped.begin(), grouped.end());

        // The unit of work is one row of one cluster, not a whole cluster:
        // a member of a few large clusters is the common case, and chunking
        // by cluster would leave all of its quadratic work on one worker.
        // Row `a` owns the pairs (a, b) for b after it in the same cluster,
        // so rows still write distinct slots.
        cluster_end.assign(grouped.size(), 0);
        size_t run_start = 0;
        for (size_t index = 1; index <= grouped.size(); ++index) {
            if (index == grouped.size()
                    || grouped[index].first != grouped[run_start].first) {
                for (size_t slot = run_start; slot < index; ++slot) {
                    cluster_end[slot] = index;
                }
                run_start = index;
            }
        }

        const size_t rows = grouped.size();
        // More workers than rows do nothing; capping before the arithmetic
        // keeps 4 * workers from wrapping for any count the size_t interface
        // admits, as take_pairs does.
        const size_t workers =
            std::max<size_t>(1, std::min(pool.NumThreads(), rows));
        // chunk_size is the ceiling, not the chunk: a member small enough to
        // fit one chunk still spreads over the pool, as take_pairs does.
        const size_t spread = (rows + 4 * workers - 1) / (4 * workers);
        const size_t rows_per_chunk =
            std::max<size_t>(1, std::min(options.chunk_size, spread));
        pool.ParallelFor(0, rows, rows_per_chunk,
            [&](size_t chunk_begin, size_t chunk_end) {
                for (size_t a = chunk_begin; a < chunk_end; ++a) {
                    const size_t last = cluster_end[a];
                    const size_t i = grouped[a].second;
                    for (size_t b = a + 1; b < last; ++b) {
                        const size_t j = grouped[b].second;
                        const size_t index = (i < j)
                            ? pair_to_condensed(i, j, num_items)
                            : pair_to_condensed(j, i, num_items);
                        data[index] += 1.0;
                    }
                }
            });
    }

    std::atomic<size_t> unobserved{0};
    // More workers than rows do nothing; capping before the arithmetic keeps
    // 4 * workers from wrapping for any count the size_t interface admits,
    // and the final min keeps a chunk near SIZE_MAX from wrapping
    // ThreadPool's own (range + chunk_size - 1) chunk count.
    const size_t finalize_workers =
        std::max<size_t>(1, std::min(pool.NumThreads(), num_items));
    // chunk_size is the ceiling, not the chunk: ThreadPool starts at most one
    // worker per chunk, so a run small enough to fit one chunk would finalize
    // single-threaded, and the shrinking upper-triangular rows stay within a
    // bounded imbalance only when each worker takes several chunks.
    const size_t finalize_spread =
        (num_items + 4 * finalize_workers - 1) / (4 * finalize_workers);
    const size_t items_per_chunk =
        std::max<size_t>(1, std::min(options.chunk_size, finalize_spread));
    pool.ParallelFor(0, num_items, items_per_chunk,
        [&](size_t chunk_begin, size_t chunk_end) {
            size_t local_unobserved = 0;
            for (size_t i = chunk_begin; i < chunk_end; ++i) {
                const uint64_t* left = masks.data() + i * words;
                for (size_t j = i + 1; j < num_items; ++j) {
                    const uint64_t* right = masks.data() + j * words;
                    uint32_t observed = 0;
                    for (size_t word = 0; word < words; ++word) {
                        observed += popcount64(left[word] & right[word]);
                    }
                    const size_t index = pair_to_condensed(i, j, num_items);
                    if (observed == 0u) {
                        data[index] = 1.0;
                        ++local_unobserved;
                    } else {
                        data[index] = 1.0 - data[index]
                            / static_cast<double>(observed);
                    }
                }
            }
            unobserved.fetch_add(local_unobserved, std::memory_order_relaxed);
        });

    destination.Finalize();

    ConsensusMatrixSummary summary;
    summary.num_partitions = num_members;
    summary.unobserved_pairs = unobserved.load(std::memory_order_relaxed);
    return summary;
}

std::vector<int> consensus_components(const StorageBackend& matrix,
                                      double threshold,
                                      const ConsensusOptions& options) {
    (void)options;
    const double* data = contiguous_data(matrix, "consensus_components");
    const size_t num_items = matrix.NumSamples();
    if (num_items < 2) {
        throw std::invalid_argument(
            "consensus_components requires at least 2 items");
    }
    if (!std::isfinite(threshold) || threshold < 0.0 || threshold > 1.0) {
        throw std::invalid_argument(
            "consensus_components: threshold must be a finite fraction in "
            "[0, 1]");
    }

    std::vector<size_t> parent(num_items);
    std::iota(parent.begin(), parent.end(), size_t{0});
    std::vector<size_t> rank(num_items, 0);

    // Iterative find with path halving: a recursive find would risk the stack
    // on a chain of a million items.
    auto find = [&parent](size_t item) {
        while (parent[item] != item) {
            parent[item] = parent[parent[item]];
            item = parent[item];
        }
        return item;
    };

    // Both sides of the comparison lose the same rounding this way, so a pair
    // whose support equals the threshold merges.
    const double cutoff = 1.0 - threshold;
    size_t index = 0;
    for (size_t i = 0; i < num_items; ++i) {
        for (size_t j = i + 1; j < num_items; ++j, ++index) {
            if (data[index] > cutoff) {
                continue;
            }
            size_t left = find(i);
            size_t right = find(j);
            if (left == right) {
                continue;
            }
            if (rank[left] < rank[right]) {
                std::swap(left, right);
            }
            parent[right] = left;
            if (rank[left] == rank[right]) {
                ++rank[left];
            }
        }
    }

    // Numbered by smallest member position, so the labels depend only on the
    // matrix and the threshold.
    std::vector<int> labels(num_items, -1);
    std::vector<int> label_of_root(num_items, -1);
    int next_label = 0;
    for (size_t i = 0; i < num_items; ++i) {
        const size_t root = find(i);
        if (label_of_root[root] < 0) {
            label_of_root[root] = next_label++;
        }
        labels[i] = label_of_root[root];
    }
    return labels;
}

ConsensusStrength consensus_strength(const StorageBackend& matrix,
                                     const std::vector<int>& labels,
                                     const ConsensusOptions& options) {
    const double* data = contiguous_data(matrix, "consensus_strength");
    const size_t num_items = matrix.NumSamples();
    if (labels.size() != num_items) {
        throw std::invalid_argument(
            "consensus_strength: " + std::to_string(labels.size())
            + " labels for " + std::to_string(num_items) + " items");
    }
    if (options.chunk_size == 0) {
        throw std::invalid_argument(
            "consensus_strength: chunk_size must be positive");
    }

    std::vector<int> distinct;
    distinct.reserve(labels.size());
    for (const int label : labels) {
        if (label >= 0) {
            distinct.push_back(label);
        }
    }
    std::sort(distinct.begin(), distinct.end());
    distinct.erase(std::unique(distinct.begin(), distinct.end()),
                   distinct.end());

    std::vector<size_t> cluster_of_item(num_items,
                                        std::numeric_limits<size_t>::max());
    std::vector<size_t> sizes(distinct.size(), 0);
    for (size_t i = 0; i < num_items; ++i) {
        if (labels[i] < 0) {
            continue;
        }
        const size_t cluster = static_cast<size_t>(
            std::lower_bound(distinct.begin(), distinct.end(), labels[i])
            - distinct.begin());
        cluster_of_item[i] = cluster;
        ++sizes[cluster];
    }

    const double nan_value = std::numeric_limits<double>::quiet_NaN();
    std::vector<double> sums(num_items, 0.0);
    ThreadPool pool(options.num_threads);
    // More workers than rows do nothing; capping before the arithmetic keeps
    // 4 * workers from wrapping for any count the size_t interface admits,
    // and the final min keeps a chunk near SIZE_MAX from wrapping
    // ThreadPool's own (range + chunk_size - 1) chunk count.
    const size_t workers =
        std::max<size_t>(1, std::min(pool.NumThreads(), num_items));
    // chunk_size is the ceiling, not the chunk: ThreadPool starts at most one
    // worker per chunk, so a matrix small enough to fit one chunk would be
    // scanned single-threaded however many threads the caller asked for.
    const size_t spread = (num_items + 4 * workers - 1) / (4 * workers);
    const size_t items_per_chunk =
        std::max<size_t>(1, std::min(options.chunk_size, spread));
    pool.ParallelFor(0, num_items, items_per_chunk,
        [&](size_t chunk_begin, size_t chunk_end) {
            for (size_t i = chunk_begin; i < chunk_end; ++i) {
                const size_t cluster = cluster_of_item[i];
                if (cluster == std::numeric_limits<size_t>::max()) {
                    continue;
                }
                double total = 0.0;
                for (size_t j = 0; j < num_items; ++j) {
                    if (j == i || cluster_of_item[j] != cluster) {
                        continue;
                    }
                    const size_t index = (i < j)
                        ? pair_to_condensed(i, j, num_items)
                        : pair_to_condensed(j, i, num_items);
                    total += 1.0 - data[index];
                }
                sums[i] = total;
            }
        });

    ConsensusStrength strength;
    strength.item_consensus.assign(num_items, nan_value);
    strength.cluster_consensus.assign(distinct.size(), nan_value);
    std::vector<double> cluster_totals(distinct.size(), 0.0);
    for (size_t i = 0; i < num_items; ++i) {
        const size_t cluster = cluster_of_item[i];
        if (cluster == std::numeric_limits<size_t>::max()) {
            continue;
        }
        const size_t size = sizes[cluster];
        if (size > 1) {
            strength.item_consensus[i] = sums[i] / static_cast<double>(size - 1);
            cluster_totals[cluster] += sums[i];
        }
    }
    for (size_t cluster = 0; cluster < distinct.size(); ++cluster) {
        const size_t size = sizes[cluster];
        if (size > 1) {
            strength.cluster_consensus[cluster] = cluster_totals[cluster]
                / (static_cast<double>(size) * static_cast<double>(size - 1));
        }
    }
    return strength;
}

}  // namespace OECluster
