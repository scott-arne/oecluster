/**
 * @file SNNWeights.cpp
 * @brief Shared-nearest-neighbor Jaccard weights over a KNNGraph.
 */

#include "SNNWeights.h"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

#include "oecluster/ThreadPool.h"

#include "DiversityValidation.h"

namespace OECluster::detail {

namespace {

constexpr const char* SNN_WEIGHTS_NAME = "snn_weights";

// A row costs about k * (k + 1) merge steps; a unit of roughly this many
// steps amortizes the claim without starving threads on small graphs.
constexpr size_t SNN_STEPS_PER_UNIT = size_t{1} << 16;

size_t shared_count(const uint32_t* a, const uint32_t* b, size_t width) {
    size_t shared = 0;
    const uint32_t* a_end = a + width;
    const uint32_t* b_end = b + width;
    while (a != a_end && b != b_end) {
        if (*a < *b) {
            ++a;
        } else if (*b < *a) {
            ++b;
        } else {
            ++shared;
            ++a;
            ++b;
        }
    }
    return shared;
}

}  // namespace

WeightedGraph snn_weights(const KNNGraph& graph, double prune,
                          size_t num_threads) {
    const size_t n = graph.NumItems();
    validate_leiden_item_count(n, SNN_WEIGHTS_NAME);
    WeightedGraph result;
    result.num_nodes = n;
    result.offsets.assign(n + 1, 0);
    if (n == 0) {
        return result;
    }
    const size_t k = graph.K();
    const size_t width = k + 1;
    const std::vector<size_t>& indices = graph.Indices();

    // N+(i), row i plus i itself, sorted by index for the merge and for the
    // membership test.
    std::vector<uint32_t> hood(n * width);
    for (size_t i = 0; i < n; ++i) {
        uint32_t* row = hood.data() + i * width;
        for (size_t m = 0; m < k; ++m) {
            row[m] = static_cast<uint32_t>(indices[i * k + m]);
        }
        row[k] = static_cast<uint32_t>(i);
        std::sort(row, row + width);
    }

    // A negative slot marks an arc that is not owned or was pruned. Each
    // unit writes only its own rows' slots, so the workers share no mutable
    // state, and the slots do not depend on how rows are split.
    std::vector<double> slot_weight(n * k, -1.0);
    ThreadPool pool(capped_threads(num_threads, n));
    const size_t rows_per_unit =
        std::max<size_t>(1, SNN_STEPS_PER_UNIT / (k * width));
    pool.ParallelFor(0, n, rows_per_unit, [&](size_t begin, size_t end) {
        for (size_t i = begin; i < end; ++i) {
            const uint32_t* hood_i = hood.data() + i * width;
            for (size_t m = 0; m < k; ++m) {
                const size_t j = indices[i * k + m];
                const uint32_t* hood_j = hood.data() + j * width;
                // A mutual pair is owned by its smaller index; a one-way arc
                // by the only row that names the other item.
                const bool mutual = std::binary_search(
                    hood_j, hood_j + width, static_cast<uint32_t>(i));
                if (mutual && j < i) {
                    continue;
                }
                const size_t shared = shared_count(hood_i, hood_j, width);
                const double weight = static_cast<double>(shared) /
                                      static_cast<double>(2 * width - shared);
                if (weight >= prune) {
                    slot_weight[i * k + m] = weight;
                }
            }
        }
    });
    hood = std::vector<uint32_t>();

    std::vector<size_t>& offsets = result.offsets;
    for (size_t i = 0; i < n; ++i) {
        for (size_t m = 0; m < k; ++m) {
            if (slot_weight[i * k + m] >= 0.0) {
                ++offsets[i + 1];
                ++offsets[indices[i * k + m] + 1];
            }
        }
    }
    for (size_t i = 0; i < n; ++i) {
        offsets[i + 1] += offsets[i];
    }
    result.neighbors.resize(offsets[n]);
    result.weights.resize(offsets[n]);
    std::vector<size_t> cursor(offsets.begin(), offsets.end() - 1);
    for (size_t i = 0; i < n; ++i) {
        for (size_t m = 0; m < k; ++m) {
            const double weight = slot_weight[i * k + m];
            if (weight < 0.0) {
                continue;
            }
            const size_t j = indices[i * k + m];
            result.neighbors[cursor[i]] = static_cast<uint32_t>(j);
            result.weights[cursor[i]++] = weight;
            result.neighbors[cursor[j]] = static_cast<uint32_t>(i);
            result.weights[cursor[j]++] = weight;
        }
    }

    // Each pair is owned once, so a row names each neighbor once and sorting
    // by neighbor alone is a total order.
    std::vector<std::pair<uint32_t, double>> row;
    for (size_t i = 0; i < n; ++i) {
        row.clear();
        for (size_t e = offsets[i]; e < offsets[i + 1]; ++e) {
            row.emplace_back(result.neighbors[e], result.weights[e]);
        }
        std::sort(row.begin(), row.end(),
                  [](const auto& a, const auto& b) { return a.first < b.first; });
        for (size_t e = offsets[i]; e < offsets[i + 1]; ++e) {
            result.neighbors[e] = row[e - offsets[i]].first;
            result.weights[e] = row[e - offsets[i]].second;
        }
    }
    return result;
}

}  // namespace OECluster::detail
