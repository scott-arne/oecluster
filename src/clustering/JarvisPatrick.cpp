/**
 * @file JarvisPatrick.cpp
 * @brief Jarvis-Patrick clustering over a KNNGraph and over raw input.
 */

#include "oecluster/clustering/JarvisPatrick.h"

#include <algorithm>
#include <cstddef>
#include <iterator>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "DiversityValidation.h"
#include "KNNGraphBuild.h"

namespace OECluster {

namespace {

constexpr const char* JARVIS_PATRICK_NAME = "jarvis_patrick";

class JarvisPatrickUnionFind {
public:
    explicit JarvisPatrickUnionFind(size_t n) : parent_(n) {
        std::iota(parent_.begin(), parent_.end(), size_t{0});
    }

    size_t Find(size_t x) {
        while (parent_[x] != x) {
            parent_[x] = parent_[parent_[x]];
            x = parent_[x];
        }
        return x;
    }

    void Union(size_t a, size_t b) {
        const size_t root_a = Find(a);
        const size_t root_b = Find(b);
        if (root_a != root_b) {
            parent_[std::max(root_a, root_b)] = std::min(root_a, root_b);
        }
    }

private:
    std::vector<size_t> parent_;
};

void validate_kmin(size_t kmin, size_t k) {
    if (kmin >= k) {
        throw std::invalid_argument(
            std::string(JARVIS_PATRICK_NAME) + " kmin must be less than k = " +
            std::to_string(k) + ", got " + std::to_string(kmin) +
            "; a mutual pair shares at most k - 1 neighbors");
    }
}

// Raw-input validation in the shared order: chunk_size, zero items, k, then
// kmin, so a bad kmin costs no comparison. Returns false for zero items.
bool validate_raw(size_t n, const JarvisPatrickOptions& options) {
    detail::validate_chunk_size(options.chunk_size, JARVIS_PATRICK_NAME);
    if (n == 0) {
        return false;
    }
    detail::validate_knn_k(n, options.k, JARVIS_PATRICK_NAME);
    validate_kmin(options.kmin, options.k);
    return true;
}

KNNGraphOptions graph_options(const JarvisPatrickOptions& options) {
    KNNGraphOptions graph;
    graph.k = options.k;
    graph.num_threads = options.num_threads;
    graph.chunk_size = options.chunk_size;
    return graph;
}

}  // namespace

JarvisPatrickResult jarvis_patrick(const KNNGraph& graph, size_t kmin) {
    const size_t n = graph.NumItems();
    const size_t k = graph.K();
    if (n == 0) {
        return JarvisPatrickResult({}, {}, k, kmin);
    }
    validate_kmin(kmin, k);

    // Rows are stored by (distance, index); membership tests and the
    // intersection need them by index.
    std::vector<size_t> sorted = graph.Indices();
    for (size_t i = 0; i < n; ++i) {
        std::sort(sorted.begin() + i * k, sorted.begin() + (i + 1) * k);
    }

    JarvisPatrickUnionFind sets(n);
    for (size_t i = 0; i < n; ++i) {
        const auto row_i = sorted.begin() + i * k;
        for (size_t m = 0; m < k; ++m) {
            const size_t j = row_i[m];
            // A link needs both directions, so j < i is decided from j's row.
            if (j < i) {
                continue;
            }
            const auto row_j = sorted.begin() + j * k;
            if (!std::binary_search(row_j, row_j + k, i)) {
                continue;
            }
            size_t shared = 0;
            auto a = row_i;
            auto b = row_j;
            while (a != row_i + k && b != row_j + k) {
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
            if (shared >= kmin) {
                sets.Union(i, j);
            }
        }
    }

    // Ascending items meet each component first at its smallest member, which
    // orders clusters by smallest member and lists members ascending.
    const size_t unassigned = std::numeric_limits<size_t>::max();
    std::vector<size_t> label_of_root(n, unassigned);
    std::vector<ClusterLabel> labels(n);
    Clusters clusters;
    for (size_t i = 0; i < n; ++i) {
        const size_t root = sets.Find(i);
        if (label_of_root[root] == unassigned) {
            label_of_root[root] = clusters.size();
            clusters.emplace_back();
        }
        labels[i] = static_cast<ClusterLabel>(label_of_root[root]);
        clusters[label_of_root[root]].push_back(i);
    }
    return JarvisPatrickResult(std::move(labels), std::move(clusters), k, kmin);
}

JarvisPatrickResult jarvis_patrick(const StorageBackend& storage,
                                   const JarvisPatrickOptions& options) {
    if (!validate_raw(storage.NumSamples(), options)) {
        return JarvisPatrickResult({}, {}, options.k, options.kmin);
    }
    return jarvis_patrick(
        detail::build_knn_graph(storage, graph_options(options), JARVIS_PATRICK_NAME),
        options.kmin);
}

JarvisPatrickResult jarvis_patrick(PairwiseComparison& comparison,
                                   const JarvisPatrickOptions& options) {
    if (!validate_raw(comparison.Size(), options)) {
        return JarvisPatrickResult({}, {}, options.k, options.kmin);
    }
    return jarvis_patrick(
        detail::build_knn_graph(comparison, graph_options(options),
                                JARVIS_PATRICK_NAME),
        options.kmin);
}

}  // namespace OECluster
