/**
 * @file KNNGraph.h
 * @brief k-nearest-neighbor graph over distance matrices and comparisons.
 */

#ifndef OECLUSTER_CLUSTERING_KNNGRAPH_H
#define OECLUSTER_CLUSTERING_KNNGRAPH_H

#include <cstddef>
#include <vector>

#include "oecluster/PairwiseComparison.h"
#include "oecluster/StorageBackend.h"

namespace OECluster {

/** @brief Options for knn_graph. */
struct KNNGraphOptions {
    /// Neighbors per item; 1 <= k <= n - 1 for n >= 1 items.
    size_t k = 0;
    /// Worker threads; 0 selects the hardware concurrency.
    size_t num_threads = 0;
    /// Pairwise distances per work unit; at least one.
    size_t chunk_size = 4096;
};

/**
 * @brief Each item's k nearest other items, with their raw distances.
 *
 * Indices() and Distances() are row-major with NumItems() * K() entries. Row
 * i occupies [i * K(), (i + 1) * K()), never names item i, and is ordered by
 * ascending (distance, index), so equal distances go to the lower index. The
 * values are distances, not affinities: a larger value is a weaker tie.
 */
class KNNGraph {
public:
    /// An empty graph: zero items and K() == 0.
    KNNGraph() = default;

    /**
     * @brief Adopt row-major neighbor arrays after validating them.
     *
     * :param num_items: Number of items (rows).
     * :param k: Neighbors per row; stored as given when num_items is zero.
     * :param indices: num_items * k neighbor indices.
     * :param distances: num_items * k distances, aligned with indices.
     * :raises std::invalid_argument: If num_items * k overflows size_t, an
     *     array has the wrong length, k is outside 1..num_items - 1 for
     *     num_items >= 1, or a row names an out-of-range item, its own item
     *     or one item twice, holds a non-finite distance, or is not ordered
     *     by ascending (distance, index).
     */
    KNNGraph(size_t num_items, size_t k, std::vector<size_t> indices,
             std::vector<double> distances);

    size_t NumItems() const { return num_items_; }
    size_t K() const { return k_; }
    const std::vector<size_t>& Indices() const { return indices_; }
    const std::vector<double>& Distances() const { return distances_; }

private:
    size_t num_items_ = 0;
    size_t k_ = 0;
    std::vector<size_t> indices_;
    std::vector<double> distances_;
};

/**
 * @brief The k-nearest-neighbor graph of a precomputed distance matrix.
 *
 * Each row is an independent bounded selection over every other item, so
 * the result does not depend on num_threads or chunk_size. Sparse storage
 * must hold every pair at or within its cutoff, as pdist() with a cutoff
 * writes it: the builder cannot verify that, and a missing nearer pair would
 * silently change a row. Duplicate sparse entries count once, with the value
 * Get() reports.
 *
 * :param storage: Dense, memory-mapped or finalized sparse storage.
 * :param options: k and threading options.
 * :returns: The graph; empty, keeping options.k, for zero items.
 * :raises std::invalid_argument: On a zero chunk_size, k outside
 *     1..n - 1 (so one item is refused), data-less non-sparse storage, or a
 *     sparse item with fewer than k distinct stored neighbors.
 * :raises std::runtime_error: If a distance read is NaN or infinite (sparse:
 *     stored entries only).
 */
KNNGraph knn_graph(const StorageBackend& storage, const KNNGraphOptions& options);

/**
 * @brief The k-nearest-neighbor graph over a comparison, evaluated lazily.
 *
 * Row i evaluates Compare(min(i, j), max(i, j)) for every j != i, N(N-1)
 * calls in total, and keeps O(N * k) results plus one comparison clone per
 * running work unit. Compare must return the same value on every call and
 * every clone; every built-in comparison does.
 *
 * :param comparison: Distance comparison; cloned once per running unit.
 * :param options: k and threading options.
 * :returns: The graph; empty, keeping options.k, for zero items.
 * :raises std::invalid_argument: On a zero chunk_size or k outside 1..n - 1.
 * :raises ComparisonError: If the comparison's facts rule out ranking its
 *     distances.
 * :raises std::runtime_error: If a comparison returns NaN or infinity.
 */
KNNGraph knn_graph(PairwiseComparison& comparison, const KNNGraphOptions& options);

}  // namespace OECluster

#endif  // OECLUSTER_CLUSTERING_KNNGRAPH_H
