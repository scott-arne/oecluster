/**
 * @file KNNGraph.cpp
 * @brief k-nearest-neighbor graph: validation and row-owned construction.
 */

#include "oecluster/clustering/KNNGraph.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "DistanceAccess.h"
#include "DiversityValidation.h"
#include "KNNGraphBuild.h"
#include "oecluster/ThreadPool.h"

namespace OECluster {

namespace {

constexpr const char* KNN_NAME = "knn_graph";

std::runtime_error knn_non_finite_error(const std::string& caller, size_t a,
                                        size_t b) {
    return std::runtime_error(caller + " read a non-finite distance between items " +
                              std::to_string(std::min(a, b)) + " and " +
                              std::to_string(std::max(a, b)));
}

// Checked before the product is formed, so no later offset can wrap.
size_t checked_graph_size(size_t num_items, size_t k) {
    if (k != 0 && num_items > std::numeric_limits<size_t>::max() / k) {
        throw std::invalid_argument("KNNGraph num_items * k overflows size_t");
    }
    return num_items * k;
}

// chunk_size counts pairwise distances, as in the threshold graph, and a row
// reads n - 1 of them. The cap at n keeps a near-SIZE_MAX chunk_size from
// overflowing ParallelFor's ceiling arithmetic.
size_t rows_per_unit(size_t n, size_t chunk_size) {
    return std::min(n, std::max<size_t>(1, chunk_size / (n - 1)));
}

// A bounded max-heap of (distance, index). The lexicographic order is total,
// so the kept set is unique whatever order candidates arrive in, which is
// what makes the result independent of threads and chunking.
class RowSelector {
public:
    explicit RowSelector(size_t k) : k_(k) { heap_.reserve(k); }

    void Offer(double distance, size_t index) {
        const Candidate candidate(distance, index);
        if (heap_.size() < k_) {
            heap_.push_back(candidate);
            std::push_heap(heap_.begin(), heap_.end());
        } else if (candidate < heap_.front()) {
            std::pop_heap(heap_.begin(), heap_.end());
            heap_.back() = candidate;
            std::push_heap(heap_.begin(), heap_.end());
        }
    }

    // Writes the kept candidates in ascending order and empties the heap.
    void Drain(size_t* indices, double* distances) {
        std::sort_heap(heap_.begin(), heap_.end());
        for (size_t m = 0; m < heap_.size(); ++m) {
            distances[m] = heap_[m].first;
            indices[m] = heap_[m].second;
        }
        heap_.clear();
    }

private:
    using Candidate = std::pair<double, size_t>;
    size_t k_;
    std::vector<Candidate> heap_;
};

KNNGraph matrix_graph(const double* data, size_t n,
                      const KNNGraphOptions& options,
                      const std::string& caller) {
    const size_t k = options.k;
    std::vector<size_t> indices(checked_graph_size(n, k));
    std::vector<double> distances(indices.size());
    ThreadPool pool(detail::capped_threads(options.num_threads, n));
    // Each unit writes only its own rows' slots, which are sized before the
    // pool starts, so the workers share no mutable state.
    pool.ParallelFor(0, n, rows_per_unit(n, options.chunk_size),
                     [&](size_t begin, size_t end) {
        RowSelector row(k);
        for (size_t i = begin; i < end; ++i) {
            for (size_t j = 0; j < n; ++j) {
                if (j == i) {
                    continue;
                }
                const double distance = detail::dense_distance(data, n, i, j);
                if (!std::isfinite(distance)) {
                    throw knn_non_finite_error(caller, i, j);
                }
                row.Offer(distance, j);
            }
            row.Drain(&indices[i * k], &distances[i * k]);
        }
    });
    return KNNGraph(n, k, std::move(indices), std::move(distances));
}

}  // namespace

KNNGraph::KNNGraph(size_t num_items, size_t k, std::vector<size_t> indices,
                   std::vector<double> distances)
    : num_items_(num_items),
      k_(k),
      indices_(std::move(indices)),
      distances_(std::move(distances)) {
    const size_t total = checked_graph_size(num_items_, k_);
    if (indices_.size() != total || distances_.size() != total) {
        throw std::invalid_argument(
            "KNNGraph needs " + std::to_string(total) +
            " indices and distances for " + std::to_string(num_items_) +
            " items at k = " + std::to_string(k_) + ", got " +
            std::to_string(indices_.size()) + " indices and " +
            std::to_string(distances_.size()) + " distances");
    }
    if (num_items_ == 0) {
        return;
    }
    detail::validate_knn_k(num_items_, k_, "KNNGraph");

    // last_row[j] == r marks j as already seen in row r; rows are sorted by
    // distance, so a repeated index need not be adjacent to its first copy.
    std::vector<size_t> last_row(num_items_, std::numeric_limits<size_t>::max());
    for (size_t r = 0; r < num_items_; ++r) {
        const std::string row = "KNNGraph row " + std::to_string(r);
        for (size_t m = 0; m < k_; ++m) {
            const size_t j = indices_[r * k_ + m];
            const double distance = distances_[r * k_ + m];
            if (j >= num_items_) {
                throw std::invalid_argument(row + " names item " + std::to_string(j) +
                                            ", outside the " +
                                            std::to_string(num_items_) + " items");
            }
            if (j == r) {
                throw std::invalid_argument(row + " contains its own item");
            }
            if (last_row[j] == r) {
                throw std::invalid_argument(row + " repeats item " + std::to_string(j));
            }
            last_row[j] = r;
            if (!std::isfinite(distance)) {
                throw std::invalid_argument(row + " has a non-finite distance");
            }
            if (m > 0) {
                const std::pair<double, size_t> previous(distances_[r * k_ + m - 1],
                                                         indices_[r * k_ + m - 1]);
                if (!(previous < std::make_pair(distance, j))) {
                    throw std::invalid_argument(
                        row + " is not ordered by ascending (distance, index)");
                }
            }
        }
    }
}

namespace detail {

KNNGraph build_knn_graph(const StorageBackend& storage,
                         const KNNGraphOptions& options,
                         const std::string& caller) {
    validate_chunk_size(options.chunk_size, caller);
    const size_t n = storage.NumSamples();
    if (n == 0) {
        return KNNGraph(0, options.k, {}, {});
    }
    validate_knn_k(n, options.k, caller);
    const double* data = storage.Data();
    if (data == nullptr) {
        throw std::invalid_argument(
            caller + " requires dense, memory-mapped or sparse storage; this "
            "storage has no data array");
    }
    return matrix_graph(data, n, options, caller);
}

KNNGraph build_knn_graph(PairwiseComparison&, const KNNGraphOptions&,
                         const std::string&) {
    throw std::logic_error("knn_graph over a comparison is not implemented yet");
}

}  // namespace detail

KNNGraph knn_graph(const StorageBackend& storage, const KNNGraphOptions& options) {
    return detail::build_knn_graph(storage, options, KNN_NAME);
}

KNNGraph knn_graph(PairwiseComparison& comparison, const KNNGraphOptions& options) {
    return detail::build_knn_graph(comparison, options, KNN_NAME);
}

}  // namespace OECluster
