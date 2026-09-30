/**
 * @file KNNGraphBuild.h
 * @brief knn_graph builders shared with jarvis_patrick, parameterized by the
 * public entry point's name so that each caller's messages name itself.
 */

#ifndef OECLUSTER_SRC_CLUSTERING_KNNGRAPHBUILD_H
#define OECLUSTER_SRC_CLUSTERING_KNNGRAPHBUILD_H

#include <cstddef>
#include <stdexcept>
#include <string>

#include "oecluster/PairwiseComparison.h"
#include "oecluster/StorageBackend.h"
#include "oecluster/clustering/KNNGraph.h"

namespace OECluster::detail {

// For n >= 1. One item has no other item to be near, so no k fits it.
inline void validate_knn_k(size_t n, size_t k, const std::string& caller) {
    if (n == 1) {
        throw std::invalid_argument(
            caller + " needs at least two items: a single item has no neighbors");
    }
    if (k == 0 || k > n - 1) {
        throw std::invalid_argument(caller + " k must be between 1 and " +
                                    std::to_string(n - 1) + " for " +
                                    std::to_string(n) + " items, got " +
                                    std::to_string(k));
    }
}

KNNGraph build_knn_graph(const StorageBackend& storage,
                         const KNNGraphOptions& options,
                         const std::string& caller);

KNNGraph build_knn_graph(PairwiseComparison& comparison,
                         const KNNGraphOptions& options,
                         const std::string& caller);

}  // namespace OECluster::detail

#endif  // OECLUSTER_SRC_CLUSTERING_KNNGRAPHBUILD_H
