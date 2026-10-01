/**
 * @file Leiden.h
 * @brief Leiden community detection over a shared-nearest-neighbor graph.
 */

#ifndef OECLUSTER_CLUSTERING_LEIDEN_H
#define OECLUSTER_CLUSTERING_LEIDEN_H

#include <cstddef>
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

#include "oecluster/PairwiseComparison.h"
#include "oecluster/StorageBackend.h"
#include "oecluster/clustering/ClusterTypes.h"
#include "oecluster/clustering/KNNGraph.h"

namespace OECluster {

/** @brief Quality function optimized by leiden. */
enum class LeidenObjective {
    Modularity,  ///< Reichardt-Bornholdt modularity with a resolution.
    CPM,         ///< Constant Potts Model.
};

/** @brief Options for leiden. */
struct LeidenOptions {
    /// Neighbors per item for the graph (raw overloads only); 1 <= k <= n - 1.
    size_t k = 0;
    LeidenObjective objective = LeidenObjective::Modularity;
    /// Resolution; finite and >= 0.
    double resolution = 1.0;
    /// SNN Jaccard weights below this are dropped; 0 <= prune < 1.
    double prune = 1.0 / 15.0;
    /// Refinement randomness; finite and > 0.
    double theta = 0.01;
    /// -1 iterates until a pass changes nothing; otherwise >= 0 passes.
    int64_t n_iterations = -1;
    /// Seed for the random stream.
    uint64_t seed = 0;
    /// Worker threads for the graph and the SNN weights; 0 selects the
    /// hardware concurrency.
    size_t num_threads = 0;
    /// Pairwise distances per work unit for the graph; at least one.
    size_t chunk_size = 4096;
};

/**
 * @brief Leiden result: labels, members, the partition's quality, and the
 * options that shaped it.
 *
 * Clusters are ordered by their smallest member and list members in
 * ascending index. Every item is assigned; an item with no retained edge is
 * a singleton.
 */
class LeidenResult : public ClusteringResult {
public:
    LeidenResult() = default;
    LeidenResult(std::vector<ClusterLabel> labels, Clusters clusters,
                 double quality, size_t iterations, LeidenObjective objective,
                 double resolution, size_t k)
        : ClusteringResult(std::move(labels), std::move(clusters)),
          quality_(quality),
          iterations_(iterations),
          objective_(objective),
          resolution_(resolution),
          k_(k) {}

    /// Objective value of the returned partition on the SNN graph.
    double Quality() const { return quality_; }
    /// Full passes run, including a final pass that changed nothing.
    size_t Iterations() const { return iterations_; }
    LeidenObjective Objective() const { return objective_; }
    double Resolution() const { return resolution_; }
    /// Neighbors per item in the graph; 0 on a default-constructed result.
    size_t K() const { return k_; }
    std::string Method() const override { return "leiden"; }

private:
    double quality_ = 0.0;
    size_t iterations_ = 0;
    LeidenObjective objective_ = LeidenObjective::Modularity;
    double resolution_ = 0.0;
    size_t k_ = 0;
};

}  // namespace OECluster

#endif  // OECLUSTER_CLUSTERING_LEIDEN_H
