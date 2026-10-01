/**
 * @file LeidenEngine.h
 * @brief The serial Leiden optimizer behind leiden, phase by phase.
 *
 * Every phase is a separate function so that tests can drive it on a
 * hand-built level. The functions are pure apart from the LeidenRng they
 * take by reference.
 */

#ifndef OECLUSTER_SRC_CLUSTERING_LEIDENENGINE_H
#define OECLUSTER_SRC_CLUSTERING_LEIDENENGINE_H

#include <cstddef>
#include <cstdint>
#include <random>
#include <vector>

#include "oecluster/clustering/ClusterTypes.h"
#include "oecluster/clustering/Leiden.h"

#include "SNNWeights.h"

namespace OECluster::detail {

/**
 * @brief One level of the optimization: a symmetric CSR graph without
 * self-loops in its rows, plus per-node self-loop weight, strength and item
 * count.
 *
 * strength[v] is the sum of v's row weights plus 2 * self[v]; m is the sum of
 * the row weights with each edge counted once, plus the sum of self.
 */
struct LeidenLevel {
    size_t num_nodes = 0;
    std::vector<size_t> offsets;      // num_nodes + 1 entries
    std::vector<uint32_t> neighbors;  // ascending within each node
    std::vector<double> weights;      // aligned with neighbors
    std::vector<double> self;
    std::vector<double> strength;
    std::vector<size_t> size;
    double m = 0.0;
};

/// The base level of a weighted graph: no self-loops, one item per node. The
/// graph is taken by value and its arrays moved in, so a caller that moves
/// its graph in keeps one copy of the CSR, not two.
LeidenLevel make_base_level(WeightedGraph graph);

/**
 * @brief The engine's only source of randomness.
 *
 * std::mt19937_64's output sequence is fixed by the standard; the bounded,
 * uniform and shuffle operations are written here because the standard
 * distributions and std::shuffle differ between standard libraries.
 */
class LeidenRng {
public:
    explicit LeidenRng(uint64_t seed) : engine_(seed) {}

    uint64_t Next() { return engine_(); }
    /// Uniform integer in [0, bound) by rejection; bound >= 1.
    uint64_t Below(uint64_t bound);
    /// Uniform double in [0, 1) from the top 53 bits.
    double Uniform();
    /// Fisher-Yates from the last position down.
    void Shuffle(std::vector<uint32_t>& values);

private:
    std::mt19937_64 engine_;
};

/** @brief The optimization options the engine reads. */
struct LeidenParams {
    LeidenObjective objective = LeidenObjective::Modularity;
    double resolution = 1.0;
    double theta = 0.01;
};

/**
 * @brief Reported quality of a partition of a level.
 *
 * Modularity: sum_c [e_c / m - resolution * (K_c / 2m)^2], 0 when m == 0.
 * CPM: sum_c [e_c - resolution * N_c (N_c - 1) / 2]. Summed in ascending
 * community id.
 *
 * :param community: Community id per node, each below num_nodes.
 */
double leiden_quality(const LeidenLevel& level,
                      const std::vector<uint32_t>& community,
                      const LeidenParams& params);

/**
 * @brief Unnormalized gain of adding an isolated node to a community.
 *
 * :param weight_to: Weight from the node to the community.
 * :param node_strength: The node's strength.
 * :param node_size: The node's item count.
 * :param community_strength: The community's strength, without the node.
 * :param community_size: The community's item count, without the node.
 * :param m: The level's total weight.
 */
double leiden_gain(double weight_to, double node_strength, double node_size,
                   double community_strength, double community_size,
                   double m, const LeidenParams& params);

/**
 * @brief Whether a set S inside community C is well connected to C.
 *
 * :param weight_to_rest: Weight between S and C - S.
 * :param set_strength: Strength of S.
 * :param set_size: Item count of S.
 * :param community_strength: Strength of C.
 * :param community_size: Item count of C.
 * :param m: The level's total weight.
 */
bool leiden_well_connected(double weight_to_rest, double set_strength,
                           double set_size, double community_strength,
                           double community_size, double m,
                           const LeidenParams& params);

/**
 * @brief Draw one refinement candidate.
 *
 * Candidates with a negative gain are never chosen; the others are drawn
 * with probability proportional to exp((gain - max gain) / theta), the
 * first whose cumulative weight exceeds u times the total.
 *
 * :param gains: Gain per candidate.
 * :param u: A uniform double in [0, 1).
 * :returns: The chosen index, or gains.size() if no gain is >= 0.
 */
size_t leiden_select_candidate(const std::vector<double>& gains, double theta,
                               double u);

/**
 * @brief Local moving: move nodes between communities until none gains.
 *
 * :param community: Community id per node, each below num_nodes; updated in
 *     place.
 */
void leiden_local_move(const LeidenLevel& level, std::vector<uint32_t>& community,
                       const LeidenParams& params, LeidenRng& rng);

/**
 * @brief Refinement: merge singletons within each community.
 *
 * :param community: The local-moving partition.
 * :returns: Refined community per node, named by a member's node id; every
 *     refined community lies inside one community.
 */
std::vector<uint32_t> leiden_refine(const LeidenLevel& level,
                                    const std::vector<uint32_t>& community,
                                    const LeidenParams& params, LeidenRng& rng);

/** @brief The next level and how the current level maps onto it. */
struct LeidenAggregate {
    LeidenLevel level;
    /// Aggregate node per current node.
    std::vector<uint32_t> node_of;
    /// Initial community per aggregate node.
    std::vector<uint32_t> community;
};

/**
 * @brief Aggregation: refined communities become the next level's nodes.
 *
 * :param refined: Refined community per node.
 * :param community: The local-moving partition each refined community lies
 *     inside.
 */
LeidenAggregate leiden_aggregate(const LeidenLevel& level,
                                 const std::vector<uint32_t>& refined,
                                 const std::vector<uint32_t>& community);

/**
 * @brief One full pass from a base-level partition.
 *
 * :param partition: Community id per base node, each below num_nodes.
 * :returns: Community id per base node after the pass (not renumbered).
 */
std::vector<uint32_t> leiden_pass(const LeidenLevel& base,
                                  const std::vector<uint32_t>& partition,
                                  const LeidenParams& params, LeidenRng& rng);

/** @brief The driver's output. */
struct LeidenRun {
    /// Canonical labels: clusters numbered by their smallest member.
    std::vector<ClusterLabel> labels;
    double quality = 0.0;
    size_t iterations = 0;
};

/**
 * @brief Run passes from singletons.
 *
 * :param graph: The weighted graph, moved into the base level.
 * :param n_iterations: -1 repeats passes until one leaves the canonical
 *     labels unchanged, counting that pass; otherwise exactly that many.
 */
LeidenRun run_leiden(WeightedGraph graph, const LeidenParams& params,
                     int64_t n_iterations, uint64_t seed);

}  // namespace OECluster::detail

#endif  // OECLUSTER_SRC_CLUSTERING_LEIDENENGINE_H
