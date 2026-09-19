/**
 * @file KMedoids.h
 * @brief k-medoids (PAM) clustering over precomputed distances.
 */

#ifndef OECLUSTER_CLUSTERING_KMEDOIDS_H
#define OECLUSTER_CLUSTERING_KMEDOIDS_H

#include <cstddef>
#include <string>
#include <utility>
#include <vector>

#include "oecluster/StorageBackend.h"
#include "oecluster/clustering/ClusterTypes.h"

namespace OECluster {

/**
 * @brief Initial medoid selection strategy for k-medoids clustering.
 */
enum class KMedoidsInit {
    Build,           ///< Greedy PAM BUILD; best seeds, O(k N^2).
    FarthestFirst,   ///< Deterministic MaxMin; spread-out seeds, O(N^2 + k N).
    Explicit         ///< Caller-supplied medoid indices.
};

/**
 * @brief Options for k-medoids clustering.
 */
struct KMedoidsOptions {
    size_t n_clusters = 2;                        ///< Medoids to place; must be in [1, item count].
    KMedoidsInit init = KMedoidsInit::Build;      ///< Initial medoid selection strategy.
    std::vector<size_t> initial_medoids;          ///< Starting medoids; required iff init == Explicit.
    size_t max_iterations = 100;                  ///< Swap iterations before giving up.
    size_t num_threads = 0;                       ///< Worker threads; 0 auto-detects hardware concurrency.
    size_t chunk_size = 4096;                     ///< Chunk size for parallelized scans (build, swap, assignment, global medoid).
};

/**
 * @brief k-medoids clustering result with the chosen medoids and objective.
 */
class KMedoidsResult : public ClusteringResult {
public:
    KMedoidsResult() = default;
    KMedoidsResult(std::vector<ClusterLabel> labels, Clusters members,
                   std::vector<size_t> medoids, double cost,
                   size_t num_iterations, bool converged)
        : ClusteringResult(std::move(labels), std::move(members)),
          medoids_(std::move(medoids)), cost_(cost),
          num_iterations_(num_iterations), converged_(converged) {}

    /** @brief Medoid item index per cluster; entry i is the medoid of label i. */
    const std::vector<size_t>& Medoids() const { return medoids_; }
    /** @brief Sum of every item's distance to its assigned medoid. */
    double Cost() const { return cost_; }
    /** @brief Swap iterations performed. */
    size_t NumIterations() const { return num_iterations_; }
    /** @brief False when the iteration cap was reached before convergence. */
    bool Converged() const { return converged_; }

    std::string Method() const override { return "k_medoids"; }

private:
    std::vector<size_t> medoids_;
    double cost_ = 0.0;
    size_t num_iterations_ = 0;
    bool converged_ = false;
};

/**
 * @brief Cluster a complete precomputed distance matrix with k-medoids (PAM).
 *
 * Places exactly ``options.n_clusters`` medoids, each a real member of the
 * input, minimizing the sum of every item's distance to its assigned medoid.
 * Initialization is PAM BUILD, deterministic farthest-first, or caller
 * supplied; the swap phase is FastPAM1, an exact algebraic reformulation of
 * textbook PAM rather than an approximation of it.
 *
 * Four guarantees hold for every returned result:
 *
 * - Every item gets a label. ``NOISE_LABEL`` is never emitted.
 * - There are exactly ``options.n_clusters`` non-empty clusters, and
 *   ``Medoids()`` is sorted ascending with label ``i`` belonging to the
 *   ``i``-th smallest medoid index.
 * - ``SparseStorage`` is refused: a cutoff-filtered store cannot answer the
 *   arbitrary pair lookups the swap scan makes.
 * - ``Converged() == true`` means a *verified* local optimum: no single medoid
 *   swap lowers the cost the result reports, checked by recomputation rather
 *   than inferred from the incremental algebra. It is never returned without
 *   that check having run and found nothing. ``Converged() == false`` means
 *   the iteration cap was reached and no such claim is made.
 *
 * Output is byte-identical across runs, ``num_threads`` values and
 * ``chunk_size`` values; those two options change how long the call takes and
 * nothing else.
 *
 * Two preconditions are assumed rather than checked, matching every other
 * algorithm in the library. Distances must be finite -- a NaN produces
 * undefined clusters -- and non-negative, which the exactly-k guarantee leans
 * on because it makes a medoid's zero self-distance the smallest value any
 * item can see. The Python entry point's gate enforces the first by scanning
 * the stored data.
 *
 * No step of this algorithm appeals to the triangle inequality, so a
 * non-metric dissimilarity is a legitimate input.
 *
 * :param storage: Complete pairwise distance storage.
 * :param options: k-medoids clustering options.
 * :returns: Labels, clusters, medoids, cost, iteration count and convergence.
 * :raises std::invalid_argument: If options are invalid or storage is incomplete.
 * :raises std::out_of_range: If an explicit medoid index is outside the storage.
 */
KMedoidsResult k_medoids_cluster(const StorageBackend& storage,
                                 const KMedoidsOptions& options);

}  // namespace OECluster

#endif  // OECLUSTER_CLUSTERING_KMEDOIDS_H
