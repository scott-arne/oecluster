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

#include "oecluster/PairwiseComparison.h"
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
    /// Comparison forms only: the most pair comparisons a streamed run may
    /// make. 0 uses the default limit. The storage overload ignores it; the
    /// Python entry point refuses it there rather than ignoring it.
    size_t max_comparisons = 0;
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
 * ``chunk_size`` values; those two options change runtime and memory use but
 * not the result.
 *
 * Three preconditions are assumed rather than checked, matching every other
 * algorithm in the library. Distances must be finite -- a NaN produces
 * undefined clusters -- and non-negative, which the exactly-k guarantee leans
 * on because it makes a medoid's zero self-distance the smallest value any
 * item can see. Finally, the sums that compute the objective must be
 * representable; finite per-entry distances whose accumulated cost overflows
 * to infinity silently change the chosen medoid while reporting
 * ``cost=inf, converged=True``. The Python entry point's gate enforces the
 * finiteness precondition by scanning the stored data.
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

/**
 * @brief Cluster from a comparison with k-medoids, without a stored matrix.
 *
 * Scores pairs on demand instead of reading a matrix, which trades memory for
 * repeated work: the result is bit-identical to clustering a matrix filled
 * through the same ``Compare``, but a run costs one pass over every pair per
 * BUILD step and per swap iteration, plus about ``n_clusters`` passes for the
 * verification that ends it, where the storage overload pays one pass in
 * total. Peak auxiliary memory is ``O(n_clusters * item count)``.
 *
 * ``options.max_comparisons`` bounds that work. The unavoidable part -- the
 * initializer, the medoid column cache, and one verification pass, or
 * ``(init passes + n_clusters) * N * (N - 1) + N * n_clusters`` -- is priced
 * before any pair is scored and refused if it exceeds the budget. The
 * initializer is the smaller part of that: ``KMedoidsInit::Build`` costs
 * ``n_clusters`` passes over every pair where ``KMedoidsInit::FarthestFirst``
 * costs one, making the floors ``2 * n_clusters`` passes against
 * ``n_clusters + 1``. Switching therefore saves nothing at
 * ``n_clusters == 1``, a quarter at two, and under half however large
 * ``n_clusters`` grows; the verification pass underneath is what neither can
 * shed. A refusal prices both floors for the ``n_clusters`` actually asked
 * for and says whether the cheaper one would clear the limit, rather than
 * leaving a rule of thumb to be applied. The swap loop then stops
 * between passes once the budget is spent -- before the predicted-delta scan
 * or before the verification pass, each priced for what it costs -- returning
 * ``Converged() == false`` with ``NumIterations() < max_iterations``, which is
 * how a budget stop is told apart from the iteration cap.
 *
 * **``max_comparisons`` is accurate to within an additive
 * ``O(n_clusters * N)``; it is not a hard ceiling, and a run can finish over
 * it.** What is priced is the passes, which is where a streamed run's time
 * goes, and they are priced conservatively -- the verification pass is
 * costed once in the unavoidable part and charged again when it runs. What is
 * not priced is the row bookkeeping each initializer does around its passes:
 * BUILD refreshes its nearest-distance table once per step, and
 * farthest-first folds one row per pick. That remainder is deliberately left
 * out, being a rounding term against the ``O(n_clusters * N^2)`` the limit
 * does price -- 8% of a single pass at ``N = 50, n_clusters = 4``, and 0.045%
 * at ``N = 45000, n_clusters = 20``. It cannot grow into a multiple of the
 * limit, but a caller who needs a guarantee rather than a brake should leave
 * that slack when choosing a value.
 *
 * Unlike the storage overload, a non-finite distance is refused as it is
 * scored rather than assumed absent: a streamed run has no stored matrix for
 * the Python gate to scan.
 *
 * :param comparison: Pairwise comparison over the items to cluster.
 * :param options: k-medoids clustering options.
 * :returns: Labels, clusters, medoids, cost, iteration count and convergence.
 * :raises std::invalid_argument: If options are invalid, or the unavoidable
 *     comparison count exceeds ``max_comparisons``.
 * :raises std::out_of_range: If an explicit medoid index is outside the items.
 * :raises std::runtime_error: If a scored distance is not finite.
 * :raises ComparisonError: If the comparison reports similarities, a non-zero
 *     self-distance, possibly non-finite values, per-pair feature subsets, or
 *     is not repeatable.
 */
KMedoidsResult k_medoids_cluster(PairwiseComparison& comparison,
                                 const KMedoidsOptions& options);

}  // namespace OECluster

#endif  // OECLUSTER_CLUSTERING_KMEDOIDS_H
