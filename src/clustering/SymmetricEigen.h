/**
 * @file SymmetricEigen.h
 * @brief Eigenvalues of a dense symmetric matrix, for the set diversity
 *        scores.
 *
 * In-tree rather than a LAPACK or Eigen dependency: FetchContent cannot reach
 * either from behind the build firewall, and the scores need eigenvalues only.
 */

#ifndef OECLUSTER_SRC_CLUSTERING_SYMMETRICEIGEN_H
#define OECLUSTER_SRC_CLUSTERING_SYMMETRICEIGEN_H

#include <cstddef>
#include <vector>

namespace OECluster::detail {

/**
 * @brief The element count of a dense n x n double matrix, or a refusal.
 *
 * A caller can raise the exact-score ceiling, so the ceiling alone does not
 * rule out an n whose square, or whose square's byte count, wraps size_t.
 *
 * :param n: Matrix order.
 * :returns: n * n.
 * :raises std::invalid_argument: If n * n or n * n * sizeof(double) is not
 *     representable in size_t, or n * n exceeds vector<double>::max_size().
 */
size_t dense_kernel_elements(size_t n);

/**
 * @brief Eigenvalues of a dense symmetric matrix, in ascending order.
 *
 * Householder reduction to tridiagonal form without accumulating the
 * transformation, then implicit QL with Wilkinson shifts. Single-threaded with
 * a fixed operation order, so the same input gives bit-identical output.
 *
 * :param a: Row-major n x n matrix; consumed. Only the lower triangle,
 *     diagonal included, is read.
 * :param n: Matrix order; 0 returns an empty vector.
 * :param max_iterations: QL iterations allowed per eigenvalue. A test seam:
 *     the public entry points always pass the default.
 * :returns: The n eigenvalues, ascending.
 * :raises std::invalid_argument: If a.size() is not n * n, or n * n is not
 *     representable (see dense_kernel_elements).
 * :raises std::runtime_error: If an eigenvalue has not converged after
 *     max_iterations iterations.
 */
std::vector<double> symmetric_eigenvalues(std::vector<double> a, size_t n,
                                          unsigned max_iterations = 30);

}  // namespace OECluster::detail

#endif  // OECLUSTER_SRC_CLUSTERING_SYMMETRICEIGEN_H
