/**
 * @file SetDiversity.h
 * @brief The Vendi score and log-determinant diversity of a set.
 *
 * Both read the set's similarity kernel, built from distances through a
 * DiversityKernel, and take either a precomputed distance matrix or a
 * comparison evaluated lazily. The exact scores decompose the kernel with an
 * in-tree symmetric eigenvalue solver and are capped at max_exact items;
 * order-2 Vendi needs no spectrum and runs at any size.
 */

#ifndef OECLUSTER_CLUSTERING_SETDIVERSITY_H
#define OECLUSTER_CLUSTERING_SETDIVERSITY_H

#include <cstddef>
#include <limits>

#include "oecluster/PairwiseComparison.h"
#include "oecluster/StorageBackend.h"

namespace OECluster {

/** @brief How a distance d becomes a kernel entry. */
enum class DiversityKernel {
    Complement,  ///< 1 - d; distances must lie in [0, 1].
    Laplacian,   ///< exp(-d / bandwidth); distances must be non-negative.
};

/** @brief Options for vendi_score. */
struct VendiOptions {
    /// 1 (Shannon, from the spectrum) or 2 (from the Frobenius norm).
    unsigned order = 1;
    /// Distance-to-kernel transform.
    DiversityKernel kernel = DiversityKernel::Complement;
    /// Laplacian bandwidth, finite and positive; NaN (unset) for Complement.
    double bandwidth = std::numeric_limits<double>::quiet_NaN();
    /// Largest input order 1 decomposes; at least one. Order 2 ignores it.
    size_t max_exact = 2048;
    /// Worker threads for the comparison overload; 0 auto-detects hardware
    /// concurrency. An explicit value is capped at the item count.
    size_t num_threads = 0;
    /// Rows per work unit, with the same scope as num_threads; at least one.
    size_t chunk_size = 256;
};

/** @brief The result of vendi_score. */
struct VendiResult {
    /// The Vendi score, between 1 and the item count for a PSD kernel.
    double score = std::numeric_limits<double>::quiet_NaN();
    /// The order it was computed at.
    unsigned order = 1;
    /// Number of items scored.
    size_t size = 0;
    /// The kernel it was computed with.
    DiversityKernel kernel = DiversityKernel::Complement;
    /// The smallest kernel eigenvalue, on the scale of K; NaN for order 2.
    double min_eigenvalue = std::numeric_limits<double>::quiet_NaN();
    /// Sum of |lambda| / n over the eigenvalues below -tol, the mass order 1
    /// dropped; NaN for order 2.
    double negative_mass = std::numeric_limits<double>::quiet_NaN();
};

/** @brief Options for logdet_diversity. */
struct LogDetOptions {
    /// Added to every eigenvalue; finite and non-negative.
    double ridge = 0.0;
    /// Distance-to-kernel transform.
    DiversityKernel kernel = DiversityKernel::Complement;
    /// Laplacian bandwidth, finite and positive; NaN (unset) for Complement.
    double bandwidth = std::numeric_limits<double>::quiet_NaN();
    /// Largest input decomposed; at least one.
    size_t max_exact = 2048;
    /// Worker threads for the comparison overload; 0 auto-detects hardware
    /// concurrency. An explicit value is capped at the item count.
    size_t num_threads = 0;
    /// Rows per work unit, with the same scope as num_threads; at least one.
    size_t chunk_size = 256;
};

/** @brief The result of logdet_diversity. */
struct LogDetResult {
    /// log det(K + ridge I), or -infinity when K + ridge I is not numerically
    /// positive definite.
    double score = std::numeric_limits<double>::quiet_NaN();
    /// The ridge it was computed with.
    double ridge = 0.0;
    /// Number of items scored.
    size_t size = 0;
    /// The kernel it was computed with.
    DiversityKernel kernel = DiversityKernel::Complement;
    /// The smallest kernel eigenvalue, before the ridge.
    double min_eigenvalue = std::numeric_limits<double>::quiet_NaN();
    /// Number of ridged eigenvalues at or below the tolerance.
    size_t nonpositive_count = 0;
};

/**
 * @brief The Vendi score (Friedman and Dieng, TMLR 2023) of a precomputed
 * distance matrix.
 *
 * Order 1 is exp(-sum p log p) over p = lambda / n for the kernel eigenvalues
 * above tol = n * DBL_EPSILON * max|lambda|, dropping the rest without
 * renormalizing, as the reference implementation does. Order 2 is
 * n^2 / ||K||_F^2, which equals the reference's value for a PSD kernel.
 *
 * :param storage: Complete dense or memory-mapped distances.
 * :param options: Order, kernel, ceiling and threading.
 * :returns: The score and, for order 1, the spectrum's negative part.
 * :raises std::invalid_argument: On an invalid option, incomplete storage,
 *     no items, more than max_exact items at order 1, a dense kernel too
 *     large to allocate, or a distance the kernel cannot use (non-finite,
 *     outside [0, 1] for Complement, negative for Laplacian).
 * :raises std::runtime_error: If the eigenvalue solver does not converge.
 */
VendiResult vendi_score(const StorageBackend& storage,
                        const VendiOptions& options);

/**
 * @brief The Vendi score of a lazily evaluated comparison.
 *
 * Every pair is compared once, as ``Compare(i, j)`` with i < j, and never a
 * self-pair.
 *
 * :param comparison: Distances to score; cloned once per concurrent worker.
 * :param options: Order, kernel, ceiling and threading.
 * :returns: As the matrix overload.
 * :raises ComparisonError: If the comparison's facts report a similarity, a
 *     nonzero self-distance, NaN-present data or subset-scored data.
 * :raises std::invalid_argument: As the matrix overload.
 * :raises std::runtime_error: As the matrix overload.
 */
VendiResult vendi_score(PairwiseComparison& comparison,
                        const VendiOptions& options);

/**
 * @brief Positive-definite log-determinant diversity of a precomputed
 * distance matrix.
 *
 * log det(K + ridge I) from the kernel's spectrum. A ridged eigenvalue at or
 * below n * DBL_EPSILON * max|mu| scores -infinity, including for an
 * indefinite kernel whose determinant happens to be positive.
 *
 * :param storage: Complete dense or memory-mapped distances.
 * :param options: Ridge, kernel, ceiling and threading.
 * :returns: The score, the smallest eigenvalue and the nonpositive count.
 * :raises std::invalid_argument: As vendi_score's matrix overload, plus a
 *     negative or non-finite ridge.
 * :raises std::runtime_error: If the eigenvalue solver does not converge.
 */
LogDetResult logdet_diversity(const StorageBackend& storage,
                              const LogDetOptions& options);

/**
 * @brief Positive-definite log-determinant diversity of a lazily evaluated
 * comparison.
 *
 * :param comparison: Distances to score; cloned once per concurrent worker.
 * :param options: Ridge, kernel, ceiling and threading.
 * :returns: As the matrix overload.
 * :raises ComparisonError: As vendi_score's comparison overload.
 * :raises std::invalid_argument: As the matrix overload.
 * :raises std::runtime_error: As the matrix overload.
 */
LogDetResult logdet_diversity(PairwiseComparison& comparison,
                              const LogDetOptions& options);

}  // namespace OECluster

#endif  // OECLUSTER_CLUSTERING_SETDIVERSITY_H
