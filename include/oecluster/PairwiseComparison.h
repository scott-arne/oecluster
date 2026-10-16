/**
 * @file PairwiseComparison.h
 * @brief Pairwise comparison interface.
 */

#ifndef OECLUSTER_PAIRWISECOMPARISON_H
#define OECLUSTER_PAIRWISECOMPARISON_H

#include <cstddef>
#include <memory>
#include <string>
#include "oecluster/GateFacts.h"

namespace OECluster {

struct CDistOptions;
struct PDistOptions;
class StorageBackend;
class PairwiseComparison;

#ifndef SWIG
namespace detail {
/// Declared here so the class below can name it a friend, and documented in
/// src/CompareUnordered.h, which is where its contract is written.
///
/// Hidden from SWIG, and not merely for tidiness: measured on 2026-10-16 with
/// the guard removed, SWIG flattens the namespace and generates
/// ``oecluster.compare_unordered``, a public Python function that reaches the
/// protected seam on any comparison. The guard is what keeps the bindings
/// from offering a second way to compare two items.
double compare_unordered(PairwiseComparison& comparison, size_t i, size_t j);
}  // namespace detail
#endif

/**
 * @brief Abstract base class for pairwise comparison methods.
 *
 * Implementations own the data and scorer state needed to compare indexed
 * items. They must be cloneable so each worker thread can evaluate pairs
 * without synchronizing scorer internals.
 */
class PairwiseComparison {
public:
    virtual ~PairwiseComparison() = default;

    /**
     * @brief Compare items i and j.
     *
     * The returned scalar follows the comparison's configured output mode,
     * usually distance or similarity.
     *
     * Both indices must name items the comparison holds. Every concrete
     * comparison in this library refuses an index at or beyond ``Size()`` with
     * ``ComparisonError`` rather than reading past its container.
     *
     * :param i: Index of first item.
     * :param j: Index of second item.
     * :returns: Comparison value between items i and j.
     */
    virtual double Compare(size_t i, size_t j) = 0;

    /**
     * @brief Report the metric facts the capability gate checks.
     *
     * The base implementation reports ``Capability::Unknown`` for all three
     * capability facts -- ``is_distance``, ``zero_self`` and ``triangle`` --
     * and ``DataIntegrity::Complete``. Subclasses that can prove a
     * property should override this; the gate reads facts off the comparison
     * object rather than off the Python builder so that directly constructed
     * objects cannot bypass it.
     *
     * :returns: The facts describing this comparison's scores.
     */
    virtual GateFacts Facts() const { return GateFacts(); }

    /**
     * @brief Optionally compute all pairwise values in bulk.
     *
     * Implementations that own a more efficient batch kernel can override this
     * method. Returning ``false`` asks the generic engine to fall back to
     * per-pair ``Compare()`` calls.
     *
     * :param storage: Storage backend to write results into.
     * :param options: Pairwise computation options.
     * :returns: ``true`` when the implementation handled the computation.
     */
    virtual bool TryPDist(StorageBackend& storage, const PDistOptions& options) {
        (void)storage;
        (void)options;
        return false;
    }

    /**
     * @brief Optionally compute cross-distance values in bulk.
     *
     * Implementations that own a more efficient batch kernel can override this
     * method. Returning ``false`` asks the generic engine to fall back to
     * per-pair ``Compare()`` calls.
     *
     * :param n_a: Number of items in the first set.
     * :param output: Row-major output array of length ``n_a * (Size() - n_a)``.
     * :param options: Cross-distance computation options.
     * :returns: ``true`` when the implementation handled the computation.
     */
    virtual bool TryCDist(size_t n_a, double* output, const CDistOptions& options) {
        (void)n_a;
        (void)output;
        (void)options;
        return false;
    }

    /**
     * @brief Create a thread-local copy.
     *
     * :returns: A new comparison instance that can be used independently.
     */
    virtual std::unique_ptr<PairwiseComparison> Clone() const = 0;

    /**
     * @brief Number of items loaded into this comparison.
     *
     * :returns: Number of items.
     */
    virtual size_t Size() const = 0;

    /**
     * @brief Human-readable comparison name.
     *
     * :returns: Name of the comparison.
     */
    virtual std::string ComparisonName() const = 0;

protected:
    /**
     * @brief The value every caller reads for the unordered pair {i, j}.
     *
     * Defaults to ``Compare(i, j)``. A comparison whose ``Compare`` is
     * directional overrides this so that a reader presenting d(i, j) as
     * symmetric does not get a value that depends on which index it happened
     * to pass first. Direction is not the only way that can happen: what an
     * override owes its callers is a value that is a function of the unordered
     * pair, so a comparison carrying scorer state from one call to the next
     * has to deal with that too. ``ROCSComparison`` is both, and does both.
     *
     * ``pdist`` fills its single triangle through here, and so does every
     * other reader of a symmetric distance -- ``cluster_report``,
     * ``activity_landscape``, ``modelability``, ``maxmin_select``, ``circles``,
     * ``knn_graph``, ``sphere_exclusion`` away from the neighbor order,
     * ``vendi_score`` and ``logdet_diversity`` -- so that the lazy route and a
     * matrix from ``pdist`` answer with the same number. The paths that still
     * read ``Compare(min(i, j), max(i, j))`` instead are the ones that refuse
     * ROCS outright; ``src/clustering/DiversityValidation.h`` is where that
     * refusal and the consequence of lifting it are recorded.
     *
     * It is protected because it is not a second public way to compare two
     * items: ``Compare`` keeps its own meaning, which ``cdist`` and direct
     * callers rely on. Callers reach it through
     * ``detail::compare_unordered``, declared above, which passes the pair in
     * ascending index order -- that, and not this default, is what makes the
     * unoverridden forward to ``Compare`` a function of the unordered pair.
     */
    virtual double CompareUnordered(size_t i, size_t j) { return Compare(i, j); }

#ifndef SWIG
    /// The one friend, and the reason the accessor exists: every caller of
    /// the seam is a free function rather than a member, so without a shared
    /// accessor each would need a friend declaration of its own. ``pdist``
    /// had one until 2026-10-16 and no longer does; it goes through the
    /// accessor like everything else.
    ///
    /// Guarded for the reason given on the declaration above, and the guard
    /// has to cover this line too: the friend names ``detail::`` and will not
    /// compile where that declaration is hidden. A second, smaller reason
    /// retired with the ``friend void pdist`` this replaces -- measured on
    /// 2026-10-15, SWIG matched that declaration to the ``pdist`` it already
    /// wraps and attached the comment above it to the generated
    /// ``oecluster.pdist`` as that function's Python docstring.
    friend double detail::compare_unordered(PairwiseComparison& comparison,
                                            size_t i, size_t j);
#endif
};

}  // namespace OECluster

#endif  // OECLUSTER_PAIRWISECOMPARISON_H
