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
     * @brief The value ``pdist`` stores for the unordered pair {i, j}.
     *
     * Defaults to ``Compare(i, j)``. A comparison whose ``Compare`` is
     * directional overrides this so that the single stored triangle does not
     * depend on which index ``pdist`` happened to pass first. It is protected
     * because it is not a second public way to compare two items: ``Compare``
     * keeps its own meaning, which ``cdist`` and direct callers rely on.
     */
    virtual double CompareUnordered(size_t i, size_t j) { return Compare(i, j); }

#ifndef SWIG
    /// ``pdist`` is the one caller of the seam above, and it is a free function
    /// rather than a member, so it needs friendship to reach a protected
    /// member. The call sits inside ``pdist``'s ``ParallelFor`` lambda, whose
    /// body is access-checked in the context of the enclosing function, and it
    /// goes through a ``std::unique_ptr<PairwiseComparison>``; neither changes
    /// what this declaration has to say. No default argument here: ``PDist.h``
    /// owns the one for ``options``, and a friend declaration that is also the
    /// first declaration of a function may not supply one.
    ///
    /// Hidden from SWIG. ``swig/oecluster.i`` ``%include``s this header before
    /// it declares ``pdist`` itself, and SWIG wraps a friend declaration as a
    /// namespace-scope function, so leaving this visible would have SWIG wrap
    /// ``pdist`` here -- against parameter types that are still incomplete at
    /// this point -- and then redeclare it later.
    friend void pdist(PairwiseComparison& comparison, StorageBackend& storage,
                      const PDistOptions& options);
#endif
};

}  // namespace OECluster

#endif  // OECLUSTER_PAIRWISECOMPARISON_H
