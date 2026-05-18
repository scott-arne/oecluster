/**
 * @file PairwiseComparison.h
 * @brief Pairwise comparison interface.
 */

#ifndef OECLUSTER_PAIRWISECOMPARISON_H
#define OECLUSTER_PAIRWISECOMPARISON_H

#include <cstddef>
#include <memory>
#include <string>

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
     * :param i: Index of first item.
     * :param j: Index of second item.
     * :returns: Comparison value between items i and j.
     */
    virtual double Compare(size_t i, size_t j) = 0;

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
};

}  // namespace OECluster

#endif  // OECLUSTER_PAIRWISECOMPARISON_H
