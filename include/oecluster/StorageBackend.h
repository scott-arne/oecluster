/**
 * @file StorageBackend.h
 * @brief Storage backend interface for distance matrices.
 */

#ifndef OECLUSTER_STORAGEBACKEND_H
#define OECLUSTER_STORAGEBACKEND_H

#include <cstddef>
#include <vector>

namespace OECluster {

/**
 * @brief Abstract base class for storing pairwise distance matrices.
 *
 * StorageBackend defines the interface for storing and retrieving pairwise
 * distances in various formats (dense, sparse, memory-mapped, etc.).
 * All implementations use condensed distance matrix indexing where only
 * the upper triangle (i < j) is stored.
 */
class StorageBackend {
public:
    virtual ~StorageBackend() = default;

    /**
     * @brief Store a distance value for a pair of items.
     *
     * :param i: Index of first item (must be less than j).
     * :param j: Index of second item (must be greater than i).
     * :param value: Distance value to store.
     */
    virtual void Set(size_t i, size_t j, double value) = 0;

    /**
     * @brief Retrieve a distance value for a pair of items.
     *
     * Handles symmetry automatically: Get(i, j) == Get(j, i).
     * Diagonal elements (i == j) always return 0.0.
     *
     * :param i: Index of first item.
     * :param j: Index of second item.
     * :returns: Distance value between items i and j.
     */
    virtual double Get(size_t i, size_t j) const = 0;

    /**
     * @brief Get the number of items in the dataset.
     *
     * :returns: Number of items N.
     */
    virtual size_t NumItems() const = 0;

    /**
     * @brief Get the number of pairwise distances stored.
     *
     * :returns: Number of pairs, equal to N * (N - 1) / 2.
     */
    virtual size_t NumPairs() const = 0;

    /**
     * @brief Get a pointer to the raw distance data.
     *
     * For backends that store data contiguously in memory, this returns
     * a pointer to the condensed distance array. For backends that don't
     * support this (e.g., sparse), returns nullptr.
     *
     * :returns: Pointer to raw data, or nullptr if not supported.
     */
    virtual double* Data() = 0;

    /**
     * @brief Get a const pointer to the raw distance data.
     *
     * :returns: Const pointer to raw data, or nullptr if not supported.
     */
    virtual const double* Data() const = 0;

    /**
     * @brief Finalize the storage backend after all data has been written.
     *
     * Some backends (e.g., sparse) may need post-processing after all
     * distances have been stored. Default implementation is a no-op.
     */
    virtual void Finalize() {}
};

/**
 * @brief Dense storage backend using condensed distance matrix indexing.
 *
 * DenseStorage stores all N*(N-1)/2 pairwise distances in a contiguous
 * array using scipy-style condensed indexing. For items i and j where i < j,
 * the distance is stored at index: n * i - i * (i + 1) / 2 + j - i - 1.
 *
 * Example::
 *
 *     DenseStorage storage(4);  // 4 items, 6 pairs
 *     storage.Set(0, 1, 1.5);
 *     storage.Set(1, 2, 2.5);
 *     double dist = storage.Get(0, 1);  // Returns 1.5
 */
class DenseStorage : public StorageBackend {
public:
    /**
     * @brief Construct a DenseStorage for n items.
     *
     * :param n: Number of items in the dataset.
     */
    explicit DenseStorage(size_t n);

    void Set(size_t i, size_t j, double value) override;
    double Get(size_t i, size_t j) const override;
    size_t NumItems() const override;
    size_t NumPairs() const override;
    double* Data() override;
    const double* Data() const override;

private:
    /**
     * @brief Calculate the condensed index for a pair (i, j) where i < j.
     *
     * :param n: Total number of items.
     * :param i: Index of first item (must be less than j).
     * :param j: Index of second item (must be greater than i).
     * :returns: Condensed index in the range [0, n*(n-1)/2).
     */
    static size_t CondensedIndex(size_t n, size_t i, size_t j);

    size_t n_;                    ///< Number of items
    std::vector<double> data_;    ///< Condensed distance array
};

}  // namespace OECluster

#endif  // OECLUSTER_STORAGEBACKEND_H
