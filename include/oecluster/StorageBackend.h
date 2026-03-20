/**
 * @file StorageBackend.h
 * @brief Storage backend interface for distance matrices.
 */

#ifndef OECLUSTER_STORAGEBACKEND_H
#define OECLUSTER_STORAGEBACKEND_H

#include <cstddef>
#include <string>
#include <vector>
#include <unordered_map>
#include <tuple>
#include <thread>
#include <shared_mutex>

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

/**
 * @brief Memory-mapped storage backend using condensed distance matrix indexing.
 *
 * MMapStorage uses the same condensed indexing as DenseStorage but backs the
 * data with a memory-mapped file for out-of-core operation. The file persists
 * after the object is destroyed, allowing data to be reused across sessions.
 *
 * Example::
 *
 *     MMapStorage storage("distances.bin", 4);  // 4 items, 6 pairs
 *     storage.Set(0, 1, 1.5);
 *     storage.Set(1, 2, 2.5);
 *     double dist = storage.Get(0, 1);  // Returns 1.5
 */
class MMapStorage : public StorageBackend {
public:
    /**
     * @brief Construct a memory-mapped storage backend.
     *
     * Creates or opens a file at path and maps it into memory.
     * If the file exists and has the correct size, it is reused.
     * Otherwise, a new file is created and zero-filled.
     *
     * :param path: File path for the memory-mapped file.
     * :param n: Number of items.
     * :raises StorageError: If file operations fail.
     */
    MMapStorage(const std::string& path, size_t n);
    ~MMapStorage() override;

    MMapStorage(const MMapStorage&) = delete;
    MMapStorage& operator=(const MMapStorage&) = delete;

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

    size_t n_;              ///< Number of items
    size_t num_pairs_;      ///< Number of pairs
    size_t file_size_;      ///< Size of the mapped file in bytes
    std::string path_;      ///< File path
    int fd_;                ///< File descriptor
    double* data_;          ///< Pointer to mapped memory
};

/**
 * @brief Sparse storage backend that stores only distances below a cutoff.
 *
 * SparseStorage filters distances during Set() operations, storing only values
 * at or below the specified cutoff threshold. It uses per-thread buffers for
 * thread-safe concurrent writes, which are merged during Finalize() for efficient
 * lookups.
 *
 * Example::
 *
 *     SparseStorage storage(4, 0.5);  // 4 items, cutoff = 0.5
 *     storage.Set(0, 1, 0.3);  // Stored (below cutoff)
 *     storage.Set(0, 2, 0.8);  // Ignored (above cutoff)
 *     storage.Finalize();
 *     double dist = storage.Get(0, 1);  // Returns 0.3
 *     double dist2 = storage.Get(0, 2);  // Returns 0.0 (not stored)
 */
class SparseStorage : public StorageBackend {
public:
    /**
     * @brief Construct a SparseStorage for n items with cutoff threshold.
     *
     * :param n: Number of items in the dataset.
     * :param cutoff: Maximum distance value to store (inclusive).
     */
    SparseStorage(size_t n, double cutoff);
    ~SparseStorage() override = default;

    void Set(size_t i, size_t j, double value) override;
    double Get(size_t i, size_t j) const override;
    size_t NumItems() const override;
    size_t NumPairs() const override;
    double* Data() override;
    const double* Data() const override;
    void Finalize() override;

    /**
     * @brief Access merged entries after Finalize().
     *
     * :returns: Vector of tuples (i, j, value) for all stored distances.
     */
    const std::vector<std::tuple<size_t, size_t, double>>& Entries() const;

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

    size_t n_;                ///< Number of items
    double cutoff_;           ///< Maximum distance value to store

    /// Per-thread buffers for concurrent Set() calls
    std::shared_mutex buffers_mutex_;
    std::unordered_map<std::thread::id,
                       std::vector<std::tuple<size_t, size_t, double>>> thread_buffers_;

    /// Merged results after Finalize()
    std::vector<std::tuple<size_t, size_t, double>> merged_entries_;
    std::unordered_map<size_t, double> lookup_;  ///< condensed index -> value
};

}  // namespace OECluster

#endif  // OECLUSTER_STORAGEBACKEND_H
