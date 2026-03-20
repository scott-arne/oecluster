/**
 * @file StorageBackend.cpp
 * @brief Implementation of storage backend classes.
 */

#include "oecluster/StorageBackend.h"
#include <cassert>
#include <algorithm>

namespace OECluster {

// DenseStorage implementation

DenseStorage::DenseStorage(size_t n)
    : n_(n), data_(n * (n - 1) / 2, 0.0) {
}

void DenseStorage::Set(size_t i, size_t j, double value) {
    // Ensure i < j for upper triangle storage
    if (i > j) {
        std::swap(i, j);
    }
    assert(i < j && "Set requires i < j after swap");

    size_t index = CondensedIndex(n_, i, j);
    data_[index] = value;
}

double DenseStorage::Get(size_t i, size_t j) const {
    // Diagonal is always zero
    if (i == j) {
        return 0.0;
    }

    // Ensure i < j for upper triangle storage
    if (i > j) {
        std::swap(i, j);
    }

    size_t index = CondensedIndex(n_, i, j);
    return data_[index];
}

size_t DenseStorage::NumItems() const {
    return n_;
}

size_t DenseStorage::NumPairs() const {
    return data_.size();
}

double* DenseStorage::Data() {
    return data_.data();
}

const double* DenseStorage::Data() const {
    return data_.data();
}

size_t DenseStorage::CondensedIndex(size_t n, size_t i, size_t j) {
    // Formula: n * i - i * (i + 1) / 2 + j - i - 1
    // This matches scipy.spatial.distance.squareform condensed indexing
    assert(i < j && "CondensedIndex requires i < j");
    return n * i - i * (i + 1) / 2 + j - i - 1;
}

}  // namespace OECluster
