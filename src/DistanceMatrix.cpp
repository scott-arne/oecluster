/**
 * @file DistanceMatrix.cpp
 * @brief Distance matrix container implementation.
 */

#include "oecluster/DistanceMatrix.h"

#include <utility>

#include "oecluster/StorageBackend.h"

namespace OECluster {

DistanceMatrix::DistanceMatrix(std::unique_ptr<StorageBackend> storage,
                               std::string metric_name,
                               std::vector<std::string> labels)
    : storage_(std::move(storage)),
      metric_name_(std::move(metric_name)),
      labels_(std::move(labels)) {}

StorageBackend& DistanceMatrix::Storage() {
    return *storage_;
}

const StorageBackend& DistanceMatrix::Storage() const {
    return *storage_;
}

const std::string& DistanceMatrix::MetricName() const {
    return metric_name_;
}

const std::vector<std::string>& DistanceMatrix::Labels() const {
    return labels_;
}

size_t DistanceMatrix::NumItems() const {
    return storage_->NumItems();
}

size_t DistanceMatrix::NumPairs() const {
    return storage_->NumPairs();
}

}  // namespace OECluster
