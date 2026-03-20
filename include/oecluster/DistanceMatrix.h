/**
 * @file DistanceMatrix.h
 * @brief Distance matrix container and operations.
 */

#ifndef OECLUSTER_DISTANCEMATRIX_H
#define OECLUSTER_DISTANCEMATRIX_H

#include <memory>
#include <string>
#include <vector>

namespace OECluster {

class StorageBackend;

/**
 * @brief Container for a distance matrix with metadata.
 *
 * Wraps a StorageBackend and associates it with a metric name and
 * optional item labels.
 */
class DistanceMatrix {
public:
    /**
     * @brief Construct a DistanceMatrix.
     *
     * :param storage: Owned storage backend containing the distances.
     * :param metric_name: Name of the metric used to compute distances.
     * :param labels: Optional labels for each item.
     */
    DistanceMatrix(std::unique_ptr<StorageBackend> storage,
                   std::string metric_name,
                   std::vector<std::string> labels = {});

    /**
     * @brief Access the underlying storage backend.
     *
     * :returns: Reference to the storage backend.
     */
    StorageBackend& Storage();

    /**
     * @brief Access the underlying storage backend (const).
     *
     * :returns: Const reference to the storage backend.
     */
    const StorageBackend& Storage() const;

    /**
     * @brief Get the metric name.
     *
     * :returns: Name of the distance metric.
     */
    const std::string& MetricName() const;

    /**
     * @brief Get the item labels.
     *
     * :returns: Vector of labels (may be empty).
     */
    const std::vector<std::string>& Labels() const;

    /**
     * @brief Get the number of items.
     *
     * :returns: Number of items N.
     */
    size_t NumItems() const;

    /**
     * @brief Get the number of pairwise distances.
     *
     * :returns: Number of pairs N*(N-1)/2.
     */
    size_t NumPairs() const;

private:
    std::unique_ptr<StorageBackend> storage_;
    std::string metric_name_;
    std::vector<std::string> labels_;
};

}  // namespace OECluster

#endif  // OECLUSTER_DISTANCEMATRIX_H
