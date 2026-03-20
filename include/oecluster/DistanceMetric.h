/**
 * @file DistanceMetric.h
 * @brief Distance metric interface.
 */

#ifndef OECLUSTER_DISTANCEMETRIC_H
#define OECLUSTER_DISTANCEMETRIC_H

#include <cstddef>
#include <memory>
#include <string>

namespace OECluster {

/**
 * @brief Abstract base class for distance metrics.
 *
 * Defines the interface for computing pairwise distances between items.
 * Implementations must be cloneable so that each thread in a parallel
 * computation can own its own scorer objects, avoiding synchronization.
 */
class DistanceMetric {
public:
    virtual ~DistanceMetric() = default;

    /**
     * @brief Compute distance between items i and j.
     *
     * Must be thread-safe when called on a Clone()'d instance.
     *
     * :param i: Index of first item.
     * :param j: Index of second item.
     * :returns: Distance value between items i and j.
     */
    virtual double Distance(size_t i, size_t j) = 0;

    /**
     * @brief Create a thread-local copy.
     *
     * Each clone owns its own scorer objects so Distance() can be
     * called without synchronization.
     *
     * :returns: A new instance that is independent of this one.
     */
    virtual std::unique_ptr<DistanceMetric> Clone() const = 0;

    /**
     * @brief Number of items loaded into this metric.
     *
     * :returns: Number of items.
     */
    virtual size_t Size() const = 0;

    /**
     * @brief Human-readable metric name.
     *
     * :returns: Name of the metric.
     */
    virtual std::string Name() const = 0;
};

}  // namespace OECluster

#endif  // OECLUSTER_DISTANCEMETRIC_H
