/**
 * @file Error.h
 * @brief Exception hierarchy for OECluster.
 */

#ifndef OECLUSTER_ERROR_H
#define OECLUSTER_ERROR_H

#include <stdexcept>
#include <string>

namespace OECluster {

/**
 * @brief Base exception for all OECluster errors.
 */
class OEClusterError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

/**
 * @brief Exception for metric-related failures (Distance, Clone).
 */
class MetricError : public OEClusterError {
public:
    using OEClusterError::OEClusterError;
};

/**
 * @brief Exception for storage-related failures (MMap, file I/O).
 */
class StorageError : public OEClusterError {
public:
    using OEClusterError::OEClusterError;
};

}  // namespace OECluster

#endif  // OECLUSTER_ERROR_H
