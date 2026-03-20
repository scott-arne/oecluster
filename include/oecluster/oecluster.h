/**
 * @file oecluster.h
 * @brief Main umbrella header for the OECluster library.
 *
 * Include this header to access all OECluster functionality.
 */

#ifndef OECLUSTER_OECLUSTER_H
#define OECLUSTER_OECLUSTER_H

#define OECLUSTER_VERSION_MAJOR 0
#define OECLUSTER_VERSION_MINOR 1
#define OECLUSTER_VERSION_PATCH 0

namespace OECluster {

// Forward declarations
class DistanceMetric;
class StorageBackend;
class DenseStorage;
class MMapStorage;
class SparseStorage;
class ThreadPool;
class DistanceMatrix;

}  // namespace OECluster

#include "oecluster/Error.h"
#include "oecluster/DistanceMetric.h"
#include "oecluster/StorageBackend.h"
#include "oecluster/ThreadPool.h"
#include "oecluster/PDist.h"
#include "oecluster/DistanceMatrix.h"

#endif  // OECLUSTER_OECLUSTER_H
