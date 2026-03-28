/**
 * @file oecluster.h
 * @brief Main umbrella header for the OECluster library.
 *
 * Include this header to access all OECluster functionality.
 */

#ifndef OECLUSTER_OECLUSTER_H
#define OECLUSTER_OECLUSTER_H

#define OECLUSTER_VERSION_MAJOR 3
#define OECLUSTER_VERSION_MINOR 1
#define OECLUSTER_VERSION_PATCH 5

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
#include "oecluster/CDist.h"
#include "oecluster/DistanceMatrix.h"

#include "oecluster/metrics/FingerprintMetric.h"
#include "oecluster/metrics/ROCSMetric.h"
#include "oecluster/metrics/SuperposeMetric.h"

#endif  // OECLUSTER_OECLUSTER_H
