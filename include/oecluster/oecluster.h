/**
 * @file oecluster.h
 * @brief Main umbrella header for the OECluster library.
 *
 * Include this header to access all OECluster functionality.
 */

#ifndef OECLUSTER_OECLUSTER_H
#define OECLUSTER_OECLUSTER_H

#define OECLUSTER_VERSION_MAJOR 3
#define OECLUSTER_VERSION_MINOR 3
#define OECLUSTER_VERSION_PATCH 0

namespace OECluster {

// Forward declarations
class PairwiseComparison;
class StorageBackend;
class DenseStorage;
class MMapStorage;
class SparseStorage;
class ThreadPool;
class DistanceMatrix;

}  // namespace OECluster

#include "oecluster/Error.h"
#include "oecluster/PairwiseComparison.h"
#include "oecluster/StorageBackend.h"
#include "oecluster/ThreadPool.h"
#include "oecluster/PDist.h"
#include "oecluster/CDist.h"
#include "oecluster/DistanceMatrix.h"

#include "oecluster/comparisons/FingerprintComparison.h"
#include "oecluster/comparisons/ROCSComparison.h"
#include "oecluster/comparisons/SuperposeComparison.h"

#include "oecluster/clustering/ClusterTypes.h"
#include "oecluster/clustering/Butina.h"
#include "oecluster/clustering/Centroid.h"
#include "oecluster/clustering/DBSCAN.h"
#include "oecluster/clustering/HDBSCAN.h"
#include "oecluster/clustering/Agglomerative.h"
#include "oecluster/clustering/BitBirch.h"

#endif  // OECLUSTER_OECLUSTER_H
