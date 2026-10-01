/**
 * @file oecluster.h
 * @brief Main umbrella header for the OECluster library.
 *
 * Include this header to access all OECluster functionality.
 */

#ifndef OECLUSTER_OECLUSTER_H
#define OECLUSTER_OECLUSTER_H

#define OECLUSTER_VERSION_MAJOR 5
#define OECLUSTER_VERSION_MINOR 11
#define OECLUSTER_VERSION_PATCH 1

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
#include "oecluster/GateFacts.h"
#include "oecluster/PairwiseComparison.h"
#include "oecluster/CondensedIndex.h"
#include "oecluster/StorageBackend.h"
#include "oecluster/ThreadPool.h"
#include "oecluster/PDist.h"
#include "oecluster/CDist.h"
#include "oecluster/DistanceMatrix.h"
#include "oecluster/DescriptorStatistics.h"

#include "oecluster/comparisons/DescriptorComparison.h"
#include "oecluster/comparisons/FingerprintComparison.h"
#include "oecluster/comparisons/MCSComparison.h"
#include "oecluster/comparisons/RMSDComparison.h"
#include "oecluster/comparisons/ROCSComparison.h"
#include "oecluster/comparisons/SuperposeComparison.h"

#include "oecluster/clustering/ClusterTypes.h"
#include "oecluster/clustering/Butina.h"
#include "oecluster/clustering/Representative.h"
#include "oecluster/clustering/DBSCAN.h"
#include "oecluster/clustering/HDBSCAN.h"
#include "oecluster/clustering/Agglomerative.h"
#include "oecluster/clustering/BitBirch.h"
#include "oecluster/clustering/KMedoids.h"
#include "oecluster/clustering/ClusterReport.h"
#include "oecluster/clustering/PartitionAgreement.h"
#include "oecluster/clustering/SARCoherence.h"
#include "oecluster/clustering/DiversitySelection.h"
#include "oecluster/clustering/SetDiversity.h"
#include "oecluster/clustering/SphereExclusion.h"
#include "oecluster/clustering/KNNGraph.h"
#include "oecluster/clustering/JarvisPatrick.h"
#include "oecluster/clustering/Leiden.h"
#include "oecluster/clustering/MurckoScaffold.h"

#endif  // OECLUSTER_OECLUSTER_H
