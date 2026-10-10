/**
 * @file agglomerative_oracle.h
 * @brief The 5.20.0 heap linkage algorithm, frozen as a differential oracle.
 *
 * Production complete, average and weighted linkage moved to the row-cache
 * kernel. This file keeps the algorithm they moved from, verbatim, so the
 * kernel can be checked against it bit for bit. It is test-only code and must
 * not be changed to track the kernel: a divergence between the two is the
 * finding the differential test exists to produce.
 */

#ifndef OECLUSTER_TESTS_AGGLOMERATIVE_ORACLE_H
#define OECLUSTER_TESTS_AGGLOMERATIVE_ORACLE_H

#include "oecluster/StorageBackend.h"
#include "oecluster/clustering/Agglomerative.h"

namespace agglomerative_oracle {

/// Cluster with the 5.20.0 heap algorithm, for every linkage including single.
OECluster::AgglomerativeResult heap_cluster(
    const OECluster::StorageBackend& storage,
    const OECluster::AgglomerativeOptions& options);

}  // namespace agglomerative_oracle

#endif  // OECLUSTER_TESTS_AGGLOMERATIVE_ORACLE_H
