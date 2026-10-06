/**
 * @file AgglomerativeHeap.h
 * @brief The generic Lance-Williams heap algorithm behind agglomerative clustering.
 */

#ifndef OECLUSTER_SRC_CLUSTERING_AGGLOMERATIVEHEAP_H
#define OECLUSTER_SRC_CLUSTERING_AGGLOMERATIVEHEAP_H

#include "oecluster/StorageBackend.h"
#include "oecluster/clustering/Agglomerative.h"

namespace OECluster::detail {

/**
 * @brief Agglomerative clustering by the generic heap algorithm, any linkage.
 *
 * agglomerative_cluster() uses it for complete, average and weighted linkage.
 * It still accepts Single, which is how the tests compare the spanning-tree
 * path with the algorithm single linkage used through 5.19.0.
 *
 * :raises std::invalid_argument: As agglomerative_cluster().
 */
AgglomerativeResult agglomerative_heap(const StorageBackend& storage,
                                       const AgglomerativeOptions& options);

}  // namespace OECluster::detail

#endif  // OECLUSTER_SRC_CLUSTERING_AGGLOMERATIVEHEAP_H
