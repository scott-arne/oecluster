/**
 * @file Subset.h
 * @brief Row-subset gathers over distance storage and fingerprint batches.
 */
#ifndef OECLUSTER_SUBSET_H
#define OECLUSTER_SUBSET_H

#include <cstddef>
#include <vector>

#include "oecluster/StorageBackend.h"
#include "oefp/batch.h"

namespace OECluster {

/**
 * @brief Copy the pairwise distances among selected items into a smaller storage.
 *
 * ``destination`` is caller-allocated with ``NumSamples() == indices.size()``
 * (the pdist arrangement, so no owned pointer crosses the Python boundary)
 * and receives ``source.Get(indices[a], indices[b])`` at ``(a, b)`` for every
 * ``a < b``, so the subset follows the given order. A dense or memory-mapped
 * source needs a ``DenseStorage`` destination and is gathered in parallel
 * over destination rows; a sparse source needs a ``SparseStorage``
 * destination with the same cutoff and is gathered in one pass over its
 * merged entries (so the source must have been finalized), copying only
 * stored pairs. A memory-mapped destination is refused: two mappings of one
 * file would alias the source unseen by the identity check. ``Finalize()``
 * is called on the destination in every case.
 *
 * :param source: Storage holding the full matrix.
 * :param indices: Distinct positions in ``[0, source.NumSamples())``, in the
 *     order the subset should take.
 * :param destination: Storage for the subset, sized to ``indices.size()``.
 * :param num_threads: Worker threads for the dense gather; 0 auto-detects.
 * :param chunk_size: Upper bound on destination rows per work unit; the
 *     chunk used is ``max(1, min(chunk_size, ceil(m / (4 * workers))))``.
 * :raises std::invalid_argument: For empty indices, an index out of range, a
 *     repeated index, a destination of another size, a zero chunk, a
 *     destination aliasing the source, a memory-mapped destination, a
 *     storage-kind mismatch, unequal sparse cutoffs, or a sparse destination
 *     that already holds entries (finalized or not).
 */
void take_pairs(const StorageBackend& source,
                const std::vector<size_t>& indices,
                StorageBackend& destination,
                size_t num_threads = 0,
                size_t chunk_size = 4096);

/**
 * @brief Build a batch holding the selected rows of another, in the given order.
 *
 * :param batch: Source batch.
 * :param indices: Distinct row positions in ``[0, batch.Size())``.
 * :returns: A new batch with the source's fingerprint spec and the selected
 *     rows.
 * :raises std::invalid_argument: For empty indices, an index out of range or
 *     a repeated index.
 */
OEFP::OEFPBatch take_fingerprints(const OEFP::OEFPBatch& batch,
                                  const std::vector<size_t>& indices);

}  // namespace OECluster

#endif  // OECLUSTER_SUBSET_H
