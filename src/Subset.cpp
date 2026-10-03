/**
 * @file Subset.cpp
 * @brief Row-subset gathers over distance storage and fingerprint batches.
 */
#include "oecluster/Subset.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_set>
#include <vector>

#include "oecluster/ThreadPool.h"
#include "oefp/fingerprint.h"

namespace OECluster {

namespace {

constexpr size_t ABSENT = static_cast<size_t>(-1);

/// The index rules are the same for both gathers: non-empty, in range,
/// distinct. The duplicate set is sized by the subset, not by the source, so
/// a small subset of a huge batch stays at subset-scaled cost.
void check_indices(const std::vector<size_t>& indices, size_t num_items,
                   const char* function_name) {
    if (indices.empty()) {
        throw std::invalid_argument(
            std::string(function_name) + " requires at least one index");
    }
    std::unordered_set<size_t> seen;
    seen.reserve(indices.size());
    for (size_t index : indices) {
        if (index >= num_items) {
            throw std::invalid_argument(
                std::string(function_name) + ": index " + std::to_string(index)
                + " is out of range for " + std::to_string(num_items)
                + " items");
        }
        if (!seen.insert(index).second) {
            throw std::invalid_argument(
                std::string(function_name) + ": index " + std::to_string(index)
                + " appears more than once");
        }
    }
}

void gather_dense(const StorageBackend& source,
                  const std::vector<size_t>& indices,
                  StorageBackend& destination,
                  size_t num_threads, size_t chunk_size) {
    const size_t m = indices.size();
    if (m < 2) {
        return;  // one item owns no pair
    }
    ThreadPool pool(num_threads);
    // More workers than rows do nothing; capping before the arithmetic keeps
    // 4 * workers from wrapping for any count the size_t interface admits.
    const size_t workers = std::max<size_t>(1, std::min(pool.NumThreads(), m));
    // chunk_size is a ceiling, not the chunk: a subset small enough to fit
    // one chunk still spreads over the pool, and the shrinking
    // upper-triangular rows stay within a bounded imbalance.
    const size_t spread = (m + 4 * workers - 1) / (4 * workers);
    const size_t rows_per_chunk =
        std::max<size_t>(1, std::min(chunk_size, spread));
    pool.ParallelFor(0, m, rows_per_chunk, [&](size_t begin, size_t end) {
        for (size_t a = begin; a < end; ++a) {
            const size_t source_a = indices[a];
            for (size_t b = a + 1; b < m; ++b) {
                // Distinct pairs land in distinct slots, which is how pdist
                // already writes from several threads.
                destination.Set(a, b, source.Get(source_a, indices[b]));
            }
        }
    });
}

void gather_sparse(const SparseStorage& source,
                   const std::vector<size_t>& indices,
                   SparseStorage& destination) {
    std::vector<size_t> position(source.NumSamples(), ABSENT);
    for (size_t a = 0; a < indices.size(); ++a) {
        position[indices[a]] = a;
    }
    // Entries() is the merged list, so its order and any superseded
    // duplicates carry over verbatim, as to_file keeps them.
    for (const auto& entry : source.Entries()) {
        const size_t a = position[std::get<0>(entry)];
        const size_t b = position[std::get<1>(entry)];
        if (a == ABSENT || b == ABSENT) {
            continue;
        }
        destination.Set(a, b, std::get<2>(entry));
    }
}

}  // namespace

void take_pairs(const StorageBackend& source,
                const std::vector<size_t>& indices,
                StorageBackend& destination,
                size_t num_threads,
                size_t chunk_size) {
    check_indices(indices, source.NumSamples(), "take_pairs");
    if (destination.NumSamples() != indices.size()) {
        throw std::invalid_argument(
            "take_pairs: destination holds "
            + std::to_string(destination.NumSamples()) + " items but "
            + std::to_string(indices.size()) + " indices were given");
    }
    if (chunk_size == 0) {
        throw std::invalid_argument("take_pairs: chunk_size must be positive");
    }
    if (&destination == &source) {
        throw std::invalid_argument(
            "take_pairs: destination must not be the source; a gather into "
            "itself would race its own reads");
    }
    // A memory-mapped destination is refused with every other backend: two
    // MMapStorage objects can map one file, which the identity check above
    // cannot see, and an unknown backend may not honour concurrent Set calls.
    // The subset is small enough to own in memory, which is what the Python
    // layer always allocates.
    const auto* sparse_source = dynamic_cast<const SparseStorage*>(&source);
    auto* sparse_destination = dynamic_cast<SparseStorage*>(&destination);
    if (sparse_source != nullptr) {
        if (sparse_destination == nullptr) {
            throw std::invalid_argument(
                "take_pairs: a sparse source needs a SparseStorage destination");
        }
    } else if (dynamic_cast<DenseStorage*>(&destination) == nullptr) {
        throw std::invalid_argument(
            "take_pairs: a dense or memory-mapped source needs a DenseStorage "
            "destination; a memory-mapped destination could share the "
            "source's file and other backends are not supported");
    }
    if (sparse_source != nullptr) {
        const double source_cutoff = sparse_source->Cutoff();
        const double destination_cutoff = sparse_destination->Cutoff();
        // A NaN cutoff is legal (nothing is ever above it) and must match itself.
        const bool same_cutoff =
            source_cutoff == destination_cutoff
            || (std::isnan(source_cutoff) && std::isnan(destination_cutoff));
        if (!same_cutoff) {
            throw std::invalid_argument(
                "take_pairs: the sparse destination's cutoff must equal the "
                "source's; a lower one would drop entries and a higher one "
                "would claim coverage the source never had");
        }
        // Finalize() merges rather than clears, so stale entries would
        // survive the gather, and Entries() shows only merged ones: a write
        // still sitting in a thread buffer is invisible until Finalize()
        // folds it in. Finalizing first makes every earlier write visible to
        // the check; on a fresh destination it is a no-op. A dense
        // destination needs no such check: every one of its pairs is
        // overwritten.
        sparse_destination->Finalize();
        if (!sparse_destination->Entries().empty()) {
            throw std::invalid_argument(
                "take_pairs: the sparse destination already holds entries");
        }
        gather_sparse(*sparse_source, indices, *sparse_destination);
    } else {
        gather_dense(source, indices, destination, num_threads, chunk_size);
    }
    destination.Finalize();
}

OEFP::OEFPBatch take_fingerprints(const OEFP::OEFPBatch& batch,
                                  const std::vector<size_t>& indices) {
    check_indices(indices, batch.Size(), "take_fingerprints");
    const size_t words = batch.WordsPerFingerprint();
    // Append reserves exactly one more row per call, so m appends reallocate
    // m times; FromFingerprints reserves once. check_indices guarantees the
    // vector is non-empty.
    std::vector<OEFP::OEFP> rows;
    rows.reserve(indices.size());
    for (size_t index : indices) {
        const std::uint64_t* row = batch.RowWords(index);
        rows.emplace_back(batch.Spec(),
                          std::vector<std::uint64_t>(row, row + words));
    }
    return OEFP::OEFPBatch::FromFingerprints(rows);
}

}  // namespace OECluster
