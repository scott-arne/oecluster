/**
 * @file BitBirchTree.h
 * @brief Internal BitBirch tree state.
 */

#ifndef OECLUSTER_CLUSTERING_BITBIRCH_TREE_H
#define OECLUSTER_CLUSTERING_BITBIRCH_TREE_H

#include <memory>
#include <utility>
#include <vector>

#include "oefp/batch.h"
#include "oecluster/clustering/BitBirch.h"
#include "BitBirchKernels.h"

namespace OECluster::detail {

class BitBirchNode;

struct BitBirchSubcluster {
    size_t n_samples = 0;
    BitBirchLinearSum linear_sum;
    std::vector<uint64_t> centroid_words;
    uint32_t centroid_popcount = 0;
    Cluster members;
    BitBirchNode* child = nullptr;
    BitBirchSubcluster* parent = nullptr;

    BitBirchSubcluster() = default;
    BitBirchSubcluster(
        BitBirchLinearSum linear_sum,
        std::vector<uint64_t> centroid_words,
        uint32_t centroid_popcount,
        Cluster members);

    void UpdateFrom(const BitBirchSubcluster& other);
    void Subtract(const BitBirchSubcluster& removed);
    bool TryMerge(
        const BitBirchSubcluster& nominee,
        const BitBirchOptions& options);
};

struct BitBirchPruneSelection {
    BitBirchNode* leaf = nullptr;
    size_t subcluster_index = 0;
    BitBirchSubcluster* subcluster = nullptr;
};

class BitBirchNode {
public:
    BitBirchNode(size_t branching_factor, bool is_leaf);

    bool InsertSubcluster(
        std::unique_ptr<BitBirchSubcluster> subcluster,
        const BitBirchOptions& options,
        BitBirchSubcluster* parent,
        std::vector<std::unique_ptr<BitBirchNode>>& node_owner);

    void AppendSubcluster(std::unique_ptr<BitBirchSubcluster> subcluster);
    void ReplaceWithSplit(
        const BitBirchSubcluster* existing,
        std::unique_ptr<BitBirchSubcluster> left,
        std::unique_ptr<BitBirchSubcluster> right);

    bool IsLeaf() const;
    void SetLeaf(bool is_leaf);
    bool HasChildSubclusters() const;
    size_t Size() const;
    const std::vector<std::unique_ptr<BitBirchSubcluster>>& Subclusters() const;
    std::vector<std::unique_ptr<BitBirchSubcluster>>& MutableSubclusters();

    BitBirchNode* prev_leaf = nullptr;
    BitBirchNode* next_leaf = nullptr;

private:
    size_t branching_factor_;
    bool is_leaf_;
    std::vector<std::unique_ptr<BitBirchSubcluster>> subclusters_;
};

class BitBirchTree {
public:
    explicit BitBirchTree(const BitBirchOptions& options);

    void Fit(const OEFP::OEFPBatch& fingerprints);
    /**
     * @brief Fit rows [begin, end) of the batch, storing global member ids.
     *
     * Identical to fitting the whole batch when begin==0 and end==Size(). Used by
     * the fast-mode parallel engine so each chunk tree carries original row
     * indices into the merge phase.
     */
    void Fit(const OEFP::OEFPBatch& fingerprints, size_t begin, size_t end);
    void FitSubclusters(std::vector<std::unique_ptr<BitBirchSubcluster>> subclusters);
    std::pair<
        std::vector<std::unique_ptr<BitBirchSubcluster>>,
        std::vector<std::unique_ptr<BitBirchSubcluster>>>
    PrepareReclusteringSubclusters(const OEFP::OEFPBatch& fingerprints) const;
    void ReassignTopClusters(
        const OEFP::OEFPBatch& fingerprints,
        size_t top_clusters,
        size_t num_threads);
    Cluster RedistributeLargestCluster(const OEFP::OEFPBatch& fingerprints);
    std::vector<const BitBirchSubcluster*> LeafSubclusters() const;
    BitBirchResult Result(const OEFP::FingerprintSpec& spec, size_t n_items) const;

    /**
     * @brief Fill out_tree with the fast-mode (partition + parallel fit + merge)
     *        clustering of the batch.
     *
     * Deterministic in (fingerprints, options) alone; options.num_threads affects
     * only concurrency, not the result. For partition_count(N)==1 this is exactly
     * the sequential Fit. out_tree must be freshly constructed with options.
     */
    static void BuildFastTree(const OEFP::OEFPBatch& fingerprints,
                              const BitBirchOptions& options,
                              BitBirchTree& out_tree);

private:
    BitBirchOptions options_;
    size_t index_tracker_ = 0;
    size_t size_bits_ = 0;
    size_t words_per_fingerprint_ = 0;
    std::vector<std::unique_ptr<BitBirchNode>> nodes_;
    std::unique_ptr<BitBirchNode> dummy_leaf_;
    BitBirchNode* root_ = nullptr;

    void EnsureInitialized(size_t size_bits, size_t words_per_fingerprint);
    BitBirchNode* MakeNode(bool is_leaf);
    std::unique_ptr<BitBirchSubcluster> MakeSubclusterFromRow(
        const OEFP::OEFPBatch& fingerprints,
        size_t row) const;
    std::unique_ptr<BitBirchSubcluster> MakeSubclusterForMember(
        const OEFP::OEFPBatch& fingerprints,
        size_t member) const;
    void SplitRootIfNeeded(bool split);
    BitBirchPruneSelection SelectLargestLeafSubcluster();
    Cluster RemoveSelectedSubcluster(const BitBirchPruneSelection& selection);
    void RelinkAfterLeafRemoval(BitBirchNode* leaf, BitBirchSubcluster* parent);
    void SubtractFromAncestors(
        BitBirchSubcluster* parent,
        const BitBirchSubcluster& removed);
    void ReinsertMembers(const OEFP::OEFPBatch& fingerprints, const Cluster& members);
};

std::pair<std::unique_ptr<BitBirchSubcluster>, std::unique_ptr<BitBirchSubcluster>>
SplitNode(
    BitBirchNode& node,
    const BitBirchOptions& options,
    std::vector<std::unique_ptr<BitBirchNode>>& node_owner);

/**
 * @brief Number of contiguous fast-mode partitions for an input of size n.
 *
 * Deterministic in n alone (independent of thread count) so fast-mode output
 * does not depend on num_threads. Returns 1 for inputs at or below the chunk
 * target, which makes fast mode byte-identical to strict parity for small n.
 */
size_t partition_count(size_t n);

}  // namespace OECluster::detail

#endif  // OECLUSTER_CLUSTERING_BITBIRCH_TREE_H
