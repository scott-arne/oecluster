/**
 * @file BitBirchTree.cpp
 * @brief Internal BitBirch tree state.
 */

#include "BitBirchTree.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <map>
#include <stdexcept>

#include "oecluster/ThreadPool.h"

namespace OECluster::detail {
namespace {

BitBirchLinearSum add_linear_sums(
    const BitBirchLinearSum& lhs,
    const BitBirchLinearSum& rhs) {
    if (lhs.size() != rhs.size()) {
        throw std::invalid_argument("BitBirch linear sums must have the same width");
    }
    BitBirchLinearSum out(lhs.size(), 0);
    for (size_t i = 0; i < lhs.size(); ++i) {
        out[i] = lhs[i] + rhs[i];
    }
    return out;
}

bool numpy_argmax_greater(const double candidate, const double best) {
    if (std::isnan(candidate)) {
        return !std::isnan(best);
    }
    if (std::isnan(best)) {
        return false;
    }
    return candidate > best;
}

bool numpy_argmin_less(const double candidate, const double best) {
    if (std::isnan(candidate)) {
        return !std::isnan(best);
    }
    if (std::isnan(best)) {
        return false;
    }
    return candidate < best;
}

size_t argmax_numpy(const std::vector<double>& values) {
    size_t best_index = 0;
    double best_value = values.front();
    for (size_t index = 1; index < values.size(); ++index) {
        if (numpy_argmax_greater(values[index], best_value)) {
            best_index = index;
            best_value = values[index];
        }
    }
    return best_index;
}

size_t argmin_numpy(const std::vector<double>& values) {
    size_t best_index = 0;
    double best_value = values.front();
    for (size_t index = 1; index < values.size(); ++index) {
        if (numpy_argmin_less(values[index], best_value)) {
            best_index = index;
            best_value = values[index];
        }
    }
    return best_index;
}

std::vector<double> similarities_to_centroid(
    const std::vector<std::unique_ptr<BitBirchSubcluster>>& subclusters,
    const std::vector<uint64_t>& centroid_words,
    const uint32_t centroid_popcount) {
    std::vector<double> similarities;
    similarities.reserve(subclusters.size());
    for (const auto& subcluster : subclusters) {
        similarities.push_back(TanimotoPackedVectors(
            subcluster->centroid_words,
            subcluster->centroid_popcount,
            centroid_words,
            centroid_popcount));
    }
    return similarities;
}

size_t closest_subcluster_index(
    const std::vector<std::unique_ptr<BitBirchSubcluster>>& subclusters,
    const BitBirchSubcluster& nominee) {
    size_t best_index = 0;
    double best_similarity = TanimotoPackedVectors(
        subclusters.front()->centroid_words,
        subclusters.front()->centroid_popcount,
        nominee.centroid_words,
        nominee.centroid_popcount);

    for (size_t index = 1; index < subclusters.size(); ++index) {
        const double similarity = TanimotoPackedVectors(
            subclusters[index]->centroid_words,
            subclusters[index]->centroid_popcount,
            nominee.centroid_words,
            nominee.centroid_popcount);
        if (numpy_argmax_greater(similarity, best_similarity)) {
            best_index = index;
            best_similarity = similarity;
        }
    }
    return best_index;
}

std::pair<size_t, size_t> max_separation(
    const std::vector<std::unique_ptr<BitBirchSubcluster>>& subclusters) {
    if (subclusters.size() < 2) {
        throw std::invalid_argument("BitBirch split requires at least two subclusters");
    }

    BitBirchLinearSum centroid_sum(subclusters.front()->linear_sum.size(), 0);
    for (const auto& subcluster : subclusters) {
        UpdateLinearSumFromWords(
            centroid_sum,
            subcluster->centroid_words.data(),
            centroid_sum.size());
    }
    const std::vector<uint64_t> centroid_words =
        BinaryCentroid(centroid_sum, subclusters.size());
    const uint32_t centroid_popcount = PopCountWords(centroid_words);
    const std::vector<double> sims_to_centroid =
        similarities_to_centroid(subclusters, centroid_words, centroid_popcount);
    const size_t first_seed = argmin_numpy(sims_to_centroid);

    const auto& first = subclusters[first_seed];
    const std::vector<double> sims_to_first = similarities_to_centroid(
        subclusters,
        first->centroid_words,
        first->centroid_popcount);
    const size_t second_seed = argmin_numpy(sims_to_first);
    return {first_seed, second_seed};
}

std::unique_ptr<BitBirchNode> make_node_like(
    const BitBirchNode& node,
    const BitBirchOptions& options) {
    return std::make_unique<BitBirchNode>(options.branching_factor, node.IsLeaf());
}

std::unique_ptr<BitBirchSubcluster> clone_leaf_subcluster(
    const BitBirchSubcluster& subcluster) {
    return std::make_unique<BitBirchSubcluster>(
        subcluster.linear_sum,
        subcluster.centroid_words,
        subcluster.centroid_popcount,
        subcluster.members);
}

}  // namespace

BitBirchSubcluster::BitBirchSubcluster(
    BitBirchLinearSum linear_sum_in,
    std::vector<uint64_t> centroid_words_in,
    const uint32_t centroid_popcount_in,
    Cluster members_in)
    : n_samples(members_in.size()),
      linear_sum(std::move(linear_sum_in)),
      centroid_words(std::move(centroid_words_in)),
      centroid_popcount(centroid_popcount_in),
      members(std::move(members_in)) {}

void BitBirchSubcluster::UpdateFrom(const BitBirchSubcluster& other) {
    if (n_samples == 0) {
        n_samples = other.n_samples;
        linear_sum = other.linear_sum;
        centroid_words = other.centroid_words;
        centroid_popcount = other.centroid_popcount;
        members = other.members;
        return;
    }

    linear_sum = add_linear_sums(linear_sum, other.linear_sum);
    n_samples += other.n_samples;
    members.insert(members.end(), other.members.begin(), other.members.end());
    centroid_words = BinaryCentroid(linear_sum, n_samples);
    centroid_popcount = PopCountWords(centroid_words);
}

void BitBirchSubcluster::Subtract(const BitBirchSubcluster& removed) {
    if (removed.n_samples > n_samples) {
        throw std::invalid_argument(
            "BitBirch prune cannot subtract more samples than a subcluster contains");
    }
    if (linear_sum.size() != removed.linear_sum.size()) {
        throw std::invalid_argument("BitBirch prune linear-sum widths differ");
    }
    for (size_t index = 0; index < linear_sum.size(); ++index) {
        if (removed.linear_sum[index] > linear_sum[index]) {
            throw std::invalid_argument(
                "BitBirch prune would make a linear sum negative");
        }
        linear_sum[index] -= removed.linear_sum[index];
    }
    n_samples -= removed.n_samples;
    members.erase(
        std::remove_if(
            members.begin(),
            members.end(),
            [&removed](const size_t member) {
                return std::find(
                           removed.members.begin(),
                           removed.members.end(),
                           member) != removed.members.end();
            }),
        members.end());
    centroid_words = BinaryCentroid(linear_sum, n_samples);
    centroid_popcount = PopCountWords(centroid_words);
}

bool BitBirchSubcluster::TryMerge(
    const BitBirchSubcluster& nominee,
    const BitBirchOptions& options) {
    BitBirchLinearSum new_linear_sum = add_linear_sums(linear_sum, nominee.linear_sum);
    const size_t new_n = n_samples + nominee.n_samples;
    const std::vector<uint64_t> new_centroid = BinaryCentroid(new_linear_sum, new_n);

    if (!AcceptBitBirchMerge(
            options.merge_criterion,
            options.threshold,
            options.tolerance,
            new_linear_sum,
            new_centroid,
            new_n,
            linear_sum,
            nominee.linear_sum,
            n_samples,
            nominee.n_samples)) {
        return false;
    }

    n_samples = new_n;
    linear_sum = std::move(new_linear_sum);
    centroid_words = new_centroid;
    centroid_popcount = PopCountWords(centroid_words);
    members.insert(members.end(), nominee.members.begin(), nominee.members.end());
    return true;
}

BitBirchNode::BitBirchNode(const size_t branching_factor, const bool is_leaf)
    : branching_factor_(branching_factor), is_leaf_(is_leaf) {}

bool BitBirchNode::InsertSubcluster(
    std::unique_ptr<BitBirchSubcluster> subcluster,
    const BitBirchOptions& options,
    BitBirchSubcluster* parent,
    std::vector<std::unique_ptr<BitBirchNode>>& node_owner) {
    if (subclusters_.empty()) {
        subcluster->parent = parent;
        AppendSubcluster(std::move(subcluster));
        return false;
    }

    const size_t closest_index = closest_subcluster_index(subclusters_, *subcluster);
    BitBirchSubcluster* closest = subclusters_[closest_index].get();

    if (closest->child != nullptr) {
        const BitBirchSubcluster inserted_snapshot = *subcluster;
        const bool split_child = closest->child->InsertSubcluster(
            std::move(subcluster),
            options,
            closest,
            node_owner);

        if (!split_child) {
            closest->UpdateFrom(inserted_snapshot);
            return false;
        }

        auto [left, right] = SplitNode(*closest->child, options, node_owner);
        if (!options.singly) {
            left->parent = closest->parent;
            right->parent = closest->parent;
        }
        ReplaceWithSplit(closest, std::move(left), std::move(right));
        return subclusters_.size() > branching_factor_;
    }

    if (closest->TryMerge(*subcluster, options)) {
        closest->parent = options.singly ? closest->parent : parent;
        return false;
    }

    if (subclusters_.size() < branching_factor_) {
        subcluster->parent = options.singly ? subcluster->parent : parent;
        AppendSubcluster(std::move(subcluster));
        return false;
    }

    subcluster->parent = options.singly ? subcluster->parent : parent;
    AppendSubcluster(std::move(subcluster));
    return true;
}

void BitBirchNode::AppendSubcluster(std::unique_ptr<BitBirchSubcluster> subcluster) {
    subclusters_.push_back(std::move(subcluster));
}

void BitBirchNode::ReplaceWithSplit(
    const BitBirchSubcluster* existing,
    std::unique_ptr<BitBirchSubcluster> left,
    std::unique_ptr<BitBirchSubcluster> right) {
    for (size_t index = 0; index < subclusters_.size(); ++index) {
        if (subclusters_[index].get() == existing) {
            subclusters_[index] = std::move(left);
            AppendSubcluster(std::move(right));
            return;
        }
    }
    throw std::invalid_argument("BitBirch split target was not found");
}

bool BitBirchNode::IsLeaf() const {
    return is_leaf_;
}

void BitBirchNode::SetLeaf(const bool is_leaf) {
    is_leaf_ = is_leaf;
}

bool BitBirchNode::HasChildSubclusters() const {
    return std::any_of(
        subclusters_.begin(),
        subclusters_.end(),
        [](const std::unique_ptr<BitBirchSubcluster>& subcluster) {
            return subcluster->child != nullptr;
        });
}

size_t BitBirchNode::Size() const {
    return subclusters_.size();
}

const std::vector<std::unique_ptr<BitBirchSubcluster>>& BitBirchNode::Subclusters() const {
    return subclusters_;
}

std::vector<std::unique_ptr<BitBirchSubcluster>>& BitBirchNode::MutableSubclusters() {
    return subclusters_;
}

BitBirchTree::BitBirchTree(const BitBirchOptions& options) : options_(options) {}

void BitBirchTree::Fit(const OEFP::OEFPBatch& fingerprints) {
    EnsureInitialized(fingerprints.SizeBits(), fingerprints.WordsPerFingerprint());
    for (size_t row = 0; row < fingerprints.Size(); ++row) {
        std::unique_ptr<BitBirchSubcluster> subcluster =
            MakeSubclusterFromRow(fingerprints, row);
        const bool split = root_->InsertSubcluster(
            std::move(subcluster),
            options_,
            nullptr,
            nodes_);
        SplitRootIfNeeded(split);
        ++index_tracker_;
    }
}

void BitBirchTree::FitSubclusters(
    std::vector<std::unique_ptr<BitBirchSubcluster>> subclusters) {
    if (subclusters.empty()) {
        return;
    }
    EnsureInitialized(subclusters.front()->linear_sum.size(), subclusters.front()->centroid_words.size());
    for (auto& subcluster : subclusters) {
        const bool split = root_->InsertSubcluster(
            std::move(subcluster),
            options_,
            nullptr,
            nodes_);
        SplitRootIfNeeded(split);
        ++index_tracker_;
    }
}

std::pair<
    std::vector<std::unique_ptr<BitBirchSubcluster>>,
    std::vector<std::unique_ptr<BitBirchSubcluster>>>
BitBirchTree::PrepareReclusteringSubclusters(
    const OEFP::OEFPBatch& fingerprints) const {
    std::vector<const BitBirchSubcluster*> ordered = LeafSubclusters();
    std::stable_sort(
        ordered.begin(),
        ordered.end(),
        [](const BitBirchSubcluster* lhs, const BitBirchSubcluster* rhs) {
            return lhs->n_samples > rhs->n_samples;
        });
    if (ordered.empty()) {
        return {};
    }

    std::vector<std::unique_ptr<BitBirchSubcluster>> rest;
    rest.reserve(ordered.size() - 1);
    for (size_t index = 1; index < ordered.size(); ++index) {
        rest.push_back(clone_leaf_subcluster(*ordered[index]));
    }

    std::vector<std::unique_ptr<BitBirchSubcluster>> largest_singletons;
    largest_singletons.reserve(ordered.front()->members.size());
    for (const size_t member : ordered.front()->members) {
        largest_singletons.push_back(MakeSubclusterForMember(fingerprints, member));
    }
    return {std::move(rest), std::move(largest_singletons)};
}

void BitBirchTree::ReassignTopClusters(
    const OEFP::OEFPBatch& fingerprints,
    const size_t top_clusters,
    const size_t num_threads) {
    if (top_clusters == 0) {
        return;
    }

    std::vector<BitBirchSubcluster*> ordered;
    if (dummy_leaf_) {
        BitBirchNode* leaf = dummy_leaf_->next_leaf;
        while (leaf != nullptr) {
            for (auto& subcluster : leaf->MutableSubclusters()) {
                ordered.push_back(subcluster.get());
            }
            leaf = leaf->next_leaf;
        }
    }
    std::stable_sort(
        ordered.begin(),
        ordered.end(),
        [](const BitBirchSubcluster* lhs, const BitBirchSubcluster* rhs) {
            return lhs->members.size() > rhs->members.size();
        });
    if (ordered.empty()) {
        return;
    }

    const size_t selected_count = std::min(top_clusters, ordered.size());
    std::vector<std::vector<uint64_t>> centroids;
    std::vector<uint32_t> centroid_popcounts;
    centroids.reserve(selected_count);
    centroid_popcounts.reserve(selected_count);

    Cluster selected_members;
    for (size_t index = 0; index < selected_count; ++index) {
        BitBirchLinearSum linear_sum(size_bits_, 0);
        for (const size_t member : ordered[index]->members) {
            UpdateLinearSumFromWords(linear_sum, fingerprints.RowWords(member), size_bits_);
            selected_members.push_back(member);
        }
        centroids.push_back(BinaryCentroid(linear_sum, ordered[index]->members.size()));
        centroid_popcounts.push_back(PopCountWords(centroids.back()));
    }

    auto assign_member = [&](const size_t member) {
        size_t best_index = 0;
        double best_similarity = TanimotoWords(
            fingerprints.RowWords(member),
            fingerprints.PopCount(member),
            centroids.front().data(),
            centroid_popcounts.front(),
            words_per_fingerprint_);
        for (size_t centroid_index = 1; centroid_index < selected_count; ++centroid_index) {
            const double similarity = TanimotoWords(
                fingerprints.RowWords(member),
                fingerprints.PopCount(member),
                centroids[centroid_index].data(),
                centroid_popcounts[centroid_index],
                words_per_fingerprint_);
            if (numpy_argmax_greater(similarity, best_similarity)) {
                best_index = centroid_index;
                best_similarity = similarity;
            }
        }
        return best_index;
    };

    std::vector<size_t> assigned_centroids(selected_members.size(), 0);
    if (selected_members.size() >= 256 && (num_threads == 0 || num_threads > 1)) {
        ThreadPool pool(num_threads);
        pool.ParallelFor(
            0,
            selected_members.size(),
            256,
            [&](const size_t begin, const size_t end) {
                for (size_t index = begin; index < end; ++index) {
                    assigned_centroids[index] = assign_member(selected_members[index]);
                }
            });
    } else {
        for (size_t index = 0; index < selected_members.size(); ++index) {
            assigned_centroids[index] = assign_member(selected_members[index]);
        }
    }

    std::map<size_t, Cluster> reassigned;
    for (size_t index = 0; index < selected_members.size(); ++index) {
        reassigned[assigned_centroids[index]].push_back(selected_members[index]);
    }

    size_t subcluster_index = 0;
    for (const auto& [_, members] : reassigned) {
        if (subcluster_index >= selected_count) {
            break;
        }
        BitBirchSubcluster* subcluster = ordered[subcluster_index];
        BitBirchLinearSum linear_sum(size_bits_, 0);
        for (const size_t member : members) {
            UpdateLinearSumFromWords(linear_sum, fingerprints.RowWords(member), size_bits_);
        }
        subcluster->members = members;
        subcluster->n_samples = members.size();
        subcluster->linear_sum = std::move(linear_sum);
        subcluster->centroid_words =
            BinaryCentroid(subcluster->linear_sum, subcluster->n_samples);
        subcluster->centroid_popcount = PopCountWords(subcluster->centroid_words);
        ++subcluster_index;
    }
}

Cluster BitBirchTree::RedistributeLargestCluster(const OEFP::OEFPBatch& fingerprints) {
    const BitBirchPruneSelection selection = SelectLargestLeafSubcluster();
    if (selection.subcluster == nullptr) {
        return {};
    }
    if (selection.subcluster->parent == nullptr) {
        throw std::invalid_argument(
            "BitBirch redistribute_largest_cluster requires a split tree");
    }

    const BitBirchSubcluster removed_snapshot = *selection.subcluster;
    BitBirchSubcluster* parent = selection.subcluster->parent;
    Cluster removed_members = RemoveSelectedSubcluster(selection);
    SubtractFromAncestors(parent, removed_snapshot);
    ReinsertMembers(fingerprints, removed_members);
    return removed_members;
}

std::vector<const BitBirchSubcluster*> BitBirchTree::LeafSubclusters() const {
    std::vector<const BitBirchSubcluster*> subclusters;
    if (!dummy_leaf_) {
        return subclusters;
    }

    const BitBirchNode* leaf = dummy_leaf_->next_leaf;
    while (leaf != nullptr) {
        for (const auto& subcluster : leaf->Subclusters()) {
            subclusters.push_back(subcluster.get());
        }
        leaf = leaf->next_leaf;
    }
    return subclusters;
}

BitBirchResult BitBirchTree::Result(
    const OEFP::FingerprintSpec& spec,
    const size_t n_items) const {
    BitBirchResult result;
    result.labels.assign(n_items, NOISE_LABEL);

    std::vector<const BitBirchSubcluster*> ordered = LeafSubclusters();
    std::stable_sort(
        ordered.begin(),
        ordered.end(),
        [](const BitBirchSubcluster* lhs, const BitBirchSubcluster* rhs) {
            return lhs->members.size() > rhs->members.size();
        });

    std::vector<OEFP::OEFP> centroid_fingerprints;
    centroid_fingerprints.reserve(ordered.size());
    for (size_t label = 0; label < ordered.size(); ++label) {
        const BitBirchSubcluster& subcluster = *ordered[label];
        result.clusters.push_back(subcluster.members);
        result.cluster_sizes.push_back(subcluster.members.size());
        for (const size_t member : subcluster.members) {
            result.labels.at(member) = static_cast<ClusterLabel>(label);
        }
        centroid_fingerprints.emplace_back(spec, subcluster.centroid_words);
    }

    result.centroids = centroid_fingerprints.empty()
        ? OEFP::OEFPBatch(spec)
        : OEFP::OEFPBatch::FromFingerprints(centroid_fingerprints);
    return result;
}

void BitBirchTree::EnsureInitialized(
    const size_t size_bits,
    const size_t words_per_fingerprint) {
    if (root_ != nullptr) {
        return;
    }
    size_bits_ = size_bits;
    words_per_fingerprint_ = words_per_fingerprint;
    root_ = MakeNode(true);
    dummy_leaf_ = std::make_unique<BitBirchNode>(options_.branching_factor, true);
    dummy_leaf_->next_leaf = root_;
    root_->prev_leaf = dummy_leaf_.get();
}

BitBirchNode* BitBirchTree::MakeNode(const bool is_leaf) {
    nodes_.push_back(std::make_unique<BitBirchNode>(options_.branching_factor, is_leaf));
    return nodes_.back().get();
}

std::unique_ptr<BitBirchSubcluster> BitBirchTree::MakeSubclusterFromRow(
    const OEFP::OEFPBatch& fingerprints,
    const size_t row) const {
    return MakeSubclusterForMember(fingerprints, row);
}

std::unique_ptr<BitBirchSubcluster> BitBirchTree::MakeSubclusterForMember(
    const OEFP::OEFPBatch& fingerprints,
    const size_t member) const {
    BitBirchLinearSum linear_sum(size_bits_, 0);
    UpdateLinearSumFromWords(linear_sum, fingerprints.RowWords(member), size_bits_);

    const uint64_t* row_words = fingerprints.RowWords(member);
    std::vector<uint64_t> centroid_words(
        row_words,
        row_words + words_per_fingerprint_);
    Cluster members{member};
    return std::make_unique<BitBirchSubcluster>(
        std::move(linear_sum),
        std::move(centroid_words),
        fingerprints.PopCount(member),
        std::move(members));
}

void BitBirchTree::SplitRootIfNeeded(const bool split) {
    if (!split) {
        return;
    }

    auto [left, right] = SplitNode(*root_, options_, nodes_);
    root_ = MakeNode(false);
    root_->AppendSubcluster(std::move(left));
    root_->AppendSubcluster(std::move(right));
}

BitBirchPruneSelection BitBirchTree::SelectLargestLeafSubcluster() {
    BitBirchPruneSelection best;
    if (!dummy_leaf_) {
        return best;
    }

    BitBirchNode* leaf = dummy_leaf_->next_leaf;
    while (leaf != nullptr) {
        auto& subclusters = leaf->MutableSubclusters();
        for (size_t index = 0; index < subclusters.size(); ++index) {
            BitBirchSubcluster* candidate = subclusters[index].get();
            if (best.subcluster == nullptr ||
                candidate->members.size() > best.subcluster->members.size()) {
                best.leaf = leaf;
                best.subcluster_index = index;
                best.subcluster = candidate;
            }
        }
        leaf = leaf->next_leaf;
    }
    return best;
}

Cluster BitBirchTree::RemoveSelectedSubcluster(
    const BitBirchPruneSelection& selection) {
    if (selection.leaf == nullptr || selection.subcluster == nullptr) {
        return {};
    }

    auto& subclusters = selection.leaf->MutableSubclusters();
    if (selection.subcluster_index >= subclusters.size() ||
        subclusters[selection.subcluster_index].get() != selection.subcluster) {
        throw std::invalid_argument("BitBirch prune selection is stale");
    }

    BitBirchSubcluster* parent = selection.subcluster->parent;
    Cluster removed_members = selection.subcluster->members;
    subclusters.erase(
        subclusters.begin() +
        static_cast<std::ptrdiff_t>(selection.subcluster_index));

    if (subclusters.empty()) {
        RelinkAfterLeafRemoval(selection.leaf, parent);
    }
    return removed_members;
}

void BitBirchTree::RelinkAfterLeafRemoval(
    BitBirchNode* leaf,
    BitBirchSubcluster* parent) {
    if (parent == nullptr) {
        if (leaf->prev_leaf != nullptr) {
            leaf->prev_leaf->next_leaf = leaf->next_leaf;
        }
        if (leaf->next_leaf != nullptr) {
            leaf->next_leaf->prev_leaf = leaf->prev_leaf;
        }
        return;
    }

    parent->child = nullptr;
    BitBirchNode* parent_node =
        parent->parent == nullptr ? root_ : parent->parent->child;
    if (parent_node != nullptr && !parent_node->HasChildSubclusters()) {
        parent_node->SetLeaf(true);
        if (leaf->prev_leaf != nullptr) {
            leaf->prev_leaf->next_leaf = parent_node;
        }
        parent_node->prev_leaf = leaf->prev_leaf;
        parent_node->next_leaf = leaf->next_leaf;
        if (leaf->next_leaf != nullptr) {
            leaf->next_leaf->prev_leaf = parent_node;
        }
        return;
    }

    if (leaf->prev_leaf != nullptr) {
        leaf->prev_leaf->next_leaf = leaf->next_leaf;
    }
    if (leaf->next_leaf != nullptr) {
        leaf->next_leaf->prev_leaf = leaf->prev_leaf;
    }
}

void BitBirchTree::SubtractFromAncestors(
    BitBirchSubcluster* parent,
    const BitBirchSubcluster& removed) {
    while (parent != nullptr) {
        parent->Subtract(removed);
        parent = parent->parent;
    }
}

void BitBirchTree::ReinsertMembers(
    const OEFP::OEFPBatch& fingerprints,
    const Cluster& members) {
    for (const size_t member : members) {
        std::unique_ptr<BitBirchSubcluster> subcluster =
            MakeSubclusterForMember(fingerprints, member);
        const bool split = root_->InsertSubcluster(
            std::move(subcluster),
            options_,
            nullptr,
            nodes_);
        SplitRootIfNeeded(split);
        ++index_tracker_;
    }
}

std::pair<std::unique_ptr<BitBirchSubcluster>, std::unique_ptr<BitBirchSubcluster>>
SplitNode(
    BitBirchNode& node,
    const BitBirchOptions& options,
    std::vector<std::unique_ptr<BitBirchNode>>& node_owner) {
    auto left_node = make_node_like(node, options);
    auto right_node = make_node_like(node, options);
    BitBirchNode* left_node_ptr = left_node.get();
    BitBirchNode* right_node_ptr = right_node.get();
    node_owner.push_back(std::move(left_node));
    node_owner.push_back(std::move(right_node));

    if (node.IsLeaf()) {
        if (node.prev_leaf != nullptr) {
            node.prev_leaf->next_leaf = left_node_ptr;
        }
        left_node_ptr->prev_leaf = node.prev_leaf;
        left_node_ptr->next_leaf = right_node_ptr;
        right_node_ptr->prev_leaf = left_node_ptr;
        right_node_ptr->next_leaf = node.next_leaf;
        if (node.next_leaf != nullptr) {
            node.next_leaf->prev_leaf = right_node_ptr;
        }
    }

    const auto [first_seed, second_seed] = max_separation(node.Subclusters());
    const auto& first = node.Subclusters()[first_seed];
    const auto& second = node.Subclusters()[second_seed];
    const std::vector<double> sims_to_first = similarities_to_centroid(
        node.Subclusters(),
        first->centroid_words,
        first->centroid_popcount);
    const std::vector<double> sims_to_second = similarities_to_centroid(
        node.Subclusters(),
        second->centroid_words,
        second->centroid_popcount);

    auto left_summary = std::make_unique<BitBirchSubcluster>();
    auto right_summary = std::make_unique<BitBirchSubcluster>();
    left_summary->child = left_node_ptr;
    right_summary->child = right_node_ptr;
    BitBirchSubcluster* left_summary_ptr = left_summary.get();
    BitBirchSubcluster* right_summary_ptr = right_summary.get();

    std::vector<std::unique_ptr<BitBirchSubcluster>> moved;
    moved.swap(node.MutableSubclusters());
    for (size_t index = 0; index < moved.size(); ++index) {
        const bool node_one_closer =
            (index == first_seed) || (sims_to_first[index] > sims_to_second[index]);
        if (node_one_closer) {
            left_summary->UpdateFrom(*moved[index]);
            if (!options.singly) {
                moved[index]->parent = left_summary_ptr;
            }
            left_node_ptr->AppendSubcluster(std::move(moved[index]));
        } else {
            right_summary->UpdateFrom(*moved[index]);
            if (!options.singly) {
                moved[index]->parent = right_summary_ptr;
            }
            right_node_ptr->AppendSubcluster(std::move(moved[index]));
        }
    }

    return {std::move(left_summary), std::move(right_summary)};
}

}  // namespace OECluster::detail
