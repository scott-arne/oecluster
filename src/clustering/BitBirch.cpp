/**
 * @file BitBirch.cpp
 * @brief Public BitBirch clustering entry points.
 */

#include "oecluster/clustering/BitBirch.h"

#include <stdexcept>

#include "BitBirchTree.h"

namespace OECluster {
namespace {

void validate_options(const BitBirchOptions& options) {
    if (options.threshold < 0.0) {
        throw std::invalid_argument("BitBirch threshold must be non-negative");
    }
    if (options.branching_factor < 1) {
        throw std::invalid_argument("BitBirch branching_factor must be at least one");
    }
    if (options.tolerance < 0.0) {
        throw std::invalid_argument("BitBirch tolerance must be non-negative");
    }
}

void validate_batch(const OEFP::OEFPBatch& fingerprints) {
    if (fingerprints.Size() > 0 && fingerprints.SizeBits() == 0) {
        throw std::invalid_argument("BitBirch requires non-zero-width fingerprints");
    }
}

}  // namespace

BitBirchResult bitbirch_cluster(
    const OEFP::OEFPBatch& fingerprints,
    const BitBirchOptions& options) {
    validate_options(options);
    validate_batch(fingerprints);

    detail::BitBirchTree tree(options);
    tree.Fit(fingerprints);
    return tree.Result(fingerprints.Spec(), fingerprints.Size());
}

BitBirchResult bitbirch_recluster(
    const OEFP::OEFPBatch& fingerprints,
    const BitBirchReclusteringOptions& options) {
    BitBirchOptions first_options;
    first_options.threshold = options.initial_threshold;
    first_options.branching_factor = options.branching_factor;
    first_options.merge_criterion = BitBirchMergeCriterion::Diameter;
    first_options.mode = options.mode;
    first_options.num_threads = options.num_threads;
    validate_options(first_options);
    validate_batch(fingerprints);

    detail::BitBirchTree first_tree(first_options);
    first_tree.Fit(fingerprints);
    auto [rest, largest_singletons] =
        first_tree.PrepareReclusteringSubclusters(fingerprints);

    BitBirchOptions second_options;
    second_options.threshold = options.second_threshold;
    second_options.branching_factor = options.branching_factor;
    second_options.merge_criterion = BitBirchMergeCriterion::Tolerance;
    second_options.tolerance = options.second_tolerance;
    second_options.singly = true;
    second_options.mode = options.mode;
    second_options.num_threads = options.num_threads;
    validate_options(second_options);

    detail::BitBirchTree second_tree(second_options);
    second_tree.FitSubclusters(std::move(rest));
    second_tree.FitSubclusters(std::move(largest_singletons));
    return second_tree.Result(fingerprints.Spec(), fingerprints.Size());
}

BitBirchResult bitbirch_refine(
    const OEFP::OEFPBatch& fingerprints,
    const BitBirchRefinementOptions& options) {
    validate_options(options.fit_options);
    validate_batch(fingerprints);
    if (options.redistribute_largest_cluster && options.fit_options.singly) {
        throw std::invalid_argument(
            "BitBirch redistribute_largest_cluster requires singly=false");
    }
    if (options.reassign_top_clusters == 1) {
        throw std::invalid_argument(
            "BitBirch reassign_top_clusters must be zero or at least two");
    }

    detail::BitBirchTree tree(options.fit_options);
    tree.Fit(fingerprints);
    if (options.redistribute_largest_cluster) {
        tree.RedistributeLargestCluster(fingerprints);
    }
    if (options.reassign_top_clusters > 0) {
        tree.ReassignTopClusters(
            fingerprints,
            options.reassign_top_clusters,
            options.num_threads);
    }
    return tree.Result(fingerprints.Spec(), fingerprints.Size());
}

}  // namespace OECluster
