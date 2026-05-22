/**
 * @file BitBirch.h
 * @brief BitBirch clustering over OEFP dense binary fingerprints.
 */

#ifndef OECLUSTER_CLUSTERING_BITBIRCH_H
#define OECLUSTER_CLUSTERING_BITBIRCH_H

#include <cstddef>
#include <vector>

#include "oefp/batch.h"
#include "oecluster/clustering/ClusterTypes.h"

namespace OECluster {

/**
 * @brief Merge criterion used when inserting BitBirch subclusters.
 */
enum class BitBirchMergeCriterion {
    Radius,
    Diameter,
    Tolerance,
    ToleranceTough
};

/**
 * @brief Execution mode for BitBirch.
 */
enum class BitBirchMode {
    StrictParity,
    Fast
};

/**
 * @brief Options for one-pass BitBirch clustering.
 */
struct BitBirchOptions {
    double threshold = 0.65;
    size_t branching_factor = 50;
    BitBirchMergeCriterion merge_criterion = BitBirchMergeCriterion::Diameter;
    double tolerance = 0.05;
    bool singly = true;
    BitBirchMode mode = BitBirchMode::StrictParity;
    size_t num_threads = 0;
};

/**
 * @brief Options for two-stage BitBirch reclustering.
 */
struct BitBirchReclusteringOptions {
    double initial_threshold = 0.65;
    double second_threshold = 0.7;
    double second_tolerance = 0.0;
    size_t branching_factor = 50;
    BitBirchMode mode = BitBirchMode::StrictParity;
    size_t num_threads = 0;
};

/**
 * @brief Options for refinement passes over a BitBirch fit.
 */
struct BitBirchRefinementOptions {
    BitBirchOptions fit_options;
    bool redistribute_largest_cluster = false;
    size_t reassign_top_clusters = 0;
    size_t num_threads = 0;
};

/**
 * @brief BitBirch clustering result with labels, clusters, and centroid fingerprints.
 */
struct BitBirchResult : public ClusteringResult {
    OEFP::OEFPBatch centroids;
    std::vector<size_t> cluster_sizes;
};

/**
 * @brief Cluster dense binary OEFP fingerprints with BitBirch.
 *
 * :param fingerprints: Dense binary fingerprint batch.
 * :param options: BitBirch clustering options.
 * :returns: Labels, clusters, centroids, and cluster sizes.
 * :raises std::invalid_argument: If options or fingerprints are invalid.
 */
BitBirchResult bitbirch_cluster(
    const OEFP::OEFPBatch& fingerprints,
    const BitBirchOptions& options = BitBirchOptions());

/**
 * @brief Cluster dense binary OEFP fingerprints with two-stage BitBirch reclustering.
 *
 * :param fingerprints: Dense binary fingerprint batch.
 * :param options: Reclustering options.
 * :returns: Labels, clusters, centroids, and cluster sizes.
 * :raises std::invalid_argument: If options or fingerprints are invalid.
 */
BitBirchResult bitbirch_recluster(
    const OEFP::OEFPBatch& fingerprints,
    const BitBirchReclusteringOptions& options = BitBirchReclusteringOptions());

/**
 * @brief Fit BitBirch and apply requested refinement passes.
 *
 * :param fingerprints: Dense binary fingerprint batch.
 * :param options: Refinement options.
 * :returns: Labels, clusters, centroids, and cluster sizes.
 * :raises std::invalid_argument: If options or fingerprints are invalid.
 */
BitBirchResult bitbirch_refine(
    const OEFP::OEFPBatch& fingerprints,
    const BitBirchRefinementOptions& options = BitBirchRefinementOptions());

}  // namespace OECluster

#endif  // OECLUSTER_CLUSTERING_BITBIRCH_H
