/**
 * @file ClusterReport.h
 * @brief Method-agnostic clustering-quality report and comparison.
 */

#ifndef OECLUSTER_CLUSTERING_CLUSTERREPORT_H
#define OECLUSTER_CLUSTERING_CLUSTERREPORT_H

#include <cstddef>
#include <limits>
#include <vector>

#include "oecluster/PairwiseComparison.h"
#include "oecluster/StorageBackend.h"
#include "oecluster/clustering/ClusterTypes.h"
#include "oecluster/clustering/Representative.h"

namespace OECluster {

/**
 * @brief Threshold preset for clustering-quality reports.
 *
 * Distances are Tanimoto/Jaccard (distance = 1 - similarity).
 */
enum class ClusterThreshold {
    Default,    ///< coverage {0.25,0.35,0.45}, boundary 0.30.
    Tight,      ///< coverage {0.20,0.30,0.40}, boundary 0.25.
    Diversity   ///< coverage {0.40,0.50,0.60}, boundary 0.40.
};

/**
 * @brief No other cluster exists, so ClusterRecord::nearest_cluster is undefined.
 *
 * Shares the value of NOISE_LABEL (ClusterTypes.h) deliberately -- both mean
 * "not a cluster" -- but the two are separate names because they answer
 * different questions. NOISE_LABEL says a sample belongs to no cluster;
 * NO_NEAREST_CLUSTER says a cluster has no neighbour to be near.
 */
constexpr ClusterLabel NO_NEAREST_CLUSTER = -1;

/**
 * @brief Which optional computations the caller asked for.
 *
 * Records the request, not the outcome. A requested metric whose value is
 * undefined still reports true here, so that NaN can be read unambiguously:
 * false means nobody asked, true with NaN means asked and undefined.
 */
struct ClusterReportRequested {
    bool pair_rank_indices = false;    ///< compute_pair_rank_indices was set.
    bool per_cluster_records = false;  ///< compute_per_cluster_records was set.
};

/**
 * @brief One row of the per-cluster table, in ClusteringResult::Members() order.
 *
 * label and nearest_cluster are ordinals into Members(), not values read out of
 * Labels(). The two are the same integer under cluster_report's partition
 * precondition; stating which one is authoritative removes the ambiguity rather
 * than relying on the coincidence.
 */
struct ClusterRecord {
    ClusterLabel label = 0;     ///< Ordinal in Members(); equals the Labels() value.
    size_t size = 0;            ///< Member count.
    size_t representative = 0;  ///< Sample index of the CONFIGURED representative.
    /// Mean within-pair distance; NaN for a singleton.
    double mean_intra_distance = std::numeric_limits<double>::quiet_NaN();
    /// Median within-pair distance; NaN for a singleton.
    double median_intra_distance = std::numeric_limits<double>::quiet_NaN();
    double radius = 0.0;    ///< Max representative-to-member distance; 0.0 for a singleton.
    double diameter = 0.0;  ///< Max pairwise distance within the cluster; 0.0 for a singleton.
    /// Mean distance from the representative to the OTHER n_k - 1 members. This
    /// denominator is n_k - 1, not the n_k the medoid scatter of the
    /// Davies-Bouldin and medoid-Dunn indices uses; the two quantities are
    /// different by design and differ by a factor of two at n_k == 2.
    double mean_representative_distance = 0.0;
    /// Ordinal in Members(), same space as label; NO_NEAREST_CLUSTER when K < 2.
    /// Ties resolve to the lowest ordinal.
    ClusterLabel nearest_cluster = NO_NEAREST_CLUSTER;
    /// Min single-linkage distance to nearest_cluster; NaN when K < 2. Not 0.0,
    /// which would read as another cluster sitting at zero distance.
    double nearest_cluster_distance = std::numeric_limits<double>::quiet_NaN();
    /// Mean over this cluster's members; NaN when K < 2. A singleton cluster
    /// scores 0 by Rousseeuw's convention: its member has no own-cluster pair,
    /// so the term is undefined rather than perfect.
    double silhouette = std::numeric_limits<double>::quiet_NaN();
    /// Pairs involving this cluster within boundary_threshold. Each violating
    /// pair is counted by both of its endpoints, so the sum over records is
    /// twice ClusterReport::boundary_violations, which counts each pair once.
    size_t boundary_violations = 0;
};

/**
 * @brief Options controlling clustering-quality report computation.
 */
struct ClusterReportOptions {
    /// Distances at which coverage_at is evaluated. A NaN entry is refused;
    /// infinity is accepted and covers every sample. See cluster_report for the
    /// partitions the NaN refusal reaches.
    std::vector<double> coverage_thresholds;
    /// Cross-cluster pairs at or below this distance are boundary violations.
    /// NaN is refused; infinity is accepted and counts every cross pair. See
    /// cluster_report for the partitions the NaN refusal reaches.
    double boundary_threshold = 0.30;
    RepresentativeMethod representative_method = RepresentativeMethod::Medoid;
    bool treat_noise_as_singletons = true;
    size_t num_threads = 0;
    /// Pairwise distances per work unit on the comparison overload, at least
    /// one; ignored by the storage overload. A unit is whole rows, and a fill
    /// block never holds more than 2^20 distances (or one row, if longer).
    size_t chunk_size = 4096;
    /// Enables c_index and baker_hubert_gamma. The two sorted arrays the
    /// indices are read off hold every pairwise distance among clustered
    /// points, Nc(Nc-1)/2 doubles in total -- roughly 400 MB at Nc = 10,000
    /// and 10 GB at Nc = 50,000. Both are built only under this flag; without
    /// it the report holds no pair-sized array.
    /// Storage overload only: the comparison overload refuses it.
    bool compute_pair_rank_indices = false;
    /// Enables ClusterReport::records. Off by default. A cluster with at most
    /// 2^20 pairs has its median taken from a buffer of its distances; a
    /// larger cluster's median is selected by re-walking its pairs, so the
    /// stage buffers at most 2^20 doubles (8 MB). Under
    /// compute_pair_rank_indices, which already holds every pair, each
    /// cluster is buffered whole.
    bool compute_per_cluster_records = false;

    /** @brief Seed coverage_thresholds and boundary_threshold from a preset. */
    explicit ClusterReportOptions(ClusterThreshold preset = ClusterThreshold::Default);
};

/**
 * @brief Clustering-quality scorecard. NaN marks an undefined metric.
 */
struct ClusterReport {
    // Basic profile.
    size_t num_samples = 0;
    size_t num_clusters = 0;
    size_t num_noise = 0;
    size_t num_singletons = 0;
    double noise_fraction = 0.0;
    double singleton_fraction = 0.0;
    double largest_cluster_fraction = 0.0;
    double cluster_size_median = 0.0;
    double cluster_size_p90 = 0.0;
    double size_gini = 0.0;
    double size_entropy = 0.0;
    // Compactness / separation.
    double mean_intra_distance = 0.0;
    double median_intra_distance = 0.0;
    double median_radius = 0.0;
    double p95_diameter = 0.0;
    /// Mean over every clustered point, singletons contributing 0 each
    /// (see ClusterRecord::silhouette). Matches Rousseeuw's definition and
    /// scikit-learn's silhouette_score.
    double silhouette = 0.0;
    double dunn_index = 0.0;
    size_t boundary_violations = 0;
    // Representative / coverage.
    double median_medoid_member_distance = 0.0;
    double representative_redundancy = 0.0;
    std::vector<double> coverage_thresholds;
    std::vector<double> coverage_at;  ///< Coverage fraction at each coverage_thresholds entry (parallel-indexed).

    // Internal indices, always computed. Exactly three fields below --
    // calinski_harabasz_medoid, davies_bouldin_medoid and
    // dunn_medoid_separation_medoid_spread -- use the true medoid and ignore
    // representative_method. The pre-existing median_medoid_member_distance is
    // medoid-named but keeps its configured-representative meaning; see the
    // class comment on cluster_report below.
    //
    // These seven default to NaN, unlike the 0.0 of every field above them in
    // this struct. The split is deliberate and is not to be tidied away: the
    // fields above predate this branch and ship in 5.0.0 with 0.0 defaults,
    // and changing an established default is a behaviour change nobody asked
    // for. New fields get the honest default -- an unpopulated metric reads as
    // undefined rather than as a measurement of zero.
    /// Higher is better. NaN when K < 2, Nc == K, or the denominator is zero.
    /// Medoid-substituted: the published Calinski-Harabasz definition uses
    /// centroids, which do not exist for a distance matrix, so each cluster's
    /// medoid stands in for its centroid and the global medoid stands in for
    /// the grand mean. The value is therefore not comparable with published
    /// figures or with scikit-learn's.
    double calinski_harabasz_medoid = std::numeric_limits<double>::quiet_NaN();
    /// Lower is better. NaN when K < 2; inf when two medoids coincide.
    /// Medoid-substituted: the published Davies-Bouldin definition uses
    /// centroids, which do not exist for a distance matrix, so each cluster's
    /// medoid stands in for its centroid. The value is therefore not
    /// comparable with published figures or with scikit-learn's.
    double davies_bouldin_medoid = std::numeric_limits<double>::quiet_NaN();
    /// Higher is better. NaN when K < 2 or the max mean within-pair distance is 0.
    double dunn_mean_separation_mean_diameter = std::numeric_limits<double>::quiet_NaN();
    /// Higher is better. NaN when K < 2 or the max medoid spread is 0.
    double dunn_medoid_separation_medoid_spread = std::numeric_limits<double>::quiet_NaN();
    /// Higher is better: positive means between-cluster distances exceed
    /// within-cluster ones. The sign convention is stated because published
    /// sources differ on it. NaN when there are no within-pairs, no
    /// between-pairs, or zero distance spread.
    double point_biserial = std::numeric_limits<double>::quiet_NaN();

    // Pair-rank indices, computed only under compute_pair_rank_indices. Read
    // `requested` to tell "nobody asked" apart from "asked and undefined".
    /// Lower is better. NaN when not requested, or when S_max == S_min.
    double c_index = std::numeric_limits<double>::quiet_NaN();
    /// Higher is better. NaN when not requested, or when s+ + s- == 0.
    double baker_hubert_gamma = std::numeric_limits<double>::quiet_NaN();

    /// Coverage over noise points only, against the configured representatives.
    /// Same length as coverage_at: coverage_thresholds.size() when K >= 1 and 0
    /// when K == 0. Every entry is NaN when the clustering has no noise.
    std::vector<double> noise_coverage_at;

    std::vector<ClusterRecord> records;  ///< Empty unless compute_per_cluster_records.
    ClusterReportRequested requested;    ///< What the caller asked for, not what was defined.
};

/**
 * @brief Two reports aligned for side-by-side comparison.
 */
struct ClusterReportComparison {
    ClusterReport a;
    ClusterReport b;
};

/**
 * @brief Compute a clustering-quality report.
 *
 * :param result: Any algorithm's clustering result (base ClusteringResult).
 * :param storage: Complete pairwise distance storage.
 * :param options: Report options (thresholds, representative method, flags).
 * :returns: A ClusterReport scorecard.
 * :raises std::invalid_argument: If storage cannot provide complete distances;
 *     if a cluster in result is empty or repeats a member; if result has at
 *     least one cluster and representative_method is HighestNeighborhood, whose
 *     neighbor threshold ClusterReportOptions has no field to supply; if result
 *     has at least one cluster and boundary_threshold or a coverage_thresholds
 *     entry is NaN, against which every comparison fails silently; if
 *     Labels() and Members() do not describe the same partition -- a sample in
 *     two clusters, a member whose Labels() entry disagrees with the cluster
 *     holding it, a clustered sample omitted from every member list, a
 *     non-noise label naming no cluster, or a noise-labelled sample sitting
 *     inside a cluster; or if any distance that reaches a reported value is
 *     not finite.
 * :raises std::out_of_range: If a cluster member is at or beyond
 *     storage.NumSamples(); if a cluster member is at or beyond Labels().size();
 *     or if result has at least one cluster and labels more samples than storage
 *     holds.
 */
ClusterReport cluster_report(
    const ClusteringResult& result,
    const StorageBackend& storage,
    const ClusterReportOptions& options = ClusterReportOptions());

/**
 * @brief Compute a clustering-quality report over a comparison, evaluated lazily.
 *
 * The result is bit-identical to the storage overload over a matrix holding
 * Compare(min(i, j), max(i, j)) for every pair, except that a zero median is
 * always +0.0. Memory does not grow with the number of pairs: distances are
 * computed in bounded blocks on worker threads and consumed in a fixed order,
 * and medians are selected without storing every pair. Only the pairs a
 * reported value reads are compared; a non-finite value elsewhere is never
 * seen.
 *
 * Precondition: Compare(i, j) returns the same value on every call and every
 * clone, and is called only with i < j.
 *
 * :param result: Any algorithm's clustering result (base ClusteringResult).
 * :param comparison: Distance comparison; cloned once per running unit.
 * :param options: Report options; num_threads and chunk_size drive the fill.
 * :returns: A ClusterReport scorecard.
 * :raises ComparisonError: If the comparison reports similarities, a nonzero
 *     self-distance, or possibly non-finite values; checked after chunk_size
 *     and before any comparison runs.
 * :raises std::invalid_argument: If chunk_size is zero (checked first) or
 *     compute_pair_rank_indices is set (checked after the facts), both before
 *     any comparison runs; otherwise on the storage overload's refusals, with
 *     a non-finite distance reported at the earliest pair in reading order
 *     whatever the thread count.
 * :raises std::out_of_range: As the storage overload, against
 *     comparison.Size().
 */
ClusterReport cluster_report(
    const ClusteringResult& result,
    PairwiseComparison& comparison,
    const ClusterReportOptions& options = ClusterReportOptions());

/**
 * @brief Pair two reports for side-by-side reading. No agreement math.
 *
 * :param a: First clustering-quality report.
 * :param b: Second clustering-quality report.
 * :returns: A comparison structure containing both reports.
 */
ClusterReportComparison compare_reports(const ClusterReport& a, const ClusterReport& b);

}  // namespace OECluster

#endif  // OECLUSTER_CLUSTERING_CLUSTERREPORT_H
