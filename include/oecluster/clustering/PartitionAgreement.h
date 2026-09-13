/**
 * @file PartitionAgreement.h
 * @brief Agreement between two labelings of the same samples.
 *
 * The input is labels and nothing else -- no distance matrix and no storage
 * backend, unlike ClusterReport.h. That is deliberate: a caller comparing two
 * clustering methods should not have to own a distance matrix to learn whether
 * the two found the same structure.
 */

#ifndef OECLUSTER_CLUSTERING_PARTITIONAGREEMENT_H
#define OECLUSTER_CLUSTERING_PARTITIONAGREEMENT_H

#include <cstddef>
#include <limits>
#include <string>
#include <vector>

#include "oecluster/clustering/ClusterTypes.h"

namespace OECluster {

/**
 * @brief How noise samples enter the contingency table.
 *
 * A sample is noise on a side when its label there is negative -- not only
 * NOISE_LABEL, matching how ClusterReport and labels_to_clusters already read
 * labels. In ``scaffold_agreement`` an empty scaffold string is noise on the
 * string side and follows the same three readings.
 */
enum class NoiseHandling {
    Singletons,  ///< Each noise sample becomes its own cluster. The default.
    Grouped,     ///< All of one side's noise samples form one cluster.
    Excluded     ///< Noise samples are dropped from both partitions.
};

/**
 * @brief Records the request, not the outcome.
 *
 * A requested metric whose value is undefined still reports true here, so NaN
 * reads unambiguously: false means nobody asked, true with NaN means asked and
 * undefined.
 */
struct PartitionAgreementRequested {
    bool adjusted_mutual_information = false;
};

/// Options for partition_agreement and scaffold_agreement.
struct PartitionAgreementOptions {
    /// How negatively-labelled samples enter the table. See NoiseHandling.
    NoiseHandling noise_handling = NoiseHandling::Singletons;

    /// Adjusted mutual information. Off by default: its expected-MI correction
    /// is the one term whose cost grows with the cluster count rather than the
    /// sample count.
    bool compute_adjusted_mutual_information = false;
};

/**
 * @brief The agreement scorecard.
 *
 * Every field is NaN when fewer than two samples survive noise handling, and
 * every field is 1.0 when the two partitions are identical after noise
 * handling; the per-field notes describe the remaining case. Pair counts are
 * accumulated in uint64_t, which is exact for any sample count below 2^32.
 */
struct PartitionAgreement {
    /// Samples entering the table. Equals the input length except under
    /// NoiseHandling::Excluded, which can shrink it.
    size_t num_samples = 0;
    /// Distinct clusters on each side after noise handling. Under Singletons
    /// these include one entry per noise sample.
    size_t num_clusters_a = 0;
    size_t num_clusters_b = 0;

    /// Hubert-Arabie adjusted Rand index. 1.0 is exact agreement, 0.0 is the
    /// value expected by chance, and negative values are worse than chance.
    /// Never NaN once two or more samples survive: the only inputs whose
    /// denominator vanishes are identical partitions, which report 1.0.
    double adjusted_rand_index = std::numeric_limits<double>::quiet_NaN();
    /// Geometric mean of pair precision and pair recall. Range [0, 1].
    /// NaN when either side is all singletons and the partitions differ, so
    /// that no coincident pair exists to normalize against; scikit-learn
    /// reports 0.0 there instead.
    double fowlkes_mallows = std::numeric_limits<double>::quiet_NaN();
    /// Mutual information over the arithmetic mean of the two entropies.
    /// Range [0, 1]. Assigned from the same computed value as v_measure, so
    /// the two are bitwise equal. Never NaN once two or more samples survive.
    double normalized_mutual_information =
        std::numeric_limits<double>::quiet_NaN();
    /// MI / H(a), the fraction of side A's information that side B explains.
    /// Asymmetric: homogeneity and completeness are the only fields that
    /// change when the arguments are swapped. NaN when H(a) is zero, that is
    /// when side A is a single cluster and the partitions differ;
    /// scikit-learn reports 1.0 there instead.
    double homogeneity = std::numeric_limits<double>::quiet_NaN();
    /// MI / H(b). Asymmetric; see homogeneity. NaN when H(b) is zero, that is
    /// when side B is a single cluster and the partitions differ;
    /// scikit-learn reports 1.0 there instead.
    double completeness = std::numeric_limits<double>::quiet_NaN();
    /// The beta = 1 V-measure, assigned from 2*MI/(H(a)+H(b)) -- the same
    /// value as normalized_mutual_information, and equal to the harmonic mean
    /// of homogeneity and completeness wherever that harmonic mean is
    /// defined. It stays defined where the harmonic form does not: a
    /// zero-entropy side gives v_measure == 0.0 while homogeneity or
    /// completeness is NaN. Both names are reported because both are in
    /// common use.
    double v_measure = std::numeric_limits<double>::quiet_NaN();
    /// Mutual information corrected for chance, normalized by the arithmetic
    /// mean of the two entropies. NaN unless
    /// PartitionAgreementOptions::compute_adjusted_mutual_information was set;
    /// consult requested to tell "not asked" from "asked and undefined". This
    /// is the one metric whose denominator is clamped rather than reported as
    /// NaN, matching scikit-learn.
    double adjusted_mutual_information =
        std::numeric_limits<double>::quiet_NaN();

    PartitionAgreementRequested requested;
};

/**
 * @brief Agreement between two labelings of the same samples.
 *
 * Reads Labels() only. Unlike cluster_report this never cross-checks Labels()
 * against Members(), so a result whose two views disagree is not rejected here.
 *
 * :param a: The reference labeling. Drives homogeneity.
 * :param b: The candidate labeling. Drives completeness.
 * :param options: Noise handling and the AMI opt-in.
 * :returns: The agreement scorecard.
 * :raises std::invalid_argument: If the labelings differ in length, or are
 *     empty.
 */
PartitionAgreement partition_agreement(
    const ClusteringResult& a, const ClusteringResult& b,
    const PartitionAgreementOptions& options = PartitionAgreementOptions());

/// Overload taking raw label vectors, for callers holding labels that never
/// came from a ClusteringResult.
PartitionAgreement partition_agreement(
    const std::vector<ClusterLabel>& a, const std::vector<ClusterLabel>& b,
    const PartitionAgreementOptions& options = PartitionAgreementOptions());

/**
 * @brief Agreement between a clustering and a per-sample scaffold annotation.
 *
 * The clustering is side A and the scaffold annotation is side B, following
 * the positional rule. So completeness carries the scaffold-purity reading --
 * whether each cluster's members share a single scaffold, the same question
 * RepresentativeMetrics::scaffold_purity asks per cluster -- and homogeneity
 * carries its transpose, whether each scaffold landed in a single cluster.
 * The symmetric metrics are unaffected.
 *
 * An empty scaffold string is missing data and follows options.noise_handling.
 *
 * :param result: The clustering to score.
 * :param scaffold_labels: One scaffold string per sample.
 * :param options: Noise handling and the AMI opt-in.
 * :returns: The agreement scorecard.
 * :raises std::invalid_argument: If scaffold_labels.size() does not equal
 *     result.NumSamples(), or either is empty.
 */
PartitionAgreement scaffold_agreement(
    const ClusteringResult& result,
    const std::vector<std::string>& scaffold_labels,
    const PartitionAgreementOptions& options = PartitionAgreementOptions());

/// Overload taking a raw label vector, the scaffold counterpart of the
/// raw-label partition_agreement. Side A is `labels`, side B is
/// `scaffold_labels`.
PartitionAgreement scaffold_agreement(
    const std::vector<ClusterLabel>& labels,
    const std::vector<std::string>& scaffold_labels,
    const PartitionAgreementOptions& options = PartitionAgreementOptions());

}  // namespace OECluster

#endif  // OECLUSTER_CLUSTERING_PARTITIONAGREEMENT_H
