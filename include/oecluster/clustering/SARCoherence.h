/**
 * @file SARCoherence.h
 * @brief Whether structural neighbourhood predicts activity.
 *
 * Three entry points, split by the shape of their input rather than by the
 * metric they report. ``sar_coherence`` takes cluster labels and a continuous
 * activity. ``activity_landscape`` and ``modelability`` take pairwise
 * distances, from a precomputed matrix or from a comparison evaluated lazily
 * -- the first with a continuous activity, the second with class annotations
 * -- because a neighbourhood question cannot be answered from labels alone.
 */

#ifndef OECLUSTER_CLUSTERING_SARCOHERENCE_H
#define OECLUSTER_CLUSTERING_SARCOHERENCE_H

#include <cstddef>
#include <limits>
#include <string>
#include <vector>

#include "oecluster/PairwiseComparison.h"
#include "oecluster/StorageBackend.h"
#include "oecluster/clustering/ClusterTypes.h"
#include "oecluster/clustering/PartitionAgreement.h"

namespace OECluster {

/**
 * @brief Per-cluster activity summary.
 *
 * One record per cluster that retained at least one scored member. Present
 * unconditionally: the aggregate effect size is uninterpretable without the
 * group means behind it.
 */
struct ClusterActivity {
    /// The cluster's label, as it appears in the input. Under
    /// NoiseHandling::Grouped the merged noise row reports -1; under
    /// NoiseHandling::Singletons each noise sample keeps its own label, so
    /// several rows may share a label and are distinguished by position.
    ClusterLabel label = 0;
    /// Members with a finite activity value.
    std::size_t num_scored = 0;
    /// Mean activity over the scored members.
    double mean_activity = std::numeric_limits<double>::quiet_NaN();
    /// Population standard deviation over the scored members. NaN when
    /// num_scored < 2, where a spread is not defined.
    double stddev_activity = std::numeric_limits<double>::quiet_NaN();
};

/** @brief Options for sar_coherence. */
struct SARCoherenceOptions {
    /// How negatively-labelled samples are treated. Defaults to Excluded,
    /// unlike PartitionAgreementOptions: a noise point promoted to a singleton
    /// cluster has a group mean equal to its own value and so reads as
    /// perfectly explained variance.
    NoiseHandling noise_handling = NoiseHandling::Excluded;
};

/**
 * @brief How much of the activity variance the clustering explains.
 *
 * Reports the raw effect size and its chance-corrected counterpart. Compare
 * clusterings with different cluster counts on omega_squared; eta_squared
 * rises with the cluster count even when activity is independent of the
 * labels, by roughly (K - 1) / (n - 1), where n is num_scored and K is
 * num_clusters. Both are post-admission counts: they are taken after noise
 * handling and missing activities have removed samples, so the inflation term
 * is the one for the data that actually entered the decomposition, not for the
 * input as supplied.
 *
 * Use this when the question is whether a clustering groups molecules that
 * behave alike. It needs no distance matrix, so it answers that question for a
 * labelling of any provenance, including one that did not come from this
 * library.
 *
 * No p-value is reported. The one-way ANOVA F test behind these effect sizes
 * assumes the groups were fixed in advance, which is exactly what a clustering
 * violates.
 */
struct SARCoherence {
    /// Length of the input activity vector, always, whatever noise_handling
    /// is. This differs from PartitionAgreement::num_samples, which reports the
    /// samples that entered the table and so shrinks under Excluded. Read
    /// num_scored for the surviving count.
    std::size_t num_samples = 0;
    /// Samples with a finite activity that survived noise handling.
    std::size_t num_scored = 0;
    /// Clusters holding at least one scored sample.
    std::size_t num_clusters = 0;
    /// SS_between / SS_total, in [0, 1]. NaN when fewer than two samples were
    /// scored or when the activity has no variance.
    double eta_squared = std::numeric_limits<double>::quiet_NaN();
    /// Chance-corrected effect size. May be negative, which means the
    /// separation is weaker than chance would produce. Not clamped. NaN in the
    /// eta_squared cases and when every scored sample is its own cluster.
    double omega_squared = std::numeric_limits<double>::quiet_NaN();
    /// Ordered by first appearance of the label among the scored samples, not
    /// among the raw input: a cluster whose earliest member was dropped takes
    /// its position from its earliest surviving member instead.
    std::vector<ClusterActivity> clusters;
};

/**
 * @brief Activity variance explained by a clustering result.
 *
 * :param result: Any clustering result; only its labels are read.
 * :param activity: One value per sample, NaN for a missing measurement.
 * :param options: See SARCoherenceOptions.
 * :returns: The effect sizes and the per-cluster table.
 * :raises std::invalid_argument: If activity is empty, its length differs from
 *     the label count, it holds an infinity, or its magnitudes are large
 *     enough that a mean or a sum of squares overflows.
 */
SARCoherence sar_coherence(const ClusteringResult& result,
                           const std::vector<double>& activity,
                           const SARCoherenceOptions& options =
                               SARCoherenceOptions());

/**
 * @brief Activity variance explained by a labeling.
 *
 * :param labels: Per-sample labels; negative values are noise.
 * :param activity: One value per sample, NaN for a missing measurement.
 * :param options: See SARCoherenceOptions.
 * :returns: The effect sizes and the per-cluster table.
 * :raises std::invalid_argument: As the ClusteringResult overload.
 */
SARCoherence sar_coherence(const std::vector<ClusterLabel>& labels,
                           const std::vector<double>& activity,
                           const SARCoherenceOptions& options =
                               SARCoherenceOptions());

/** @brief Options for activity_landscape. */
struct ActivityLandscapeOptions {
    /// Pairs at or below this distance are structurally near. Shares
    /// ClusterReportOptions::boundary_threshold's default because it encodes
    /// the same judgement.
    double distance_threshold = 0.30;
    /// Activity differences at or above this are sharp. One log unit.
    double activity_threshold = 1.0;
    /// RMODI band half-width, in activity standard deviations. The published
    /// default from Ruiz and Gomez-Nieto (2018).
    double rmodi_delta = 0.625;
    /// 0 selects the hardware concurrency. Results do not depend on this value.
    std::size_t num_threads = 0;
    /// Pairwise distances per work unit on the comparison overload, at least
    /// one; ignored by the storage overload. A unit is whole rows and never
    /// fewer than 64, the floor the storage overload also applies, because
    /// each unit merges two length-n buffers. Results do not depend on this
    /// value.
    std::size_t chunk_size = 4096;
};

/**
 * @brief Continuous activity against descriptor-space distance.
 *
 * Cliff density and SALI describe how abruptly activity changes between near
 * neighbours. RMODI describes how often a molecule's nearest neighbour shares
 * its activity band.
 *
 * Use this when the question is whether the descriptor space is smooth enough
 * for a regression model to learn, or where its cliffs are. Every reported
 * value is independent of num_threads, bit for bit. That is thread-count
 * invariance, not permutation invariance: reordering the samples reorders the
 * summation and may move the last digit.
 */
struct ActivityLandscape {
    /// Rows in the distance matrix, or the comparison's Size().
    std::size_t num_samples = 0;
    /// Samples with a finite activity value.
    std::size_t num_scored = 0;
    /// Pairs with both endpoints scored: num_scored * (num_scored - 1) / 2.
    std::size_t num_pairs_scored = 0;
    /// Scored pairs that are both near and sharply different.
    std::size_t num_cliffs = 0;
    /// num_cliffs / num_pairs_scored. NaN when num_pairs_scored is 0.
    double cliff_density = std::numeric_limits<double>::quiet_NaN();
    /// Scored pairs at distance exactly 0, which SALI cannot score. Excluded
    /// from max_sali and mean_sali and reported here so the exclusion is
    /// visible.
    std::size_t num_zero_distance_pairs = 0;
    /// Largest |activity difference| / distance over scored pairs with a
    /// nonzero distance. NaN when there are none.
    double max_sali = std::numeric_limits<double>::quiet_NaN();
    /// Mean of the same quantity over the same pairs. NaN when there are none.
    double mean_sali = std::numeric_limits<double>::quiet_NaN();
    /// Fraction of scored molecules whose nearest same-band neighbour is
    /// strictly closer than their nearest different-band neighbour. NaN when
    /// num_scored < 2.
    ///
    /// Ruiz and Gomez-Nieto (2018) define this as a ratio of the two nearest
    /// distances thresholded at one; comparing the two distances directly is
    /// equivalent wherever both are finite, and avoids forming a quotient that
    /// the zero-distance case cannot express. A molecule with no neighbour on
    /// one side has that side's minimum left at +infinity, so the comparison
    /// still decides: this is A3's extension, since the publication does not
    /// cover the case.
    double rmodi = std::numeric_limits<double>::quiet_NaN();
    /// Population standard deviation of the scored activities: the sigma that
    /// rmodi_delta multiplies. NaN when num_scored < 2.
    double activity_stddev = std::numeric_limits<double>::quiet_NaN();
};

/**
 * @brief Cliff density, SALI and RMODI over a precomputed distance matrix.
 *
 * :param storage: Complete pairwise distances. SparseStorage is rejected.
 * :param activity: One value per sample, NaN for a missing measurement.
 * :param options: See ActivityLandscapeOptions.
 * :returns: The landscape summary.
 * :raises std::invalid_argument: If storage is incomplete, activity is empty
 *     or mismatched or holds an infinity, an option is negative or non-finite,
 *     a distance is negative or non-finite, or an accumulator overflows. The
 *     distance check is part of the pair sweep, which is skipped when fewer
 *     than two samples are scored, so below that a corrupt matrix is reported
 *     as undefined metrics rather than refused.
 */
ActivityLandscape activity_landscape(const StorageBackend& storage,
                                     const std::vector<double>& activity,
                                     const ActivityLandscapeOptions& options =
                                         ActivityLandscapeOptions());

/**
 * @brief Cliff density, SALI and RMODI over a comparison, evaluated lazily.
 *
 * Bit-identical to the storage overload over a matrix holding
 * Compare(min(i, j), max(i, j)) for every pair. Each pair of scored samples is
 * compared once, on worker threads, and memory does not grow with the number
 * of pairs. Only pairs of scored samples are compared; a non-finite value on
 * a pair with a missing activity is never seen.
 *
 * Precondition: Compare(i, j) returns the same value on every call and every
 * clone, and is called only with i < j.
 *
 * :param comparison: Distance comparison; cloned once per running unit.
 * :param activity: One value per sample, NaN for a missing measurement.
 * :param options: See ActivityLandscapeOptions; chunk_size sizes the units.
 * :returns: The landscape summary.
 * :raises ComparisonError: If the comparison reports similarities, a nonzero
 *     self-distance, possibly non-finite values, or per-pair feature subsets;
 *     checked after chunk_size and before any comparison runs.
 * :raises std::invalid_argument: If chunk_size is zero (checked first), on
 *     the storage overload's activity and option refusals (checked after the
 *     facts, against comparison.Size()), or on its distance and overflow
 *     refusals during the sweep. With more than one bad distance, which one is
 *     reported depends on the thread schedule, as on the storage overload.
 */
ActivityLandscape activity_landscape(PairwiseComparison& comparison,
                                     const std::vector<double>& activity,
                                     const ActivityLandscapeOptions& options =
                                         ActivityLandscapeOptions());

/** @brief Per-class nearest-neighbour concordance. */
struct ClassConcordance {
    /// The class string, as it appears in the input.
    std::string label;
    /// Members of this class with a non-empty class string.
    std::size_t num_members = 0;
    /// Fraction of those members whose nearest scored neighbour shares the
    /// class. NaN when the class is the only scored class, where no molecule
    /// has a neighbour that could differ.
    double fraction_same_class = std::numeric_limits<double>::quiet_NaN();
};

/** @brief Options for modelability. */
struct ModelabilityOptions {
    /// 0 selects the hardware concurrency. Results do not depend on this value.
    std::size_t num_threads = 0;
    /// Pairwise distances per work unit on the comparison overload, at least
    /// one; ignored by the storage overload. A unit is whole rows, at least
    /// one. Results do not depend on this value.
    std::size_t chunk_size = 4096;
};

/**
 * @brief The modelability index (Golbraikh et al., 2014), for K classes.
 *
 * The mean over classes of the fraction of members whose nearest neighbour
 * shares the class. Low values say a classifier is unlikely to learn this
 * dataset from this descriptor, whatever the model.
 *
 * Use this before fitting a classifier, to find out whether the descriptor
 * carries the signal at all. Nearest-neighbour ties resolve to the lowest
 * scored index, so the result is independent of num_threads, bit for bit. That
 * is thread-count invariance, not permutation invariance: the tie rule is
 * defined on the scored positions, so reordering the samples can pick a
 * different tied neighbour.
 */
struct Modelability {
    /// Length of the input class vector.
    std::size_t num_samples = 0;
    /// Samples with a non-empty class string.
    std::size_t num_scored = 0;
    /// Distinct non-empty class strings.
    std::size_t num_classes = 0;
    /// Unweighted mean of fraction_same_class over the classes. NaN when
    /// num_classes < 2.
    double modi = std::numeric_limits<double>::quiet_NaN();
    /// Ordered by first appearance of the class string among the scored
    /// samples. The only thing dropped here is the empty class string, which
    /// is no class at all, so this coincides with input order over the classes
    /// that exist.
    std::vector<ClassConcordance> classes;
};

/**
 * @brief Nearest-neighbour class concordance over a precomputed matrix.
 *
 * :param storage: Complete pairwise distances. SparseStorage is rejected.
 * :param activity_classes: One class string per sample; empty means missing.
 * :param options: See ModelabilityOptions.
 * :returns: MODI and the per-class table.
 * :raises std::invalid_argument: If storage is incomplete, activity_classes is
 *     empty or mismatched, or a distance is negative or non-finite. The
 *     distance check needs a pair to read, so it runs only when at least two
 *     samples are scored; below that a corrupt matrix is reported as undefined
 *     metrics rather than refused. One class with two or more scored samples
 *     is still checked.
 */
Modelability modelability(const StorageBackend& storage,
                          const std::vector<std::string>& activity_classes,
                          const ModelabilityOptions& options =
                              ModelabilityOptions());

/**
 * @brief Nearest-neighbour class concordance over a comparison, evaluated lazily.
 *
 * Bit-identical to the storage overload over a matrix holding
 * Compare(min(i, j), max(i, j)) for every pair. Each scored row is scanned in
 * full, as on the storage overload, so a run with n scored samples and at
 * least two scored classes makes n * (n - 1) comparisons: every pair twice,
 * both times as Compare(min, max). With one scored class nothing is scored
 * and each pair is only validated, n * (n - 1) / 2 comparisons.
 * Memory does not grow with the number of pairs. Only pairs of scored samples
 * are compared; a non-finite value on a pair with an empty class is never
 * seen.
 *
 * Precondition: Compare(i, j) returns the same value on every call and every
 * clone, and is called only with i < j.
 *
 * :param comparison: Distance comparison; cloned once per running unit.
 * :param activity_classes: One class string per sample; empty means missing.
 * :param options: See ModelabilityOptions; chunk_size sizes the units.
 * :returns: MODI and the per-class table.
 * :raises ComparisonError: If the comparison reports similarities, a nonzero
 *     self-distance, possibly non-finite values, or per-pair feature subsets;
 *     checked after chunk_size and before any comparison runs.
 * :raises std::invalid_argument: If chunk_size is zero (checked first), on
 *     the storage overload's class refusals (checked after the facts, against
 *     comparison.Size()), or on a negative or non-finite distance during the
 *     sweep.
 */
Modelability modelability(PairwiseComparison& comparison,
                          const std::vector<std::string>& activity_classes,
                          const ModelabilityOptions& options =
                              ModelabilityOptions());

}  // namespace OECluster

#endif  // OECLUSTER_CLUSTERING_SARCOHERENCE_H
