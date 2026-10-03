/**
 * @file test_umbrella_header.cpp
 * @brief The umbrella header must reach every public surface it advertises.
 *
 * One symbol is named per header included by ``include/oecluster/oecluster.h``,
 * always in a form that needs the definition rather than a declaration. The
 * limit of the technique is transitive reach: a header dropped from the
 * umbrella that some remaining umbrella header includes on its own -- GateFacts.h
 * through PairwiseComparison.h, for one -- still arrives, and nothing here
 * notices. What is caught is a surface that becomes unreachable through the
 * umbrella entirely.
 */

#include <gtest/gtest.h>
#include <stdexcept>
#include <oecluster/oecluster.h>

namespace {

/// Instantiate to force a complete type. A forward declaration does not satisfy
/// ``sizeof``, so this compiles only where the umbrella header carried the
/// definition in; the value it returns is incidental.
template <typename T>
size_t definition_size() {
    return sizeof(T);
}

}  // namespace

// No other OECluster include: naming these types is the assertion, and it is
// made when this file compiles rather than when the tests run.
TEST(UmbrellaHeaderTest, ReachesDescriptorStatistics) {
    const OECluster::DescriptorStatisticsOptions options;
    EXPECT_TRUE(options.sources.empty());
    EXPECT_FALSE(options.inverse_covariance);

    const OECluster::DescriptorStatisticsResult result;
    EXPECT_EQ(result.num_rows, 0u);
    EXPECT_TRUE(result.columns.empty());
}

TEST(UmbrellaHeaderTest, ReachesDescriptorComparison) {
    const OECluster::DescriptorOptions options;
    EXPECT_TRUE(options.sources.empty());
    EXPECT_EQ(options.metric, "standardized_euclidean");
    EXPECT_EQ(options.missing, "complete_case");

    // The type name is the assertion; cannot construct without molecules.
    const OECluster::DescriptorComparison* ptr = nullptr;
    EXPECT_EQ(ptr, nullptr);
}

TEST(UmbrellaHeaderTest, ReachesRMSDComparison) {
    const OECluster::RMSDOptions options;
    EXPECT_FALSE(options.overlay);
    EXPECT_TRUE(options.automorph);
    EXPECT_TRUE(options.heavy_only);

    // The type name is the assertion; cannot construct without molecules.
    const OECluster::RMSDComparison* ptr = nullptr;
    EXPECT_EQ(ptr, nullptr);
}

TEST(UmbrellaHeaderTest, ReachesMCSComparison) {
    const OECluster::MCSOptions options;
    EXPECT_EQ(options.search_mode, OECluster::MCSSearchMode::Approximate);
    EXPECT_EQ(options.match_level, OECluster::MCSMatchLevel::Default);
    EXPECT_EQ(options.max_matches, 1024u);
    EXPECT_FALSE(options.similarity);

    // The type name is the assertion; cannot construct without molecules.
    const OECluster::MCSComparison* ptr = nullptr;
    EXPECT_EQ(ptr, nullptr);
}

TEST(UmbrellaHeaderTest, ReachesTheCoreHeaders) {
    EXPECT_NE(definition_size<OECluster::ComparisonError>(), 0u);
    EXPECT_NE(definition_size<OECluster::GateFacts>(), 0u);
    EXPECT_NE(definition_size<OECluster::PairwiseComparison>(), 0u);
    EXPECT_NE(definition_size<OECluster::DenseStorage>(), 0u);
    EXPECT_NE(definition_size<OECluster::ThreadPool>(), 0u);
    EXPECT_NE(definition_size<OECluster::PDistOptions>(), 0u);
    EXPECT_NE(definition_size<OECluster::CDistOptions>(), 0u);
    EXPECT_NE(definition_size<OECluster::DistanceMatrix>(), 0u);

    // CondensedIndex.h publishes free functions rather than a type, so calling
    // one is what needs the header.
    EXPECT_EQ(OECluster::pair_to_condensed(0, 1, 3), 0u);
}

TEST(UmbrellaHeaderTest, ReachesTheRemainingComparisons) {
    EXPECT_NE(definition_size<OECluster::FingerprintOptions>(), 0u);
    EXPECT_NE(definition_size<OECluster::ROCSOptions>(), 0u);
    EXPECT_NE(definition_size<OECluster::SuperposeOptions>(), 0u);
}

TEST(UmbrellaHeaderTest, ReachesTheClusteringHeaders) {
    EXPECT_NE(definition_size<OECluster::ButinaOptions>(), 0u);
    EXPECT_NE(definition_size<OECluster::RepresentativeOptions>(), 0u);
    EXPECT_NE(definition_size<OECluster::DBSCANOptions>(), 0u);
    EXPECT_NE(definition_size<OECluster::HDBSCANOptions>(), 0u);
    EXPECT_NE(definition_size<OECluster::AgglomerativeOptions>(), 0u);
    EXPECT_NE(definition_size<OECluster::BitBirchOptions>(), 0u);
    EXPECT_NE(definition_size<OECluster::KMedoidsOptions>(), 0u);
    EXPECT_NE(definition_size<OECluster::MurckoOptions>(), 0u);
    EXPECT_NE(definition_size<OECluster::MurckoResult>(), 0u);

    // ClusterTypes.h is reached through the free function rather than through
    // ClusteringResult, which every header above would drag in anyway.
    EXPECT_TRUE(OECluster::labels_to_clusters({}).empty());
}

// ClusterReport.h was absent from the umbrella header until 5.3.0, so A1's
// entire public surface was unreachable through the entry point the
// documentation names. All three quality headers are covered here so the
// omission cannot recur for any of them.
TEST(UmbrellaHeaderTest, ReachesTheClusterQualityHeaders) {
    const OECluster::ClusterReportOptions report_options;
    EXPECT_TRUE(report_options.treat_noise_as_singletons);

    const OECluster::PartitionAgreementOptions agreement_options;
    EXPECT_EQ(agreement_options.noise_handling,
              OECluster::NoiseHandling::Singletons);
    EXPECT_FALSE(agreement_options.compute_adjusted_mutual_information);

    const OECluster::SARCoherenceOptions coherence_options;
    EXPECT_EQ(coherence_options.noise_handling,
              OECluster::NoiseHandling::Excluded);

    const OECluster::ActivityLandscapeOptions landscape_options;
    EXPECT_DOUBLE_EQ(landscape_options.distance_threshold, 0.30);
    EXPECT_DOUBLE_EQ(landscape_options.activity_threshold, 1.0);
    EXPECT_DOUBLE_EQ(landscape_options.rmodi_delta, 0.625);

    const OECluster::ModelabilityOptions model_options;
    EXPECT_EQ(model_options.num_threads, 0u);
}

// DiversitySelection.h joined the umbrella in 5.8.0; naming both options
// structs keeps it from silently dropping out, as ClusterReport.h once did.
TEST(UmbrellaHeaderTest, ReachesTheDiversitySelectionHeader) {
    const OECluster::MaxMinOptions maxmin_options;
    EXPECT_EQ(maxmin_options.count, 0u);
    EXPECT_EQ(maxmin_options.seed_mode, OECluster::MaxMinSeed::Index);
    EXPECT_EQ(maxmin_options.chunk_size, 256u);

    const OECluster::CirclesOptions circles_options;
    EXPECT_EQ(circles_options.method, OECluster::CirclesMethod::MaxMin);
}

// SetDiversity.h joined the umbrella in 5.9.0; naming both options structs
// keeps it from silently dropping out, as ClusterReport.h once did.
TEST(UmbrellaHeaderTest, ReachesTheSetDiversityHeader) {
    const OECluster::VendiOptions vendi_options;
    EXPECT_EQ(vendi_options.order, 1u);
    EXPECT_EQ(vendi_options.kernel, OECluster::DiversityKernel::Complement);
    EXPECT_EQ(vendi_options.max_exact, 2048u);

    const OECluster::LogDetOptions logdet_options;
    EXPECT_DOUBLE_EQ(logdet_options.ridge, 0.0);
}

// SphereExclusion.h joined the umbrella in 5.10.0.
TEST(UmbrellaHeaderTest, ReachesTheSphereExclusionHeader) {
    const OECluster::SphereExclusionOptions options;
    EXPECT_EQ(options.order, OECluster::SphereOrder::Input);
    EXPECT_EQ(options.assignment, OECluster::SphereAssignment::First);
    EXPECT_EQ(options.chunk_size, 4096u);
    EXPECT_EQ(OECluster::SphereExclusionResult().Method(), "sphere_exclusion");
}

// KNNGraph.h joined the umbrella in 5.11.0.
TEST(UmbrellaHeaderTest, ReachesTheKNNGraphHeader) {
    const OECluster::KNNGraphOptions options;
    EXPECT_EQ(options.k, 0u);
    EXPECT_EQ(options.num_threads, 0u);
    EXPECT_EQ(options.chunk_size, 4096u);
    EXPECT_EQ(OECluster::KNNGraph().NumItems(), 0u);
}

// JarvisPatrick.h joined the umbrella in 5.11.0.
TEST(UmbrellaHeaderTest, ReachesTheJarvisPatrickHeader) {
    const OECluster::JarvisPatrickOptions options;
    EXPECT_EQ(options.k, 0u);
    EXPECT_EQ(options.kmin, 0u);
    EXPECT_EQ(options.chunk_size, 4096u);
    EXPECT_EQ(OECluster::JarvisPatrickResult().Method(), "jarvis_patrick");
}

// Leiden.h joined the umbrella in 5.12.0.
TEST(UmbrellaHeaderTest, ReachesTheLeidenHeader) {
    const OECluster::LeidenOptions options;
    EXPECT_EQ(options.k, 0u);
    EXPECT_TRUE(options.objective == OECluster::LeidenObjective::Modularity);
    EXPECT_EQ(options.resolution, 1.0);
    EXPECT_EQ(options.n_iterations, -1);
    EXPECT_EQ(options.seed, 0u);
    EXPECT_EQ(options.chunk_size, 4096u);
    EXPECT_EQ(OECluster::LeidenResult().Method(), "leiden");
}

// ISimReport.h joined the umbrella in 5.14.0; naming its options structs keeps
// it from silently dropping out, as ClusterReport.h once did.
TEST(UmbrellaHeaderTest, ReachesTheISimReportHeader) {
    const OECluster::ISimOptions isim_options;
    EXPECT_EQ(isim_options.metric, "tanimoto");
    const OECluster::ISimReportOptions report_options;
    EXPECT_FALSE(report_options.compute_centroid_indices);
    EXPECT_NE(definition_size<OECluster::ISimReport>(), 0u);
}

// Subset.h joined the umbrella in 5.16.0; both gathers are called so the
// definitions, not only the declarations, must be reachable.
TEST(UmbrellaHeaderTest, ReachesTheSubsetHeader) {
    OECluster::DenseStorage source(2);
    source.Set(0, 1, 0.25);
    OECluster::DenseStorage destination(2);
    OECluster::take_pairs(source, {1, 0}, destination);
    EXPECT_DOUBLE_EQ(destination.Get(0, 1), 0.25);
    EXPECT_THROW(OECluster::take_fingerprints(OEFP::OEFPBatch(), {}),
                 std::invalid_argument);
}

// Consensus.h joined the umbrella in 5.17.0; all three kernels are called so
// the definitions, not only the declarations, must be reachable.
TEST(UmbrellaHeaderTest, ReachesTheConsensusHeader) {
    OECluster::DenseStorage destination(2);
    const OECluster::ConsensusMatrixSummary summary =
        OECluster::coassociation_distances(2, {0, 2}, {0, 1}, {0, 0},
                                           destination);
    EXPECT_EQ(summary.num_partitions, 1u);
    EXPECT_DOUBLE_EQ(destination.Get(0, 1), 0.0);

    const std::vector<int> labels =
        OECluster::consensus_components(destination, 0.5);
    EXPECT_EQ(labels, (std::vector<int>{0, 0}));

    const OECluster::ConsensusStrength strength =
        OECluster::consensus_strength(destination, labels);
    EXPECT_DOUBLE_EQ(strength.item_consensus[0], 1.0);
}
