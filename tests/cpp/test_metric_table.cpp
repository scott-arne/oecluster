#include <gtest/gtest.h>
#include <oefp/oefp.h>
#include <limits>
#include <string>
#include <utility>
#include <vector>
#include "../../src/comparisons/MetricTable.h"
#include "oecluster/Error.h"

using namespace OECluster;

// The 17 metrics probed against OEFP 0.3.0. These rows are the oracle: if OEFP
// changes a capability, this test must fail loudly rather than let the gate
// silently start admitting or refusing different matrices.
TEST(MetricTableTest, OEFPCapabilityTableIsUnchanged) {
    const std::vector<double> variances{1.0, 2.0, 3.0};
    const std::vector<double> inverse_covariance{1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};

    struct Row {
        const char* label;
        OEFP::Metric metric;
        bool zero_self;
        bool triangle;
    };

    const std::vector<Row> rows{
        {"jaccard", OEFP::Metric::Jaccard(), true, true},
        {"tanimoto", OEFP::Metric::Tanimoto(), false, false},
        {"dice", OEFP::Metric::Dice(), true, false},
        {"manhattan", OEFP::Metric::Manhattan(), true, true},
        {"euclidean", OEFP::Metric::Euclidean(), true, true},
        {"bray_curtis", OEFP::Metric::BrayCurtis(), true, false},
        {"canberra", OEFP::Metric::Canberra(), true, true},
        {"hamming", OEFP::Metric::Hamming(), true, true},
        {"sokal_sneath", OEFP::Metric::SokalSneath(), true, true},
        {"matching", OEFP::Metric::Matching(), true, true},
        {"rogers_tanimoto", OEFP::Metric::RogersTanimoto(), true, true},
        {"russell_rao", OEFP::Metric::RussellRao(), false, true},
        {"kulsinski", OEFP::Metric::Kulsinski(), false, true},
        {"sokal_michener", OEFP::Metric::SokalMichener(), true, true},
        {"seuclidean", OEFP::Metric::StandardizedEuclidean(variances), true, true},
        {"mahalanobis", OEFP::Metric::Mahalanobis(inverse_covariance), true, true},
        {"chebyshev", OEFP::Metric::Chebyshev(), true, true},
    };

    ASSERT_EQ(rows.size(), 17u);
    for (const Row& row : rows) {
        EXPECT_EQ(row.metric.HasZeroSelfDistance(), row.zero_self) << row.label;
        EXPECT_EQ(row.metric.SatisfiesTriangleInequality(), row.triangle) << row.label;
    }
}

TEST(MetricTableTest, DefaultFingerprintPathResolvesToJaccard) {
    MetricParams params;
    const OEFP::Metric metric = resolve_metric("tanimoto", false, params, MetricSurface::Fingerprint);
    EXPECT_EQ(metric.Name(), OEFP::MetricName::Jaccard);
    EXPECT_TRUE(metric.HasZeroSelfDistance());
    EXPECT_TRUE(metric.SatisfiesTriangleInequality());
}

TEST(MetricTableTest, TanimotoSimilarityResolvesToTanimoto) {
    MetricParams params;
    const OEFP::Metric metric = resolve_metric("tanimoto", true, params, MetricSurface::Fingerprint);
    EXPECT_EQ(metric.Name(), OEFP::MetricName::Tanimoto);
}

// Only ASCII case is folded. Camel case such as "BrayCurtis" is not a
// supported spelling and resolves as an unknown metric.
TEST(MetricTableTest, UpperCaseSnakeNameIsAccepted) {
    MetricParams params;
    EXPECT_EQ(resolve_metric("BRAY_CURTIS", false, params, MetricSurface::Fingerprint).Name(),
              OEFP::MetricName::BrayCurtis);
}

TEST(MetricTableTest, EuclideanIsNowSupportedOnFingerprints) {
    MetricParams params;
    EXPECT_NO_THROW(resolve_metric("euclidean", false, params, MetricSurface::Fingerprint));
}

TEST(MetricTableTest, EveryTableRowResolvesToItsOwnMetric) {
    // One shared params object with every field populated to a valid value, so a
    // single params satisfies all 29 rows.
    MetricParams params;
    params.p = 2.0;
    params.tversky_alpha = 0.5;
    params.tversky_beta = 0.5;
    params.variances = {1.0, 2.0};
    params.inverse_covariance = {1.0, 0.0, 0.0, 1.0};

    struct Row {
        const char* name;
        MetricSurface surface;
        bool similarity;
        OEFP::MetricName expected;
    };

    const std::vector<Row> rows{
        // Fingerprint surface, similarity=false (17 rows)
        {"jaccard", MetricSurface::Fingerprint, false, OEFP::MetricName::Jaccard},
        {"tanimoto", MetricSurface::Fingerprint, false, OEFP::MetricName::Jaccard},
        {"dice", MetricSurface::Fingerprint, false, OEFP::MetricName::Dice},
        {"sokal_sneath", MetricSurface::Fingerprint, false, OEFP::MetricName::SokalSneath},
        {"matching", MetricSurface::Fingerprint, false, OEFP::MetricName::Matching},
        {"rogers_tanimoto", MetricSurface::Fingerprint, false, OEFP::MetricName::RogersTanimoto},
        {"russell_rao", MetricSurface::Fingerprint, false, OEFP::MetricName::RussellRao},
        {"kulsinski", MetricSurface::Fingerprint, false, OEFP::MetricName::Kulsinski},
        {"sokal_michener", MetricSurface::Fingerprint, false, OEFP::MetricName::SokalMichener},
        {"euclidean", MetricSurface::Fingerprint, false, OEFP::MetricName::Euclidean},
        {"manhattan", MetricSurface::Fingerprint, false, OEFP::MetricName::Manhattan},
        {"chebyshev", MetricSurface::Fingerprint, false, OEFP::MetricName::Chebyshev},
        {"hamming", MetricSurface::Fingerprint, false, OEFP::MetricName::Hamming},
        {"canberra", MetricSurface::Fingerprint, false, OEFP::MetricName::Canberra},
        {"bray_curtis", MetricSurface::Fingerprint, false, OEFP::MetricName::BrayCurtis},
        {"minkowski", MetricSurface::Fingerprint, false, OEFP::MetricName::Minkowski},
        {"tversky", MetricSurface::Fingerprint, false, OEFP::MetricName::Tversky},
        // Fingerprint surface, similarity=true (2 rows)
        {"tanimoto", MetricSurface::Fingerprint, true, OEFP::MetricName::Tanimoto},
        {"tversky", MetricSurface::Fingerprint, true, OEFP::MetricName::Tversky},
        // Descriptor surface, similarity=false (10 rows)
        {"euclidean", MetricSurface::Descriptor, false, OEFP::MetricName::Euclidean},
        {"manhattan", MetricSurface::Descriptor, false, OEFP::MetricName::Manhattan},
        {"chebyshev", MetricSurface::Descriptor, false, OEFP::MetricName::Chebyshev},
        {"hamming", MetricSurface::Descriptor, false, OEFP::MetricName::Hamming},
        {"canberra", MetricSurface::Descriptor, false, OEFP::MetricName::Canberra},
        {"bray_curtis", MetricSurface::Descriptor, false, OEFP::MetricName::BrayCurtis},
        {"minkowski", MetricSurface::Descriptor, false, OEFP::MetricName::Minkowski},
        {"standardized_euclidean", MetricSurface::Descriptor, false, OEFP::MetricName::StandardizedEuclidean},
        {"seuclidean", MetricSurface::Descriptor, false, OEFP::MetricName::StandardizedEuclidean},
        {"mahalanobis", MetricSurface::Descriptor, false, OEFP::MetricName::Mahalanobis},
    };

    ASSERT_EQ(rows.size(), 29u);
    for (const Row& row : rows) {
        SCOPED_TRACE(std::string(row.name) + " on " +
                     (row.surface == MetricSurface::Fingerprint ? "Fingerprint" : "Descriptor"));
        const OEFP::Metric metric = resolve_metric(row.name, row.similarity, params, row.surface);
        EXPECT_EQ(metric.Name(), row.expected);
    }
}

TEST(MetricTableTest, MinkowskiUsesTheSuppliedExponent) {
    MetricParams params;
    params.p = 0.5;
    const OEFP::Metric metric = resolve_metric("minkowski", false, params, MetricSurface::Descriptor);
    EXPECT_EQ(metric.Name(), OEFP::MetricName::Minkowski);
    EXPECT_DOUBLE_EQ(metric.P(), 0.5);
    EXPECT_FALSE(metric.SatisfiesTriangleInequality());
}

TEST(MetricTableTest, MinkowskiAboveOneSatisfiesTheTriangleInequality) {
    MetricParams params;
    params.p = 3.0;
    const OEFP::Metric metric = resolve_metric("minkowski", false, params, MetricSurface::Descriptor);
    EXPECT_DOUBLE_EQ(metric.P(), 3.0);
    EXPECT_TRUE(metric.SatisfiesTriangleInequality());
}

TEST(MetricTableTest, MinkowskiRejectsNonPositiveOrNaNExponent) {
    const std::vector<double> invalid_exponents{
        0.0,
        -1.0,
        std::numeric_limits<double>::quiet_NaN(),
    };

    for (double p : invalid_exponents) {
        SCOPED_TRACE(p);
        MetricParams params;
        params.p = p;
        EXPECT_THROW(resolve_metric("minkowski", false, params, MetricSurface::Descriptor),
                     ComparisonError);
    }
}

TEST(MetricTableTest, TverskyIsNotAMetricSpace) {
    MetricParams params;
    params.tversky_alpha = 0.3;
    params.tversky_beta = 0.7;
    const OEFP::Metric metric = resolve_metric("tversky", false, params, MetricSurface::Fingerprint);
    EXPECT_EQ(metric.Name(), OEFP::MetricName::Tversky);
    EXPECT_FALSE(metric.HasZeroSelfDistance());
    // Tversky is asymmetric. Swapping the weights changes which side of the
    // comparison is penalised while leaving every capability flag identical,
    // so the flags alone cannot pin the forwarding.
    EXPECT_DOUBLE_EQ(metric.Alpha(), 0.3);
    EXPECT_DOUBLE_EQ(metric.Beta(), 0.7);
}

TEST(MetricTableTest, TverskyRejectsWeightsOutsideTheUnitInterval) {
    const std::vector<std::pair<double, double>> rejected{
        {-1.0, 0.5},
        {0.5, -1.0},
        {1.5, 0.5},
        {0.5, 42.0},
        {std::numeric_limits<double>::quiet_NaN(), 0.5},
        {0.5, std::numeric_limits<double>::quiet_NaN()},
    };

    for (const auto& weights : rejected) {
        MetricParams params;
        params.tversky_alpha = weights.first;
        params.tversky_beta = weights.second;
        EXPECT_THROW(resolve_metric("tversky", false, params, MetricSurface::Fingerprint),
                     ComparisonError)
            << "alpha=" << weights.first << " beta=" << weights.second;
    }

    // The bounds are closed, so both endpoints must still resolve.
    MetricParams bounds;
    bounds.tversky_alpha = 0.0;
    bounds.tversky_beta = 1.0;
    EXPECT_NO_THROW(resolve_metric("tversky", false, bounds, MetricSurface::Fingerprint));
}

TEST(MetricTableTest, SeuclideanIsAnAliasOfStandardizedEuclidean) {
    MetricParams params;
    params.variances = {1.0, 2.0};
    const OEFP::Metric alias =
        resolve_metric("seuclidean", false, params, MetricSurface::Descriptor);
    const OEFP::Metric canonical =
        resolve_metric("standardized_euclidean", false, params, MetricSurface::Descriptor);
    EXPECT_EQ(alias.Name(), canonical.Name());
    // Name() is StandardizedEuclidean whatever was forwarded, so the variances
    // are the only thing that actually pins the alias to the same construction.
    EXPECT_EQ(alias.Variances(), params.variances);
    EXPECT_EQ(canonical.Variances(), params.variances);
}

TEST(MetricTableTest, MahalanobisForwardsTheInverseCovariance) {
    // variances is deliberately populated and different: the mahalanobis row
    // sits directly beneath two rows that forward params.variances, and
    // forwarding the wrong member is otherwise invisible.
    MetricParams params;
    params.variances = {9.0, 9.0};
    params.inverse_covariance = {1.0, 0.0, 0.0, 1.0};
    const OEFP::Metric metric =
        resolve_metric("mahalanobis", false, params, MetricSurface::Descriptor);
    EXPECT_EQ(metric.Name(), OEFP::MetricName::Mahalanobis);
    EXPECT_EQ(metric.InverseCovariance(), params.inverse_covariance);
    EXPECT_TRUE(metric.Variances().empty());
}

TEST(MetricTableTest, DescriptorMetricOnFingerprintSurfaceIsRejected) {
    MetricParams params;
    params.variances = {1.0};
    try {
        resolve_metric("mahalanobis", false, params, MetricSurface::Fingerprint);
        FAIL() << "expected ComparisonError";
    } catch (const ComparisonError& error) {
        EXPECT_NE(std::string(error.what()).find("comparison=\"descriptor\""), std::string::npos);
    }
}

TEST(MetricTableTest, BitSetMetricOnDescriptorSurfaceIsRejected) {
    MetricParams params;
    try {
        resolve_metric("jaccard", false, params, MetricSurface::Descriptor);
        FAIL() << "expected ComparisonError";
    } catch (const ComparisonError& error) {
        EXPECT_NE(std::string(error.what()).find("comparison=\"fingerprint\""), std::string::npos);
    }
}

TEST(MetricTableTest, SimilarityOnAMetricWithoutASimilarityFormIsRejected) {
    MetricParams params;
    try {
        resolve_metric("dice", true, params, MetricSurface::Fingerprint);
        FAIL() << "expected ComparisonError";
    } catch (const ComparisonError& error) {
        const std::string message(error.what());
        EXPECT_NE(message.find("'tanimoto'"), std::string::npos);
        EXPECT_NE(message.find("'tversky'"), std::string::npos);
    }
}

TEST(MetricTableTest, HaversineIsRejectedWithItsOwnRationale) {
    MetricParams params;
    for (MetricSurface surface : {MetricSurface::Fingerprint, MetricSurface::Descriptor}) {
        try {
            resolve_metric("haversine", false, params, surface);
            FAIL() << "expected ComparisonError";
        } catch (const ComparisonError& error) {
            EXPECT_NE(std::string(error.what()).find("latitude"), std::string::npos);
        }
    }
}

TEST(MetricTableTest, UnknownMetricListsTheSupportedNames) {
    MetricParams params;
    try {
        resolve_metric("cosine", false, params, MetricSurface::Fingerprint);
        FAIL() << "expected ComparisonError";
    } catch (const ComparisonError& error) {
        const std::string message(error.what());
        EXPECT_NE(message.find("cosine"), std::string::npos);
        EXPECT_NE(message.find(supported_metric_names(MetricSurface::Fingerprint)),
                  std::string::npos);
    }
}

// Pinned in full rather than by sampled needles. A single flipped surface flag
// routes a metric to an OEFP path that refuses it with a bare
// std::invalid_argument, bypassing this layer's ComparisonError contract. A
// later task adding a metric is expected to update this expectation
// deliberately.
TEST(MetricTableTest, SupportedNamesAreTheFullPerSurfaceAcceptSets) {
    EXPECT_EQ(supported_metric_names(MetricSurface::Fingerprint),
              "jaccard, tanimoto, dice, sokal_sneath, matching, rogers_tanimoto, russell_rao, "
              "kulsinski, sokal_michener, euclidean, manhattan, chebyshev, hamming, canberra, "
              "bray_curtis, minkowski, tversky");
    EXPECT_EQ(supported_metric_names(MetricSurface::Descriptor),
              "euclidean, manhattan, chebyshev, hamming, canberra, bray_curtis, minkowski, "
              "standardized_euclidean, seuclidean, mahalanobis");
}
