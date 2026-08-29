#include <gtest/gtest.h>
#include <oefp/oefp.h>
#include <string>
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

TEST(MetricTableTest, MinkowskiUsesTheSuppliedExponent) {
    MetricParams params;
    params.p = 0.5;
    const OEFP::Metric metric = resolve_metric("minkowski", false, params, MetricSurface::Descriptor);
    EXPECT_EQ(metric.Name(), OEFP::MetricName::Minkowski);
    EXPECT_FALSE(metric.SatisfiesTriangleInequality());
}

TEST(MetricTableTest, MinkowskiRejectsNonPositiveExponent) {
    MetricParams params;
    params.p = 0.0;
    EXPECT_THROW(resolve_metric("minkowski", false, params, MetricSurface::Descriptor),
                 ComparisonError);
}

TEST(MetricTableTest, TverskyIsNotAMetricSpace) {
    MetricParams params;
    params.tversky_alpha = 0.3;
    params.tversky_beta = 0.7;
    const OEFP::Metric metric = resolve_metric("tversky", false, params, MetricSurface::Fingerprint);
    EXPECT_EQ(metric.Name(), OEFP::MetricName::Tversky);
    EXPECT_FALSE(metric.HasZeroSelfDistance());
}

TEST(MetricTableTest, TverskyRejectsNegativeWeights) {
    MetricParams params;
    params.tversky_alpha = -1.0;
    EXPECT_THROW(resolve_metric("tversky", false, params, MetricSurface::Fingerprint),
                 ComparisonError);
}

TEST(MetricTableTest, SeuclideanIsAnAliasOfStandardizedEuclidean) {
    MetricParams params;
    params.variances = {1.0, 2.0};
    EXPECT_EQ(resolve_metric("seuclidean", false, params, MetricSurface::Descriptor).Name(),
              resolve_metric("standardized_euclidean", false, params, MetricSurface::Descriptor).Name());
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
        EXPECT_NE(message.find("tanimoto"), std::string::npos);
    }
}

TEST(MetricTableTest, SupportedNamesDifferBySurface) {
    const std::string fingerprint = supported_metric_names(MetricSurface::Fingerprint);
    const std::string descriptor = supported_metric_names(MetricSurface::Descriptor);
    EXPECT_NE(fingerprint.find("jaccard"), std::string::npos);
    EXPECT_EQ(fingerprint.find("mahalanobis"), std::string::npos);
    EXPECT_NE(descriptor.find("mahalanobis"), std::string::npos);
    EXPECT_EQ(descriptor.find("jaccard"), std::string::npos);
}
