/**
 * @file test_set_diversity.cpp
 * @brief The Vendi score and log-determinant diversity over distance matrices
 *        and comparisons.
 */

#include <gtest/gtest.h>

#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "oecluster/Error.h"
#include "oecluster/GateFacts.h"
#include "oecluster/StorageBackend.h"
#include "oecluster/clustering/SetDiversity.h"

#include "diversity_test_support.h"

using namespace OECluster;
using namespace diversity_test;

namespace {

const double NaN = std::numeric_limits<double>::quiet_NaN();
const double INF = std::numeric_limits<double>::infinity();

VendiOptions Order(unsigned order) {
    VendiOptions options;
    options.order = order;
    return options;
}

LogDetOptions Ridge(double ridge) {
    LogDetOptions options;
    options.ridge = ridge;
    return options;
}

// Distances in [0, 1) from a fixed formula, varied enough that the kernel has
// no special structure.
std::vector<double> UnitDistances(size_t n) {
    return Condensed(n, [](size_t i, size_t j) {
        return static_cast<double>((i * 7919 + j * 104729) % 1000) / 1000.0;
    });
}

std::vector<double> Constant(size_t n, double distance) {
    return Condensed(n, [=](size_t, size_t) { return distance; });
}

// d(0, 1) = d(0, 2) = 0 and d(1, 2) = 1: the complement kernel is
// [[1, 1, 1], [1, 1, 0], [1, 0, 1]], with eigenvalues 1 - sqrt 2, 1, 1 + sqrt 2.
std::vector<double> IndefiniteDistances() { return {0.0, 0.0, 1.0}; }

void ExpectInvalidArgument(const std::function<void()>& call,
                           const std::string& message) {
    try {
        call();
        FAIL() << "expected std::invalid_argument: " << message;
    } catch (const std::invalid_argument& error) {
        EXPECT_EQ(std::string(error.what()), message);
    }
}

// -inf == -inf, but EXPECT_NEAR on two infinities computes inf - inf = NaN.
void ExpectSameScore(double actual, double expected, double tolerance) {
    if (std::isinf(expected)) {
        EXPECT_EQ(actual, expected);
    } else {
        EXPECT_NEAR(actual, expected, tolerance);
    }
}

}  // namespace

TEST(SetDiversityTest, AnIdentityKernelScoresTheItemCount) {
    const DenseStorage storage = MakeStorage(5, Constant(5, 1.0));

    const VendiResult one = vendi_score(storage, Order(1));
    EXPECT_NEAR(one.score, 5.0, 1e-12);
    EXPECT_EQ(one.order, 1u);
    EXPECT_EQ(one.size, 5u);
    EXPECT_EQ(one.kernel, DiversityKernel::Complement);
    EXPECT_NEAR(one.min_eigenvalue, 1.0, 1e-12);
    EXPECT_EQ(one.negative_mass, 0.0);
    EXPECT_DOUBLE_EQ(vendi_score(storage, Order(2)).score, 5.0);

    const LogDetResult plain = logdet_diversity(storage, Ridge(0.0));
    EXPECT_NEAR(plain.score, 0.0, 1e-12);
    EXPECT_EQ(plain.nonpositive_count, 0u);
    EXPECT_EQ(plain.size, 5u);
    EXPECT_NEAR(logdet_diversity(storage, Ridge(0.5)).score,
                5.0 * std::log(1.5), 1e-12);
}

TEST(SetDiversityTest, IdenticalItemsScoreOne) {
    const DenseStorage storage = MakeStorage(4, Constant(4, 0.0));

    EXPECT_NEAR(vendi_score(storage, Order(1)).score, 1.0, 1e-12);
    EXPECT_DOUBLE_EQ(vendi_score(storage, Order(2)).score, 1.0);

    const LogDetResult singular = logdet_diversity(storage, Ridge(0.0));
    EXPECT_EQ(singular.score, -INF);
    EXPECT_EQ(singular.nonpositive_count, 3u);
    EXPECT_NEAR(singular.min_eigenvalue, 0.0, 1e-12);

    const LogDetResult ridged = logdet_diversity(storage, Ridge(0.5));
    EXPECT_EQ(ridged.nonpositive_count, 0u);
    EXPECT_EQ(ridged.ridge, 0.5);
    EXPECT_NEAR(ridged.score, 3.0 * std::log(0.5) + std::log(4.5), 1e-12);
}

TEST(SetDiversityTest, EqualBlocksScoreTheBlockCount) {
    // Three blocks of four identical items, the blocks at distance 1.
    const DenseStorage storage = MakeStorage(12, Condensed(12, [](size_t i, size_t j) {
        return i / 4 == j / 4 ? 0.0 : 1.0;
    }));
    EXPECT_NEAR(vendi_score(storage, Order(1)).score, 3.0, 1e-12);
    EXPECT_DOUBLE_EQ(vendi_score(storage, Order(2)).score, 3.0);
}

TEST(SetDiversityTest, AnIndefiniteKernelDropsWithoutRenormalizing) {
    const VendiResult result =
        vendi_score(MakeStorage(3, IndefiniteDistances()), Order(1));
    const double root2 = std::sqrt(2.0);
    EXPECT_NEAR(result.min_eigenvalue, 1.0 - root2, 1e-12);
    EXPECT_NEAR(result.negative_mass, (root2 - 1.0) / 3.0, 1e-12);

    const double p1 = 1.0 / 3.0;
    const double p2 = (1.0 + root2) / 3.0;
    EXPECT_NEAR(result.score, std::exp(-(p1 * std::log(p1) + p2 * std::log(p2))),
                1e-12);
    EXPECT_NEAR(result.score, 1.7177654733746037, 1e-12);

    // Renormalizing the kept mass would give about 1.83 instead.
    const double kept = p1 + p2;
    const double q1 = p1 / kept;
    const double q2 = p2 / kept;
    const double renormalized = std::exp(-(q1 * std::log(q1) + q2 * std::log(q2)));
    EXPECT_GT(std::fabs(result.score - renormalized), 1e-3);
}

TEST(SetDiversityTest, OrderTwoIsTheFrobeniusForm) {
    const VendiResult result =
        vendi_score(MakeStorage(3, IndefiniteDistances()), Order(2));
    // n^2 / (n + 2 * sum_{i<j} K_ij^2) = 9 / (3 + 2 * 2).
    EXPECT_DOUBLE_EQ(result.score, 9.0 / 7.0);
    EXPECT_EQ(result.order, 2u);
    EXPECT_TRUE(std::isnan(result.min_eigenvalue));
    EXPECT_TRUE(std::isnan(result.negative_mass));
}

TEST(SetDiversityTest, TheLaplacianKernelMatchesTheHandComputedCase) {
    // Points at 0, 1 and 2 with bandwidth 2: K = [[1, a, b], [a, 1, a], [b, a, 1]]
    // with a = exp(-1/2) and b = exp(-1). (1, 0, -1) gives 1 - b; the
    // symmetric vectors (x, y, x) give (2 + b +- sqrt(b^2 + 8a^2)) / 2.
    const DenseStorage storage = MakeStorage(3, Positions({0.0, 1.0, 2.0}));
    const double a = std::exp(-0.5);
    const double b = std::exp(-1.0);
    const double root = std::sqrt(b * b + 8.0 * a * a);
    const std::vector<double> lambdas = {(2.0 + b - root) / 2.0, 1.0 - b,
                                         (2.0 + b + root) / 2.0};

    VendiOptions vendi = Order(1);
    vendi.kernel = DiversityKernel::Laplacian;
    vendi.bandwidth = 2.0;
    double entropy = 0.0;
    double log_det = 0.0;
    for (const double lambda : lambdas) {
        entropy -= lambda / 3.0 * std::log(lambda / 3.0);
        log_det += std::log(lambda);
    }
    const VendiResult one = vendi_score(storage, vendi);
    EXPECT_NEAR(one.score, std::exp(entropy), 1e-12);
    EXPECT_NEAR(one.min_eigenvalue, lambdas[0], 1e-12);
    EXPECT_EQ(one.kernel, DiversityKernel::Laplacian);

    vendi.order = 2;
    EXPECT_NEAR(vendi_score(storage, vendi).score,
                9.0 / (3.0 + 2.0 * (2.0 * a * a + b * b)), 1e-14);

    LogDetOptions logdet;
    logdet.kernel = DiversityKernel::Laplacian;
    logdet.bandwidth = 2.0;
    EXPECT_NEAR(logdet_diversity(storage, logdet).score, log_det, 1e-12);
}

TEST(SetDiversityTest, ANonsingularIndefiniteKernelScoresMinusInfinity) {
    // Two copies of the indefinite triple, the copies at distance 1: the
    // kernel is block-diagonal with eigenvalues (1 - sqrt 2) twice, 1 twice
    // and (1 + sqrt 2) twice, so its determinant is (-1)^2 = 1 > 0. The
    // positive-definite rule still refuses it.
    const DenseStorage storage = MakeStorage(6, Condensed(6, [](size_t i, size_t j) {
        if (i / 3 != j / 3) {
            return 1.0;
        }
        return (i % 3 == 1 && j % 3 == 2) ? 1.0 : 0.0;
    }));
    const LogDetResult result = logdet_diversity(storage, Ridge(0.0));
    EXPECT_EQ(result.score, -INF);
    EXPECT_EQ(result.nonpositive_count, 2u);
    EXPECT_NEAR(result.min_eigenvalue, 1.0 - std::sqrt(2.0), 1e-12);
}

TEST(SetDiversityTest, TheRidgeToleranceIsTheBoundary) {
    // All ones, n = 4: the zero modes survive only when ridge clears
    // tol = n * eps * (n + ridge), about 16 eps.
    const DenseStorage storage = MakeStorage(4, Constant(4, 0.0));
    const double tolerance = 16.0 * std::numeric_limits<double>::epsilon();

    const LogDetResult below = logdet_diversity(storage, Ridge(0.1 * tolerance));
    EXPECT_EQ(below.score, -INF);
    EXPECT_EQ(below.nonpositive_count, 3u);

    const LogDetResult above = logdet_diversity(storage, Ridge(10.0 * tolerance));
    EXPECT_TRUE(std::isfinite(above.score));
    EXPECT_EQ(above.nonpositive_count, 0u);
}

TEST(SetDiversityComparisonTest, MatchesTheMatrix) {
    const size_t n = 12;
    const std::vector<double> distances = UnitDistances(n);
    const DenseStorage storage = MakeStorage(n, distances);
    for (const DiversityKernel kernel :
         {DiversityKernel::Complement, DiversityKernel::Laplacian}) {
        const double bandwidth = kernel == DiversityKernel::Laplacian ? 0.5 : NaN;
        for (const unsigned order : {1u, 2u}) {
            VendiOptions options = Order(order);
            options.kernel = kernel;
            options.bandwidth = bandwidth;
            TableComparison table(n, distances);
            const double expected = vendi_score(storage, options).score;
            const double actual = vendi_score(table, options).score;
            if (order == 2) {
                EXPECT_EQ(actual, expected) << "order 2 must agree bit-for-bit";
            } else {
                EXPECT_NEAR(actual, expected, 1e-12);
            }
        }
        LogDetOptions logdet = Ridge(0.1);
        logdet.kernel = kernel;
        logdet.bandwidth = bandwidth;
        TableComparison table(n, distances);
        ExpectSameScore(logdet_diversity(table, logdet).score,
                        logdet_diversity(storage, logdet).score, 1e-12);
    }
}

TEST(SetDiversityComparisonTest, IsBitIdenticalAtEveryThreadCountAndChunkSize) {
    const size_t n = 23;
    const std::vector<double> distances = UnitDistances(n);
    const double order_two = vendi_score(MakeStorage(n, distances), Order(2)).score;
    TableComparison baseline_table(n, distances);
    VendiOptions baseline_options = Order(1);
    baseline_options.num_threads = 1;
    const double order_one = vendi_score(baseline_table, baseline_options).score;
    LogDetOptions baseline_logdet = Ridge(0.1);
    baseline_logdet.num_threads = 1;
    TableComparison baseline_logdet_table(n, distances);
    const double logdet =
        logdet_diversity(baseline_logdet_table, baseline_logdet).score;

    for (const size_t threads : {size_t{1}, size_t{2}, size_t{7}}) {
        for (const size_t chunk : {size_t{1}, size_t{3}, size_t{256}}) {
            SCOPED_TRACE("threads " + std::to_string(threads) + ", chunk " +
                         std::to_string(chunk));
            for (const unsigned order : {1u, 2u}) {
                VendiOptions options = Order(order);
                options.num_threads = threads;
                options.chunk_size = chunk;
                TableComparison table(n, distances);
                EXPECT_EQ(vendi_score(table, options).score,
                          order == 2 ? order_two : order_one);
                EXPECT_LE(table.NumClones(), threads);
            }
            LogDetOptions options = Ridge(0.1);
            options.num_threads = threads;
            options.chunk_size = chunk;
            TableComparison table(n, distances);
            EXPECT_EQ(logdet_diversity(table, options).score, logdet);
            EXPECT_LE(table.NumClones(), threads);
        }
    }
}

// chunk_size 1 keeps the comparisons on the ThreadPool path, as in
// CirclesComparisonTest.CapsAnAbsurdThreadCount.
TEST(SetDiversityComparisonTest, CapsAnAbsurdThreadCount) {
    const size_t n = 9;
    const std::vector<double> distances = UnitDistances(n);
    VendiOptions options = Order(2);
    options.num_threads = std::size_t{1} << 61;
    options.chunk_size = 1;
    TableComparison table(n, distances);
    EXPECT_EQ(vendi_score(table, options).score,
              vendi_score(MakeStorage(n, distances), Order(2)).score);
}

TEST(SetDiversityComparisonTest, ProvesCloneIsolation) {
    const size_t n = 40;
    const std::vector<double> distances = UnitDistances(n);
    const DenseStorage storage = MakeStorage(n, distances);
    for (const unsigned order : {1u, 2u}) {
        SCOPED_TRACE("order " + std::to_string(order));
        VendiOptions options = Order(order);
        options.num_threads = 4;
        options.chunk_size = 1;
        IsolationComparison comparison(n, distances);

        const VendiResult result = vendi_score(comparison, options);

        EXPECT_EQ(comparison.Violations(), 0u);
        EXPECT_TRUE(comparison.OverlapObserved())
            << "Overlap not observed; test may be flaky on this machine";
        EXPECT_EQ(result.score, vendi_score(storage, Order(order)).score);
    }
}

TEST(SetDiversityComparisonTest, AMaximalChunkSizeRunsAsOneChunk) {
    const size_t n = 10;
    const std::vector<double> distances = UnitDistances(n);
    const DenseStorage storage = MakeStorage(n, distances);
    for (const unsigned order : {1u, 2u}) {
        VendiOptions options = Order(order);
        options.chunk_size = std::numeric_limits<size_t>::max();
        TableComparison table(n, distances);
        EXPECT_EQ(vendi_score(table, options).score,
                  vendi_score(storage, Order(order)).score);
    }
    LogDetOptions logdet = Ridge(0.1);
    logdet.chunk_size = std::numeric_limits<size_t>::max();
    TableComparison table(n, distances);
    ExpectSameScore(logdet_diversity(table, logdet).score,
                    logdet_diversity(storage, Ridge(0.1)).score, 0.0);
}

TEST(SetDiversityValidationTest, RefusesEachInvalidOption) {
    const DenseStorage storage = MakeStorage(3, UnitDistances(3));

    for (const unsigned order : {0u, 3u}) {
        ExpectInvalidArgument([&] { vendi_score(storage, Order(order)); },
                              "Vendi score order must be 1 or 2");
    }

    VendiOptions bad_kernel;
    bad_kernel.kernel = static_cast<DiversityKernel>(7);
    ExpectInvalidArgument([&] { vendi_score(storage, bad_kernel); },
                          "Vendi score kernel is not a known DiversityKernel");

    for (const double bandwidth : {NaN, INF, 0.0, -1.0}) {
        VendiOptions laplacian;
        laplacian.kernel = DiversityKernel::Laplacian;
        laplacian.bandwidth = bandwidth;
        ExpectInvalidArgument(
            [&] { vendi_score(storage, laplacian); },
            "Vendi score Laplacian kernel requires a finite, positive bandwidth");
    }

    VendiOptions stray_bandwidth;
    stray_bandwidth.bandwidth = 1.0;
    ExpectInvalidArgument(
        [&] { vendi_score(storage, stray_bandwidth); },
        "Vendi score bandwidth applies only to the Laplacian kernel");

    VendiOptions zero_ceiling;
    zero_ceiling.max_exact = 0;
    ExpectInvalidArgument([&] { vendi_score(storage, zero_ceiling); },
                          "Vendi score max_exact must be at least one");

    VendiOptions zero_chunk;
    zero_chunk.chunk_size = 0;
    ExpectInvalidArgument([&] { vendi_score(storage, zero_chunk); },
                          "Vendi score chunk_size must be at least one");

    LogDetOptions logdet_kernel;
    logdet_kernel.kernel = static_cast<DiversityKernel>(7);
    ExpectInvalidArgument(
        [&] { logdet_diversity(storage, logdet_kernel); },
        "Log-determinant diversity kernel is not a known DiversityKernel");

    for (const double ridge : {-1.0, INF, NaN}) {
        ExpectInvalidArgument(
            [&] { logdet_diversity(storage, Ridge(ridge)); },
            "Log-determinant diversity ridge must be finite and non-negative");
    }

    LogDetOptions logdet_ceiling;
    logdet_ceiling.max_exact = 0;
    ExpectInvalidArgument(
        [&] { logdet_diversity(storage, logdet_ceiling); },
        "Log-determinant diversity max_exact must be at least one");

    LogDetOptions logdet_chunk;
    logdet_chunk.chunk_size = 0;
    ExpectInvalidArgument(
        [&] { logdet_diversity(storage, logdet_chunk); },
        "Log-determinant diversity chunk_size must be at least one");
}

TEST(SetDiversityValidationTest, RefusesStorageItCannotRead) {
    ExpectInvalidArgument(
        [] { vendi_score(SparseStorage(4, 0.5), VendiOptions()); },
        "Vendi score requires complete pairwise distances; SparseStorage is "
        "not supported");
    ExpectInvalidArgument(
        [] { vendi_score(NullDataStorage(4), VendiOptions()); },
        "Vendi score requires contiguous dense or memory-mapped storage");
    ExpectInvalidArgument(
        [] { vendi_score(DenseStorage(0), VendiOptions()); },
        "Vendi score requires at least one item");
    ExpectInvalidArgument(
        [] { logdet_diversity(DenseStorage(0), LogDetOptions()); },
        "Log-determinant diversity requires at least one item");
}

TEST(SetDiversityValidationTest, ValidatesOptionsBeforeTheInput) {
    // Options (step 1) before storage (step 2) and before the item count
    // (step 3): the chunk check is split from the size check for this.
    ExpectInvalidArgument([] { vendi_score(NullDataStorage(4), Order(3)); },
                          "Vendi score order must be 1 or 2");
    VendiOptions zero_chunk;
    zero_chunk.chunk_size = 0;
    ExpectInvalidArgument([&] { vendi_score(DenseStorage(0), zero_chunk); },
                          "Vendi score chunk_size must be at least one");
    TableComparison empty(0, {});
    ExpectInvalidArgument([&] { vendi_score(empty, zero_chunk); },
                          "Vendi score chunk_size must be at least one");
}

TEST(SetDiversityValidationTest, RefusesAnInputAboveTheExactCeiling) {
    CountingComparison counter(5, GateFacts());
    VendiOptions vendi;
    vendi.max_exact = 4;
    ExpectInvalidArgument(
        [&] { vendi_score(counter, vendi); },
        "Vendi score computes an exact spectrum of at most max_exact (4) items, "
        "but the input has 5; raise max_exact (memory grows as 8n^2 bytes and "
        "time as n^3) or use the Vendi score with order=2, which needs no "
        "spectrum");
    LogDetOptions logdet;
    logdet.max_exact = 4;
    ExpectInvalidArgument(
        [&] { logdet_diversity(counter, logdet); },
        "Log-determinant diversity computes an exact spectrum of at most "
        "max_exact (4) items, but the input has 5; raise max_exact (memory "
        "grows as 8n^2 bytes and time as n^3) or use the Vendi score with "
        "order=2, which needs no spectrum");
    EXPECT_EQ(counter.Count(), 0u) << "Compare called before the ceiling";

    // Order 2 needs no spectrum, so the ceiling does not apply to it.
    vendi.order = 2;
    EXPECT_DOUBLE_EQ(vendi_score(counter, vendi).score, 1.0);
    EXPECT_EQ(counter.Count(), 10u);
}

TEST(SetDiversityValidationTest, RefusesADenseKernelTooLargeToAllocate) {
    CountingComparison huge(size_t{1} << 33, GateFacts());
    VendiOptions vendi;
    vendi.max_exact = std::numeric_limits<size_t>::max();
    ExpectInvalidArgument(
        [&] { vendi_score(huge, vendi); },
        "A dense 8589934592 x 8589934592 kernel is too large to allocate");
    LogDetOptions logdet;
    logdet.max_exact = std::numeric_limits<size_t>::max();
    ExpectInvalidArgument(
        [&] { logdet_diversity(huge, logdet); },
        "A dense 8589934592 x 8589934592 kernel is too large to allocate");
    EXPECT_EQ(huge.Count(), 0u);
}

TEST(SetDiversityValidationTest, RefusesDistancesTheKernelCannotUse) {
    const std::string above =
        "complement kernel requires distances in [0, 1], but d(0, 2) = 1.5; use "
        "the Laplacian kernel (kernel=\"laplacian\") for other distances";
    const std::string below =
        "complement kernel requires distances in [0, 1], but d(0, 2) = -0.25; "
        "use the Laplacian kernel (kernel=\"laplacian\") for other distances";
    const std::vector<double> too_far = {0.2, 1.5, 0.3};
    const std::vector<double> negative = {0.2, -0.25, 0.3};
    const std::vector<double> missing = {0.2, NaN, 0.3};

    for (const unsigned order : {1u, 2u}) {
        SCOPED_TRACE("order " + std::to_string(order));
        ExpectInvalidArgument(
            [&] { vendi_score(MakeStorage(3, too_far), Order(order)); },
            "Vendi score " + above);
        TableComparison table(3, too_far);
        ExpectInvalidArgument([&] { vendi_score(table, Order(order)); },
                              "Vendi score " + above);
    }
    ExpectInvalidArgument(
        [&] { vendi_score(MakeStorage(3, negative), Order(1)); },
        "Vendi score " + below);
    ExpectInvalidArgument(
        [&] { logdet_diversity(MakeStorage(3, too_far), LogDetOptions()); },
        "Log-determinant diversity " + above);

    VendiOptions laplacian;
    laplacian.kernel = DiversityKernel::Laplacian;
    laplacian.bandwidth = 1.0;
    ExpectInvalidArgument(
        [&] { vendi_score(MakeStorage(3, negative), laplacian); },
        "Vendi score Laplacian kernel requires non-negative distances, but "
        "d(0, 2) = -0.25");
    TableComparison negative_table(3, negative);
    ExpectInvalidArgument(
        [&] { vendi_score(negative_table, laplacian); },
        "Vendi score Laplacian kernel requires non-negative distances, but "
        "d(0, 2) = -0.25");

    ExpectInvalidArgument(
        [&] { vendi_score(MakeStorage(3, missing), Order(1)); },
        "Vendi score read a non-finite distance between items 0 and 2");
    TableComparison missing_table(3, missing);
    ExpectInvalidArgument(
        [&] { vendi_score(missing_table, Order(2)); },
        "Vendi score read a non-finite distance between items 0 and 2");
}

TEST(SetDiversityComparisonTest, RefusesComparisonsItsFactsRuleOut) {
    const auto expect = [](GateFacts facts, const std::string& message) {
        SCOPED_TRACE(message);
        CountingComparison counter(6, facts);
        try {
            vendi_score(counter, VendiOptions());
            FAIL() << "expected ComparisonError: " << message;
        } catch (const ComparisonError& error) {
            EXPECT_EQ(std::string(error.what()), "Vendi score " + message);
        }
        try {
            logdet_diversity(counter, LogDetOptions());
            FAIL() << "expected ComparisonError: " << message;
        } catch (const ComparisonError& error) {
            EXPECT_EQ(std::string(error.what()),
                      "Log-determinant diversity " + message);
        }
        EXPECT_EQ(counter.Count(), 0u) << "Compare called before facts refusal";
    };

    GateFacts similarity;
    similarity.is_distance = Capability::No;
    expect(similarity,
           "requires distances, but the comparison reports similarities");

    GateFacts nonzero_self;
    nonzero_self.is_distance = Capability::Yes;
    nonzero_self.zero_self = Capability::No;
    expect(nonzero_self,
           "requires a zero self-distance, but the comparison reports that "
           "d(x, x) is not zero");

    GateFacts nan_present;
    nan_present.is_distance = Capability::Yes;
    nan_present.zero_self = Capability::Yes;
    nan_present.data_integrity = DataIntegrity::NaNPresent;
    expect(nan_present,
           "cannot rank distances the comparison declares may be non-finite "
           "(missing='propagate')");

    GateFacts subset_scored;
    subset_scored.is_distance = Capability::Yes;
    subset_scored.zero_self = Capability::Yes;
    subset_scored.data_integrity = DataIntegrity::SubsetScored;
    expect(subset_scored,
           "cannot rank distances scored on per-pair feature subsets "
           "(missing='ignore'); they are not mutually comparable");

    TableComparison empty(0, {});
    ExpectInvalidArgument([&] { vendi_score(empty, VendiOptions()); },
                          "Vendi score requires at least one item");
}

TEST(SetDiversityTest, OneItemScoresOne) {
    const DenseStorage storage = MakeStorage(1, {});
    EXPECT_EQ(vendi_score(storage, Order(1)).score, 1.0);
    EXPECT_EQ(vendi_score(storage, Order(2)).score, 1.0);
    EXPECT_NEAR(logdet_diversity(storage, Ridge(0.5)).score, std::log(1.5),
                1e-15);
}
