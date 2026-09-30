/**
 * @file test_sphere_exclusion.cpp
 * @brief Sphere exclusion over distance matrices and comparisons.
 */

#include <gtest/gtest.h>

#include <algorithm>
#include <cstddef>
#include <functional>
#include <limits>
#include <numeric>
#include <random>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

#include "oecluster/Error.h"
#include "oecluster/GateFacts.h"
#include "oecluster/StorageBackend.h"
#include "oecluster/clustering/Butina.h"
#include "oecluster/clustering/DiversitySelection.h"
#include "oecluster/clustering/SphereExclusion.h"

#include "diversity_test_support.h"

using namespace OECluster;
using namespace diversity_test;

namespace {

const double NaN = std::numeric_limits<double>::quiet_NaN();
const double INF = std::numeric_limits<double>::infinity();

SphereExclusionOptions Options(double threshold,
                               SphereOrder order = SphereOrder::Input) {
    SphereExclusionOptions options;
    options.distance_threshold = threshold;
    options.order = order;
    return options;
}

SphereExclusionOptions PermutationOptions(double threshold,
                                          std::vector<size_t> permutation) {
    SphereExclusionOptions options = Options(threshold, SphereOrder::Permutation);
    options.permutation = std::move(permutation);
    return options;
}

double At(size_t n, const std::vector<double>& condensed, size_t a, size_t b) {
    if (a == b) {
        return 0.0;
    }
    const size_t i = std::min(a, b);
    const size_t j = std::max(a, b);
    return condensed[n * i - i * (i + 1) / 2 + j - i - 1];
}

std::vector<size_t> Shuffled(size_t n, unsigned seed) {
    std::vector<size_t> order(n);
    std::iota(order.begin(), order.end(), size_t{0});
    std::mt19937 generator(seed);
    std::shuffle(order.begin(), order.end(), generator);
    return order;
}

// Leader clustering from its item-centric definition: in the given order,
// each item joins the first existing center within the threshold, or becomes
// a center itself.
Clusters OracleLeader(size_t n, const std::vector<double>& condensed,
                      double threshold, const std::vector<size_t>& order) {
    Clusters clusters;
    for (const size_t item : order) {
        bool placed = false;
        for (Cluster& cluster : clusters) {
            if (At(n, condensed, cluster.front(), item) <= threshold) {
                cluster.push_back(item);
                placed = true;
                break;
            }
        }
        if (!placed) {
            clusters.push_back(Cluster{item});
        }
    }
    for (Cluster& cluster : clusters) {
        std::sort(cluster.begin() + 1, cluster.end());
    }
    return clusters;
}

void ExpectClusters(const SphereExclusionResult& result,
                    const Clusters& expected) {
    EXPECT_EQ(result.Members(), expected);
    ASSERT_EQ(result.Centers().size(), expected.size());
    for (size_t c = 0; c < expected.size(); ++c) {
        EXPECT_EQ(result.Centers()[c], expected[c].front());
        for (const size_t member : expected[c]) {
            EXPECT_EQ(result.Labels()[member], static_cast<ClusterLabel>(c));
        }
    }
}

void ExpectInvariants(const SphereExclusionResult& result, size_t n,
                      const std::vector<double>& condensed, double threshold) {
    ASSERT_EQ(result.NumSamples(), n);
    ASSERT_EQ(result.Centers().size(), result.NumClusters());
    std::vector<int> seen(n, 0);
    for (size_t c = 0; c < result.NumClusters(); ++c) {
        const Cluster& cluster = result.Members()[c];
        ASSERT_FALSE(cluster.empty());
        EXPECT_EQ(result.Centers()[c], cluster.front());
        EXPECT_TRUE(std::is_sorted(cluster.begin() + 1, cluster.end()));
        for (const size_t member : cluster) {
            ++seen[member];
            EXPECT_EQ(result.Labels()[member], static_cast<ClusterLabel>(c));
            EXPECT_LE(At(n, condensed, cluster.front(), member), threshold);
        }
    }
    for (size_t item = 0; item < n; ++item) {
        EXPECT_EQ(seen[item], 1) << "item " << item;
    }
    for (size_t a = 0; a < result.Centers().size(); ++a) {
        for (size_t b = a + 1; b < result.Centers().size(); ++b) {
            EXPECT_GT(At(n, condensed, result.Centers()[a], result.Centers()[b]),
                      threshold);
        }
    }
}

void ExpectInvalidArgument(const std::function<void()>& call,
                           const std::string& message) {
    try {
        call();
        FAIL() << "expected std::invalid_argument: " << message;
    } catch (const std::invalid_argument& error) {
        EXPECT_EQ(std::string(error.what()), message);
    }
}

void ExpectRuntimeError(const std::function<void()>& call,
                        const std::string& message) {
    try {
        call();
        FAIL() << "expected std::runtime_error: " << message;
    } catch (const std::runtime_error& error) {
        EXPECT_EQ(std::string(error.what()), message);
    }
}

}  // namespace

TEST(SphereExclusionTest, InputOrderCentersAreTheSequentialPacking) {
    for (const size_t n : {size_t{2}, size_t{7}, size_t{12}, size_t{25}}) {
        const std::vector<double> condensed = Scrambled(n);
        const DenseStorage storage = MakeStorage(n, condensed);
        std::vector<size_t> identity(n);
        std::iota(identity.begin(), identity.end(), size_t{0});
        for (const double threshold : {0.0, 2.5, 4.0}) {
            SCOPED_TRACE(::testing::Message() << "n " << n << ", threshold "
                                              << threshold);
            CirclesOptions sequential;
            sequential.method = CirclesMethod::Sequential;
            const SphereExclusionResult result =
                sphere_exclusion(storage, Options(threshold));
            EXPECT_EQ(result.Centers(),
                      circles(storage, threshold, sequential).members);
            ExpectClusters(result,
                           OracleLeader(n, condensed, threshold, identity));
            EXPECT_EQ(result.Method(), "sphere_exclusion");
        }
    }
}

TEST(SphereExclusionTest, PermutationOrderTakesCentersInTheGivenOrder) {
    const DenseStorage line = MakeStorage(5, Line(5));
    const SphereExclusionResult input = sphere_exclusion(line, Options(1.0));
    ExpectClusters(input, {{0, 1}, {2, 3}, {4}});
    EXPECT_EQ(input.Labels(), std::vector<ClusterLabel>({0, 0, 1, 1, 2}));

    const SphereExclusionResult identity =
        sphere_exclusion(line, PermutationOptions(1.0, {0, 1, 2, 3, 4}));
    EXPECT_EQ(identity.Members(), input.Members());
    EXPECT_EQ(identity.Labels(), input.Labels());

    const SphereExclusionResult reversed =
        sphere_exclusion(line, PermutationOptions(1.0, {4, 3, 2, 1, 0}));
    ExpectClusters(reversed, {{4, 3}, {2, 1}, {0}});
    EXPECT_EQ(reversed.Labels(), std::vector<ClusterLabel>({2, 1, 1, 0, 0}));
    EXPECT_EQ(reversed.Centers(), std::vector<size_t>({4, 2, 0}));

    const SphereExclusionResult mixed =
        sphere_exclusion(line, PermutationOptions(1.0, {2, 0, 4, 1, 3}));
    ExpectClusters(mixed, {{2, 1, 3}, {0}, {4}});
    EXPECT_EQ(mixed.Labels(), std::vector<ClusterLabel>({1, 0, 0, 0, 2}));

    const size_t n = 17;
    const std::vector<double> condensed = Scrambled(n);
    const DenseStorage storage = MakeStorage(n, condensed);
    for (unsigned seed = 1; seed <= 5; ++seed) {
        const std::vector<size_t> order = Shuffled(n, seed);
        ExpectClusters(sphere_exclusion(storage, PermutationOptions(2.5, order)),
                       OracleLeader(n, condensed, 2.5, order));
    }
}

TEST(SphereExclusionTest, NeighborsOrderEqualsButina) {
    struct Case {
        size_t n;
        std::vector<double> condensed;
        double threshold;
    };
    std::vector<Case> cases{{12, ScrambledSixths(12), 2.0 / 6.0},
                            {12, ScrambledSixths(12), 0.5},
                            {20, Hashed(20), 0.2},
                            {20, Hashed(20), 0.35}};
    for (unsigned seed = 1; seed <= 4; ++seed) {
        cases.push_back({40, Quantized(40, seed, 5), 0.2});
        cases.push_back({17, Quantized(17, seed, 5), 0.5});
    }
    for (const Case& fixture : cases) {
        const DenseStorage storage = MakeStorage(fixture.n, fixture.condensed);
        for (const bool reordering : {false, true}) {
            for (const size_t threads : {size_t{1}, size_t{4}}) {
                SCOPED_TRACE(::testing::Message()
                             << "n " << fixture.n << ", threshold "
                             << fixture.threshold << ", reordering "
                             << reordering << ", threads " << threads);
                ButinaOptions butina_options;
                butina_options.distance_threshold = fixture.threshold;
                butina_options.reordering = reordering;
                butina_options.num_threads = threads;
                butina_options.chunk_size = 7;
                const ButinaResult butina =
                    butina_cluster(storage, butina_options);

                SphereExclusionOptions options =
                    Options(fixture.threshold, SphereOrder::Neighbors);
                options.reordering = reordering;
                options.num_threads = threads;
                options.chunk_size = 7;
                const SphereExclusionResult result =
                    sphere_exclusion(storage, options);
                EXPECT_EQ(result.Members(), butina.Members());
                EXPECT_EQ(result.Labels(), butina.Labels());
                ExpectClusters(result, butina.Members());
            }
        }
    }
}

// Item 1 is 1.5 from both centers (a tie, kept by the earlier center), and
// item 3 is nearer the later center than the one that claimed it first.
TEST(SphereExclusionTest, NearestAssignmentMovesItemsToTheirClosestCenter) {
    const DenseStorage storage = MakeStorage(4, Positions({0.0, 1.5, 3.0, 1.6}));

    const SphereExclusionResult first = sphere_exclusion(storage, Options(2.0));
    ExpectClusters(first, {{0, 1, 3}, {2}});
    EXPECT_EQ(first.Labels(), std::vector<ClusterLabel>({0, 0, 1, 0}));

    SphereExclusionOptions options = Options(2.0);
    options.assignment = SphereAssignment::Nearest;
    const SphereExclusionResult nearest = sphere_exclusion(storage, options);
    ExpectClusters(nearest, {{0, 1}, {2, 3}});
    EXPECT_EQ(nearest.Labels(), std::vector<ClusterLabel>({0, 0, 1, 1}));
    EXPECT_EQ(nearest.Centers(), first.Centers());
}

// Butina's first-claim clusters are [[3, 1, 2], [0]]; item 1 is nearer the
// later center 0, so nearest assignment keeps the centers and moves it.
TEST(SphereExclusionTest, NeighborsWithNearestKeepsButinasCentersOnly) {
    const DenseStorage storage =
        MakeStorage(4, {0.3, 2.1, 2.0, 2.8, 0.8, 0.8});
    ButinaOptions butina_options;
    butina_options.distance_threshold = 1.0;
    EXPECT_EQ(butina_cluster(storage, butina_options).Members(),
              Clusters({{3, 1, 2}, {0}}));

    SphereExclusionOptions options = Options(1.0, SphereOrder::Neighbors);
    ExpectClusters(sphere_exclusion(storage, options), {{3, 1, 2}, {0}});
    options.assignment = SphereAssignment::Nearest;
    const SphereExclusionResult nearest = sphere_exclusion(storage, options);
    ExpectClusters(nearest, {{3, 2}, {0, 1}});
    EXPECT_EQ(nearest.Centers(), std::vector<size_t>({3, 0}));
}

// ThreadPool's ceiling division wraps to zero chunks for a chunk size near
// SIZE_MAX; the graph builder clamps to its pair count so no pair is skipped.
TEST(SphereExclusionTest, NeighborsIgnoresAHugeChunkSize) {
    const DenseStorage storage = MakeStorage(5, Positions({0, 1, 3, 7, 8}));
    const Clusters expected{{4, 3}, {1, 0}, {2}};
    SphereExclusionOptions options = Options(1.5, SphereOrder::Neighbors);
    ButinaOptions butina_options;
    butina_options.distance_threshold = 1.5;
    for (const size_t chunk : {size_t{1}, size_t{4096},
                               std::numeric_limits<size_t>::max()}) {
        SCOPED_TRACE(::testing::Message() << "chunk " << chunk);
        options.chunk_size = chunk;
        butina_options.chunk_size = chunk;
        for (const size_t threads : {size_t{1}, size_t{4}}) {
            options.num_threads = threads;
            butina_options.num_threads = threads;
            ExpectClusters(sphere_exclusion(storage, options), expected);
            EXPECT_EQ(butina_cluster(storage, butina_options).Members(),
                      expected);
        }
    }
}

// The neighbor order hands num_threads to the threshold-graph builder, which
// parallelizes over pairs; the request is capped at the item count.
TEST(SphereExclusionTest, NeighborsCapsAnAbsurdThreadCount) {
    const size_t n = 9;
    const DenseStorage storage =
        MakeStorage(n, Positions({0, 1, 3, 7, 8, 12, 13, 20, 21}));
    SphereExclusionOptions options = Options(1.5, SphereOrder::Neighbors);
    const SphereExclusionResult expected = sphere_exclusion(storage, options);
    options.num_threads = std::size_t{1} << 61;
    options.chunk_size = 1;
    EXPECT_EQ(sphere_exclusion(storage, options).Members(), expected.Members());
}

TEST(SphereExclusionTest, EveryOrderAndAssignmentKeepsTheInvariants) {
    for (const size_t n : {size_t{2}, size_t{5}, size_t{17}, size_t{40}}) {
        for (unsigned seed = 1; seed <= 4; ++seed) {
            const std::vector<double> condensed = Quantized(n, seed, 5);
            const DenseStorage storage = MakeStorage(n, condensed);
            std::vector<SphereExclusionOptions> configurations;
            for (const double threshold : {0.0, 0.2, 0.5, 1.0}) {
                configurations.push_back(Options(threshold));
                configurations.push_back(
                    PermutationOptions(threshold, Shuffled(n, seed)));
                SphereExclusionOptions neighbors =
                    Options(threshold, SphereOrder::Neighbors);
                configurations.push_back(neighbors);
                neighbors.reordering = true;
                configurations.push_back(neighbors);
            }
            for (const SphereExclusionOptions& first_options : configurations) {
                SCOPED_TRACE(::testing::Message()
                             << "n " << n << ", seed " << seed
                             << ", threshold " << first_options.distance_threshold
                             << ", order " << static_cast<int>(first_options.order)
                             << ", reordering " << first_options.reordering);
                const SphereExclusionResult first =
                    sphere_exclusion(storage, first_options);
                ExpectInvariants(first, n, condensed,
                                 first_options.distance_threshold);

                SphereExclusionOptions nearest_options = first_options;
                nearest_options.assignment = SphereAssignment::Nearest;
                const SphereExclusionResult nearest =
                    sphere_exclusion(storage, nearest_options);
                ExpectInvariants(nearest, n, condensed,
                                 first_options.distance_threshold);
                EXPECT_EQ(nearest.Centers(), first.Centers());
                const std::vector<size_t>& centers = nearest.Centers();
                for (size_t item = 0; item < n; ++item) {
                    const auto label =
                        static_cast<size_t>(nearest.Labels()[item]);
                    const double own = At(n, condensed, centers[label], item);
                    for (size_t k = 0; k < centers.size(); ++k) {
                        const double other = At(n, condensed, centers[k], item);
                        if (k < label) {
                            EXPECT_GT(other, own) << "item " << item;
                        } else {
                            EXPECT_GE(other, own) << "item " << item;
                        }
                    }
                }
            }
        }
    }
}

TEST(SphereExclusionTest, RefusesANonFiniteDistance) {
    std::vector<double> condensed = Positions({0.0, 5.0, 10.0});
    condensed[2] = NaN;  // (1, 2)
    const DenseStorage storage = MakeStorage(3, condensed);
    const std::string message =
        "sphere_exclusion read a non-finite distance between items 1 and 2";
    ExpectRuntimeError([&] { sphere_exclusion(storage, Options(1.0)); }, message);
    ExpectRuntimeError(
        [&] { sphere_exclusion(storage, PermutationOptions(1.0, {0, 1, 2})); },
        message);
    ExpectRuntimeError(
        [&] {
            sphere_exclusion(storage, Options(1.0, SphereOrder::Neighbors));
        },
        message);

    // The neighbor order scans before building the graph, so an entry the
    // graph would read as "not a neighbor" still raises.
    std::vector<double> infinite = Positions({0.0, 5.0, 10.0});
    infinite[0] = INF;  // (0, 1)
    ExpectRuntimeError(
        [&] {
            sphere_exclusion(MakeStorage(3, infinite),
                             Options(1.0, SphereOrder::Neighbors));
        },
        "sphere_exclusion read a non-finite distance between items 0 and 1");
}

// The first-claim pass never reads (1, 2): item 1 is claimed by item 0 before
// item 2 becomes a center. Only the nearest pass reads it.
TEST(SphereExclusionTest, NearestAssignmentRefusesANonFiniteDistance) {
    const DenseStorage storage = MakeStorage(3, {1.0, 5.0, NaN});
    ExpectClusters(sphere_exclusion(storage, Options(1.0)), {{0, 1}, {2}});

    SphereExclusionOptions options = Options(1.0);
    options.assignment = SphereAssignment::Nearest;
    ExpectRuntimeError(
        [&] { sphere_exclusion(storage, options); },
        "sphere_exclusion read a non-finite distance between items 2 and 1");
}

TEST(SphereExclusionValidationTest, RefusesEachInvalidRequest) {
    const DenseStorage storage = MakeStorage(4, Line(4));

    for (const double threshold : {NaN, INF}) {
        ExpectInvalidArgument(
            [&] { sphere_exclusion(storage, Options(threshold)); },
            "sphere_exclusion distance_threshold must be finite");
    }
    ExpectInvalidArgument([&] { sphere_exclusion(storage, Options(-1.0)); },
                          "sphere_exclusion distance_threshold must be "
                          "non-negative");

    ExpectInvalidArgument(
        [&] {
            sphere_exclusion(storage, Options(1.0, static_cast<SphereOrder>(7)));
        },
        "Unknown sphere_exclusion order");
    SphereExclusionOptions bad_assignment = Options(1.0);
    bad_assignment.assignment = static_cast<SphereAssignment>(7);
    ExpectInvalidArgument([&] { sphere_exclusion(storage, bad_assignment); },
                          "Unknown sphere_exclusion assignment");

    SphereExclusionOptions input_reordering = Options(1.0);
    input_reordering.reordering = true;
    ExpectInvalidArgument([&] { sphere_exclusion(storage, input_reordering); },
                          "sphere_exclusion reordering requires the Neighbors "
                          "order");
    SphereExclusionOptions permutation_reordering =
        PermutationOptions(1.0, {0, 1, 2, 3});
    permutation_reordering.reordering = true;
    ExpectInvalidArgument(
        [&] { sphere_exclusion(storage, permutation_reordering); },
        "sphere_exclusion reordering requires the Neighbors order");

    for (const SphereOrder order : {SphereOrder::Input, SphereOrder::Neighbors}) {
        SphereExclusionOptions stray = Options(1.0, order);
        stray.permutation = {0, 1, 2, 3};
        ExpectInvalidArgument([&] { sphere_exclusion(storage, stray); },
                              "sphere_exclusion permutation requires the "
                              "Permutation order");
    }

    SphereExclusionOptions zero_chunk = Options(1.0);
    zero_chunk.chunk_size = 0;
    ExpectInvalidArgument([&] { sphere_exclusion(storage, zero_chunk); },
                          "sphere_exclusion chunk_size must be at least one");

    ExpectInvalidArgument(
        [&] { sphere_exclusion(storage, PermutationOptions(1.0, {})); },
        "sphere_exclusion permutation has 0 entries for 4 items");
    ExpectInvalidArgument(
        [&] { sphere_exclusion(storage, PermutationOptions(1.0, {0, 1, 2})); },
        "sphere_exclusion permutation has 3 entries for 4 items");
    ExpectInvalidArgument(
        [&] {
            sphere_exclusion(storage, PermutationOptions(1.0, {0, 1, 2, 4}));
        },
        "sphere_exclusion permutation entry 4 is outside the item range");
    ExpectInvalidArgument(
        [&] {
            sphere_exclusion(storage, PermutationOptions(1.0, {0, 1, 1, 3}));
        },
        "sphere_exclusion permutation repeats item 1");
}

TEST(SphereExclusionValidationTest, RefusesStorageItCannotRead) {
    ExpectInvalidArgument(
        [] { sphere_exclusion(SparseStorage(4, 0.5), Options(1.0)); },
        "sphere_exclusion requires complete pairwise distances; SparseStorage "
        "is not supported");
    ExpectInvalidArgument(
        [] { sphere_exclusion(NullDataStorage(4), Options(1.0)); },
        "sphere_exclusion requires contiguous dense or memory-mapped storage");

    // Options are checked before the storage, and the storage before the
    // permutation's length.
    ExpectInvalidArgument(
        [] { sphere_exclusion(SparseStorage(4, 0.5), Options(-1.0)); },
        "sphere_exclusion distance_threshold must be non-negative");
    ExpectInvalidArgument(
        [] { sphere_exclusion(NullDataStorage(4), PermutationOptions(1.0, {0})); },
        "sphere_exclusion requires contiguous dense or memory-mapped storage");
}

TEST(SphereExclusionTest, ZeroAndOneItems) {
    const DenseStorage empty(0);
    const DenseStorage single(1);
    std::vector<SphereExclusionOptions> configurations{
        Options(1.0), Options(1.0, SphereOrder::Neighbors)};
    SphereExclusionOptions nearest = Options(1.0);
    nearest.assignment = SphereAssignment::Nearest;
    configurations.push_back(nearest);
    for (const SphereExclusionOptions& options : configurations) {
        const SphereExclusionResult none = sphere_exclusion(empty, options);
        EXPECT_EQ(none.NumSamples(), 0u);
        EXPECT_EQ(none.NumClusters(), 0u);
        EXPECT_TRUE(none.Centers().empty());

        ExpectClusters(sphere_exclusion(single, options), {{0}});
    }
    EXPECT_EQ(sphere_exclusion(empty, PermutationOptions(1.0, {})).NumClusters(),
              0u);
    ExpectClusters(sphere_exclusion(single, PermutationOptions(1.0, {0})), {{0}});
    ExpectInvalidArgument([&] { sphere_exclusion(empty, Options(-1.0)); },
                          "sphere_exclusion distance_threshold must be "
                          "non-negative");
}
