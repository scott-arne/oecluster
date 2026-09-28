/**
 * @file test_circles.cpp
 * @brief The #Circles packing over distance matrices and comparisons.
 */

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "oecluster/StorageBackend.h"
#include "oecluster/clustering/DiversitySelection.h"

#include "diversity_test_support.h"

using namespace OECluster;
using namespace diversity_test;

namespace {

const double NaN = std::numeric_limits<double>::quiet_NaN();

CirclesOptions MethodOptions(CirclesMethod method) {
    CirclesOptions options;
    options.method = method;
    return options;
}

double CondensedAt(size_t n, const std::vector<double>& condensed, size_t a,
                   size_t b) {
    const size_t i = std::min(a, b);
    const size_t j = std::max(a, b);
    return condensed[n * i - i * (i + 1) / 2 + j - i - 1];
}

// The reference greedy pass, written from its definition.
std::vector<size_t> OracleSequential(size_t n,
                                     const std::vector<double>& condensed,
                                     double threshold) {
    std::vector<size_t> members;
    for (size_t candidate = 0; candidate < n; ++candidate) {
        bool accept = true;
        for (const size_t member : members) {
            accept = accept && CondensedAt(n, condensed, member, candidate) > threshold;
        }
        if (accept) {
            members.push_back(candidate);
        }
    }
    return members;
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

}  // namespace

TEST(CirclesTest, BothMethodsFindTheKnownPacking) {
    const DenseStorage storage =
        MakeStorage(6, Positions({0, 0.5, 1, 3, 3.2, 10}));

    const CirclesResult maxmin =
        circles(storage, 1.0, MethodOptions(CirclesMethod::MaxMin));
    const CirclesResult sequential =
        circles(storage, 1.0, MethodOptions(CirclesMethod::Sequential));

    EXPECT_EQ(maxmin.count, 3u);
    EXPECT_EQ(maxmin.members, std::vector<size_t>({0, 5, 4}));
    EXPECT_EQ(maxmin.threshold, 1.0);
    EXPECT_EQ(maxmin.method, CirclesMethod::MaxMin);
    EXPECT_EQ(sequential.count, 3u);
    EXPECT_EQ(sequential.members, std::vector<size_t>({0, 3, 5}));
    EXPECT_EQ(sequential.method, CirclesMethod::Sequential);
}

// Either method yields a valid packing, hence a lower bound; neither dominates.
TEST(CirclesTest, TheMethodsCanDisagree) {
    const DenseStorage storage = MakeStorage(5, Positions({0, 2, 3, 4, 6}));

    EXPECT_EQ(circles(storage, 1.5, MethodOptions(CirclesMethod::MaxMin)).members,
              std::vector<size_t>({0, 4, 2}));
    EXPECT_EQ(
        circles(storage, 1.5, MethodOptions(CirclesMethod::Sequential)).members,
        std::vector<size_t>({0, 1, 3, 4}));
}

TEST(CirclesTest, MaxMinIsTheThresholdSelectionFromItemZero) {
    for (const size_t n : {size_t{1}, size_t{2}, size_t{7}, size_t{12}}) {
        const DenseStorage storage = MakeStorage(n, Scrambled(n));
        for (const double threshold : {0.0, 1.0, 2.5, 4.0, 6.0}) {
            MaxMinOptions options;
            options.threshold = threshold;
            EXPECT_EQ(
                circles(storage, threshold, MethodOptions(CirclesMethod::MaxMin))
                    .members,
                maxmin_select(storage, options).indices)
                << "n " << n << ", threshold " << threshold;
        }
    }
}

TEST(CirclesTest, SequentialMatchesTheReferencePass) {
    for (const size_t n : {size_t{1}, size_t{2}, size_t{7}, size_t{12}}) {
        const std::vector<double> condensed = Scrambled(n);
        const DenseStorage storage = MakeStorage(n, condensed);
        for (const double threshold : {0.0, 1.0, 2.5, 4.0, 6.0}) {
            EXPECT_EQ(circles(storage, threshold,
                              MethodOptions(CirclesMethod::Sequential))
                          .members,
                      OracleSequential(n, condensed, threshold))
                << "n " << n << ", threshold " << threshold;
        }
    }
}

TEST(CirclesTest, EveryPairOfMembersIsBeyondTheThreshold) {
    const size_t n = 15;
    const std::vector<double> condensed = Scrambled(n);
    const DenseStorage storage = MakeStorage(n, condensed);
    for (const CirclesMethod method :
         {CirclesMethod::MaxMin, CirclesMethod::Sequential}) {
        for (const double threshold : {0.0, 1.0, 3.0, 5.0}) {
            const CirclesResult result =
                circles(storage, threshold, MethodOptions(method));
            EXPECT_EQ(result.count, result.members.size());
            for (size_t a = 0; a < result.members.size(); ++a) {
                for (size_t b = a + 1; b < result.members.size(); ++b) {
                    EXPECT_GT(CondensedAt(n, condensed, result.members[a],
                                          result.members[b]),
                              threshold);
                }
            }
        }
    }
}

TEST(CirclesTest, DuplicatesAtAZeroThresholdPackIntoOneCircle) {
    const DenseStorage storage = MakeStorage(4, std::vector<double>(6, 0.0));
    for (const CirclesMethod method :
         {CirclesMethod::MaxMin, CirclesMethod::Sequential}) {
        EXPECT_EQ(circles(storage, 0.0, MethodOptions(method)).members,
                  std::vector<size_t>({0}));
    }
}

TEST(CirclesTest, OneItemIsOneCircle) {
    const DenseStorage storage = MakeStorage(1, {});
    for (const CirclesMethod method :
         {CirclesMethod::MaxMin, CirclesMethod::Sequential}) {
        EXPECT_EQ(circles(storage, 0.5, MethodOptions(method)).count, 1u);
    }
}

TEST(CirclesTest, MaxMinRefusesANonFiniteDistanceItReads) {
    std::vector<double> condensed = Line(4);
    condensed[1] = NaN;  // (0, 2)
    const DenseStorage storage = MakeStorage(4, condensed);

    ExpectInvalidArgument(
        [&] { circles(storage, 0.5, MethodOptions(CirclesMethod::MaxMin)); },
        "Diversity selection read a non-finite distance between items 0 and 2");
}

// Candidate 2 is already rejected by member 0, and the NaN to member 1 still
// refuses: an answer must not depend on which member was compared first.
TEST(CirclesTest, SequentialRefusesANaNEvenAfterAnEarlyRejection) {
    std::vector<double> condensed = Positions({0, 5, 0.5});
    condensed[2] = NaN;  // (1, 2)
    const DenseStorage storage = MakeStorage(3, condensed);

    ExpectInvalidArgument(
        [&] { circles(storage, 1.0, MethodOptions(CirclesMethod::Sequential)); },
        "Diversity selection read a non-finite distance between items 1 and 2");
}

TEST(CirclesValidationTest, RefusesEachInvalidRequest) {
    const DenseStorage storage = MakeStorage(4, Line(4));

    ExpectInvalidArgument(
        [&] {
            circles(storage, 1.0, MethodOptions(static_cast<CirclesMethod>(7)));
        },
        "Unknown #Circles method");

    CirclesOptions zero_chunk;
    zero_chunk.chunk_size = 0;
    ExpectInvalidArgument([&] { circles(storage, 1.0, zero_chunk); },
                          "#Circles chunk_size must be at least one");

    for (const double threshold :
         {NaN, std::numeric_limits<double>::infinity(), -1.0}) {
        ExpectInvalidArgument(
            [&] { circles(storage, threshold, CirclesOptions()); },
            "#Circles threshold must be finite and non-negative");
    }
}

TEST(CirclesValidationTest, RefusesStorageItCannotRead) {
    ExpectInvalidArgument(
        [] { circles(SparseStorage(4, 0.5), 1.0, CirclesOptions()); },
        "#Circles requires complete pairwise distances; SparseStorage is not "
        "supported");
    ExpectInvalidArgument(
        [] { circles(NullDataStorage(4), 1.0, CirclesOptions()); },
        "#Circles requires contiguous dense or memory-mapped storage");
    ExpectInvalidArgument(
        [] { circles(DenseStorage(0), 1.0, CirclesOptions()); },
        "#Circles requires at least one item");
}

// Every packing is an independent set of the graph joining pairs at or within
// the threshold, so neither method can exceed the exhaustive maximum.
TEST(CirclesTest, NeitherMethodExceedsTheMaximumIndependentSet) {
    for (const size_t n : {size_t{4}, size_t{7}, size_t{10}}) {
        const std::vector<double> condensed = Scrambled(n);
        const DenseStorage storage = MakeStorage(n, condensed);
        for (const double threshold : {1.0, 2.5, 4.0}) {
            size_t maximum = 0;
            for (size_t subset = 1; subset < (size_t{1} << n); ++subset) {
                bool independent = true;
                size_t size = 0;
                for (size_t a = 0; a < n && independent; ++a) {
                    if ((subset >> a & 1u) == 0) {
                        continue;
                    }
                    ++size;
                    for (size_t b = a + 1; b < n; ++b) {
                        if ((subset >> b & 1u) != 0 &&
                            CondensedAt(n, condensed, a, b) <= threshold) {
                            independent = false;
                            break;
                        }
                    }
                }
                if (independent) {
                    maximum = std::max(maximum, size);
                }
            }
            for (const CirclesMethod method :
                 {CirclesMethod::MaxMin, CirclesMethod::Sequential}) {
                EXPECT_LE(circles(storage, threshold, MethodOptions(method)).count,
                          maximum)
                    << "n " << n << ", threshold " << threshold;
            }
        }
    }
}
