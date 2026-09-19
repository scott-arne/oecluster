#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <random>
#include <stdexcept>
#include <vector>

#include "oecluster/StorageBackend.h"
#include "oecluster/clustering/ClusterTypes.h"
#include "oecluster/clustering/KMedoids.h"

#include "../../src/clustering/MaxMinKernel.h"

using namespace OECluster;

namespace {

// Four items on a line at 0, 0.5, 10, 11: farthest-first from item 0 must
// reach the far pair before the near one.
//
// Item 1 sits at 0.5 rather than 1.0 so the third selection is decided by
// distance instead of by the tie rule. At 1.0, items 1 and 2 would both be 1.0
// from the selected set {0, 3} and the smaller-index rule would pick item 1 --
// a correct answer, but one that tests the tie rule rather than farthest-first
// coverage. At 0.5 item 2 wins outright at 1.0 against item 1's 0.5.
// Every position is dyadic, so every distance below is exact.
DenseStorage MakeLineStorage() {
    DenseStorage storage(4);
    storage.Set(0, 1, 0.5);
    storage.Set(0, 2, 10.0);
    storage.Set(0, 3, 11.0);
    storage.Set(1, 2, 9.5);
    storage.Set(1, 3, 10.5);
    storage.Set(2, 3, 1.0);
    return storage;
}

DenseStorage MakeAllZeroStorage(size_t n) {
    DenseStorage storage(n);
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i + 1; j < n; ++j) {
            storage.Set(i, j, 0.0);
        }
    }
    return storage;
}

// Four items at 0, 4, 9, 10: the pair (2, 3) is clustered together (distance 1.0)
// while item 1 is midway between them. This fixture separates the MaxMin recurrence
// from "farthest from the seed" — the nearest-distance update after selecting item 3
// lowers item 2's distance from 9.0 to 1.0, which causes item 1 (at 4.0) to be picked
// next instead of item 2. Without that update item 2 would be chosen. Every position
// is an integer, so every distance is exact.
DenseStorage MakeClusteredLineStorage() {
    DenseStorage storage(4);
    storage.Set(0, 1, 4.0);
    storage.Set(0, 2, 9.0);
    storage.Set(0, 3, 10.0);
    storage.Set(1, 2, 5.0);
    storage.Set(1, 3, 6.0);
    storage.Set(2, 3, 1.0);
    return storage;
}

// Two tight triples separated by a wide gap. Items 0-2 sit within 0.5 of one
// another around 0.0; items 3-5 sit within 0.5 of one another around 10.0.
// Every position is a dyadic rational, so every distance, partial sum and
// total is exactly representable: the interior items 2 and 3 tie at 30.0 on
// total distance with no rounding for the tie rule to be at the mercy of.
DenseStorage MakeTwoTriplesStorage() {
    const double positions[6] = {0.0, 0.25, 0.5, 10.0, 10.25, 10.5};
    DenseStorage storage(6);
    for (size_t i = 0; i < 6; ++i) {
        for (size_t j = i + 1; j < 6; ++j) {
            storage.Set(i, j, std::abs(positions[i] - positions[j]));
        }
    }
    return storage;
}

// A backend that is not SparseStorage and reports a nonempty pair count, yet
// hands back a null buffer. No production backend behaves this way, which is
// precisely why the second storage guard needs a test double: without one the
// guard is unreachable and a future backend could trip it unnoticed. Only the
// six pure virtuals are overridden; Finalize() keeps its base no-op.
class NullDataStorage : public StorageBackend {
public:
    explicit NullDataStorage(size_t n) : n_(n) {}

    void Set(size_t, size_t, double) override {}
    double Get(size_t, size_t) const override { return 0.0; }
    size_t NumSamples() const override { return n_; }
    size_t NumPairs() const override { return n_ * (n_ - 1) / 2; }
    double* Data() override { return nullptr; }
    const double* Data() const override { return nullptr; }

private:
    size_t n_;
};

// Three distinct items, each duplicated once, so three pairs sit at distance 0.
DenseStorage MakeDuplicateRowStorage() {
    const double positions[6] = {0.0, 0.0, 5.0, 5.0, 9.0, 9.0};
    DenseStorage storage(6);
    for (size_t i = 0; i < 6; ++i) {
        for (size_t j = i + 1; j < 6; ++j) {
            storage.Set(i, j, std::abs(positions[i] - positions[j]));
        }
    }
    return storage;
}

// Textbook PAM, carried by the test so the production loop has something
// independent to match: for every (medoid, non-medoid) pair, reassign every
// item from scratch and recompute the total. No caching of any kind.
double NaiveTotalCost(const double* data, size_t n,
                      const std::vector<size_t>& medoids) {
    double total = 0.0;
    for (size_t j = 0; j < n; ++j) {
        double nearest = std::numeric_limits<double>::infinity();
        for (const size_t medoid : medoids) {
            const double distance = detail::dense_distance(data, n, j, medoid);
            if (distance < nearest) {
                nearest = distance;
            }
        }
        total += nearest;
    }
    return total;
}

// The labels the medoid set implies, written out independently of the
// production assignment cache. ``medoids`` is ascending, as the production
// result's is, so "the first slot at the minimum" is the smaller-medoid-item
// tie rule and the two label vectors are comparable element by element.
std::vector<ClusterLabel> NaiveLabels(const double* data, size_t n,
                                      const std::vector<size_t>& medoids) {
    const size_t k = medoids.size();
    std::vector<ClusterLabel> labels(n, 0);
    for (size_t j = 0; j < n; ++j) {
        size_t chosen = k;
        double best = 0.0;
        for (size_t slot = 0; slot < k; ++slot) {
            // The self-assignment rule: a medoid belongs to its own slot even
            // when a duplicate medoid sits at distance zero from it.
            if (medoids[slot] == j) {
                chosen = slot;
                break;
            }
            const double distance =
                detail::dense_distance(data, n, j, medoids[slot]);
            if (chosen == k || distance < best) {
                chosen = slot;
                best = distance;
            }
        }
        labels[j] = static_cast<ClusterLabel>(chosen);
    }
    return labels;
}

std::vector<size_t> NaivePamSwapPhase(const double* data, size_t n,
                                      std::vector<size_t> medoids) {
    const size_t k = medoids.size();
    double current = NaiveTotalCost(data, n, medoids);

    for (size_t iteration = 0; iteration < 1000; ++iteration) {
        std::vector<bool> is_medoid(n, false);
        for (const size_t medoid : medoids) {
            is_medoid[medoid] = true;
        }

        bool found = false;
        double best_total = 0.0;
        size_t best_entering = 0;
        size_t best_leaving_item = 0;
        size_t best_slot = 0;

        for (size_t h = 0; h < n; ++h) {
            if (is_medoid[h]) {
                continue;
            }
            for (size_t slot = 0; slot < k; ++slot) {
                std::vector<size_t> trial = medoids;
                trial[slot] = h;
                const double total = NaiveTotalCost(data, n, trial);
                if (!(total < current)) {
                    continue;
                }
                // The same lexicographic key the production loop uses:
                // (score, entering item, leaving medoid item).
                const bool wins =
                    !found || total < best_total ||
                    (total == best_total &&
                     (h < best_entering ||
                      (h == best_entering && medoids[slot] < best_leaving_item)));
                if (wins) {
                    found = true;
                    best_total = total;
                    best_entering = h;
                    best_leaving_item = medoids[slot];
                    best_slot = slot;
                }
            }
        }

        if (!found) {
            break;
        }
        medoids[best_slot] = best_entering;
        current = best_total;
    }

    std::sort(medoids.begin(), medoids.end());
    return medoids;
}

// Euclidean distances over random 2-D points: metric and continuous.
DenseStorage MakeEuclideanStorage(size_t n, uint32_t seed, double scale) {
    std::mt19937 engine(seed);
    std::uniform_real_distribution<double> coordinate(0.0, 1.0);
    std::vector<double> x(n);
    std::vector<double> y(n);
    for (size_t i = 0; i < n; ++i) {
        x[i] = coordinate(engine);
        y[i] = coordinate(engine);
    }

    DenseStorage storage(n);
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i + 1; j < n; ++j) {
            const double dx = x[i] - x[j];
            const double dy = y[i] - y[j];
            storage.Set(i, j, scale * std::sqrt(dx * dx + dy * dy));
        }
    }
    return storage;
}

// Random symmetric values with a zeroed diagonal, violating the triangle
// inequality outright.
DenseStorage MakeNonMetricStorage(size_t n, uint32_t seed, double scale) {
    std::mt19937 engine(seed);
    std::uniform_real_distribution<double> value(0.0, 1.0);

    DenseStorage storage(n);
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i + 1; j < n; ++j) {
            storage.Set(i, j, scale * value(engine));
        }
    }
    return storage;
}

// The exactly representable family's metric half: Manhattan distances between
// integer lattice points are integers and satisfy the triangle inequality, so
// this one is a theorem *and* a metric. MakeSmallIntegerStorage below is the
// non-metric half -- unconstrained symmetric integers violate the triangle
// inequality freely -- and the spec asks for both.
DenseStorage MakeIntegerLatticeStorage(size_t n, uint32_t seed) {
    std::mt19937 engine(seed);
    std::uniform_int_distribution<int> coordinate(0, 20);
    std::vector<int> x(n);
    std::vector<int> y(n);
    for (size_t i = 0; i < n; ++i) {
        x[i] = coordinate(engine);
        y[i] = coordinate(engine);
    }

    DenseStorage storage(n);
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i + 1; j < n; ++j) {
            // Differenced as doubles so <cmath>'s std::abs applies; the values
            // are small integers, so every step is exact.
            const double distance =
                std::abs(static_cast<double>(x[i] - x[j])) +
                std::abs(static_cast<double>(y[i] - y[j]));
            storage.Set(i, j, distance);
        }
    }
    return storage;
}

// Small integers in a double: every distance, partial sum and total is exactly
// representable, so both cost expressions agree bit for bit and parity is a
// theorem rather than an observation. Unconstrained, so this is the non-metric
// instance of the family.
DenseStorage MakeSmallIntegerStorage(size_t n, uint32_t seed) {
    std::mt19937 engine(seed);
    std::uniform_int_distribution<int> value(0, 9);

    DenseStorage storage(n);
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i + 1; j < n; ++j) {
            storage.Set(i, j, static_cast<double>(value(engine)));
        }
    }
    return storage;
}

// Parity isolates the swap phase, so both sides must start from the same
// medoids: PAM's answer depends on its seeds, and comparing a BUILD-seeded run
// against a reference seeded some other way would test the initializer rather
// than the kernel. The seeds are therefore handed in explicitly, and the
// production call runs with init = Explicit so it performs no initialization
// of its own.
void ExpectParityFromSeeds(const DenseStorage& storage,
                           const std::vector<size_t>& seeds) {
    const size_t n = storage.NumSamples();
    const double* data = storage.Data();

    KMedoidsOptions options;
    options.n_clusters = seeds.size();
    options.init = KMedoidsInit::Explicit;
    options.initial_medoids = seeds;

    const KMedoidsResult result = k_medoids_cluster(storage, options);
    const std::vector<size_t> reference = NaivePamSwapPhase(data, n, seeds);

    EXPECT_EQ(result.Medoids(), reference);
    // Labels as well as medoids: two runs can agree on the medoid set and still
    // disagree on which cluster a tied item landed in, and the spec asserts all
    // three of medoids, labels and cost.
    EXPECT_EQ(result.Labels(), NaiveLabels(data, n, reference));
    // EXPECT_EQ, not EXPECT_DOUBLE_EQ: the latter accepts a four-ULP
    // difference, and the spec asks this test for a bit-identical cost. Both
    // sides sum the same summands in ascending item order, so the result is
    // exact and a tolerant comparison would only hide a real divergence in the
    // delta algebra. This is the one assertion in the suite that pins the
    // FastPAM1 arithmetic against a brute-force total; every other cost check
    // recomputes its expectation differently and stays on EXPECT_DOUBLE_EQ.
    EXPECT_EQ(result.Cost(), NaiveTotalCost(data, n, reference));
}

// Two independent seed sets per matrix: the first k items, and the test's own
// MaxMin selection. Neither borrows anything from the production initializers,
// so a bug in BUILD cannot hide a bug in the swap phase or vice versa.
void ExpectParityWithNaivePam(const DenseStorage& storage, size_t k) {
    std::vector<size_t> leading(k);
    for (size_t slot = 0; slot < k; ++slot) {
        leading[slot] = slot;
    }
    ExpectParityFromSeeds(storage, leading);

    ExpectParityFromSeeds(
        storage,
        detail::maxmin_select_from(storage.Data(), storage.NumSamples(), k, 0));
}

}  // namespace

TEST(MaxMinKernelTest, SelectsFarthestFirstFromTheSeed) {
    const DenseStorage storage = MakeLineStorage();

    const std::vector<size_t> selection =
        detail::maxmin_select_from(storage.Data(), storage.NumSamples(), 3, 0);

    EXPECT_EQ(selection, std::vector<size_t>({0, 3, 2}));
}

TEST(MaxMinKernelTest, LowersTheNearestDistanceAfterEachSelection) {
    const DenseStorage storage = MakeClusteredLineStorage();

    const std::vector<size_t> selection =
        detail::maxmin_select_from(storage.Data(), storage.NumSamples(), 3, 0);

    EXPECT_EQ(selection, std::vector<size_t>({0, 3, 1}));
}

TEST(MaxMinKernelTest, HonorsTheSeed) {
    const DenseStorage storage = MakeLineStorage();

    const std::vector<size_t> selection =
        detail::maxmin_select_from(storage.Data(), storage.NumSamples(), 2, 2);

    EXPECT_EQ(selection, std::vector<size_t>({2, 0}));
}

// On an all-zero matrix every candidate ties at every step. Without the
// selected mask the smaller-index rule would return the seed `count` times.
TEST(MaxMinKernelTest, ExcludesSelectedItemsWhenEveryDistanceTies) {
    const DenseStorage storage = MakeAllZeroStorage(5);

    const std::vector<size_t> selection =
        detail::maxmin_select_from(storage.Data(), storage.NumSamples(), 3, 2);

    EXPECT_EQ(selection, std::vector<size_t>({2, 0, 1}));
}

TEST(MaxMinKernelTest, SelectsEveryItemWhenCountEqualsTheItemCount) {
    const DenseStorage storage = MakeAllZeroStorage(3);

    const std::vector<size_t> selection =
        detail::maxmin_select_from(storage.Data(), storage.NumSamples(), 3, 1);

    EXPECT_EQ(selection, std::vector<size_t>({1, 0, 2}));
}

TEST(MaxMinKernelTest, RejectsAnOutOfRangeCount) {
    const DenseStorage storage = MakeLineStorage();

    EXPECT_THROW(
        detail::maxmin_select_from(storage.Data(), storage.NumSamples(), 0, 0),
        std::invalid_argument);
    EXPECT_THROW(
        detail::maxmin_select_from(storage.Data(), storage.NumSamples(), 5, 0),
        std::invalid_argument);
}

TEST(MaxMinKernelTest, RejectsAnOutOfRangeSeed) {
    const DenseStorage storage = MakeLineStorage();

    EXPECT_THROW(
        detail::maxmin_select_from(storage.Data(), storage.NumSamples(), 2, 4),
        std::invalid_argument);
}

TEST(KMedoidsValidationTest, RefusesSparseStorage) {
    SparseStorage storage(4, 0.5);
    KMedoidsOptions options;
    options.n_clusters = 2;

    try {
        k_medoids_cluster(storage, options);
        FAIL() << "expected std::invalid_argument";
    } catch (const std::invalid_argument& error) {
        // The shared text from detail::validate_complete_distance_storage,
        // named with this algorithm. Asserted rather than merely caught: the
        // message is what tells a caller which of the two storage failures
        // they hit.
        EXPECT_STREQ(error.what(),
                     "K-medoids clustering requires complete pairwise "
                     "distances; SparseStorage is not supported");
    }
}

// The second half of the first validation row: a backend that is not
// SparseStorage, reports pairs, and still hands back a null buffer. No
// production backend does this, so the test supplies one -- without it the row
// is only half covered and a caller would reach dense_distance() with nullptr.
TEST(KMedoidsValidationTest, RefusesAPairCountingBackendWithNoBuffer) {
    NullDataStorage storage(6);
    KMedoidsOptions options;
    options.n_clusters = 2;

    try {
        k_medoids_cluster(storage, options);
        FAIL() << "expected std::invalid_argument";
    } catch (const std::invalid_argument& error) {
        EXPECT_STREQ(error.what(),
                     "K-medoids clustering requires contiguous dense or "
                     "memory-mapped storage");
    }
}

TEST(KMedoidsValidationTest, RefusesAZeroChunkSize) {
    KMedoidsOptions options;
    options.n_clusters = 2;
    options.chunk_size = 0;

    try {
        k_medoids_cluster(MakeTwoTriplesStorage(), options);
        FAIL() << "expected std::invalid_argument";
    } catch (const std::invalid_argument& error) {
        EXPECT_STREQ(error.what(), "K-medoids chunk_size must be at least one");
    }
}

TEST(KMedoidsValidationTest, RefusesZeroMaxIterations) {
    KMedoidsOptions options;
    options.n_clusters = 2;
    options.max_iterations = 0;

    try {
        k_medoids_cluster(MakeTwoTriplesStorage(), options);
        FAIL() << "expected std::invalid_argument";
    } catch (const std::invalid_argument& error) {
        EXPECT_STREQ(error.what(),
                     "K-medoids max_iterations must be at least one");
    }
}

TEST(KMedoidsValidationTest, RefusesZeroClusters) {
    KMedoidsOptions options;
    options.n_clusters = 0;

    try {
        k_medoids_cluster(MakeTwoTriplesStorage(), options);
        FAIL() << "expected std::invalid_argument";
    } catch (const std::invalid_argument& error) {
        EXPECT_STREQ(error.what(), "K-medoids n_clusters must be at least one");
    }
}

TEST(KMedoidsValidationTest, RefusesMoreClustersThanItems) {
    KMedoidsOptions options;
    options.n_clusters = 7;

    try {
        k_medoids_cluster(MakeTwoTriplesStorage(), options);
        FAIL() << "expected std::invalid_argument";
    } catch (const std::invalid_argument& error) {
        EXPECT_STREQ(error.what(),
                     "K-medoids n_clusters must be at most the item count");
    }
}

TEST(KMedoidsValidationTest, RefusesSeedsWithoutExplicitInitialization) {
    KMedoidsOptions options;
    options.n_clusters = 2;
    options.init = KMedoidsInit::Build;
    options.initial_medoids = {0, 3};

    try {
        k_medoids_cluster(MakeTwoTriplesStorage(), options);
        FAIL() << "expected std::invalid_argument";
    } catch (const std::invalid_argument& error) {
        EXPECT_STREQ(
            error.what(),
            "K-medoids initial_medoids requires an explicit initialization");
    }
}

TEST(KMedoidsValidationTest, RefusesTheWrongNumberOfSeeds) {
    KMedoidsOptions options;
    options.n_clusters = 2;
    options.init = KMedoidsInit::Explicit;
    options.initial_medoids = {0};

    try {
        k_medoids_cluster(MakeTwoTriplesStorage(), options);
        FAIL() << "expected std::invalid_argument";
    } catch (const std::invalid_argument& error) {
        EXPECT_STREQ(
            error.what(),
            "K-medoids initial_medoids must hold exactly n_clusters indices");
    }
}

TEST(KMedoidsValidationTest, RefusesAnOutOfRangeSeed) {
    KMedoidsOptions options;
    options.n_clusters = 2;
    options.init = KMedoidsInit::Explicit;
    options.initial_medoids = {0, 6};

    try {
        k_medoids_cluster(MakeTwoTriplesStorage(), options);
        FAIL() << "expected std::out_of_range";
    } catch (const std::out_of_range& error) {
        EXPECT_STREQ(
            error.what(),
            "K-medoids initial_medoids index is outside the storage range");
    }

    // Verify that range is reported before uniqueness regardless of which
    // problem appears first in the list. The duplicate appears at positions 0
    // and 1, while the out-of-range index is at position 2.
    options.n_clusters = 3;
    options.initial_medoids = {3, 3, 9};

    try {
        k_medoids_cluster(MakeTwoTriplesStorage(), options);
        FAIL() << "expected std::out_of_range";
    } catch (const std::out_of_range& error) {
        EXPECT_STREQ(
            error.what(),
            "K-medoids initial_medoids index is outside the storage range");
    }
}

TEST(KMedoidsValidationTest, RefusesDuplicateSeeds) {
    KMedoidsOptions options;
    options.n_clusters = 2;
    options.init = KMedoidsInit::Explicit;
    options.initial_medoids = {3, 3};

    try {
        k_medoids_cluster(MakeTwoTriplesStorage(), options);
        FAIL() << "expected std::invalid_argument";
    } catch (const std::invalid_argument& error) {
        EXPECT_STREQ(error.what(), "K-medoids initial_medoids must be unique");
    }
}

TEST(KMedoidsValidationTest, RefusesAnUnknownInitializationMethod) {
    KMedoidsOptions options;
    options.n_clusters = 2;
    options.init = static_cast<KMedoidsInit>(99);

    try {
        k_medoids_cluster(MakeTwoTriplesStorage(), options);
        FAIL() << "expected std::invalid_argument";
    } catch (const std::invalid_argument& error) {
        EXPECT_STREQ(error.what(), "Unknown k-medoids initialization method");
    }
}

// An empty matrix needs no rule of its own: n_clusters >= 1 fails the item
// count bound and n_clusters == 0 fails the row above it.
TEST(KMedoidsValidationTest, RefusesAnEmptyMatrix) {
    DenseStorage storage(0);
    KMedoidsOptions options;
    options.n_clusters = 1;

    EXPECT_THROW(k_medoids_cluster(storage, options), std::invalid_argument);
}

TEST(KMedoidsDegenerateTest, ShortCircuitsWhenEveryItemIsAMedoid) {
    KMedoidsOptions options;
    options.n_clusters = 6;

    const KMedoidsResult result = k_medoids_cluster(MakeTwoTriplesStorage(), options);

    EXPECT_EQ(result.Medoids(), std::vector<size_t>({0, 1, 2, 3, 4, 5}));
    EXPECT_EQ(result.Labels(), std::vector<ClusterLabel>({0, 1, 2, 3, 4, 5}));
    EXPECT_EQ(result.NumClusters(), 6u);
    EXPECT_DOUBLE_EQ(result.Cost(), 0.0);
    EXPECT_EQ(result.NumIterations(), 0u);
    EXPECT_TRUE(result.Converged());
    EXPECT_EQ(result.Method(), "k_medoids");
}

TEST(KMedoidsInitializationTest, BothInitializersAgreeOnTheGlobalMedoid) {
    KMedoidsOptions build_options;
    build_options.n_clusters = 1;
    build_options.init = KMedoidsInit::Build;

    KMedoidsOptions maxmin_options;
    maxmin_options.n_clusters = 1;
    maxmin_options.init = KMedoidsInit::FarthestFirst;

    const KMedoidsResult from_build =
        k_medoids_cluster(MakeTwoTriplesStorage(), build_options);
    const KMedoidsResult from_maxmin =
        k_medoids_cluster(MakeTwoTriplesStorage(), maxmin_options);

    EXPECT_EQ(from_build.Medoids(), from_maxmin.Medoids());
    // Items 2 and 3 are the two interior points and tie on total distance;
    // the smaller index wins.
    EXPECT_EQ(from_build.Medoids(), std::vector<size_t>({2}));
}

// The k == 1 case above cannot tell the two initializers apart: BUILD never
// enters its gain loop and farthest-first never makes a MaxMin selection, so
// both collapse to the global medoid. At k == 2 they diverge, and the divergence
// is what pins each one's own rule.
//
// Both start at the global medoid, item 1 (items 1 and 2 tie at 20.5 total
// distance and the smaller index wins). BUILD then scores candidates by the
// total distance they remove: item 2 and item 3 both score exactly 19.0, so the
// smaller index wins and BUILD takes item 2. Farthest-first instead takes the
// item that is furthest from the selected set, which is item 3 at 10.5.
//
// Both answers cost 1.5, and both are global optima -- every pair straddling the
// gap costs 1.5 -- so the swap loop Task 3 adds cannot move either one and these
// expectations survive it.
TEST(KMedoidsInitializationTest, TheTwoInitializersDivergeOnTheSecondMedoid) {
    const DenseStorage storage = MakeLineStorage();

    KMedoidsOptions build_options;
    build_options.n_clusters = 2;
    build_options.init = KMedoidsInit::Build;

    KMedoidsOptions maxmin_options;
    maxmin_options.n_clusters = 2;
    maxmin_options.init = KMedoidsInit::FarthestFirst;

    const KMedoidsResult from_build = k_medoids_cluster(storage, build_options);
    const KMedoidsResult from_maxmin = k_medoids_cluster(storage, maxmin_options);

    EXPECT_EQ(from_build.Medoids(), std::vector<size_t>({1, 2}));
    EXPECT_EQ(from_maxmin.Medoids(), std::vector<size_t>({1, 3}));
    EXPECT_NE(from_build.Medoids(), from_maxmin.Medoids());

    // Equal cost, different medoids: the cost alone would not separate them.
    EXPECT_DOUBLE_EQ(from_build.Cost(), 1.5);
    EXPECT_DOUBLE_EQ(from_maxmin.Cost(), 1.5);
    EXPECT_EQ(from_build.Labels(), std::vector<ClusterLabel>({0, 0, 1, 1}));
    EXPECT_EQ(from_maxmin.Labels(), std::vector<ClusterLabel>({0, 0, 1, 1}));
}

// The all-zero matrix ties every gain, every MaxMin distance and every delta
// at once. It is the single input that fails if the selected mask is dropped
// from either initializer: the medoid list comes back with duplicates and one
// cluster is empty.
TEST(KMedoidsDegenerateTest, AllZeroMatrixStillProducesDistinctMedoids) {
    for (const KMedoidsInit init :
         {KMedoidsInit::Build, KMedoidsInit::FarthestFirst}) {
        for (const size_t k : {size_t{2}, size_t{3}, size_t{5}}) {
            KMedoidsOptions options;
            options.n_clusters = k;
            options.init = init;

            const KMedoidsResult result =
                k_medoids_cluster(MakeAllZeroStorage(8), options);

            std::vector<size_t> medoids = result.Medoids();
            ASSERT_EQ(medoids.size(), k);
            EXPECT_TRUE(std::is_sorted(medoids.begin(), medoids.end()));
            EXPECT_EQ(std::unique(medoids.begin(), medoids.end()) - medoids.begin(),
                      static_cast<long>(k));
            ASSERT_EQ(result.NumClusters(), k);
            for (const Cluster& cluster : result.Members()) {
                EXPECT_FALSE(cluster.empty());
            }
            EXPECT_DOUBLE_EQ(result.Cost(), 0.0);
        }
    }
}

// The nearest-medoid scan visits slots in order, so a tie is only resolved
// correctly if an equal distance fails to displace the incumbent. Inverting that
// comparison passes every other test in this file: no other test reads Labels()
// on an input where a non-medoid is equidistant from two medoids. On an all-zero
// matrix every non-medoid ties against every medoid, so each one must land in
// the slot of the smallest medoid index; the inverted rule would send them all
// to the largest instead.
//
// Note the scope: because assemble() sorts the medoids before assigning, medoid
// item order and slot order always agree at this boundary, so this pins the
// direction of the tie rule and not its key.
TEST(KMedoidsTieRuleTest, EquidistantItemsTakeTheSmallestMedoidIndex) {
    KMedoidsOptions two_options;
    two_options.n_clusters = 2;
    const KMedoidsResult two =
        k_medoids_cluster(MakeAllZeroStorage(8), two_options);
    EXPECT_EQ(two.Medoids(), std::vector<size_t>({0, 1}));
    EXPECT_EQ(two.Labels(),
              std::vector<ClusterLabel>({0, 1, 0, 0, 0, 0, 0, 0}));

    KMedoidsOptions three_options;
    three_options.n_clusters = 3;
    const KMedoidsResult three =
        k_medoids_cluster(MakeAllZeroStorage(8), three_options);
    EXPECT_EQ(three.Medoids(), std::vector<size_t>({0, 1, 2}));
    EXPECT_EQ(three.Labels(),
              std::vector<ClusterLabel>({0, 1, 2, 0, 0, 0, 0, 0}));
}

// Duplicates at distance 0 are what the self-assignment rule exists for:
// without it both duplicate medoids land in the smaller-index slot and the
// other cluster comes back empty.
//
// The seeds are chosen so the configuration is already globally optimal: with
// medoids {0, 1, 2, 4} on positions {0, 0, 5, 5, 9, 9} every item sits at
// distance 0 from its nearest medoid, so the total cost is 0 and no swap can
// lower it. That matters because Task 3 adds the swap loop and re-runs this
// file: a seed set the optimizer would improve on -- {0, 1} at k == 2 costs
// 28, and swapping item 1 for item 2 drops it to 8 -- would pass here and then
// fail the moment the swap phase lands, and the surviving result would no
// longer hold two duplicate medoids at all. Optimal seeds keep the assertion
// about duplicates rather than about the optimizer.
TEST(KMedoidsDegenerateTest, DuplicateMedoidsStillYieldExactlyKNonEmptyClusters) {
    KMedoidsOptions options;
    options.n_clusters = 4;
    options.init = KMedoidsInit::Explicit;
    options.initial_medoids = {0, 1, 2, 4};

    const KMedoidsResult result =
        k_medoids_cluster(MakeDuplicateRowStorage(), options);

    ASSERT_EQ(result.NumClusters(), 4u);
    for (const Cluster& cluster : result.Members()) {
        EXPECT_FALSE(cluster.empty());
    }
    // Items 0 and 1 are the duplicate pair, and both are medoids: without the
    // self-assignment rule item 1 would be absorbed into slot 0 at distance 0
    // and slot 1 would come back empty.
    EXPECT_EQ(result.Medoids(), std::vector<size_t>({0, 1, 2, 4}));
    EXPECT_DOUBLE_EQ(result.Cost(), 0.0);
}

TEST(KMedoidsAssemblyTest, MedoidsAreSortedAndLabelsFollowThem) {
    KMedoidsOptions options;
    options.n_clusters = 2;
    options.init = KMedoidsInit::Explicit;
    options.initial_medoids = {4, 1};

    const KMedoidsResult result =
        k_medoids_cluster(MakeTwoTriplesStorage(), options);

    EXPECT_EQ(result.Medoids(), std::vector<size_t>({1, 4}));
    EXPECT_EQ(result.Labels(), std::vector<ClusterLabel>({0, 0, 0, 1, 1, 1}));
    EXPECT_EQ(result.Members()[0], Cluster({0, 1, 2}));
    EXPECT_EQ(result.Members()[1], Cluster({3, 4, 5}));
}

TEST(KMedoidsAssemblyTest, CostMatchesTheReturnedAssignment) {
    KMedoidsOptions options;
    options.n_clusters = 2;
    options.init = KMedoidsInit::Explicit;
    options.initial_medoids = {1, 4};

    const DenseStorage storage = MakeTwoTriplesStorage();
    const KMedoidsResult result = k_medoids_cluster(storage, options);

    double expected = 0.0;
    for (size_t j = 0; j < storage.NumSamples(); ++j) {
        const size_t medoid = result.Medoids()[static_cast<size_t>(result.Labels()[j])];
        expected += detail::dense_distance(storage.Data(), storage.NumSamples(),
                                           j, medoid);
    }
    EXPECT_DOUBLE_EQ(result.Cost(), expected);
}

TEST(KMedoidsParityTest, MatchesNaivePamOnMetricContinuousMatrices) {
    for (const size_t n : {size_t{20}, size_t{31}, size_t{40}}) {
        for (const size_t k : {size_t{2}, size_t{3}, size_t{5}}) {
            ExpectParityWithNaivePam(
                MakeEuclideanStorage(n, static_cast<uint32_t>(n * 31 + k), 1.0),
                k);
        }
    }
}

// Random symmetric values violate the triangle inequality outright. Parity here
// is the test that a metric assumption has not been smuggled into the kernel:
// any pruning rule that appeals to the triangle inequality skips a candidate
// naive PAM evaluates, and the two answers part company.
TEST(KMedoidsParityTest, MatchesNaivePamOnNonMetricContinuousMatrices) {
    for (const size_t n : {size_t{20}, size_t{31}, size_t{40}}) {
        for (const size_t k : {size_t{2}, size_t{3}, size_t{5}}) {
            ExpectParityWithNaivePam(
                MakeNonMetricStorage(n, static_cast<uint32_t>(n * 17 + k), 1.0),
                k);
        }
    }
}

// Exactly representable: no rounding exists, so identical output is a theorem.
// Both instances of the family are covered -- lattice distances are metric,
// unconstrained symmetric integers are not.
TEST(KMedoidsParityTest, MatchesNaivePamOnSmallIntegerMetricMatrices) {
    for (const size_t n : {size_t{20}, size_t{31}}) {
        for (const size_t k : {size_t{2}, size_t{4}}) {
            ExpectParityWithNaivePam(
                MakeIntegerLatticeStorage(n, static_cast<uint32_t>(n * 13 + k)),
                k);
        }
    }
}

TEST(KMedoidsParityTest, MatchesNaivePamOnSmallIntegerNonMetricMatrices) {
    for (const size_t n : {size_t{20}, size_t{31}}) {
        for (const size_t k : {size_t{2}, size_t{4}}) {
            ExpectParityWithNaivePam(
                MakeSmallIntegerStorage(n, static_cast<uint32_t>(n * 7 + k)), k);
        }
    }
}

TEST(KMedoidsParityTest, MatchesNaivePamOnDuplicateRowsAndAllZeros) {
    ExpectParityWithNaivePam(MakeDuplicateRowStorage(), 2);
    ExpectParityWithNaivePam(MakeDuplicateRowStorage(), 3);
    ExpectParityWithNaivePam(MakeAllZeroStorage(8), 2);
    ExpectParityWithNaivePam(MakeAllZeroStorage(8), 3);
}

TEST(KMedoidsSwapTest, ReachesTheOptimumFromADeliberatelyBadSeed) {
    const DenseStorage storage = MakeTwoTriplesStorage();

    KMedoidsOptions bad_seed;
    bad_seed.n_clusters = 2;
    bad_seed.init = KMedoidsInit::Explicit;
    bad_seed.initial_medoids = {0, 1};

    KMedoidsOptions from_build;
    from_build.n_clusters = 2;

    const KMedoidsResult recovered = k_medoids_cluster(storage, bad_seed);
    const KMedoidsResult expected = k_medoids_cluster(storage, from_build);

    EXPECT_EQ(recovered.Medoids(), expected.Medoids());
    EXPECT_TRUE(recovered.Converged());
    EXPECT_GT(recovered.NumIterations(), 0u);
}

TEST(KMedoidsSwapTest, PerformsNoIterationsWhenSeededAtTheOptimum) {
    const DenseStorage storage = MakeTwoTriplesStorage();

    KMedoidsOptions options;
    options.n_clusters = 2;
    options.init = KMedoidsInit::Explicit;
    options.initial_medoids = {1, 4};

    const KMedoidsResult result = k_medoids_cluster(storage, options);

    EXPECT_EQ(result.Medoids(), std::vector<size_t>({1, 4}));
    EXPECT_EQ(result.NumIterations(), 0u);
    EXPECT_TRUE(result.Converged());
}
