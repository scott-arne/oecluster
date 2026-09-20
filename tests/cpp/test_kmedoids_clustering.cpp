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

#include "../../src/clustering/KMedoidsSwapKernel.h"
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

    // The seed list is what distinguishes one call from the next within a
    // single matrix, and these are EXPECT_ rather than ASSERT_ assertions, so
    // without it a regression prints unlabeled vector diffs and then keeps
    // going into unrelated cases. Traced here rather than in the test bodies so
    // every caller inherits it.
    testing::Message seed_trace;
    seed_trace << "seeds = {";
    for (size_t slot = 0; slot < seeds.size(); ++slot) {
        seed_trace << (slot == 0 ? "" : ", ") << seeds[slot];
    }
    seed_trace << "}";
    SCOPED_TRACE(seed_trace);

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
    //
    // How much the assertion proves depends on the matrix, which matters when
    // one day it fails. On the small-integer families every distance and
    // partial sum is exactly representable, so ranking candidates by the
    // predicted delta and ranking them by a brute-force total are the same
    // ordering: agreement is a theorem there, and a failure is a kernel bug
    // with nothing else to blame. On the continuous families it is a regression
    // pin instead. The two orderings agree at these seeds because no near-tie
    // happens to straddle a rounding boundary, not because they must, so a
    // different libm, a different floating-point contraction setting or a new
    // seed can part them with no code change at all. Diagnose a
    // continuous-family failure as toolchain drift first and a kernel bug
    // second; a small-integer-family failure admits no such excuse.
    //
    // The 2^53 fixtures are the exception on both counts. Their distances are
    // integers, but they are placed exactly where partial sums stop being
    // exactly representable, so agreement is not a theorem there either and the
    // toolchain-drift diagnosis applies -- see MakeCancellingSwapStorage.
    EXPECT_EQ(result.Cost(), NaiveTotalCost(data, n, reference));
}

// Two seed sets per matrix: the first k items, and a MaxMin selection. The
// second does call the production farthest-first kernel -- maxmin_select_from
// is exactly what farthest_first_initialize delegates to -- but that does not
// weaken the comparison, because the seed list is an *input* handed identically
// to the production call and to NaivePamSwapPhase rather than an oracle for
// either. A bug in maxmin_select_from moves where both sides start and still
// cannot make a broken swap phase look correct; only the two swap phases are
// under test. Two seed sets rather than one so the parity claim does not rest
// on a single starting configuration.
void ExpectParityWithNaivePam(const DenseStorage& storage, size_t k) {
    SCOPED_TRACE(testing::Message()
                 << "n = " << storage.NumSamples() << ", k = " << k);

    std::vector<size_t> leading(k);
    for (size_t slot = 0; slot < k; ++slot) {
        leading[slot] = slot;
    }
    ExpectParityFromSeeds(storage, leading);

    ExpectParityFromSeeds(
        storage,
        detail::maxmin_select_from(storage.Data(), storage.NumSamples(), k, 0));
}

// The KMedoidsOptions defaults, restated so a white-box call exercises the same
// chunking the public entry point does.
constexpr size_t KERNEL_THREADS = 0;
constexpr size_t KERNEL_CHUNK = 4096;

// Five items on a line at 0, 1, 2, 20, 21. The three-item group has an odd
// size, so its median is the single item 1 rather than a tied pair -- which is
// what lets the winning swap below be unique. An even-sized group on a line has
// two optimal medians and every such swap ties, which is useful for the tie-rule
// test and useless for a test that wants one winner.
DenseStorage MakeSkewedLineStorage() {
    const double positions[5] = {0.0, 1.0, 2.0, 20.0, 21.0};
    DenseStorage storage(5);
    for (size_t i = 0; i < 5; ++i) {
        for (size_t j = i + 1; j < 5; ++j) {
            storage.Set(i, j, std::abs(positions[i] - positions[j]));
        }
    }
    return storage;
}

// Five items on a line at 0, 1, 100, 106, 103. Items 0 and 1 are the intended
// medoid pair for this fixture and sit together at the origin; the far group is
// items 2, 3 and 4, and its unique median is deliberately item 4 -- the
// highest-index non-medoid. Any candidate loop that stops early, or that prunes
// by distance from the current medoids, misses precisely the winner here.
DenseStorage MakeFarMedianStorage() {
    const double positions[5] = {0.0, 1.0, 100.0, 106.0, 103.0};
    DenseStorage storage(5);
    for (size_t i = 0; i < 5; ++i) {
        for (size_t j = i + 1; j < 5; ++j) {
            storage.Set(i, j, std::abs(positions[i] - positions[j]));
        }
    }
    return storage;
}

// The mirror of MakeFarMedianStorage: five items on a line at 100, 0, 1, 99, 102,
// where the far group is items 0, 3 and 4 and its unique median is item 0 -- the
// LOWEST-index non-medoid under the intended seeding. MakeFarMedianStorage guards
// the end of the candidate range, since its winner is the highest-index item; this
// one guards the start, which a loop beginning one past its chunk would skip.
DenseStorage MakeLowIndexMedianStorage() {
    const double positions[5] = {100.0, 0.0, 1.0, 99.0, 102.0};
    DenseStorage storage(5);
    for (size_t i = 0; i < 5; ++i) {
        for (size_t j = i + 1; j < 5; ++j) {
            storage.Set(i, j, std::abs(positions[i] - positions[j]));
        }
    }
    return storage;
}

// MakeFarMedianStorage with its last two items exchanged, so the far group's
// unique median moves from item 4 to item 3. The two existing median fixtures
// both place their winner at an EVEN index, which leaves a candidate loop that
// stepped two at a time indistinguishable from one that stepped by one. Here
// such a loop never evaluates the winner at all.
DenseStorage MakeOddIndexMedianStorage() {
    const double positions[5] = {0.0, 1.0, 100.0, 103.0, 106.0};
    DenseStorage storage(5);
    for (size_t i = 0; i < 5; ++i) {
        for (size_t j = i + 1; j < 5; ++j) {
            storage.Set(i, j, std::abs(positions[i] - positions[j]));
        }
    }
    return storage;
}

constexpr double TWO_POW_53 = 9007199254740992.0;
constexpr double TWO_POW_54 = 18014398509481984.0;

// Five unit-spaced items and one item a full 2^53 away from all of them. The
// spacing is the point: ulp(2^53) is 2, so a unit addend lands exactly halfway
// between two representable doubles and rounds to even. Accumulating the small
// distances first and the large one last keeps them (they reach 2^53 + 4);
// accumulating the large one first absorbs every one of them (the sum stays at
// 2^53). No fixture built from small integers can tell those two orders apart,
// which is what this one is built to do: witness the ascending-order contract
// that total_cost and recomputed_total both promise. It is the deliberate
// witness rather than the only one -- MakeSpeculativeUndoStorage below turns
// out to witness the same contract on one of its trial configurations, as a
// side effect of the arithmetic it was built for.
DenseStorage MakeOrderSensitiveSumStorage() {
    const size_t n = 6;
    DenseStorage storage(n);
    for (size_t i = 0; i < n - 1; ++i) {
        for (size_t j = i + 1; j < n - 1; ++j) {
            storage.Set(i, j, 1.0);
        }
        storage.Set(i, n - 1, TWO_POW_53);
    }
    return storage;
}

// Seven items, laid out so the predicted delta and the from-scratch
// recomputation of the same swap disagree by exactly two units.
//
// Items 0 and 6 are the medoids and items 2, 3 and 4 are their near neighbours;
// item 5 is the entrant under test and item 1 is a distant item that pins the
// running sum at 2^53, where ulp is 2 and a unit addend rounds to even. Which
// of the unit-sized terms survive therefore depends on the order they are added
// in, and the two summations here -- a delta accumulated over per-item
// corrections, and a total accumulated over per-item minima -- add them in
// different orders. The 2^54 entries are sentinels that keep every candidate
// other than the one under test unattractive, sometimes by becoming those
// candidates' minima: entering item 2 while leaving slot 0, for instance,
// gives medoids {2, 6} and leaves item 1 paying 2^54, which is
// min(d(1, 2), d(1, 6)). They are load-bearing, not padding.
//
// Both fixtures below assume the default rounding mode and no reassociation of
// the kernel's sums. Neither accumulation contains a multiply, so FMA
// contraction cannot reach them, and the build sets no fast-math flag; but a
// failure here should still be diagnosed as a floating-point-flag change before
// it is diagnosed as a kernel bug.
//
// :param medoid_gap: d(0, 6). Moves the second nearest distance of the two
//     medoids, which is what decides whether the disagreement lands on the
//     favourable or the unfavourable side.
// :param entrant_distance: d(2, 5) and d(3, 5), the entrant's pull on the two
//     items it would capture.
DenseStorage MakeCancellingSwapStorage(double medoid_gap,
                                       double entrant_distance) {
    DenseStorage storage(7);
    storage.Set(0, 1, TWO_POW_53);
    storage.Set(0, 2, 4.0);
    storage.Set(0, 3, 1000.0);
    storage.Set(0, 4, 2.0);
    storage.Set(0, 5, 0.0);
    storage.Set(0, 6, medoid_gap);
    storage.Set(1, 2, TWO_POW_54);
    storage.Set(1, 3, TWO_POW_54);
    storage.Set(1, 4, TWO_POW_54);
    storage.Set(1, 5, 0.0);
    storage.Set(1, 6, TWO_POW_54);
    storage.Set(2, 3, TWO_POW_54);
    storage.Set(2, 4, TWO_POW_54);
    storage.Set(2, 5, entrant_distance);
    storage.Set(2, 6, 1000.0);
    storage.Set(3, 4, TWO_POW_54);
    storage.Set(3, 5, entrant_distance);
    storage.Set(3, 6, 4.0);
    storage.Set(4, 5, TWO_POW_54);
    storage.Set(4, 6, TWO_POW_53 + 102.0);
    storage.Set(5, 6, TWO_POW_54);
    return storage;
}

// From medoids {0, 6} the predicted delta for (enter 5, leave slot 1) is -2,
// but rebuilding the assignment from scratch returns the unchanged cost, so the
// speculative swap has to be undone.
DenseStorage MakeSpeculativeUndoStorage() {
    return MakeCancellingSwapStorage(TWO_POW_53 + 6.0, 1.0);
}

// The mirror image: from medoids {0, 6} no predicted delta is negative, yet
// (enter 5, leave slot 1) does lower the recomputed total by two, so only the
// verification pass can find it.
DenseStorage MakeMissedImprovementStorage() {
    return MakeCancellingSwapStorage(TWO_POW_53, 3.0);
}

// The three inputs the swap kernels take together, bundled so a white-box test
// states its starting configuration once instead of repeating the build ritual.
struct KernelState {
    std::vector<size_t> medoids;
    std::vector<size_t> slot_of;
    std::vector<detail::Assignment> assignments;
    double cost = 0.0;
};

// Built through the production helpers on purpose: a hand-rolled Assignment
// cache would be testing the test's own arithmetic, and the d2 field in
// particular is easy to get subtly wrong.
KernelState MakeKernelState(const DenseStorage& storage,
                            std::vector<size_t> medoids) {
    const size_t n = storage.NumSamples();

    KernelState state;
    state.medoids = std::move(medoids);
    state.slot_of.assign(n, n);
    detail::refresh_slot_map(state.slot_of, state.medoids);
    state.assignments = detail::build_assignments(
        storage.Data(), n, state.medoids, KERNEL_THREADS, KERNEL_CHUNK);
    state.cost = detail::total_cost(state.assignments);
    return state;
}

// The end-to-end statement that the fast kernel selects what an exhaustive
// recomputation would: its predicted delta is the recomputed difference, and no
// other pair recomputes to a smaller total. EXPECT_EQ on the delta is only
// legitimate on the exactly representable families -- see the long note in
// ExpectParityFromSeeds -- so this helper takes integer fixtures only.
void ExpectPredictedDeltaMatchesRecomputation(
    const DenseStorage& storage, const std::vector<size_t>& medoids) {
    const size_t n = storage.NumSamples();
    const double* data = storage.Data();

    SCOPED_TRACE(testing::Message() << "n = " << n << ", k = " << medoids.size());

    const KernelState state = MakeKernelState(storage, medoids);
    const detail::SwapCandidate winner = detail::best_predicted_swap(
        data, n, state.medoids, state.assignments, state.slot_of, KERNEL_THREADS,
        KERNEL_CHUNK);
    ASSERT_TRUE(winner.valid);

    const double winning_total = detail::recomputed_total(
        data, n, state.assignments, winner.leaving_slot, winner.entering_item);
    EXPECT_EQ(winner.score, winning_total - state.cost);

    for (size_t slot = 0; slot < state.medoids.size(); ++slot) {
        for (size_t h = 0; h < n; ++h) {
            if (state.slot_of[h] != n) {
                continue;
            }
            SCOPED_TRACE(testing::Message()
                         << "leaving slot = " << slot << ", entering item = " << h);
            EXPECT_GE(detail::recomputed_total(data, n, state.assignments, slot, h),
                      winning_total);
        }
    }
}

// The continuous counterpart of ExpectPredictedDeltaMatchesRecomputation. The
// predicted delta and a recomputed total sum the same quantities in different
// orders, so on a continuous matrix they may differ in the last bits and the
// bit-exact assertion above would measure the toolchain rather than the kernel.
// What survives is the selection: whichever pair the kernel returns recomputes
// to the exhaustive minimum, so a near-tie resolved the other way still passes.
void ExpectSelectedPairIsAnExhaustiveMinimum(
    const DenseStorage& storage, const std::vector<size_t>& medoids) {
    const size_t n = storage.NumSamples();
    const double* data = storage.Data();

    SCOPED_TRACE(testing::Message() << "n = " << n << ", k = " << medoids.size());

    const KernelState state = MakeKernelState(storage, medoids);
    const detail::SwapCandidate winner = detail::best_predicted_swap(
        data, n, state.medoids, state.assignments, state.slot_of, KERNEL_THREADS,
        KERNEL_CHUNK);
    ASSERT_TRUE(winner.valid);

    double exhaustive = std::numeric_limits<double>::infinity();
    for (size_t slot = 0; slot < state.medoids.size(); ++slot) {
        for (size_t h = 0; h < n; ++h) {
            if (state.slot_of[h] != n) {
                continue;
            }
            exhaustive = std::min(
                exhaustive,
                detail::recomputed_total(data, n, state.assignments, slot, h));
        }
    }

    EXPECT_DOUBLE_EQ(
        detail::recomputed_total(data, n, state.assignments, winner.leaving_slot,
                                 winner.entering_item),
        exhaustive);
}

// Every (leaving slot, entering item) pair on one fixture, checking that the
// cached shortcut and a full rebuild agree bit for bit.
void ExpectRecomputedTotalMatchesARebuild(const DenseStorage& storage,
                                          const std::vector<size_t>& medoids) {
    const size_t n = storage.NumSamples();
    const double* data = storage.Data();

    SCOPED_TRACE(testing::Message() << "n = " << n << ", k = " << medoids.size());

    const KernelState state = MakeKernelState(storage, medoids);

    for (size_t slot = 0; slot < state.medoids.size(); ++slot) {
        for (size_t h = 0; h < n; ++h) {
            // Skipping existing medoids is not a shortcut: a trial set holding
            // the same item twice has fewer than k distinct medoids, which is
            // not a configuration the swap phase can ever reach.
            if (state.slot_of[h] != n) {
                continue;
            }
            SCOPED_TRACE(testing::Message()
                         << "leaving slot = " << slot << ", entering item = " << h);

            std::vector<size_t> trial = state.medoids;
            trial[slot] = h;
            const double rebuilt = detail::total_cost(detail::build_assignments(
                data, n, trial, KERNEL_THREADS, KERNEL_CHUNK));

            EXPECT_EQ(detail::recomputed_total(data, n, state.assignments, slot, h),
                      rebuilt);
        }
    }
}

// No (leaving slot, entering item) pair may have a from-scratch recomputed
// total below the reported cost. This asserts the convergence contract
// directly rather than through a proxy.
void ExpectVerifiedLocalOptimum(const DenseStorage& storage,
                                const KMedoidsResult& result) {
    const size_t n = storage.NumSamples();
    const double* data = storage.Data();
    const std::vector<size_t>& medoids = result.Medoids();

    std::vector<bool> is_medoid(n, false);
    for (const size_t medoid : medoids) {
        is_medoid[medoid] = true;
    }

    for (size_t h = 0; h < n; ++h) {
        if (is_medoid[h]) {
            continue;
        }
        for (size_t slot = 0; slot < medoids.size(); ++slot) {
            std::vector<size_t> trial = medoids;
            trial[slot] = h;
            EXPECT_GE(NaiveTotalCost(data, n, trial), result.Cost())
                << "swap (slot " << slot << ", item " << h << ") lowers the cost";
        }
    }
}

DenseStorage PermuteStorage(const DenseStorage& storage,
                            const std::vector<size_t>& order) {
    const size_t n = storage.NumSamples();
    DenseStorage permuted(n);
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i + 1; j < n; ++j) {
            permuted.Set(i, j,
                         detail::dense_distance(storage.Data(), n,
                                                order[i], order[j]));
        }
    }
    return permuted;
}

// The whole thread-count by chunk-size cross product against a single-threaded
// baseline. Held as a helper so more than one matrix can be put through it: the
// kernel's reductions are documented as order-fixed, and how loudly a schedule
// that broke that would show up depends on the matrix.
void ExpectIdenticalAcrossSchedules(const DenseStorage& storage, size_t k) {
    SCOPED_TRACE(testing::Message()
                 << "n = " << storage.NumSamples() << ", k = " << k);

    KMedoidsOptions baseline_options;
    baseline_options.n_clusters = k;
    baseline_options.num_threads = 1;
    baseline_options.chunk_size = 4096;
    const KMedoidsResult baseline = k_medoids_cluster(storage, baseline_options);

    // SIZE_MAX is a legal chunk size -- validation only rejects zero -- and it
    // is the value that overflows an unguarded (n + chunk_size - 1) ceiling to
    // zero chunks. A scan that runs no chunk at all returns its
    // default-initialized output: a single cluster at cost 0, reported as
    // converged. This row is what fails if effective_chunk_size() is dropped.
    for (const size_t threads : {size_t{1}, size_t{2}, size_t{4}, size_t{8}}) {
        for (const size_t chunk : {size_t{1}, size_t{7}, size_t{4096},
                                   std::numeric_limits<size_t>::max()}) {
            SCOPED_TRACE(testing::Message() << "num_threads = " << threads
                                            << ", chunk_size = " << chunk);
            KMedoidsOptions options;
            options.n_clusters = k;
            options.num_threads = threads;
            options.chunk_size = chunk;

            const KMedoidsResult result = k_medoids_cluster(storage, options);

            EXPECT_EQ(result.Medoids(), baseline.Medoids());
            EXPECT_EQ(result.Labels(), baseline.Labels());
            EXPECT_EQ(result.NumIterations(), baseline.NumIterations());
            EXPECT_EQ(result.Converged(), baseline.Converged());
            // Bit-identical, not merely close.
            EXPECT_EQ(result.Cost(), baseline.Cost());
        }
    }
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

// This row pins the RESULT contract at k == n, not the short circuit that
// produces it. The n_clusters == n branch in k_medoids_cluster is a performance
// optimization whose output is bit-identical to what the general path returns:
// with it deleted, build_initialize masks already-selected items and so
// necessarily ends holding all n of them, both swap kernels then see an empty
// candidate set and converge at zero iterations, and assemble sorts the medoids
// to the same identity partition at the same cost. The only difference is O(n^3)
// work against O(1), and no observable the public API exposes can see it -- the
// kernels take storage.Data() as a raw pointer, so a counting StorageBackend
// cannot measure distance reads either. Measured: with the short circuit
// removed the whole suite still passes. Treat that as known, not as a gap
// this file can close; catching it needs a complexity guard rather than a
// behavioral assertion.
TEST(KMedoidsDegenerateTest, ReturnsTheIdentityPartitionWhenEveryItemIsAMedoid) {
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

// The rows below call the swap kernels directly. k_medoids_cluster runs the
// fast prediction and the verification pass one after the other, and either
// alone reproduces the naive reference, so an end-to-end test cannot tell
// whether the fast kernel ran at all -- disabling it costs a factor of k in
// speed and nothing in output. Reaching the kernels from here is what makes
// that observable.

TEST(KMedoidsSwapKernelTest, PredictsTheExactDeltaForTheWinningSwap) {
    const DenseStorage storage = MakeSkewedLineStorage();
    const KernelState state = MakeKernelState(storage, {0, 3});

    // Items sit at 0, 1, 2, 20, 21 and the medoids are items 0 and 3, so every
    // item pays its distance to the nearer of positions 0 and 20:
    // 0 + 1 + 2 + 0 + 1 = 4.
    EXPECT_EQ(state.cost, 4.0);

    const detail::SwapCandidate candidate = detail::best_predicted_swap(
        storage.Data(), 5, state.medoids, state.assignments, state.slot_of,
        KERNEL_THREADS, KERNEL_CHUNK);

    ASSERT_TRUE(candidate.valid);
    EXPECT_EQ(candidate.entering_item, 1u);
    EXPECT_EQ(candidate.leaving_slot, 0u);
    EXPECT_EQ(candidate.leaving_item, 0u);

    // Derived from the definition of the swap rather than from the kernel's
    // algebra. Moving slot 0 from item 0 to item 1 puts the medoids at
    // positions 1 and 20, so the five items pay 1 + 0 + 1 + 0 + 1 = 3 against
    // the 4 above, a delta of -1. The five rival pairs are worse by hand too:
    // {2, 20} and {0, 21} both total 4 (delta 0), {0, 2} totals 38, {0, 1}
    // totals 40 and {21, 20} totals 57, so the winner is unique.
    EXPECT_EQ(candidate.score, -1.0);
}

TEST(KMedoidsSwapKernelTest, PredictedDeltaEqualsTheExactRecomputedDifference) {
    ExpectPredictedDeltaMatchesRecomputation(MakeSkewedLineStorage(), {0, 3});
    ExpectPredictedDeltaMatchesRecomputation(MakeFarMedianStorage(), {0, 1});
    ExpectPredictedDeltaMatchesRecomputation(MakeSmallIntegerStorage(12, 3),
                                             {0, 1, 2});
    ExpectPredictedDeltaMatchesRecomputation(MakeIntegerLatticeStorage(12, 5),
                                             {0, 1, 2});

    const DenseStorage wide_integers = MakeSmallIntegerStorage(20, 9);
    ExpectPredictedDeltaMatchesRecomputation(
        wide_integers, detail::maxmin_select_from(wide_integers.Data(), 20, 4, 0));

    const DenseStorage wide_lattice = MakeIntegerLatticeStorage(20, 21);
    ExpectPredictedDeltaMatchesRecomputation(
        wide_lattice, detail::maxmin_select_from(wide_lattice.Data(), 20, 4, 0));

    // Continuous matrices get the selection claim without the bit-exact delta.
    ExpectSelectedPairIsAnExhaustiveMinimum(MakeEuclideanStorage(16, 4, 1.0),
                                            {0, 1, 2});
    ExpectSelectedPairIsAnExhaustiveMinimum(MakeNonMetricStorage(16, 8, 1.0),
                                            {0, 1, 2});
}

TEST(KMedoidsSwapKernelTest, ScoresEveryNonMedoidCandidate) {
    const DenseStorage storage = MakeFarMedianStorage();
    const KernelState state = MakeKernelState(storage, {0, 1});

    // Item 4 is the last item the candidate loop reaches and is also the only
    // winner: entering it costs 7 against 10 for items 2 and 3. Chunk sizes 1
    // and 2 split the candidate range across several chunks, so the answer also
    // has to survive the serial reduction over chunk winners.
    for (const size_t chunk : {size_t{1}, size_t{2}, KERNEL_CHUNK}) {
        SCOPED_TRACE(testing::Message() << "chunk_size = " << chunk);
        const detail::SwapCandidate candidate = detail::best_predicted_swap(
            storage.Data(), 5, state.medoids, state.assignments, state.slot_of,
            KERNEL_THREADS, chunk);

        ASSERT_TRUE(candidate.valid);
        EXPECT_EQ(candidate.entering_item, 4u);
        EXPECT_EQ(state.slot_of[candidate.entering_item], 5u);
    }

    // The entering item is never one of the current medoids, on fixtures large
    // enough that an off-by-one in the candidate filter could hide.
    const DenseStorage integers = MakeSmallIntegerStorage(20, 9);
    const std::vector<size_t> seeds =
        detail::maxmin_select_from(integers.Data(), 20, 4, 0);
    const KernelState wide = MakeKernelState(integers, seeds);
    const detail::SwapCandidate candidate = detail::best_predicted_swap(
        integers.Data(), 20, wide.medoids, wide.assignments, wide.slot_of,
        KERNEL_THREADS, KERNEL_CHUNK);

    ASSERT_TRUE(candidate.valid);
    EXPECT_EQ(wide.slot_of[candidate.entering_item], 20u);
}

TEST(KMedoidsSwapKernelTest, ResolvesTiedDeltasToTheSmallestEnteringThenLeavingItem) {
    {
        // Two symmetric triples: shifting either medoid one step inward saves
        // exactly the same 0.25, so entering item 1 and entering item 4 tie and
        // the smaller entering item has to win.
        const DenseStorage storage = MakeTwoTriplesStorage();
        const KernelState state = MakeKernelState(storage, {0, 3});
        const double* data = storage.Data();

        EXPECT_EQ(detail::recomputed_total(data, 6, state.assignments, 0, 1),
                  detail::recomputed_total(data, 6, state.assignments, 1, 4));

        const detail::SwapCandidate candidate = detail::best_predicted_swap(
            data, 6, state.medoids, state.assignments, state.slot_of,
            KERNEL_THREADS, KERNEL_CHUNK);

        ASSERT_TRUE(candidate.valid);
        EXPECT_EQ(candidate.score, -0.25);
        EXPECT_EQ(candidate.entering_item, 1u);
        EXPECT_EQ(candidate.leaving_slot, 0u);
    }

    {
        // Both medoids sit at the origin, so entering item 4 saves the same
        // amount whichever of them leaves and the leaving key decides. The
        // medoid list is descending on purpose: slot 0 holds item 1 and slot 1
        // holds item 0, so a tie rule written on the slot index would answer
        // slot 0 while the required item rule answers slot 1.
        const DenseStorage storage = MakeFarMedianStorage();
        const KernelState state = MakeKernelState(storage, {1, 0});
        const double* data = storage.Data();

        EXPECT_EQ(detail::recomputed_total(data, 5, state.assignments, 0, 4),
                  detail::recomputed_total(data, 5, state.assignments, 1, 4));

        const detail::SwapCandidate candidate = detail::best_predicted_swap(
            data, 5, state.medoids, state.assignments, state.slot_of,
            KERNEL_THREADS, KERNEL_CHUNK);

        ASSERT_TRUE(candidate.valid);
        EXPECT_EQ(candidate.entering_item, 4u);
        EXPECT_EQ(candidate.leaving_item, 0u);
        EXPECT_EQ(candidate.leaving_slot, 1u);
    }
}

TEST(KMedoidsSwapKernelTest, VerificationPassSelectsTheExhaustiveMinimumTotal) {
    {
        // The medoid list ascending, which fixes the selected total and the
        // winning entering item.
        const DenseStorage storage = MakeFarMedianStorage();
        const KernelState state = MakeKernelState(storage, {0, 1});

        // Both medoids sit at the origin while items 2, 3 and 4 sit near 100, so
        // the starting total is 0 + 0 + 99 + 105 + 102 = 306.
        EXPECT_EQ(state.cost, 306.0);

        const detail::SwapCandidate verified = detail::verification_pass(
            storage.Data(), 5, state.medoids, state.assignments, state.slot_of,
            state.cost, KERNEL_THREADS, KERNEL_CHUNK);

        ASSERT_TRUE(verified.valid);
        EXPECT_EQ(verified.entering_item, 4u);
        EXPECT_EQ(verified.leaving_slot, 0u);
        EXPECT_EQ(verified.leaving_item, 0u);

        // This score is a TOTAL, not a delta: verification_pass ranks candidates by
        // the recomputed cost of the whole configuration, which is what makes its
        // comparison against the current cost bit-exact. Do not "fix" it to -299.
        // Moving item 0 out for item 4 leaves medoids at positions 1 and 103, so
        // the five items pay 1 + 0 + 3 + 3 + 0 = 7. Entering item 4 at the other
        // slot totals 7 as well, and the four remaining pairs all total 10, so 7 is
        // the exhaustive minimum and the smaller leaving item breaks the tie.
        EXPECT_EQ(verified.score, 7.0);
    }

    {
        // The same fixture and the same medoid set, listed descending: slot 0
        // holds item 1 and slot 1 holds item 0. Entering item 4 recomputes to 7
        // whichever medoid leaves, so the leaving key decides, and a pass that
        // ranked or recorded the leaving SLOT would answer slot 0 where the
        // required item rule answers slot 1. The block above cannot tell those
        // two apart, because {0, 1} makes every item index equal its slot index.
        const DenseStorage storage = MakeFarMedianStorage();
        const KernelState state = MakeKernelState(storage, {1, 0});

        EXPECT_EQ(state.cost, 306.0);

        const detail::SwapCandidate verified = detail::verification_pass(
            storage.Data(), 5, state.medoids, state.assignments, state.slot_of,
            state.cost, KERNEL_THREADS, KERNEL_CHUNK);

        ASSERT_TRUE(verified.valid);
        EXPECT_EQ(verified.entering_item, 4u);
        EXPECT_EQ(verified.leaving_item, 0u);
        EXPECT_EQ(verified.leaving_slot, 1u);
        EXPECT_EQ(verified.score, 7.0);
    }
}

TEST(KMedoidsSwapKernelTest, VerificationPassConsidersTheLowestIndexedItemAsAnEntrant) {
    // Both medoids sit at the origin while items 0, 3 and 4 sit near 100, so the
    // starting total is 99 + 0 + 0 + 98 + 101 = 298.
    const DenseStorage storage = MakeLowIndexMedianStorage();
    const KernelState state = MakeKernelState(storage, {1, 2});

    EXPECT_EQ(state.cost, 298.0);

    const detail::SwapCandidate verified = detail::verification_pass(
        storage.Data(), 5, state.medoids, state.assignments, state.slot_of,
        state.cost, KERNEL_THREADS, KERNEL_CHUNK);

    // Item 0 is the unique best entrant, and it is item zero on purpose: a
    // candidate loop that began one past the start of its chunk would skip it,
    // answer item 3 at a total of 5, and still look like a successful
    // verification. Every other verification-pass fixture seeds item 0 as a
    // medoid, so this row is the only one where that mistake is visible.
    ASSERT_TRUE(verified.valid);
    EXPECT_EQ(verified.entering_item, 0u);
    EXPECT_EQ(verified.score, 4.0);

    // Removing either near medoid leaves the other one item away from its
    // partner, so both slots tie at 4 and the smaller leaving ITEM decides.
    EXPECT_EQ(verified.leaving_item, 1u);
    EXPECT_EQ(verified.leaving_slot, 0u);
}

TEST(KMedoidsSwapKernelTest, VerificationPassConsidersAnOddIndexedEntrantInEveryChunk) {
    // Medoids at 0 and 1 with items 2, 3 and 4 near 100, so the starting total
    // is 0 + 0 + 99 + 102 + 105 = 306.
    const DenseStorage storage = MakeOddIndexMedianStorage();
    const KernelState state = MakeKernelState(storage, {0, 1});

    EXPECT_EQ(state.cost, 306.0);

    // This loop carries two unrelated claims, so prune it with care.
    //
    // Per-chunk staging: num_threads = 1 fixes the order chunks are claimed in,
    // so with one worker the chunks run in ascending order and a pass that wrote
    // every chunk's local winner to the same slot would report the LAST chunk's
    // winner rather than the best one. Chunk sizes 1 and 2 put the winner in
    // chunk 3 of 5 and chunk 1 of 3 respectively -- never in chunk zero, and
    // never in the final chunk. KERNEL_CHUNK is a single chunk and says nothing
    // about staging.
    //
    // Candidate stride: chunk size 1 holds one candidate per chunk, so a loop
    // stepping two at a time is invisible there; only chunk size 2 and
    // KERNEL_CHUNK witness it.
    for (const size_t chunk : {size_t{1}, size_t{2}, KERNEL_CHUNK}) {
        SCOPED_TRACE(testing::Message() << "chunk_size = " << chunk);

        const detail::SwapCandidate verified = detail::verification_pass(
            storage.Data(), 5, state.medoids, state.assignments, state.slot_of,
            state.cost, 1, chunk);

        // Item 3 is the unique best entrant, and its index is odd on purpose:
        // both other median fixtures put their winner at an even index, so a
        // candidate loop advancing two at a time would answer item 2 at a total
        // of 10 on all of them and still look like a successful verification.
        ASSERT_TRUE(verified.valid);
        EXPECT_EQ(verified.entering_item, 3u);
        EXPECT_EQ(verified.score, 7.0);
        EXPECT_EQ(verified.leaving_item, 0u);
        EXPECT_EQ(verified.leaving_slot, 0u);
    }
}

TEST(KMedoidsSwapKernelTest, VerificationPassReturnsNoCandidateAtALocalOptimum) {
    const DenseStorage storage = MakeFarMedianStorage();
    const KernelState state = MakeKernelState(storage, {0, 4});
    const double* data = storage.Data();

    // {0, 4} is the exhaustive optimum over all ten pairs of these five items.
    EXPECT_EQ(state.cost, 7.0);

    // {1, 4} recomputes to exactly 7 as well, which is what makes the strict
    // inequality observable: accepting an equal-cost swap here would send the
    // loop between two configurations of the same cost forever.
    EXPECT_EQ(detail::recomputed_total(data, 5, state.assignments, 0, 1), 7.0);

    const detail::SwapCandidate verified = detail::verification_pass(
        data, 5, state.medoids, state.assignments, state.slot_of, state.cost,
        KERNEL_THREADS, KERNEL_CHUNK);

    EXPECT_FALSE(verified.valid);
}

TEST(KMedoidsSwapKernelTest, RecomputedTotalIsBitIdenticalToARebuiltAssignmentCost) {
    // The speculative undo in the optimization loop and the convergence
    // argument both rest on this: the cached shortcut must agree with a full
    // rebuild to the bit, or an accepted swap could report a cost the returned
    // labels do not produce.
    ExpectRecomputedTotalMatchesARebuild(MakeFarMedianStorage(), {0, 1});
    ExpectRecomputedTotalMatchesARebuild(MakeSmallIntegerStorage(12, 31),
                                         {0, 1, 2});
    ExpectRecomputedTotalMatchesARebuild(MakeIntegerLatticeStorage(12, 17),
                                         {0, 1, 2, 3});

    const DenseStorage integers = MakeSmallIntegerStorage(14, 47);
    ExpectRecomputedTotalMatchesARebuild(
        integers, detail::maxmin_select_from(integers.Data(), 14, 4, 0));

    // Degenerate matrices exercise the self-assignment rule inside the rebuild:
    // a medoid keeps its own slot even when a duplicate sits at distance zero,
    // and the shortcut has to reach the same totals anyway.
    ExpectRecomputedTotalMatchesARebuild(MakeDuplicateRowStorage(), {0, 2});
    ExpectRecomputedTotalMatchesARebuild(MakeAllZeroStorage(6), {0, 1});

    // The two fixtures below are the ones whose sums are order-sensitive here.
    // Everything above is small integers or a degenerate matrix, exact in any
    // accumulation order, and so cannot witness the ascending-order contract at
    // all. Measured by reversing recomputed_total's accumulation loop: the
    // focused suite then loses exactly two rows, this one and
    // KMedoidsDynamicRangeTest.HoldsParityWhenAPredictedImprovementIsUndone,
    // and inside this one only the two fixtures below report a pair.
    //
    // MakeOrderSensitiveSumStorage is built for that and fails six pairs under
    // the reversal, each reporting 2^53 from the shortcut against 2^53 + 4 from
    // the rebuild.
    ExpectRecomputedTotalMatchesARebuild(MakeOrderSensitiveSumStorage(), {0, 1});

    // The undo fixture fails exactly one pair -- leaving slot 1, entering item
    // 5, reporting 2^53 + 8 against a rebuilt 2^53 + 10 -- and that pair is the
    // swap the speculative undo branch turns on. Note which sum is the delicate
    // one: at the base medoids {0, 6} the per-item minima are 0, 2^53, 4, 4, 2,
    // 0, 0, every small term even against ulp(2^53) = 2, so the base total is
    // exact in any order and this helper never asserts it anyway. At the trial
    // set {0, 5} the minima are 0, 0, 1, 1, 2, 0, 2^53 + 6, two odd units
    // against the 2^53 term, and ascending and reversed genuinely differ. That
    // is the total the helper compares, and pinning that the shortcut and the
    // rebuild agree on it is what lets
    // KMedoidsDynamicRangeTest.HoldsParityWhenAPredictedImprovementIsUndone
    // state its precondition through the shortcut instead of through the rebuild
    // the swap loop actually performs.
    ExpectRecomputedTotalMatchesARebuild(MakeSpeculativeUndoStorage(), {0, 6});
}

TEST(KMedoidsDeterminismTest, IsIdenticalAcrossThreadCountsAndChunkSizes) {
    ExpectIdenticalAcrossSchedules(MakeEuclideanStorage(40, 4242, 1.0), 4);

    // A continuous matrix on its own understates what this row claims. Its
    // distances are all of one magnitude, so a reduction whose order depended on
    // the schedule would move the answer by a few last bits and might not move
    // the selected medoids at all.
    //
    // k = 1 on MakeOrderSensitiveSumStorage is the one configuration of that
    // fixture where the 2^53 term reaches the reported cost: the single medoid
    // is item 0, item 5 stays a non-medoid, and the total is four unit-sized
    // terms plus 2^53. Ascending that sums to 2^53 + 4; adding the far term
    // first absorbs all four units and leaves 2^53, so a reduction that summed
    // in a schedule-dependent order would report different whole-unit costs at
    // different chunk sizes rather than differing in the last bits. At any
    // larger k the far item is itself selected as a medoid and the cost is a sum
    // of small integers, exact in every order and blind to this.
    //
    // A reduction that is merely in the wrong FIXED order moves every run in
    // this cross product together and cannot be caught here at all; that is
    // pinned by the EXPECT_EQ(predicted.score, ...) assertions in the two
    // KMedoidsDynamicRangeTest rows.
    ExpectIdenticalAcrossSchedules(MakeOrderSensitiveSumStorage(), 1);
}

// Continuous distances only: a tied input is legitimately permutation-variant
// under any index-based tie rule, so asserting invariance there would assert
// something false.
TEST(KMedoidsPermutationTest, TheMedoidSetSurvivesARowPermutation) {
    const DenseStorage storage = MakeEuclideanStorage(30, 909, 1.0);
    const std::vector<size_t> order = {
        7,  22, 3,  18, 11, 0,  25, 14, 5,  29,
        1,  16, 9,  23, 12, 27, 4,  19, 8,  21,
        2,  26, 15, 6,  28, 10, 24, 13, 20, 17};

    KMedoidsOptions options;
    options.n_clusters = 4;
    options.init = KMedoidsInit::FarthestFirst;

    const KMedoidsResult original = k_medoids_cluster(storage, options);
    const KMedoidsResult permuted =
        k_medoids_cluster(PermuteStorage(storage, order), options);

    std::vector<size_t> mapped;
    mapped.reserve(permuted.Medoids().size());
    for (const size_t medoid : permuted.Medoids()) {
        mapped.push_back(order[medoid]);
    }
    std::sort(mapped.begin(), mapped.end());

    EXPECT_EQ(mapped, original.Medoids());
}

// PAM guarantees only a local optimum in general, so this targets data where
// local and global coincide.
TEST(KMedoidsGlobalOptimumTest, MatchesExhaustiveSearchOnWellSeparatedData) {
    const DenseStorage storage = MakeTwoTriplesStorage();
    const size_t n = storage.NumSamples();
    const double* data = storage.Data();

    KMedoidsOptions options;
    options.n_clusters = 2;
    const KMedoidsResult result = k_medoids_cluster(storage, options);

    double best = std::numeric_limits<double>::infinity();
    for (size_t a = 0; a < n; ++a) {
        for (size_t b = a + 1; b < n; ++b) {
            const double total = NaiveTotalCost(data, n, {a, b});
            if (total < best) {
                best = total;
            }
        }
    }

    EXPECT_DOUBLE_EQ(result.Cost(), best);
}

// The expected medoids are written literally, so the suite cannot be
// self-consistently wrong.
TEST(KMedoidsKnownAnswerTest, FindsTheMedoidOfEachTightTriple) {
    KMedoidsOptions options;
    options.n_clusters = 2;

    const KMedoidsResult result = k_medoids_cluster(MakeTwoTriplesStorage(), options);

    EXPECT_EQ(result.Medoids(), std::vector<size_t>({1, 4}));
    EXPECT_EQ(result.Labels(), std::vector<ClusterLabel>({0, 0, 0, 1, 1, 1}));
    // 0.25 + 0 + 0.25 per triple, and every term is exactly representable.
    EXPECT_DOUBLE_EQ(result.Cost(), 1.0);
    EXPECT_TRUE(result.Converged());
}

TEST(KMedoidsCostTest, CostMatchesAnIndependentSumOverTheReturnedLabels) {
    for (const uint32_t seed : {uint32_t{11}, uint32_t{12}, uint32_t{13}}) {
        const DenseStorage storage = MakeNonMetricStorage(25, seed, 1.0);
        KMedoidsOptions options;
        options.n_clusters = 3;

        const KMedoidsResult result = k_medoids_cluster(storage, options);

        double expected = 0.0;
        for (size_t j = 0; j < storage.NumSamples(); ++j) {
            const size_t medoid =
                result.Medoids()[static_cast<size_t>(result.Labels()[j])];
            expected += detail::dense_distance(
                storage.Data(), storage.NumSamples(), j, medoid);
        }
        EXPECT_DOUBLE_EQ(result.Cost(), expected);
    }
}

TEST(KMedoidsCostTest, CostStrictlyDecreasesAcrossAcceptedIterations) {
    const DenseStorage storage = MakeEuclideanStorage(30, 777, 1.0);

    KMedoidsOptions options;
    options.n_clusters = 3;
    options.init = KMedoidsInit::Explicit;
    options.initial_medoids = {0, 1, 2};

    const KMedoidsResult converged_result = k_medoids_cluster(storage, options);
    ASSERT_GT(converged_result.NumIterations(), 0u);

    // Seeded with the cost of the untouched seeds rather than with infinity, so
    // the cap == 1 comparison asserts that the first accepted iteration
    // decreased the cost instead of passing vacuously.
    double previous = NaiveTotalCost(storage.Data(), storage.NumSamples(),
                                     {0, 1, 2});
    for (size_t cap = 1; cap <= converged_result.NumIterations(); ++cap) {
        KMedoidsOptions capped = options;
        capped.max_iterations = cap;
        const double cost = k_medoids_cluster(storage, capped).Cost();
        EXPECT_LT(cost, previous);
        previous = cost;
    }
}

// Also the initialization-convergence row: both initializers must converge
// inside the default cap on ordinary data. They may land on different local
// optima -- PAM guarantees no more than that -- so the assertion is that each
// result is a verified optimum, not that the two agree.
TEST(KMedoidsConvergenceTest, EveryConvergedResultIsAVerifiedLocalOptimum) {
    for (const KMedoidsInit init :
         {KMedoidsInit::Build, KMedoidsInit::FarthestFirst}) {
        for (const uint32_t seed : {uint32_t{31}, uint32_t{32}, uint32_t{33}}) {
            for (const size_t k : {size_t{2}, size_t{4}}) {
                const DenseStorage metric = MakeEuclideanStorage(28, seed, 1.0);
                const DenseStorage non_metric =
                    MakeNonMetricStorage(28, seed, 1.0);

                KMedoidsOptions options;
                options.n_clusters = k;
                options.init = init;

                const KMedoidsResult from_metric =
                    k_medoids_cluster(metric, options);
                ASSERT_TRUE(from_metric.Converged());
                ExpectVerifiedLocalOptimum(metric, from_metric);

                const KMedoidsResult from_non_metric =
                    k_medoids_cluster(non_metric, options);
                ASSERT_TRUE(from_non_metric.Converged());
                ExpectVerifiedLocalOptimum(non_metric, from_non_metric);
            }
        }
    }
}

// Parity and verified local optimality at the two extremes of absolute scale,
// fifteen orders of magnitude apart. That is the whole claim: these matrices
// are randomly generated, and on them the predicted delta and the recomputed
// total agree about which swaps improve, so no improvement decision here turns
// on the disagreement between the two. An implementation that terminated on
// predicted deltas alone still passes this row -- measured: with the swap
// loop's verification result replaced by a default-constructed candidate, so
// that the loop converges on the fast path alone, the only failure in the
// focused suite is
// KMedoidsDynamicRangeTest.HoldsParityWhenTheVerificationPassAppliesTheSwap,
// which is handcrafted for exactly that purpose.
TEST(KMedoidsDynamicRangeTest, HoldsParityAndOptimalityAtExtremeScales) {
    for (const double scale : {1e-9, 1e6}) {
        const DenseStorage metric = MakeEuclideanStorage(24, 5150, scale);
        const DenseStorage non_metric = MakeNonMetricStorage(24, 5151, scale);

        ExpectParityWithNaivePam(metric, 3);
        ExpectParityWithNaivePam(non_metric, 3);

        KMedoidsOptions options;
        options.n_clusters = 3;
        const KMedoidsResult result = k_medoids_cluster(metric, options);
        ASSERT_TRUE(result.Converged());
        ExpectVerifiedLocalOptimum(metric, result);
    }
}

TEST(KMedoidsDynamicRangeTest, HoldsParityOnAMatrixMixingBothScales) {
    DenseStorage storage(24);
    const DenseStorage small = MakeEuclideanStorage(24, 6160, 1e-9);
    const DenseStorage large = MakeEuclideanStorage(24, 6161, 1e6);
    for (size_t i = 0; i < 24; ++i) {
        for (size_t j = i + 1; j < 24; ++j) {
            const DenseStorage& source = (i + j) % 2 == 0 ? small : large;
            storage.Set(i, j, detail::dense_distance(source.Data(), 24, i, j));
        }
    }

    ExpectParityWithNaivePam(storage, 3);

    KMedoidsOptions options;
    options.n_clusters = 3;
    const KMedoidsResult result = k_medoids_cluster(storage, options);
    ASSERT_TRUE(result.Converged());
    ExpectVerifiedLocalOptimum(storage, result);
}

// The randomly generated dynamic-range matrices never actually separate the
// predicted delta from the recomputed total: on every one of them the two agree
// on which swaps improve, so the swap loop's two correction branches stay
// unvisited. These last two rows are built to separate them, one in each
// direction. They are the only rows in the suite that reach either branch.

// A predicted improvement that is not one. What the undo protects is the
// strict-improvement invariant the termination argument rests on: an accepted
// configuration must have a strictly smaller recomputed cost than its
// predecessor, which is what makes a configuration impossible to revisit.
// Deleting the trial rebuild -- taking the prediction at its word -- bills a
// zero-improvement move as progress and reports {0, 5} at 2^53 + 10 after one
// iteration where the loop should report {0, 6} at 2^53 + 10 after none. The
// two partitions cost exactly the same; it is the accounting that breaks.
TEST(KMedoidsDynamicRangeTest, HoldsParityWhenAPredictedImprovementIsUndone) {
    const DenseStorage storage = MakeSpeculativeUndoStorage();
    const std::vector<size_t> seeds = {0, 6};

    const KernelState state = MakeKernelState(storage, seeds);
    const detail::SwapCandidate predicted = detail::best_predicted_swap(
        storage.Data(), storage.NumSamples(), state.medoids, state.assignments,
        state.slot_of, KERNEL_THREADS, KERNEL_CHUNK);
    ASSERT_TRUE(predicted.valid);
    // The precondition for the branch under test: the fast path is entered
    // because the prediction is negative, and the rebuild then declines it.
    //
    // These two are also load-bearing on their own, and not redundant with the
    // end-to-end assertions below. Together with the counterpart score
    // assertion in the row after this one, they are all that pins the ascending
    // accumulation order of `shared` and `correction` in best_predicted_swap:
    // reversing that inner loop changes the answer here to score 0.0 entering
    // item 3 while leaving every end-to-end result in this file unchanged.
    EXPECT_EQ(predicted.score, -2.0);
    EXPECT_EQ(predicted.entering_item, 5u);
    // Stated through the cached shortcut rather than through the rebuild the
    // loop itself performs. The two agree on this matrix, and
    // RecomputedTotalIsBitIdenticalToARebuiltAssignmentCost pins that agreement
    // on this exact fixture rather than leaving it to the general claim.
    EXPECT_EQ(detail::recomputed_total(storage.Data(), storage.NumSamples(),
                                       state.assignments,
                                       predicted.leaving_slot,
                                       predicted.entering_item),
              state.cost);

    ExpectParityFromSeeds(storage, seeds);

    KMedoidsOptions options;
    options.n_clusters = 2;
    options.init = KMedoidsInit::Explicit;
    options.initial_medoids = seeds;
    const KMedoidsResult result = k_medoids_cluster(storage, options);

    EXPECT_EQ(result.Medoids(), std::vector<size_t>({0, 6}));
    EXPECT_EQ(result.Cost(), TWO_POW_53 + 10.0);
    // The undone swap must not be billed as progress.
    EXPECT_EQ(result.NumIterations(), 0u);
    EXPECT_TRUE(result.Converged());
    ExpectVerifiedLocalOptimum(storage, result);
}

// The reverse: no predicted delta is negative, so the fast path declines, and
// only the terminal verification pass sees the improvement. Deleting that pass
// reports medoids {0, 6} at cost 2^53 + 10 and calls it converged.
TEST(KMedoidsDynamicRangeTest, HoldsParityWhenTheVerificationPassAppliesTheSwap) {
    const DenseStorage storage = MakeMissedImprovementStorage();
    const std::vector<size_t> seeds = {0, 6};

    const KernelState state = MakeKernelState(storage, seeds);
    const detail::SwapCandidate predicted = detail::best_predicted_swap(
        storage.Data(), storage.NumSamples(), state.medoids, state.assignments,
        state.slot_of, KERNEL_THREADS, KERNEL_CHUNK);
    ASSERT_TRUE(predicted.valid);
    // The counterpart of the score assertion in the row above, and the other
    // half of what pins best_predicted_swap's accumulation order: reversing that
    // inner loop turns this 0.0 into -2.0, and outside these two rows nothing in
    // the suite notices.
    EXPECT_EQ(predicted.score, 0.0);
    EXPECT_EQ(detail::recomputed_total(storage.Data(), storage.NumSamples(),
                                       state.assignments, 1, 5),
              TWO_POW_53 + 8.0);

    ExpectParityFromSeeds(storage, seeds);

    KMedoidsOptions options;
    options.n_clusters = 2;
    options.init = KMedoidsInit::Explicit;
    options.initial_medoids = seeds;
    const KMedoidsResult result = k_medoids_cluster(storage, options);

    EXPECT_EQ(result.Medoids(), std::vector<size_t>({0, 5}));
    EXPECT_EQ(result.Cost(), TWO_POW_53 + 8.0);
    EXPECT_EQ(result.NumIterations(), 1u);
    EXPECT_TRUE(result.Converged());
    ExpectVerifiedLocalOptimum(storage, result);
}

// Complements the exactly-k assertions already pinned in
// KMedoidsDegenerateTest.AllZeroMatrixStillProducesDistinctMedoids: on a
// matrix where every candidate swap has a delta of exactly zero, no swap
// improves, so the run must stop at its seeds and still claim convergence.
TEST(KMedoidsDegenerateTest, AllZeroMatrixConvergesAtTheInitializationSeeds) {
    for (const KMedoidsInit init :
         {KMedoidsInit::Build, KMedoidsInit::FarthestFirst}) {
        for (const size_t k : {size_t{2}, size_t{3}, size_t{5}}) {
            KMedoidsOptions options;
            options.n_clusters = k;
            options.init = init;

            const KMedoidsResult result =
                k_medoids_cluster(MakeAllZeroStorage(8), options);

            EXPECT_EQ(result.NumIterations(), 0u);
            EXPECT_TRUE(result.Converged());
        }
    }
}

TEST(KMedoidsCapTest, ReachingTheCapReportsAValidUnconvergedPartition) {
    const DenseStorage storage = MakeEuclideanStorage(30, 777, 1.0);

    KMedoidsOptions options;
    options.n_clusters = 3;
    options.init = KMedoidsInit::Explicit;
    options.initial_medoids = {0, 1, 2};
    options.max_iterations = 1;

    const KMedoidsResult result = k_medoids_cluster(storage, options);

    EXPECT_FALSE(result.Converged());
    EXPECT_EQ(result.NumIterations(), 1u);
    EXPECT_EQ(result.NumClusters(), 3u);
    EXPECT_EQ(result.Medoids().size(), 3u);

    double expected = 0.0;
    for (size_t j = 0; j < storage.NumSamples(); ++j) {
        const size_t medoid =
            result.Medoids()[static_cast<size_t>(result.Labels()[j])];
        expected += detail::dense_distance(storage.Data(),
                                           storage.NumSamples(), j, medoid);
    }
    EXPECT_DOUBLE_EQ(result.Cost(), expected);

    // The cap skips the verification pass, so this input is deliberately one
    // that is not yet a local optimum: an improving swap still exists.
    const std::vector<size_t>& medoids = result.Medoids();
    bool improving_swap_exists = false;
    for (size_t h = 0; h < storage.NumSamples(); ++h) {
        for (size_t slot = 0; slot < medoids.size(); ++slot) {
            std::vector<size_t> trial = medoids;
            trial[slot] = h;
            if (NaiveTotalCost(storage.Data(), storage.NumSamples(), trial) <
                result.Cost()) {
                improving_swap_exists = true;
            }
        }
    }
    EXPECT_TRUE(improving_swap_exists);
}

TEST(KMedoidsDegenerateTest, SingleClusterReturnsTheGlobalMedoid) {
    const DenseStorage storage = MakeEuclideanStorage(25, 3131, 1.0);
    const size_t n = storage.NumSamples();
    const double* data = storage.Data();

    KMedoidsOptions options;
    options.n_clusters = 1;
    const KMedoidsResult result = k_medoids_cluster(storage, options);

    size_t expected = 0;
    double best = std::numeric_limits<double>::infinity();
    for (size_t candidate = 0; candidate < n; ++candidate) {
        double sum = 0.0;
        for (size_t j = 0; j < n; ++j) {
            sum += detail::dense_distance(data, n, candidate, j);
        }
        if (sum < best) {
            best = sum;
            expected = candidate;
        }
    }

    EXPECT_EQ(result.Medoids(), std::vector<size_t>({expected}));
    EXPECT_EQ(result.NumIterations(), 0u);
    EXPECT_TRUE(result.Converged());
    EXPECT_EQ(result.NumClusters(), 1u);
}
