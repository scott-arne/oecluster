#include <gtest/gtest.h>

#include <oechem.h>

#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <string>
#include <utility>
#include <vector>

#include "oecluster/Error.h"
#include "oecluster/GateFacts.h"
#include "oecluster/PairwiseComparison.h"
#include "oecluster/StorageBackend.h"
#include "oecluster/clustering/SARCoherence.h"
#include "oecluster/comparisons/FingerprintComparison.h"
#include "diversity_test_support.h"
#include "report_test_support.h"

using namespace OECluster;
using diversity_test::CountingComparison;
using diversity_test::Hashed;
using diversity_test::IsolationComparison;
using diversity_test::MakeStorage;
using diversity_test::Quantized;
using diversity_test::TableComparison;
using report_test::SameDouble;

namespace {

constexpr size_t UNBOUNDED = std::numeric_limits<size_t>::max();
constexpr double NOT_A_NUMBER = std::numeric_limits<double>::quiet_NaN();

// 150 rows is enough for activity_landscape's 64-row floor to leave several
// units; 20000 is the one chunk_size here that lifts a unit above the floor.
constexpr size_t N = 150;
const size_t CHUNKS[] = {1, 7, 4096, 20000, UNBOUNDED};
const size_t THREADS[] = {1, 2, 4, 0};

size_t CondensedIndex(size_t n, size_t i, size_t j) {
    return n * i - i * (i + 1) / 2 + j - i - 1;
}

// Missing at 5, 18, 31, ...: those rows and their pairs are never read.
std::vector<double> Activity(size_t n) {
    std::vector<double> activity(n);
    for (size_t i = 0; i < n; ++i) {
        activity[i] = i % 13 == 5 ? NOT_A_NUMBER : static_cast<double>((i * 7) % 11) * 0.5;
    }
    return activity;
}

// Three classes, missing at 4, 13, 22, ...
std::vector<std::string> Classes(size_t n) {
    const char* names[] = {"a", "b", "c"};
    std::vector<std::string> classes(n);
    for (size_t i = 0; i < n; ++i) {
        classes[i] = i % 9 == 4 ? "" : names[(i * 7 + i / 5) % 3];
    }
    return classes;
}

// One class, so modelability takes its serial validation scan.
std::vector<std::string> OneClass(size_t n) {
    std::vector<std::string> classes(n);
    for (size_t i = 0; i < n; ++i) {
        classes[i] = i % 9 == 4 ? "" : "a";
    }
    return classes;
}

// Hashed has no zero distances; Quantized on three levels has zeros and
// heavy ties, which drive num_zero_distance_pairs and the tie rules.
std::vector<std::vector<double>> Fixtures(size_t n) {
    return {Hashed(n), Quantized(n, 2, 3)};
}

::testing::AssertionResult SameLandscape(const ActivityLandscape& a,
                                         const ActivityLandscape& b) {
    if (a.num_samples != b.num_samples || a.num_scored != b.num_scored ||
        a.num_pairs_scored != b.num_pairs_scored || a.num_cliffs != b.num_cliffs ||
        a.num_zero_distance_pairs != b.num_zero_distance_pairs) {
        return ::testing::AssertionFailure() << "a count differs";
    }
    const std::pair<const char*, std::pair<double, double>> reals[] = {
        {"cliff_density", {a.cliff_density, b.cliff_density}},
        {"max_sali", {a.max_sali, b.max_sali}},
        {"mean_sali", {a.mean_sali, b.mean_sali}},
        {"rmodi", {a.rmodi, b.rmodi}},
        {"activity_stddev", {a.activity_stddev, b.activity_stddev}},
    };
    for (const auto& real : reals) {
        if (!SameDouble(real.second.first, real.second.second)) {
            return ::testing::AssertionFailure()
                   << real.first << ": " << report_test::Hex(real.second.first) << " vs "
                   << report_test::Hex(real.second.second);
        }
    }
    return ::testing::AssertionSuccess();
}

::testing::AssertionResult SameModelability(const Modelability& a, const Modelability& b) {
    if (a.num_samples != b.num_samples || a.num_scored != b.num_scored ||
        a.num_classes != b.num_classes || a.classes.size() != b.classes.size()) {
        return ::testing::AssertionFailure() << "a count differs";
    }
    if (!SameDouble(a.modi, b.modi)) {
        return ::testing::AssertionFailure()
               << "modi: " << report_test::Hex(a.modi) << " vs " << report_test::Hex(b.modi);
    }
    for (size_t k = 0; k < a.classes.size(); ++k) {
        if (a.classes[k].label != b.classes[k].label ||
            a.classes[k].num_members != b.classes[k].num_members ||
            !SameDouble(a.classes[k].fraction_same_class, b.classes[k].fraction_same_class)) {
            return ::testing::AssertionFailure() << "class row " << k << " differs";
        }
    }
    return ::testing::AssertionSuccess();
}

std::string MessageOf(const std::function<void()>& call) {
    try {
        call();
    } catch (const std::exception& error) {
        return error.what();
    }
    return "<no exception>";
}

}  // namespace

// TableComparison throws on a reversed pair or a self-pair, so every parity
// test below also proves the engines read only Compare(min, max).
TEST(SARComparisonTest, LandscapeMatchesTheMatrixForEveryThreadAndChunk) {
    const std::vector<double> activity = Activity(N);
    for (const std::vector<double>& condensed : Fixtures(N)) {
        const DenseStorage storage = MakeStorage(N, condensed);
        ActivityLandscapeOptions options;
        const ActivityLandscape expected = activity_landscape(storage, activity, options);
        ASSERT_GT(expected.num_cliffs, 0u);
        for (const size_t threads : THREADS) {
            for (const size_t chunk : CHUNKS) {
                TableComparison comparison(N, condensed);
                options.num_threads = threads;
                options.chunk_size = chunk;
                EXPECT_TRUE(SameLandscape(expected,
                                          activity_landscape(comparison, activity, options)))
                    << "threads " << threads << " chunk " << chunk;
            }
        }
    }
}

TEST(SARComparisonTest, ModelabilityMatchesTheMatrixForEveryThreadAndChunk) {
    for (const std::vector<std::string>& classes : {Classes(N), OneClass(N)}) {
        for (const std::vector<double>& condensed : Fixtures(N)) {
            const DenseStorage storage = MakeStorage(N, condensed);
            ModelabilityOptions options;
            const Modelability expected = modelability(storage, classes, options);
            for (const size_t threads : THREADS) {
                for (const size_t chunk : CHUNKS) {
                    TableComparison comparison(N, condensed);
                    options.num_threads = threads;
                    options.chunk_size = chunk;
                    EXPECT_TRUE(SameModelability(expected,
                                                 modelability(comparison, classes, options)))
                        << "classes " << expected.num_classes << " threads " << threads
                        << " chunk " << chunk;
                }
            }
        }
    }
}

TEST(SARComparisonTest, MatchesAMatrixFilledFromAFingerprintComparison) {
    const std::vector<const char*> smiles = {
        "c1ccccc1", "c1ccc(O)cc1", "c1ccc(N)cc1", "c1ccc(Cl)cc1",
        "CCCCCCCC", "CCCCCCO", "CCCCCCN", "CC(C)CCO",
        "C1CCCCC1", "C1CCCCC1O", "c1ccncc1", "CC(=O)Oc1ccccc1C(=O)O",
    };
    std::vector<OEChem::OEGraphMol> graph_mols(smiles.size());
    std::vector<OEChem::OEMolBase*> mols;
    for (size_t i = 0; i < smiles.size(); ++i) {
        ASSERT_TRUE(OEChem::OESmilesToMol(graph_mols[i], smiles[i]));
        mols.push_back(&static_cast<OEChem::OEMolBase&>(graph_mols[i]));
    }
    FingerprintComparison comparison(mols);
    const size_t n = mols.size();
    DenseStorage storage(n);
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i + 1; j < n; ++j) {
            storage.Set(i, j, comparison.Compare(i, j));
        }
    }
    const std::vector<double> activity = {5.0, 6.2, 5.1, 7.4, 4.0, 4.3,
                                          6.9, 4.1, 5.5, 5.6, 8.0, NOT_A_NUMBER};
    const std::vector<std::string> classes = {"a", "a", "b", "a", "c", "c",
                                              "b", "c", "a", "a", "b", ""};
    for (const size_t threads : {size_t{1}, size_t{4}}) {
        ActivityLandscapeOptions landscape_options;
        landscape_options.num_threads = threads;
        EXPECT_TRUE(SameLandscape(activity_landscape(storage, activity, landscape_options),
                                  activity_landscape(comparison, activity, landscape_options)))
            << "threads " << threads;
        ModelabilityOptions model_options;
        model_options.num_threads = threads;
        model_options.chunk_size = 1;
        EXPECT_TRUE(SameModelability(modelability(storage, classes, model_options),
                                     modelability(comparison, classes, model_options)))
            << "threads " << threads;
    }
}

TEST(SARComparisonTest, LandscapeReportsTheMatrixMessageForASingleBadPair) {
    const std::vector<double> activity = Activity(N);
    for (const double bad : {NOT_A_NUMBER, std::numeric_limits<double>::infinity(), -0.25}) {
        std::vector<double> condensed = Hashed(N);
        condensed[CondensedIndex(N, 3, 16)] = bad;
        const DenseStorage storage = MakeStorage(N, condensed);
        ActivityLandscapeOptions options;
        const std::string expected =
            MessageOf([&] { activity_landscape(storage, activity, options); });
        EXPECT_EQ(expected,
                  "activity_landscape: the distance between samples 3 and 16 must be "
                  "finite and non-negative");
        for (const size_t threads : {size_t{1}, size_t{4}}) {
            for (const size_t chunk : {size_t{1}, size_t{20000}}) {
                TableComparison comparison(N, condensed);
                options.num_threads = threads;
                options.chunk_size = chunk;
                EXPECT_EQ(MessageOf([&] { activity_landscape(comparison, activity, options); }),
                          expected)
                    << "bad " << bad << " threads " << threads << " chunk " << chunk;
            }
        }
    }
}

TEST(SARComparisonTest, ModelabilityReportsTheMatrixMessageForASingleBadPair) {
    for (const std::vector<std::string>& classes : {Classes(N), OneClass(N)}) {
        for (const double bad : {NOT_A_NUMBER, std::numeric_limits<double>::infinity(), -0.25}) {
            std::vector<double> condensed = Hashed(N);
            condensed[CondensedIndex(N, 3, 16)] = bad;
            const DenseStorage storage = MakeStorage(N, condensed);
            ModelabilityOptions options;
            const std::string expected =
                MessageOf([&] { modelability(storage, classes, options); });
            EXPECT_EQ(expected,
                      "modelability: the distance between samples 3 and 16 must be finite "
                      "and non-negative");
            for (const size_t threads : {size_t{1}, size_t{4}}) {
                for (const size_t chunk : {size_t{1}, size_t{4096}}) {
                    TableComparison comparison(N, condensed);
                    options.num_threads = threads;
                    options.chunk_size = chunk;
                    EXPECT_EQ(MessageOf([&] { modelability(comparison, classes, options); }),
                              expected)
                        << "bad " << bad << " threads " << threads << " chunk " << chunk;
                }
            }
        }
    }
}

TEST(SARComparisonTest, ANonFiniteValueOnAnUnreadPairIsNeverSeen) {
    // 5 and 18 have no activity, so their pair is never read by the landscape.
    // 4 and 13 have no class, so their pair is never read by modelability. The
    // two guarantees are independent -- the landscape still reads (4, 13) (its
    // endpoints have activity) and modelability still reads (5, 18) (its
    // endpoints have a class) -- so each poisoned pair is kept on its own
    // comparison rather than shared, or the other function would see it.
    const DenseStorage clean = MakeStorage(N, Hashed(N));

    std::vector<double> landscape_poisoned = Hashed(N);
    landscape_poisoned[CondensedIndex(N, 5, 18)] = NOT_A_NUMBER;
    TableComparison landscape_comparison(N, landscape_poisoned);
    EXPECT_TRUE(SameLandscape(activity_landscape(clean, Activity(N)),
                              activity_landscape(landscape_comparison, Activity(N))));

    std::vector<double> model_poisoned = Hashed(N);
    model_poisoned[CondensedIndex(N, 4, 13)] = NOT_A_NUMBER;
    TableComparison model_comparison(N, model_poisoned);
    EXPECT_TRUE(SameModelability(modelability(clean, Classes(N)),
                                 modelability(model_comparison, Classes(N))));
}

TEST(SARComparisonTest, ClonesAreIsolatedAndBoundedByTheItems) {
    const std::vector<double> condensed = Hashed(N);
    const DenseStorage storage = MakeStorage(N, condensed);
    {
        IsolationComparison isolation(N, condensed);
        ActivityLandscapeOptions options;
        options.num_threads = 4;
        options.chunk_size = 1;
        EXPECT_TRUE(SameLandscape(activity_landscape(storage, Activity(N), options),
                                  activity_landscape(isolation, Activity(N), options)));
        EXPECT_EQ(isolation.Violations(), 0u);
        EXPECT_TRUE(isolation.OverlapObserved())
            << "Overlap not observed; test may be flaky on this machine";
    }
    {
        IsolationComparison isolation(N, condensed);
        ModelabilityOptions options;
        options.num_threads = 4;
        options.chunk_size = 1;
        EXPECT_TRUE(SameModelability(modelability(storage, Classes(N), options),
                                     modelability(isolation, Classes(N), options)));
        EXPECT_EQ(isolation.Violations(), 0u);
        EXPECT_TRUE(isolation.OverlapObserved())
            << "Overlap not observed; test may be flaky on this machine";
    }

    // More workers requested than items, the largest count, and the automatic
    // count: whatever is asked for, no more clones than items are made.
    const size_t n = 10;
    const std::vector<double> small = Hashed(n);
    for (const size_t threads : {size_t{64}, UNBOUNDED, size_t{0}}) {
        TableComparison landscape_table(n, small);
        ActivityLandscapeOptions landscape_options;
        landscape_options.num_threads = threads;
        landscape_options.chunk_size = 1;
        activity_landscape(landscape_table, Activity(n), landscape_options);
        EXPECT_GE(landscape_table.NumClones(), 1u);
        EXPECT_LE(landscape_table.NumClones(), n) << "threads " << threads;

        TableComparison model_table(n, small);
        ModelabilityOptions model_options;
        model_options.num_threads = threads;
        model_options.chunk_size = 1;
        modelability(model_table, Classes(n), model_options);
        EXPECT_GE(model_table.NumClones(), 1u);
        EXPECT_LE(model_table.NumClones(), n) << "threads " << threads;
    }
}

TEST(SARComparisonTest, CountsTheDocumentedComparisons) {
    // 24 samples: Activity leaves 22 scored (5 and 18 missing), Classes and
    // OneClass leave 21 (4, 13 and 22 missing). Only scored pairs are read.
    const size_t n = 24;
    for (const size_t threads : {size_t{1}, size_t{4}}) {
        CountingComparison landscape_counter(n, GateFacts());
        ActivityLandscapeOptions landscape_options;
        landscape_options.num_threads = threads;
        activity_landscape(landscape_counter, Activity(n), landscape_options);
        EXPECT_EQ(landscape_counter.Count(), 22u * 21u / 2u) << "threads " << threads;

        // Two or more classes: every row visits every other scored sample,
        // so each pair is compared twice.
        CountingComparison classes_counter(n, GateFacts());
        ModelabilityOptions model_options;
        model_options.num_threads = threads;
        modelability(classes_counter, Classes(n), model_options);
        EXPECT_EQ(classes_counter.Count(), 21u * 20u) << "threads " << threads;

        // One class: nothing to score, so each pair is only validated once.
        CountingComparison one_class_counter(n, GateFacts());
        modelability(one_class_counter, OneClass(n), model_options);
        EXPECT_EQ(one_class_counter.Count(), 21u * 20u / 2u) << "threads " << threads;
    }
}

TEST(SARComparisonTest, RefusesFactsBeforeAnyComparison) {
    GateFacts similarity;
    similarity.is_distance = Capability::No;
    GateFacts nonzero_self;
    nonzero_self.zero_self = Capability::No;
    GateFacts nan_present;
    nan_present.data_integrity = DataIntegrity::NaNPresent;
    GateFacts subset;
    subset.data_integrity = DataIntegrity::SubsetScored;
    for (const GateFacts& facts : {similarity, nonzero_self, nan_present, subset}) {
        CountingComparison comparison(24, facts);
        // A mismatched length is checked after the facts, so the facts win.
        EXPECT_THROW(activity_landscape(comparison, Activity(20)), ComparisonError);
        EXPECT_THROW(modelability(comparison, Classes(20)), ComparisonError);
        EXPECT_EQ(comparison.Count(), 0u);
    }
}

TEST(SARComparisonTest, RefusesAZeroChunkFirst) {
    GateFacts similarity;
    similarity.is_distance = Capability::No;
    CountingComparison comparison(24, similarity);
    ActivityLandscapeOptions landscape_options;
    landscape_options.chunk_size = 0;
    EXPECT_EQ(MessageOf([&] {
                  activity_landscape(comparison, Activity(24), landscape_options);
              }),
              "activity_landscape chunk_size must be at least one");
    ModelabilityOptions model_options;
    model_options.chunk_size = 0;
    EXPECT_EQ(MessageOf([&] { modelability(comparison, Classes(24), model_options); }),
              "modelability chunk_size must be at least one");
    EXPECT_EQ(comparison.Count(), 0u);
}

TEST(SARComparisonTest, NamesTheComparisonInTheLengthRefusals) {
    TableComparison comparison(20, Hashed(20));
    EXPECT_EQ(MessageOf([&] { activity_landscape(comparison, Activity(24)); }),
              "activity_landscape: activity has 24 entries but the comparison has 20 samples");
    EXPECT_EQ(MessageOf([&] { modelability(comparison, Classes(24)); }),
              "modelability: activity_classes has 24 entries but the comparison has 20 "
              "samples");
    const DenseStorage storage = MakeStorage(20, Hashed(20));
    EXPECT_EQ(MessageOf([&] { modelability(storage, Classes(24)); }),
              "modelability: activity_classes has 24 entries but the storage has 20 samples");
}

TEST(SARComparisonTest, RefusesBadOptionsBeforeAnyComparison) {
    CountingComparison comparison(24, GateFacts());
    ActivityLandscapeOptions options;
    options.rmodi_delta = -1.0;
    EXPECT_EQ(MessageOf([&] { activity_landscape(comparison, Activity(24), options); }),
              "activity_landscape: rmodi_delta must be non-negative");
    EXPECT_EQ(comparison.Count(), 0u);
}
