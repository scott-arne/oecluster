/**
 * @file test_comparison_repeatability.cpp
 * @brief Every built-in comparison family returns the same bits for a pair,
 * whatever its clone has scored before.
 *
 * The comparison-built threshold graph scores every pair twice, once to size
 * each row and once to fill it, on whichever clone the scheduler hands the
 * chunk. Its exactness rests on Compare being repeatable across differing
 * clone histories, which PairwiseComparison does not promise and ROCS clones
 * hold mutable toolkit state for. These tests are the evidence for that
 * precondition, family by family.
 */

#include <gtest/gtest-spi.h>
#include <gtest/gtest.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <memory>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <oebio.h>
#include <oechem.h>
#include <oeomega2.h>

#include "oecluster/StorageBackend.h"
#include "oecluster/comparisons/DescriptorComparison.h"
#include "oecluster/comparisons/FingerprintComparison.h"
#include "oecluster/comparisons/MCSComparison.h"
#include "oecluster/comparisons/RMSDComparison.h"
#include "oecluster/comparisons/ROCSComparison.h"
#include "oecluster/comparisons/SuperposeComparison.h"

#include "streaming_test_support.h"

using namespace OECluster;
using namespace streaming_test;

namespace {

struct Pair {
    size_t i;
    size_t j;
};

uint64_t Bits(double value) {
    uint64_t bits = 0;
    std::memcpy(&bits, &value, sizeof bits);
    return bits;
}

std::vector<Pair> AllPairs(size_t n) {
    std::vector<Pair> pairs;
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i + 1; j < n; ++j) {
            pairs.push_back({i, j});
        }
    }
    return pairs;
}

std::vector<Pair> Shuffled(std::vector<Pair> pairs, unsigned seed) {
    std::mt19937 generator(seed);
    std::shuffle(pairs.begin(), pairs.end(), generator);
    return pairs;
}

// Every pair scored on `clone` in `order`; the values come back in condensed
// order whatever the order of scoring.
std::vector<double> Score(PairwiseComparison& clone, size_t n,
                          const std::vector<Pair>& order) {
    std::vector<double> values(n * (n - 1) / 2);
    for (const Pair& pair : order) {
        values[n * pair.i - pair.i * (pair.i + 1) / 2 + pair.j - pair.i - 1] =
            clone.Compare(pair.i, pair.j);
    }
    return values;
}

void ExpectSameBits(const std::vector<double>& actual,
                    const std::vector<double>& expected) {
    ASSERT_EQ(actual.size(), expected.size());
    for (size_t k = 0; k < actual.size(); ++k) {
        EXPECT_EQ(Bits(actual[k]), Bits(expected[k]))
            << "pair " << k << ": " << actual[k] << " vs " << expected[k];
    }
}

// The three checks the design names. (a) One clone reads every pair a
// second time after a shuffled history of the others. (b) Fresh clones,
// each first given its own shuffled history, read every pair again; so does
// a clone with no history at all. (c) Whole graphs built under several
// thread and chunk settings, where the scheduler decides which clone scores
// which chunk on each pass, equal the graph of the baseline values.
void ExpectRepeatable(PairwiseComparison& prototype,
                      const std::vector<size_t>& thread_counts,
                      const std::vector<size_t>& chunk_sizes) {
    const size_t n = prototype.Size();
    ASSERT_GE(n, 2u);
    const std::vector<Pair> pairs = AllPairs(n);

    std::unique_ptr<PairwiseComparison> first = prototype.Clone();
    const std::vector<double> baseline = Score(*first, n, Shuffled(pairs, 1));
    ExpectSameBits(Score(*first, n, Shuffled(pairs, 2)), baseline);

    for (const unsigned seed : {3u, 4u}) {
        std::unique_ptr<PairwiseComparison> clone = prototype.Clone();
        Score(*clone, n, Shuffled(pairs, seed + 100));
        ExpectSameBits(Score(*clone, n, Shuffled(pairs, seed)), baseline);
    }
    std::unique_ptr<PairwiseComparison> fresh = prototype.Clone();
    ExpectSameBits(Score(*fresh, n, pairs), baseline);

    DenseStorage storage(n);
    for (const Pair& pair : pairs) {
        storage.Set(pair.i, pair.j,
                    baseline[n * pair.i - pair.i * (pair.i + 1) / 2 + pair.j -
                             pair.i - 1]);
    }
    // The median admits about half the pairs, so the graph is neither empty
    // nor complete; clamped because the builder refuses a negative threshold.
    std::vector<double> sorted = baseline;
    std::sort(sorted.begin(), sorted.end());
    ThresholdGraphOptions options;
    options.threshold = std::max(0.0, sorted[sorted.size() / 2]);
    options.num_threads = 1;
    const auto expected = Rows(BuildThresholdNeighborGraph(storage, options));
    for (const size_t threads : thread_counts) {
        for (const size_t chunk : chunk_sizes) {
            SCOPED_TRACE(::testing::Message()
                         << "threads " << threads << ", chunk " << chunk);
            options.num_threads = threads;
            options.chunk_size = chunk;
            // A pass mismatch is the builder's own net for this defect;
            // reported as a failure so every family still runs.
            try {
                EXPECT_EQ(Rows(BuildThresholdNeighborGraph(prototype, options)),
                          expected);
            } catch (const std::logic_error& error) {
                ADD_FAILURE() << error.what();
            }
        }
    }
}

const std::vector<size_t> THREADS{1, 2, 4};
const std::vector<size_t> CHUNKS{1, 3, 4096};

std::vector<OEChem::OEMolBase*> Pointers(std::vector<OEChem::OEGraphMol>& mols) {
    std::vector<OEChem::OEMolBase*> pointers;
    for (auto& mol : mols) {
        pointers.push_back(&static_cast<OEChem::OEMolBase&>(mol));
    }
    return pointers;
}

std::vector<OEChem::OEGraphMol> GraphMols(const std::vector<const char*>& smiles) {
    std::vector<OEChem::OEGraphMol> mols;
    for (const char* smi : smiles) {
        mols.emplace_back();
        EXPECT_TRUE(OEChem::OESmilesToMol(mols.back(), smi)) << smi;
    }
    return mols;
}

// One Omega conformer per molecule: OEShape finds nothing to overlay on
// planar input (test_rocs_comparison.cpp explains), so ROCS needs 3D.
std::vector<std::shared_ptr<OEChem::OEMol>> Conformers(
    const std::vector<const char*>& smiles) {
    OEConfGen::OEOmega omega;
    omega.SetMaxConfs(1);
    omega.SetStrictStereo(false);
    std::vector<std::shared_ptr<OEChem::OEMol>> mols;
    for (const char* smi : smiles) {
        auto mol = std::make_shared<OEChem::OEMol>();
        EXPECT_TRUE(OEChem::OESmilesToMol(*mol, smi)) << smi;
        EXPECT_TRUE(omega(*mol)) << smi;
        mols.push_back(mol);
    }
    return mols;
}

// Several poses of one molecule, each its own single-conformer item, so RMSD
// has a shared topology and genuinely different coordinates.
std::vector<std::shared_ptr<OEChem::OEMol>> Poses(const char* smiles,
                                                  unsigned int count) {
    OEConfGen::OEOmega omega;
    omega.SetMaxConfs(count);
    omega.SetStrictStereo(false);
    OEChem::OEMol ensemble;
    EXPECT_TRUE(OEChem::OESmilesToMol(ensemble, smiles)) << smiles;
    EXPECT_TRUE(omega(ensemble)) << smiles;
    std::vector<std::shared_ptr<OEChem::OEMol>> poses;
    for (OESystem::OEIter<OEChem::OEConfBase> conf = ensemble.GetConfs(); conf;
         ++conf) {
        poses.push_back(std::make_shared<OEChem::OEMol>(*conf));
    }
    return poses;
}

std::vector<std::shared_ptr<OEChem::OEMol>> Topological(
    const std::vector<const char*>& smiles) {
    std::vector<std::shared_ptr<OEChem::OEMol>> mols;
    for (const char* smi : smiles) {
        auto mol = std::make_shared<OEChem::OEMol>();
        EXPECT_TRUE(OEChem::OESmilesToMol(*mol, smi)) << smi;
        mols.push_back(mol);
    }
    return mols;
}

const std::vector<const char*> FINGERPRINT_SMILES{
    "CCO",      "CCCO",       "CCCCO",   "c1ccccc1", "Cc1ccccc1", "CCc1ccccc1",
    "CC(=O)O",  "CC(=O)OC",   "CCN",     "CCCN",     "C1CCCCC1",  "c1ccncc1"};

const std::vector<const char*> SHAPE_SMILES{
    "CC(=O)Oc1ccccc1C(=O)O", "CCN(CC)CCOC(=O)c1ccccc1", "c1ccc2ccccc2c1",
    "CCCCCCO", "OC(=O)c1ccccc1", "CC(C)Cc1ccc(cc1)C(C)C(=O)O"};

const std::vector<const char*> MCS_SMILES{
    "c1ccccc1", "Cc1ccccc1", "C1CCCCC1", "Oc1ccccc1",
    "CN1CC[C@]23c4c5ccc(O)c4O[C@H]2[C@@H](O)C=C[C@H]3[C@H]1C5",
    "CC1(C)S[C@@H]2[C@H](NC(=O)Cc3ccccc3)C(=O)N2[C@H]1C(=O)O"};

// Drifts with its own clone's history: the defect the checks exist to catch.
class HistoryComparison : public PairwiseComparison {
public:
    explicit HistoryComparison(size_t n) : n_(n) {}

    double Compare(size_t, size_t) override {
        return 1.0 + 1e-9 * static_cast<double>(calls_++);
    }
    std::unique_ptr<PairwiseComparison> Clone() const override {
        return std::make_unique<HistoryComparison>(n_);
    }
    size_t Size() const override { return n_; }
    std::string ComparisonName() const override { return "history"; }

private:
    size_t n_;
    size_t calls_ = 0;
};

}  // namespace

// Without this, a harness that compared nothing would pass every family.
TEST(ComparisonRepeatabilityTest, TheChecksCatchAHistoryDependentComparison) {
    HistoryComparison drifting(4);
    ::testing::TestPartResultArray failures;
    {
        ::testing::ScopedFakeTestPartResultReporter reporter(
            ::testing::ScopedFakeTestPartResultReporter::INTERCEPT_ONLY_CURRENT_THREAD,
            &failures);
        ExpectRepeatable(drifting, {1}, {4096});
    }
    EXPECT_GT(failures.size(), 0);
}

TEST(ComparisonRepeatabilityTest, Fingerprint) {
    std::vector<OEChem::OEGraphMol> mols = GraphMols(FINGERPRINT_SMILES);
    FingerprintComparison comparison(Pointers(mols));
    ExpectRepeatable(comparison, THREADS, CHUNKS);
}

TEST(ComparisonRepeatabilityTest, Descriptor) {
    std::vector<OEChem::OEGraphMol> mols = GraphMols(FINGERPRINT_SMILES);
    DescriptorComparison comparison(Pointers(mols));
    ExpectRepeatable(comparison, THREADS, CHUNKS);
}

TEST(ComparisonRepeatabilityTest, ROCS) {
    ROCSComparison comparison(Conformers(SHAPE_SMILES));
    ExpectRepeatable(comparison, THREADS, CHUNKS);
}

TEST(ComparisonRepeatabilityTest, RMSDInFrameAndOverlaid) {
    const std::vector<std::shared_ptr<OEChem::OEMol>> poses =
        Poses("CCCCCCOc1ccccc1", 6);
    ASSERT_GE(poses.size(), 3u);
    RMSDComparison in_frame(poses);
    ExpectRepeatable(in_frame, THREADS, CHUNKS);
    RMSDOptions overlay;
    overlay.overlay = true;
    RMSDComparison overlaid(poses, overlay);
    ExpectRepeatable(overlaid, THREADS, CHUNKS);
}

TEST(ComparisonRepeatabilityTest, MCS) {
    MCSComparison comparison(Topological(MCS_SMILES));
    ExpectRepeatable(comparison, THREADS, CHUNKS);
}

// Only two design units ship as assets, so the first is read a second time
// as a third item: three pairs give each reading a history of other pairs,
// and two threads at one pair per chunk hand chunks to different clones.
// Kept that small because a SiteHopper comparison costs almost a second.
TEST(ComparisonRepeatabilityTest, SuperposeAndSiteHopper) {
    std::vector<std::shared_ptr<OEBio::OEDesignUnit>> dus;
    for (const char* name :
         {"spruce_5FQD_1_5FQD_1-ALIGNED_BC__DU__LVY_B-1438.oedu",
          "spruce_8G66_1_8G66_1-ALIGNED_BC__DU__YOT_B-502.oedu",
          "spruce_5FQD_1_5FQD_1-ALIGNED_BC__DU__LVY_B-1438.oedu"}) {
        auto du = std::make_shared<OEBio::OEDesignUnit>();
        const std::string path = std::string(TEST_ASSETS_DIR) + "/" + name;
        if (!OEBio::OEReadDesignUnit(path, *du)) {
            GTEST_SKIP() << "Cannot read test asset: " << path;
        }
        dus.push_back(du);
    }
    for (const SuperposeMethod method :
         {SuperposeMethod::GlobalCarbonAlpha, SuperposeMethod::SiteHopper}) {
        SCOPED_TRACE(::testing::Message() << "method " << static_cast<int>(method));
        SuperposeOptions options;
        options.method = method;
        SuperposeComparison comparison(dus, options);
        ExpectRepeatable(comparison, {1, 2}, {1});
    }
}
