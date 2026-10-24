/**
 * @file test_worker_gate.cpp
 * @brief Every site that starts comparison workers honours the pool-mode gate.
 *
 * Under an unsafe memory-pool mode, a comparison that runs OpenEye code must
 * be scored by one thread at a time. The measure is the peak number of Compare
 * calls in flight: a gated call can legitimately compare on a worker thread in
 * one phase and inline on the calling thread in another, so a count of
 * distinct thread ids would not be the right measure.
 */

#include <gtest/gtest.h>

#include <algorithm>
#include <atomic>
#include <chrono>
#include <memory>
#include <mutex>
#include <set>
#include <string>
#include <thread>
#include <vector>

#include <oesystem.h>

#include "oecluster/CDist.h"
#include "oecluster/PDist.h"
#include "oecluster/PoolMode.h"
#include "oecluster/StorageBackend.h"
#include "oecluster/clustering/Agglomerative.h"
#include "oecluster/clustering/DiversitySelection.h"
#include "oecluster/comparisons/FingerprintComparison.h"

#include <oechem.h>

using namespace OECluster;

namespace {

const unsigned int SAFE = OESystem::OEMemPoolMode::Default;
const unsigned int UNSAFE = OESystem::OEMemPoolMode::UnboundedCache;

struct Overlap {
    std::atomic<int> in_flight{0};
    std::atomic<int> peak{0};
};

/// Distances |i - j|, a metric with zero self-distance. Compare counts the
/// calls in flight across every clone and keeps the peak, and sleeps briefly
/// so that workers which run at all run at once.
class OverlapComparison : public PairwiseComparison {
public:
    OverlapComparison(size_t n, bool runs_openeye)
        : n_(n), runs_openeye_(runs_openeye), overlap_(std::make_shared<Overlap>()) {}

    double Compare(size_t i, size_t j) override {
        const int now = overlap_->in_flight.fetch_add(1) + 1;
        int seen = overlap_->peak.load();
        while (now > seen && !overlap_->peak.compare_exchange_weak(seen, now)) {
        }
        std::this_thread::sleep_for(std::chrono::microseconds(500));
        overlap_->in_flight.fetch_sub(1);
        return static_cast<double>(i > j ? i - j : j - i);
    }

    GateFacts Facts() const override {
        GateFacts facts;
        facts.is_distance = Capability::Yes;
        facts.zero_self = Capability::Yes;
        facts.triangle = Capability::Yes;
        return facts;
    }

    bool RunsOpenEyeCode() const override { return runs_openeye_; }

    std::unique_ptr<PairwiseComparison> Clone() const override {
        return std::make_unique<OverlapComparison>(*this);
    }

    size_t Size() const override { return n_; }
    std::string ComparisonName() const override { return "overlap"; }
    int Peak() const { return overlap_->peak.load(); }

private:
    size_t n_;
    bool runs_openeye_;
    std::shared_ptr<Overlap> overlap_;
};

/// Records the worker count each bulk hook receives, then declines.
class BulkRecordingComparison : public OverlapComparison {
public:
    BulkRecordingComparison() : OverlapComparison(4, true) {}

    bool TryPDist(StorageBackend&, const PDistOptions& options) override {
        pdist_threads = options.num_threads;
        return false;
    }

    bool TryCDist(size_t, double*, const CDistOptions& options) override {
        cdist_threads = options.num_threads;
        return false;
    }

    size_t pdist_threads = 99;
    size_t cdist_threads = 99;
};

int pdist_peak(unsigned int mode, bool runs_openeye) {
    CallerPoolModeScope scope(mode);
    OverlapComparison comparison(24, runs_openeye);
    DenseStorage storage(24);
    PDistOptions options;
    options.num_threads = 4;
    options.chunk_size = 8;
    pdist(comparison, storage, options);
    return comparison.Peak();
}

int cdist_peak(unsigned int mode) {
    CallerPoolModeScope scope(mode);
    OverlapComparison comparison(24, true);
    std::vector<double> output(12 * 12);
    CDistOptions options;
    options.num_threads = 4;
    options.chunk_size = 8;
    cdist(comparison, 12, output.data(), options);
    return comparison.Peak();
}

int maxmin_peak(unsigned int mode) {
    CallerPoolModeScope scope(mode);
    OverlapComparison comparison(24, true);
    MaxMinOptions options;
    options.count = 4;
    options.num_threads = 4;
    options.chunk_size = 4;
    maxmin_select(comparison, options);
    return comparison.Peak();
}

int prim_peak(unsigned int mode) {
    CallerPoolModeScope scope(mode);
    OverlapComparison comparison(48, true);
    AgglomerativeOptions options;
    options.linkage = AgglomerativeLinkageMethod::Single;
    options.n_clusters = 2;
    options.num_threads = 4;
    agglomerative_cluster(comparison, options);
    return comparison.Peak();
}

}  // namespace

TEST(WorkerGateTest, PdistComparesOneAtATimeUnderAnUnsafePool) {
    EXPECT_EQ(pdist_peak(UNSAFE, true), 1);
    EXPECT_GT(pdist_peak(SAFE, true), 1);
}

TEST(WorkerGateTest, CdistComparesOneAtATimeUnderAnUnsafePool) {
    EXPECT_EQ(cdist_peak(UNSAFE), 1);
    EXPECT_GT(cdist_peak(SAFE), 1);
}

TEST(WorkerGateTest, ChunkedComparisonsCompareOneAtATimeUnderAnUnsafePool) {
    EXPECT_EQ(maxmin_peak(UNSAFE), 1);
    EXPECT_GT(maxmin_peak(SAFE), 1);
}

TEST(WorkerGateTest, PrimComparesOneAtATimeUnderAnUnsafePool) {
    EXPECT_EQ(prim_peak(UNSAFE), 1);
    EXPECT_GT(prim_peak(SAFE), 1);
}

TEST(WorkerGateTest, ScoringThatRunsNoOpenEyeCodeStaysParallel) {
    EXPECT_GT(pdist_peak(UNSAFE, false), 1);
}

/// A real FingerprintComparison that records the worker count its own bulk
/// hook receives, then runs that hook.
class RecordingFingerprint : public FingerprintComparison {
public:
    using FingerprintComparison::FingerprintComparison;

    bool TryPDist(StorageBackend& storage, const PDistOptions& options) override {
        pdist_threads = options.num_threads;
        return FingerprintComparison::TryPDist(storage, options);
    }

    bool TryCDist(size_t n_a, double* output, const CDistOptions& options) override {
        cdist_threads = options.num_threads;
        return FingerprintComparison::TryCDist(n_a, output, options);
    }

    size_t pdist_threads = 99;
    size_t cdist_threads = 99;
};

TEST(WorkerGateTest, AFingerprintComparisonKeepsItsRequestUnderAnUnsafePool) {
    // Fingerprint scoring runs no OpenEye code, so its OEFP kernel gets the
    // caller's request whatever the mode. How OEFP then threads that request
    // is OEFP's contract; what oecluster owes is the count it hands over.
    std::vector<OEChem::OEGraphMol> owned(6);
    const char* smiles[] = {"c1ccccc1", "Cc1ccccc1", "CCO", "c1ccncc1",
                            "CC(=O)O", "c1ccc2ccccc2c1"};
    std::vector<OEChem::OEMolBase*> mols;
    for (size_t i = 0; i < owned.size(); ++i) {
        OEChem::OESmilesToMol(owned[i], smiles[i]);
        mols.push_back(&static_cast<OEChem::OEMolBase&>(owned[i]));
    }
    CallerPoolModeScope scope(UNSAFE);
    RecordingFingerprint fingerprint(mols);

    DenseStorage storage(6);
    PDistOptions pdist_options;
    pdist_options.num_threads = 8;
    pdist(fingerprint, storage, pdist_options);
    EXPECT_EQ(fingerprint.pdist_threads, 8u);

    std::vector<double> output(3 * 3);
    CDistOptions cdist_options;
    cdist_options.num_threads = 8;
    cdist(fingerprint, 3, output.data(), cdist_options);
    EXPECT_EQ(fingerprint.cdist_threads, 8u);
}

TEST(WorkerGateTest, TheBulkHooksReceiveTheCappedCount) {
    for (const auto& [mode, expected] :
         std::vector<std::pair<unsigned int, size_t>>{{UNSAFE, 1}, {SAFE, 8}}) {
        SCOPED_TRACE(testing::Message() << "mode = " << mode);
        CallerPoolModeScope scope(mode);
        BulkRecordingComparison comparison;
        DenseStorage storage(4);
        PDistOptions pdist_options;
        pdist_options.num_threads = 8;
        pdist(comparison, storage, pdist_options);
        EXPECT_EQ(comparison.pdist_threads, expected);

        std::vector<double> output(2 * 2);
        CDistOptions cdist_options;
        cdist_options.num_threads = 8;
        cdist(comparison, 2, output.data(), cdist_options);
        EXPECT_EQ(comparison.cdist_threads, expected);
    }
}

TEST(WorkerGateTest, ProgressStaysOnTheCallingThreadWhenTheCallerAsks) {
    for (const bool on_caller : {true, false}) {
        SCOPED_TRACE(testing::Message() << "callbacks_on_calling_thread = " << on_caller);
        CallerPoolModeScope scope(SAFE, on_caller);
        const std::thread::id caller = std::this_thread::get_id();
        std::mutex mutex;
        std::set<std::thread::id> seen;
        auto record = [&](size_t, size_t) {
            std::lock_guard<std::mutex> lock(mutex);
            seen.insert(std::this_thread::get_id());
        };

        OverlapComparison comparison(24, false);
        DenseStorage storage(24);
        PDistOptions pdist_options;
        pdist_options.num_threads = 4;
        pdist_options.chunk_size = 8;
        pdist_options.progress = record;
        pdist(comparison, storage, pdist_options);

        std::vector<double> output(12 * 12);
        CDistOptions cdist_options;
        cdist_options.num_threads = 4;
        cdist_options.chunk_size = 8;
        cdist_options.progress = record;
        cdist(comparison, 12, output.data(), cdist_options);

        if (on_caller) {
            EXPECT_EQ(seen, std::set<std::thread::id>{caller});
        } else {
            EXPECT_EQ(seen.count(caller), 0u);
        }
    }
}
