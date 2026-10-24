/**
 * @file test_pool_mode.cpp
 * @brief The memory-pool mode a call runs under, and the worker counts it allows.
 */

#include <gtest/gtest.h>

#include <memory>
#include <optional>
#include <string>
#include <thread>
#include <vector>

#include <oechem.h>
#include <oesystem.h>

#include "oecluster/PairwiseComparison.h"
#include "oecluster/PoolMode.h"
#include "oecluster/comparisons/DescriptorComparison.h"
#include "oecluster/comparisons/FingerprintComparison.h"

using OECluster::callbacks_on_calling_thread;
using OECluster::CallerPoolModeScope;
using OECluster::effective_pool_mode;
using OECluster::pool_is_thread_safe;
using OECluster::resolve_pool_mode;
using OECluster::safe_worker_count;

namespace {

const unsigned int SAFE = OESystem::OEMemPoolMode::Default;
const unsigned int UNSAFE = OESystem::OEMemPoolMode::UnboundedCache;

/// A comparison that answers RunsOpenEyeCode as told.
class FixedComparison : public OECluster::PairwiseComparison {
public:
    explicit FixedComparison(bool runs_openeye) : runs_openeye_(runs_openeye) {}
    double Compare(size_t, size_t) override { return 0.0; }
    std::unique_ptr<PairwiseComparison> Clone() const override {
        return std::make_unique<FixedComparison>(*this);
    }
    size_t Size() const override { return 2; }
    std::string ComparisonName() const override { return "fixed"; }
    bool RunsOpenEyeCode() const override { return runs_openeye_; }

private:
    bool runs_openeye_;
};

/// A comparison that does not override RunsOpenEyeCode.
class SilentComparison : public OECluster::PairwiseComparison {
public:
    double Compare(size_t, size_t) override { return 0.0; }
    std::unique_ptr<PairwiseComparison> Clone() const override {
        return std::make_unique<SilentComparison>(*this);
    }
    size_t Size() const override { return 2; }
    std::string ComparisonName() const override { return "silent"; }
};

}  // namespace

TEST(PoolModeTest, TheSuiteRunsUnderAThreadSafeNativePool) {
    // Every safe-mode expectation below assumes the native default.
    EXPECT_TRUE(pool_is_thread_safe(OESystem::OEGetMemPoolMode()));
}

TEST(PoolModeTest, RecognizesOnlyTheThreadSafeModes) {
    // Every single bit OEMemPoolMode defines as of 2026.1.0. The members are
    // namespace-scoped const unsigned int rather than enumerators, so a vendor
    // addition compiles silently and simply goes unlisted -- this table is a
    // record of what was checked, not a compile-time guard.
    EXPECT_TRUE(pool_is_thread_safe(OESystem::OEMemPoolMode::Mutexed));
    EXPECT_TRUE(pool_is_thread_safe(OESystem::OEMemPoolMode::ThreadLocal));
    EXPECT_TRUE(pool_is_thread_safe(OESystem::OEMemPoolMode::Default));

    // SingleThreaded is 0, so it is the absence of a safe bit rather than a
    // flag: asserting it is unsafe is asserting the predicate's default.
    EXPECT_FALSE(pool_is_thread_safe(OESystem::OEMemPoolMode::SingleThreaded));
    // The cache and allocator-strategy bits say nothing about thread safety on
    // their own. Spinlocked is treated as unsafe deliberately: the SDK header
    // describes it as "Spinlock - for Java threading", which names an intended
    // consumer rather than a guarantee for this caller, and the cost of being
    // wrong in this direction is lost parallelism rather than corruption.
    EXPECT_FALSE(pool_is_thread_safe(OESystem::OEMemPoolMode::BoundedCache));
    EXPECT_FALSE(pool_is_thread_safe(OESystem::OEMemPoolMode::UnboundedCache));
    EXPECT_FALSE(pool_is_thread_safe(OESystem::OEMemPoolMode::System));
    EXPECT_FALSE(pool_is_thread_safe(OESystem::OEMemPoolMode::Spinlocked));

    // A safe bit stays safe when a cache bit joins it, which is the only
    // combination Default itself relies on.
    EXPECT_TRUE(pool_is_thread_safe(OESystem::OEMemPoolMode::Mutexed |
                                    OESystem::OEMemPoolMode::BoundedCache));
    EXPECT_FALSE(pool_is_thread_safe(OESystem::OEMemPoolMode::System |
                                     OESystem::OEMemPoolMode::UnboundedCache));
}

TEST(PoolModeTest, EveryCombinationOfTheModeBitsIsClassifiedByItsSafeBits) {
    // The six documented bits: Mutexed 1, ThreadLocal 2, BoundedCache 4,
    // UnboundedCache 8, System 16, Spinlocked 32.
    const unsigned int safe_bits = OESystem::OEMemPoolMode::Mutexed |
                                   OESystem::OEMemPoolMode::ThreadLocal;
    for (unsigned int mode = 0; mode < 64u; ++mode) {
        SCOPED_TRACE(testing::Message() << "mode = " << mode);
        EXPECT_EQ(pool_is_thread_safe(mode), (mode & safe_bits) != 0u);
    }
}

TEST(PoolModeTest, WithNoScopeTheEffectiveModeIsTheNativeOne) {
    EXPECT_EQ(effective_pool_mode(), OESystem::OEGetMemPoolMode());
    EXPECT_FALSE(callbacks_on_calling_thread());
}

TEST(PoolModeTest, ScopesNestAndRestore) {
    {
        CallerPoolModeScope outer(UNSAFE);
        EXPECT_EQ(effective_pool_mode(), UNSAFE);
        EXPECT_FALSE(callbacks_on_calling_thread());
        {
            CallerPoolModeScope inner(SAFE, true);
            EXPECT_EQ(effective_pool_mode(), SAFE);
            EXPECT_TRUE(callbacks_on_calling_thread());
        }
        EXPECT_EQ(effective_pool_mode(), UNSAFE);
        EXPECT_FALSE(callbacks_on_calling_thread());
    }
    EXPECT_EQ(effective_pool_mode(), OESystem::OEGetMemPoolMode());
}

TEST(PoolModeTest, AnEmptyScopeLeavesTheOuterModeInForce) {
    CallerPoolModeScope outer(UNSAFE);
    {
        CallerPoolModeScope inner(std::nullopt, true);
        EXPECT_EQ(effective_pool_mode(), UNSAFE);
        EXPECT_TRUE(callbacks_on_calling_thread());
    }
    EXPECT_EQ(effective_pool_mode(), UNSAFE);
}

TEST(PoolModeTest, AnEmptyScopeAloneLeavesTheNativeMode) {
    CallerPoolModeScope scope(std::nullopt);
    EXPECT_EQ(effective_pool_mode(), OESystem::OEGetMemPoolMode());
}

TEST(PoolModeTest, AScopeBelongsToItsThread) {
    CallerPoolModeScope scope(UNSAFE, true);
    unsigned int seen_mode = 0;
    bool seen_callbacks = true;
    std::thread other([&] {
        seen_mode = effective_pool_mode();
        seen_callbacks = callbacks_on_calling_thread();
    });
    other.join();
    EXPECT_EQ(seen_mode, OESystem::OEGetMemPoolMode());
    EXPECT_FALSE(seen_callbacks);
}

TEST(PoolModeTest, ResolvePrefersTheCallerModeThenTheEffectiveOne) {
    EXPECT_EQ(resolve_pool_mode(UNSAFE), UNSAFE);
    EXPECT_EQ(resolve_pool_mode(SAFE), SAFE);
    EXPECT_EQ(resolve_pool_mode(std::nullopt), OESystem::OEGetMemPoolMode());
    // With no Python value, an unsafe outer scope (standing in for an unsafe
    // native mode, which one process cannot set twice) decides.
    CallerPoolModeScope outer(UNSAFE);
    EXPECT_EQ(resolve_pool_mode(std::nullopt), UNSAFE);
    EXPECT_EQ(resolve_pool_mode(SAFE), SAFE);
}

TEST(PoolModeTest, SafeWorkerCountCapsOnlyOpenEyeWorkUnderAnUnsafeMode) {
    const FixedComparison openeye(true);
    const FixedComparison plain(false);
    for (const size_t requested : {size_t{0}, size_t{1}, size_t{8}}) {
        SCOPED_TRACE(testing::Message() << "requested = " << requested);
        {
            CallerPoolModeScope scope(SAFE);
            EXPECT_EQ(safe_worker_count(openeye, requested), requested);
            EXPECT_EQ(safe_worker_count(plain, requested), requested);
        }
        {
            CallerPoolModeScope scope(UNSAFE);
            EXPECT_EQ(safe_worker_count(openeye, requested), 1u);
            EXPECT_EQ(safe_worker_count(plain, requested), requested);
        }
    }
}

TEST(PoolModeTest, AComparisonIsAssumedToRunOpenEyeCode) {
    EXPECT_TRUE(SilentComparison().RunsOpenEyeCode());
}

TEST(PoolModeTest, FingerprintAndDescriptorScoringRunNoOpenEyeCode) {
    std::vector<OEChem::OEGraphMol> owned(2);
    OEChem::OESmilesToMol(owned[0], "c1ccccc1");
    OEChem::OESmilesToMol(owned[1], "Cc1ccccc1");
    // OEGraphMol reaches OEMolBase through a conversion, not inheritance, as
    // test_descriptor_comparison.cpp:155 spells it.
    std::vector<OEChem::OEMolBase*> mols = {
        &static_cast<OEChem::OEMolBase&>(owned[0]),
        &static_cast<OEChem::OEMolBase&>(owned[1])};

    const OECluster::FingerprintComparison fingerprint(mols);
    EXPECT_FALSE(fingerprint.RunsOpenEyeCode());

    OECluster::DescriptorOptions options;
    options.metric = "euclidean";
    const OECluster::DescriptorComparison descriptor(mols, options);
    EXPECT_FALSE(descriptor.RunsOpenEyeCode());
}
