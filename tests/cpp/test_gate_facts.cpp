#include <gtest/gtest.h>
#include "oecluster/GateFacts.h"
#include "oecluster/PairwiseComparison.h"

using namespace OECluster;

/// Minimal comparison that does not override Facts(), exercising the base default.
class SilentComparison : public PairwiseComparison {
public:
    double Compare(size_t, size_t) override { return 0.0; }
    std::unique_ptr<PairwiseComparison> Clone() const override {
        return std::make_unique<SilentComparison>();
    }
    size_t Size() const override { return 0; }
    std::string ComparisonName() const override { return "silent"; }
};

/// Comparison that stamps concrete facts, exercising the override path.
class StampedComparison : public SilentComparison {
public:
    GateFacts Facts() const override {
        GateFacts facts;
        facts.zero_self = Capability::Yes;
        facts.triangle = Capability::No;
        facts.data_integrity = DataIntegrity::SubsetScored;
        return facts;
    }
};

TEST(GateFactsTest, DefaultsAreUnknownAndComplete) {
    GateFacts facts;
    EXPECT_EQ(facts.zero_self, Capability::Unknown);
    EXPECT_EQ(facts.triangle, Capability::Unknown);
    EXPECT_EQ(facts.data_integrity, DataIntegrity::Complete);
}

TEST(GateFactsTest, BaseComparisonReportsUnknown) {
    SilentComparison comparison;
    const GateFacts facts = comparison.Facts();
    EXPECT_EQ(facts.zero_self, Capability::Unknown);
    EXPECT_EQ(facts.triangle, Capability::Unknown);
    EXPECT_EQ(facts.data_integrity, DataIntegrity::Complete);
}

TEST(GateFactsTest, OverrideIsVisibleThroughBasePointer) {
    StampedComparison stamped;
    PairwiseComparison* base = &stamped;
    const GateFacts facts = base->Facts();
    EXPECT_EQ(facts.zero_self, Capability::Yes);
    EXPECT_EQ(facts.triangle, Capability::No);
    EXPECT_EQ(facts.data_integrity, DataIntegrity::SubsetScored);
}
