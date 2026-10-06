/**
 * @file mst_test_support.h
 * @brief Fixtures shared by the spanning-tree clustering tests.
 */

#ifndef OECLUSTER_TESTS_CPP_MST_TEST_SUPPORT_H
#define OECLUSTER_TESTS_CPP_MST_TEST_SUPPORT_H

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <oechem.h>

#include "oecluster/GateFacts.h"
#include "oecluster/PairwiseComparison.h"
#include "oecluster/StorageBackend.h"

#include "../../src/clustering/HDBSCANLinkage.h"

namespace mst_test {

/**
 * @brief The sequential Prim's algorithm HDBSCAN shipped through 5.19.0.
 *
 * Kept verbatim, apart from reading a condensed vector instead of storage, as
 * the oracle the shared kernel must reproduce edge for edge.
 */
inline std::vector<OECluster::detail::HDBSCANMSTEdge> LegacyPrim(
    size_t n, const std::vector<double>& condensed,
    const std::vector<double>& core_distances, double alpha) {
    auto distance = [&](size_t i, size_t j) {
        if (i > j) {
            std::swap(i, j);
        }
        return condensed[n * i - i * (i + 1) / 2 + j - i - 1];
    };
    if (n < 2) {
        return {};
    }
    std::vector<OECluster::detail::HDBSCANMSTEdge> mst;
    mst.reserve(n - 1);
    std::vector<bool> in_tree(n, false);
    std::vector<double> min_reachability(n, std::numeric_limits<double>::infinity());
    std::vector<size_t> current_sources(n, 0);
    size_t current_node = 0;

    for (size_t step = 0; step < n - 1; ++step) {
        in_tree[current_node] = true;

        double new_reachability = std::numeric_limits<double>::infinity();
        size_t source_node = 0;
        size_t new_node = 0;

        for (size_t candidate = 0; candidate < n; ++candidate) {
            if (in_tree[candidate]) {
                continue;
            }

            const double next_min_reach = min_reachability[candidate];
            const size_t next_source = current_sources[candidate];
            const double reach = std::max({core_distances[current_node],
                                           core_distances[candidate],
                                           distance(current_node, candidate) / alpha});

            if (reach < next_min_reach) {
                min_reachability[candidate] = reach;
                current_sources[candidate] = current_node;
                if (reach < new_reachability) {
                    new_reachability = reach;
                    source_node = current_node;
                    new_node = candidate;
                }
            } else if (next_min_reach < new_reachability) {
                new_reachability = next_min_reach;
                source_node = next_source;
                new_node = candidate;
            }
        }

        mst.push_back(OECluster::detail::HDBSCANMSTEdge{source_node, new_node,
                                                        new_reachability});
        current_node = new_node;
    }

    return mst;
}

inline bool SameEdges(const std::vector<OECluster::detail::HDBSCANMSTEdge>& a,
                      const std::vector<OECluster::detail::HDBSCANMSTEdge>& b) {
    if (a.size() != b.size()) {
        return false;
    }
    for (size_t k = 0; k < a.size(); ++k) {
        if (a[k].current_node != b[k].current_node || a[k].next_node != b[k].next_node ||
            a[k].distance != b[k].distance) {
            return false;
        }
    }
    return true;
}

/**
 * @brief A condensed table that counts every evaluation, per unordered pair.
 *
 * Clones share the counts, so a test can assert that a pass compared each
 * pair exactly once whatever the thread schedule. A reversed pair throws, which
 * pins the Compare(min, max) rule.
 */
class PairCountingComparison : public OECluster::PairwiseComparison {
public:
    PairCountingComparison(size_t n, std::vector<double> condensed,
                           OECluster::GateFacts facts = OECluster::GateFacts())
        : n_(n),
          condensed_(std::make_shared<const std::vector<double>>(std::move(condensed))),
          facts_(facts),
          counts_(std::make_shared<std::vector<std::atomic<size_t>>>(condensed_->size())) {}

    double Compare(size_t i, size_t j) override {
        if (i >= j) {
            throw std::logic_error("PairCountingComparison read (" + std::to_string(i) +
                                   ", " + std::to_string(j) + ")");
        }
        const size_t index = n_ * i - i * (i + 1) / 2 + j - i - 1;
        (*counts_)[index].fetch_add(1, std::memory_order_relaxed);
        return (*condensed_)[index];
    }

    OECluster::GateFacts Facts() const override { return facts_; }
    std::unique_ptr<OECluster::PairwiseComparison> Clone() const override {
        return std::make_unique<PairCountingComparison>(*this);
    }
    size_t Size() const override { return n_; }
    std::string ComparisonName() const override { return "pair-counting"; }

    // Evaluations of each pair, in condensed order.
    std::vector<size_t> Counts() const {
        std::vector<size_t> counts;
        for (const auto& count : *counts_) {
            counts.push_back(count.load());
        }
        return counts;
    }

    size_t Total() const {
        size_t total = 0;
        for (const size_t count : Counts()) {
            total += count;
        }
        return total;
    }

    void Reset() {
        for (auto& count : *counts_) {
            count.store(0);
        }
    }

private:
    size_t n_;
    std::shared_ptr<const std::vector<double>> condensed_;
    OECluster::GateFacts facts_;
    std::shared_ptr<std::vector<std::atomic<size_t>>> counts_;
};

/**
 * @brief An item count too large to allocate for, and no pair to read.
 *
 * Reaches a size check that must fire before any allocation or Compare call.
 */
class HugeComparison : public OECluster::PairwiseComparison {
public:
    explicit HugeComparison(size_t n) : n_(n) {}
    double Compare(size_t, size_t) override {
        throw std::logic_error("HugeComparison was compared");
    }
    std::unique_ptr<OECluster::PairwiseComparison> Clone() const override {
        return std::make_unique<HugeComparison>(*this);
    }
    size_t Size() const override { return n_; }
    std::string ComparisonName() const override { return "huge"; }

private:
    size_t n_;
};

/**
 * @brief A condensed table whose Compare throws on one pair.
 *
 * Proves that an exception raised inside a parallel pass reaches the caller,
 * and that a later run over a well-behaved comparison is unaffected.
 */
class ThrowingComparison : public PairCountingComparison {
public:
    ThrowingComparison(size_t n, std::vector<double> condensed, size_t i, size_t j)
        : PairCountingComparison(n, std::move(condensed)), i_(i), j_(j) {}

    double Compare(size_t i, size_t j) override {
        if (i == i_ && j == j_) {
            throw std::runtime_error("ThrowingComparison refused (" + std::to_string(i) +
                                     ", " + std::to_string(j) + ")");
        }
        return PairCountingComparison::Compare(i, j);
    }
    std::unique_ptr<OECluster::PairwiseComparison> Clone() const override {
        return std::make_unique<ThrowingComparison>(*this);
    }

private:
    size_t i_;
    size_t j_;
};

/**
 * @brief A comparison named "rocs" over a table, to reach the refusal cheaply.
 */
class NamedRocsComparison : public PairCountingComparison {
public:
    using PairCountingComparison::PairCountingComparison;
    std::unique_ptr<OECluster::PairwiseComparison> Clone() const override {
        return std::make_unique<NamedRocsComparison>(*this);
    }
    std::string ComparisonName() const override { return "rocs"; }
};

// Small molecules whose fingerprint distances tie often.
inline const std::vector<const char*>& FingerprintSmiles() {
    static const std::vector<const char*> smiles{
        "CCO",       "CCCO",     "CCCCO",       "c1ccccc1",    "Cc1ccccc1", "CCc1ccccc1",
        "CC(=O)O",   "CC(=O)OC", "CCN",         "CCCN",        "C1CCCCC1",  "c1ccncc1",
        "CCCCCO",    "c1ccc2ccccc2c1",          "CC(C)O",      "OCCO",      "CCOC(=O)C",
        "c1ccc(O)cc1", "c1ccc(N)cc1", "Clc1ccccc1", "CC(C)(C)O", "CCCCCC"};
    return smiles;
}

inline std::vector<OEChem::OEGraphMol> FingerprintMolecules() {
    std::vector<OEChem::OEGraphMol> mols;
    for (const char* smi : FingerprintSmiles()) {
        mols.emplace_back();
        if (!OEChem::OESmilesToMol(mols.back(), smi)) {
            throw std::runtime_error(std::string("unparseable SMILES ") + smi);
        }
    }
    return mols;
}

inline std::vector<OEChem::OEMolBase*> Pointers(std::vector<OEChem::OEGraphMol>& mols) {
    std::vector<OEChem::OEMolBase*> pointers;
    for (auto& mol : mols) {
        pointers.push_back(&static_cast<OEChem::OEMolBase&>(mol));
    }
    return pointers;
}

}  // namespace mst_test

#endif  // OECLUSTER_TESTS_CPP_MST_TEST_SUPPORT_H
