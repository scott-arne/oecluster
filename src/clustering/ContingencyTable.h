/**
 * @file ContingencyTable.h
 * @brief Contingency-table construction for the partition-agreement metrics.
 *
 * Header-only and inline, matching InternalIndices.h. These live outside
 * PartitionAgreement.cpp's anonymous namespace so the C++ tests can drive them
 * directly: an interning mistake, an unsorted cell list, or a multiply that
 * wraps at 32 bits all fail by returning a plausible number rather than by
 * crashing.
 */

#ifndef OECLUSTER_SRC_CLUSTERING_CONTINGENCYTABLE_H
#define OECLUSTER_SRC_CLUSTERING_CONTINGENCYTABLE_H

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

#include "oecluster/clustering/ClusterTypes.h"
#include "oecluster/clustering/PartitionAgreement.h"

namespace OECluster::detail {

/// The sparse contingency table of two interned labelings.
struct ContingencyTable {
    /// Sorted by (row, col), so the cell list is identical between two runs on
    /// identical input. This does not make accumulation order independent of
    /// the *sample* order — the ids follow first appearance — so `score()`
    /// sorts its own terms before summing them.
    struct Cell {
        uint32_t row;
        uint32_t col;
        uint64_t count;
    };

    uint64_t num_samples = 0;
    std::vector<uint64_t> marginals_a;
    std::vector<uint64_t> marginals_b;
    std::vector<Cell> cells;
};

/// True when this label means "no cluster" on its own side. Any negative
/// integer is noise, matching ClusterReport and labels_to_clusters; the empty
/// string is noise on the scaffold side.
inline bool is_noise(ClusterLabel label) { return label < 0; }
inline bool is_noise(const std::string& label) { return label.empty(); }

/// Pass one, run for each side before any id is assigned: ORs this side's
/// noise samples into `drop`. Only NoiseHandling::Excluded sets any bit; the
/// other two readings keep every sample. Because both sides OR into the same
/// mask, `drop` holds the union of the two sides' exclusions once both calls
/// have returned -- a sample noisy on either side is dropped. Discovery is
/// separated from interning precisely because the combined mask must exist
/// before the first id is handed out.
template <typename Label>
inline void mark_excluded(const std::vector<Label>& labels,
                          NoiseHandling noise_handling,
                          std::vector<bool>& drop) {
    if (noise_handling != NoiseHandling::Excluded) {
        return;
    }
    for (size_t i = 0; i < labels.size(); ++i) {
        if (is_noise(labels[i])) {
            drop[i] = true;
        }
    }
}

/// Pass two, run for each side with the finalized mask: maps the surviving
/// labels to dense ids 0..K-1 in first-appearance order, skipping every sample
/// `drop` marks. Returns one id per surviving sample in input order and
/// reports K through `num_ids`. Templated on the label type because the two
/// sides key on different types -- ClusterLabel for a labeling, std::string
/// for a scaffold annotation -- and a single key type cannot serve both.
template <typename Label>
inline std::vector<uint32_t> intern_side(const std::vector<Label>& labels,
                                         NoiseHandling noise_handling,
                                         const std::vector<bool>& drop,
                                         uint32_t& num_ids) {
    std::vector<uint32_t> ids;
    ids.reserve(labels.size());
    std::unordered_map<Label, uint32_t> seen;
    uint32_t next_id = 0;
    uint32_t grouped_id = 0;
    bool grouped_assigned = false;

    for (size_t i = 0; i < labels.size(); ++i) {
        if (drop[i]) {
            continue;
        }
        const Label& label = labels[i];
        if (is_noise(label)) {
            // Excluded has already removed every noise sample through `drop`,
            // so only Singletons and Grouped reach here.
            if (noise_handling == NoiseHandling::Grouped) {
                if (!grouped_assigned) {
                    grouped_id = next_id++;
                    grouped_assigned = true;
                }
                ids.push_back(grouped_id);
            } else {
                ids.push_back(next_id++);
            }
            continue;
        }
        auto it = seen.find(label);
        if (it == seen.end()) {
            it = seen.emplace(label, next_id++).first;
        }
        ids.push_back(it->second);
    }

    num_ids = next_id;
    return ids;
}

/// Shared body of the two build_contingency overloads; they differ only in
/// which is_noise the templates above resolve to.
template <typename LabelA, typename LabelB>
inline ContingencyTable build_contingency_impl(const std::vector<LabelA>& a,
                                               const std::vector<LabelB>& b,
                                               NoiseHandling noise_handling) {
    // The public entry points validate this before calling, so the invariant
    // holds in production; the C++ tests drive the build_contingency overloads
    // directly and would index `drop` and `ids_b` out of bounds without it.
    // Both the library and the test binary build Release with -DNDEBUG, so
    // this documents the precondition and traps it under a debug build rather
    // than guarding the suite that normally runs.
    assert(a.size() == b.size() && "build_contingency requires equal lengths");

    std::vector<bool> drop(a.size(), false);
    mark_excluded(a, noise_handling, drop);
    mark_excluded(b, noise_handling, drop);

    uint32_t num_a = 0;
    uint32_t num_b = 0;
    const std::vector<uint32_t> ids_a =
        intern_side(a, noise_handling, drop, num_a);
    const std::vector<uint32_t> ids_b =
        intern_side(b, noise_handling, drop, num_b);

    ContingencyTable table;
    table.num_samples = static_cast<uint64_t>(ids_a.size());
    table.marginals_a.assign(num_a, 0);
    table.marginals_b.assign(num_b, 0);

    // The key packing is exact because intern_side yields uint32_t ids.
    std::unordered_map<uint64_t, uint64_t> counts;
    counts.reserve(ids_a.size());
    for (size_t i = 0; i < ids_a.size(); ++i) {
        const uint32_t row = ids_a[i];
        const uint32_t col = ids_b[i];
        ++table.marginals_a[row];
        ++table.marginals_b[col];
        ++counts[(static_cast<uint64_t>(row) << 32) | col];
    }

    table.cells.reserve(counts.size());
    for (const auto& entry : counts) {
        table.cells.push_back(ContingencyTable::Cell{
            static_cast<uint32_t>(entry.first >> 32),
            static_cast<uint32_t>(entry.first & 0xFFFFFFFFu), entry.second});
    }

    // An unordered_map's iteration order is not reproducible across runs and
    // floating-point addition is not associative, so without this sort two runs
    // on identical input could differ in their last bits. It fixes the order
    // for one sample order only: interned ids follow first appearance, so a
    // permutation of the samples renumbers the clusters and reorders these
    // cells. Permutation invariance is established in `score()`, which sorts
    // its summation terms by the values being added rather than by id.
    std::sort(table.cells.begin(), table.cells.end(),
              [](const ContingencyTable::Cell& lhs,
                 const ContingencyTable::Cell& rhs) {
                  return lhs.row != rhs.row ? lhs.row < rhs.row
                                            : lhs.col < rhs.col;
              });
    return table;
}

inline ContingencyTable build_contingency(const std::vector<ClusterLabel>& a,
                                          const std::vector<ClusterLabel>& b,
                                          NoiseHandling noise_handling) {
    return build_contingency_impl(a, b, noise_handling);
}

inline ContingencyTable build_contingency(const std::vector<ClusterLabel>& a,
                                          const std::vector<std::string>& b,
                                          NoiseHandling noise_handling) {
    return build_contingency_impl(a, b, noise_handling);
}

/// One distinct cluster size and the number of clusters that have it.
struct SizeMultiplicity {
    uint64_t size;
    uint32_t multiplicity;
};

/// Collapse one side's marginals into its distinct-size histogram, strictly
/// ascending by `size` and with no repeated size. This is the cnt_a / cnt_b of
/// the expected-mutual-information sum, and the ascending order is what fixes
/// that sum's accumulation order.
inline std::vector<SizeMultiplicity> marginal_histogram(
    const std::vector<uint64_t>& marginals) {
    std::vector<uint64_t> sizes(marginals);
    std::sort(sizes.begin(), sizes.end());

    std::vector<SizeMultiplicity> histogram;
    for (uint64_t size : sizes) {
        if (!histogram.empty() && histogram.back().size == size) {
            ++histogram.back().multiplicity;
        } else {
            histogram.push_back(SizeMultiplicity{size, 1});
        }
    }
    return histogram;
}

/// Visit every (row size, column size) pair exactly once, outer loop over
/// `hist_a` and inner over `hist_b`, both ascending, calling
/// `visit(uint64_t u, uint64_t w, uint64_t multiplicity)` with `multiplicity`
/// equal to cnt_a[u] * cnt_b[w]. Factored out rather than written inline in
/// the AMI routine so the traversal itself is observable: the tests record the
/// callbacks and assert the whole sequence, which is the only exact check
/// available on a computation whose total is floating-point.
template <typename Visit>
inline void for_each_marginal_pair(const std::vector<SizeMultiplicity>& hist_a,
                                   const std::vector<SizeMultiplicity>& hist_b,
                                   Visit&& visit) {
    for (const SizeMultiplicity& entry_a : hist_a) {
        for (const SizeMultiplicity& entry_b : hist_b) {
            // Both factors are uint32_t and both reach 70000 on a supported
            // input, where the product is 4900000000 -- it wraps to 605032704
            // in 32-bit arithmetic and silently corrupts E[MI]. Widen first.
            const uint64_t multiplicity =
                static_cast<uint64_t>(entry_a.multiplicity) *
                static_cast<uint64_t>(entry_b.multiplicity);
            visit(entry_a.size, entry_b.size, multiplicity);
        }
    }
}

/// log(k!) for k in [0, n], by prefix-summing log(k). Costs 8*(n+1) bytes --
/// 800 KB at n = 100k -- and turns each inner term into nine table lookups and
/// eight adds instead of nine lgamma calls.
inline std::vector<double> log_factorials(uint64_t n) {
    std::vector<double> table(static_cast<size_t>(n) + 1, 0.0);
    for (uint64_t k = 2; k <= n; ++k) {
        table[static_cast<size_t>(k)] =
            table[static_cast<size_t>(k - 1)] + std::log(static_cast<double>(k));
    }
    return table;
}

/// The inner sum of Vinh et al. (2010) for one pair of marginal values. It
/// depends on (i, j) only through (u, w), which is what makes the grouping in
/// expected_mutual_information exact rather than approximate.
inline double expected_term(uint64_t u, uint64_t w, uint64_t n,
                            const std::vector<double>& logfact) {
    const auto lf = [&logfact](uint64_t k) {
        return logfact[static_cast<size_t>(k)];
    };
    const double total = static_cast<double>(n);
    // The hypergeometric support: a cell cannot be emptier than u + w - n, and
    // count = 0 contributes nothing because 0 * log(...) is 0.
    const uint64_t lower = (u + w > n) ? (u + w - n) : 1;
    const uint64_t upper = std::min(u, w);
    double sum = 0.0;
    for (uint64_t count = lower; count <= upper; ++count) {
        const double log_p = lf(u) + lf(w) + lf(n - u) + lf(n - w) - lf(n) -
                             lf(count) - lf(u - count) - lf(w - count) -
                             lf(n - u - w + count);
        sum += (static_cast<double>(count) / total) *
               std::log(total * static_cast<double>(count) /
                        (static_cast<double>(u) * static_cast<double>(w))) *
               std::exp(log_p);
    }
    return sum;
}

/**
 * @brief Expected mutual information under the hypergeometric model.
 *
 * The outer sums run over every (i, j) pair of marginals, including pairs
 * whose observed cell count is zero, so this cannot ride the sparse cell list.
 * Grouping equal marginal values is exact, not an approximation, and bounds the
 * work by the number of distinct cluster sizes rather than the cluster count --
 * two all-singleton sides collapse from N^2 pairs to one. This sum needs no
 * extra canonicalization: the histograms are built from the multiset of
 * cluster sizes and traversed in ascending order, both of which are unchanged
 * by a permutation of the samples.
 */
inline double expected_mutual_information(const ContingencyTable& table,
                                          const std::vector<double>& logfact) {
    double expected = 0.0;
    for_each_marginal_pair(
        marginal_histogram(table.marginals_a),
        marginal_histogram(table.marginals_b),
        [&](uint64_t u, uint64_t w, uint64_t multiplicity) {
            expected += static_cast<double>(multiplicity) *
                        expected_term(u, w, table.num_samples, logfact);
        });
    return expected;
}

}  // namespace OECluster::detail

#endif  // OECLUSTER_SRC_CLUSTERING_CONTINGENCYTABLE_H
