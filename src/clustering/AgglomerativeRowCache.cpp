/**
 * @file AgglomerativeRowCache.cpp
 * @brief Active-slot condensed matrix with cached row minima, on a step team.
 */

#include "AgglomerativeRowCache.h"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <functional>
#include <limits>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <thread>
#include <utility>

#include "StepTeam.h"
#include "ThresholdGraph.h"
#include "oecluster/ThreadPool.h"

namespace OECluster::detail {

namespace {

constexpr double INFINITE_DISTANCE = std::numeric_limits<double>::infinity();
constexpr size_t NO_SLOT = std::numeric_limits<size_t>::max();
constexpr size_t NO_NODE = std::numeric_limits<size_t>::max();

// Average linkage weighs by cluster sizes; Weighted uses unweighted 0.5 factor
// per the reference definition. Byte for byte the 5.20.0 expression: the
// operand order and the rounding of each product, the sum and the division are
// what make a merge height bit-identical, so this must not be rearranged.
double update_linkage_distance(
    AgglomerativeLinkageMethod linkage,
    double left_distance,
    double right_distance,
    size_t left_size,
    size_t right_size) {
    switch (linkage) {
        case AgglomerativeLinkageMethod::Single:
            return std::min(left_distance, right_distance);
        case AgglomerativeLinkageMethod::Complete:
            return std::max(left_distance, right_distance);
        case AgglomerativeLinkageMethod::Average:
            return ((static_cast<double>(left_size) * left_distance) +
                    (static_cast<double>(right_size) * right_distance)) /
                   static_cast<double>(left_size + right_size);
        case AgglomerativeLinkageMethod::Weighted:
            return 0.5 * (left_distance + right_distance);
    }

    throw std::invalid_argument("Unknown agglomerative linkage method");
}

// One active pair under the selection key: ascending (distance, lower node id,
// higher node id). The node ids are public ids, never slot indices: slots are
// recycled and ids are not, so a slot comparison would reorder ties the moment
// a merged cluster took a retired slot.
struct Candidate {
    double distance = INFINITE_DISTANCE;
    size_t left = NO_NODE;
    size_t right = NO_NODE;
    size_t left_slot = NO_SLOT;
    size_t right_slot = NO_SLOT;
};

bool precedes(const Candidate& lhs, const Candidate& rhs) {
    if (lhs.distance != rhs.distance) {
        return lhs.distance < rhs.distance;
    }
    if (lhs.left != rhs.left) {
        return lhs.left < rhs.left;
    }
    return lhs.right < rhs.right;
}

Candidate make_candidate(double distance, size_t slot_a, size_t node_a,
                         size_t slot_b, size_t node_b) {
    if (node_b < node_a) {
        std::swap(slot_a, slot_b);
        std::swap(node_a, node_b);
    }
    return Candidate{distance, node_a, node_b, slot_a, slot_b};
}

// Within one row the canonical key collapses to (distance, partner node id).
// Row x compares (min(nx, ny), max(nx, ny)) across partners y with nx fixed:
// a partner below nx wins on the left element, a partner above nx loses on it,
// and within either group the smaller ny wins on whichever element differs. So
// the smaller partner id wins in every case, and a row scan needs no swap.
bool row_precedes(double distance, size_t partner_node, double best_distance,
                  size_t best_partner_node) {
    return distance < best_distance ||
           (distance == best_distance && partner_node < best_partner_node);
}

/** @brief A row's cached minimum: the nearest active cluster and its distance. */
struct RowBest {
    double distance = INFINITE_DISTANCE;
    size_t partner_slot = NO_SLOT;
};

// What one participant contributes to a merge step. Padded so that
// participants publishing at the same time do not share a cache line.
struct alignas(64) StepResult {
    // Row `kept`'s minimum over the rows this participant updated, which is
    // where every one of that row's new distances was produced.
    double kept_distance = INFINITE_DISTANCE;
    size_t kept_partner = NO_SLOT;
    size_t kept_partner_node = NO_NODE;
    // The smallest row minimum among those rows, under the full key.
    Candidate global;
    uint64_t rescans = 0;
    uint64_t rescanned_slots = 0;
};

/**
 * @brief Condensed indexing over the N slots, with the row bases precomputed.
 *
 * detail::condensed_index() validates and throws, which a loop running once
 * per active pair per merge cannot afford; the row base table also removes the
 * multiply and divide from the inner loop.
 */
class SlotIndex {
public:
    explicit SlotIndex(size_t n)
        : base_(n) {
        for (size_t i = 0; i < n; ++i) {
            base_[i] = n * i - i * (i + 1) / 2;
        }
    }

    size_t Base(size_t row) const { return base_[row]; }

    size_t operator()(size_t i, size_t j) const {
        const size_t low = i < j ? i : j;
        const size_t high = i < j ? j : i;
        return base_[low] + high - low - 1;
    }

private:
    std::vector<size_t> base_;
};

// The input and the workspace share the condensed layout over N, so the copy
// is index-for-index. Reading each value once, in order, is what lets a
// memory-mapped input stay evictable afterwards.
void copy_input(const StorageBackend& storage, std::vector<double>& work,
                const RowCacheOptions& options, ThreadPool& pool,
                const SlotIndex& index) {
    const size_t n = storage.NumSamples();
    const double* data = storage.Data();
    // Capped at the row count so a near-SIZE_MAX chunk_size cannot overflow
    // ParallelFor's ceiling arithmetic and silently skip the copy.
    const size_t chunk = std::min(options.chunk_size, n);
    pool.ParallelFor(0, n, chunk, [&](size_t begin, size_t end) {
        for (size_t i = begin; i < end; ++i) {
            if (i + 1 >= n) {
                continue;
            }
            size_t k = index.Base(i);
            for (size_t j = i + 1; j < n; ++j, ++k) {
                const double distance = data[k];
                if (!std::isfinite(distance)) {
                    throw non_finite_distance_error(options.caller, i, j);
                }
                work[k] = distance;
            }
        }
    });
}

// Row minima before the first merge, when slot i still holds leaf node i and
// the partner node id is the partner slot.
void initial_row_minima(size_t n, const std::vector<double>& work,
                        std::vector<RowBest>& best, size_t chunk,
                        ThreadPool& pool, const SlotIndex& index) {
    pool.ParallelFor(0, n, chunk, [&](size_t begin, size_t end) {
        for (size_t x = begin; x < end; ++x) {
            double row_distance = INFINITE_DISTANCE;
            size_t row_partner = NO_SLOT;
            for (size_t y = 0; y < x; ++y) {
                const double distance = work[index.Base(y) + x - y - 1];
                if (row_precedes(distance, y, row_distance, row_partner)) {
                    row_distance = distance;
                    row_partner = y;
                }
            }
            size_t k = index.Base(x);
            for (size_t y = x + 1; y < n; ++y, ++k) {
                if (row_precedes(work[k], y, row_distance, row_partner)) {
                    row_distance = work[k];
                    row_partner = y;
                }
            }
            best[x] = RowBest{row_distance, row_partner};
        }
    });
}

}  // namespace

size_t row_cache_participants(size_t num_threads, size_t n, size_t hardware) {
    return resolve_participants(num_threads, n,
                                std::min(hardware, ROW_CACHE_DEFAULT_PARTICIPANTS));
}

LinkageTree agglomerative_row_cache(const StorageBackend& storage,
                                    const RowCacheOptions& options) {
    LinkageTree tree;
    const size_t n = storage.NumSamples();
    if (n < 2) {
        return tree;
    }

    const SlotIndex index(n);
    std::vector<double> work(n * (n - 1) / 2);

    std::vector<size_t> slot_node(n);
    std::iota(slot_node.begin(), slot_node.end(), size_t{0});
    std::vector<size_t> cluster_size(n, 1);
    std::vector<RowBest> best(n);
    // The active slots, compacted so a step divides a contiguous range, and
    // the inverse map that retires one in constant time. Their order is free:
    // every reduction over them runs under a total order, so the winner does
    // not depend on where a slot sits.
    std::vector<size_t> live(n);
    std::iota(live.begin(), live.end(), size_t{0});
    std::vector<size_t> position(n);
    std::iota(position.begin(), position.end(), size_t{0});

    {
        ThreadPool pool(options.num_threads);
        copy_input(storage, work, options, pool, index);
        initial_row_minima(n, work, best, std::min(options.chunk_size, n), pool,
                           index);
    }

    Candidate current;
    for (size_t x = 0; x < n; ++x) {
        const Candidate candidate =
            make_candidate(best[x].distance, x, slot_node[x], best[x].partner_slot,
                           slot_node[best[x].partner_slot]);
        if (precedes(candidate, current)) {
            current = candidate;
        }
    }

    tree.children_left.reserve(n - 1);
    tree.children_right.reserve(n - 1);
    tree.distances.reserve(n - 1);
    tree.cluster_sizes.reserve(n - 1);

    // Captured by the step body and rewritten before each Run().
    size_t slot_of_low_node = 0;
    size_t slot_of_high_node = 0;
    size_t size_low = 0;
    size_t size_high = 0;
    size_t kept = 0;
    size_t kept_node = 0;
    size_t live_size = 0;
    const AgglomerativeLinkageMethod linkage = options.linkage;

    auto step = [&](size_t begin, size_t end, StepResult& out) {
        for (size_t p = begin; p < end; ++p) {
            const size_t x = live[p];
            const double distance_low = work[index(x, slot_of_low_node)];
            const double distance_high = work[index(x, slot_of_high_node)];
            const double merged = update_linkage_distance(
                linkage, distance_low, distance_high, size_low, size_high);
            work[index(x, kept)] = merged;

            const size_t cached = best[x].partner_slot;
            if (cached == slot_of_low_node || cached == slot_of_high_node) {
                // The cached minimum left with a child. Every update returns a
                // value in [min(d_low, d_high), max(d_low, d_high)], so the
                // exact merged distance is at least the old row minimum and
                // can never win the row outright; but average linkage rounds
                // its two products, the sum and the division separately and
                // can land one ULP below it, so there is no shortcut and the
                // row is rescanned unconditionally.
                double row_distance = INFINITE_DISTANCE;
                size_t row_partner = NO_SLOT;
                size_t row_partner_node = NO_NODE;
                for (size_t q = 0; q < live_size; ++q) {
                    const size_t y = live[q];
                    if (y == x) {
                        continue;
                    }
                    const double candidate_distance = work[index(x, y)];
                    const size_t candidate_node = slot_node[y];
                    if (row_precedes(candidate_distance, candidate_node, row_distance,
                                     row_partner_node)) {
                        row_distance = candidate_distance;
                        row_partner = y;
                        row_partner_node = candidate_node;
                    }
                }
                best[x] = RowBest{row_distance, row_partner};
                ++out.rescans;
                out.rescanned_slots += live_size - 1;
            } else if (row_precedes(merged, kept_node, best[x].distance,
                                    slot_node[cached])) {
                // The cached partner survived, so the cache is still valid for
                // every old candidate and only the merged cluster can beat it.
                best[x] = RowBest{merged, kept};
            }

            const size_t node = slot_node[x];
            if (row_precedes(merged, node, out.kept_distance, out.kept_partner_node)) {
                out.kept_distance = merged;
                out.kept_partner = x;
                out.kept_partner_node = node;
            }
            const Candidate candidate =
                make_candidate(best[x].distance, x, node, best[x].partner_slot,
                               slot_node[best[x].partner_slot]);
            if (precedes(candidate, out.global)) {
                out.global = candidate;
            }
        }
    };

    const size_t participants = row_cache_participants(
        options.num_threads, n, std::thread::hardware_concurrency());
    const size_t cutoff = options.serial_cutoff.value_or(ROW_CACHE_SERIAL_CUTOFF);

    std::unique_ptr<StepTeam> team;
    std::vector<StepResult> published(participants);
    // Claimed by every participant; kept off the line holding the step's
    // read-only size and unit.
    struct StepShape {
        size_t size = 0;
        size_t unit = 1;
    };
    alignas(64) std::atomic<size_t> next_position{0};
    alignas(64) StepShape shape;
    std::function<void(size_t)> step_body;
    // The widest step updates every slot but the kept one, n - 2 of them, so
    // below the cutoff no step can take the team path and every thread started
    // here would be joined without having updated a row. Starting them costs
    // more than the whole pass at small item counts.
    if (row_cache_team_is_useful(participants, n, cutoff)) {
        team = std::make_unique<StepTeam>(participants);
        step_body = [&](size_t participant) {
            StepResult local;
            while (!team->Failed()) {
                const size_t begin =
                    next_position.fetch_add(shape.unit, std::memory_order_relaxed);
                if (begin >= shape.size) {
                    break;
                }
                step(begin, std::min(begin + shape.unit, shape.size), local);
            }
            published[participant] = local;
        };
    }

    size_t next_node = n;
    while (live.size() > 1 && tree.distances.size() < options.target_merges) {
        // Snapshot both children before anything is mutated: the update of the
        // kept row weighs by the sizes they had, and the kept slot is about to
        // take a new id and the combined size.
        slot_of_low_node = current.left_slot;
        slot_of_high_node = current.right_slot;
        size_low = cluster_size[slot_of_low_node];
        size_high = cluster_size[slot_of_high_node];
        const size_t merged_size = size_low + size_high;
        kept = std::min(slot_of_low_node, slot_of_high_node);
        const size_t retired = std::max(slot_of_low_node, slot_of_high_node);

        tree.children_left.push_back(current.left);
        tree.children_right.push_back(current.right);
        tree.distances.push_back(current.distance);
        tree.cluster_sizes.push_back(merged_size);

        kept_node = next_node++;
        slot_node[kept] = kept_node;
        cluster_size[kept] = merged_size;

        const size_t retired_position = position[retired];
        live[retired_position] = live.back();
        position[live[retired_position]] = retired_position;
        live.pop_back();
        // The kept slot sits last so the step's range is a prefix that does
        // not contain it; its own row is reduced from the values the step
        // writes, which is the only place they all exist at once.
        const size_t kept_position = position[kept];
        const size_t last = live.size() - 1;
        std::swap(live[kept_position], live[last]);
        position[live[kept_position]] = kept_position;
        position[live[last]] = last;

        live_size = live.size();
        const size_t range = live_size - 1;
        StepResult combined;
        if (participants == 1 || !team || range < cutoff) {
            step(0, range, combined);
        } else {
            shape.size = range;
            shape.unit = work_unit(range, participants);
            next_position.store(0, std::memory_order_relaxed);
            team->Run(step_body);
            for (const StepResult& result : published) {
                if (row_precedes(result.kept_distance, result.kept_partner_node,
                                 combined.kept_distance, combined.kept_partner_node)) {
                    combined.kept_distance = result.kept_distance;
                    combined.kept_partner = result.kept_partner;
                    combined.kept_partner_node = result.kept_partner_node;
                }
                if (precedes(result.global, combined.global)) {
                    combined.global = result.global;
                }
                combined.rescans += result.rescans;
                combined.rescanned_slots += result.rescanned_slots;
            }
        }

        best[kept] = RowBest{combined.kept_distance, combined.kept_partner};
        current = combined.global;
        if (combined.kept_partner != NO_SLOT) {
            const Candidate kept_candidate = make_candidate(
                combined.kept_distance, kept, kept_node, combined.kept_partner,
                slot_node[combined.kept_partner]);
            if (precedes(kept_candidate, current)) {
                current = kept_candidate;
            }
        }
        if (options.stats != nullptr) {
            options.stats->rescans += combined.rescans;
            options.stats->rescanned_slots += combined.rescanned_slots;
            options.stats->updates += range;
        }
    }

    return tree;
}

}  // namespace OECluster::detail
