/**
 * @file PrimMST.cpp
 * @brief Dense Prim's algorithm on a persistent team of threads.
 */

#include "PrimMST.h"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <functional>
#include <limits>
#include <memory>
#include <numeric>
#include <sstream>
#include <thread>
#include <utility>

#include "DistanceAccess.h"
#include "StepTeam.h"
#include "ThresholdGraph.h"

namespace OECluster::detail {

namespace {

constexpr double INFINITE_REACH = std::numeric_limits<double>::infinity();
constexpr size_t NO_NODE = std::numeric_limits<size_t>::max();

// One participant's best candidate for a step. Padded so that participants
// publishing at the same time do not share a cache line.
struct alignas(64) BestCandidate {
    double reach = INFINITE_REACH;
    size_t node = NO_NODE;
    size_t position = NO_NODE;
};

// The total order the 5.19.0 scan applied implicitly: the smallest reach, and
// among equal reaches the lowest index, which is the first one it met.
bool precedes(double reach, size_t node, const BestCandidate& best) {
    return reach < best.reach || (reach == best.reach && node < best.node);
}

// The domain every read passes through. -0.0 == 0.0, so the canonicalization
// changes no comparison, only the sign bit of a zero that reaches the output.
double checked_distance(double distance, bool non_negative,
                        const std::string& caller, size_t a, size_t b) {
    if (!std::isfinite(distance)) {
        throw non_finite_distance_error(caller, std::min(a, b), std::max(a, b));
    }
    if (non_negative && distance < 0.0) {
        throw negative_distance_error(caller, std::min(a, b), std::max(a, b));
    }
    return distance == 0.0 ? 0.0 : distance;
}

// read(participant, current, candidate) returns a checked distance.
template <typename Read>
std::vector<HDBSCANMSTEdge> prim_kernel(size_t n, const PrimWeights& weights,
                                        size_t participants, size_t cutoff,
                                        const std::string& caller, Read& read) {
    std::vector<HDBSCANMSTEdge> mst;
    if (n < 2) {
        return mst;
    }
    mst.reserve(n - 1);

    const bool mutual = !weights.core.empty();
    const double* core = mutual ? weights.core.data() : nullptr;
    std::vector<double> reach(n, INFINITE_REACH);
    std::vector<size_t> source(n, 0);
    std::vector<size_t> remaining(n - 1);
    std::iota(remaining.begin(), remaining.end(), size_t{1});
    size_t current = 0;

    auto scan = [&](size_t participant, size_t begin, size_t end, BestCandidate& best) {
        for (size_t position = begin; position < end; ++position) {
            const size_t candidate = remaining[position];
            double& candidate_reach = reach[candidate];
            // The weight is at least max(core_current, core_candidate), and
            // candidate_reach >= core_candidate always holds once set, so in
            // either case the strict-< update below could not fire.
            const bool pruned = weights.prune &&
                                (core[current] >= candidate_reach ||
                                 core[candidate] == candidate_reach);
            if (!pruned) {
                const double distance = read(participant, current, candidate);
                double weight = distance;
                if (mutual) {
                    const double quotient = distance / weights.alpha;
                    if (!std::isfinite(quotient)) {
                        throw alpha_overflow_error(caller, weights.alpha);
                    }
                    weight = std::max({core[current], core[candidate], quotient});
                }
                if (weight < candidate_reach) {
                    candidate_reach = weight;
                    source[candidate] = current;
                }
            }
            if (precedes(candidate_reach, candidate, best)) {
                best = BestCandidate{candidate_reach, candidate, position};
            }
        }
    };

    std::unique_ptr<StepTeam> team;
    std::vector<BestCandidate> published(participants);
    // Claimed by every participant; kept off the line holding the step's
    // read-only size and unit.
    struct StepShape {
        size_t size = 0;
        size_t unit = 1;
    };
    alignas(64) std::atomic<size_t> next_position{0};
    alignas(64) StepShape shape;
    std::function<void(size_t)> step_body;
    // The largest candidate list a step ever sees is n - 1, so below the cutoff
    // no step can take the team path and every thread started here would be
    // joined without having scanned anything. Starting them costs more than the
    // whole pass at small item counts.
    if (participants > 1 && n - 1 >= cutoff) {
        team = std::make_unique<StepTeam>(participants);
        step_body = [&](size_t participant) {
            BestCandidate local;
            while (!team->Failed()) {
                const size_t begin =
                    next_position.fetch_add(shape.unit, std::memory_order_relaxed);
                if (begin >= shape.size) {
                    break;
                }
                scan(participant, begin, std::min(begin + shape.unit, shape.size), local);
            }
            published[participant] = local;
        };
    }

    for (size_t step = 0; step + 1 < n; ++step) {
        const size_t remaining_count = remaining.size();
        BestCandidate best;
        if (participants == 1 || remaining_count < cutoff) {
            scan(0, 0, remaining_count, best);
        } else {
            shape.size = remaining_count;
            shape.unit = work_unit(remaining_count, participants);
            next_position.store(0, std::memory_order_relaxed);
            team->Run(step_body);
            for (const BestCandidate& candidate : published) {
                if (precedes(candidate.reach, candidate.node, best)) {
                    best = candidate;
                }
            }
        }
        mst.push_back(HDBSCANMSTEdge{source[best.node], best.node, best.reach});
        // Order does not matter: the next step's choice depends only on the keys.
        remaining[best.position] = remaining.back();
        remaining.pop_back();
        current = best.node;
    }
    return mst;
}

}  // namespace

size_t prim_participants(size_t num_threads, size_t n, size_t hardware) {
    return resolve_participants(num_threads, n,
                                std::min(hardware, PRIM_DEFAULT_PARTICIPANTS));
}

std::invalid_argument alpha_overflow_error(const std::string& caller, double alpha) {
    std::ostringstream message;
    message << caller << " alpha=" << alpha
            << " is too small: a distance divided by alpha is not finite";
    return std::invalid_argument(message.str());
}

std::vector<HDBSCANMSTEdge> prim_mst(const StorageBackend& storage,
                                     const PrimWeights& weights,
                                     const PrimOptions& options) {
    const size_t n = storage.NumSamples();
    const double* data = storage.Data();
    const bool non_negative = !weights.core.empty();
    auto read = [&](size_t, size_t a, size_t b) {
        return checked_distance(dense_distance(data, n, a, b), non_negative,
                                options.caller, a, b);
    };
    const size_t participants = prim_participants(
        options.num_threads, n, std::thread::hardware_concurrency());
    return prim_kernel(n, weights, participants,
                       options.serial_cutoff.value_or(PRIM_SERIAL_CUTOFF_MATRIX),
                       options.caller, read);
}

std::vector<HDBSCANMSTEdge> prim_mst(PairwiseComparison& comparison,
                                     const PrimWeights& weights,
                                     const PrimOptions& options) {
    const size_t n = comparison.Size();
    const size_t participants = prim_participants(
        options.num_threads, n, std::thread::hardware_concurrency());
    // Built here, serially, before the team starts: Clone() is const, but
    // nothing promises that a const call is safe to make concurrently.
    std::vector<std::unique_ptr<PairwiseComparison>> clones;
    clones.reserve(participants);
    for (size_t participant = 0; participant < participants; ++participant) {
        clones.push_back(comparison.Clone());
    }
    const bool non_negative = !weights.core.empty();
    auto read = [&](size_t participant, size_t a, size_t b) {
        return checked_distance(
            clones[participant]->Compare(std::min(a, b), std::max(a, b)),
            non_negative, options.caller, a, b);
    };
    return prim_kernel(n, weights, participants,
                       options.serial_cutoff.value_or(2 * participants),
                       options.caller, read);
}

}  // namespace OECluster::detail
