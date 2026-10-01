/**
 * @file LeidenEngine.cpp
 * @brief The serial Leiden optimizer: local moving, refinement, aggregation.
 */

#include "LeidenEngine.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numeric>
#include <utility>
#include <vector>

namespace OECluster::detail {

namespace {

constexpr uint32_t NO_COMMUNITY = std::numeric_limits<uint32_t>::max();

// Per-target weight sums that reset in O(1): a slot is live only while its
// stamp matches the current epoch, so no O(n) clear runs per node.
class WeightAccumulator {
public:
    explicit WeightAccumulator(size_t n) : weight_(n, 0.0), stamp_(n, 0) {}

    void Begin() {
        ++epoch_;
        touched_.clear();
    }

    void Add(uint32_t target, double weight) {
        if (stamp_[target] != epoch_) {
            stamp_[target] = epoch_;
            weight_[target] = 0.0;
            touched_.push_back(target);
        }
        weight_[target] += weight;
    }

    double Get(uint32_t target) const {
        return stamp_[target] == epoch_ ? weight_[target] : 0.0;
    }

    /// Targets in order of first Add since Begin.
    std::vector<uint32_t>& Touched() { return touched_; }

private:
    std::vector<double> weight_;
    std::vector<size_t> stamp_;
    std::vector<uint32_t> touched_;
    size_t epoch_ = 0;
};

}  // namespace

LeidenLevel make_base_level(WeightedGraph graph) {
    LeidenLevel level;
    level.num_nodes = graph.num_nodes;
    level.offsets = std::move(graph.offsets);
    level.neighbors = std::move(graph.neighbors);
    level.weights = std::move(graph.weights);
    level.self.assign(level.num_nodes, 0.0);
    level.strength.assign(level.num_nodes, 0.0);
    level.size.assign(level.num_nodes, 1);
    for (size_t v = 0; v < level.num_nodes; ++v) {
        for (size_t e = level.offsets[v]; e < level.offsets[v + 1]; ++e) {
            level.strength[v] += level.weights[e];
            if (level.neighbors[e] > v) {
                level.m += level.weights[e];
            }
        }
    }
    return level;
}

uint64_t LeidenRng::Below(uint64_t bound) {
    // 2^64 mod bound: rejecting outputs below it leaves a range whose length
    // is a multiple of bound, so the remainder is unbiased.
    const uint64_t threshold = (uint64_t{0} - bound) % bound;
    for (;;) {
        const uint64_t x = Next();
        if (x >= threshold) {
            return x % bound;
        }
    }
}

double LeidenRng::Uniform() {
    return static_cast<double>(Next() >> 11) * 0x1.0p-53;
}

void LeidenRng::Shuffle(std::vector<uint32_t>& values) {
    for (size_t i = values.size(); i > 1; --i) {
        const size_t j = static_cast<size_t>(Below(i));
        std::swap(values[i - 1], values[j]);
    }
}

double leiden_quality(const LeidenLevel& level,
                      const std::vector<uint32_t>& community,
                      const LeidenParams& params) {
    const size_t n = level.num_nodes;
    std::vector<double> internal(n, 0.0);
    std::vector<double> strength(n, 0.0);
    std::vector<double> size(n, 0.0);
    for (size_t v = 0; v < n; ++v) {
        const uint32_t c = community[v];
        internal[c] += level.self[v];
        strength[c] += level.strength[v];
        size[c] += static_cast<double>(level.size[v]);
        for (size_t e = level.offsets[v]; e < level.offsets[v + 1]; ++e) {
            const uint32_t u = level.neighbors[e];
            if (u > v && community[u] == c) {
                internal[c] += level.weights[e];
            }
        }
    }
    double quality = 0.0;
    if (params.objective == LeidenObjective::Modularity) {
        if (level.m == 0.0) {
            return 0.0;
        }
        for (size_t c = 0; c < n; ++c) {
            const double share = strength[c] / (2.0 * level.m);
            quality += internal[c] / level.m - params.resolution * share * share;
        }
        return quality;
    }
    for (size_t c = 0; c < n; ++c) {
        quality += internal[c] - params.resolution * size[c] * (size[c] - 1.0) / 2.0;
    }
    return quality;
}

double leiden_gain(double weight_to, double node_strength, double node_size,
                   double community_strength, double community_size,
                   double m, const LeidenParams& params) {
    if (params.objective == LeidenObjective::Modularity) {
        if (m == 0.0) {
            return weight_to;
        }
        return weight_to -
               params.resolution * node_strength * community_strength / (2.0 * m);
    }
    return weight_to - params.resolution * node_size * community_size;
}

bool leiden_well_connected(double weight_to_rest, double set_strength,
                           double set_size, double community_strength,
                           double community_size, double m,
                           const LeidenParams& params) {
    if (params.objective == LeidenObjective::Modularity) {
        if (m == 0.0) {
            return weight_to_rest >= 0.0;
        }
        return weight_to_rest >= params.resolution * set_strength *
                                     (community_strength - set_strength) /
                                     (2.0 * m);
    }
    return weight_to_rest >=
           params.resolution * set_size * (community_size - set_size);
}

size_t leiden_select_candidate(const std::vector<double>& gains, double theta,
                               double u) {
    double best = -std::numeric_limits<double>::infinity();
    for (const double gain : gains) {
        if (gain >= 0.0 && gain > best) {
            best = gain;
        }
    }
    if (best < 0.0) {
        return gains.size();
    }
    double total = 0.0;
    for (const double gain : gains) {
        if (gain >= 0.0) {
            total += std::exp((gain - best) / theta);
        }
    }
    const double target = u * total;
    double cumulative = 0.0;
    size_t last = gains.size();
    for (size_t i = 0; i < gains.size(); ++i) {
        if (gains[i] < 0.0) {
            continue;
        }
        cumulative += std::exp((gains[i] - best) / theta);
        last = i;
        if (cumulative > target) {
            return i;
        }
    }
    // Rounding can leave the cumulative sum at or below target for u near 1.
    return last;
}

void leiden_local_move(const LeidenLevel& level, std::vector<uint32_t>& community,
                       const LeidenParams& params, LeidenRng& rng) {
    const size_t n = level.num_nodes;
    if (n == 0 || level.m == 0.0) {
        return;
    }
    std::vector<double> total_strength(n, 0.0);
    std::vector<double> total_size(n, 0.0);
    std::vector<size_t> members(n, 0);
    for (size_t v = 0; v < n; ++v) {
        total_strength[community[v]] += level.strength[v];
        total_size[community[v]] += static_cast<double>(level.size[v]);
        ++members[community[v]];
    }
    // Pushed in descending order so the smallest empty id is on top.
    std::vector<uint32_t> empty;
    for (size_t id = n; id-- > 0;) {
        if (members[id] == 0) {
            empty.push_back(static_cast<uint32_t>(id));
        }
    }

    std::vector<uint32_t> queue(n);
    std::iota(queue.begin(), queue.end(), uint32_t{0});
    rng.Shuffle(queue);
    // A ring buffer of capacity n: in_queue keeps each node in it at most
    // once.
    std::vector<char> in_queue(n, 1);
    size_t head = 0;
    size_t count = n;
    WeightAccumulator weight_to(n);

    while (count > 0) {
        const uint32_t v = queue[head];
        head = (head + 1) % n;
        --count;
        in_queue[v] = 0;

        const uint32_t d = community[v];
        const double node_strength = level.strength[v];
        const double node_size = static_cast<double>(level.size[v]);
        total_strength[d] -= node_strength;
        total_size[d] -= node_size;
        --members[d];

        weight_to.Begin();
        for (size_t e = level.offsets[v]; e < level.offsets[v + 1]; ++e) {
            weight_to.Add(community[level.neighbors[e]], level.weights[e]);
        }
        auto gain_of = [&](uint32_t c) {
            return leiden_gain(weight_to.Get(c), node_strength, node_size,
                               total_strength[c], total_size[c], level.m, params);
        };

        // Non-empty candidates: the neighbors' communities, and d unless v
        // was its sole member. Equal gains go to the lowest id.
        uint32_t best = NO_COMMUNITY;
        double best_gain = -std::numeric_limits<double>::infinity();
        auto consider = [&](uint32_t c) {
            const double gain = gain_of(c);
            if (gain > best_gain || (gain == best_gain && c < best)) {
                best = c;
                best_gain = gain;
            }
        };
        for (const uint32_t c : weight_to.Touched()) {
            consider(c);
        }
        if (members[d] > 0) {
            consider(d);
        }
        // The empty candidate's gain is 0 under both objectives, and it
        // ranks after every non-empty candidate with an equal gain.
        const uint32_t empty_candidate =
            members[d] == 0 ? d : (empty.empty() ? NO_COMMUNITY : empty.back());
        if (empty_candidate != NO_COMMUNITY &&
            (best == NO_COMMUNITY || 0.0 > best_gain)) {
            best = empty_candidate;
            best_gain = 0.0;
        }
        const uint32_t target = best_gain > gain_of(d) ? best : d;

        if (target != d) {
            if (target == empty_candidate && target != d) {
                empty.pop_back();
            }
            if (members[d] == 0) {
                empty.push_back(d);
            }
        }
        community[v] = target;
        total_strength[target] += node_strength;
        total_size[target] += node_size;
        ++members[target];

        if (target != d) {
            for (size_t e = level.offsets[v]; e < level.offsets[v + 1]; ++e) {
                const uint32_t u = level.neighbors[e];
                if (community[u] != target && !in_queue[u]) {
                    in_queue[u] = 1;
                    queue[(head + count) % n] = u;
                    ++count;
                }
            }
        }
    }
}

}  // namespace OECluster::detail
