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

// Communities renumbered by first appearance in node order, so numbering by
// smallest member.
std::vector<uint32_t> canonical(const std::vector<uint32_t>& community) {
    std::vector<uint32_t> label_of(community.size(), NO_COMMUNITY);
    std::vector<uint32_t> labels(community.size());
    uint32_t next = 0;
    for (size_t v = 0; v < community.size(); ++v) {
        uint32_t& label = label_of[community[v]];
        if (label == NO_COMMUNITY) {
            label = next++;
        }
        labels[v] = label;
    }
    return labels;
}

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
    // With m == 0 every pass starts from singletons (aggregation preserves
    // m), and a singleton's only empty candidate is its own id, so no node
    // can move under either objective.
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
    // Pushed in descending order so the smallest empty id is on top. Ids
    // emptied during the phase are then pushed LIFO, as the spec requires,
    // not re-sorted.
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
        // Exact zeros for an emptied community: subtraction can leave a
        // residue that would give d a nonzero gain as the empty candidate.
        if (--members[d] == 0) {
            total_strength[d] = 0.0;
            total_size[d] = 0.0;
        }

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
            if (target == empty_candidate) {
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

std::vector<uint32_t> leiden_refine(const LeidenLevel& level,
                                    const std::vector<uint32_t>& community,
                                    const LeidenParams& params, LeidenRng& rng) {
    const size_t n = level.num_nodes;
    std::vector<uint32_t> refined(n);
    std::iota(refined.begin(), refined.end(), uint32_t{0});
    if (n == 0 || level.m == 0.0) {
        return refined;
    }

    // Members of each community in ascending node order.
    std::vector<size_t> start(n + 1, 0);
    for (size_t v = 0; v < n; ++v) {
        ++start[community[v] + 1];
    }
    for (size_t c = 0; c < n; ++c) {
        start[c + 1] += start[c];
    }
    std::vector<uint32_t> members(n);
    {
        std::vector<size_t> cursor(start.begin(), start.end() - 1);
        for (size_t v = 0; v < n; ++v) {
            members[cursor[community[v]]++] = static_cast<uint32_t>(v);
        }
    }

    std::vector<double> community_strength(n, 0.0);
    std::vector<double> community_size(n, 0.0);
    // Refined-community state, indexed by refined id.
    std::vector<double> set_strength(n);
    std::vector<double> set_size(n);
    std::vector<size_t> set_count(n, 1);
    std::vector<double> set_external(n, 0.0);
    for (size_t v = 0; v < n; ++v) {
        community_strength[community[v]] += level.strength[v];
        community_size[community[v]] += static_cast<double>(level.size[v]);
        set_strength[v] = level.strength[v];
        set_size[v] = static_cast<double>(level.size[v]);
        for (size_t e = level.offsets[v]; e < level.offsets[v + 1]; ++e) {
            if (community[level.neighbors[e]] == community[v]) {
                set_external[v] += level.weights[e];
            }
        }
    }

    WeightAccumulator weight_to(n);
    std::vector<uint32_t> order;
    std::vector<uint32_t> candidates;
    std::vector<double> gains;
    for (size_t c = 0; c < n; ++c) {
        if (start[c] == start[c + 1]) {
            continue;
        }
        order.assign(members.begin() + start[c], members.begin() + start[c + 1]);
        rng.Shuffle(order);
        const double c_strength = community_strength[c];
        const double c_size = community_size[c];
        for (const uint32_t v : order) {
            const uint32_t own = refined[v];
            if (set_count[own] != 1) {
                continue;
            }
            if (!leiden_well_connected(set_external[own], set_strength[own],
                                       set_size[own], c_strength, c_size,
                                       level.m, params)) {
                continue;
            }
            weight_to.Begin();
            for (size_t e = level.offsets[v]; e < level.offsets[v + 1]; ++e) {
                const uint32_t u = level.neighbors[e];
                if (community[u] == c) {
                    weight_to.Add(refined[u], level.weights[e]);
                }
            }
            std::vector<uint32_t>& touched = weight_to.Touched();
            std::sort(touched.begin(), touched.end());
            candidates.clear();
            gains.clear();
            for (const uint32_t t : touched) {
                if (!leiden_well_connected(set_external[t], set_strength[t],
                                           set_size[t], c_strength, c_size,
                                           level.m, params)) {
                    continue;
                }
                candidates.push_back(t);
                gains.push_back(leiden_gain(weight_to.Get(t), set_strength[own],
                                            set_size[own], set_strength[t],
                                            set_size[t], level.m, params));
            }
            if (candidates.empty()) {
                continue;
            }
            const size_t chosen =
                leiden_select_candidate(gains, params.theta, rng.Uniform());
            if (chosen == candidates.size()) {
                continue;
            }
            const uint32_t t = candidates[chosen];
            // The pair's mutual weight leaves both external sums.
            set_external[t] += set_external[own] - 2.0 * weight_to.Get(t);
            set_strength[t] += set_strength[own];
            set_size[t] += set_size[own];
            set_count[t] += 1;
            set_count[own] = 0;
            refined[v] = t;
        }
    }
    return refined;
}

LeidenAggregate leiden_aggregate(const LeidenLevel& level,
                                 const std::vector<uint32_t>& refined,
                                 const std::vector<uint32_t>& community) {
    const size_t n = level.num_nodes;
    LeidenAggregate result;
    result.node_of.resize(n);
    // Ascending nodes meet each refined community first at its smallest
    // member, which numbers the new nodes by smallest member.
    std::vector<uint32_t> id_of(n, NO_COMMUNITY);
    std::vector<uint32_t> first_member;
    for (size_t v = 0; v < n; ++v) {
        uint32_t& id = id_of[refined[v]];
        if (id == NO_COMMUNITY) {
            id = static_cast<uint32_t>(first_member.size());
            first_member.push_back(static_cast<uint32_t>(v));
        }
        result.node_of[v] = id;
    }
    const size_t count = first_member.size();

    std::vector<uint32_t> initial_of(n, NO_COMMUNITY);
    result.community.resize(count);
    uint32_t next_community = 0;
    for (size_t a = 0; a < count; ++a) {
        uint32_t& initial = initial_of[community[first_member[a]]];
        if (initial == NO_COMMUNITY) {
            initial = next_community++;
        }
        result.community[a] = initial;
    }

    LeidenLevel& next = result.level;
    next.num_nodes = count;
    next.self.assign(count, 0.0);
    next.strength.assign(count, 0.0);
    next.size.assign(count, 0);
    next.m = level.m;

    std::vector<size_t> start(count + 1, 0);
    for (size_t v = 0; v < n; ++v) {
        const uint32_t a = result.node_of[v];
        ++start[a + 1];
        next.self[a] += level.self[v];
        next.strength[a] += level.strength[v];
        next.size[a] += level.size[v];
    }
    for (size_t a = 0; a < count; ++a) {
        start[a + 1] += start[a];
    }
    std::vector<uint32_t> members(n);
    {
        std::vector<size_t> cursor(start.begin(), start.end() - 1);
        for (size_t v = 0; v < n; ++v) {
            members[cursor[result.node_of[v]]++] = static_cast<uint32_t>(v);
        }
    }

    // Each cross pair is summed from its smaller node, once to count the
    // rows and again, in the same order, to fill them, so both directions
    // hold the same double and no O(E) edge table is ever alive beside the
    // two levels.
    WeightAccumulator weight_to(count);
    const auto sum_cross = [&](size_t a, bool add_self) {
        weight_to.Begin();
        for (size_t p = start[a]; p < start[a + 1]; ++p) {
            const uint32_t v = members[p];
            for (size_t e = level.offsets[v]; e < level.offsets[v + 1]; ++e) {
                const uint32_t u = level.neighbors[e];
                const uint32_t b = result.node_of[u];
                if (b == a) {
                    if (add_self && u > v) {
                        next.self[a] += level.weights[e];
                    }
                } else if (b > a) {
                    weight_to.Add(b, level.weights[e]);
                }
            }
        }
    };

    next.offsets.assign(count + 1, 0);
    for (size_t a = 0; a < count; ++a) {
        sum_cross(a, true);
        next.offsets[a + 1] += weight_to.Touched().size();
        for (const uint32_t b : weight_to.Touched()) {
            ++next.offsets[b + 1];
        }
    }
    for (size_t a = 0; a < count; ++a) {
        next.offsets[a + 1] += next.offsets[a];
    }
    next.neighbors.resize(next.offsets[count]);
    next.weights.resize(next.offsets[count]);
    std::vector<size_t> cursor(next.offsets.begin(), next.offsets.end() - 1);
    for (size_t a = 0; a < count; ++a) {
        sum_cross(a, false);
        for (const uint32_t b : weight_to.Touched()) {
            const double weight = weight_to.Get(b);
            next.neighbors[cursor[a]] = b;
            next.weights[cursor[a]++] = weight;
            next.neighbors[cursor[b]] = static_cast<uint32_t>(a);
            next.weights[cursor[b]++] = weight;
        }
    }
    std::vector<std::pair<uint32_t, double>> row;
    for (size_t a = 0; a < count; ++a) {
        row.clear();
        for (size_t e = next.offsets[a]; e < next.offsets[a + 1]; ++e) {
            row.emplace_back(next.neighbors[e], next.weights[e]);
        }
        std::sort(row.begin(), row.end(),
                  [](const auto& x, const auto& y) { return x.first < y.first; });
        for (size_t e = next.offsets[a]; e < next.offsets[a + 1]; ++e) {
            next.neighbors[e] = row[e - next.offsets[a]].first;
            next.weights[e] = row[e - next.offsets[a]].second;
        }
    }
    return result;
}

std::vector<uint32_t> leiden_pass(const LeidenLevel& base,
                                  const std::vector<uint32_t>& partition,
                                  const LeidenParams& params, LeidenRng& rng) {
    const size_t n = base.num_nodes;
    std::vector<uint32_t> community = partition;
    std::vector<uint32_t> node_of_base(n);
    std::iota(node_of_base.begin(), node_of_base.end(), uint32_t{0});
    // Aggregates only shrink; the base level is read in place, not copied.
    const LeidenLevel* current = &base;
    LeidenLevel owned;
    for (;;) {
        leiden_local_move(*current, community, params, rng);
        const std::vector<uint32_t> refined =
            leiden_refine(*current, community, params, rng);
        LeidenAggregate aggregate = leiden_aggregate(*current, refined, community);
        if (aggregate.level.num_nodes == current->num_nodes) {
            break;
        }
        for (size_t i = 0; i < n; ++i) {
            node_of_base[i] = aggregate.node_of[node_of_base[i]];
        }
        community = std::move(aggregate.community);
        owned = std::move(aggregate.level);
        current = &owned;
    }
    std::vector<uint32_t> result(n);
    for (size_t i = 0; i < n; ++i) {
        result[i] = community[node_of_base[i]];
    }
    return result;
}

LeidenRun run_leiden(WeightedGraph graph, const LeidenParams& params,
                     int64_t n_iterations, uint64_t seed) {
    const LeidenLevel base = make_base_level(std::move(graph));
    LeidenRng rng(seed);
    std::vector<uint32_t> labels(base.num_nodes);
    std::iota(labels.begin(), labels.end(), uint32_t{0});
    LeidenRun run;
    if (n_iterations >= 0) {
        for (int64_t pass = 0; pass < n_iterations; ++pass) {
            labels = canonical(leiden_pass(base, labels, params, rng));
            ++run.iterations;
        }
    } else {
        for (;;) {
            std::vector<uint32_t> next =
                canonical(leiden_pass(base, labels, params, rng));
            ++run.iterations;
            if (next == labels) {
                break;
            }
            labels = std::move(next);
        }
    }
    run.quality = leiden_quality(base, labels, params);
    run.labels.assign(labels.begin(), labels.end());
    return run;
}

}  // namespace OECluster::detail
