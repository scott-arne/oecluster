/**
 * @file isim_test_support.h
 * @brief Brute-force iSIM oracle and fingerprint builders for the iSIM tests.
 *
 * Every quantity is recomputed pair by pair with bit-at-a-time counting, so
 * the oracle shares no code with the kernels under test. Integer pair sums
 * equal the kernels' count-vector identities exactly; the doubles they turn
 * into round the same way, which is why agreement is required at rel 1e-14.
 */

#ifndef OECLUSTER_TESTS_CPP_ISIM_TEST_SUPPORT_H
#define OECLUSTER_TESTS_CPP_ISIM_TEST_SUPPORT_H

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <initializer_list>
#include <limits>
#include <string>
#include <vector>

#include "oefp/batch.h"
#include "oefp/fingerprint.h"

namespace isim_test {

inline OEFP::OEFP make_fp(const size_t size_bits, std::initializer_list<size_t> on_bits) {
    OEFP::FingerprintSpec spec;
    spec.size_bits = size_bits;
    spec.value_type = OEFP::FingerprintValueType::Binary;
    spec.source_name = "test";
    OEFP::OEFP fp(spec);
    for (const size_t bit : on_bits) {
        fp.SetBit(bit);
    }
    return fp;
}

inline OEFP::OEFPBatch make_batch(const std::vector<OEFP::OEFP>& fps) {
    return OEFP::OEFPBatch::FromFingerprints(fps);
}

// xorshift64 at roughly 12.5% density. No bit is forced on, so a row can be
// all-zero; callers that need the zero-union paths build them explicitly.
inline OEFP::OEFPBatch make_random_batch(const size_t rows, const size_t bits,
                                         uint64_t seed) {
    std::vector<OEFP::OEFP> fps;
    fps.reserve(rows);
    auto next = [&seed]() {
        seed ^= seed << 13; seed ^= seed >> 7; seed ^= seed << 17;
        return seed;
    };
    for (size_t i = 0; i < rows; ++i) {
        OEFP::FingerprintSpec spec;
        spec.size_bits = bits;
        spec.value_type = OEFP::FingerprintValueType::Binary;
        spec.source_name = "test";
        OEFP::OEFP fp(spec);
        for (size_t b = 0; b < bits; ++b) {
            if ((next() & 7u) == 0u) fp.SetBit(b);
        }
        fps.push_back(fp);
    }
    return OEFP::OEFPBatch::FromFingerprints(fps);
}

inline bool bit_set(const OEFP::OEFPBatch& batch, size_t row, size_t bit) {
    return ((batch.RowWords(row)[bit / 64u] >> (bit % 64u)) & 1u) != 0u;
}

inline uint64_t popcount(const OEFP::OEFPBatch& batch, size_t row) {
    uint64_t count = 0;
    for (size_t b = 0; b < batch.SizeBits(); ++b) count += bit_set(batch, row, b) ? 1u : 0u;
    return count;
}

inline uint64_t intersection(const OEFP::OEFPBatch& batch, size_t i, size_t j) {
    uint64_t count = 0;
    for (size_t b = 0; b < batch.SizeBits(); ++b) {
        count += (bit_set(batch, i, b) && bit_set(batch, j, b)) ? 1u : 0u;
    }
    return count;
}

inline uint64_t union_count(const OEFP::OEFPBatch& batch, size_t i, size_t j) {
    return popcount(batch, i) + popcount(batch, j) - intersection(batch, i, j);
}

inline double ratio(uint64_t intersections, uint64_t unions) {
    return unions == 0u ? 1.0
                        : static_cast<double>(intersections) / static_cast<double>(unions);
}

inline double distance(const OEFP::OEFPBatch& batch, size_t i, size_t j) {
    return 1.0 - ratio(intersection(batch, i, j), union_count(batch, i, j));
}

inline bool same_double(double a, double b) {
    if (std::isnan(a) || std::isnan(b)) return std::isnan(a) && std::isnan(b);
    return a == b && std::signbit(a) == std::signbit(b);
}

inline void expect_rel(double actual, double expected, const std::string& what) {
    if (std::isnan(expected) || std::isinf(expected)) {
        EXPECT_TRUE(same_double(actual, expected)) << what << ": " << actual << " vs " << expected;
        return;
    }
    EXPECT_LE(std::fabs(actual - expected), 1e-14 * std::fabs(expected))
        << what << ": " << actual << " vs " << expected;
}

inline double median(std::vector<double> values) {
    if (values.empty()) return 0.0;
    std::sort(values.begin(), values.end());
    const size_t mid = values.size() / 2;
    return values.size() % 2 == 1 ? values[mid] : (values[mid - 1] + values[mid]) / 2.0;
}

/// Everything isim_report computes beyond the exact profile, by brute force.
struct OracleReport {
    double isim_intra_distance = std::numeric_limits<double>::quiet_NaN();
    double isim_inter_distance = std::numeric_limits<double>::quiet_NaN();
    double median_radius = std::numeric_limits<double>::quiet_NaN();
    double median_medoid_member_distance = std::numeric_limits<double>::quiet_NaN();
    double calinski_harabasz_medoid = std::numeric_limits<double>::quiet_NaN();
    double isim_silhouette = std::numeric_limits<double>::quiet_NaN();
    double davies_bouldin_medoid = std::numeric_limits<double>::quiet_NaN();
    double dunn_medoid_separation_medoid_spread = std::numeric_limits<double>::quiet_NaN();
    std::vector<double> coverage_at;
    std::vector<double> noise_coverage_at;
    std::vector<size_t> medoids;
    std::vector<double> record_intra;
    std::vector<double> separation;
    std::vector<double> radius;
    std::vector<double> mean_medoid_distance;
    std::vector<double> record_silhouette;
    std::vector<int> nearest_cluster;
    std::vector<double> nearest_similarity;
};

// Best-by-exact-fraction with ties to the lowest index; a zero denominator
// scores 1. Small test sums keep the cross products inside uint64_t.
inline size_t oracle_medoid(const OEFP::OEFPBatch& batch, const std::vector<size_t>& set) {
    size_t best = set.front();
    uint64_t best_num = 0, best_den = 0;
    bool first = true;
    for (const size_t i : set) {
        uint64_t num = 0, den = 0;
        for (const size_t j : set) {
            if (j == i) continue;
            num += intersection(batch, i, j);
            den += union_count(batch, i, j);
        }
        if (den == 0u) { num = 1u; den = 1u; }
        const bool better = first || num * best_den > best_num * den ||
                            (num * best_den == best_num * den && i < best);
        if (better) { best = i; best_num = num; best_den = den; }
        first = false;
    }
    return best;
}

inline OracleReport oracle_report(const OEFP::OEFPBatch& batch,
                                  const std::vector<int>& labels,
                                  const std::vector<double>& thresholds) {
    OracleReport out;
    int max_label = -1;
    for (const int label : labels) max_label = std::max(max_label, label);
    const size_t k_count = static_cast<size_t>(max_label + 1);
    std::vector<std::vector<size_t>> clusters(k_count);
    std::vector<size_t> clustered;
    size_t num_noise = 0;
    for (size_t i = 0; i < labels.size(); ++i) {
        if (labels[i] >= 0) {
            clusters[static_cast<size_t>(labels[i])].push_back(i);
            clustered.push_back(i);
        } else {
            ++num_noise;
        }
    }
    if (k_count == 0) return out;

    uint64_t intra_i = 0, intra_u = 0;
    bool has_pair = false;
    out.record_intra.assign(k_count, std::numeric_limits<double>::quiet_NaN());
    for (size_t k = 0; k < k_count; ++k) {
        uint64_t ki = 0, ku = 0;
        const auto& c = clusters[k];
        for (size_t a = 0; a < c.size(); ++a)
            for (size_t b = a + 1; b < c.size(); ++b) {
                ki += intersection(batch, c[a], c[b]);
                ku += union_count(batch, c[a], c[b]);
            }
        if (c.size() >= 2) { has_pair = true; out.record_intra[k] = 1.0 - ratio(ki, ku); }
        intra_i += ki; intra_u += ku;
        out.medoids.push_back(oracle_medoid(batch, c));
    }
    if (has_pair) out.isim_intra_distance = 1.0 - ratio(intra_i, intra_u);

    std::vector<double> scatter(k_count), square(k_count);
    for (size_t k = 0; k < k_count; ++k) {
        double r = 0.0, total = 0.0, s = 0.0, sq = 0.0;
        for (const size_t m : clusters[k]) {
            const double d = distance(batch, out.medoids[k], m);
            r = std::max(r, d);
            if (m != out.medoids[k]) total += d;
            s += d; sq += d * d;
        }
        out.radius.push_back(r);
        out.mean_medoid_distance.push_back(
            clusters[k].size() < 2 ? 0.0 : total / static_cast<double>(clusters[k].size() - 1));
        scatter[k] = s / static_cast<double>(clusters[k].size());
        square[k] = sq;
    }
    out.median_radius = median(out.radius);
    out.median_medoid_member_distance = median(out.mean_medoid_distance);

    out.separation.assign(k_count, std::numeric_limits<double>::quiet_NaN());
    out.record_silhouette.assign(k_count, std::numeric_limits<double>::quiet_NaN());
    out.nearest_cluster.assign(k_count, -1);
    out.nearest_similarity.assign(k_count, std::numeric_limits<double>::quiet_NaN());
    if (k_count >= 2) {
        uint64_t xi = 0, xu = 0;
        for (size_t k = 0; k < k_count; ++k)
            for (size_t l = k + 1; l < k_count; ++l)
                for (const size_t i : clusters[k])
                    for (const size_t j : clusters[l]) {
                        xi += intersection(batch, i, j);
                        xu += union_count(batch, i, j);
                    }
        out.isim_inter_distance = 1.0 - ratio(xi, xu);

        for (size_t k = 0; k < k_count; ++k) {
            uint64_t si = 0, su = 0;
            for (const size_t i : clusters[k])
                for (const size_t j : clustered)
                    if (labels[j] != static_cast<int>(k)) {
                        si += intersection(batch, i, j);
                        su += union_count(batch, i, j);
                    }
            out.separation[k] = 1.0 - ratio(si, su);
        }

        const size_t global = oracle_medoid(batch, clustered);
        double between = 0.0, within = 0.0;
        for (size_t k = 0; k < k_count; ++k) {
            const double d = distance(batch, out.medoids[k], global);
            between += static_cast<double>(clusters[k].size()) * d * d;
            within += square[k];
        }
        if (clustered.size() > k_count) {
            const double num = between / static_cast<double>(k_count - 1);
            const double den = within / static_cast<double>(clustered.size() - k_count);
            out.calinski_harabasz_medoid =
                den == 0.0 ? std::numeric_limits<double>::quiet_NaN() : num / den;
        }

        // Centroid stage.
        double sil_total = 0.0;
        for (size_t k = 0; k < k_count; ++k) {
            double sil_sum = 0.0;
            for (const size_t i : clusters[k]) {
                double s = 0.0;
                if (clusters[k].size() >= 2) {
                    uint64_t oi = 0, ou = 0;
                    for (const size_t j : clusters[k])
                        if (j != i) { oi += intersection(batch, i, j); ou += union_count(batch, i, j); }
                    const double a = 1.0 - ratio(oi, ou);
                    double b = std::numeric_limits<double>::infinity();
                    for (size_t l = 0; l < k_count; ++l) {
                        if (l == k) continue;
                        uint64_t li = 0, lu = 0;
                        for (const size_t j : clusters[l]) { li += intersection(batch, i, j); lu += union_count(batch, i, j); }
                        b = std::min(b, 1.0 - ratio(li, lu));
                    }
                    const double m = std::max(a, b);
                    s = m == 0.0 ? 0.0 : (b - a) / m;
                }
                sil_sum += s;
            }
            sil_total += sil_sum;
            out.record_silhouette[k] = sil_sum / static_cast<double>(clusters[k].size());

            double best = -std::numeric_limits<double>::infinity();
            for (size_t l = 0; l < k_count; ++l) {
                if (l == k) continue;
                uint64_t ci = 0, cu = 0;
                for (const size_t i : clusters[k])
                    for (const size_t j : clusters[l]) { ci += intersection(batch, i, j); cu += union_count(batch, i, j); }
                const double sim = ratio(ci, cu);
                if (sim > best) { best = sim; out.nearest_cluster[k] = static_cast<int>(l); }
            }
            out.nearest_similarity[k] = best;
        }
        out.isim_silhouette = sil_total / static_cast<double>(clustered.size());

        double db_total = 0.0, min_sep = std::numeric_limits<double>::infinity(), max_spread = 0.0;
        for (size_t a = 0; a < k_count; ++a) {
            double worst = 0.0;
            for (size_t b = 0; b < k_count; ++b) {
                if (a == b) continue;
                const double sep = distance(batch, out.medoids[a], out.medoids[b]);
                min_sep = std::min(min_sep, sep);
                worst = std::max(worst, sep == 0.0 ? std::numeric_limits<double>::infinity()
                                                   : (scatter[a] + scatter[b]) / sep);
            }
            db_total += worst;
            max_spread = std::max(max_spread, 2.0 * scatter[a]);
        }
        out.davies_bouldin_medoid = db_total / static_cast<double>(k_count);
        out.dunn_medoid_separation_medoid_spread =
            max_spread == 0.0 ? std::numeric_limits<double>::quiet_NaN() : min_sep / max_spread;
    }

    if (!labels.empty() && !thresholds.empty()) {
        out.coverage_at.assign(thresholds.size(), 0.0);
        out.noise_coverage_at.assign(thresholds.size(), std::numeric_limits<double>::quiet_NaN());
        for (size_t t = 0; t < thresholds.size(); ++t) {
            size_t covered = 0, covered_noise = 0;
            for (size_t i = 0; i < labels.size(); ++i) {
                double nearest = std::numeric_limits<double>::infinity();
                for (const size_t m : out.medoids) nearest = std::min(nearest, distance(batch, i, m));
                if (nearest <= thresholds[t]) { ++covered; if (labels[i] < 0) ++covered_noise; }
            }
            out.coverage_at[t] = static_cast<double>(covered) / static_cast<double>(labels.size());
            if (num_noise > 0)
                out.noise_coverage_at[t] = static_cast<double>(covered_noise) / static_cast<double>(num_noise);
        }
    }
    return out;
}

}  // namespace isim_test

#endif  // OECLUSTER_TESTS_CPP_ISIM_TEST_SUPPORT_H
