/**
 * @file DiversitySelection.cpp
 * @brief Farthest-first subset selection and the #Circles coverage measure.
 */

#include "oecluster/clustering/DiversitySelection.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "ChunkedComparisons.h"
#include "DistanceAccess.h"
#include "DiversityValidation.h"
#include "MaxMinKernel.h"

namespace OECluster {

namespace {

using detail::ChunkedComparisons;
using detail::validate_comparison_facts;
using detail::validate_size_and_chunk;

constexpr const char* SELECTION_NAME = "MaxMin selection";
constexpr const char* CIRCLES_NAME = "#Circles";

void validate_seed_mode(MaxMinSeed seed_mode) {
    switch (seed_mode) {
        case MaxMinSeed::Index:
        case MaxMinSeed::Medoid:
        case MaxMinSeed::Farthest:
            return;
    }
    throw std::invalid_argument("Unknown MaxMin selection seed mode");
}

void validate_circles_method(CirclesMethod method) {
    switch (method) {
        case CirclesMethod::MaxMin:
        case CirclesMethod::Sequential:
            return;
    }
    throw std::invalid_argument("Unknown #Circles method");
}

void validate_selection(const MaxMinOptions& options, size_t n) {
    const bool has_threshold = !std::isnan(options.threshold);
    if (options.count == 0 && !has_threshold) {
        throw std::invalid_argument(
            "MaxMin selection requires a count, a threshold, or both");
    }
    if (options.count > n) {
        throw std::invalid_argument(
            "MaxMin selection count must be at most the item count (" +
            std::to_string(n) + ")");
    }
    if (has_threshold &&
        (std::isinf(options.threshold) || options.threshold < 0.0)) {
        throw std::invalid_argument(
            "MaxMin selection threshold must be finite and non-negative");
    }
    if (options.seed_mode == MaxMinSeed::Index && options.seed >= n) {
        throw std::invalid_argument(
            "MaxMin selection seed is outside the item range");
    }
    if (options.initial.empty()) {
        return;
    }
    if (options.seed_mode != MaxMinSeed::Index || options.seed != 0) {
        throw std::invalid_argument(
            "MaxMin selection initial cannot be combined with a seed");
    }
    if (options.count != 0 && options.initial.size() > options.count) {
        throw std::invalid_argument(
            "MaxMin selection initial holds more entries than count");
    }
    std::vector<bool> seen(n, false);
    for (const size_t index : options.initial) {
        if (index >= n) {
            throw std::invalid_argument(
                "MaxMin selection initial index is outside the item range");
        }
        if (seen[index]) {
            throw std::invalid_argument(
                "MaxMin selection initial entries must be unique");
        }
        seen[index] = true;
    }
}

void validate_circles_threshold(double threshold) {
    if (!std::isfinite(threshold) || threshold < 0.0) {
        throw std::invalid_argument(
            "#Circles threshold must be finite and non-negative");
    }
}

MaxMinStop to_public_stop(detail::MaxMinKernelStop stop) {
    switch (stop) {
        case detail::MaxMinKernelStop::Count:
            return MaxMinStop::Count;
        case detail::MaxMinKernelStop::Threshold:
            return MaxMinStop::Threshold;
        case detail::MaxMinKernelStop::Exhausted:
            break;
    }
    return MaxMinStop::Exhausted;
}

MaxMinSelection to_selection(detail::MaxMinKernelResult&& result) {
    MaxMinSelection selection;
    selection.indices = std::move(result.indices);
    selection.pick_distances = std::move(result.pick_distances);
    selection.stop = to_public_stop(result.stop);
    return selection;
}

// Fold item 0's row with item 0 masked and take the argmax over the rest.
// Masking item 0 keeps the diagonal unread, which is what lets a comparison
// with a nonzero self-distance choose the same seed as the matrix.
template <typename RowProvider>
size_t farthest_seed(RowProvider& rows, size_t n) {
    if (n == 1) {
        return 0;
    }
    std::vector<bool> selected(n, false);
    selected[0] = true;
    std::vector<double> nearest(n, 0.0);
    rows.FoldRow(0, selected, nearest, true);

    size_t best = 1;
    for (size_t j = 2; j < n; ++j) {
        if (nearest[j] > nearest[best]) {
            best = j;
        }
    }
    return best;
}

// global_medoid validates neither finiteness nor overflow, and k-medoids keeps
// it that way. Several rows saturating to infinity would hand the medoid to
// the smallest such index rather than the true minimum, so this entry point
// refuses both. One serial pass over the condensed array does the validation
// and the sums together: the pair it names is always the first non-finite one
// in condensed order, and each row still receives its terms in ascending j
// order, so the sums are bitwise those medoid_row_sums would produce.
size_t medoid_seed(const double* data, size_t n) {
    std::vector<double> sums(n, 0.0);
    size_t k = 0;
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i + 1; j < n; ++j, ++k) {
            const double distance = data[k];
            if (!std::isfinite(distance)) {
                throw detail::non_finite_distance_error(i, j);
            }
            sums[i] += distance;
            sums[j] += distance;
        }
    }
    for (size_t i = 0; i < n; ++i) {
        if (!std::isfinite(sums[i])) {
            throw std::invalid_argument(
                "MaxMin selection medoid seed: the distance sum of item " +
                std::to_string(i) + " overflows");
        }
    }
    return detail::argmin_row_sum(sums);
}

detail::MaxMinKernelRequest selection_request(const MaxMinOptions& options,
                                              size_t seed) {
    detail::MaxMinKernelRequest request;
    request.count = options.count;
    request.threshold = options.threshold;
    request.initial = options.initial.empty() ? std::vector<size_t>{seed}
                                              : options.initial;
    return request;
}

CirclesResult circles_result(std::vector<size_t> members, double threshold,
                             CirclesMethod method) {
    CirclesResult result;
    result.count = members.size();
    result.members = std::move(members);
    result.threshold = threshold;
    result.method = method;
    return result;
}

// The reference pass's decision for one candidate, from the distances to every
// current member in member order. The whole buffer is always filled first, so
// the decision and any refusal depend only on the distances, never on which
// member happened to be compared first.
bool accept_candidate(const std::vector<double>& distances,
                      const std::vector<size_t>& members, size_t candidate,
                      double threshold) {
    bool accept = true;
    for (size_t position = 0; position < distances.size(); ++position) {
        if (!std::isfinite(distances[position])) {
            throw detail::non_finite_distance_error(members[position], candidate);
        }
        if (distances[position] <= threshold) {
            accept = false;
        }
    }
    return accept;
}

detail::MaxMinKernelRequest circles_request(double threshold) {
    detail::MaxMinKernelRequest request;
    request.threshold = threshold;
    request.initial = {0};
    return request;
}

// Row provider over a comparison. Always refuses non-finite reads: every
// entry point that builds one is new, so there is no legacy to preserve.
// Compare(min, max) reads the orientation pdist stores, so an asymmetric
// comparison gives the lazy and matrix paths the same numbers.
class ComparisonRowProvider {
public:
    explicit ComparisonRowProvider(ChunkedComparisons& work) : work_(work) {}

    void FoldRow(size_t p, const std::vector<bool>& selected,
                 std::vector<double>& nearest, bool first) {
        work_.Run(selected.size(),
                  [&](PairwiseComparison& local, size_t begin, size_t end) {
            for (size_t j = begin; j < end; ++j) {
                if (selected[j]) {
                    continue;
                }
                const double distance =
                    local.Compare(std::min(p, j), std::max(p, j));
                if (!std::isfinite(distance)) {
                    throw detail::non_finite_distance_error(p, j);
                }
                if (first || distance < nearest[j]) {
                    nearest[j] = distance;
                }
            }
        });
    }

private:
    ChunkedComparisons& work_;
};

}  // namespace

MaxMinSelection maxmin_select(const StorageBackend& storage,
                              const MaxMinOptions& options) {
    validate_seed_mode(options.seed_mode);
    detail::validate_complete_distance_storage(storage, SELECTION_NAME);
    const size_t n = storage.NumSamples();
    validate_size_and_chunk(n, options.chunk_size, SELECTION_NAME);
    validate_selection(options, n);

    const double* data = storage.Data();
    detail::MatrixRowProvider rows(data, n, /*refuse_non_finite=*/true);
    size_t seed = options.seed;
    if (options.initial.empty()) {
        if (options.seed_mode == MaxMinSeed::Farthest) {
            seed = farthest_seed(rows, n);
        } else if (options.seed_mode == MaxMinSeed::Medoid) {
            seed = medoid_seed(data, n);
        }
    }
    return to_selection(
        detail::maxmin_run(rows, n, selection_request(options, seed)));
}

MaxMinSelection maxmin_select(PairwiseComparison& comparison,
                              const MaxMinOptions& options) {
    validate_seed_mode(options.seed_mode);
    validate_comparison_facts(comparison, SELECTION_NAME);
    const size_t n = comparison.Size();
    validate_size_and_chunk(n, options.chunk_size, SELECTION_NAME);
    validate_selection(options, n);
    // Refused rather than computed: a medoid needs all O(N^2) comparisons,
    // which is exactly what the lazy path exists to avoid.
    if (options.seed_mode == MaxMinSeed::Medoid) {
        throw std::invalid_argument(
            "MaxMin selection seed mode Medoid requires a distance matrix");
    }

    ChunkedComparisons work(comparison, n, options.num_threads,
                            options.chunk_size);
    ComparisonRowProvider rows(work);
    size_t seed = options.seed;
    if (options.initial.empty() && options.seed_mode == MaxMinSeed::Farthest) {
        seed = farthest_seed(rows, n);
    }
    return to_selection(
        detail::maxmin_run(rows, n, selection_request(options, seed)));
}

CirclesResult circles(const StorageBackend& storage, double threshold,
                      const CirclesOptions& options) {
    validate_circles_method(options.method);
    detail::validate_complete_distance_storage(storage, CIRCLES_NAME);
    const size_t n = storage.NumSamples();
    validate_size_and_chunk(n, options.chunk_size, CIRCLES_NAME);
    validate_circles_threshold(threshold);

    const double* data = storage.Data();
    if (options.method == CirclesMethod::MaxMin) {
        detail::MatrixRowProvider rows(data, n, /*refuse_non_finite=*/true);
        return circles_result(
            detail::maxmin_run(rows, n, circles_request(threshold)).indices,
            threshold, options.method);
    }

    std::vector<size_t> members;
    std::vector<double> distances;
    for (size_t candidate = 0; candidate < n; ++candidate) {
        distances.resize(members.size());
        for (size_t position = 0; position < members.size(); ++position) {
            distances[position] =
                detail::dense_distance(data, n, members[position], candidate);
        }
        if (accept_candidate(distances, members, candidate, threshold)) {
            members.push_back(candidate);
        }
    }
    return circles_result(std::move(members), threshold, options.method);
}

CirclesResult circles(PairwiseComparison& comparison, double threshold,
                      const CirclesOptions& options) {
    validate_circles_method(options.method);
    validate_comparison_facts(comparison, CIRCLES_NAME);
    const size_t n = comparison.Size();
    validate_size_and_chunk(n, options.chunk_size, CIRCLES_NAME);
    validate_circles_threshold(threshold);

    ChunkedComparisons work(comparison, n, options.num_threads,
                            options.chunk_size);
    if (options.method == CirclesMethod::MaxMin) {
        ComparisonRowProvider rows(work);
        return circles_result(
            detail::maxmin_run(rows, n, circles_request(threshold)).indices,
            threshold, options.method);
    }

    std::vector<size_t> members;
    std::vector<double> distances;
    for (size_t candidate = 0; candidate < n; ++candidate) {
        distances.resize(members.size());
        // Members precede the candidate in input order, so member < candidate
        // is already the (min, max) orientation.
        work.Run(members.size(),
                 [&](PairwiseComparison& local, size_t begin, size_t end) {
            for (size_t position = begin; position < end; ++position) {
                distances[position] = local.Compare(members[position], candidate);
            }
        });
        if (accept_candidate(distances, members, candidate, threshold)) {
            members.push_back(candidate);
        }
    }
    return circles_result(std::move(members), threshold, options.method);
}

}  // namespace OECluster
