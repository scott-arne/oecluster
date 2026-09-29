/**
 * @file SetDiversity.cpp
 * @brief The Vendi score and log-determinant diversity of a set.
 */

#include "oecluster/clustering/SetDiversity.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "ChunkedComparisons.h"
#include "DistanceAccess.h"
#include "DiversityValidation.h"
#include "SymmetricEigen.h"

namespace OECluster {

namespace {

constexpr const char* VENDI_NAME = "Vendi score";
constexpr const char* LOGDET_NAME = "Log-determinant diversity";

void validate_kernel(DiversityKernel kernel, double bandwidth,
                     const std::string& caller) {
    switch (kernel) {
        case DiversityKernel::Complement:
            if (!std::isnan(bandwidth)) {
                throw std::invalid_argument(
                    caller + " bandwidth applies only to the Laplacian kernel");
            }
            return;
        case DiversityKernel::Laplacian:
            if (!std::isfinite(bandwidth) || bandwidth <= 0.0) {
                throw std::invalid_argument(
                    caller +
                    " Laplacian kernel requires a finite, positive bandwidth");
            }
            return;
    }
    throw std::invalid_argument(caller +
                                " kernel is not a known DiversityKernel");
}

void validate_max_exact(size_t max_exact, const std::string& caller) {
    if (max_exact == 0) {
        throw std::invalid_argument(caller + " max_exact must be at least one");
    }
}

void validate_vendi_options(const VendiOptions& options) {
    if (options.order != 1 && options.order != 2) {
        throw std::invalid_argument("Vendi score order must be 1 or 2");
    }
    validate_kernel(options.kernel, options.bandwidth, VENDI_NAME);
    validate_max_exact(options.max_exact, VENDI_NAME);
    detail::validate_chunk_size(options.chunk_size, VENDI_NAME);
}

void validate_logdet_options(const LogDetOptions& options) {
    validate_kernel(options.kernel, options.bandwidth, LOGDET_NAME);
    if (!std::isfinite(options.ridge) || options.ridge < 0.0) {
        throw std::invalid_argument(
            "Log-determinant diversity ridge must be finite and non-negative");
    }
    validate_max_exact(options.max_exact, LOGDET_NAME);
    detail::validate_chunk_size(options.chunk_size, LOGDET_NAME);
}

// The spectrum costs 8n^2 bytes and O(n^3) time, so the default ceiling keeps
// an accidental whole-library call from exhausting either. Order 2 is the way
// past it because it needs no spectrum.
void validate_ceiling(size_t n, size_t max_exact, const std::string& caller) {
    if (n > max_exact) {
        throw std::invalid_argument(
            caller + " computes an exact spectrum of at most max_exact (" +
            std::to_string(max_exact) + ") items, but the input has " +
            std::to_string(n) +
            "; raise max_exact (memory grows as 8n^2 bytes and time as n^3) "
            "or use the Vendi score with order=2, which needs no spectrum");
    }
}

std::string format_distance(double distance) {
    std::ostringstream out;
    out << std::setprecision(std::numeric_limits<double>::max_digits10)
        << distance;
    return out.str();
}

std::string pair_name(size_t i, size_t j) {
    return "d(" + std::to_string(i) + ", " + std::to_string(j) + ")";
}

// Maps a distance to a kernel entry, refusing a distance the kernel cannot
// use. Neither the matrix gate nor the comparison facts promise d in [0, 1]
// or d >= 0, so the check happens on every read. Stateless, so workers share
// one instance.
class KernelTransform {
public:
    KernelTransform(DiversityKernel kernel, double bandwidth, const char* caller)
        : kernel_(kernel), bandwidth_(bandwidth), caller_(caller) {}

    double operator()(size_t i, size_t j, double distance) const {
        if (!std::isfinite(distance)) {
            throw std::invalid_argument(
                std::string(caller_) + " read a non-finite distance between items " +
                std::to_string(i) + " and " + std::to_string(j));
        }
        if (kernel_ == DiversityKernel::Complement) {
            if (distance < 0.0 || distance > 1.0) {
                throw std::invalid_argument(
                    std::string(caller_) +
                    " complement kernel requires distances in [0, 1], but " +
                    pair_name(i, j) + " = " + format_distance(distance) +
                    "; use the Laplacian kernel (kernel=\"laplacian\") for "
                    "other distances");
            }
            return 1.0 - distance;
        }
        // A negative distance would push the entry above 1, or overflow it.
        if (distance < 0.0) {
            throw std::invalid_argument(
                std::string(caller_) +
                " Laplacian kernel requires non-negative distances, but " +
                pair_name(i, j) + " = " + format_distance(distance));
        }
        return std::exp(-distance / bandwidth_);
    }

private:
    DiversityKernel kernel_;
    double bandwidth_;
    const char* caller_;
};

// Visits every pair i < j of a complete distance matrix, serially: reading a
// stored distance costs far less than the scoring around it.
class MatrixPairs {
public:
    MatrixPairs(const double* data, size_t n, KernelTransform transform)
        : data_(data), n_(n), transform_(transform) {}

    // body(i, j, entry), rows ascending and j ascending within a row.
    template <typename Body>
    void Run(Body&& body) {
        for (size_t i = 0; i < n_; ++i) {
            for (size_t j = i + 1; j < n_; ++j) {
                body(i, j,
                     transform_(i, j, detail::dense_distance(data_, n_, i, j)));
            }
        }
    }

private:
    const double* data_;
    size_t n_;
    KernelTransform transform_;
};

// Visits every pair i < j of a comparison. A chunk owns whole rows, so one
// worker visits all of row i's pairs, in ascending j, and body may write any
// state keyed by i without synchronization.
class ComparisonPairs {
public:
    ComparisonPairs(PairwiseComparison& comparison, size_t n, size_t num_threads,
                    size_t chunk_size, KernelTransform transform)
        : work_(comparison, n, num_threads, chunk_size),
          n_(n),
          transform_(transform) {}

    template <typename Body>
    void Run(Body&& body) {
        work_.Run(n_, [&](PairwiseComparison& local, size_t begin, size_t end) {
            for (size_t i = begin; i < end; ++i) {
                for (size_t j = i + 1; j < n_; ++j) {
                    body(i, j, transform_(i, j, local.Compare(i, j)));
                }
            }
        });
    }

private:
    detail::ChunkedComparisons work_;
    size_t n_;
    KernelTransform transform_;
};

template <typename Pairs>
std::vector<double> kernel_spectrum(Pairs& pairs, size_t n) {
    std::vector<double> kernel(detail::dense_kernel_elements(n), 0.0);
    for (size_t i = 0; i < n; ++i) {
        kernel[i * n + i] = 1.0;
    }
    // Entry (j, i) is written only by the worker that owns row i, so no two
    // workers write one element, and the kernel does not depend on how rows
    // were chunked.
    pairs.Run([&](size_t i, size_t j, double entry) {
        kernel[i * n + j] = entry;
        kernel[j * n + i] = entry;
    });
    return detail::symmetric_eigenvalues(std::move(kernel), n);
}

template <typename Pairs>
double order_two_score(Pairs& pairs, size_t n) {
    // One slot per row, summed in ascending j by the row's single owner, then
    // across rows in ascending i. The reduction order is fixed, so the score
    // is bit-identical for every thread count, chunk size and distance source.
    std::vector<double> row_sums(n, 0.0);
    pairs.Run([&](size_t i, size_t, double entry) {
        row_sums[i] += entry * entry;
    });
    double off_diagonal = 0.0;
    for (const double row_sum : row_sums) {
        off_diagonal += row_sum;
    }
    const double size = static_cast<double>(n);
    return size * size / (size + 2.0 * off_diagonal);
}

// numpy's matrix_rank rule. Duplicate items make the kernel exactly singular,
// but round-off leaves its zero eigenvalues near +-1e-16; without a tolerance
// a tiny positive one would score about -36 rather than -infinity, and differ
// by platform.
double rank_tolerance(const std::vector<double>& values) {
    double largest = 0.0;
    for (const double value : values) {
        largest = std::max(largest, std::fabs(value));
    }
    return static_cast<double>(values.size()) *
           std::numeric_limits<double>::epsilon() * largest;
}

template <typename Pairs>
VendiResult vendi(Pairs& pairs, size_t n, const VendiOptions& options) {
    VendiResult result;
    result.order = options.order;
    result.size = n;
    result.kernel = options.kernel;
    if (options.order == 2) {
        result.score = order_two_score(pairs, n);
        return result;
    }
    validate_ceiling(n, options.max_exact, VENDI_NAME);
    const std::vector<double> eigenvalues = kernel_spectrum(pairs, n);
    const double tolerance = rank_tolerance(eigenvalues);
    const double size = static_cast<double>(n);
    double entropy = 0.0;
    double negative = 0.0;
    // Eigenvalues at or below the tolerance are dropped and the rest are not
    // renormalized, as in the reference score_K; negative_mass reports what
    // the negative ones carried.
    for (const double lambda : eigenvalues) {
        if (lambda > tolerance) {
            const double p = lambda / size;
            entropy -= p * std::log(p);
        } else if (lambda < -tolerance) {
            negative += -lambda;
        }
    }
    result.score = std::exp(entropy);
    result.min_eigenvalue = eigenvalues.front();
    result.negative_mass = negative / size;
    return result;
}

template <typename Pairs>
LogDetResult logdet(Pairs& pairs, size_t n, const LogDetOptions& options) {
    validate_ceiling(n, options.max_exact, LOGDET_NAME);
    std::vector<double> shifted = kernel_spectrum(pairs, n);
    LogDetResult result;
    result.ridge = options.ridge;
    result.size = n;
    result.kernel = options.kernel;
    result.min_eigenvalue = shifted.front();
    for (double& mu : shifted) {
        mu += options.ridge;
    }
    const double tolerance = rank_tolerance(shifted);
    double sum = 0.0;
    for (const double mu : shifted) {
        if (mu <= tolerance) {
            ++result.nonpositive_count;
        } else {
            sum += std::log(mu);
        }
    }
    // Positive-definite, not log|det|: a mode at or below zero means the
    // ridged kernel is not a valid similarity kernel, even when an even number
    // of negative modes leaves the determinant positive.
    result.score = result.nonpositive_count == 0
                       ? sum
                       : -std::numeric_limits<double>::infinity();
    return result;
}

}  // namespace

VendiResult vendi_score(const StorageBackend& storage,
                        const VendiOptions& options) {
    validate_vendi_options(options);
    detail::validate_complete_distance_storage(storage, VENDI_NAME);
    const size_t n = storage.NumSamples();
    detail::validate_item_count(n, VENDI_NAME);
    MatrixPairs pairs(storage.Data(), n,
                      KernelTransform(options.kernel, options.bandwidth,
                                      VENDI_NAME));
    return vendi(pairs, n, options);
}

VendiResult vendi_score(PairwiseComparison& comparison,
                        const VendiOptions& options) {
    validate_vendi_options(options);
    detail::validate_comparison_facts(comparison, VENDI_NAME);
    const size_t n = comparison.Size();
    detail::validate_item_count(n, VENDI_NAME);
    ComparisonPairs pairs(comparison, n, options.num_threads, options.chunk_size,
                          KernelTransform(options.kernel, options.bandwidth,
                                          VENDI_NAME));
    return vendi(pairs, n, options);
}

LogDetResult logdet_diversity(const StorageBackend& storage,
                              const LogDetOptions& options) {
    validate_logdet_options(options);
    detail::validate_complete_distance_storage(storage, LOGDET_NAME);
    const size_t n = storage.NumSamples();
    detail::validate_item_count(n, LOGDET_NAME);
    MatrixPairs pairs(storage.Data(), n,
                      KernelTransform(options.kernel, options.bandwidth,
                                      LOGDET_NAME));
    return logdet(pairs, n, options);
}

LogDetResult logdet_diversity(PairwiseComparison& comparison,
                              const LogDetOptions& options) {
    validate_logdet_options(options);
    detail::validate_comparison_facts(comparison, LOGDET_NAME);
    const size_t n = comparison.Size();
    detail::validate_item_count(n, LOGDET_NAME);
    ComparisonPairs pairs(comparison, n, options.num_threads, options.chunk_size,
                          KernelTransform(options.kernel, options.bandwidth,
                                          LOGDET_NAME));
    return logdet(pairs, n, options);
}

}  // namespace OECluster
