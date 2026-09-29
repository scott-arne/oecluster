/**
 * @file SymmetricEigen.cpp
 * @brief Householder tridiagonalization and implicit QL, eigenvalues only.
 *
 * The eigenvalue-only forms of tred2 and tqli (Numerical Recipes, 3rd ed.,
 * section 11.4), on a row-major matrix.
 */

#include "SymmetricEigen.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace OECluster::detail {

size_t dense_kernel_elements(size_t n) {
    const size_t max = std::numeric_limits<size_t>::max();
    const bool elements_overflow = n != 0 && n > max / n;
    if (elements_overflow || n * n > max / sizeof(double) ||
        n * n > std::vector<double>().max_size()) {
        throw std::invalid_argument(
            "A dense " + std::to_string(n) + " x " + std::to_string(n) +
            " kernel is too large to allocate");
    }
    return n * n;
}

namespace {

// Reduces the lower triangle of z to tridiagonal form in place. On return
// d holds the diagonal, and e[1..n-1] the subdiagonal, with e[0] = 0. Every
// read and write is at (row, column) with column <= row.
void tridiagonalize(std::vector<double>& z, size_t n, std::vector<double>& d,
                    std::vector<double>& e) {
    const auto at = [&](size_t row, size_t column) -> double& {
        return z[row * n + column];
    };
    for (size_t i = n - 1; i > 0; --i) {
        const size_t l = i - 1;
        double h = 0.0;
        if (l > 0) {
            double scale = 0.0;
            for (size_t k = 0; k < i; ++k) {
                scale += std::fabs(at(i, k));
            }
            if (scale == 0.0) {
                e[i] = at(i, l);
            } else {
                for (size_t k = 0; k < i; ++k) {
                    at(i, k) /= scale;
                    h += at(i, k) * at(i, k);
                }
                double f = at(i, l);
                double g = f >= 0.0 ? -std::sqrt(h) : std::sqrt(h);
                e[i] = scale * g;
                h -= f * g;
                at(i, l) = f - g;
                f = 0.0;
                for (size_t j = 0; j < i; ++j) {
                    g = 0.0;
                    for (size_t k = 0; k <= j; ++k) {
                        g += at(j, k) * at(i, k);
                    }
                    for (size_t k = j + 1; k < i; ++k) {
                        g += at(k, j) * at(i, k);
                    }
                    e[j] = g / h;
                    f += e[j] * at(i, j);
                }
                const double hh = f / (h + h);
                for (size_t j = 0; j < i; ++j) {
                    f = at(i, j);
                    g = e[j] - hh * f;
                    e[j] = g;
                    for (size_t k = 0; k <= j; ++k) {
                        at(j, k) -= f * e[k] + g * at(i, k);
                    }
                }
            }
        } else {
            e[i] = at(i, l);
        }
    }
    e[0] = 0.0;
    for (size_t i = 0; i < n; ++i) {
        d[i] = at(i, i);
    }
}

// Implicit QL with Wilkinson shifts on the tridiagonal (d, e), where e[i]
// couples d[i] and d[i + 1]. Signed indices: the inner sweep runs down to l,
// which may be 0.
void tridiagonal_ql(std::vector<double>& d, std::vector<double>& e,
                    unsigned max_iterations) {
    const auto n = static_cast<std::ptrdiff_t>(d.size());
    const double epsilon = std::numeric_limits<double>::epsilon();
    for (std::ptrdiff_t l = 0; l < n; ++l) {
        unsigned iterations = 0;
        std::ptrdiff_t m = l;
        do {
            for (m = l; m < n - 1; ++m) {
                const double dd = std::fabs(d[m]) + std::fabs(d[m + 1]);
                if (std::fabs(e[m]) <= epsilon * dd) {
                    break;
                }
            }
            if (m == l) {
                break;
            }
            if (iterations++ == max_iterations) {
                throw std::runtime_error(
                    "symmetric_eigenvalues did not converge within " +
                    std::to_string(max_iterations) + " iterations");
            }
            double g = (d[l + 1] - d[l]) / (2.0 * e[l]);
            double r = std::hypot(g, 1.0);
            g = d[m] - d[l] + e[l] / (g + (g >= 0.0 ? r : -r));
            double s = 1.0;
            double c = 1.0;
            double p = 0.0;
            bool deflated = false;
            for (std::ptrdiff_t i = m - 1; i >= l; --i) {
                const double f = s * e[i];
                const double b = c * e[i];
                r = std::hypot(f, g);
                e[i + 1] = r;
                if (r == 0.0) {
                    // Underflow: the subdiagonal split early; restart this l.
                    d[i + 1] -= p;
                    e[m] = 0.0;
                    deflated = true;
                    break;
                }
                s = f / r;
                c = g / r;
                g = d[i + 1] - p;
                r = (d[i] - g) * s + 2.0 * c * b;
                p = s * r;
                d[i + 1] = g + p;
                g = c * r - b;
            }
            if (deflated) {
                continue;
            }
            d[l] -= p;
            e[l] = g;
            e[m] = 0.0;
        } while (m != l);
    }
}

}  // namespace

std::vector<double> symmetric_eigenvalues(std::vector<double> a, size_t n,
                                          unsigned max_iterations) {
    const size_t elements = dense_kernel_elements(n);
    if (a.size() != elements) {
        throw std::invalid_argument(
            "symmetric_eigenvalues expected " + std::to_string(elements) +
            " entries for n = " + std::to_string(n) + ", got " +
            std::to_string(a.size()));
    }
    if (n == 0) {
        return {};
    }
    // Power-of-two scaling guards against overflow in QL's deflation threshold.
    double max_abs = 0.0;
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j <= i; ++j) {
            max_abs = std::max(max_abs, std::fabs(a[i * n + j]));
        }
    }
    int exponent = 0;
    if (max_abs > 0.0 && std::isfinite(max_abs)) {
        std::frexp(max_abs, &exponent);
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j <= i; ++j) {
                a[i * n + j] = std::ldexp(a[i * n + j], -exponent);
            }
        }
    }
    std::vector<double> d(n, 0.0);
    std::vector<double> e(n, 0.0);
    tridiagonalize(a, n, d, e);
    // tridiagonalize leaves e[i] coupling d[i - 1] and d[i]; QL wants it
    // shifted down one place.
    for (size_t i = 1; i < n; ++i) {
        e[i - 1] = e[i];
    }
    e[n - 1] = 0.0;
    tridiagonal_ql(d, e, max_iterations);
    std::sort(d.begin(), d.end());
    if (exponent != 0) {
        for (size_t k = 0; k < n; ++k) {
            d[k] = std::ldexp(d[k], exponent);
        }
    }
    return d;
}

}  // namespace OECluster::detail
