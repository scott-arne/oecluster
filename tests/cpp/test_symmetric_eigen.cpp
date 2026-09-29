/**
 * @file test_symmetric_eigen.cpp
 * @brief The private symmetric eigenvalue solver behind the set diversity
 *        scores.
 */

#include <gtest/gtest.h>

#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include "../../src/clustering/SymmetricEigen.h"

using OECluster::detail::dense_kernel_elements;
using OECluster::detail::symmetric_eigenvalues;

namespace {

// A full row-major n x n matrix from an entry function.
std::vector<double> Square(size_t n,
                           const std::function<double(size_t, size_t)>& entry) {
    std::vector<double> a(n * n);
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            a[i * n + j] = entry(i, j);
        }
    }
    return a;
}

std::vector<double> RandomSymmetric(size_t n, unsigned seed) {
    std::mt19937_64 engine(seed);
    std::uniform_real_distribution<double> uniform(-1.0, 1.0);
    std::vector<double> a(n * n);
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j <= i; ++j) {
            const double value = uniform(engine);
            a[i * n + j] = value;
            a[j * n + i] = value;
        }
    }
    return a;
}

void ExpectEigenvalues(const std::vector<double>& actual,
                       const std::vector<double>& expected, double tolerance) {
    ASSERT_EQ(actual.size(), expected.size());
    for (size_t k = 0; k < expected.size(); ++k) {
        EXPECT_NEAR(actual[k], expected[k], tolerance) << "eigenvalue " << k;
    }
}

}  // namespace

TEST(SymmetricEigenTest, AnEmptyMatrixHasNoEigenvalues) {
    EXPECT_TRUE(symmetric_eigenvalues({}, 0).empty());
}

TEST(SymmetricEigenTest, AOneByOneMatrixIsItsEntry) {
    EXPECT_EQ(symmetric_eigenvalues({-2.5}, 1), std::vector<double>{-2.5});
}

TEST(SymmetricEigenTest, TwoByTwoMatchesTheClosedForm) {
    const double root5 = std::sqrt(5.0);
    ExpectEigenvalues(symmetric_eigenvalues({2.0, 1.0, 1.0, 3.0}, 2),
                      {(5.0 - root5) / 2.0, (5.0 + root5) / 2.0}, 1e-14);
}

TEST(SymmetricEigenTest, DiagonalEntriesComeBackAscending) {
    const std::vector<double> diagonal = {3.0, -1.0, 2.0, 0.5};
    const std::vector<double> a = Square(4, [&](size_t i, size_t j) {
        return i == j ? diagonal[i] : 0.0;
    });
    EXPECT_EQ(symmetric_eigenvalues(a, 4),
              (std::vector<double>{-1.0, 0.5, 2.0, 3.0}));
}

TEST(SymmetricEigenTest, AllOnesHasOneNonzeroEigenvalue) {
    const size_t n = 6;
    ExpectEigenvalues(
        symmetric_eigenvalues(Square(n, [](size_t, size_t) { return 1.0; }), n),
        {0.0, 0.0, 0.0, 0.0, 0.0, 6.0}, 1e-12);
}

TEST(SymmetricEigenTest, TridiagonalToeplitzMatchesItsKnownSpectrum) {
    const size_t n = 50;
    const std::vector<double> a = Square(n, [](size_t i, size_t j) {
        if (i == j) {
            return 2.0;
        }
        return (i + 1 == j || j + 1 == i) ? -1.0 : 0.0;
    });
    const double pi = std::acos(-1.0);
    std::vector<double> expected;
    for (size_t k = 1; k <= n; ++k) {
        expected.push_back(2.0 - 2.0 * std::cos(static_cast<double>(k) * pi /
                                                static_cast<double>(n + 1)));
    }
    ExpectEigenvalues(symmetric_eigenvalues(a, n), expected, 1e-12);
}

TEST(SymmetricEigenTest, AnIndefiniteMatrixKeepsItsNegativeEigenvalue) {
    const double root2 = std::sqrt(2.0);
    ExpectEigenvalues(
        symmetric_eigenvalues({1, 1, 1, 1, 1, 0, 1, 0, 1}, 3),
        {1.0 - root2, 1.0, 1.0 + root2}, 1e-14);
}

TEST(SymmetricEigenTest, RepeatedEigenvaluesAreAllReturned) {
    // 3I - J has the eigenvalue 3 three times and -1 once.
    const std::vector<double> a = Square(4, [](size_t i, size_t j) {
        return (i == j ? 3.0 : 0.0) - 1.0;
    });
    ExpectEigenvalues(symmetric_eigenvalues(a, 4), {-1.0, 3.0, 3.0, 3.0},
                      1e-13);
}

TEST(SymmetricEigenTest, RandomMatricesPreserveTraceAndFrobeniusNorm) {
    for (const size_t n : {size_t{2}, size_t{3}, size_t{10}, size_t{57},
                           size_t{150}, size_t{300}}) {
        for (const unsigned seed : {1u, 2u}) {
            const std::vector<double> a = RandomSymmetric(n, seed);
            double trace = 0.0;
            double squared_norm = 0.0;
            for (size_t i = 0; i < n; ++i) {
                trace += a[i * n + i];
                for (size_t j = 0; j < n; ++j) {
                    squared_norm += a[i * n + j] * a[i * n + j];
                }
            }
            const std::vector<double> eigenvalues = symmetric_eigenvalues(a, n);
            double sum = 0.0;
            double sum_squares = 0.0;
            for (size_t k = 0; k < n; ++k) {
                sum += eigenvalues[k];
                sum_squares += eigenvalues[k] * eigenvalues[k];
                if (k > 0) {
                    EXPECT_LE(eigenvalues[k - 1], eigenvalues[k]);
                }
            }
            const double norm = std::sqrt(squared_norm);
            EXPECT_NEAR(sum, trace, 1e-10 * norm) << "n " << n;
            EXPECT_NEAR(sum_squares, squared_norm, 1e-10 * squared_norm)
                << "n " << n;
        }
    }
}

TEST(SymmetricEigenTest, TheSameInputGivesBitIdenticalOutput) {
    const std::vector<double> a = RandomSymmetric(120, 7);
    EXPECT_EQ(symmetric_eigenvalues(a, 120), symmetric_eigenvalues(a, 120));
}

TEST(SymmetricEigenTest, OnlyTheLowerTriangleIsRead) {
    const size_t n = 20;
    const std::vector<double> a = RandomSymmetric(n, 3);
    std::vector<double> garbage = a;
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i + 1; j < n; ++j) {
            garbage[i * n + j] = std::numeric_limits<double>::quiet_NaN();
        }
    }
    EXPECT_EQ(symmetric_eigenvalues(garbage, n), symmetric_eigenvalues(a, n));
}

TEST(SymmetricEigenTest, ZeroIterationsThrowsWhenQLHasWorkToDo) {
    EXPECT_THROW(symmetric_eigenvalues({2.0, 1.0, 1.0, 3.0}, 2, 0),
                 std::runtime_error);
    // A diagonal matrix needs no QL iteration, so zero is enough for it.
    EXPECT_EQ(symmetric_eigenvalues({2.0, 0.0, 0.0, 1.0}, 2, 0),
              (std::vector<double>{1.0, 2.0}));
}

TEST(SymmetricEigenTest, ASizeMismatchThrows) {
    EXPECT_THROW(symmetric_eigenvalues(std::vector<double>(8, 0.0), 3),
                 std::invalid_argument);
    EXPECT_THROW(symmetric_eigenvalues({1.0}, 0), std::invalid_argument);
}

TEST(SymmetricEigenTest, AnOverflowingSizeThrowsWithoutAllocating) {
    const size_t huge = size_t{1} << 33;  // n * n overflows size_t.
    EXPECT_THROW(dense_kernel_elements(huge), std::invalid_argument);
    EXPECT_THROW(symmetric_eigenvalues({}, huge), std::invalid_argument);
    // n * n fits, but its byte count does not.
    EXPECT_THROW(dense_kernel_elements(size_t{1} << 31), std::invalid_argument);
    try {
        dense_kernel_elements(huge);
    } catch (const std::invalid_argument& error) {
        EXPECT_EQ(std::string(error.what()),
                  "A dense 8589934592 x 8589934592 kernel is too large to "
                  "allocate");
    }
    EXPECT_EQ(dense_kernel_elements(2048), size_t{2048} * 2048);
    EXPECT_EQ(dense_kernel_elements(0), size_t{0});
}
