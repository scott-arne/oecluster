/**
 * @file MetricTable.h
 * @brief Resolution of user-facing metric names to OEFP metrics.
 *
 * Private to the comparison implementations: this header is not installed and
 * is not exposed through SWIG.
 */

#ifndef OECLUSTER_METRICTABLE_H
#define OECLUSTER_METRICTABLE_H

#include <oefp/oefp.h>
#include <string>
#include <vector>

namespace OECluster {

/**
 * @brief Which comparison surface a metric name is being resolved for.
 *
 * Bit-set metrics are meaningless over real-valued descriptor columns, and
 * fitted metrics such as Mahalanobis are meaningless over fingerprint bits, so
 * the two surfaces expose disjoint-but-overlapping name sets.
 */
enum class MetricSurface { Fingerprint, Descriptor };

/**
 * @brief Parameters consumed by the parameterized metrics.
 *
 * Unused fields are ignored, so a caller only populates what its chosen metric
 * needs.
 */
struct MetricParams {
    double p = 2.0;                          ///< Minkowski exponent.
    double tversky_alpha = 0.5;              ///< Tversky weight on the reference side.
    double tversky_beta = 0.5;               ///< Tversky weight on the fit side.
    std::vector<double> variances;           ///< Standardized Euclidean per-column variances.
    std::vector<double> inverse_covariance;  ///< Row-major Mahalanobis inverse covariance.
};

/**
 * @brief Resolve a user-facing metric name to a concrete OEFP metric.
 *
 * The resolved metric is what the capability gate reads its facts from; the
 * user-facing string is never consulted again after this call.
 *
 * :param name: User-facing metric name; ASCII case is folded.
 * :param similarity: When true, resolve the metric's similarity form.
 * :param params: Parameters for the parameterized metrics.
 * :param surface: The comparison surface requesting the metric.
 * :returns: The resolved OEFP metric.
 * :raises ComparisonError: When the name is unknown, is not available on the
 *     requested surface, has no similarity form, or carries invalid parameters.
 */
OEFP::Metric resolve_metric(const std::string& name, bool similarity, const MetricParams& params,
                            MetricSurface surface);

/**
 * @brief List the metric names available on a surface, for error messages.
 *
 * :param surface: The comparison surface.
 * :returns: A comma-separated list in table order.
 */
std::string supported_metric_names(MetricSurface surface);

}  // namespace OECluster

#endif  // OECLUSTER_METRICTABLE_H
