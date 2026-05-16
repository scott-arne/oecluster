/**
 * @file FingerprintMetric.h
 * @brief Distance metric based on molecular fingerprint similarity.
 */

#ifndef OECLUSTER_METRICS_FINGERPRINTMETRIC_H
#define OECLUSTER_METRICS_FINGERPRINTMETRIC_H

#include <memory>
#include <string>
#include <vector>
#include "oecluster/DistanceMetric.h"

namespace OEChem { class OEMolBase; }

namespace OECluster {

/**
 * @brief Configuration options for fingerprint generation.
 */
struct FingerprintOptions {
    std::string fp_type = "morgan";
    unsigned int numbits = 2048;
    unsigned int min_distance = 0;
    unsigned int max_distance = 2;
    std::string similarity_func = "tanimoto";  ///< Similarity function name
    bool similarity = false;  ///< Return raw similarity instead of distance
};

/**
 * @brief Fingerprint-based distance metric using OEFP fingerprints.
 *
 * Computes all molecular fingerprints upfront during construction, then
 * returns a distance value based on the configured similarity metric.
 * The fingerprint data is immutable and shared across clones via
 * ``std::shared_ptr``.
 */
class FingerprintMetric : public DistanceMetric {
public:
    using Options = FingerprintOptions;

    /**
     * @brief Construct a FingerprintMetric from a set of molecules.
     *
     * Generates fingerprints for every molecule in *mols* using the
     * fingerprint type and parameters specified by *opts*.
     *
     * :param mols: Pointers to molecules (not owned).
     * :param opts: Fingerprint options.
     * :raises MetricError: If fingerprint generation fails for any molecule.
     */
    explicit FingerprintMetric(const std::vector<OEChem::OEMolBase*>& mols,
                               const Options& opts = Options());

    double Distance(size_t i, size_t j) override;
    std::unique_ptr<DistanceMetric> Clone() const override;
    size_t Size() const override;
    std::string Name() const override;

private:
    struct Impl;
    std::shared_ptr<const Impl> pimpl_;

    /// Private clone constructor -- shares immutable fingerprint data.
    explicit FingerprintMetric(std::shared_ptr<const Impl> impl);
};

}  // namespace OECluster

#endif  // OECLUSTER_METRICS_FINGERPRINTMETRIC_H
