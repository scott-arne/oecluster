/**
 * @file FingerprintMetric.h
 * @brief Distance metric based on molecular fingerprint Tanimoto similarity.
 */

#ifndef OECLUSTER_METRICS_FINGERPRINTMETRIC_H
#define OECLUSTER_METRICS_FINGERPRINTMETRIC_H

#ifdef OECLUSTER_HAS_GRAPHSIM

#include <memory>
#include <vector>
#include "oecluster/DistanceMetric.h"

namespace OEChem { class OEMolBase; }

namespace OECluster {

/**
 * @brief Configuration options for fingerprint generation.
 */
struct FingerprintOptions {
    unsigned int fp_type = 105;  ///< OEFPType::Tree
};

/**
 * @brief Fingerprint-based distance metric using OEGraphSim Tanimoto similarity.
 *
 * Computes all molecular fingerprints upfront during construction, then
 * returns ``1 - Tanimoto(fp_i, fp_j)`` as the distance. The fingerprint
 * data is immutable and shared across clones via ``std::shared_ptr``.
 */
class FingerprintMetric : public DistanceMetric {
public:
    using Options = FingerprintOptions;

    /**
     * @brief Construct a FingerprintMetric from a set of molecules.
     *
     * Generates fingerprints for every molecule in *mols* using the
     * fingerprint type specified by *opts*.
     *
     * :param mols: Pointers to molecules (not owned).
     * :param opts: Fingerprint options.
     * :raises MetricError: If fingerprint generation fails for any molecule.
     */
    explicit FingerprintMetric(const std::vector<OEChem::OEMolBase*>& mols,
                               const Options& opts = Options{});

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

#endif  // OECLUSTER_HAS_GRAPHSIM
#endif  // OECLUSTER_METRICS_FINGERPRINTMETRIC_H
