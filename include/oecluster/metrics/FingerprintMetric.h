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
    std::string fp_type = "circular";
    unsigned int numbits = 2048;
    unsigned int min_distance = 0;
    unsigned int max_distance = 2;
    unsigned int atom_type_mask = 0;  ///< 0 = use method default
    unsigned int bond_type_mask = 0;  ///< 0 = use method default
    std::string similarity_func = "tanimoto";  ///< Similarity function name
    bool similarity = false;  ///< Return raw similarity instead of distance
};

/**
 * @brief Fingerprint-based distance metric using OEGraphSim similarity functions.
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

    /**
     * @brief Parse a pipe-delimited string of atom type names into a bitmask.
     *
     * :param pipe_delimited: Pipe-delimited atom type names (e.g. "AtomicNumber|Aromaticity").
     * :returns: Combined bitmask of atom type flags.
     * :raises MetricError: If any token is unrecognized.
     */
    static unsigned int ParseAtomTypeMask(const std::string& pipe_delimited);

    /**
     * @brief Parse a pipe-delimited string of bond type names into a bitmask.
     *
     * :param pipe_delimited: Pipe-delimited bond type names (e.g. "BondOrder|InRing").
     * :returns: Combined bitmask of bond type flags.
     * :raises MetricError: If any token is unrecognized.
     */
    static unsigned int ParseBondTypeMask(const std::string& pipe_delimited);

private:
    struct Impl;
    std::shared_ptr<const Impl> pimpl_;

    /// Private clone constructor -- shares immutable fingerprint data.
    explicit FingerprintMetric(std::shared_ptr<const Impl> impl);
};

}  // namespace OECluster

#endif  // OECLUSTER_METRICS_FINGERPRINTMETRIC_H
