/**
 * @file FingerprintComparison.h
 * @brief Pairwise comparison based on molecular fingerprint similarity.
 */

#ifndef OECLUSTER_COMPARISONS_FINGERPRINTCOMPARISON_H
#define OECLUSTER_COMPARISONS_FINGERPRINTCOMPARISON_H

#include <memory>
#include <string>
#include <vector>
#include "oecluster/PairwiseComparison.h"

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
    std::string metric = "tanimoto";  ///< OEFP scalar metric name.
    bool similarity = false;  ///< Return raw similarity instead of distance
};

/**
 * @brief Fingerprint-based pairwise comparison using OEFP fingerprints.
 *
 * Computes all molecular fingerprints upfront during construction, then
 * delegates scalar and batch comparisons to OEFP. The fingerprint data is
 * immutable and shared across clones via ``std::shared_ptr``.
 */
class FingerprintComparison : public PairwiseComparison {
public:
    using Options = FingerprintOptions;

    /**
     * @brief Construct a FingerprintComparison from a set of molecules.
     *
     * Generates fingerprints for every molecule in *mols* using the
     * fingerprint type and parameters specified by *opts*.
     *
     * :param mols: Pointers to molecules (not owned).
     * :param opts: Fingerprint options.
     * :raises ComparisonError: If fingerprint generation fails for any molecule.
     */
    explicit FingerprintComparison(const std::vector<OEChem::OEMolBase*>& mols,
                               const Options& opts = Options());

    double Compare(size_t i, size_t j) override;
    bool TryPDist(StorageBackend& storage, const PDistOptions& options) override;
    bool TryCDist(size_t n_a, double* output, const CDistOptions& options) override;
    std::unique_ptr<PairwiseComparison> Clone() const override;
    size_t Size() const override;
    std::string ComparisonName() const override;

private:
    struct Impl;
    std::shared_ptr<const Impl> pimpl_;

    /// Private clone constructor -- shares immutable fingerprint data.
    explicit FingerprintComparison(std::shared_ptr<const Impl> impl);
};

}  // namespace OECluster

#endif  // OECLUSTER_COMPARISONS_FINGERPRINTCOMPARISON_H
