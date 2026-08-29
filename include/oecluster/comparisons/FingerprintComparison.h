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
#include "oecluster/GateFacts.h"

namespace OEChem { class OEMolBase; }

namespace OECluster {

/**
 * @brief Configuration for fingerprint-based comparison.
 *
 * Fields are grouped by the axis they control. Per-family fields are read only
 * when the matching ``fp_type`` is selected, and per-metric fields only when the
 * matching ``metric`` is selected; an unused field is ignored rather than
 * rejected, so a default-constructed struct is a working Morgan configuration.
 * The Python layer is stricter, because it can see which arguments the caller
 * actually named.
 */
struct FingerprintOptions {
    std::string fp_type = "morgan";   ///< morgan, atom_pair, topological_atom_pair, topological_torsions.
    std::string storage = "binary";   ///< binary, count, sparse, sparse_count.
    /// Folding width. The sparse storages do not fold, so this value does not
    /// shape their output -- but OEFP still validates it for the Morgan family,
    /// which rejects zero even there. Leave it nonzero.
    unsigned int numbits = 2048;
    std::string metric = "tanimoto";  ///< User-facing metric name.
    bool similarity = false;          ///< Return raw similarity instead of distance.

    unsigned int radius = 2;              ///< Morgan only.
    unsigned int min_distance = 1;        ///< Atom-pair family only (OEFP default).
    unsigned int max_distance = 30;       ///< Atom-pair family only (OEFP default).
    unsigned int torsion_atom_count = 4;  ///< Topological torsions only.
    bool use_chirality = false;           ///< Shared across families.

    double p = 2.0;              ///< Minkowski exponent; minkowski only.
    double tversky_alpha = 0.5;  ///< tversky only.
    double tversky_beta = 0.5;   ///< tversky only.
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
    GateFacts Facts() const override;

private:
    struct Impl;
    std::shared_ptr<const Impl> pimpl_;

    /// Private clone constructor -- shares immutable fingerprint data.
    explicit FingerprintComparison(std::shared_ptr<const Impl> impl);
};

}  // namespace OECluster

#endif  // OECLUSTER_COMPARISONS_FINGERPRINTCOMPARISON_H
