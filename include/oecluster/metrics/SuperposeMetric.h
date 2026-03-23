/**
 * @file SuperposeMetric.h
 * @brief Distance metric based on protein superposition using oespruce OESuperpose.
 */

#ifndef OECLUSTER_METRICS_SUPERPOSEMETRIC_H
#define OECLUSTER_METRICS_SUPERPOSEMETRIC_H

#include <memory>
#include <string>
#include <vector>
#include "oecluster/DistanceMetric.h"

namespace OEBio { class OEDesignUnit; }
namespace OEChem { class OEMolBase; }

namespace OECluster {

/**
 * @brief Superposition method for protein overlay.
 */
enum class SuperposeMethod {
    GlobalCarbonAlpha,  ///< All matched alpha carbon atoms (default)
    Global,             ///< Global superposition
    DDM,                ///< Distance Difference Matrix
    Weighted,           ///< Weighted DDM
    SSE,                ///< Secondary Structure Elements (Tanimoto score)
    SiteHopper          ///< SiteHopper patch score
};

/**
 * @brief Score type selection for superposition results.
 */
enum class SuperposeScoreType {
    Auto,       ///< Natural score for the method
    RMSD,       ///< Root mean square deviation
    Tanimoto,   ///< Tanimoto coefficient [0,1]
    PatchScore  ///< SiteHopper patch score [0,4]
};

/**
 * @brief Configuration options for protein superposition.
 */
struct SuperposeOptions {
    SuperposeMethod method = SuperposeMethod::GlobalCarbonAlpha;
    SuperposeScoreType score_type = SuperposeScoreType::Auto;
    bool similarity = false;
    std::string predicate;      ///< oeselect expression for both ref and fit
    std::string ref_predicate;  ///< Override predicate for ref
    std::string fit_predicate;  ///< Override predicate for fit
};

/**
 * @brief Protein superposition distance metric using oespruce OESuperpose.
 *
 * Supports multiple superposition methods (GlobalCarbonAlpha, Global, DDM,
 * Weighted, SSE, SiteHopper) with configurable score types and atom predicates
 * via oeselect expressions.
 *
 * Score semantics per method:
 * - RMSD group (Global, GlobalCarbonAlpha, DDM, Weighted): raw RMSD
 * - Tanimoto group (SSE): 1.0 - tanimoto (distance), raw tanimoto (similarity)
 * - PatchScore group (SiteHopper): 4.0 - patch_score (distance), raw (similarity)
 *
 * Each Clone() creates a new thread-local OESuperpose instance.
 */
class SuperposeMetric : public DistanceMetric {
public:
    using Options = SuperposeOptions;

    /**
     * @brief Construct from design units.
     *
     * :param dus: Shared pointers to design units.
     * :param opts: Superposition options.
     */
    explicit SuperposeMetric(
        const std::vector<std::shared_ptr<OEBio::OEDesignUnit>>& dus,
        const Options& opts = Options());

    /**
     * @brief Construct from molecules directly.
     *
     * Molecules are copied internally.
     *
     * :param mols: Pointers to molecules with 3D coordinates.
     * :param opts: Superposition options.
     */
    explicit SuperposeMetric(
        const std::vector<OEChem::OEMolBase*>& mols,
        const Options& opts = Options());

    ~SuperposeMetric() override;

    double Distance(size_t i, size_t j) override;
    std::unique_ptr<DistanceMetric> Clone() const override;
    size_t Size() const override;
    std::string Name() const override;

private:
    struct SharedData;
    struct ThreadLocalData;
    std::shared_ptr<const SharedData> shared_;
    std::unique_ptr<ThreadLocalData> local_;
    Options opts_;

    /// Private clone constructor -- shares structure data, creates new OESuperpose.
    SuperposeMetric(std::shared_ptr<const SharedData> shared, const Options& opts);

    /// Initialize thread-local OESuperpose with configured options.
    void InitSuperpose();
};

}  // namespace OECluster

#endif  // OECLUSTER_METRICS_SUPERPOSEMETRIC_H
