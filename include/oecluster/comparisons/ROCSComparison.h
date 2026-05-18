/**
 * @file ROCSComparison.h
 * @brief Pairwise comparison based on OEShape overlay (ROCS-style).
 */

#ifndef OECLUSTER_COMPARISONS_ROCSCOMPARISON_H
#define OECLUSTER_COMPARISONS_ROCSCOMPARISON_H

#include <memory>
#include <string>
#include <vector>
#include "oecluster/PairwiseComparison.h"

namespace OEChem { class OEMol; }

namespace OECluster {

/**
 * @brief Score type for ROCS overlay distance computation.
 */
enum class ROCSScoreType {
    ComboNorm,  ///< TanimotoCombo normalized to [0,1]: distance = 1.0 - combo/2.0
    Combo,      ///< TanimotoCombo with range [0,2]: distance = 2.0 - combo
    Shape,      ///< Shape Tanimoto with range [0,1]: distance = 1.0 - shape
    Color       ///< Color Tanimoto with range [0,1]: distance = 1.0 - color
};

/**
 * @brief Configuration options for ROCS overlay scoring.
 */
struct ROCSOptions {
    ROCSScoreType score_type = ROCSScoreType::ComboNorm;  ///< Scoring method
    unsigned int color_ff_type = 1;  ///< OEColorFFType (1=ImplicitMillsDean)
    bool similarity = false;         ///< Return raw similarity instead of distance
};

/**
 * @brief ROCS-style shape/color pairwise comparison using OEShape.
 *
 * Stores shared references to molecules and uses ``OEOverlay`` to compute
 * pairwise overlay scores. Comparison output depends on score_type and mode:
 *   - ComboNorm: ``1.0 - TanimotoCombo/2.0`` (range [0,1])
 *   - Combo: ``2.0 - TanimotoCombo`` (range [0,2])
 *   - Shape: ``1.0 - ShapeTanimoto`` (range [0,1])
 *   - Color: ``1.0 - ColorTanimoto`` (range [0,1])
 *
 * Each Clone() creates a new ``OEOverlay`` instance so that Compare()
 * can be called concurrently from different threads without locking.
 */
class ROCSComparison : public PairwiseComparison {
public:
    using Options = ROCSOptions;

    /**
     * @brief Construct a ROCSComparison from a set of molecules.
     *
     * Molecules are stored by shared_ptr and shared across clones.
     * Each molecule should have 3D coordinates for meaningful results.
     *
     * :param mols: Shared pointers to molecules.
     * :param opts: Scoring options.
     */
    explicit ROCSComparison(const std::vector<std::shared_ptr<OEChem::OEMol>>& mols,
                        const Options& opts = Options());

    ~ROCSComparison() override;

    double Compare(size_t i, size_t j) override;
    std::unique_ptr<PairwiseComparison> Clone() const override;
    size_t Size() const override;
    std::string ComparisonName() const override;

private:
    struct SharedData;
    struct ThreadLocalData;
    std::shared_ptr<const SharedData> shared_;
    std::unique_ptr<ThreadLocalData> local_;
    Options opts_;

    /// Private clone constructor -- shares molecule data, creates new overlay.
    ROCSComparison(std::shared_ptr<const SharedData> shared,
               const Options& opts);

    /// Initialize OEOverlay with configured options.
    void InitOverlay();
};

}  // namespace OECluster

#endif  // OECLUSTER_COMPARISONS_ROCSCOMPARISON_H
