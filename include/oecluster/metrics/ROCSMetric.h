/**
 * @file ROCSMetric.h
 * @brief Distance metric based on OEShape overlay (ROCS-style).
 */

#ifndef OECLUSTER_METRICS_ROCSMETRIC_H
#define OECLUSTER_METRICS_ROCSMETRIC_H

#ifdef OECLUSTER_HAS_SHAPE

#include <memory>
#include <string>
#include <vector>
#include "oecluster/DistanceMetric.h"

namespace OEChem { class OEMol; }

namespace OECluster {

/**
 * @brief Configuration options for ROCS overlay scoring.
 */
struct ROCSOptions {
    bool color_score = false;   ///< Use color Tanimoto only (range [0,1])
    bool combo_score = true;    ///< Use TanimotoCombo (shape + color, range [0,2])
};

/**
 * @brief ROCS-style shape/color overlay distance metric using OEShape.
 *
 * Stores shared references to molecules and uses ``OEOverlay`` to compute
 * pairwise overlay scores. Distance is ``1.0 - best_score`` where
 * ``best_score`` is the TanimotoCombo, color Tanimoto, or shape Tanimoto
 * depending on options.
 *
 * Each Clone() creates a new ``OEOverlay`` instance so that Distance()
 * can be called concurrently from different threads without locking.
 */
class ROCSMetric : public DistanceMetric {
public:
    using Options = ROCSOptions;

    /**
     * @brief Construct a ROCSMetric from a set of molecules.
     *
     * Molecules are stored by shared_ptr and shared across clones.
     * Each molecule should have 3D coordinates for meaningful results.
     *
     * :param mols: Shared pointers to molecules.
     * :param opts: Scoring options.
     */
    explicit ROCSMetric(const std::vector<std::shared_ptr<OEChem::OEMol>>& mols,
                        const Options& opts = Options{});

    ~ROCSMetric() override;

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

    /// Private clone constructor -- shares molecule data, creates new overlay.
    ROCSMetric(std::shared_ptr<const SharedData> shared,
               const Options& opts);
};

}  // namespace OECluster

#endif  // OECLUSTER_HAS_SHAPE
#endif  // OECLUSTER_METRICS_ROCSMETRIC_H
