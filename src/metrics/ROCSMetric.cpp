/**
 * @file ROCSMetric.cpp
 * @brief Implementation of ROCS-style shape overlay distance metric.
 */

#include "oecluster/metrics/ROCSMetric.h"

#ifdef OECLUSTER_HAS_SHAPE

#include <oechem.h>
#include <oeshape.h>
#include "oecluster/Error.h"

namespace OECluster {

struct ROCSMetric::SharedData {
    std::vector<std::shared_ptr<OEChem::OEMol>> mols;
};

struct ROCSMetric::ThreadLocalData {
    OEShape::OEOverlay overlay;
};

ROCSMetric::~ROCSMetric() = default;

ROCSMetric::ROCSMetric(const std::vector<std::shared_ptr<OEChem::OEMol>>& mols,
                       const Options& opts)
    : opts_(opts) {
    auto shared = std::make_shared<SharedData>();
    shared->mols = mols;
    shared_ = std::move(shared);
    local_ = std::make_unique<ThreadLocalData>();
}

ROCSMetric::ROCSMetric(std::shared_ptr<const SharedData> shared,
                       const Options& opts)
    : shared_(std::move(shared)),
      local_(std::make_unique<ThreadLocalData>()),
      opts_(opts) {}

double ROCSMetric::Distance(size_t i, size_t j) {
    local_->overlay.SetupRef(*shared_->mols[i]);

    OEShape::OEBestOverlayScore score;
    local_->overlay.BestOverlay(score, *shared_->mols[j]);

    float sim;
    if (opts_.color_score) {
        sim = score.GetColorTanimoto();
    } else if (opts_.combo_score) {
        sim = score.GetTanimotoCombo();
    } else {
        sim = score.GetTanimoto();
    }

    return 1.0 - static_cast<double>(sim);
}

std::unique_ptr<DistanceMetric> ROCSMetric::Clone() const {
    return std::unique_ptr<DistanceMetric>(new ROCSMetric(shared_, opts_));
}

size_t ROCSMetric::Size() const {
    return shared_->mols.size();
}

std::string ROCSMetric::Name() const {
    return "rocs";
}

}  // namespace OECluster

#endif  // OECLUSTER_HAS_SHAPE
