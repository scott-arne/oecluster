/**
 * @file ROCSMetric.cpp
 * @brief Implementation of ROCS-style shape overlay distance metric.
 */

#include "oecluster/metrics/ROCSMetric.h"

#include <oechem.h>
#include <oeshape.h>
#include "oecluster/Error.h"

namespace OECluster {

struct ROCSMetric::SharedData {
    std::vector<std::shared_ptr<OEChem::OEMol>> mols;
};

struct ROCSMetric::ThreadLocalData {
    OEShape::OEOverlayOptions overlay_opts;
    OEShape::OEOverlay overlay;

    explicit ThreadLocalData(const OEShape::OEOverlayOptions& opts)
        : overlay_opts(opts), overlay(opts) {}
};

ROCSMetric::~ROCSMetric() = default;

void ROCSMetric::InitOverlay() {
    OEShape::OEOverlayOptions overlay_opts;

    // Configure color force field if needed
    if (opts_.score_type == ROCSScoreType::Color ||
        opts_.score_type == ROCSScoreType::Combo ||
        opts_.score_type == ROCSScoreType::ComboNorm) {
        OEShape::OEColorOptions color_opts;
        color_opts.SetColorForceField(opts_.color_ff_type);
        overlay_opts.SetColorOptions(color_opts);
    }

    local_ = std::make_unique<ThreadLocalData>(overlay_opts);
}

ROCSMetric::ROCSMetric(const std::vector<std::shared_ptr<OEChem::OEMol>>& mols,
                       const Options& opts)
    : opts_(opts) {
    auto shared = std::make_shared<SharedData>();
    shared->mols = mols;
    shared_ = std::move(shared);
    InitOverlay();
}

ROCSMetric::ROCSMetric(std::shared_ptr<const SharedData> shared,
                       const Options& opts)
    : shared_(std::move(shared)),
      opts_(opts) {
    InitOverlay();
}

double ROCSMetric::Distance(size_t i, size_t j) {
    local_->overlay.SetupRef(*shared_->mols[i]);

    OEShape::OEBestOverlayScore score;
    local_->overlay.BestOverlay(score, *shared_->mols[j]);

    if (opts_.similarity) {
        switch (opts_.score_type) {
            case ROCSScoreType::ComboNorm:
                return static_cast<double>(score.GetTanimotoCombo()) / 2.0;
            case ROCSScoreType::Combo:
                return static_cast<double>(score.GetTanimotoCombo());
            case ROCSScoreType::Shape:
                return static_cast<double>(score.GetTanimoto());
            case ROCSScoreType::Color:
                return static_cast<double>(score.GetColorTanimoto());
        }
    }

    switch (opts_.score_type) {
        case ROCSScoreType::ComboNorm:
            return 1.0 - static_cast<double>(score.GetTanimotoCombo()) / 2.0;
        case ROCSScoreType::Combo:
            return 2.0 - static_cast<double>(score.GetTanimotoCombo());
        case ROCSScoreType::Shape:
            return 1.0 - static_cast<double>(score.GetTanimoto());
        case ROCSScoreType::Color:
            return 1.0 - static_cast<double>(score.GetColorTanimoto());
    }

    return 0.0;
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
