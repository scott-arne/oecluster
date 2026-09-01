/**
 * @file ROCSComparison.cpp
 * @brief Implementation of ROCS-style shape overlay comparison.
 */

#include "oecluster/comparisons/ROCSComparison.h"

#include <oechem.h>
#include <oeshape.h>
#include "oecluster/Error.h"

namespace OECluster {

struct ROCSComparison::SharedData {
    /// The comparison's own copies of the caller's molecules. Color preparation
    /// adds color atoms to a molecule, so the comparison cannot work on pointers
    /// the caller still holds; see the constructor for the full argument.
    std::vector<std::shared_ptr<OEChem::OEMol>> mols;
};

struct ROCSComparison::ThreadLocalData {
    OEShape::OEOverlayOptions overlay_opts;
    OEShape::OEOverlay overlay;

    explicit ThreadLocalData(const OEShape::OEOverlayOptions& opts)
        : overlay_opts(opts), overlay(opts) {}
};

ROCSComparison::~ROCSComparison() = default;

void ROCSComparison::InitOverlay(SharedData* prep_target) {
    OEShape::OEOverlayOptions overlay_opts;

    // Configure color force field if needed
    if (opts_.score_type == ROCSScoreType::Color ||
        opts_.score_type == ROCSScoreType::Combo ||
        opts_.score_type == ROCSScoreType::ComboNorm) {
        OEShape::OEColorOptions color_opts;
        color_opts.SetColorForceField(opts_.color_ff_type);
        overlay_opts.SetColorOptions(color_opts);

        // Setting the force field on the overlay options is not enough on its own:
        // it tells the overlay how to score color atoms, but nothing puts color
        // atoms on the molecules, so GetColorTanimoto() reads 0.0 for every pair
        // including a molecule against itself. OEOverlapPrep is what assigns them.
        //
        // SetUseHydrogens(true) overrides a prep default of false. Stripping
        // hydrogens changes the shape term as well as the color term, and it moves
        // a self-overlay off 1.0 -- phenol scores 0.98975 against itself with
        // hydrogens off -- which would make the diagonal that Facts() stamps
        // approximate rather than exact. Keeping hydrogens also reproduces the
        // shape numbers this comparison returned before color prep existed, so it
        // is the narrower behavior change of the two.
        //
        // Only the constructor that owns the molecules preps. Clones share an
        // already-prepped SharedData, and Compare() is called O(n^2) times against
        // it from several threads, so prepping anywhere but here would be both
        // quadratic and a data race.
        if (prep_target != nullptr) {
            OEShape::OEOverlapPrep prep;
            prep.SetAssignColor(true);
            prep.SetUseHydrogens(true);
            prep.Initialize(color_opts);
            for (const auto& mol : prep_target->mols) {
                prep.Prep(*mol);
            }
        }
    }

    local_ = std::make_unique<ThreadLocalData>(overlay_opts);
}

ROCSComparison::ROCSComparison(const std::vector<std::shared_ptr<OEChem::OEMol>>& mols,
                       const Options& opts)
    : opts_(opts) {
    for (size_t i = 0; i < mols.size(); ++i) {
        if (!mols[i]) {
            throw ComparisonError("ROCSComparison received null molecule pointer at index " +
                                  std::to_string(i));
        }
    }

    // Snapshot the input set rather than aliasing the caller's shared_ptrs. Color
    // preparation writes color atoms into the molecule, so sharing the caller's
    // pointers would mean constructing this object silently edits objects the
    // caller still owns. The copy is O(n) against the O(n^2) overlay matrix this
    // class exists to fill.
    auto shared = std::make_shared<SharedData>();
    shared->mols.reserve(mols.size());
    for (const auto& mol : mols) {
        shared->mols.push_back(std::make_shared<OEChem::OEMol>(*mol));
    }

    InitOverlay(shared.get());
    shared_ = std::move(shared);
}

ROCSComparison::ROCSComparison(std::shared_ptr<const SharedData> shared,
                       const Options& opts)
    : shared_(std::move(shared)),
      opts_(opts) {
    InitOverlay(nullptr);
}

double ROCSComparison::Compare(size_t i, size_t j) {
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

std::unique_ptr<PairwiseComparison> ROCSComparison::Clone() const {
    return std::unique_ptr<PairwiseComparison>(new ROCSComparison(shared_, opts_));
}

size_t ROCSComparison::Size() const {
    return shared_->mols.size();
}

std::string ROCSComparison::ComparisonName() const {
    return "rocs";
}

GateFacts ROCSComparison::Facts() const {
    GateFacts facts;
    facts.is_distance = opts_.similarity ? Capability::No : Capability::Yes;

    // Every score type now saturates on the diagonal -- shape and color Tanimoto
    // both reach 1.0 for a molecule against itself, so combo reaches 2.0 -- and
    // each distance form subtracts exactly that saturation value. The diagonal is
    // therefore zero in the distance direction for all four score types, and is
    // the saturation value rather than zero in the similarity direction.
    //
    // This rests on the prep InitOverlay applies, not on a mathematical guarantee.
    // BestOverlay optimizes from inertial-frame starting poses and returns the
    // best it finds; it is not obliged to find the identity transform even when a
    // molecule is overlaid on itself. Under the current prep it does recover it,
    // which is what makes the diagonal exact -- with hydrogens stripped it
    // measurably does not, and phenol self-scores 0.98975 on shape. Change the
    // prep and this stamp has to be re-measured.
    facts.zero_self = opts_.similarity ? Capability::No : Capability::Yes;

    facts.triangle = Capability::Unknown;
    facts.data_integrity = DataIntegrity::Complete;
    return facts;
}

}  // namespace OECluster
