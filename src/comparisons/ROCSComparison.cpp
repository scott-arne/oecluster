/**
 * @file ROCSComparison.cpp
 * @brief Implementation of ROCS-style shape overlay comparison.
 */

#include "oecluster/comparisons/ROCSComparison.h"

#include <cmath>
#include <oechem.h>
#include <oeshape.h>
#include "oecluster/Error.h"

namespace OECluster {

namespace {
/// A recovered self-overlay currently lands on its saturation value exactly, so
/// this tolerance is not load-bearing on any input measured so far. It exists
/// because OEBestOverlayScore returns float and nothing in the API promises
/// that exactness; float epsilon near 1.0 is about 1.2e-7. Flipping a hard gate
/// fact on a one-ULP drift would be a worse failure than tolerating one, and at
/// 1e-6 this still sits three orders of magnitude below the smallest genuine
/// shortfall measured -- 0.0137, for diatomic chlorine.
constexpr double SELF_SCORE_TOLERANCE = 1e-6;
}  // namespace

struct ROCSComparison::SharedData {
    /// The comparison's own copies of the caller's molecules. Color preparation
    /// adds color atoms to a molecule, so the comparison cannot work on pointers
    /// the caller still holds; see the constructor for the full argument.
    std::vector<std::shared_ptr<OEChem::OEMol>> mols;

    /// The measured diagonal, cached because clones must not repeat the O(n)
    /// measurement. This value depends on the options it was measured under.
    /// That is sound only because Clone() is the sole caller of the constructor
    /// that accepts an existing SharedData, and it always passes the same
    /// options. If SharedData is ever shared across differing options, this
    /// field has to move onto the object.
    Capability zero_self = Capability::Unknown;
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
        // ``< 3``, deliberately, and not the ``== 0`` that RMSDComparison uses:
        // OEGenerate2DCoordinates leaves the dimension at 2, and a planar molecule
        // has no shape volume, so every Tanimoto degenerates to 0/0 and every score
        // -- including a self-comparison -- comes back 1.0. Aligning this with the
        // sibling check would let that through silently.
        if (mols[i]->GetDimension() < 3) {
            throw ComparisonError(
                "ROCSComparison requires 3D coordinates: molecule at index " +
                std::to_string(i) + " has dimension " +
                std::to_string(mols[i]->GetDimension()) +
                ". Generate conformers first, for example with OEOmega.");
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
    shared_ = shared;  // Copy, not move: the measurement below writes through
                       // the mutable handle, and shared_ is a view onto const.
    // Writing through ``shared`` is safe because the object is not yet published:
    // construction is single-threaded and no clone can exist.
    shared->zero_self = MeasureZeroSelf();
}

ROCSComparison::ROCSComparison(std::shared_ptr<const SharedData> shared,
                       const Options& opts)
    : shared_(std::move(shared)),
      opts_(opts) {
    InitOverlay(nullptr);
}

Capability ROCSComparison::MeasureZeroSelf() {
    for (size_t i = 0; i < shared_->mols.size(); ++i) {
        if (std::abs(Compare(i, i)) > SELF_SCORE_TOLERANCE) {
            return Capability::No;
        }
    }
    // An empty molecule set stamps Yes vacuously, which is correct: there is no
    // diagonal to violate.
    return Capability::Yes;
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

    // Reported from the measurement the constructor took on this molecule set
    // under these options, not asserted from the score type or the direction
    // flag. No general claim about the diagonal would be true: the color force
    // field gives a small molecule such as methane no color atom at all, so its
    // color self-Tanimoto is 0.0 and its ComboNorm diagonal sits at 0.5, and the
    // overlay optimizer cannot seat a near-spherical diatomic exactly on itself,
    // so chlorine falls 0.0137 short on shape. Both stamp No, and correctly so.
    // Measuring also means the stamp tracks the prep and the toolkit rather than
    // having to be re-derived whenever either moves.
    facts.zero_self = shared_->zero_self;

    facts.triangle = Capability::Unknown;
    facts.data_integrity = DataIntegrity::Complete;
    return facts;
}

}  // namespace OECluster
