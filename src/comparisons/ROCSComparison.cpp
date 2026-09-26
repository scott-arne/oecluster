/**
 * @file ROCSComparison.cpp
 * @brief Implementation of ROCS-style shape overlay comparison.
 */

#include "oecluster/comparisons/ROCSComparison.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <oechem.h>
#include <oeshape.h>
#include "IndexRange.h"
#include "oecluster/Error.h"

namespace OECluster {

namespace {
/// A recovered self-overlay currently lands on its saturation value exactly, so
/// this tolerance is not load-bearing on any input measured so far. It exists
/// because OEBestOverlayScore returns float and nothing in the API promises
/// that exactness; float epsilon near 1.0 is about 1.2e-7. Flipping a hard gate
/// fact on a one-ULP drift would be a worse failure than tolerating one, and at
/// 1e-6 this still sits four orders of magnitude below the smallest genuine
/// shortfall measured -- 0.0137, for diatomic chlorine.
constexpr double SELF_SCORE_TOLERANCE = 1e-6;

/// Largest per-axis span, in angstroms, that a conformer may have and still be
/// overlaid. Chosen inside a measured band rather than at its edge: 3e4 still
/// self-overlays exactly, 1e5 was the first hard failure seen (std::bad_alloc),
/// and 1e6 and above come back saturated at ~0.999, which reads like a score.
/// Real chemistry sits far below -- a large protein spans roughly 200 angstroms
/// and a virus capsid roughly 1000 -- so 1e4 clears any plausible input by an
/// order of magnitude and stays an order of magnitude under the first observed
/// failure. OEShape does not fail at 1e4, and where it does fail is
/// allocation-dependent and moved between probes; the refusal says no meaningful
/// score exists that far out, not that the toolkit would crash there.
///
/// An upper bound on catastrophe, not a certificate of quality. From somewhere
/// between 100 and 200 angstroms a stretched phenol already scores 0.546 against
/// benzene on an overlay being used for the first time and exactly 1.0 on one
/// that has scored anything before it, so cross-pair values well inside this
/// limit can still be unreliable. That split is a property of OEOverlay reuse
/// rather than of the input, and is tracked separately from this guard.
constexpr double MAX_COORDINATE_EXTENT = 1e4;  // angstroms

/// Snapshot raw molecule pointers so the shared_ptr<OEMol> constructor can own
/// the validation. Null checking happens here because the conversion below
/// dereferences before that constructor ever sees the input.
std::vector<std::shared_ptr<OEChem::OEMol>> to_oemol_snapshots(
    const std::vector<OEChem::OEMolBase*>& mols) {
    std::vector<std::shared_ptr<OEChem::OEMol>> out;
    out.reserve(mols.size());
    for (size_t i = 0; i < mols.size(); ++i) {
        if (!mols[i]) {
            throw ComparisonError("ROCSComparison received null molecule pointer at index " +
                                  std::to_string(i));
        }
        out.push_back(std::make_shared<OEChem::OEMol>(*mols[i]));
    }
    return out;
}
}  // namespace

struct ROCSComparison::SharedData {
    /// The comparison's own copies of the caller's molecules. Color preparation
    /// adds color atoms to a molecule, so the comparison cannot work on pointers
    /// the caller still holds; see the constructor for the full argument.
    std::vector<std::shared_ptr<OEChem::OEMol>> mols;

    /// The two stamps taken from the measured diagonal, cached because clones
    /// must not repeat the O(n) measurement. Both values depend on the options
    /// they were measured under. That is sound only because Clone() is the sole
    /// caller of the constructor that accepts an existing SharedData, and it
    /// always passes the same options. If SharedData is ever shared across
    /// differing options, these fields have to move onto the object.
    Capability zero_self = Capability::Unknown;
    DataIntegrity data_integrity = DataIntegrity::Complete;
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
    // The null check has to come first: the snapshot below cannot copy through a
    // null pointer.
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

    for (size_t i = 0; i < shared->mols.size(); ++i) {
        // GetDimension() on its own reports a metadata attribute that OEShape
        // never reads, and the two disagree in both directions: SetCoords leaves
        // the attribute stale, so a sound conformer can report 0, while an SDF of
        // all-zero coordinates round-trips to 3. Recomputing from the coordinates
        // makes the check test what the overlay will actually see. It mutates the
        // molecule, hence the snapshot above rather than the caller's copy.
        //
        // The refresh reads the active conformer and nothing else:
        // OESetDimensionFromCoords takes an OEMolBase&, and the OEMolBase view
        // of a multiconformer OEMol is its active conformer. So the guard is
        // permissive about the rest of the ensemble. On a six-conformer Omega
        // octane, flattening any of the five non-active conformers still
        // measures 3 and is admitted, which is right: BestOverlay picks a sound
        // conformer and the molecule still seats exactly on itself. Flattening
        // the active conformer measures 2 and is refused, which is also right,
        // but not conformer by conformer -- OEShape scores that molecule 0.0
        // against everything, as reference and as fit, with the five sound
        // conformers no help at all. A degenerate active conformer poisons the
        // whole molecule.
        OEChem::OESetDimensionFromCoords(*shared->mols[i]);

        // ``< 3``, deliberately, and not the ``== 0`` that RMSDComparison uses.
        // Planar input is as unusable here as absent input, and it fails in two
        // different ways depending on the company it keeps.
        // ``OEShape::OEOverlay::SetupRef`` refuses a degenerate reference: it
        // warns, returns ``false``, and leaves the previous reference installed.
        // Compare() reuses one overlay per thread and discards that return
        // value, so when one molecule is degenerate among sound ones -- the
        // likelier accident, one failed conformer generation in an otherwise
        // good set -- BestOverlay goes on to score the stale reference against
        // the requested fit. Measured on one reused overlay, a refused flat
        // octane reference followed by phenol returns a shape Tanimoto of
        // 0.959704, which is exactly the benzene-against-phenol score from the
        // step before it. That is not a saturation value a reader might
        // question but a plausible score belonging to a different pair, and
        // which pair depends on what that thread compared last. Each worker
        // gets its own clone and its own chunk of the pair list, so the
        // inherited pair should move with the thread count as well, but that
        // part is a deduction from the partitioning, not a measurement. Only
        // when every molecule is degenerate is there no stale reference to
        // inherit; then the toolkit finds no coordinates it can overlay and
        // returns a Tanimoto of exactly 0.0, so every pair -- the diagonal
        // included -- comes back at whichever value saturation puts it at for
        // the configured score type. The four distance forms saturate at
        // ComboNorm 1.0, Combo 2.0, Shape 1.0 and Color 1.0, which
        // MeasureDiagonal stamps No on. The four similarity forms saturate at
        // 0.0 and stamp ``zero_self = Yes``, which is a pass
        // -- but the gate refuses a similarity on ``is_distance`` before it
        // reads the diagonal at all, so both halves are refused there and the
        // similarity half is refused earlier.
        //
        // Refusing at construction is still the only way to catch this. The
        // facts a saturated similarity reports -- ``is_distance = No,
        // zero_self = Yes, triangle = Unknown, data_integrity = Complete`` --
        // are a clean bill of health with nothing in them marking the matrix as
        // pure saturation. The gated entry points refuse it for being a
        // similarity, not for being empty, and every other reader of Facts() is
        // told the numbers mean something.
        if (shared->mols[i]->GetDimension() < 3) {
            throw ComparisonError(
                "ROCSComparison requires 3D coordinates: molecule at index " +
                std::to_string(i) + " has dimension " +
                std::to_string(shared->mols[i]->GetDimension()) +
                ". Generate conformers first, for example with OEOmega.");
        }

        // Coordinates the dimension attribute cannot speak for. Finite but
        // extreme geometry reports dimension 3, reaches the self-overlay in
        // MeasureDiagonal, and either kills the process or returns a saturated
        // score; a stretched molecule still seats exactly on itself, so the
        // diagonal invariant never fires and this has to run first.
        //
        // Extent, not magnitude. A frame rigidly translated to 1e6 scores within
        // 0.002 of where it started, while a stretched one either dies or
        // saturates, so a guard on ``abs(coordinate)`` would refuse input that
        // scores perfectly well.
        //
        // Every conformer, unlike the dimension refresh above, which reads only
        // the active one. Both readings are right for their own guard:
        // BestOverlay can route around a flat non-active conformer, but it grids
        // a stretched one before it chooses, and a six-conformer octane
        // stretched to 1e10 on one non-active conformer took the process down
        // with SIGBUS.
        //
        // Finiteness before extent, and not merely for the message: ``max - min``
        // over a NaN is NaN and ``NaN > MAX_COORDINATE_EXTENT`` is false, the
        // trap MeasureDiagonal documents below, so an unchecked coordinate would
        // slip past the very test that exists to refuse it.
        for (OESystem::OEIter<OEChem::OEConfBase> conf = shared->mols[i]->GetConfs();
             conf; ++conf) {
            double lo[3] = {0.0, 0.0, 0.0};
            double hi[3] = {0.0, 0.0, 0.0};
            bool have_bounds = false;
            for (OESystem::OEIter<OEChem::OEAtomBase> atom = conf->GetAtoms(); atom; ++atom) {
                double xyz[3] = {0.0, 0.0, 0.0};
                if (!conf->GetCoords(&*atom, xyz)) {
                    throw ComparisonError(
                        "ROCSComparison could not read coordinates for molecule at index " +
                        std::to_string(i) + ".");
                }
                for (int axis = 0; axis < 3; ++axis) {
                    if (!std::isfinite(xyz[axis])) {
                        throw ComparisonError(
                            "ROCSComparison requires finite coordinates: molecule at index " +
                            std::to_string(i) +
                            " has a non-finite (NaN or infinite) coordinate. Every score "
                            "involving it, the diagonal included, would be NaN.");
                    }
                    lo[axis] = have_bounds ? std::min(lo[axis], xyz[axis]) : xyz[axis];
                    hi[axis] = have_bounds ? std::max(hi[axis], xyz[axis]) : xyz[axis];
                }
                have_bounds = true;
            }
            if (!have_bounds) {
                continue;
            }
            for (int axis = 0; axis < 3; ++axis) {
                const double extent = hi[axis] - lo[axis];
                if (extent > MAX_COORDINATE_EXTENT) {
                    throw ComparisonError(
                        "ROCSComparison requires coordinates of bounded extent: molecule at "
                        "index " + std::to_string(i) + " has a conformer spanning " +
                        std::to_string(extent) + " angstroms on one axis, above the limit of " +
                        std::to_string(MAX_COORDINATE_EXTENT) +
                        ". No meaningful overlay score exists at that scale.");
                }
            }
        }
    }

    InitOverlay(shared.get());
    shared_ = shared;  // Copy, not move: the measurement below writes through
                       // the mutable handle, and shared_ is a view onto const.
    // Writing through ``shared`` is safe because the object is not yet published:
    // construction is single-threaded and no clone can exist.
    MeasureDiagonal(*shared);
}

// Delegates rather than duplicating the null scan, the snapshot, the dimension
// and coordinate guards and the diagonal measurement above, so the two
// construction paths cannot drift. The resulting double copy -- raw pointer to
// OEMol here, then the delegate's own snapshot -- is O(n) against the O(n^2)
// overlay matrix this class exists to fill.
ROCSComparison::ROCSComparison(const std::vector<OEChem::OEMolBase*>& mols,
                               const Options& opts)
    : ROCSComparison(to_oemol_snapshots(mols), opts) {}

ROCSComparison::ROCSComparison(std::shared_ptr<const SharedData> shared,
                       const Options& opts)
    : shared_(std::move(shared)),
      opts_(opts) {
    InitOverlay(nullptr);
}

void ROCSComparison::MeasureDiagonal(SharedData& target) {
    // Reads Compare(), which is a view onto shared_; writes only ``target``.
    // The two have to be the same object or this measures one and stamps
    // another. The sole call site satisfies that, having just assigned shared_.
    assert(&target == shared_.get() &&
           "MeasureDiagonal must stamp the SharedData that shared_ views");

    bool diagonal_vanishes = true;
    // Bounded by the container Compare() indexes, not by ``target``. The assert
    // above is what documents that the two are the same object; bounding on
    // ``target`` instead would turn a divergence from a wrong stamp into an
    // out-of-range read, and the assert is compiled out in a release build.
    for (size_t i = 0; i < shared_->mols.size(); ++i) {
        const double self_score = Compare(i, i);
        if (!std::isfinite(self_score)) {
            // Tested explicitly, and before the threshold, because
            // ``std::abs(NaN) > tol`` is false: a naive threshold reads a NaN
            // diagonal as vanished and stamps a tier-1 Yes on a matrix that is
            // not a number. d(x, x) == 0 is definitively false when d(x, x) is
            // not a number, so No is provable here rather than merely cautious.
            // Returning is safe only here: both stamps are already at their
            // strongest and no later entry can change either.
            target.zero_self = Capability::No;
            target.data_integrity = DataIntegrity::NaNPresent;
            return;
        }
        if (std::abs(self_score) > SELF_SCORE_TOLERANCE) {
            // Remembered rather than returned on. A later entry cannot undo
            // this, but it can still be a NaN that has to escalate
            // data_integrity, and returning here would make that stamp depend
            // on the caller's molecule order.
            diagonal_vanishes = false;
        }
    }
    // An empty molecule set stamps Yes vacuously, which is correct: there is no
    // diagonal to violate.
    target.zero_self = diagonal_vanishes ? Capability::Yes : Capability::No;
}

double ROCSComparison::Compare(size_t i, size_t j) {
    detail::check_compare_index_range("ROCSComparison", i, j, shared_->mols.size());

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
    // overlay optimizer does not always seat a molecule exactly on itself --
    // chlorine falls 0.0137 short on shape and bromine 0.0150. All stamp No, and
    // correctly so. Measuring also means the stamp tracks the prep and the
    // toolkit rather than having to be re-derived whenever either moves.
    facts.zero_self = shared_->zero_self;

    facts.triangle = Capability::Unknown;

    // Escalated to NaNPresent only when the diagonal measurement actually
    // produced a non-finite score, which GateFacts requires of any comparison
    // that produces one. The measurement sees the diagonal and nothing else, so
    // Complete here means no NaN was observed among the n self-scores, not that
    // the O(n^2) off-diagonal pairs are proven finite.
    facts.data_integrity = shared_->data_integrity;
    return facts;
}

}  // namespace OECluster
