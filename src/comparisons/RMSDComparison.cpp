/**
 * @file RMSDComparison.cpp
 * @brief Implementation of coordinate RMSD comparison.
 */

#include "oecluster/comparisons/RMSDComparison.h"

#include <cmath>
#include <oechem.h>
#include "oecluster/Error.h"

namespace OECluster {

struct RMSDComparison::SharedData {
    /**
     * @brief The structural properties one molecule had when it was validated.
     *
     * The constructor's guards are properties of the input set as it stood at
     * construction. Recording what was seen lets Compare notice that a caller
     * has since mutated a molecule out of that set.
     */
    struct ValidatedShape {
        unsigned int dimension = 0;
        unsigned int num_atoms = 0;
    };

    std::vector<std::shared_ptr<OEChem::OEMol>> mols;
    std::vector<ValidatedShape> shapes;
};

RMSDComparison::~RMSDComparison() = default;

RMSDComparison::RMSDComparison(const std::vector<std::shared_ptr<OEChem::OEMol>>& mols,
                               const Options& opts)
    : opts_(opts) {
    for (size_t i = 0; i < mols.size(); ++i) {
        if (!mols[i]) {
            throw ComparisonError("RMSDComparison received null molecule pointer at index " +
                                  std::to_string(i));
        }
        if (mols[i]->GetDimension() == 0) {
            throw ComparisonError(
                "RMSDComparison received molecule at index " + std::to_string(i) +
                " with no coordinates (GetDimension() == 0); molecules built from " +
                "SMILES without an embedding step carry no coordinates and RMSD needs them");
        }
    }

    // Require dimensional uniformity: mixing a 2D depiction with a 3D conformer is
    // not a meaningful comparison.
    if (!mols.empty()) {
        const unsigned int reference_dim = mols[0]->GetDimension();
        for (size_t i = 1; i < mols.size(); ++i) {
            if (mols[i]->GetDimension() != reference_dim) {
                throw ComparisonError(
                    "RMSDComparison received molecules with differing dimensions: item 0 has " +
                    std::to_string(reference_dim) + "D coordinates while item " +
                    std::to_string(i) + " has " + std::to_string(mols[i]->GetDimension()) +
                    "D coordinates; mixing 2D depictions with 3D conformers is not meaningful");
            }
        }
    }

    // One upfront pass rather than a per-pair check: the matrix is O(n^2) and a
    // mismatch is a property of the input set, not of a pair.
    if (!mols.empty()) {
        const std::string reference = OEChem::OEMolToSmiles(*mols[0]);
        for (size_t i = 1; i < mols.size(); ++i) {
            const std::string candidate = OEChem::OEMolToSmiles(*mols[i]);
            if (candidate != reference) {
                throw ComparisonError(
                    "RMSD requires every molecule to share one topology, but item " +
                    std::to_string(i) + " ('" + candidate + "') differs from item 0 ('" +
                    reference + "'). Use the 'rocs' comparison for cross-ligand 3D similarity");
            }
        }

        // Atom-count uniformity is required both with automorph=false (atoms matched
        // positionally) and with heavy_only=false (hydrogens counted). The topology
        // check has already established equal heavy-atom composition, so equal total
        // counts imply equal hydrogen counts.
        if (!opts.automorph || !opts.heavy_only) {
            const unsigned int reference_count = mols[0]->NumAtoms();
            for (size_t i = 1; i < mols.size(); ++i) {
                if (mols[i]->NumAtoms() != reference_count) {
                    if (!opts.automorph) {
                        throw ComparisonError(
                            "With automorph=false, atoms are matched by index, but item " +
                            std::to_string(i) + " has " + std::to_string(mols[i]->NumAtoms()) +
                            " atoms while item 0 has " + std::to_string(reference_count) +
                            "; the items must share one atom ordering");
                    } else {
                        // The automorph=true, heavy_only=false path.
                        throw ComparisonError(
                            "With heavy_only=false, hydrogens are counted, but item " +
                            std::to_string(i) + " has " + std::to_string(mols[i]->NumAtoms()) +
                            " atoms while item 0 has " + std::to_string(reference_count) +
                            "; every item must use the same hydrogen representation " +
                            "(either OEAddExplicitHydrogens on all, or heavy_only=true)");
                    }
                }
            }

            // With automorph=false, additionally require element-at-index uniformity.
            // Counts are already known equal, so indexing to reference_count is safe.
            if (!opts.automorph) {
                for (size_t i = 1; i < mols.size(); ++i) {
                    OESystem::OEIter<OEChem::OEAtomBase> ref_atom = mols[0]->GetAtoms();
                    OESystem::OEIter<OEChem::OEAtomBase> cand_atom = mols[i]->GetAtoms();
                    for (unsigned int idx = 0; idx < reference_count; ++idx, ++ref_atom, ++cand_atom) {
                        if (ref_atom->GetAtomicNum() != cand_atom->GetAtomicNum()) {
                            throw ComparisonError(
                                "With automorph=false, atoms are matched by index, but item " +
                                std::to_string(i) + " has atomic number " +
                                std::to_string(cand_atom->GetAtomicNum()) + " at index " +
                                std::to_string(idx) + " while item 0 has " +
                                std::to_string(ref_atom->GetAtomicNum()) +
                                "; the items must share one atom ordering");
                        }
                    }
                }
            }
        }
    }

    auto shared = std::make_shared<SharedData>();
    shared->mols = mols;
    shared->shapes.reserve(mols.size());
    for (const auto& mol : mols) {
        SharedData::ValidatedShape shape;
        shape.dimension = mol->GetDimension();
        shape.num_atoms = mol->NumAtoms();
        shared->shapes.push_back(shape);
    }
    shared_ = std::move(shared);
}

RMSDComparison::RMSDComparison(std::shared_ptr<const SharedData> shared, const Options& opts)
    : shared_(std::move(shared)), opts_(opts) {}

void RMSDComparison::CheckRecordedShape(size_t index) const {
    const SharedData::ValidatedShape& recorded = shared_->shapes[index];
    const OEChem::OEMol& mol = *shared_->mols[index];

    const unsigned int dimension = mol.GetDimension();
    if (dimension != recorded.dimension) {
        throw ComparisonError(
            "RMSDComparison item " + std::to_string(index) +
            " was modified after the comparison was constructed: its coordinate dimension was " +
            std::to_string(recorded.dimension) + "D when validated and is now " +
            std::to_string(dimension) +
            "D, so the comparison's validation no longer holds; rebuild the comparison from the "
            "molecules in their current state");
    }

    const unsigned int num_atoms = mol.NumAtoms();
    if (num_atoms != recorded.num_atoms) {
        throw ComparisonError(
            "RMSDComparison item " + std::to_string(index) +
            " was modified after the comparison was constructed: its atom count was " +
            std::to_string(recorded.num_atoms) + " when validated and is now " +
            std::to_string(num_atoms) +
            ", so the comparison's validation no longer holds; rebuild the comparison from the "
            "molecules in their current state");
    }
}

double RMSDComparison::Compare(size_t i, size_t j) {
    // The constructor validates the set once, but the caller keeps its own pointers
    // to these molecules. Re-check the two O(1) properties the guards rested on --
    // adding hydrogens changes the atom count, re-embedding changes the dimension --
    // so a mutated set raises rather than scoring a plausible wrong number. Both are
    // integer accessors, and measured at roughly 6ns per pair against the ~27us an
    // OERMSD automorphism search costs on benzene, so the per-pair cost is noise.
    CheckRecordedShape(i);
    CheckRecordedShape(j);

    const double value = OEChem::OERMSD(*shared_->mols[i], *shared_->mols[j],
                                        opts_.automorph, opts_.heavy_only, opts_.overlay);
    if (!std::isfinite(value)) {
        throw ComparisonError("OERMSD returned a non-finite value for items " +
                              std::to_string(i) + " and " + std::to_string(j) +
                              "; at least one of them carries a non-finite coordinate");
    }
    // -1.0 is OERMSD's documented atom-matching failure sentinel. The constructor's
    // checks are expected to prevent it: the dimension checks catch unembedded and
    // dimensionally mixed input, the topology check catches differing chemistry, and
    // the atom-count check catches differing hydrogen representation whenever counts
    // matter (automorph=false or heavy_only=false), with a per-index element check on
    // top when atoms are matched positionally. This guard is defense-in-depth: a
    // negative distance must never reach storage.
    if (value < 0.0) {
        throw ComparisonError("OERMSD could not match atoms between items " + std::to_string(i) +
                              " and " + std::to_string(j) +
                              "; the molecules share a SMILES but not a usable atom mapping");
    }
    return value;
}

std::unique_ptr<PairwiseComparison> RMSDComparison::Clone() const {
    return std::unique_ptr<PairwiseComparison>(new RMSDComparison(shared_, opts_));
}

size_t RMSDComparison::Size() const {
    return shared_->mols.size();
}

std::string RMSDComparison::ComparisonName() const {
    return "rmsd";
}

GateFacts RMSDComparison::Facts() const {
    GateFacts facts;
    // RMSD has no similarity form: every option combination returns a
    // displacement in angstroms, so the orientation is unconditional.
    facts.is_distance = Capability::Yes;
    // Zero self-distance is structural: a molecule against itself has zero
    // displacement under every combination of the options.
    facts.zero_self = Capability::Yes;
    // The triangle inequality is deliberately Unknown rather than Yes. In a
    // fixed frame RMSD is a Euclidean distance scaled by 1/sqrt(N), and a
    // minimum taken over a group acting by isometries -- the automorphism group
    // or SE(3) -- is a quotient metric, so the inequality would hold for the
    // exact minimum. But OERMSD is not documented to attain that minimum: its
    // automorphism search is not stated to be exhaustive, and an approximate
    // minimum breaks the quotient argument. Unknown is permissive at the gate,
    // so this costs no caller anything and claims only what is proven.
    facts.triangle = Capability::Unknown;
    facts.data_integrity = DataIntegrity::Complete;
    return facts;
}

}  // namespace OECluster
