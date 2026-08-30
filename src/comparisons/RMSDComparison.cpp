/**
 * @file RMSDComparison.cpp
 * @brief Implementation of coordinate RMSD comparison.
 */

#include "oecluster/comparisons/RMSDComparison.h"

#include <algorithm>
#include <cmath>
#include <map>
#include <utility>
#include <oechem.h>
#include "oecluster/Error.h"

namespace OECluster {

namespace {

/// A bond written in terms of atom iteration positions:
/// ``(lower position, higher position)``. Endpoints are normalized low-to-high
/// so that the two directions of one bond compare equal. The bond's order is
/// deliberately not part of this; see the guard that consumes it.
using IndexedBond = std::pair<unsigned int, unsigned int>;

/// The molecule's bonds as a sorted set of index-wise endpoint pairs, which two
/// molecules can be compared on exactly.
///
/// Positions come from walking ``GetAtoms()`` rather than from
/// ``OEAtomBase::GetIdx()``: index values can carry gaps once atoms have been
/// deleted, and this guard exists to protect exactly that kind of edited
/// molecule. Walking the iterator also keeps the numbering identical to the one
/// the element-at-index check uses, so the two guards speak about the same
/// positions.
std::vector<IndexedBond> indexed_bonds(const OEChem::OEMol& mol) {
    std::map<const OEChem::OEAtomBase*, unsigned int> positions;
    unsigned int position = 0;
    for (OESystem::OEIter<OEChem::OEAtomBase> atom = mol.GetAtoms(); atom; ++atom, ++position) {
        positions[&*atom] = position;
    }

    std::vector<IndexedBond> bonds;
    bonds.reserve(mol.NumBonds());
    for (OESystem::OEIter<OEChem::OEBondBase> bond = mol.GetBonds(); bond; ++bond) {
        // at() rather than operator[]: an endpoint outside the molecule's own atom
        // list would be a toolkit invariant violation, and silently folding it into
        // position 0 would make this correctness guard quietly wrong.
        const unsigned int begin = positions.at(bond->GetBgn());
        const unsigned int end = positions.at(bond->GetEnd());
        bonds.push_back(IndexedBond(std::min(begin, end), std::max(begin, end)));
    }
    std::sort(bonds.begin(), bonds.end());
    return bonds;
}

}  // namespace

struct RMSDComparison::SharedData {
    /// The comparison's own copies of the caller's molecules; see the class
    /// documentation for why the comparison owns them.
    std::vector<std::shared_ptr<const OEChem::OEMol>> mols;
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
    }

    // Snapshot the input set before validating it, so that every guard below is a
    // statement about the molecules that will actually be scored. The caller keeps
    // its own molecules and may reorder, re-embed or protonate them afterwards
    // without any of it reaching us. The copy is O(n) against the O(n^2) matrix this
    // class exists to fill, and an OEMol copy is an order of magnitude cheaper than
    // the single OERMSD call each pair already costs.
    auto shared = std::make_shared<SharedData>();
    shared->mols.reserve(mols.size());
    for (const auto& mol : mols) {
        shared->mols.push_back(std::make_shared<const OEChem::OEMol>(*mol));
    }
    const std::vector<std::shared_ptr<const OEChem::OEMol>>& owned = shared->mols;

    for (size_t i = 0; i < owned.size(); ++i) {
        if (owned[i]->GetDimension() == 0) {
            throw ComparisonError(
                "RMSDComparison received molecule at index " + std::to_string(i) +
                " with no coordinates (GetDimension() == 0); molecules built from " +
                "SMILES without an embedding step carry no coordinates and RMSD needs them");
        }
    }

    // Require dimensional uniformity: mixing a 2D depiction with a 3D conformer is
    // not a meaningful comparison.
    if (!owned.empty()) {
        const unsigned int reference_dim = owned[0]->GetDimension();
        for (size_t i = 1; i < owned.size(); ++i) {
            if (owned[i]->GetDimension() != reference_dim) {
                throw ComparisonError(
                    "RMSDComparison received molecules with differing dimensions: item 0 has " +
                    std::to_string(reference_dim) + "D coordinates while item " +
                    std::to_string(i) + " has " + std::to_string(owned[i]->GetDimension()) +
                    "D coordinates; mixing 2D depictions with 3D conformers is not meaningful");
            }
        }
    }

    // One upfront pass rather than a per-pair check: the matrix is O(n^2) and a
    // mismatch is a property of the input set, not of a pair.
    if (!owned.empty()) {
        const std::string reference = OEChem::OEMolToSmiles(*owned[0]);
        for (size_t i = 1; i < owned.size(); ++i) {
            const std::string candidate = OEChem::OEMolToSmiles(*owned[i]);
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
            const unsigned int reference_count = owned[0]->NumAtoms();
            for (size_t i = 1; i < owned.size(); ++i) {
                if (owned[i]->NumAtoms() != reference_count) {
                    if (!opts.automorph) {
                        throw ComparisonError(
                            "With automorph=false, atoms are matched by index, but item " +
                            std::to_string(i) + " has " + std::to_string(owned[i]->NumAtoms()) +
                            " atoms while item 0 has " + std::to_string(reference_count) +
                            "; the items must share one atom ordering");
                    } else {
                        // The automorph=true, heavy_only=false path.
                        throw ComparisonError(
                            "With heavy_only=false, hydrogens are counted, but item " +
                            std::to_string(i) + " has " + std::to_string(owned[i]->NumAtoms()) +
                            " atoms while item 0 has " + std::to_string(reference_count) +
                            "; every item must use the same hydrogen representation " +
                            "(either OEAddExplicitHydrogens on all, or heavy_only=true)");
                    }
                }
            }

            // With automorph=false, additionally require element-at-index uniformity.
            // Counts are already known equal, so indexing to reference_count is safe.
            if (!opts.automorph) {
                for (size_t i = 1; i < owned.size(); ++i) {
                    OESystem::OEIter<OEChem::OEAtomBase> ref_atom = owned[0]->GetAtoms();
                    OESystem::OEIter<OEChem::OEAtomBase> cand_atom = owned[i]->GetAtoms();
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

                // Matching elements at every index still does not establish a shared
                // atom ordering: explicit hydrogens can be distributed identically as
                // elements while attaching to different heavy atoms, which leaves the
                // atom count, the canonical SMILES and the element sequence all equal
                // and index-matched RMSD measuring unrelated pairs of atoms. Comparing
                // the bonds by index closes that gap. It costs O(atoms + bonds log
                // bonds) once per item, against an O(n^2) matrix of OERMSD calls.
                //
                // Bond *order* is deliberately excluded, and must stay excluded. The
                // question here is only whether index i names the same atom in every
                // item, and coordinate RMSD never reads a bond order. Chemistry is the
                // canonical-SMILES check's job, and it has already run: any order
                // difference that changes the molecule changes the SMILES and is
                // rejected there, so the only differences that survive to this point
                // are kekulizations of one aromatic system. Two kekulizations of one
                // ligand are one atom ordering, and refusing them would be a false
                // refusal on input two writers can easily produce from one molecule.
                const std::vector<IndexedBond> reference_bonds = indexed_bonds(*owned[0]);
                for (size_t i = 1; i < owned.size(); ++i) {
                    const std::vector<IndexedBond> candidate_bonds = indexed_bonds(*owned[i]);
                    if (candidate_bonds.size() != reference_bonds.size()) {
                        throw ComparisonError(
                            "With automorph=false, atoms are matched by index, but item " +
                            std::to_string(i) + " has " +
                            std::to_string(candidate_bonds.size()) + " bonds while item 0 has " +
                            std::to_string(reference_bonds.size()) +
                            "; the items must share one atom ordering. Use automorph=true "
                            "for symmetry-aware matching");
                    }
                    if (candidate_bonds != reference_bonds) {
                        throw ComparisonError(
                            "With automorph=false, atoms are matched by index, but item " +
                            std::to_string(i) +
                            " bonds different pairs of indices than item 0; matching "
                            "elements at each index is not enough, because two files whose "
                            "atom columns agree can still attach their hydrogens to "
                            "different heavy atoms, so the items do not share one atom "
                            "ordering. Use automorph=true for symmetry-aware matching");
                    }
                }
            }
        }
    }

    shared_ = std::move(shared);
}

RMSDComparison::RMSDComparison(std::shared_ptr<const SharedData> shared, const Options& opts)
    : shared_(std::move(shared)), opts_(opts) {}

double RMSDComparison::Compare(size_t i, size_t j) {
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
