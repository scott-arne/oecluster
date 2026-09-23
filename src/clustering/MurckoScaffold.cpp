/**
 * @file MurckoScaffold.cpp
 * @brief Bemis-Murcko scaffold extraction and scaffold-keyed clustering.
 */

#include "oecluster/clustering/MurckoScaffold.h"

#include <algorithm>
#include <cstddef>
#include <iterator>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <oechem.h>
#include <oemedchem.h>
#include <oesystem.h>

#include "../DescriptorBuild.h"
#include "MurckoKernels.h"
#include "oecluster/Error.h"

namespace OECluster::detail {

namespace {

// Atom map index marking a framework atom whose bonds the sidechain cut
// removed. Map indices are the only per-atom marking that survives
// OESubsetMol, and OECreateCanSmiString emits them, so scaffold_of clears
// every one of them before canonicalizing.
constexpr unsigned int CUT_ATOM_MARK = 1u;

}  // namespace

std::optional<std::string> scaffold_of(const OEChem::OEMolBase& mol,
                                       ScaffoldType type) noexcept {
    try {
        // OEGetBemisMurcko takes a non-const molecule because it perceives ring
        // membership and writes flags. Every step below therefore runs on a
        // private copy, which is both the caller-immutability guarantee and the
        // reason the parallel loop has no shared mutable state.
        OEChem::OEGraphMol work(mol);
        OEChem::OEFindRingAtomsAndBonds(work);

        bool has_ring_bond = false;
        for (OESystem::OEIter<OEChem::OEBondBase> bond = work.GetBonds(); bond;
             ++bond) {
            if (bond->IsInRing()) {
                has_ring_bond = true;
                break;
            }
        }
        // Decided before extraction, not inferred from extraction returning
        // nothing: that distinction is what keeps an unprocessable cyclic
        // molecule from being filed as noise.
        if (!has_ring_bond) {
            return std::string();
        }

        // A record with two ring-bearing components yields more than one
        // region. Union them and subset once, so the component order comes from
        // canonicalization rather than from OEMedChem's iteration order.
        OEChem::OEAtomBondSet region_union;
        size_t region_count = 0;
        for (OESystem::OEIter<OEChem::OEAtomBondSet> region =
                 OEMedChem::OEGetBemisMurcko(work, OEMedChem::OERegionType::Framework);
             region; ++region) {
            // Add's bool reports whether the member was already present, not
            // whether the union failed, so an overlap is not an error here.
            region_union.Add(*region);
            ++region_count;
        }
        if (region_count == 0u) {
            return std::nullopt;
        }

        // OESubsetMol's hydrogen-count adjustment turns the bond order lost to
        // the cut into implicit hydrogens on the atom left behind. Which atoms
        // those are is knowable only here, while the parent still holds the
        // bonds that leave the union, so they are marked now and normalized
        // after the subset. The clearing pass is not redundant: an input read
        // from mapped SMILES arrives with map indices already set.
        for (OESystem::OEIter<OEChem::OEAtomBase> atom = work.GetAtoms(); atom;
             ++atom) {
            atom->SetMapIdx(0u);
        }
        for (OESystem::OEIter<OEChem::OEAtomBase> atom = work.GetAtoms(); atom;
             ++atom) {
            if (!region_union.HasAtom(&*atom)) {
                continue;
            }
            for (OESystem::OEIter<OEChem::OEBondBase> bond = atom->GetBonds();
                 bond; ++bond) {
                if (!region_union.HasBond(&*bond)) {
                    atom->SetMapIdx(CUT_ATOM_MARK);
                    break;
                }
            }
        }

        // OEIsMemberPtr consumes the iterator it is handed, so each predicate
        // takes its own call to GetAtoms/GetBonds.
        OEChem::OEIsMemberPtr<OEChem::OEAtomBase> atom_pred(region_union.GetAtoms());
        OEChem::OEIsMemberPtr<OEChem::OEBondBase> bond_pred(region_union.GetBonds());
        OEChem::OEGraphMol framework;
        const bool adjust_h_count = true;  // Close the valences cutting sidechains opened.
        if (!OEChem::OESubsetMol(framework, work, atom_pred, bond_pred, adjust_h_count)) {
            return std::nullopt;
        }
        // A subset that reports success but produces nothing would canonicalize
        // to the empty string, which this API reserves for "acyclic". Reporting
        // the failure keeps that one meaning.
        if (framework.NumAtoms() == 0u) {
            return std::nullopt;
        }

        // OESubsetMol carries the parent's aromatic flags onto atoms whose
        // valences the cut has just changed. Caffeine's six-membered ring is
        // aromatic only because of carbonyls that are themselves sidechains, so
        // without re-perception here the Framework level emits a string OEChem
        // cannot kekulize. This runs for both levels, before anything
        // canonicalizes the result.
        OEChem::OEFindRingAtomsAndBonds(framework);
        OEChem::OEAssignAromaticFlags(framework);

        // The normalization the marks exist for. RemoveAtomProperties is the
        // library's own recomputation of formal charge and implicit hydrogen
        // count -- the same bit the Generic level already relies on inside
        // BemisMurcko -- and the two OEChem default-valence helpers are not
        // usable here: both turn a neutral sp3 secondary amine into [NH2],
        // which would corrupt the caffeine and piperazine frameworks.
        //
        // It runs on a copy and is transferred only onto the marked atoms. A
        // charged ring atom whose every bond survived the cut is part of the
        // scaffold, and stripping it yields a neutral hypervalent atom: an
        // azoniaspiro ammonium comes back with five bonds and no charge, which
        // OEChem will parse but which is not the molecule.
        OEChem::OEGraphMol normalized(framework);
        if (!OEChem::OEUncolorMol(
                normalized, OEChem::OEUncolorStrategy::RemoveAtomProperties)) {
            return std::nullopt;
        }
        // The transfer below matches atoms by index, which holds only because
        // a framework region never contains a hydrogen: RemoveAtomProperties
        // deletes explicit hydrogens and renumbers what is left, and a shifted
        // index would silently move another atom's charge. Refusing on a count
        // mismatch keeps that latent assumption from failing quietly.
        if (normalized.NumAtoms() != framework.NumAtoms()) {
            return std::nullopt;
        }
        for (OESystem::OEIter<OEChem::OEAtomBase> atom = framework.GetAtoms();
             atom; ++atom) {
            if (atom->GetMapIdx() != CUT_ATOM_MARK) {
                continue;
            }
            const OEChem::OEAtomBase* source =
                normalized.GetAtom(OEChem::OEHasAtomIdx(atom->GetIdx()));
            if (source == nullptr) {
                return std::nullopt;
            }
            atom->SetFormalCharge(source->GetFormalCharge());
            atom->SetImplicitHCount(source->GetImplicitHCount());
        }
        // OECreateCanSmiString emits atom map indices, so clearing them is a
        // correctness step and not tidiness: leaving them turns c1ccccc1 into
        // [cH:1]1ccccc1.
        for (OESystem::OEIter<OEChem::OEAtomBase> atom = framework.GetAtoms();
             atom; ++atom) {
            atom->SetMapIdx(0u);
        }
        // Charge and hydrogen counts changed, so the flags assigned just above
        // are re-derived before either level consumes them.
        OEChem::OEFindRingAtomsAndBonds(framework);
        OEChem::OEAssignAromaticFlags(framework);

        if (type == ScaffoldType::Generic) {
            // The library strategy rather than a hand-rolled rewrite: its
            // RemoveAtomProperties bit is what recomputes implicit hydrogen
            // counts, without which benzene and pyridine reduce to two
            // different strings.
            //
            // A failed reduction must not fall through: the framework-level
            // string still in `framework` would answer a Generic request with
            // Framework output, which no caller could detect.
            if (!OEChem::OEUncolorMol(framework, OEChem::OEUncolorStrategy::BemisMurcko)) {
                return std::nullopt;
            }
            // Converting every bond to single invalidates the flags the
            // re-perception above just assigned, so it is run once more.
            OEChem::OEFindRingAtomsAndBonds(framework);
            OEChem::OEAssignAromaticFlags(framework);
        }

        // After the hydrogen-count adjustment above, never instead of it:
        // adjustment closes cut valences, suppression then normalizes how those
        // now-correct hydrogens are represented so an SD-file molecule and a
        // SMILES molecule agree. Its bool reports whether anything was
        // suppressed — it is false for a molecule that only ever had implicit
        // hydrogens — so it is not a failure signal.
        OEChem::OESuppressHydrogens(framework);

        std::string smiles;
        OEChem::OECreateCanSmiString(smiles, framework);
        return smiles;
    } catch (...) {
        // Task 4 calls this from a thread pool and reports the *first* failing
        // index in input order. An exception escaping a worker would replace
        // that deterministic report with whichever worker happened to unwind
        // first, so every failure leaves by the return value.
        return std::nullopt;
    }
}

}  // namespace OECluster::detail

namespace OECluster {

namespace {

// Per-molecule extraction cost varies little, so there is nothing for a caller
// to tune; the constant is internal and capped at the item count at the use
// site.
constexpr size_t MURCKO_CHUNK = 64;

std::vector<std::string> scaffolds_for(const std::vector<OEChem::OEMolBase*>& mols,
                                       const MurckoOptions& options,
                                       const char* caller) {
    // Options before inputs, so a caller wrong in two ways learns about the
    // enumerator rather than about their data.
    if (options.scaffold != ScaffoldType::Framework &&
        options.scaffold != ScaffoldType::Generic) {
        throw std::invalid_argument("Unknown Murcko scaffold type");
    }
    if (mols.empty()) {
        throw ComparisonError(std::string(caller) +
                              " requires at least one molecule");
    }
    const std::vector<const OEChem::OEMolBase*> inputs = checked_inputs(mols, caller);

    // The hazard in calling the toolkit concurrently is not this function's
    // data -- each worker builds its own molecule copies -- but the toolkit's
    // shared molecule memory pool. A caller who selected SingleThreaded gets
    // serial extraction rather than a data race.
    const size_t num_threads = detail::dispatch_thread_count(
        OESystem::OEGetMemPoolMode(), options.num_threads, inputs.size());
    const size_t chunk = std::min(MURCKO_CHUNK, inputs.size());
    const ScaffoldType type = options.scaffold;

    std::vector<std::optional<std::string>> raw = detail::extract_all(
        inputs.size(), num_threads, chunk,
        [&](size_t i) { return detail::scaffold_of(*inputs[i], type); });

    // The verdict is a post-join scan in input order, not an exception thrown
    // through ParallelFor: that one is captured under call_once, so the
    // survivor is a race winner and cancellation may leave later molecules
    // unexamined.
    return detail::finish_extraction(std::move(raw), caller);
}

}  // namespace

std::vector<std::string> murcko_scaffolds(const std::vector<OEChem::OEMolBase*>& mols,
                                          const MurckoOptions& options) {
    return scaffolds_for(mols, options, "murcko_scaffolds");
}

MurckoResult murcko_cluster(const std::vector<OEChem::OEMolBase*>& mols,
                            const MurckoOptions& options) {
    // Shares the labeler's implementation rather than re-extracting, so the
    // two public entry points cannot drift apart about what a scaffold is.
    // The caller name is the only thing that differs: a murcko_cluster user
    // should not be told that murcko_scaffolds refused their input.
    std::vector<std::string> scaffolds = scaffolds_for(mols, options, "murcko_cluster");

    std::vector<std::string> cluster_scaffolds;
    for (const std::string& scaffold : scaffolds) {
        if (!scaffold.empty()) {
            cluster_scaffolds.push_back(scaffold);
        }
    }
    // Sorting before assigning ranks is what makes the labeling canonical:
    // permuted inputs give the same label to the same scaffold.
    std::sort(cluster_scaffolds.begin(), cluster_scaffolds.end());
    cluster_scaffolds.erase(
        std::unique(cluster_scaffolds.begin(), cluster_scaffolds.end()),
        cluster_scaffolds.end());

    std::vector<ClusterLabel> labels(scaffolds.size(), NOISE_LABEL);
    for (size_t i = 0; i < scaffolds.size(); ++i) {
        if (scaffolds[i].empty()) {
            continue;  // Acyclic: noise, and absent from every cluster.
        }
        const auto position = std::lower_bound(cluster_scaffolds.begin(),
                                               cluster_scaffolds.end(),
                                               scaffolds[i]);
        labels[i] = static_cast<ClusterLabel>(
            std::distance(cluster_scaffolds.begin(), position));
    }

    // Ranks are dense from zero, so labels_to_clusters indexes Members() by
    // label directly and ClusterScaffolds()[i] names Members()[i].
    Clusters members = labels_to_clusters(labels);
    return MurckoResult(std::move(labels), std::move(members),
                        std::move(scaffolds), std::move(cluster_scaffolds));
}

}  // namespace OECluster
