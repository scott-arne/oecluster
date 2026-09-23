/**
 * @file MurckoScaffold.cpp
 * @brief Bemis-Murcko scaffold extraction and scaffold-keyed clustering.
 */

#include "oecluster/clustering/MurckoScaffold.h"

#include <algorithm>
#include <cstddef>
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

}  // namespace

std::vector<std::string> murcko_scaffolds(const std::vector<OEChem::OEMolBase*>& mols,
                                          const MurckoOptions& options) {
    // Options before inputs, so a caller wrong in two ways learns about the
    // enumerator rather than about their data.
    if (options.scaffold != ScaffoldType::Framework &&
        options.scaffold != ScaffoldType::Generic) {
        throw std::invalid_argument("Unknown Murcko scaffold type");
    }
    if (mols.empty()) {
        throw ComparisonError("murcko_scaffolds requires at least one molecule");
    }
    const std::vector<const OEChem::OEMolBase*> inputs =
        checked_inputs(mols, "murcko_scaffolds");

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
    if (const std::optional<size_t> bad = detail::first_failure(raw)) {
        throw ComparisonError(
            "murcko_scaffolds could not extract a scaffold for molecule at index " +
            std::to_string(*bad));
    }

    std::vector<std::string> scaffolds;
    scaffolds.reserve(raw.size());
    for (std::optional<std::string>& value : raw) {
        scaffolds.push_back(std::move(*value));
    }
    return scaffolds;
}

}  // namespace OECluster
