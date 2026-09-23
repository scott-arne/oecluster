/**
 * @file MurckoKernels.h
 * @brief Internal extraction kernel for Murcko scaffold assignment.
 *
 * src-private, following BitBirchKernels.h and KMedoidsSwapKernel.h: not
 * installed and not in the umbrella header, but included directly by the C++
 * tests so extraction correctness is unit-tested without expanding the public
 * surface.
 */

#ifndef OECLUSTER_CLUSTERING_MURCKO_KERNELS_H
#define OECLUSTER_CLUSTERING_MURCKO_KERNELS_H

#include <optional>
#include <string>

#include "oecluster/clustering/MurckoScaffold.h"

namespace OEChem {
class OEMolBase;
}

namespace OECluster::detail {

/**
 * @brief Extract one molecule's Bemis-Murcko scaffold.
 *
 * Never throws for an extraction failure: the three outcomes are distinct in
 * the return type, which is what lets the threaded caller report the *first*
 * failing index in input order rather than whichever worker lost the race.
 *
 * :param mol: Molecule to read; it is copied and never modified.
 * :param type: Extraction level.
 * :returns: The scaffold SMILES, ``""`` for a molecule with no ring bonds, or
 *     ``std::nullopt`` when a ring-containing molecule yielded no framework or
 *     any SDK transformation step reported failure.
 */
std::optional<std::string> scaffold_of(const OEChem::OEMolBase& mol,
                                       ScaffoldType type);

}  // namespace OECluster::detail

#endif  // OECLUSTER_CLUSTERING_MURCKO_KERNELS_H
