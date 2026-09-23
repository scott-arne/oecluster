/**
 * @file MurckoScaffold.cpp
 * @brief Bemis-Murcko scaffold extraction and scaffold-keyed clustering.
 */

#include <oemedchem.h>

namespace OECluster::detail {

// Temporary link probe, and the only reason this translation unit exists in
// this commit. It forces a reference to OEMedChem so that adding the library to
// the link line is actually exercised before any behavior depends on it. The
// real extraction kernel replaces it.
const char* murcko_link_probe() {
    return OEMedChem::OEGetRegionTypeName(OEMedChem::OERegionType::Framework);
}

}  // namespace OECluster::detail
