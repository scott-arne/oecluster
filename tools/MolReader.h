#ifndef OEPDIST_MOLREADER_H
#define OEPDIST_MOLREADER_H

#include <algorithm>
#include <memory>
#include <string>
#include <vector>
#include <oechem.h>
#include <oebio.h>

namespace OEPDist {

/// Returns true if the path ends with a Maestro file extension.
inline bool IsMaestroFile(const std::string& path) {
    auto ends_with = [](const std::string& s, const std::string& suffix) {
        if (suffix.size() > s.size()) return false;
        return std::equal(suffix.rbegin(), suffix.rend(), s.rbegin());
    };
    return ends_with(path, ".mae") || ends_with(path, ".mae.gz")
        || ends_with(path, ".maegz");
}

struct MolSet {
    std::vector<OEChem::OEGraphMol> owned_mols;
    std::vector<OEChem::OEMolBase*> ptrs;
    std::vector<std::string> labels;
};

MolSet ReadMolecules(const std::string& path, bool verbose = false);

struct MultiConfMolSet {
    std::vector<std::shared_ptr<OEChem::OEMol>> mols;
    std::vector<std::string> labels;
};

MultiConfMolSet ReadMultiConfMolecules(const std::string& path, bool verbose = false);

struct DesignUnitSet {
    std::vector<std::shared_ptr<OEBio::OEDesignUnit>> dus;
    std::vector<std::string> labels;
};

DesignUnitSet ReadDesignUnits(const std::string& path, bool verbose = false);

}  // namespace OEPDist

#endif  // OEPDIST_MOLREADER_H
