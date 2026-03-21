#ifndef OEPDIST_MOLREADER_H
#define OEPDIST_MOLREADER_H

#include <memory>
#include <string>
#include <vector>
#include <oechem.h>

#ifdef OECLUSTER_HAS_BIO
#include <oebio.h>
#endif

namespace OEPDist {

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

#ifdef OECLUSTER_HAS_BIO
struct DesignUnitSet {
    std::vector<std::shared_ptr<OEBio::OEDesignUnit>> dus;
    std::vector<std::string> labels;
};

DesignUnitSet ReadDesignUnits(const std::string& path, bool verbose = false);
#endif

}  // namespace OEPDist

#endif  // OEPDIST_MOLREADER_H
