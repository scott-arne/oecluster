#include "MolReader.h"
#include <iostream>
#include <oemaestro/OEMaestroReader.h>

namespace OEPDist {

MolSet ReadMolecules(const std::string& path, bool verbose) {
    MolSet result;
    OEChem::OEGraphMol mol;
    size_t idx = 0;
    if (IsMaestroFile(path)) {
        OEMaestro::OEMaestroReader reader(path);
        while (reader.Read(static_cast<OEChem::OEMolBase&>(mol))) {
            result.owned_mols.push_back(mol);
            std::string title = mol.GetTitle();
            result.labels.push_back(
                title.empty() ? "mol_" + std::to_string(idx) : title);
            ++idx;
        }
    } else {
        OEChem::oemolistream ifs;
        if (!ifs.open(path)) {
            throw std::runtime_error("Failed to open: " + path);
        }
        while (OEChem::OEReadMolecule(ifs, mol)) {
            result.owned_mols.push_back(mol);
            std::string title = mol.GetTitle();
            result.labels.push_back(
                title.empty() ? "mol_" + std::to_string(idx) : title);
            ++idx;
        }
    }
    result.ptrs.reserve(result.owned_mols.size());
    for (auto& m : result.owned_mols) {
        result.ptrs.push_back(&static_cast<OEChem::OEMolBase&>(m));
    }
    if (verbose) {
        std::cerr << "Read " << result.owned_mols.size()
                  << " molecules from " << path << std::endl;
    }
    return result;
}

MultiConfMolSet ReadMultiConfMolecules(const std::string& path, bool verbose) {
    MultiConfMolSet result;
    OEChem::OEMol mol;
    size_t idx = 0;
    if (IsMaestroFile(path)) {
        OEMaestro::OEMaestroReader reader(path);
        while (reader.Read(mol)) {
            auto ptr = std::make_shared<OEChem::OEMol>(mol);
            std::string title = mol.GetTitle();
            result.labels.push_back(
                title.empty() ? "mol_" + std::to_string(idx) : title);
            result.mols.push_back(std::move(ptr));
            ++idx;
        }
    } else {
        OEChem::oemolistream ifs;
        if (!ifs.open(path)) {
            throw std::runtime_error("Failed to open: " + path);
        }
        while (OEChem::OEReadMolecule(ifs, mol)) {
            auto ptr = std::make_shared<OEChem::OEMol>(mol);
            std::string title = mol.GetTitle();
            result.labels.push_back(
                title.empty() ? "mol_" + std::to_string(idx) : title);
            result.mols.push_back(std::move(ptr));
            ++idx;
        }
    }
    if (verbose) {
        std::cerr << "Read " << result.mols.size()
                  << " molecules from " << path << std::endl;
    }
    return result;
}

DesignUnitSet ReadDesignUnits(const std::string& path, bool verbose) {
    DesignUnitSet result;
    OEChem::oemolistream ifs;
    if (!ifs.open(path)) {
        throw std::runtime_error("Failed to open: " + path);
    }
    OEBio::OEDesignUnit du;
    size_t idx = 0;
    while (OEBio::OEReadDesignUnit(ifs, du)) {
        auto ptr = std::make_shared<OEBio::OEDesignUnit>(du);
        std::string title = du.GetTitle();
        result.labels.push_back(title.empty() ? "du_" + std::to_string(idx) : title);
        result.dus.push_back(std::move(ptr));
        ++idx;
    }
    if (verbose) {
        std::cerr << "Read " << result.dus.size()
                  << " design units from " << path << std::endl;
    }
    return result;
}

}  // namespace OEPDist
