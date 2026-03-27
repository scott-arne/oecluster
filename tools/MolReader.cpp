#include "MolReader.h"
#include <iostream>
#include <oemaestro/OEMaestroReader.h>
#include <oemaestro/OEReadMaestroDesignUnit.h>
#include <oeplatform.h>

namespace OEPDist {

/// Compute fraction from stream position, returning -1 if unavailable.
static double StreamFraction(OEPlatform::oeifstream& ifs) {
    auto total = ifs.size();
    if (total <= 0) return -1.0;
    auto pos = ifs.tell();
    return static_cast<double>(pos) / static_cast<double>(total);
}

static double StreamFraction(OEChem::oemolistream& ifs) {
    auto total = ifs.size();
    if (total <= 0) return -1.0;
    auto pos = ifs.tell();
    return static_cast<double>(pos) / static_cast<double>(total);
}

MolSet ReadMolecules(const std::string& path, bool verbose,
                     ReadProgress progress) {
    MolSet result;
    OEChem::OEGraphMol mol;
    size_t idx = 0;
    if (IsMaestroFile(path)) {
        OEPlatform::oeifstream ifs(path);
        OEMaestro::OEMaestroReader reader(ifs);
        while (reader.Read(static_cast<OEChem::OEMolBase&>(mol))) {
            result.owned_mols.push_back(mol);
            std::string title = mol.GetTitle();
            result.labels.push_back(
                title.empty() ? "mol_" + std::to_string(idx) : title);
            ++idx;
            if (progress) progress(idx, StreamFraction(ifs));
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
            if (progress) progress(idx, StreamFraction(ifs));
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

MultiConfMolSet ReadMultiConfMolecules(const std::string& path, bool verbose,
                                      ReadProgress progress) {
    MultiConfMolSet result;
    OEChem::OEMol mol;
    size_t idx = 0;
    if (IsMaestroFile(path)) {
        OEPlatform::oeifstream ifs(path);
        OEMaestro::OEMaestroReader reader(ifs);
        while (reader.Read(mol)) {
            auto ptr = std::make_shared<OEChem::OEMol>(mol);
            std::string title = mol.GetTitle();
            result.labels.push_back(
                title.empty() ? "mol_" + std::to_string(idx) : title);
            result.mols.push_back(std::move(ptr));
            ++idx;
            if (progress) progress(idx, StreamFraction(ifs));
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
            if (progress) progress(idx, StreamFraction(ifs));
        }
    }
    if (verbose) {
        std::cerr << "Read " << result.mols.size()
                  << " molecules from " << path << std::endl;
    }
    return result;
}

DesignUnitSet ReadDesignUnits(const std::string& path, bool verbose,
                             ReadProgress progress) {
    DesignUnitSet result;
    OEBio::OEDesignUnit du;
    size_t idx = 0;
    if (IsMaestroFile(path)) {
        OEPlatform::oeifstream ifs(path);
        OEMaestro::OEMaestroDesignUnitReader reader(ifs);
        while (reader.Read(du)) {
            auto ptr = std::make_shared<OEBio::OEDesignUnit>(du);
            std::string title = du.GetTitle();
            result.labels.push_back(
                title.empty() ? "du_" + std::to_string(idx) : title);
            result.dus.push_back(std::move(ptr));
            ++idx;
            if (progress) progress(idx, StreamFraction(ifs));
        }
    } else {
        OEChem::oemolistream ifs;
        if (!ifs.open(path)) {
            throw std::runtime_error("Failed to open: " + path);
        }
        while (OEBio::OEReadDesignUnit(ifs, du)) {
            auto ptr = std::make_shared<OEBio::OEDesignUnit>(du);
            std::string title = du.GetTitle();
            result.labels.push_back(
                title.empty() ? "du_" + std::to_string(idx) : title);
            result.dus.push_back(std::move(ptr));
            ++idx;
            if (progress) progress(idx, StreamFraction(ifs));
        }
    }
    if (verbose) {
        std::cerr << "Read " << result.dus.size()
                  << " design units from " << path << std::endl;
    }
    return result;
}

}  // namespace OEPDist
