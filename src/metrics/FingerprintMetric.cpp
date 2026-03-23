/**
 * @file FingerprintMetric.cpp
 * @brief Implementation of fingerprint-based distance metric with configurable
 *        fingerprint types and similarity functions.
 */

#include "oecluster/metrics/FingerprintMetric.h"

#include <algorithm>
#include <cctype>
#include <sstream>
#include <unordered_map>
#include <oechem.h>
#include <oegraphsim.h>
#include "oecluster/Error.h"

namespace OECluster {

// ---------------------------------------------------------------------------
// Impl
// ---------------------------------------------------------------------------

struct FingerprintMetric::Impl {
    std::vector<OEGraphSim::OEFingerPrint> fingerprints;
    FingerprintOptions opts;
};

// ---------------------------------------------------------------------------
// Helper: lowercase a string in place
// ---------------------------------------------------------------------------

static std::string to_lower(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char c) { return std::tolower(c); });
    return s;
}

// ---------------------------------------------------------------------------
// Helper: trim whitespace from both ends
// ---------------------------------------------------------------------------

static std::string trim(const std::string& s) {
    auto start = s.find_first_not_of(" \t\r\n");
    if (start == std::string::npos) return "";
    auto end = s.find_last_not_of(" \t\r\n");
    return s.substr(start, end - start + 1);
}

// ---------------------------------------------------------------------------
// Constructor
// ---------------------------------------------------------------------------

FingerprintMetric::FingerprintMetric(const std::vector<OEChem::OEMolBase*>& mols,
                                     const Options& opts) {
    auto impl = std::make_shared<Impl>();
    impl->opts = opts;
    impl->fingerprints.resize(mols.size());

    std::string fp_lower = to_lower(opts.fp_type);

    // Resolve default atom/bond masks for the three parameterized FP types.
    unsigned int atom_mask = opts.atom_type_mask;
    unsigned int bond_mask = opts.bond_type_mask;

    if (fp_lower == "circular") {
        if (atom_mask == 0) atom_mask = OEGraphSim::OEFPAtomType::DefaultCircularAtom;
        if (bond_mask == 0) bond_mask = OEGraphSim::OEFPBondType::DefaultCircularBond;
    } else if (fp_lower == "tree") {
        if (atom_mask == 0) atom_mask = OEGraphSim::OEFPAtomType::DefaultTreeAtom;
        if (bond_mask == 0) bond_mask = OEGraphSim::OEFPBondType::DefaultTreeBond;
    } else if (fp_lower == "path") {
        if (atom_mask == 0) atom_mask = OEGraphSim::OEFPAtomType::DefaultPathAtom;
        if (bond_mask == 0) bond_mask = OEGraphSim::OEFPBondType::DefaultPathBond;
    }

    for (size_t i = 0; i < mols.size(); ++i) {
        bool ok = false;
        if (fp_lower == "circular") {
            ok = OEGraphSim::OEMakeCircularFP(impl->fingerprints[i], *mols[i],
                                               opts.numbits, opts.min_distance,
                                               opts.max_distance, atom_mask,
                                               bond_mask);
        } else if (fp_lower == "tree") {
            ok = OEGraphSim::OEMakeTreeFP(impl->fingerprints[i], *mols[i],
                                           opts.numbits, opts.min_distance,
                                           opts.max_distance, atom_mask,
                                           bond_mask);
        } else if (fp_lower == "path") {
            ok = OEGraphSim::OEMakePathFP(impl->fingerprints[i], *mols[i],
                                           opts.numbits, opts.min_distance,
                                           opts.max_distance, atom_mask,
                                           bond_mask);
        } else if (fp_lower == "maccs") {
            ok = OEGraphSim::OEMakeMACCS166FP(impl->fingerprints[i], *mols[i]);
        } else if (fp_lower == "lingo") {
            ok = OEGraphSim::OEMakeLingoFP(impl->fingerprints[i], *mols[i]);
        } else {
            throw MetricError("Unknown fingerprint type: " + opts.fp_type);
        }

        if (!ok) {
            throw MetricError("Failed to compute fingerprint for molecule " +
                              std::to_string(i));
        }
    }

    // Validate: similarity mode is incompatible with absolute distance functions
    std::string sim_check = to_lower(opts.similarity_func);
    if (opts.similarity && (sim_check == "manhattan" || sim_check == "euclidean")) {
        throw MetricError("similarity=true is incompatible with distance function: "
                          + opts.similarity_func);
    }

    pimpl_ = std::move(impl);
}

// ---------------------------------------------------------------------------
// Private clone constructor
// ---------------------------------------------------------------------------

FingerprintMetric::FingerprintMetric(std::shared_ptr<const Impl> impl)
    : pimpl_(std::move(impl)) {}

// ---------------------------------------------------------------------------
// Distance
// ---------------------------------------------------------------------------

double FingerprintMetric::Distance(size_t i, size_t j) {
    const auto& fp_i = pimpl_->fingerprints[i];
    const auto& fp_j = pimpl_->fingerprints[j];
    std::string sim_lower = to_lower(pimpl_->opts.similarity_func);
    bool sim_mode = pimpl_->opts.similarity;

    if (sim_lower == "tanimoto") {
        double sim = static_cast<double>(OEGraphSim::OETanimoto(fp_i, fp_j));
        return sim_mode ? sim : (1.0 - sim);
    } else if (sim_lower == "dice") {
        double sim = static_cast<double>(OEGraphSim::OEDice(fp_i, fp_j));
        return sim_mode ? sim : (1.0 - sim);
    } else if (sim_lower == "cosine") {
        double sim = static_cast<double>(OEGraphSim::OECosine(fp_i, fp_j));
        return sim_mode ? sim : (1.0 - sim);
    } else if (sim_lower == "manhattan") {
        return static_cast<double>(OEGraphSim::OEManhattan(fp_i, fp_j));
    } else if (sim_lower == "euclidean") {
        return static_cast<double>(OEGraphSim::OEEuclid(fp_i, fp_j));
    }

    throw MetricError("Unknown similarity function: " + pimpl_->opts.similarity_func);
}

// ---------------------------------------------------------------------------
// Clone / Size / Name
// ---------------------------------------------------------------------------

std::unique_ptr<DistanceMetric> FingerprintMetric::Clone() const {
    return std::unique_ptr<DistanceMetric>(new FingerprintMetric(pimpl_));
}

size_t FingerprintMetric::Size() const {
    return pimpl_->fingerprints.size();
}

std::string FingerprintMetric::Name() const {
    return "fingerprint";
}

// ---------------------------------------------------------------------------
// ParseAtomTypeMask
// ---------------------------------------------------------------------------

unsigned int FingerprintMetric::ParseAtomTypeMask(const std::string& pipe_delimited) {
    static const std::unordered_map<std::string, unsigned int> atom_map = {
        {"atomicnumber",      OEGraphSim::OEFPAtomType::AtomicNumber},
        {"aromaticity",       OEGraphSim::OEFPAtomType::Aromaticity},
        {"chiral",            OEGraphSim::OEFPAtomType::Chiral},
        {"formalcharge",      OEGraphSim::OEFPAtomType::FormalCharge},
        {"hvydegree",         OEGraphSim::OEFPAtomType::HvyDegree},
        {"hybridization",     OEGraphSim::OEFPAtomType::Hybridization},
        {"inring",            OEGraphSim::OEFPAtomType::InRing},
        {"hcount",            OEGraphSim::OEFPAtomType::HCount},
        {"eqhalogen",         OEGraphSim::OEFPAtomType::EqHalogen},
        {"eqaromatic",        OEGraphSim::OEFPAtomType::EqAromatic},
        {"eqhbondacceptor",   OEGraphSim::OEFPAtomType::EqHBondAcceptor},
        {"eqhbonddonor",      OEGraphSim::OEFPAtomType::EqHBondDonor},
        {"defaultatom",       OEGraphSim::OEFPAtomType::DefaultAtom},
        {"defaultpathatom",   OEGraphSim::OEFPAtomType::DefaultPathAtom},
        {"defaultcircularatom", OEGraphSim::OEFPAtomType::DefaultCircularAtom},
        {"defaulttreeatom",   OEGraphSim::OEFPAtomType::DefaultTreeAtom},
    };

    unsigned int mask = 0;
    std::istringstream stream(pipe_delimited);
    std::string token;
    while (std::getline(stream, token, '|')) {
        std::string key = to_lower(trim(token));
        // Strip "oefpatomtype::" or "oefpatomtype_" prefix if present
        for (const auto& prefix : {"oefpatomtype::", "oefpatomtype_"}) {
            std::string p(prefix);
            if (key.size() > p.size() && key.compare(0, p.size(), p) == 0) {
                key = key.substr(p.size());
                break;
            }
        }
        auto it = atom_map.find(key);
        if (it == atom_map.end()) {
            throw MetricError("Unknown atom type token: " + trim(token));
        }
        mask |= it->second;
    }
    return mask;
}

// ---------------------------------------------------------------------------
// ParseBondTypeMask
// ---------------------------------------------------------------------------

unsigned int FingerprintMetric::ParseBondTypeMask(const std::string& pipe_delimited) {
    static const std::unordered_map<std::string, unsigned int> bond_map = {
        {"bondorder",          OEGraphSim::OEFPBondType::BondOrder},
        {"chiral",             OEGraphSim::OEFPBondType::Chiral},
        {"inring",             OEGraphSim::OEFPBondType::InRing},
        {"defaultbond",        OEGraphSim::OEFPBondType::DefaultBond},
        {"defaultpathbond",    OEGraphSim::OEFPBondType::DefaultPathBond},
        {"defaultcircularbond", OEGraphSim::OEFPBondType::DefaultCircularBond},
        {"defaulttreebond",    OEGraphSim::OEFPBondType::DefaultTreeBond},
    };

    unsigned int mask = 0;
    std::istringstream stream(pipe_delimited);
    std::string token;
    while (std::getline(stream, token, '|')) {
        std::string key = to_lower(trim(token));
        // Strip "oefpbondtype::" or "oefpbondtype_" prefix if present
        for (const auto& prefix : {"oefpbondtype::", "oefpbondtype_"}) {
            std::string p(prefix);
            if (key.size() > p.size() && key.compare(0, p.size(), p) == 0) {
                key = key.substr(p.size());
                break;
            }
        }
        auto it = bond_map.find(key);
        if (it == bond_map.end()) {
            throw MetricError("Unknown bond type token: " + trim(token));
        }
        mask |= it->second;
    }
    return mask;
}

}  // namespace OECluster
