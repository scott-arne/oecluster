/**
 * @file RMSDComparison.cpp
 * @brief Implementation of coordinate RMSD comparison.
 */

#include "oecluster/comparisons/RMSDComparison.h"

#include <oechem.h>
#include "oecluster/Error.h"

namespace OECluster {

struct RMSDComparison::SharedData {
    std::vector<std::shared_ptr<OEChem::OEMol>> mols;
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

    // One upfront pass rather than a per-pair check: the matrix is O(n^2) and a
    // mismatch is a property of the input set, not of a pair.
    if (!mols.empty()) {
        const std::string reference = OEChem::OEMolToSmiles(*mols[0]);
        for (size_t i = 1; i < mols.size(); ++i) {
            const std::string candidate = OEChem::OEMolToSmiles(*mols[i]);
            if (candidate != reference) {
                throw ComparisonError(
                    "RMSD requires every molecule to share one topology, but item " +
                    std::to_string(i) + " ('" + candidate + "') differs from item 0 ('" +
                    reference + "'). Use the 'rocs' comparison for cross-ligand 3D similarity");
            }
        }
    }

    auto shared = std::make_shared<SharedData>();
    shared->mols = mols;
    shared_ = std::move(shared);
}

RMSDComparison::RMSDComparison(std::shared_ptr<const SharedData> shared, const Options& opts)
    : shared_(std::move(shared)), opts_(opts) {}

double RMSDComparison::Compare(size_t i, size_t j) {
    const double value = OEChem::OERMSD(*shared_->mols[i], *shared_->mols[j],
                                        opts_.automorph, opts_.heavy_only, opts_.overlay);
    // -1.0 is OERMSD's atom-matching failure sentinel. The upfront topology
    // check should make it unreachable, but matching can fail independently of
    // SMILES equality and a negative distance must never reach storage.
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
    //
    // The triangle inequality is deliberately Unknown rather than Yes. In a
    // fixed frame RMSD is a Euclidean distance scaled by 1/sqrt(N), and a
    // minimum taken over a group acting by isometries -- the automorphism group
    // or SE(3) -- is a quotient metric, so the inequality would hold for the
    // exact minimum. But OERMSD is not documented to attain that minimum: its
    // automorphism search is not stated to be exhaustive, and an approximate
    // minimum breaks the quotient argument. Unknown is permissive at the gate,
    // so this costs no caller anything and claims only what is proven.
    facts.zero_self = Capability::Yes;
    facts.triangle = Capability::Unknown;
    facts.data_integrity = DataIntegrity::Complete;
    return facts;
}

}  // namespace OECluster
