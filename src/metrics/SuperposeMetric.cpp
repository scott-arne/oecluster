/**
 * @file SuperposeMetric.cpp
 * @brief Implementation of protein superposition RMSD distance metric.
 */

#include "oecluster/metrics/SuperposeMetric.h"

#ifdef OECLUSTER_HAS_BIO

#include <oechem.h>
#include <oebio.h>
#include "oecluster/Error.h"

namespace OECluster {

struct SuperposeMetric::SharedData {
    std::vector<OEChem::OEGraphMol> mols;
};

SuperposeMetric::~SuperposeMetric() = default;

SuperposeMetric::SuperposeMetric(
    const std::vector<std::shared_ptr<OEBio::OEDesignUnit>>& dus,
    const Options& opts)
    : opts_(opts) {
    auto shared = std::make_shared<SharedData>();
    shared->mols.reserve(dus.size());
    for (const auto& du : dus) {
        OEChem::OEGraphMol mol;
        du->GetComponents(mol, OEBio::OEDesignUnitComponents::Protein);
        shared->mols.push_back(std::move(mol));
    }
    shared_ = std::move(shared);
}

SuperposeMetric::SuperposeMetric(
    const std::vector<OEChem::OEMolBase*>& mols,
    const Options& opts)
    : opts_(opts) {
    auto shared = std::make_shared<SharedData>();
    shared->mols.reserve(mols.size());
    for (auto* mol : mols) {
        shared->mols.emplace_back(*mol);
    }
    shared_ = std::move(shared);
}

SuperposeMetric::SuperposeMetric(
    std::shared_ptr<const SharedData> shared,
    const Options& opts)
    : shared_(std::move(shared)),
      opts_(opts) {}

double SuperposeMetric::Distance(size_t i, size_t j) {
    // Create mutable copies since OEGetAlignment takes non-const refs
    OEChem::OEGraphMol mol_i(shared_->mols[i]);
    OEChem::OEGraphMol mol_j(shared_->mols[j]);

    OEBio::OESequenceAlignment alignment =
        OEBio::OEGetAlignment(mol_i.SCMol(), mol_j.SCMol(),
                               0xF,  // OEAssumption::Default
                               opts_.alignment_method,
                               opts_.gap_penalty,
                               opts_.extend_penalty);

    double rmsd = OEBio::OERMSD(mol_i.SCMol(), mol_j.SCMol(), alignment,
                                 opts_.only_calpha, opts_.overlay);
    return rmsd;
}

std::unique_ptr<DistanceMetric> SuperposeMetric::Clone() const {
    return std::unique_ptr<DistanceMetric>(
        new SuperposeMetric(shared_, opts_));
}

size_t SuperposeMetric::Size() const {
    return shared_->mols.size();
}

std::string SuperposeMetric::Name() const {
    return "superpose";
}

}  // namespace OECluster

#endif  // OECLUSTER_HAS_BIO
