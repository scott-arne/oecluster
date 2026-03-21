/**
 * @file SiteHopperMetric.cpp
 * @brief Implementation of binding site comparison distance metric.
 */

#include "oecluster/metrics/SiteHopperMetric.h"

#ifdef OECLUSTER_HAS_BIO

#include <oechem.h>
#include <oebio.h>
#include "oecluster/Error.h"

namespace OECluster {

struct SiteHopperMetric::SharedData {
    std::vector<OEChem::OEGraphMol> sites;
};

SiteHopperMetric::~SiteHopperMetric() = default;

SiteHopperMetric::SiteHopperMetric(
    const std::vector<std::shared_ptr<OEBio::OEDesignUnit>>& dus,
    const Options& opts)
    : opts_(opts) {
    auto shared = std::make_shared<SharedData>();
    shared->sites.reserve(dus.size());
    for (const auto& du : dus) {
        OEChem::OEGraphMol site;
        // Extract the protein component from the binding site region.
        // Use the full protein if no specific binding site component
        // is available.
        if (!du->GetComponents(site, OEBio::OEDesignUnitComponents::Protein)) {
            throw MetricError(
                "SiteHopperMetric: failed to extract protein from design unit");
        }
        shared->sites.push_back(std::move(site));
    }
    shared_ = std::move(shared);
}

SiteHopperMetric::SiteHopperMetric(
    std::shared_ptr<const SharedData> shared,
    const Options& opts)
    : shared_(std::move(shared)),
      opts_(opts) {}

double SiteHopperMetric::Distance(size_t i, size_t j) {
    // Create mutable copies since OEGetAlignment takes non-const refs
    OEChem::OEGraphMol site_i(shared_->sites[i]);
    OEChem::OEGraphMol site_j(shared_->sites[j]);

    OEBio::OESequenceAlignment alignment =
        OEBio::OEGetAlignment(site_i.SCMol(), site_j.SCMol(),
                               0xF,  // OEAssumption::Default
                               opts_.alignment_method,
                               opts_.gap_penalty,
                               opts_.extend_penalty);

    double rmsd = OEBio::OERMSD(site_i.SCMol(), site_j.SCMol(), alignment,
                                 opts_.only_calpha, true);
    return rmsd;
}

std::unique_ptr<DistanceMetric> SiteHopperMetric::Clone() const {
    return std::unique_ptr<DistanceMetric>(
        new SiteHopperMetric(shared_, opts_));
}

size_t SiteHopperMetric::Size() const {
    return shared_->sites.size();
}

std::string SiteHopperMetric::Name() const {
    return "sitehopper";
}

}  // namespace OECluster

#endif  // OECLUSTER_HAS_BIO
