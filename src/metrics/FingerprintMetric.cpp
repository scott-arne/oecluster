/**
 * @file FingerprintMetric.cpp
 * @brief Implementation of fingerprint-based Tanimoto distance metric.
 */

#include "oecluster/metrics/FingerprintMetric.h"

#ifdef OECLUSTER_HAS_GRAPHSIM

#include <oechem.h>
#include <oegraphsim.h>
#include "oecluster/Error.h"

namespace OECluster {

struct FingerprintMetric::Impl {
    std::vector<OEGraphSim::OEFingerPrint> fingerprints;
};

FingerprintMetric::FingerprintMetric(const std::vector<OEChem::OEMolBase*>& mols,
                                     const Options& opts) {
    auto impl = std::make_shared<Impl>();
    impl->fingerprints.resize(mols.size());
    for (size_t i = 0; i < mols.size(); ++i) {
        if (!OEGraphSim::OEMakeFP(impl->fingerprints[i], *mols[i], opts.fp_type)) {
            throw MetricError("Failed to compute fingerprint for molecule " +
                              std::to_string(i));
        }
    }
    pimpl_ = std::move(impl);
}

FingerprintMetric::FingerprintMetric(std::shared_ptr<const Impl> impl)
    : pimpl_(std::move(impl)) {}

double FingerprintMetric::Distance(size_t i, size_t j) {
    float sim = OEGraphSim::OETanimoto(pimpl_->fingerprints[i],
                                       pimpl_->fingerprints[j]);
    return 1.0 - static_cast<double>(sim);
}

std::unique_ptr<DistanceMetric> FingerprintMetric::Clone() const {
    return std::unique_ptr<DistanceMetric>(new FingerprintMetric(pimpl_));
}

size_t FingerprintMetric::Size() const {
    return pimpl_->fingerprints.size();
}

std::string FingerprintMetric::Name() const {
    return "fingerprint";
}

}  // namespace OECluster

#endif  // OECLUSTER_HAS_GRAPHSIM
