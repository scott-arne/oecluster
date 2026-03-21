/**
 * @file SiteHopperMetric.h
 * @brief Distance metric based on binding site comparison using OEBio.
 */

#ifndef OECLUSTER_METRICS_SITEHOPPERMETRIC_H
#define OECLUSTER_METRICS_SITEHOPPERMETRIC_H

#ifdef OECLUSTER_HAS_BIO

#include <memory>
#include <string>
#include <vector>
#include "oecluster/DistanceMetric.h"

namespace OEBio { class OEDesignUnit; }

namespace OECluster {

/**
 * @brief Configuration options for binding site comparison.
 */
struct SiteHopperOptions {
    unsigned int alignment_method = 2;  ///< OESeqAlignmentMethod (2=PAM250)
    int gap_penalty = -10;              ///< Gap penalty for sequence alignment
    int extend_penalty = -2;            ///< Gap extension penalty
    bool only_calpha = true;            ///< Use only C-alpha atoms for RMSD
};

/**
 * @brief Binding site distance metric using protein RMSD from OEBio.
 *
 * Extracts binding site protein components from design units and computes
 * pairwise RMSD between them using sequence alignment and superposition.
 *
 * Distance values are unbounded (not normalized to [0,1]).
 *
 * Each Clone() shares the immutable binding site structures. Distance()
 * creates temporary copies since OEGetAlignment requires non-const refs.
 */
class SiteHopperMetric : public DistanceMetric {
public:
    using Options = SiteHopperOptions;

    /**
     * @brief Construct from design units (extracts binding site components).
     *
     * :param dus: Shared pointers to design units.
     * :param opts: Site comparison options.
     */
    explicit SiteHopperMetric(const std::vector<std::shared_ptr<OEBio::OEDesignUnit>>& dus,
                              const Options& opts = Options{});

    ~SiteHopperMetric() override;

    double Distance(size_t i, size_t j) override;
    std::unique_ptr<DistanceMetric> Clone() const override;
    size_t Size() const override;
    std::string Name() const override;

private:
    struct SharedData;
    std::shared_ptr<const SharedData> shared_;
    Options opts_;

    /// Private clone constructor -- shares site data.
    SiteHopperMetric(std::shared_ptr<const SharedData> shared, const Options& opts);
};

}  // namespace OECluster

#endif  // OECLUSTER_HAS_BIO
#endif  // OECLUSTER_METRICS_SITEHOPPERMETRIC_H
