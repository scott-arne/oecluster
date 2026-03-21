/**
 * @file SuperposeMetric.h
 * @brief Distance metric based on protein superposition RMSD using OEBio.
 */

#ifndef OECLUSTER_METRICS_SUPERPOSEMETRIC_H
#define OECLUSTER_METRICS_SUPERPOSEMETRIC_H

#ifdef OECLUSTER_HAS_BIO

#include <memory>
#include <string>
#include <vector>
#include "oecluster/DistanceMetric.h"

namespace OEBio { class OEDesignUnit; }
namespace OEChem { class OEMolBase; }

namespace OECluster {

/**
 * @brief Configuration options for protein superposition RMSD.
 */
struct SuperposeOptions {
    unsigned int alignment_method = 2;  ///< OESeqAlignmentMethod (2=PAM250)
    int gap_penalty = -10;              ///< Gap penalty for sequence alignment
    int extend_penalty = -2;            ///< Gap extension penalty
    bool only_calpha = true;            ///< Use only C-alpha atoms for RMSD
    bool overlay = true;                ///< Superpose before computing RMSD
};

/**
 * @brief Protein superposition distance metric using OEBio OERMSD.
 *
 * Computes pairwise RMSD between protein structures after sequence
 * alignment and optimal overlay. Accepts either design units (from
 * which the protein component is extracted) or molecules directly.
 *
 * Distance values are unbounded (not normalized to [0,1]).
 *
 * Each Clone() shares the immutable protein structures. Distance()
 * creates temporary copies since OEGetAlignment requires non-const refs.
 */
class SuperposeMetric : public DistanceMetric {
public:
    using Options = SuperposeOptions;

    /**
     * @brief Construct from design units (extracts protein components).
     *
     * :param dus: Shared pointers to design units.
     * :param opts: Superposition options.
     */
    explicit SuperposeMetric(const std::vector<std::shared_ptr<OEBio::OEDesignUnit>>& dus,
                             const Options& opts = Options{});

    /**
     * @brief Construct from molecules directly.
     *
     * Molecules are copied internally. The pointers are only used
     * during construction.
     *
     * :param mols: Pointers to molecules with 3D coordinates.
     * :param opts: Superposition options.
     */
    explicit SuperposeMetric(const std::vector<OEChem::OEMolBase*>& mols,
                             const Options& opts = Options{});

    ~SuperposeMetric() override;

    double Distance(size_t i, size_t j) override;
    std::unique_ptr<DistanceMetric> Clone() const override;
    size_t Size() const override;
    std::string Name() const override;

private:
    struct SharedData;
    std::shared_ptr<const SharedData> shared_;
    Options opts_;

    /// Private clone constructor -- shares molecule data.
    SuperposeMetric(std::shared_ptr<const SharedData> shared, const Options& opts);
};

}  // namespace OECluster

#endif  // OECLUSTER_HAS_BIO
#endif  // OECLUSTER_METRICS_SUPERPOSEMETRIC_H
