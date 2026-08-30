/**
 * @file RMSDComparison.h
 * @brief Pairwise comparison by root-mean-square deviation of coordinates.
 */

#ifndef OECLUSTER_COMPARISONS_RMSDCOMPARISON_H
#define OECLUSTER_COMPARISONS_RMSDCOMPARISON_H

#include <cstddef>
#include <memory>
#include <string>
#include <vector>
#include "oecluster/PairwiseComparison.h"

namespace OEChem { class OEMol; }

namespace OECluster {

/**
 * @brief Configuration options for RMSD scoring.
 *
 * The defaults are OEChem's own, so the default call measures symmetry-aware
 * heavy-atom RMSD in the input frame -- correct for docked poses that already
 * share a receptor frame. Conformer work sets ``overlay``.
 */
struct RMSDOptions {
    bool overlay = false;    ///< Superpose before measuring.
    bool automorph = true;   ///< Symmetry-aware atom matching.
    bool heavy_only = true;  ///< Skip hydrogens.
};

/**
 * @brief Coordinate RMSD comparison over molecules that share a topology.
 *
 * Every item must have the same canonical SMILES; the constructor checks this
 * once and directs callers with mixed input to the ROCS comparison. Multi-
 * conformer input is expanded to one item per conformer by the Python layer
 * before the comparison is built.
 *
 * Molecules are held by ``shared_ptr`` to const shared state and are never
 * modified: ``OEChem::OERMSD`` takes both molecules by const reference, so
 * ``overlay`` changes what is measured without writing coordinates back.
 */
class RMSDComparison : public PairwiseComparison {
public:
    using Options = RMSDOptions;

    /**
     * @brief Construct an RMSDComparison from a set of single-conformer molecules.
     *
     * :param mols: Shared pointers to molecules with 3D coordinates.
     * :param opts: Scoring options.
     * :raises ComparisonError: When a pointer is null or an item's topology
     *     differs from the first item's.
     */
    explicit RMSDComparison(const std::vector<std::shared_ptr<OEChem::OEMol>>& mols,
                            const Options& opts = Options());

    ~RMSDComparison() override;

    double Compare(size_t i, size_t j) override;
    std::unique_ptr<PairwiseComparison> Clone() const override;
    size_t Size() const override;
    std::string ComparisonName() const override;
    GateFacts Facts() const override;

private:
    struct SharedData;
    std::shared_ptr<const SharedData> shared_;
    Options opts_;

    /// Private clone constructor -- shares the immutable molecule list.
    RMSDComparison(std::shared_ptr<const SharedData> shared, const Options& opts);
};

}  // namespace OECluster

#endif  // OECLUSTER_COMPARISONS_RMSDCOMPARISON_H
