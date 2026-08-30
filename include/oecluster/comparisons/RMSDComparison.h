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
    bool overlay = false;  ///< Superpose before measuring.
    /**
     * Symmetry-aware atom matching. With it off, atoms are matched by index,
     * so every item must additionally share one atom ordering -- a shared
     * canonical SMILES does not imply that.
     */
    bool automorph = true;
    bool heavy_only = true;  ///< Skip hydrogens.
};

/**
 * @brief Coordinate RMSD comparison over molecules that share a topology.
 *
 * Every item must have the same canonical SMILES; the constructor checks this
 * once and directs callers with mixed input to the ROCS comparison. With
 * ``automorph=false`` every item must additionally share one atom ordering
 * (see ``RMSDOptions::automorph``).
 *
 * Only each molecule's active conformer is measured: a multi-conformer
 * ``OEMol`` yields one number per molecule, not per pose. Callers wanting
 * per-pose distances should expand conformers into separate items first (the
 * Python layer does this by default with ``expand_conformers=True``).
 *
 * Molecules are held by ``shared_ptr`` to const shared state and are never
 * modified: ``OEChem::OERMSD`` takes both molecules by const reference, so
 * ``overlay`` changes what is measured without writing coordinates back.
 */
class RMSDComparison : public PairwiseComparison {
public:
    using Options = RMSDOptions;

    /**
     * @brief Construct an RMSDComparison from molecules that carry coordinates.
     *
     * :param mols: Shared pointers to molecules; only each molecule's active
     *     conformer is measured. The comparison measures whatever coordinate set
     *     is present (2D or 3D).
     * :param opts: Scoring options.
     * :raises ComparisonError: When a pointer is null, an item carries no
     *     coordinates, an item's topology differs from the first item's, or
     *     (with ``automorph=false``) the items do not share one atom ordering.
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
