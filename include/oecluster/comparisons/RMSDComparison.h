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
     * and the caller is asserting that index i names the same atom in every
     * item -- a shared canonical SMILES does not imply that. The constructor
     * verifies only that the assertion is *structurally admissible*: the same
     * element at each index, and the same bonds between the same indices among
     * the atoms being scored (``heavy_only`` decides which those are), though
     * not the same bond orders.
     *
     * That check cannot confirm the correspondence is the one intended, because
     * some wrong correspondences leave no structural trace. Atoms that are
     * symmetry-equivalent -- the three hydrogens of a methyl, the two oxygens
     * of a carboxylate -- can be permuted between two files with every element,
     * bond and canonical SMILES still agreeing, and index matching will then
     * report a real nonzero displacement between what are in fact two poses of
     * one molecule. Callers who cannot guarantee their atom order should leave
     * this true.
     */
    bool automorph = true;
    bool heavy_only = true;  ///< Skip hydrogens.
};

/**
 * @brief Coordinate RMSD comparison over molecules that share a topology.
 *
 * Every item must have the same canonical SMILES; the constructor checks this
 * once and directs callers with mixed input to the ROCS comparison. Every item
 * must also carry coordinates of the same dimension, so a 2D depiction is never
 * measured against a 3D conformer. With ``automorph=false`` every item must
 * additionally share one atom ordering (see ``RMSDOptions::automorph``), and
 * with ``heavy_only=false`` every item must share one hydrogen representation,
 * since otherwise a suppressed-hydrogen molecule would score as identical to an
 * explicit-hydrogen one whose hydrogens are displaced.
 *
 * Only each molecule's active conformer is measured: a multi-conformer
 * ``OEMol`` yields one number per molecule, not per pose. Callers wanting
 * per-pose distances should expand conformers into separate items first (the
 * Python layer does this by default with ``expand_conformers=True``).
 *
 * Molecules are never modified: ``OEChem::OERMSD`` takes both by const
 * reference, so ``overlay`` changes what is measured without writing
 * coordinates back.
 *
 * The constructor copies every input molecule and measures those copies, so
 * nothing the caller does to its own molecules afterwards can change a score.
 * That makes repeated scoring of one comparison reproducible, and it removes
 * the race in which a caller mutates a molecule while a ``pdist`` reads it --
 * it is not a general thread-safety claim, only the removal of that one race.
 * A caller who wants to score modified molecules constructs a new comparison.
 * The cost is one copy of the input set: O(n) molecules against the O(n^2)
 * distance matrix this class exists to produce.
 */
class RMSDComparison : public PairwiseComparison {
public:
    using Options = RMSDOptions;

    /**
     * @brief Construct an RMSDComparison from molecules that carry coordinates.
     *
     * :param mols: Shared pointers to molecules; each is copied into the
     *     comparison, and only each molecule's active conformer is measured. The
     *     comparison measures whatever coordinate set is present (2D or 3D).
     * :param opts: Scoring options.
     * :raises ComparisonError: When a pointer is null, an item carries no
     *     coordinates, the items' coordinate dimensions differ, an item's
     *     topology differs from the first item's, (with ``automorph=false``) the
     *     items do not share one atom ordering, or (with ``heavy_only=false``)
     *     the items do not share one hydrogen representation.
     */
    explicit RMSDComparison(const std::vector<std::shared_ptr<OEChem::OEMol>>& mols,
                            const Options& opts = Options());

    ~RMSDComparison() override;

    /**
     * @brief Measure the RMSD between items i and j.
     *
     * :param i: Index of the first item.
     * :param j: Index of the second item.
     * :returns: The RMSD in the coordinate units of the input.
     * :raises ComparisonError: When ``OEChem::OERMSD`` reports a non-finite
     *     value or an atom-matching failure.
     */
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
