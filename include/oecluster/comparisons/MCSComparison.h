/**
 * @file MCSComparison.h
 * @brief Pairwise comparison by maximum common substructure.
 */

#ifndef OECLUSTER_COMPARISONS_MCSCOMPARISON_H
#define OECLUSTER_COMPARISONS_MCSCOMPARISON_H

#include <cstddef>
#include <memory>
#include <string>
#include <vector>
#include "oecluster/PairwiseComparison.h"

namespace OEChem { class OEMol; }
namespace OEChem { class OEMolBase; }

namespace OECluster {

/**
 * @brief Which OEChem MCS algorithm runs.
 *
 * The default deliberately inverts the toolkit's own: ``OEMCSType::Default`` is
 * exhaustive, but exhaustive search is one to three orders of magnitude slower
 * and is not reliably better. On a 53-bond against 54-bond macrolide pair it
 * took 16.2 s and returned 50 matched bonds where approximate took 8.5 ms and
 * returned 51.
 */
enum class MCSSearchMode {
    Approximate,  ///< Fast heuristic search. The default.
    Exhaustive    ///< Complete search, subject to the toolkit's own truncation limit.
};

/**
 * @brief How strictly atoms and bonds must correspond to count as matched.
 *
 * Three presets rather than raw ``OEExprOpts`` bitmasks, so the public surface
 * names chemistry rather than toolkit flags.
 *
 * ``Default`` already constrains bond order, so ``Exact`` does not add it:
 * ethane against ethene scores 1.0 under ``Loose`` and 0.0 under ``Default``.
 * What ``Exact`` adds is ring membership on both atoms and bonds, hydrogen
 * count, degree, and strict rather than lenient formal charge. Ring membership
 * is the one that surprises: cyclohexane against hexane scores 0.83 under
 * ``Default`` and 0.0 under ``Exact``.
 *
 * ``ExactAtoms`` also sets the isotope and chirality bits, but neither was
 * observed to change a score. Isotope matching is directional -- an unlabelled
 * pattern matches a labelled target -- and ``Compare`` takes the larger of the
 * two directions, so the unlabelled direction wins: 13C-butane against butane
 * scores 1.0 at every level. Do not rely on either to separate stereoisomers
 * or labelled analogues.
 */
enum class MCSMatchLevel {
    Default,  ///< OEChem's default atom and bond expressions.
    Exact,    ///< Adds ring membership, hydrogen count, degree, strict charge.
    Loose     ///< Atomic number only, with bonds unconstrained.
};

/**
 * @brief Configuration options for MCS scoring.
 */
struct MCSOptions {
    MCSSearchMode search_mode = MCSSearchMode::Approximate;  ///< Which algorithm runs.
    MCSMatchLevel match_level = MCSMatchLevel::Default;      ///< Matching strictness.
    /**
     * How many matches one directed search may enumerate. The toolkit's own
     * default. Quality saturates well below it -- 256 already reproduced the
     * saturated bond count on all 45 pairs measured -- so this carries fourfold
     * headroom. Small values are destructive rather than merely faster: at 1,
     * 22 of those 45 pairs scored wrong. Zero is refused at construction.
     */
    unsigned int max_matches = 1024;
    /**
     * Return the bond Tanimoto itself rather than ``1 - Tanimoto``. Unlike the
     * RMSD comparison, MCS has a natural similarity form, so both orientations
     * are supported.
     */
    bool similarity = false;
};

/**
 * @brief Maximum-common-substructure comparison scored as Tanimoto over bonds.
 *
 * With ``c`` the matched-bond count and ``|A|``, ``|B|`` the two snapshots'
 * bond counts, the similarity is ``c / (|A| + |B| - c)`` and the distance is
 * one minus that.
 *
 * Hydrogens are suppressed in the constructor's snapshot, all three
 * ``OESuppressHydrogens`` retention flags explicitly false. The toolkit
 * defaults ``retainIsotope`` to true, which would leave a deuterium in place as
 * an explicit atom and put a labelled analogue on a different bond denominator
 * from its parent -- a difference a topological score has no way to mean.
 *
 * Suppression folds a hydrogen into the implicit hydrogen count of the atom it
 * hangs off, so a hydrogen stays explicit when that fold has nowhere to go:
 * when it has no single owner, being bonded to two atoms, or when it carries a
 * formal charge, which an implicit count cannot hold. Diborane keeps both
 * bridging hydrogens and is scored on its four B-H bonds, having no heavy-atom
 * bonds at all; ``[H-][Li+]`` keeps its hydride and is accepted where neutral
 * ``[H][Li]`` loses its only bond and is refused. Molecular hydrogen is the
 * degenerate case: one atom absorbs the other and survives as its owner, and
 * with no bonds left ``[H][H]`` is refused. Isotope is not a survival cause
 * once ``retainIsotope`` is false: perdeuterated benzene suppresses to the same
 * six-bond ring as benzene. Drug-like input has none of these, so the
 * denominator is a heavy-atom bond count in practice but not by construction.
 *
 * Approximate search is asymmetric: the match found with A as the pattern need
 * not equal the one found with B as the pattern. Since ``pdist`` fills only one
 * triangle, ``Compare`` runs both directions and takes the larger count, which
 * is symmetric by construction and never worse than either direction alone.
 *
 * Molecules with no bonds after hydrogen suppression -- methane, water, argon
 * -- are refused at construction, because bond Tanimoto is undefined rather
 * than merely extreme for them.
 *
 * ``Clone()`` deep-copies the molecule snapshots rather than aliasing the
 * parent's. ``pdist`` and ``cdist`` build every clone serially on the calling
 * thread before entering their parallel loop, so private per-clone molecules
 * mean no ``OEMol`` is reachable from two threads during the parallel phase.
 * That is a statement about one ``pdist``/``cdist`` invocation; sharing a
 * single ``MCSComparison`` across threads of the caller's own and calling
 * ``Compare`` on it directly bypasses ``Clone()`` and is unsupported, exactly
 * as it is for the other comparisons. The price is about 10.4 KB per molecule
 * per clone, so memory grows as ``num_threads x n``.
 *
 * There is no metric guarantee. No triangle-inequality violation appeared in
 * 74,400 ordered triples, but the inclusion-exclusion bound that would prove
 * the Jaccard metric property was violated 66 times over 59,280 triples, so the
 * proof is unavailable rather than merely unattempted.
 */
class MCSComparison : public PairwiseComparison {
public:
    using Options = MCSOptions;

    /**
     * @brief Construct an MCSComparison from molecules.
     *
     * Coordinates are irrelevant: the score is topological, so molecules parsed
     * from SMILES need no embedding step and a multi-conformer ``OEMol`` is
     * scored once rather than once per pose.
     *
     * :param mols: Shared pointers to molecules. Each is copied into the
     *     comparison and hydrogen-suppressed, so nothing the caller does to its
     *     own molecules afterwards can change a score.
     * :param opts: Scoring options.
     * :raises ComparisonError: When a pointer is null, a molecule has no bonds
     *     after hydrogen suppression, ``max_matches`` is zero, or
     *     ``search_mode`` or ``match_level`` carries a value outside its enum.
     */
    explicit MCSComparison(const std::vector<std::shared_ptr<OEChem::OEMol>>& mols,
                           const Options& opts = Options());

    /**
     * @brief Construct from molecules that need not be multiconformer.
     *
     * Molecules are snapshotted into ``OEMol`` internally, so an
     * ``OEGraphMol`` is accepted directly. The comparison never reads
     * coordinates, so nothing is lost by the narrower input type.
     *
     * An ``OEMol`` satisfies this overload's bindings typecheck as well as the
     * strict one above, so something has to rank them. That something is
     * typemap precedence, not declaration order: SWIG ranks by argument count
     * first and precedence second, and equal precedence leaves the order to an
     * unspecified tie-break that was measured to differ between argument-count
     * groups. See the precedence note on the ``OEMolBase*`` typecheck in
     * ``swig/oecluster.i``.
     *
     * The ranking has to be right even though this class cannot suffer a wrong
     * one. Binding here rather than above costs a multiconformer ``OEMol`` its
     * non-active conformers, which is unobservable to a comparison that never
     * reads coordinates -- hence "nothing is lost" above. The same typemap
     * serves ``ROCSComparison``, where that loss changes every score, so the
     * ordering is load-bearing there and merely correct here.
     *
     * :param mols: Pointers to molecules.
     * :param opts: MCS options.
     */
    explicit MCSComparison(const std::vector<OEChem::OEMolBase*>& mols,
                           const Options& opts = Options());

    ~MCSComparison() override;

    /**
     * @brief Score the maximum common substructure between items i and j.
     *
     * ``Compare(i, i)`` short-circuits to 0.0, or to 1.0 under ``similarity``,
     * without searching.
     *
     * :param i: Index of the first item.
     * :param j: Index of the second item.
     * :returns: The bond-Tanimoto distance, or the similarity under
     *     ``MCSOptions::similarity``.
     * :raises ComparisonError: When an index is past the end, or when the
     *     toolkit refuses the search construction, the ranking functor or the
     *     match budget.
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

    /// Private clone constructor -- takes the snapshot set it is handed. Unlike
    /// the other comparisons, ``Clone()`` hands it a fresh deep copy rather than
    /// the parent's own; see the class documentation.
    MCSComparison(std::shared_ptr<const SharedData> shared, const Options& opts);
};

}  // namespace OECluster

#endif  // OECLUSTER_COMPARISONS_MCSCOMPARISON_H
