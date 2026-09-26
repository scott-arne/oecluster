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
 */
enum class MCSMatchLevel {
    Default,  ///< OEChem's default atom and bond expressions.
    Exact,    ///< Adds hydrogen count, charge, degree and bond order.
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
 * With ``c`` the matched-bond count and ``|A|``, ``|B|`` the two molecules'
 * heavy-atom bond counts, the similarity is ``c / (|A| + |B| - c)`` and the
 * distance is one minus that.
 *
 * Hydrogens are always suppressed in the constructor's snapshot, all three
 * ``OESuppressHydrogens`` retention flags explicitly false. The toolkit
 * defaults ``retainIsotope`` to true, which would leave a deuterium in place as
 * an explicit atom and put a labelled analogue on a different bond denominator
 * from its parent -- a difference a topological score has no way to mean.
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
