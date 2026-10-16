/**
 * @file ROCSComparison.h
 * @brief Pairwise comparison based on OEShape overlay (ROCS-style).
 */

#ifndef OECLUSTER_COMPARISONS_ROCSCOMPARISON_H
#define OECLUSTER_COMPARISONS_ROCSCOMPARISON_H

#include <initializer_list>
#include <memory>
#include <string>
#include <vector>
#include "oecluster/PairwiseComparison.h"

namespace OEChem { class OEMol; }
namespace OEChem { class OEMolBase; }

namespace OECluster {

/**
 * @brief Score type for ROCS overlay distance computation.
 */
enum class ROCSScoreType {
    ComboNorm,  ///< TanimotoCombo normalized to [0,1]: distance = 1.0 - combo/2.0
    Combo,      ///< TanimotoCombo with range [0,2]: distance = 2.0 - combo
    Shape,      ///< Shape Tanimoto with range [0,1]: distance = 1.0 - shape
    Color       ///< Color Tanimoto with range [0,1]: distance = 1.0 - color
};

/**
 * @brief Configuration options for ROCS overlay scoring.
 */
struct ROCSOptions {
    ROCSScoreType score_type = ROCSScoreType::ComboNorm;  ///< Scoring method
    unsigned int color_ff_type = 1;  ///< OEColorFFType (1=ImplicitMillsDean)
    bool similarity = false;         ///< Return raw similarity instead of distance
};

/**
 * @brief ROCS-style shape/color pairwise comparison using OEShape.
 *
 * Copies the molecules it is given and uses ``OEOverlay`` to compute pairwise
 * overlay scores. The copy is not incidental: preparing the color atoms the
 * color term needs modifies the molecule, so the comparison must not work on
 * objects the caller still owns. Comparison output depends on score_type and
 * mode:
 *   - ComboNorm: ``1.0 - TanimotoCombo/2.0`` (range [0,1])
 *   - Combo: ``2.0 - TanimotoCombo`` (range [0,2])
 *   - Shape: ``1.0 - ShapeTanimoto`` (range [0,1])
 *   - Color: ``1.0 - ColorTanimoto`` (range [0,1])
 *
 * Every search builds and discards its own ``OEOverlay``, so clones share no
 * scorer state and Compare() can be called concurrently from different
 * threads without locking. A clone used to own one long-lived overlay
 * instead; see below for why it does not.
 *
 * ``Compare(i, j)`` is directional: it makes ``i`` the overlay reference and
 * searches only ``j``'s conformers, so the two orientations of a pair are two
 * different searches whose answers may differ. They need not: a self-pair, a
 * duplicate and any genuinely symmetric pair agree either way, which is why
 * ``CompareUnordered`` runs one search when ``i == j``. ``pdist`` fills one
 * triangle and presents it as symmetric, so it routes through a protected
 * seam, ``CompareUnordered``, which scores both orientations and keeps the
 * better overlay; so does every other reader of a symmetric distance that
 * accepts a ROCS comparison, which is why those readers and a ``pdist``
 * matrix agree. That seam is not part of the public surface and does not
 * appear in the bindings. ``cdist`` and direct ``Compare`` callers keep the directional
 * meaning, which is the one their output documents -- a ``cdist`` row says what
 * each reference finds in the fit set, and symmetrizing it would answer a
 * different question.
 *
 * Direction was only half of it. ``OEOverlay`` also carries state between
 * calls, enough to move an accepted input by 0.45 rather than by a rounding
 * error, so every search here -- ``CompareUnordered``'s two and ``Compare``'s
 * one -- runs on an overlay built for it and discarded after it. Neither
 * function reads or writes state that outlives the call, which is what makes
 * a score a function of its arguments. ``Compare`` reused the per-thread
 * overlay until 2026-10-16, and a ``cdist`` cell moved between 0.546 and a
 * saturated 1.0 with the thread count and with where an unrelated molecule
 * sat in the target list.
 *
 * The cost is two searches per off-diagonal pair instead of one, so a
 * ``pdist`` over n molecules runs n*(n-1) searches rather than n*(n-1)/2, and
 * every other reader of the seam pays the same doubling on the pairs it
 * reads. The n self-overlays the constructor takes are unchanged. Building an overlay per
 * search costs on top of that, and the fraction falls as the search grows,
 * since the construction is fixed and the search is not:
 * ``benchmarks/rocs_overlay_reuse.py`` measures the two strategies against
 * each other in one process, and on 2026-10-16 put fresh over shared at a
 * median of 1.20 at one conformer, 1.09 to 1.11 at twenty and 1.06 at fifty.
 *
 * What this buys is a seam value that is a function of its pair, and so a
 * ``pdist`` matrix that does not depend on input order and a set of lazy
 * readers that agree with it. Measured 2026-10-15: on three
 * flexible multiconformer molecules at twenty conformers each, the orientation
 * gap reaches 0.0355 on ComboNorm, 0.0710 on Combo, 0.0660 on Shape and 0.0010
 * on Color, and the forward-against-reversed residual is 0.0 on all four; on a
 * phenol displaced to 200 angstroms with benzene and octane, all six
 * permutations give a residual of 0.0, against 0.4538 before the fresh
 * overlay. Two fixtures, not a proof: the tests assert 1e-6 rather than
 * equality, because nothing in the toolkit promises bit-reproducibility.
 */
class ROCSComparison : public PairwiseComparison {
public:
    using Options = ROCSOptions;

    /**
     * @brief Construct a ROCSComparison from a set of molecules.
     *
     * Each molecule is copied, and the copies are shared across clones; the
     * caller's molecules are left untouched. Every molecule must carry 3D
     * coordinates, and the constructor refuses the set otherwise. The predicate
     * is OEChem's own dimension attribute, recomputed from the coordinates on
     * the comparison's copy before it is read, so a stale attribute --
     * ``SetCoords`` does not refresh one -- causes no refusal. That is an
     * axis count, not a geometric rank: a planar molecule rotated out of the
     * xy-plane counts three, and is admitted. Admitting it is correct, and the
     * reason this is the right guard rather than a rank test: across every input
     * class measured, the recomputed attribute coincides exactly with whether
     * OEShape can find coordinates to overlay at all. Genuinely linear molecules
     * such as N#N, O=C=O and C#N overlay well -- self shape distances of 0.0,
     * 0.0132 and 0.0128, the latter two no further off than the diatomic
     * chlorine that ROCSComparison.cpp already documents as not seating exactly
     * -- and a rank test would refuse all three. They are legitimate input in
     * the narrow sense that matters here, which is that this constructor should
     * admit them. It is not a claim about what happens afterwards: under the
     * default ComboNorm distance their diagonals are 0.500, 0.507 and 0.0067,
     * all nonzero, so all three stamp ``zero_self = No`` and the capability gate
     * refuses them with no override available.
     *
     * Construction also measures the diagonal once -- one self-overlay per
     * molecule, O(n) against the O(n^2) matrix this class exists to fill -- so
     * that Facts() can report ``zero_self`` from evidence.
     *
     * :param mols: Shared pointers to molecules.
     * :param opts: Scoring options.
     * :raises ComparisonError: Among the reasons: a molecule pointer is null;
     *     a molecule's recomputed dimension attribute is below three; a
     *     molecule's coordinates cannot be read; a coordinate is NaN or
     *     infinite; or a conformer spans more than 1e4 angstroms on one axis.
     */
    explicit ROCSComparison(const std::vector<std::shared_ptr<OEChem::OEMol>>& mols,
                        const Options& opts = Options());

    /**
     * @brief Construct from molecules that need not be multiconformer.
     *
     * Molecules are snapshotted into ``OEMol`` internally, so an
     * ``OEGraphMol`` is accepted directly. Such an input carries one
     * conformer, so ``BestOverlay`` has a single pose to choose from rather
     * than an ensemble; pass an ``OEMol`` when conformer search matters.
     *
     * C++ callers have to spell the conversion: an ``OEGraphMol*`` does not
     * convert to an ``OEMolBase*`` implicitly, and a pointer ``static_cast``
     * is rejected as well. Write ``&static_cast<OEChem::OEMolBase&>(graph_mol)``.
     * Callers arriving through the bindings never see this; the typemap
     * converts for them.
     *
     * An ``OEMol`` satisfies this overload's bindings typecheck as well as the
     * strict one above, so something has to rank them. That something is
     * typemap precedence, not declaration order: SWIG ranks by argument count
     * first and precedence second, and equal precedence leaves the order to an
     * unspecified tie-break that was measured to differ between argument-count
     * groups. See the precedence note on the ``OEMolBase*`` typecheck in
     * ``swig/oecluster.i``, which is what keeps a multiconformer ``OEMol`` off
     * this overload's active-conformer view -- a view that would silently
     * change every score it touched.
     *
     * :param mols: Pointers to molecules with 3D coordinates.
     * :param opts: ROCS options.
     */
    explicit ROCSComparison(const std::vector<OEChem::OEMolBase*>& mols,
                            const Options& opts = Options());

#ifndef SWIG
    /**
     * @brief Disambiguate a braced initializer between the two overloads above.
     *
     * ``ROCSComparison({})`` and ``ROCSComparison({nullptr})`` compiled before
     * the ``OEMolBase*`` overload existed and are ambiguous with it: both
     * vector types can be brace-initialized from either, so neither candidate
     * wins. An initializer-list constructor is considered ahead of both and
     * restores those spellings.
     *
     * It also takes ``{sp, sp}``, which used to reach the strict vector
     * constructor, because a braced list matches an initializer list exactly
     * where a vector needs a user-defined conversion. Delegating to the strict
     * constructor is therefore not a convenience but the requirement: routing
     * those molecules through the ``OEMolBase`` view instead would collapse
     * every multiconformer ``OEMol`` to its active pose and change the scores
     * such a call already returns.
     *
     * Hidden from SWIG deliberately. ``swig/oecluster.i`` ``%include``s this
     * header, so an unguarded third constructor would add a third candidate to
     * the overload dispatch whose ranking the precedence note in that file
     * exists to pin. This is a C++ source-compatibility measure with no
     * bindings surface.
     *
     * :param mols: Shared pointers to molecules, as for the strict overload.
     * :param opts: ROCS options.
     */
    explicit ROCSComparison(std::initializer_list<std::shared_ptr<OEChem::OEMol>> mols,
                            const Options& opts = Options());
#endif

    ~ROCSComparison() override;

    double Compare(size_t i, size_t j) override;
    std::unique_ptr<PairwiseComparison> Clone() const override;
    size_t Size() const override;
    std::string ComparisonName() const override;
    GateFacts Facts() const override;

protected:
    /// Score both orientations of {i, j} on overlays built for the purpose and
    /// keep the better one, so that the single triangle ``pdist`` stores
    /// depends neither on which index it passed first nor on what this object
    /// scored before. Equal indices run one search, both orientations there
    /// being the same search.
    ///
    /// Declared here rather than beside the public overrides above, and
    /// deliberately: the access of an override is set by its own declaration,
    /// so listing it there would make it public and give callers a second,
    /// contradictory way to compare two items. The base declares it protected
    /// precisely to keep that from happening.
    ///
    /// Both searches run on overlays built here and discarded here, so the two
    /// argument orders run the same two searches and have been measured to
    /// agree exactly. The tests still assert a tolerance rather than equality,
    /// that being a measurement of this toolkit on these fixtures rather than
    /// a promise the toolkit makes.
    double CompareUnordered(size_t i, size_t j) override;

private:
    struct SharedData;
    struct ThreadLocalData;
    std::shared_ptr<const SharedData> shared_;
    std::unique_ptr<ThreadLocalData> local_;
    Options opts_;

    /// Private clone constructor -- shares molecule data, and takes its own
    /// copy of the options each of its searches builds an overlay from.
    ROCSComparison(std::shared_ptr<const SharedData> shared,
               const Options& opts);

    /// Build the ``OEOverlayOptions`` every search constructs its own
    /// ``OEOverlay`` from; the name is left from when this held the overlay
    /// itself. On the three color-bearing score types, a non-null
    /// ``prep_target`` additionally has color atoms assigned to the molecules
    /// it holds; on ``Shape``, which never reads the color term, the target is
    /// ignored. Only the constructor that owns the snapshot passes a target;
    /// clones inherit an already-prepared set.
    void InitOverlay(SharedData* prep_target);

    /// Measure the diagonal by scoring every molecule against itself, writing
    /// both the ``zero_self`` and the ``data_integrity`` stamps into ``target``.
    /// Runs once, at construction, and the result is cached in ``SharedData`` so
    /// that clones inherit it: Facts() is const and called freely, so measuring
    /// there would put n overlays behind an accessor and repeat them for every
    /// clone.
    void MeasureDiagonal(SharedData& target);
};

}  // namespace OECluster

#endif  // OECLUSTER_COMPARISONS_ROCSCOMPARISON_H
