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
 * Each Clone() creates a new ``OEOverlay`` instance so that Compare()
 * can be called concurrently from different threads without locking.
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

private:
    struct SharedData;
    struct ThreadLocalData;
    std::shared_ptr<const SharedData> shared_;
    std::unique_ptr<ThreadLocalData> local_;
    Options opts_;

    /// Private clone constructor -- shares molecule data, creates new overlay.
    ROCSComparison(std::shared_ptr<const SharedData> shared,
               const Options& opts);

    /// Initialize OEOverlay with configured options. On the three color-bearing
    /// score types, a non-null ``prep_target`` additionally has color atoms
    /// assigned to the molecules it holds; on ``Shape``, which never reads the
    /// color term, the target is ignored. Only the constructor that owns the
    /// snapshot passes a target; clones inherit an already-prepared set.
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
