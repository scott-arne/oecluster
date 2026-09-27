/**
 * @file MCSComparison.cpp
 * @brief Implementation of maximum-common-substructure comparison.
 */

#include "oecluster/comparisons/MCSComparison.h"

#include <algorithm>
#include <exception>
#include <string>
#include <utility>
#include <oechem.h>
#include "IndexRange.h"
#include "oecluster/Error.h"

namespace OECluster {

namespace {

/// The OEChem atom and bond expression bitmasks one match level selects.
struct MatchExpressions {
    unsigned int atom;
    unsigned int bond;
};

MatchExpressions expressions_for(MCSMatchLevel level) {
    switch (level) {
        case MCSMatchLevel::Exact:
            return MatchExpressions{OEChem::OEExprOpts::ExactAtoms,
                                    OEChem::OEExprOpts::ExactBonds};
        case MCSMatchLevel::Loose:
            // A bond expression of zero constrains nothing, so an aromatic ring
            // matches its saturated analogue. Atomic number alone on the atoms
            // is the rest of what "loose" means.
            return MatchExpressions{OEChem::OEExprOpts::AtomicNumber, 0u};
        case MCSMatchLevel::Default:
            break;
    }
    return MatchExpressions{OEChem::OEExprOpts::DefaultAtoms,
                            OEChem::OEExprOpts::DefaultBonds};
}

/// The toolkit search-type constant for a search mode. Spelled out rather than
/// defaulted, because ``OEMCSType::Default`` is exhaustive and this library's
/// default is not.
unsigned int mcs_type_for(MCSSearchMode mode) {
    return mode == MCSSearchMode::Exhaustive ? OEChem::OEMCSType::Exhaustive
                                             : OEChem::OEMCSType::Approximate;
}

/// Read every atom's and bond's ring and aromatic flags once.
///
/// These are reads, so they cannot alter perception. What they do is
/// materialise any lazily-computed cache while a single thread still owns the
/// molecule, which keeps the per-clone copies from having to do it later. This
/// toolkit exposes no ``OEFindRingAtomAndBond``, so a read traversal is the
/// available lever.
void warm_perception(const OEChem::OEMol& mol) {
    for (OESystem::OEIter<OEChem::OEAtomBase> atom = mol.GetAtoms(); atom; ++atom) {
        (void)atom->IsInRing();
        (void)atom->IsAromatic();
    }
    for (OESystem::OEIter<OEChem::OEBondBase> bond = mol.GetBonds(); bond; ++bond) {
        (void)bond->IsInRing();
        (void)bond->IsAromatic();
    }
}

/// The largest matched-bond count one directed search finds.
///
/// Not the first match's count: approximate search does not yield best-first,
/// so every match the iterator produces within the budget is examined. An empty
/// iterator means a genuine zero overlap, which is distinguishable from failure
/// only because the three toolkit statuses below are checked first -- each of
/// them fails by producing a plausible score rather than an error.
///
/// ``Match`` hands back a heap-allocated iterator pointer rather than a range,
/// so it cannot drive a range-``for``; the ``OEIter`` handle owns the pointer
/// and is the form the toolkit's own MaximumCommonSS example uses. ``umatch``
/// is pinned true so the budget is spent on structurally distinct matches
/// rather than on automorphic repeats of one; it was measured not to change the
/// score on any of 120 pairs.
unsigned int directed_bond_count(const OEChem::OEMol& pattern, const OEChem::OEMol& target,
                                 const MCSOptions& opts, size_t pattern_index,
                                 size_t target_index) {
    const std::string pair = "items " + std::to_string(pattern_index) + " (pattern) and " +
                             std::to_string(target_index) + " (target)";
    const MatchExpressions expr = expressions_for(opts.match_level);

    // Spec 5.2: a toolkit throw has to reach the caller as a ComparisonError
    // naming the pair, never as a bare SDK or standard exception, and never as a
    // fabricated "no match". The three status checks below already raise
    // ComparisonError, so the first handler rethrows those unchanged -- this
    // mirrors FingerprintComparison.cpp:322-328.
    try {
        OEChem::OEMCSSearch search(pattern, expr.atom, expr.bond, mcs_type_for(opts.search_mode));
        if (!search) {
            throw ComparisonError("OEMCSSearch could not be constructed for " + pair +
                                  "; an unusable search yields an empty iterator, which would "
                                  "otherwise be scored as zero overlap");
        }
        if (!search.SetMCSFunc(OEChem::OEMCSMaxBondsCompleteCycles(1.0))) {
            throw ComparisonError("OEMCSSearch::SetMCSFunc was rejected for " + pair +
                                  "; the toolkit default ranking functor would otherwise stand in "
                                  "for the one this comparison scores against");
        }
        if (!search.SetMaxMatches(opts.max_matches)) {
            throw ComparisonError("OEMCSSearch::SetMaxMatches was rejected for " + pair +
                                  "; the toolkit default budget would otherwise stand in for the "
                                  "requested one");
        }

        unsigned int best = 0;
        for (OESystem::OEIter<const OEChem::OEMatchBase> match =
                 search.Match(target, /*umatch=*/true);
             match; ++match) {
            best = std::max<unsigned int>(best, match->NumBonds());
        }
        return best;
    } catch (const ComparisonError&) {
        throw;
    } catch (const std::exception& exc) {
        throw ComparisonError("The MCS search failed for " + pair + ": " +
                              std::string(exc.what()));
    } catch (...) {
        // OEChem signals most failures through OEThrow rather than by throwing,
        // so nothing is known to land here; without it an unrecognised throw
        // would reach Python as SWIG's bare "Unknown C++ exception", naming
        // neither the comparison nor the pair.
        throw ComparisonError("The MCS search failed for " + pair +
                              " with an unrecognised toolkit exception");
    }
}

/// Snapshot raw molecule pointers so the shared_ptr<OEMol> constructor can own
/// the validation. Null checking happens here because the conversion below
/// dereferences before that constructor ever sees the input.
std::vector<std::shared_ptr<OEChem::OEMol>> to_oemol_snapshots(
    const std::vector<OEChem::OEMolBase*>& mols) {
    std::vector<std::shared_ptr<OEChem::OEMol>> out;
    out.reserve(mols.size());
    for (size_t i = 0; i < mols.size(); ++i) {
        if (!mols[i]) {
            throw ComparisonError("MCSComparison received null molecule pointer at index " +
                                  std::to_string(i));
        }
        out.push_back(std::make_shared<OEChem::OEMol>(*mols[i]));
    }
    return out;
}

}  // namespace

struct MCSComparison::SharedData {
    /// Hydrogen-suppressed copies of the caller's molecules.
    std::vector<std::shared_ptr<const OEChem::OEMol>> mols;
    /// Each snapshot's bond count, which supplies the Tanimoto denominator.
    std::vector<unsigned int> bonds;
};

std::vector<const void*>
MCSComparisonSnapshotAccess::SnapshotAddresses(const MCSComparison& cmp) {
    std::vector<const void*> out;
    out.reserve(cmp.shared_->mols.size());
    for (const auto& mol : cmp.shared_->mols) {
        out.push_back(static_cast<const void*>(mol.get()));
    }
    return out;
}

MCSComparison::~MCSComparison() = default;

MCSComparison::MCSComparison(const std::vector<std::shared_ptr<OEChem::OEMol>>& mols,
                             const Options& opts)
    : opts_(opts) {
    // The toolkit accepts a zero budget silently -- SetMaxMatches(0) returns
    // true and GetMaxMatches() then reads 0 -- and every pair would score as
    // completely dissimilar with no error anywhere. So it has to be caught here.
    if (opts.max_matches == 0) {
        throw ComparisonError(
            "MCSComparison requires max_matches >= 1; the toolkit reads zero as a budget of "
            "zero matches rather than as unlimited, so every pair would score distance 1");
    }
    // A bad integer arriving through the bindings must not reach OEMCSSearch as
    // an undefined search type or expression pair.
    const int search_mode_value = static_cast<int>(opts.search_mode);
    if (search_mode_value < static_cast<int>(MCSSearchMode::Approximate) ||
        search_mode_value > static_cast<int>(MCSSearchMode::Exhaustive)) {
        throw ComparisonError("MCSComparison received an out-of-range search_mode value: " +
                              std::to_string(search_mode_value));
    }
    const int match_level_value = static_cast<int>(opts.match_level);
    if (match_level_value < static_cast<int>(MCSMatchLevel::Default) ||
        match_level_value > static_cast<int>(MCSMatchLevel::Loose)) {
        throw ComparisonError("MCSComparison received an out-of-range match_level value: " +
                              std::to_string(match_level_value));
    }

    for (size_t i = 0; i < mols.size(); ++i) {
        if (!mols[i]) {
            throw ComparisonError("MCSComparison received null molecule pointer at index " +
                                  std::to_string(i));
        }
    }

    auto shared = std::make_shared<SharedData>();
    shared->mols.reserve(mols.size());
    shared->bonds.reserve(mols.size());
    for (size_t i = 0; i < mols.size(); ++i) {
        auto snapshot = std::make_shared<OEChem::OEMol>(*mols[i]);
        // All three retention flags explicit. The toolkit defaults
        // retainIsotope to true, which leaves a deuterium in place as an
        // explicit atom: benzene-d1 would carry 7 bonds against benzene's 6 and
        // score 0.857 similarity against its own parent.
        OEChem::OESuppressHydrogens(*snapshot, false, false, false);

        const unsigned int bond_count = snapshot->NumBonds();
        if (bond_count == 0) {
            throw ComparisonError(
                "MCSComparison received molecule at index " + std::to_string(i) + " ('" +
                std::string(snapshot->GetTitle()) +
                "') with no bonds after hydrogen suppression; bond Tanimoto has a zero "
                "denominator there, and calling two bondless molecules identical would "
                "cluster argon with water");
        }
        warm_perception(*snapshot);

        shared->mols.push_back(std::shared_ptr<const OEChem::OEMol>(std::move(snapshot)));
        shared->bonds.push_back(bond_count);
    }

    shared_ = std::move(shared);
}

// Delegates rather than duplicating the four validations and the
// snapshot-and-suppress loop above, so neither can be changed for one path and
// forgotten for the other. The two are not quite interchangeable: this path
// null-checks while snapshotting, before the delegate validates the options, so
// given both a null pointer and a bad option it reports the null pointer where
// the direct path reports the option. Unreachable from Python, where the
// bindings reject a null list element before any constructor runs, and harmless
// in C++, where both inputs are errors and either message names a real one.
//
// The resulting double copy -- raw pointer to OEMol here, then the delegate's
// own snapshot -- is O(n) against the O(n^2) matrix this class exists to fill.
MCSComparison::MCSComparison(const std::vector<OEChem::OEMolBase*>& mols,
                             const Options& opts)
    : MCSComparison(to_oemol_snapshots(mols), opts) {}

// Delegates to the strict constructor, which is the whole point of it: a
// braced list of shared_ptr used to land there directly and must still end up
// there. See the header for why this overload exists at all.
MCSComparison::MCSComparison(std::initializer_list<std::shared_ptr<OEChem::OEMol>> mols,
                             const Options& opts)
    : MCSComparison(std::vector<std::shared_ptr<OEChem::OEMol>>(mols), opts) {}

MCSComparison::MCSComparison(std::shared_ptr<const SharedData> shared, const Options& opts)
    : shared_(std::move(shared)), opts_(opts) {}

double MCSComparison::Compare(size_t i, size_t j) {
    detail::check_compare_index_range("MCSComparison", i, j, shared_->mols.size());

    // The diagonal is answered without searching. Self-distances measured
    // exactly zero anyway, but zero_self is a tier-1 gate fact in distance mode,
    // and one approximate search under-matching a large molecule against itself
    // would break clustering outright.
    if (i == j) {
        return opts_.similarity ? 1.0 : 0.0;
    }

    // Both directions, larger count wins: approximate search is asymmetric and
    // pdist fills only one triangle, so a one-direction scorer would make the
    // matrix depend on input order.
    const unsigned int matched =
        std::max(directed_bond_count(*shared_->mols[i], *shared_->mols[j], opts_, i, j),
                 directed_bond_count(*shared_->mols[j], *shared_->mols[i], opts_, j, i));

    // A common substructure cannot have more bonds than either molecule, so the
    // denominator is at least max(|A|, |B|) and the constructor's zero-bond
    // refusal makes that at least 1.
    const unsigned int denominator = shared_->bonds[i] + shared_->bonds[j] - matched;
    const double similarity = static_cast<double>(matched) / static_cast<double>(denominator);
    return opts_.similarity ? similarity : 1.0 - similarity;
}

std::unique_ptr<PairwiseComparison> MCSComparison::Clone() const {
    // Deep copy, deliberately unlike the other comparisons, which alias the
    // parent's snapshots. pdist and cdist build every clone serially on the
    // calling thread before entering ParallelFor, so giving each clone private
    // molecules means no OEMol is reachable from two threads during the parallel
    // phase -- not a narrower window, but none. Neither re-suppression nor
    // re-validation is needed: the parent already did both, and the bond counts
    // are copied verbatim.
    auto copy = std::make_shared<SharedData>();
    copy->mols.reserve(shared_->mols.size());
    for (const auto& mol : shared_->mols) {
        copy->mols.push_back(std::make_shared<const OEChem::OEMol>(*mol));
    }
    copy->bonds = shared_->bonds;
    return std::unique_ptr<PairwiseComparison>(new MCSComparison(std::move(copy), opts_));
}

size_t MCSComparison::Size() const {
    return shared_->mols.size();
}

std::string MCSComparison::ComparisonName() const {
    return "mcs";
}

GateFacts MCSComparison::Facts() const {
    GateFacts facts;
    facts.is_distance = opts_.similarity ? Capability::No : Capability::Yes;
    // zero_self tracks similarity for the same reason is_distance does: it
    // asserts that Compare(i, i) is zero, and under similarity it is 1.0.
    // Reporting Yes would be false metadata even though the gate refuses on
    // is_distance first, because GateFacts documents the two as independent.
    facts.zero_self = opts_.similarity ? Capability::No : Capability::Yes;
    // Unknown states exactly what is known: no proof and no counterexample.
    // Bond Tanimoto would be a metric if the matched-bond count were a
    // set-intersection cardinality, but it is not -- the inclusion-exclusion
    // bound genuine intersections satisfy was violated 66 times over 59,280
    // triples, so the Jaccard proof is unavailable. Against that, no actual
    // triangle violation appeared in 74,400 ordered triples. Unknown is
    // permissive at the gate, so it costs no caller anything.
    facts.triangle = Capability::Unknown;
    // Every pair is scored by one formula over one quantity, so nothing here is
    // SubsetScored: that name is for a score computed over a per-pair subset of
    // the available dimensions, which is the descriptor comparison's
    // missing='ignore' and not this.
    facts.data_integrity = DataIntegrity::Complete;
    return facts;
}

}  // namespace OECluster
