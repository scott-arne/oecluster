"""
Metric-capability gate for precomputed distance matrices.

The clustering algorithms in this package assume their input is a metric: an
item's distance to itself is zero, and the triangle inequality holds. Several
supported comparisons produce values satisfying neither. This module records
what a comparison actually guarantees (:func:`facts_from_comparison`) and
refuses to run an algorithm whose assumptions those guarantees violate
(:func:`require_metric`).

The gate has two tiers. Tier 1 -- distance orientation, zero self-distance,
and finite data -- is a hard refusal that no keyword overrides, because the
algorithms silently produce wrong output rather than failing. Finiteness is
the one tier-1 fact the gate measures off the data rather than reads off the
recorded stamp. Tier 2 -- the triangle inequality, subset-scored distances,
and proven probe violations -- is a soundness warning that
``allow_nonmetric=True`` overrides.
"""

import math

import numpy as np

from . import oecluster as _oecluster

_CAPABILITY = {
    _oecluster.Capability_Unknown: "unknown",
    _oecluster.Capability_No: False,
    _oecluster.Capability_Yes: True,
}

_DATA_INTEGRITY = {
    _oecluster.DataIntegrity_Complete: "complete",
    _oecluster.DataIntegrity_NaNPresent: "nan_present",
    _oecluster.DataIntegrity_SubsetScored: "subset_scored",
}


def default_facts():
    """
    Build the facts dictionary for a matrix that made no claims.

    Every value is the "unknown" form, and no unknown fact refuses on its own.
    This is what a distance matrix written before 5.0.0 loads as: the file
    recorded no capability claim, so the gate must not invent one. The one
    thing these facts cannot buy is a pass on non-finite data, which the gate
    measures off the stored distances rather than reads off the stamp.

    :returns: Fresh facts dictionary.
    """
    return {
        'is_distance': "unknown",
        'zero_self': "unknown",
        'triangle': "unknown",
        'data_integrity': "unknown",
        'metric_probe': "not_run",
        'probe_violations': 0,
        'probe_sampled': 0,
    }


def facts_from_comparison(comparison_obj):
    """
    Read the capability facts a comparison object stamps about itself.

    Read this *after* the distances are computed: a comparison may only
    discover a non-finite value while scoring pairs.

    :param comparison_obj: A C++ comparison exposing ``Facts()``.
    :returns: Facts dictionary; all-unknown for an object without ``Facts()``.
    """
    facts = default_facts()
    getter = getattr(comparison_obj, "Facts", None)
    if getter is None:
        return facts

    native = getter()
    facts['is_distance'] = _CAPABILITY.get(native.is_distance, "unknown")
    facts['zero_self'] = _CAPABILITY.get(native.zero_self, "unknown")
    facts['triangle'] = _CAPABILITY.get(native.triangle, "unknown")
    facts['data_integrity'] = _DATA_INTEGRITY.get(
        native.data_integrity, "unknown")
    return facts


def _has_nonfinite(distance_matrix):
    """
    Measure whether a matrix currently holds a non-finite distance.

    This reads the data, not the recorded stamp, so it still answers for a
    matrix whose values changed after the stamp was taken.

    :param distance_matrix: A :class:`SymmetricDistanceMatrix`.
    :returns: True if any distance the algorithms would read is NaN or
        infinite.
    """
    # Sparse storage is scanned entry by entry, not through ``.condensed``:
    # for a sparse matrix ``.condensed`` is a densified copy cached on first
    # access, so it answers for a snapshot the algorithms never read, while
    # costing a Python loop over every pair and seeding a cache that then
    # silently diverges from the storage. The entry list is deliberately left
    # un-deduplicated -- ``ThresholdGraph`` tests every tuple ``Entries()``
    # returns against the threshold, superseded duplicates included, so the raw
    # list is exactly what the algorithms see, while ``Get`` reports only the
    # last write for a pair. Values still sitting in an unmerged ``Set`` buffer
    # are deliberately out of scope for the same reason: until ``Finalize``
    # folds them in, they are invisible to ``Entries``, to ``Get`` and to the
    # algorithms, so refusing on them would be over-refusal.
    if isinstance(distance_matrix.storage, _oecluster.SparseStorage):
        return any(not math.isfinite(value)
                   for _, _, value in distance_matrix.storage._entries())

    # ``isfinite`` rather than ``isnan``: infinities break the algorithms the
    # same way, the native stamp already escalates on any non-finite value, and
    # the refusal below says "non-finite entries".
    return not np.isfinite(distance_matrix.condensed).all()


def check_allow_nonmetric(allow_nonmetric):
    """
    Refuse an ``allow_nonmetric`` that is not a bool or numpy.bool_.

    :param allow_nonmetric: Caller value.
    :raises TypeError: If the value is not a bool or numpy.bool_.
    """
    # A malformed override is a call the caller must fix whatever the matrix
    # looks like, so it is rejected before any fact is read. Truthiness would
    # be the wrong rule here: allow_nonmetric="False" reads to a caller as
    # "off" while switching the tier-2 checks off. numpy.bool_ is permitted
    # because it coerces faithfully and is never silently reinterpreted.
    if not isinstance(allow_nonmetric, (bool, np.bool_)):
        raise TypeError(
            "allow_nonmetric must be True or False, "
            f"not {type(allow_nonmetric).__name__} "
            f"({allow_nonmetric!r}). A truthy value would "
            "silently disable a safety check.")


def require_metric(distance_matrix, caller, *, allow_nonmetric=False):
    """
    Refuse to run a metric-assuming algorithm on a non-metric matrix.

    A malformed ``allow_nonmetric`` is rejected before any of the checks
    below, which then run in this order, the first failure raising:

    1. ``is_distance is False`` -- hard refusal.
    2. ``zero_self is False`` -- hard refusal.
    3. ``data_integrity == "nan_present"``, or the data itself holds a
       non-finite entry -- hard refusal.
    4. ``triangle is False`` -- overridable.
    5. ``data_integrity == "subset_scored"`` -- overridable.
    6. ``metric_probe == "violations_found"`` -- overridable.

    A fact of ``"unknown"`` never refuses on its own. Check 3 is the only one
    that does not decide from the stamp alone: it also measures the stored
    distances, so a matrix stamped ``"unknown"`` -- or ``"complete"`` -- that
    actually holds a non-finite value is refused on the measurement.

    The first two checks are separate because they answer separate questions.
    ``is_distance`` is about orientation -- whether a large number means "far"
    or "close" -- and a similarity clustered as a distance inverts every
    result. ``zero_self`` is about the diagonal, and a measure that scores a
    thing as unlike itself breaks the algorithms differently. Neither implies
    the other, and methane demonstrates both directions: it carries no colour
    features, so its ROCS colour *similarity* scores 0.0 against itself, while
    its ROCS ``combo_norm`` *distance* sits at 0.5 on the diagonal.

    :param distance_matrix: The matrix about to be clustered.
    :param caller: Name of the calling entry point, used in the messages.
    :param allow_nonmetric: Proceed despite a tier-2 violation.
    :raises TypeError: If allow_nonmetric is not a bool or numpy.bool_.
    :raises ValueError: If a check refuses.
    """
    check_allow_nonmetric(allow_nonmetric)

    facts = distance_matrix.facts
    name = distance_matrix.comparison_name

    if facts['is_distance'] is False:
        raise ValueError(
            f"{caller} requires distances, but this matrix holds a similarity "
            f"({name}). Recompute with similarity=False, or choose a metric "
            f"with a distance form. This cannot be overridden: every "
            f"clustering algorithm here reads a small value as 'close'.")

    if facts['zero_self'] is False:
        raise ValueError(
            f"{caller} requires a zero self-distance, but the measure behind "
            f"this matrix ({name}) does not score an item as identical to "
            f"itself. This cannot be overridden: every clustering algorithm "
            f"here treats the diagonal as the closest any pair can be.")

    # The stamp describes the data as computed. The caller can still write
    # through ``.condensed``, or through ``storage.Set`` -- plus ``Finalize``
    # on sparse storage -- and ``from_file`` trusts a file's stamp outright, so
    # the one tier-1 fact that is cheap to re-measure is re-measured, and a
    # stamp that disagrees with the data loses. Testing the stamp first
    # short-circuits only a matrix already stamped ``nan_present``; every other
    # matrix, including the common ``complete`` one, pays for the scan.
    if (facts['data_integrity'] == "nan_present"
            or _has_nonfinite(distance_matrix)):
        # ``missing`` is a descriptor-only option. Offering it to anyone else
        # sends them into a second refusal -- "Unknown kwargs for <name>
        # comparison: ['missing']" -- over a keyword their comparison never
        # took. The population this scan exists for, per the comment above, is
        # the one least likely to have a recompute available at all, so the
        # remedy that always holds is the one always named.
        remedy = ("Recompute with missing='complete_case', or remove the "
                  "offending items."
                  if name == "descriptor"
                  else "Remove the offending items.")
        raise ValueError(
            f"{caller} requires a complete distance matrix, but this matrix "
            f"contains non-finite entries. {remedy} This cannot be "
            f"overridden: every comparison against NaN is false, so the "
            f"clusters would depend on iteration order.")

    if allow_nonmetric:
        return

    if facts['triangle'] is False:
        raise ValueError(
            f"the measure behind this matrix ({name}) violates the triangle "
            f"inequality; {caller} assumes a metric. "
            f"Pass allow_nonmetric=True to proceed anyway.")

    if facts['data_integrity'] == "subset_scored":
        raise ValueError(
            f"this matrix was scored on a per-pair subset of features "
            f"(missing='ignore'), so its distances are not mutually "
            f"comparable; {caller} assumes a metric. "
            f"Pass allow_nonmetric=True to proceed anyway.")

    if facts['metric_probe'] == "violations_found":
        raise ValueError(
            f"precomputed matrix has proven triangle inequality violations "
            f"({facts['probe_violations']} of {facts['probe_sampled']} "
            f"sampled triples); {caller} assumes a metric. "
            f"Pass allow_nonmetric=True to proceed anyway.")


def require_comparable(distance_matrix, caller):
    """
    Refuse a matrix whose distances cannot be ranked or thresholded.

    The SAR-coherence metrics do not assume a metric. They rank distances
    against one another -- a nearest neighbour, a band minimum -- and compare
    them against a threshold, and a triangle-inequality violation leaves both
    operations meaningful. So this gate keeps ``require_metric``'s three tier-1
    refusals (orientation, the diagonal, non-finite entries) and waives the
    tier-2 ones, with a single exception.

    That exception is ``data_integrity == "subset_scored"``. A matrix computed
    with ``missing='ignore'`` scores each pair on whatever features that pair
    happens to share, so two of its distances answer different questions and
    the smaller one is not necessarily the nearer pair. Ranking is exactly what
    these metrics do, so the refusal is not overridable: there is no reading of
    the result that would be correct.

    :param distance_matrix: The matrix about to be scored.
    :param caller: Name of the calling entry point, used in the messages.
    :raises ValueError: If a check refuses.
    """
    require_metric(distance_matrix, caller, allow_nonmetric=True)

    if distance_matrix.facts['data_integrity'] == "subset_scored":
        raise ValueError(
            f"this matrix was scored on a per-pair subset of features "
            f"(missing='ignore'), so its distances are not mutually "
            f"comparable; {caller} ranks and thresholds distances against one "
            f"another. Recompute with missing='complete_case'. This cannot be "
            f"overridden: a nearest neighbour picked out of incomparable "
            f"distances is not a nearest neighbour.")


def condensed_lookup(condensed, n, i, j):
    """
    Look up ``d(i, j)`` in a condensed distance array, vectorized.

    ``i`` and ``j`` must be elementwise distinct. This is a precondition, not
    a check: enforcing it would cost a comparison over the ~1e5 elements the
    probe passes, and the probe already filters. Violating it does not raise.
    A diagonal entry has no place in a condensed array, so the arithmetic
    yields an index that numpy resolves against some other pair -- ``i == j
    == 0`` yields ``-1``, which wraps to the last element -- and the caller
    receives an ordinary-looking distance belonging to a different pair.

    :param condensed: 1-D condensed distance array for ``n`` items.
    :param n: Number of items.
    :param i: Array of row indices.
    :param j: Array of column indices, elementwise distinct from ``i``.
    :returns: Array of distances.
    """
    lo = np.minimum(i, j).astype(np.int64)
    hi = np.maximum(i, j).astype(np.int64)
    return condensed[n * lo + hi - ((lo + 2) * (lo + 1)) // 2]


def probe_triangle(condensed, n, *, samples=100000, seed=0):
    """
    Sample triples looking for a proven triangle-inequality violation.

    The probe can only disprove: a violating triple proves the matrix is not a
    metric, while finding none proves nothing. The caller therefore records
    the outcome separately from the ``triangle`` capability, which stays
    unknown either way.

    The comparison carries a tolerance that grows with the longest side of the
    triple, so that floating-point noise in a genuine metric does not read as
    a violation.

    :param condensed: 1-D condensed distance array.
    :param n: Number of items.
    :param samples: Number of triples to draw; a non-positive count skips the
        probe.
    :param seed: Seed for the sampler, fixed so results are reproducible.
    :returns: Dict with ``metric_probe``, ``probe_violations``, and
              ``probe_sampled``. The two counts are of *distinct* inequalities,
              not of draws.
    """
    skipped = {'metric_probe': "not_run", 'probe_violations': 0,
               'probe_sampled': 0}
    if n < 3 or samples <= 0:
        return skipped

    rng = np.random.default_rng(seed)
    i = rng.integers(0, n, size=samples)
    j = rng.integers(0, n, size=samples)
    k = rng.integers(0, n, size=samples)
    distinct = (i != j) & (j != k) & (i != k)
    i, j, k = i[distinct], j[distinct], k[distinct]
    if i.size == 0:
        return skipped

    # Drawing with replacement means the same inequality arrives many times
    # over, and ``d(i, k) <= d(i, j) + d(j, k)`` is the same inequality with
    # the two ends swapped. Counting draws would report more triples than a
    # small matrix contains -- a 4-item matrix admits 12 of these, and the
    # default draw reported 6233 violations of 37646 where 2 of 12 exist --
    # and ``require_metric`` prints both counts as its evidence.
    lo = np.minimum(i, k)
    hi = np.maximum(i, k)
    # One integer code per triple rather than ``np.unique(rows, axis=0)``,
    # which lexsorts a structured view: measured 7.6 ms against 51 ms for the
    # same 62129 triples at n = 60. Overflowing int64 would take n > 2e6,
    # whose condensed array alone would be some 17 TB.
    _, keep = np.unique((lo.astype(np.int64) * n + j) * n + hi,
                        return_index=True)
    i, j, k = lo[keep], j[keep], hi[keep]

    d_ij = condensed_lookup(condensed, n, i, j)
    d_jk = condensed_lookup(condensed, n, j, k)
    d_ik = condensed_lookup(condensed, n, i, k)
    scale = np.maximum(1.0, np.maximum(d_ij, np.maximum(d_jk, d_ik)))
    violations = int(np.count_nonzero(d_ik > d_ij + d_jk + 1e-9 * scale))

    return {
        'metric_probe': ("violations_found" if violations
                         else "no_violations_found"),
        'probe_violations': violations,
        'probe_sampled': int(i.size),
    }
