"""
Metric-capability gate for precomputed distance matrices.

The clustering algorithms in this package assume their input is a metric: an
item's distance to itself is zero, and the triangle inequality holds. Several
supported comparisons produce values satisfying neither. This module records
what a comparison actually guarantees (:func:`facts_from_comparison`) and
refuses to run an algorithm whose assumptions those guarantees violate
(:func:`require_metric`).

The gate has two tiers. Tier 1 -- distance orientation, zero self-distance,
and NaN-free data -- is a hard refusal that no keyword overrides, because the
algorithms silently produce wrong output rather than failing. Tier 2 -- the
triangle inequality, subset-scored distances, and proven probe violations --
is a soundness warning that ``allow_nonmetric=True`` overrides.
"""

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

    Every value is the "unknown" form, which the gate never refuses. This is
    what a distance matrix written before 5.0.0 loads as: the file recorded no
    capability claim, so the gate must not invent one.

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


def require_metric(distance_matrix, caller, *, allow_nonmetric=False):
    """
    Refuse to run a metric-assuming algorithm on a non-metric matrix.

    A malformed ``allow_nonmetric`` is rejected before any of the checks
    below, which then run in this order, the first failure raising:

    1. ``is_distance is False`` -- hard refusal.
    2. ``zero_self is False`` -- hard refusal.
    3. ``data_integrity == "nan_present"`` -- hard refusal.
    4. ``triangle is False`` -- overridable.
    5. ``data_integrity == "subset_scored"`` -- overridable.
    6. ``metric_probe == "violations_found"`` -- overridable.

    A fact of ``"unknown"`` never refuses.

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

    if facts['data_integrity'] == "nan_present":
        raise ValueError(
            f"{caller} requires a complete distance matrix, but this matrix "
            f"contains non-finite entries. Recompute with "
            f"missing='complete_case', or remove the offending items. This "
            f"cannot be overridden: every comparison against NaN is false, "
            f"so the clusters would depend on iteration order.")

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
