"""
Parameter selection: sweep one clustering keyword and pick a value by a
validity index.

The clustering entry points and the two scorecards (:func:`cluster_report`,
:func:`isim_report`) already exist; this module only orchestrates them. It
runs one algorithm over an explicit grid of one keyword, scores every
partition, and ranks the rows under a named criterion. :class:`ClusteringSpec`
is the reusable description of "which algorithm, with which fixed options"
that the sweep runs; it is deliberately unbound from the input so that later
workflow features can rerun it on a resample or combine several into an
ensemble.

Everything from the package is resolved at call time rather than imported at
the top: the package ``__init__`` imports this module before the clustering
functions and scorecards it needs are defined there.
"""

import math
import numbers
import types
from collections.abc import Callable, Mapping
from typing import Any, NamedTuple

import numpy as np

from . import oecluster as _oecluster

# Every public entry point that returns a ClusteringResult, by its public
# name. knn_graph returns a graph and is left out on purpose.
_ROSTER_NAMES = (
    "butina", "dbscan", "hdbscan", "agglomerative", "k_medoids", "bitbirch",
    "bitbirch_recluster", "bitbirch_refine", "sphere_exclusion",
    "jarvis_patrick", "leiden", "murcko",
)

_ROSTER: dict[str, Callable[..., Any]] | None = None


def _package():
    """The initialized package, looked up at call time (see module docstring)."""
    import oecluster
    return oecluster


def _roster():
    """Public name -> clustering function, built on first use."""
    global _ROSTER
    roster = _ROSTER
    if roster is None:
        package = _package()
        roster = {name: getattr(package, name) for name in _ROSTER_NAMES}
        _ROSTER = roster
    return roster


def _options_equal(left, right):
    """Compare option mappings, handling NumPy arrays safely.

    NumPy's == operator is elementwise and returns a boolean array for
    multi-element arrays, which cannot be converted to a single truth value.
    This helper compares values element-wise for arrays and by equality
    for other types.

    :param left: First options mapping.
    :param right: Second options mapping.
    :returns: True if both mappings have the same keys and equal values.
    """
    if set(left) != set(right):
        return False
    for key in left:
        left_val = left[key]
        right_val = right[key]
        # Identity first, as Python's own container equality does; a NaN
        # option would otherwise make a spec unequal to itself.
        if left_val is right_val:
            continue
        # Use np.array_equal for NumPy arrays to avoid the ambiguous truth value error.
        if isinstance(left_val, np.ndarray) or isinstance(right_val, np.ndarray):
            if not np.array_equal(left_val, right_val):
                return False
        else:
            if left_val != right_val:
                return False
    return True


def _roster_name(function):
    """The public name of a roster function, else the callable's own name."""
    for name, candidate in _roster().items():
        if candidate is function:
            return name
    return getattr(function, "__name__", type(function).__name__)


class ClusteringSpec:
    """A clustering entry point plus the keyword options it is run with.

    The spec is unbound from the input: :meth:`run` takes the items each
    time, so one spec can be rerun on another input unchanged. A roster
    function is recorded under its public name, which is what a later
    serialization needs; any other callable that returns a
    :class:`ClusteringResult` is accepted under its own ``__name__``.

    The options are not validated at construction; the algorithm validates
    them when it runs, exactly as a direct call would.

    Specs compare equal on the same algorithm object and equal options.
    Defining ``__eq__`` without ``__hash__`` makes them unhashable, which is
    intended: the options are compared by value.
    """

    def __init__(self, algorithm, **options):
        """
        :param algorithm: A roster name (case-insensitive), the matching
            function object, or any callable returning a
            :class:`ClusteringResult`.
        :param options: Keyword options fixed for every run.
        :raises TypeError: If ``algorithm`` is neither a string nor callable.
        :raises ValueError: If the name is not in the roster.
        """
        if isinstance(algorithm, str):
            key = algorithm.lower()
            roster = _roster()
            if key not in roster:
                raise ValueError(
                    f"Unknown clustering algorithm {algorithm!r}; choose one "
                    f"of {', '.join(_ROSTER_NAMES)}")
            self._algorithm = roster[key]
            self._name = key
        elif callable(algorithm):
            self._algorithm = algorithm
            self._name = _roster_name(algorithm)
        else:
            raise TypeError(
                "ClusteringSpec() expects a roster name or a callable, not "
                f"{type(algorithm).__name__}")
        self._options = types.MappingProxyType(dict(options))

    @property
    def name(self):
        """The public roster name, or the callable's ``__name__``."""
        return self._name

    @property
    def algorithm(self):
        """The callable that :meth:`run` invokes."""
        return self._algorithm

    @property
    def options(self):
        """Read-only view of the fixed keyword options."""
        return self._options

    def run(self, items, **overrides):
        """Cluster ``items`` with the fixed options, overrides winning.

        :param items: Passed to the algorithm untouched.
        :param overrides: Keyword options that replace fixed ones for this run.
        :returns: The algorithm's :class:`ClusteringResult`.
        :raises TypeError: If the callable returned something else.
        """
        options = {**self._options, **overrides}
        result = self._algorithm(items, **options)
        if not isinstance(result, _package().ClusteringResult):
            raise TypeError(
                f"{self!r} returned {type(result).__name__}, not a "
                "ClusteringResult")
        return result

    def __eq__(self, other):
        if not isinstance(other, ClusteringSpec):
            return NotImplemented
        return (self._algorithm is other._algorithm
                and _options_equal(self._options, other._options))

    def __repr__(self):
        parts = [repr(self._name)]
        parts.extend(f"{key}={value!r}" for key, value in self._options.items())
        return f"ClusteringSpec({', '.join(parts)})"


# Validity indices by direction: +1 maximize, -1 minimize. This is the only
# place a direction is recorded; a test checks every key against the
# scorecards' field lists. Profile metrics (counts, fractions, size
# statistics, intra distances, radius, diameter) are left out because they
# have no direction of their own.
_CRITERIA: dict[str, int] = {
    "silhouette": 1,
    "isim_silhouette": 1,
    "dunn_index": 1,
    "dunn_mean_separation_mean_diameter": 1,
    "dunn_medoid_separation_medoid_spread": 1,
    "calinski_harabasz_medoid": 1,
    "point_biserial": 1,
    "baker_hubert_gamma": 1,
    "davies_bouldin_medoid": -1,
    "c_index": -1,
}

_DEFAULT_CRITERION = {
    "cluster_report": "silhouette",
    "isim_report": "isim_silhouette",
}

# Criteria that exist only when the scorer runs an optional stage. The sweep
# switches the stage on for its criterion; a caller who switched it off and
# asked for that criterion gets a refusal rather than a table of NaN.
_STAGE_FLAGS = {
    "cluster_report": {
        "c_index": "compute_pair_rank_indices",
        "baker_hubert_gamma": "compute_pair_rank_indices",
    },
    "isim_report": {
        "isim_silhouette": "compute_centroid_indices",
        "davies_bouldin_medoid": "compute_centroid_indices",
        "dunn_medoid_separation_medoid_spread": "compute_centroid_indices",
    },
}

# The pair-rank indices sort every pairwise distance, which cluster_report
# refuses to build from a comparison, so they need a matrix input.
_MATRIX_ONLY = frozenset({"c_index", "baker_hubert_gamma"})


class _Scorer(NamedTuple):
    """The scorecard a sweep scores with, chosen from the input kind."""

    name: str
    function: Callable[..., Any]
    report_class: Any
    matrix: bool


def _is_fingerprint_batch(items):
    try:
        import oefp
    except ImportError:
        return False
    return isinstance(items, oefp.OEFPBatch)


def _scorer_for(items):
    """Pick the scorer for ``items``, refusing kinds no scorer accepts.

    A raw item sequence is refused rather than accepted with a
    ``comparison=`` option: the sweep would then name the comparison twice,
    once for the algorithm and once for the scorer, with nothing keeping the
    two the same. A prebuilt comparison is the single definition both read.

    :raises TypeError: For a kind outside the three accepted ones.
    :raises ValueError: For sparse storage, which cluster_report refuses;
        an algorithm such as leiden would otherwise cluster it and the sweep
        would fail only at scoring.
    """
    package = _package()
    if _is_fingerprint_batch(items):
        return _Scorer("isim_report", package.isim_report, package.ISimReport,
                       False)
    if isinstance(items, package.CrossDistanceMatrix):
        raise TypeError(
            "select_parameter() requires a SymmetricDistanceMatrix; a "
            "CrossDistanceMatrix is rectangular and no clustering algorithm "
            "accepts it")
    if isinstance(items, package.SymmetricDistanceMatrix):
        # ValueError, not TypeError: the argument's type is right, its
        # storage is not.
        if isinstance(items.storage, package.SparseStorage):
            raise ValueError(  # noqa: TRY004
                "select_parameter() requires complete pairwise distances to "
                "score every partition; SparseStorage is not supported by "
                "cluster_report")
        return _Scorer("cluster_report", package.cluster_report,
                       package.ClusterReport, True)
    if isinstance(items, _oecluster.PairwiseComparison):
        return _Scorer("cluster_report", package.cluster_report,
                       package.ClusterReport, False)
    raise TypeError(
        "select_parameter() expects a SymmetricDistanceMatrix, a prebuilt "
        "comparison such as FingerprintComparison(mols), or an oefp.OEFPBatch, "
        f"not {type(items).__name__}; wrap a raw item sequence in a prebuilt "
        "comparison so that clustering and scoring read the same distances")


def _resolve_criterion(criterion, scorer, report_options):
    """The criterion to rank by, defaulted per scorer and checked.

    Merges the opt-in stage flag the criterion needs into ``report_options``
    (the sweep's own copy, never the caller's mapping).

    :raises TypeError: If ``criterion`` is not a string, or a stage flag the
        caller put in ``report_options`` is not a bool.
    :raises ValueError: If it is not a validity index, the scorer does not
        produce it, it needs a matrix and the input is a comparison, or the
        caller switched its stage off.
    """
    if criterion is None:
        criterion = _DEFAULT_CRITERION[scorer.name]
    elif not isinstance(criterion, str):
        raise TypeError(
            "criterion must be a str naming a report field, not "
            f"{type(criterion).__name__}")
    if criterion not in _CRITERIA:
        raise ValueError(
            f"criterion {criterion!r} is not a validity index; choose one of "
            f"{', '.join(_CRITERIA)}")
    produced = scorer.report_class._SCALAR_FIELDS
    if criterion not in produced:
        available = ", ".join(name for name in _CRITERIA if name in produced)
        raise ValueError(
            f"criterion {criterion!r} is not produced by {scorer.name}; its "
            f"validity indices are {available}")
    if criterion in _MATRIX_ONLY and not scorer.matrix:
        raise ValueError(
            f"criterion {criterion!r} needs every pairwise distance, which "
            "cluster_report does not build from a comparison; pass a "
            "SymmetricDistanceMatrix")
    # The scorers refuse truthiness for their flags; checking here keeps a
    # non-bool from surfacing only after the first clustering has run.
    for known_flag in set(_STAGE_FLAGS[scorer.name].values()):
        if known_flag in report_options and not isinstance(
                report_options[known_flag], (bool, np.bool_)):
            raise TypeError(
                f"report_options[{known_flag!r}] must be True or False, not "
                f"{type(report_options[known_flag]).__name__}")
    flag = _STAGE_FLAGS[scorer.name].get(criterion)
    if flag is not None:
        if flag in report_options:
            if not report_options[flag]:
                raise ValueError(
                    f"criterion {criterion!r} requires "
                    f"report_options[{flag!r}]=True")
        else:
            report_options[flag] = True
    return criterion


def _validate_bounds(max_noise_fraction, min_clusters, max_clusters):
    """Check the three bounds; returns them normalized to plain Python types.

    ``bool`` is refused for every bound, as the scorecards refuse truthiness
    for their flags: ``True`` as a noise cap is a mistake, not a value of 1.
    Normalizing a NumPy scalar bound to plain ``float``/``int`` here, rather
    than leaving it as ``np.float64``/``np.int64``, keeps the rejection text
    :func:`_first_violation` builds free of the ``np.float64(...)`` repr.

    :raises TypeError: For a bound of the wrong type.
    :raises ValueError: For a bound out of range or an inverted pair.
    :returns: ``(max_noise_fraction, min_clusters, max_clusters)``, each
        None or its normalized type.
    """
    if max_noise_fraction is not None:
        if (isinstance(max_noise_fraction, bool)
                or not isinstance(max_noise_fraction, numbers.Real)):
            raise TypeError("max_noise_fraction must be a real number or None")
        max_noise_fraction = float(max_noise_fraction)
        if math.isnan(max_noise_fraction) or not 0.0 <= max_noise_fraction <= 1.0:
            raise ValueError("max_noise_fraction must be between 0 and 1")
    bounds = {"min_clusters": min_clusters, "max_clusters": max_clusters}
    for name, bound in bounds.items():
        if bound is None:
            continue
        if isinstance(bound, bool) or not isinstance(bound, numbers.Integral):
            raise TypeError(f"{name} must be an int or None")
        if bound < 1:
            raise ValueError(f"{name} must be at least 1")
        bounds[name] = int(bound)
    min_clusters, max_clusters = bounds["min_clusters"], bounds["max_clusters"]
    if (min_clusters is not None and max_clusters is not None
            and min_clusters > max_clusters):
        raise ValueError("min_clusters must not exceed max_clusters")
    return max_noise_fraction, min_clusters, max_clusters


def _first_violation(report, max_noise_fraction, min_clusters, max_clusters):
    """The first bound ``report`` violates, as text, or None."""
    # ``not (a <= b)`` rather than ``a > b`` so a NaN noise fraction counts
    # as a violation instead of slipping through.
    if (max_noise_fraction is not None
            and not report.noise_fraction <= max_noise_fraction):
        return (f"noise_fraction {report.noise_fraction!r} > "
                f"max_noise_fraction {max_noise_fraction!r}")
    if min_clusters is not None and report.num_clusters < min_clusters:
        return f"num_clusters {report.num_clusters!r} < min_clusters {min_clusters!r}"
    if max_clusters is not None and report.num_clusters > max_clusters:
        return f"num_clusters {report.num_clusters!r} > max_clusters {max_clusters!r}"
    return None


class SweepRow(NamedTuple):
    """One grid value of a sweep: its partition, its report and its rank input.

    ``score`` is the criterion read off ``report`` as a float, NaN when the
    scorecard left it undefined. ``eligible`` is False when a bound was
    violated, and ``rejection`` then names the first violated bound.
    """

    value: Any
    result: Any
    report: Any
    score: float
    eligible: bool
    rejection: str | None


class ParameterSelection:
    """Read-only outcome of :func:`select_parameter`.

    ``rows`` holds every grid value in grid order with its result and report,
    so the table can be re-ranked under any other field those reports
    already carry without rerunning anything. ``winner`` is None when no
    eligible row has a non-NaN score; the table is still here to show why.

    Assignment to any attribute, public or private, raises
    ``AttributeError``: the fields are set through ``object.__setattr__`` in
    the constructor and ``__setattr__`` rejects everything after, the same
    arrangement :class:`ClusterReport` uses.
    """

    def __init__(self, spec, parameter, criterion, rows, winner_index):
        object.__setattr__(self, "_spec", spec)
        object.__setattr__(self, "_parameter", parameter)
        object.__setattr__(self, "_criterion", criterion)
        object.__setattr__(self, "_rows", tuple(rows))
        object.__setattr__(self, "_winner_index", winner_index)

    def __setattr__(self, name, value):
        """ParameterSelection is read-only; reject attribute assignment."""
        raise AttributeError(
            f"ParameterSelection is read-only; cannot set {name!r}")

    @property
    def spec(self):
        """The :class:`ClusteringSpec` that was run."""
        return self._spec

    @property
    def parameter(self):
        """The keyword that was swept."""
        return self._parameter

    @property
    def criterion(self):
        """The report field the rows were ranked by."""
        return self._criterion

    @property
    def rows(self):
        """Tuple of :class:`SweepRow`, one per grid value, in grid order."""
        return self._rows

    @property
    def winner_index(self):
        """Index of the winning row into :attr:`rows`, or None."""
        return self._winner_index

    @property
    def winner(self):
        """The winning :class:`SweepRow`, or None."""
        if self._winner_index is None:
            return None
        return self._rows[self._winner_index]

    @property
    def columns(self):
        """Column names of :meth:`to_table`."""
        return (self._parameter, self._criterion, "num_clusters",
                "noise_fraction", "eligible", "rejection")

    def to_table(self):
        """Return a list of ``(value, score, num_clusters, noise_fraction,
        eligible, rejection)`` tuples, one per row, matching :attr:`columns`.

        A fresh list each call, the ``compare_reports().to_table()``
        convention of plain tuples a caller can hand to a table printer.
        """
        return [(row.value, row.score, row.report.num_clusters,
                 row.report.noise_fraction, row.eligible, row.rejection)
                for row in self._rows]

    def __repr__(self):
        header = (f"ParameterSelection(spec={self._spec!r}, "
                  f"parameter={self._parameter!r}, "
                  f"criterion={self._criterion!r})")
        labels = self.columns[:4]
        # A grid built with np.linspace holds np.float64 values, whose own
        # repr is "np.float64(0.05)"; .item() unwraps the plain Python
        # scalar for display without changing SweepRow.value itself.
        cells = [(repr(row.value.item() if isinstance(row.value, np.generic)
                      else row.value),
                 f"{row.score:.4f}", str(row.report.num_clusters),
                 f"{row.report.noise_fraction:.4g}")
                 for row in self._rows]
        widths = [max(len(label), *(len(cell[i]) for cell in cells))
                  for i, label in enumerate(labels)]
        lines = [header,
                 "   " + "  ".join(f"{label:<{width}}"
                                   for label, width in zip(labels, widths))]
        for index, (row, cell) in enumerate(zip(self._rows, cells)):
            marker = " * " if index == self._winner_index else "   "
            line = marker + "  ".join(f"{text:<{width}}"
                                      for text, width in zip(cell, widths))
            if not row.eligible:
                line += f"  rejected: {row.rejection}"
            lines.append(line.rstrip())
        return "\n".join(lines)


def _winner_index(rows, direction):
    """Index of the best eligible non-NaN score, earliest on ties, or None.

    Infinities take part: Davies-Bouldin reports +inf for coincident
    medoids, which is a defined worst value under minimization, not a
    missing one. Only NaN is excluded.
    """
    best = None
    best_key = None
    for index, row in enumerate(rows):
        if not row.eligible or math.isnan(row.score):
            continue
        key = direction * row.score
        if best_key is None or key > best_key:
            best = index
            best_key = key
    return best


def select_parameter(algorithm, items, parameter, values, *, criterion=None,
                     max_noise_fraction=None, min_clusters=None,
                     max_clusters=None, report_options=None):
    """Sweep one clustering keyword over a grid and pick a value by a
    validity index.

    Every grid value is clustered with the spec, scored with
    :func:`cluster_report` (a matrix or a prebuilt comparison) or
    :func:`isim_report` (a fingerprint batch), and kept as a
    :class:`SweepRow`. Argument validation runs before the first clustering,
    so a bad call fails at once rather than after a long run. Grid points run
    one after another in the caller's order.

    :param algorithm: A :class:`ClusteringSpec`, a roster name, or a callable;
        the last two mean a spec with no fixed options.
    :param items: A ``SymmetricDistanceMatrix`` with complete storage, a
        prebuilt comparison, or an ``oefp.OEFPBatch``; passed to the
        algorithm and the scorer untouched.
    :param parameter: The keyword to sweep.
    :param values: The grid, a non-empty iterable (not a string).
    :param criterion: A validity index; None means ``"silhouette"``, or
        ``"isim_silhouette"`` for fingerprints.
    :param max_noise_fraction: Rows above this noise fraction are ineligible.
    :param min_clusters: Rows with fewer clusters are ineligible.
    :param max_clusters: Rows with more clusters are ineligible.
    :param report_options: Mapping forwarded to the scorer unchanged, with
        the opt-in stage flag the criterion needs added to a copy.
    :returns: A :class:`ParameterSelection`.
    :raises TypeError: For a non-string ``parameter``, a string or
        non-iterable ``values``, a non-mapping ``report_options``, a
        non-string ``criterion``, a stage flag in ``report_options`` that
        is not a bool, a bound of the wrong type, or an unsupported
        ``items`` kind.
    :raises ValueError: For an empty ``parameter`` or grid, a bound out of
        range, sparse storage, a criterion that is not a validity index or
        that the scorer does not produce, a pair-rank criterion with a
        prebuilt comparison, or a criterion whose stage flag the caller set
        to False.
    """
    spec = (algorithm if isinstance(algorithm, ClusteringSpec)
            else ClusteringSpec(algorithm))
    if not isinstance(parameter, str):
        raise TypeError(
            f"parameter must be a str naming a keyword of {spec.name}, not "
            f"{type(parameter).__name__}")
    if not parameter:
        raise ValueError("parameter must not be empty")
    if isinstance(values, str):
        raise TypeError("values must be an iterable of parameter values, not a str")
    try:
        grid = tuple(values)
    except TypeError as exc:
        raise TypeError("values must be an iterable of parameter values") from exc
    if not grid:
        raise ValueError("values must not be empty")
    if report_options is None:
        options = {}
    elif isinstance(report_options, Mapping):
        options = dict(report_options)
    else:
        raise TypeError(
            "report_options must be a mapping or None, not "
            f"{type(report_options).__name__}")
    max_noise_fraction, min_clusters, max_clusters = _validate_bounds(
        max_noise_fraction, min_clusters, max_clusters)
    scorer = _scorer_for(items)
    criterion = _resolve_criterion(criterion, scorer, options)
    direction = _CRITERIA[criterion]

    rows = []
    for value in grid:
        result = spec.run(items, **{parameter: value})
        report = scorer.function(result, items, **options)
        score = float(getattr(report, criterion))
        rejection = _first_violation(report, max_noise_fraction, min_clusters,
                                     max_clusters)
        rows.append(SweepRow(value, result, report, score, rejection is None,
                             rejection))
    return ParameterSelection(spec, parameter, criterion, rows,
                              _winner_index(rows, direction))
