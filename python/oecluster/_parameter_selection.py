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

import types
from collections.abc import Callable
from typing import Any

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

    Nothing is validated at construction. The algorithm validates its
    options when it runs, exactly as a direct call would.

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
                and dict(self._options) == dict(other._options))

    def __repr__(self):
        parts = [repr(self._name)]
        parts.extend(f"{key}={value!r}" for key, value in self._options.items())
        return f"ClusteringSpec({', '.join(parts)})"
