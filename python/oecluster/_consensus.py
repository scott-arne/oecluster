"""
Consensus clustering over an ensemble of partitions.

:func:`consensus` accumulates how often each pair of items was placed in the
same cluster, divides by how often the pair was observed together, and reads
a partition off the resulting co-association distances. The ensemble is
whatever the earlier workflow slices produce: the rows of a
:func:`select_parameter` sweep, the resamples of a :func:`cluster_stability`
run, or a sequence of partitions built by hand.

The base class is imported eagerly because :class:`ConsensusResult` subclasses
it; the package imports this module at the end of its ``__init__`` for that
reason. Everything else from the package is resolved at call time, as
:mod:`._parameter_selection` does.
"""

import math
import numbers
import os
from typing import NamedTuple

import numpy as np

from . import ClusteringResult, _gate
from . import oecluster as _oecluster
from ._parameter_selection import ClusteringSpec, _package, _roster
from ._stability import _dense_codes, _positions, _thread_count

# Marks a member that covers every item. A private object rather than None so
# that a caller's own ``(None, labels)`` tuple reaches _positions and is
# refused instead of being credited with observing every item.
_FULL = object()

# Resolved only when no ``method`` extracts the partition; see _resolve_threshold.
_DEFAULT_THRESHOLD = 0.5

# Roster entries that cluster fingerprints or molecules and cannot read a
# distance matrix. Refusing them by name turns a TypeError from deep inside
# the entry point into a refusal that says what is wrong.
_NON_MATRIX_METHODS = frozenset({
    "bitbirch", "bitbirch_recluster", "bitbirch_refine", "murcko",
})


def _member_labels(labels, expected, what):
    """Validate one member's (or a method result's) labels.

    The one-dimensional check is not redundant with the length check:
    ``ClusteringResult`` builds its labels with ``numpy.asarray(list(labels),
    dtype=numpy.intp)``, so a list of one-element lists becomes an ``(n, 1)``
    array whose ``num_samples`` is still ``n``, and such an array would
    broadcast rather than index.

    :raises TypeError: For a non-integer label array.
    :raises ValueError: For a non-one-dimensional array or the wrong length.
    """
    values = np.asarray(labels)
    if values.dtype.kind not in "iu":
        raise TypeError(
            f"{what} labels must be integers, not {values.dtype}")
    if values.ndim != 1:
        raise ValueError(
            f"{what} labels must be one-dimensional, not shape {values.shape}")
    if values.size != expected:
        raise ValueError(
            f"{what} has {values.size} labels but {expected} were expected")
    return values


def _resolve_members(ensemble):
    """Turn an ensemble into ``(members, inferred)``.

    A member is ``(positions, labels)`` with ``positions`` _FULL for a full
    partition. ``inferred`` is the item count the ensemble implies, from the
    first full member or a :class:`ClusterStability`'s reference, or None
    when every member is partial. Reconciling it with the caller's
    ``num_items`` is the caller's job, so that the ensemble's own structure
    is settled first, as the validation order requires.

    :raises TypeError: For an ensemble kind this function does not accept.
    :raises ValueError: For an empty ensemble or members of differing size.
    """
    package = _package()

    if isinstance(ensemble, package.ClusterStability):
        if ensemble.indices is None:
            raise ValueError(
                "a ClusterStability contributes its resamples only when it "
                "kept them; rerun cluster_stability with keep_partitions=True")
        members = [(positions, labels)
                   for positions, labels in zip(ensemble.indices,
                                                ensemble.labels)]
        return members, ensemble.reference.num_samples

    if isinstance(ensemble, package.ParameterSelection):
        rows = ensemble.rows
        if not rows:
            raise ValueError("the ensemble is empty")
        return ([(_FULL, row.result.labels) for row in rows],
                rows[0].result.num_samples)

    if isinstance(ensemble, (str, bytes)) or not hasattr(ensemble, "__iter__"):
        raise TypeError(
            "consensus() expects a ClusterStability, a ParameterSelection, or "
            "a sequence of ClusteringResult and (positions, labels) members, "
            f"not {type(ensemble).__name__}")

    elements = list(ensemble)
    if not elements:
        raise ValueError("the ensemble is empty")

    members = []
    total = None
    for index, element in enumerate(elements):
        if isinstance(element, package.ClusteringResult):
            # Only the first full member sets the inferred count. Later ones
            # are checked against the resolved count in _validate_members,
            # which runs after num_items itself has been validated, so a
            # malformed num_items is reported before a member mismatch.
            if total is None:
                total = element.num_samples
            members.append((_FULL, element.labels))
        elif isinstance(element, tuple) and len(element) == 2:
            members.append(element)
        else:
            raise TypeError(
                f"member {index} must be a ClusteringResult or a "
                f"(positions, labels) tuple, not {type(element).__name__}")
    return members, total


def _validate_members(members, num_items):
    """Check every member and return ``(positions, dense codes)`` pairs."""
    checked = []
    for index, (positions, labels) in enumerate(members):
        what = f"member {index}"
        if positions is _FULL:
            values = _member_labels(labels, num_items, what)
            checked.append((_FULL, _dense_codes(values)))
            continue
        chosen = _positions(positions, num_items)
        values = _member_labels(labels, chosen.size, what)
        checked.append((chosen, _dense_codes(values)))
    return checked


def _item_total(num_items):
    if isinstance(num_items, bool) or not isinstance(num_items, numbers.Integral):
        raise TypeError("num_items must be an int or None")
    if num_items < 2:
        raise ValueError("num_items must be at least 2")
    return int(num_items)


def _resolve_threshold(threshold, method):
    """The merge threshold, or None when a spec extracts the partition.

    ``threshold`` defaults to None rather than 0.5 so that an omitted value
    stays distinguishable from an explicit one; otherwise every ``method=``
    call would look like a caller asking for both.
    """
    if method is not None:
        if threshold is not None:
            raise ValueError(
                "threshold and method are mutually exclusive: a threshold the "
                "spec cannot see would be silently ignored")
        return None
    if threshold is None:
        return _DEFAULT_THRESHOLD
    if isinstance(threshold, bool) or not isinstance(threshold, numbers.Real):
        raise TypeError("threshold must be a real number")
    value = float(threshold)
    if math.isnan(value) or not 0.0 <= value <= 1.0:
        raise ValueError("threshold must be between 0 and 1")
    return value


def _resolve_method(method):
    if method is None:
        return None
    spec = method if isinstance(method, ClusteringSpec) else ClusteringSpec(method)
    # Compared by identity rather than by name: a foreign callable takes its
    # own ``__name__``, and one that happens to be called ``murcko`` is still
    # a caller's function that may well read a matrix.
    roster = _roster()
    if any(spec.algorithm is roster[name] for name in _NON_MATRIX_METHODS):
        raise ValueError(
            f"method={spec.name!r} clusters fingerprints or molecules and "
            "cannot read the consensus distance matrix; choose an entry point "
            "that takes a SymmetricDistanceMatrix")
    return spec


def _native_options(num_threads):
    options = _oecluster.ConsensusOptions()
    options.num_threads = num_threads
    return options


def _build_matrix(members, num_items, options, output):
    """Accumulate the co-association matrix and wrap it."""
    package = _package()
    offsets = [0]
    chunks_positions = []
    chunks_labels = []
    for positions, codes in members:
        if positions is _FULL:
            positions = np.arange(num_items, dtype=np.intp)
        chunks_positions.append(positions)
        chunks_labels.append(codes)
        offsets.append(offsets[-1] + int(positions.size))
    flat_positions = np.concatenate(chunks_positions)
    flat_labels = np.concatenate(chunks_labels)

    storage = (_oecluster.MMapStorage(os.fspath(output), num_items)
               if output is not None else _oecluster.DenseStorage(num_items))
    summary = _oecluster.coassociation_distances(
        num_items,
        _oecluster.SizeTVector(offsets),
        _oecluster.SizeTVector(flat_positions.tolist()),
        _oecluster.IntVector(flat_labels.tolist()),
        storage,
        options)

    matrix = package.SymmetricDistanceMatrix(
        storage, "consensus", None,
        {"num_partitions": int(summary.num_partitions)},
        {"is_distance": True, "zero_self": True, "data_integrity": "complete"})
    # A co-association distance obeys the triangle inequality for some
    # ensembles and not others, and a method= spec runs through the ordinary
    # clustering gate, so the probe records what is true of this one.
    matrix._facts.update(_gate.probe_triangle(matrix.condensed, num_items))
    return matrix, storage, int(summary.unobserved_pairs), int(summary.num_partitions)


def _clusters_from_labels(labels, num_clusters):
    return tuple(tuple(int(item) for item in np.flatnonzero(labels == label))
                 for label in range(num_clusters))


class ConsensusRecord(NamedTuple):
    """One consensus cluster and the evidence behind it."""

    label: int
    size: int
    cluster_consensus: float


class ConsensusResult(ClusteringResult):
    """Outcome of :func:`consensus`: a partition plus its evidence.

    The partition is canonical, `0..K-1` with `-1` for noise and
    ``clusters[k]`` holding exactly the items labelled `k`, because the
    package's own consumers require it: ``cluster_report`` refuses a label
    that is not its cluster's ordinal, and both it and
    ``partition_agreement`` refuse labels beyond the signed 32-bit range.

    Attributes other than the inherited ``labels`` and ``clusters`` are set
    here and read through properties, the convention every other result
    subclass follows.
    """

    def __init__(self, labels, clusters, *, matrix, threshold, spec,
                 num_partitions, unobserved_pairs, records, item_consensus,
                 agreement):
        super().__init__(labels, clusters)
        item_consensus.setflags(write=False)
        self._matrix = matrix
        self._threshold = threshold
        self._spec = spec
        self._num_partitions = num_partitions
        self._unobserved_pairs = unobserved_pairs
        self._records = tuple(records)
        self._item_consensus = item_consensus
        self._agreement = tuple(float(value) for value in agreement)

    @property
    def method(self):
        return "consensus"

    @property
    def matrix(self):
        """The consensus :class:`SymmetricDistanceMatrix`."""
        return self._matrix

    @property
    def threshold(self):
        """The merge threshold, or None when a spec extracted the partition."""
        return self._threshold

    @property
    def spec(self):
        """The :class:`ClusteringSpec` that extracted the partition, or None."""
        return self._spec

    @property
    def num_partitions(self):
        """Number of ensemble members."""
        return self._num_partitions

    @property
    def unobserved_pairs(self):
        """Pairs no member observed together."""
        return self._unobserved_pairs

    @property
    def records(self):
        """Tuple of :class:`ConsensusRecord` in ascending label order."""
        return self._records

    @property
    def item_consensus(self):
        """Read-only mean co-association of each item with its own cluster."""
        return self._item_consensus

    @property
    def agreement(self):
        """Tuple of one adjusted Rand index per member, NaN where undefined."""
        return self._agreement

    @property
    def mean_agreement(self):
        """Mean of the defined agreement entries, NaN when none is defined."""
        defined = [value for value in self._agreement if not math.isnan(value)]
        return float(np.mean(defined)) if defined else math.nan

    @property
    def columns(self):
        """Column names of :meth:`to_table`."""
        return ("label", "size", "cluster_consensus")

    def to_table(self):
        """Return a fresh list of tuples, one per record, matching
        :attr:`columns`."""
        return [tuple(record) for record in self._records]

    def __repr__(self):
        extraction = (f"threshold={self._threshold}" if self._spec is None
                      else f"spec={self._spec!r}")
        header = (f"ConsensusResult(num_partitions={self._num_partitions}, "
                  f"{extraction}, num_clusters={self.num_clusters}, "
                  f"mean_agreement={self.mean_agreement:.4f})")
        names = self.columns
        cells = [(str(record.label), str(record.size),
                  f"{record.cluster_consensus:.4f}")
                 for record in self._records]
        widths = [max([len(name)] + [len(cell[i]) for cell in cells])
                  for i, name in enumerate(names)]
        lines = [header,
                 "   " + "  ".join(f"{name:<{width}}"
                                   for name, width in zip(names, widths))]
        for cell in cells:
            lines.append(("   " + "  ".join(f"{text:<{width}}"
                                            for text, width in zip(cell, widths))).rstrip())
        return "\n".join(lines)


def consensus(ensemble, *, num_items=None, threshold=None, method=None,
              noise="singletons", num_threads=0, output=None):
    """Combine an ensemble of partitions into one, with its evidence.

    Every pair of items gets a co-association distance, ``1 - co / obs``,
    where ``co`` counts the members that placed both items in one cluster and
    ``obs`` the members that observed both. The default extraction unions
    every pair whose co-association is at least ``threshold``; passing
    ``method`` instead runs any matrix-consuming clustering spec on the
    co-association matrix.

    :param ensemble: A ``ClusterStability`` (its resamples), a
        ``ParameterSelection`` (its rows), or a sequence of
        ``ClusteringResult`` and ``(positions, labels)`` members.
    :param num_items: The item count; inferred from a full member or a
        stability reference, required when every member is partial, and
        reconciled with the inferred value when both exist.
    :param threshold: Co-association fraction in ``[0, 1]`` for the default
        extraction; defaults to 0.5. Mutually exclusive with ``method``.
    :param method: ``None``, or a ``ClusteringSpec``, roster name or callable
        that accepts a ``SymmetricDistanceMatrix``. The matrix is passed by
        reference, not copied, so the method must not write to it: a mutation
        leaves the retained ``matrix`` disagreeing with ``num_partitions``,
        ``unobserved_pairs`` and the consensus statistics, all of which are
        computed from it. The same holds for writing through
        ``result.matrix.condensed`` after the call.
    :param noise: ``"singletons"``, ``"grouped"`` or ``"excluded"``, forwarded
        to :func:`partition_agreement`; it does not affect the matrix.
    :param num_threads: Worker threads for the native passes; 0 selects the
        hardware concurrency.
    :param output: Optional path for a memory-mapped matrix.
    :returns: A :class:`ConsensusResult`.
    :raises TypeError: For an argument of the wrong type, including ``bool``
        where an int or float is expected.
    :raises ValueError: For an argument out of range, an empty or
        unreconcilable ensemble, both ``threshold`` and ``method``, a
        non-matrix roster name, or a spec result of the wrong size or shape.
    """
    package = _package()
    members, inferred = _resolve_members(ensemble)
    if num_items is None:
        if inferred is None:
            raise ValueError(
                "num_items is required when the ensemble holds only partial "
                "members: their positions alone do not bound the item count")
        num_items = _item_total(inferred)
    else:
        num_items = _item_total(num_items)
        if inferred is not None and inferred != num_items:
            raise ValueError(
                f"num_items={num_items} disagrees with the ensemble, which "
                f"covers {inferred} items")
    members = _validate_members(members, num_items)
    threshold = _resolve_threshold(threshold, method)
    spec = _resolve_method(method)
    num_threads = _thread_count(num_threads)
    if output is not None and not isinstance(output, (str, os.PathLike)):
        raise TypeError("output must be a path or None")
    # Resolved before the matrix is built: by the time the first agreement
    # call would catch a misspelled mode, an O(N^2) matrix exists and
    # output= has left a file behind.
    package._noise_handling(noise)

    options = _native_options(num_threads)
    matrix, storage, unobserved_pairs, num_partitions = _build_matrix(
        members, num_items, options, output)

    if spec is None:
        # list() first: the kernel returns a SWIG IntVector, which numpy
        # would otherwise wrap as an object array, exactly as
        # ClusteringResult converts a native Labels() vector.
        labels = np.asarray(
            list(_oecluster.consensus_components(storage, threshold, options)),
            dtype=np.intp)
    else:
        result = spec.run(matrix)
        if not isinstance(result, package.ClusteringResult):
            raise TypeError(
                f"{spec!r} returned {type(result).__name__}, not a "
                "ClusteringResult")
        if result.num_samples != num_items:
            raise ValueError(
                f"{spec!r} returned {result.num_samples} labels for the "
                f"{num_items} consensus items")
        labels = _dense_codes(
            _member_labels(result.labels, num_items, f"{spec!r}"))

    num_clusters = int(labels.max()) + 1 if labels.size and labels.max() >= 0 else 0
    clusters = _clusters_from_labels(labels, num_clusters)

    strength = _oecluster.consensus_strength(
        storage, _oecluster.IntVector(labels.tolist()), options)
    item_consensus = np.asarray(list(strength.item_consensus),
                                dtype=np.float64)
    cluster_consensus = list(strength.cluster_consensus)
    records = tuple(
        ConsensusRecord(label, len(clusters[label]),
                        float(cluster_consensus[label]))
        for label in range(num_clusters))

    agreement = []
    for positions, codes in members:
        observed = labels if positions is _FULL else labels[positions]
        agreement.append(package.partition_agreement(
            observed, codes, noise=noise).adjusted_rand_index)

    return ConsensusResult(
        labels, clusters, matrix=matrix, threshold=threshold, spec=spec,
        num_partitions=num_partitions, unobserved_pairs=unobserved_pairs,
        records=records, item_consensus=item_consensus, agreement=agreement)
