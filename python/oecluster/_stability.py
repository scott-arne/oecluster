"""
Cluster stability by resampling, and the row-subset primitive beneath it.

:func:`take` selects items of a ``SymmetricDistanceMatrix`` or an
``oefp.OEFPBatch`` by position through the native gathers ``take_pairs`` and
``take_fingerprints``. :func:`cluster_stability` reruns a
:class:`ClusteringSpec` on repeated subsamples drawn without replacement,
matches every reference cluster to its best counterpart in each resample by
Jaccard overlap (Hennig's clusterboot statistics), and scores each resample's
partition against the reference with the adjusted Rand index.

Everything from the package is resolved at call time rather than imported at
the top, as :mod:`._parameter_selection` does: the package ``__init__``
imports this module before the entry points it needs are defined there.
"""

import math
import numbers

import numpy as np

from . import _gate
from . import oecluster as _oecluster
from ._parameter_selection import _is_fingerprint_batch, _package

# Hennig's clusterboot thresholds: a best Jaccard below the first means the
# cluster dissolved in that resample, above the second that it was recovered.
_DISSOLVED_BELOW = 0.5
_RECOVERED_ABOVE = 0.75


# --- take --------------------------------------------------------------------

def _positions(indices, num_items):
    """Validate ``indices`` for ``num_items`` items; returns an ``intp`` array.

    ``numpy.asarray`` first, then the dtype is checked rather than coerced:
    ``[0.9, 1.1]`` must not silently select rows 0 and 1, which is the
    ``operator.index`` rule the package applies to labels.

    :raises TypeError: For a non-1-D or non-integer array.
    :raises ValueError: For an empty, out-of-range or repeated position.
    """
    positions = np.asarray(indices)
    if positions.ndim != 1:
        raise TypeError(
            "indices must be a one-dimensional sequence of positions, not "
            f"{positions.ndim}-D")
    if positions.size == 0:
        raise ValueError("indices must select at least one item")
    if positions.dtype.kind not in "iu":
        raise TypeError(
            f"indices must be integers, not {positions.dtype}; float, boolean "
            "and string positions are refused rather than coerced")
    if bool(np.any(positions < 0)) or bool(np.any(positions >= num_items)):
        raise ValueError(
            f"indices must lie in [0, {num_items}) for {num_items} items")
    if np.unique(positions).size != positions.size:
        raise ValueError(
            "indices must be distinct; a repeated position would select an "
            "item twice")
    return positions.astype(np.intp)


def _thread_count(num_threads):
    if (isinstance(num_threads, bool)
            or not isinstance(num_threads, numbers.Integral)):
        raise TypeError("num_threads must be an int")
    if num_threads < 0:
        raise ValueError("num_threads must be non-negative")
    return int(num_threads)


def _data_facts(subset, source_facts):
    """The facts a proper subset re-measures: the probe and NaN presence.

    The comparison facts (``is_distance``, ``zero_self``, ``triangle``) and
    a complete or subset-scored ``data_integrity`` describe the comparison
    and are inherited by the caller; these two describe the data.
    """
    sparse = isinstance(subset.storage, _oecluster.SparseStorage)
    facts = {}
    if not sparse and source_facts["metric_probe"] != "not_run":
        # The from_condensed defaults, so a subset is probed as a fresh
        # precomputed matrix would be. A sparse subset is left alone: the
        # probe needs the condensed form, which a sparse matrix would have to
        # densify, and inheriting never claims less than the source did.
        facts.update(_gate.probe_triangle(subset.condensed, subset.num_samples))
    if source_facts["data_integrity"] == "nan_present":
        # The sparse scan reads the stored entries, as the metric gate does;
        # ``condensed`` on a sparse matrix is a densified copy.
        if sparse:
            finite = all(math.isfinite(value)
                         for _, _, value in subset.storage._entries())
        else:
            finite = bool(np.isfinite(subset.condensed).all())
        facts["data_integrity"] = "complete" if finite else "nan_present"
    return facts


def _take_matrix(matrix, positions, num_threads):
    package = _package()
    source = matrix.storage
    count = int(positions.size)
    sparse = isinstance(source, _oecluster.SparseStorage)
    destination = (_oecluster.SparseStorage(count, source.Cutoff()) if sparse
                   else _oecluster.DenseStorage(count))
    _oecluster.take_pairs(source, _oecluster.SizeTVector(positions.tolist()),
                          destination, num_threads)
    labels = matrix.labels
    sliced = [labels[int(p)] for p in positions] if len(labels) else None
    subset = package.SymmetricDistanceMatrix(
        destination, matrix.comparison_name, sliced, dict(matrix.params),
        matrix.facts)
    if count < matrix.num_samples:
        # A proper subset re-measures the data-dependent facts; a permutation
        # of every item excluded nothing and keeps them all. ``_facts`` is the
        # new matrix's own dict (``facts`` hands out copies), written once
        # here before anyone else can see the subset.
        subset._facts.update(_data_facts(subset, matrix.facts))
    return subset


def take(items, indices, *, num_threads=0):
    """Select items by position, keeping the input's kind.

    A ``SymmetricDistanceMatrix`` yields a new matrix over the selected items
    in the given order (dense and memory-mapped sources become
    ``DenseStorage``, a sparse source stays sparse with the same cutoff), with
    the comparison name, parameters, sliced labels and the comparison facts
    carried over. A permutation of every item inherits every fact; a proper
    subset of a dense or memory-mapped source re-runs the triangle probe when
    the source's had run, and a source stamped NaN-present is re-measured on
    the subset. An ``oefp.OEFPBatch`` yields a new batch with the same
    fingerprint spec and the selected rows.

    :param items: A ``SymmetricDistanceMatrix`` or an ``oefp.OEFPBatch``.
    :param indices: Integer positions, distinct and in range, in any order.
    :param num_threads: Worker threads for the dense and memory-mapped gather
        (0 selects the hardware concurrency); ignored by the sparse and
        fingerprint gathers, which are single passes.
    :returns: A matrix or batch of the same kind as ``items``.
    :raises TypeError: For another ``items`` kind, a non-1-D or non-integer
        ``indices`` array, or a non-integer ``num_threads``.
    :raises ValueError: For empty, out-of-range or repeated positions, or a
        negative ``num_threads``.
    """
    package = _package()
    if isinstance(items, package.SymmetricDistanceMatrix):
        positions = _positions(indices, items.num_samples)
        return _take_matrix(items, positions, _thread_count(num_threads))
    if _is_fingerprint_batch(items):
        import oefp
        positions = _positions(indices, items.size)
        _thread_count(num_threads)
        native = _oecluster.take_fingerprints(
            items, _oecluster.SizeTVector(positions.tolist()))
        return oefp.OEFPBatch._from_native(native)
    raise TypeError(
        "take() expects a SymmetricDistanceMatrix or an oefp.OEFPBatch, not "
        f"{type(items).__name__}")
