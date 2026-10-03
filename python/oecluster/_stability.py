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
from typing import NamedTuple

import numpy as np

from . import _gate
from . import oecluster as _oecluster
from ._parameter_selection import ClusteringSpec, _is_fingerprint_batch, _package

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


def _integrity_after_nan_exclusion(subset):
    """What a NaN-free proper subset of a NaN-present source may claim.

    A descriptor comparison under ``missing="ignore"`` scores each pair on
    the dimensions both items have, which is subset-scored; an observed
    NaN then escalates the stamp to NaN-present and hides it (the
    comparison ranks an observed NaN above the policy). Excluding the NaN
    items uncovers the subset-scored matrix underneath, not a complete
    one. Under ``propagate`` and ``complete_case`` every remaining pair
    was scored on the full data, so the subset is complete. A descriptor
    matrix that records no policy (one saved before 5.16.0) is treated as
    subset-scored, which never claims more than the data supports.
    """
    if subset.comparison_name != "descriptor":
        return "complete"
    policy = subset.params.get("missing")
    if policy in ("propagate", "complete_case"):
        return "complete"
    return "subset_scored"


def _data_facts(subset, source_facts):
    """The facts a proper subset re-measures: the probe and NaN presence.

    The comparison facts (``is_distance``, ``zero_self``, ``triangle``) and
    a complete or subset-scored ``data_integrity`` describe the comparison
    and are inherited by the caller; the probe and the NaN presence describe
    the data. The NaN re-measure defers to ``_integrity_after_nan_exclusion``
    for what a NaN-free subset may claim.
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
        facts["data_integrity"] = (_integrity_after_nan_exclusion(subset)
                                   if finite else "nan_present")
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
    the subset (a NaN-free subset of an ignore-scored descriptor matrix is
    subset-scored, not complete). An ``oefp.OEFPBatch`` yields a new batch
    with the same fingerprint spec and the selected rows.

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


# --- cluster_stability -------------------------------------------------------

class ClusterStabilityRecord(NamedTuple):
    """One reference cluster's stability over the resamples it appeared in."""

    label: int
    size: int
    mean_jaccard: float
    dissolved: float
    recovered: float
    evaluated: int


class ClusterStability:
    """Read-only outcome of :func:`cluster_stability`.

    ``records`` holds one :class:`ClusterStabilityRecord` per non-noise
    reference cluster in ascending label order; ``jaccard`` is the ``(K, R)``
    matrix behind them, NaN where a cluster had no member in a resample;
    ``agreement`` is one adjusted Rand index per resample, NaN where
    :func:`partition_agreement` leaves it undefined.

    Assignment to any attribute, public or private, raises
    ``AttributeError``: the fields are set through ``object.__setattr__`` in
    the constructor and ``__setattr__`` rejects everything after, the same
    arrangement :class:`ParameterSelection` uses.
    """

    def __init__(self, spec, reference, records, jaccard, agreement,
                 resamples, fraction, seed, indices, labels):
        jaccard.setflags(write=False)
        object.__setattr__(self, "_spec", spec)
        object.__setattr__(self, "_reference", reference)
        object.__setattr__(self, "_records", tuple(records))
        object.__setattr__(self, "_jaccard", jaccard)
        object.__setattr__(self, "_agreement",
                           tuple(float(value) for value in agreement))
        object.__setattr__(self, "_resamples", resamples)
        object.__setattr__(self, "_fraction", fraction)
        object.__setattr__(self, "_seed", seed)
        object.__setattr__(self, "_indices", indices)
        object.__setattr__(self, "_labels", labels)

    def __setattr__(self, name, value):
        """ClusterStability is read-only; reject attribute assignment."""
        raise AttributeError(
            f"ClusterStability is read-only; cannot set {name!r}")

    @property
    def spec(self):
        """The :class:`ClusteringSpec` that was resampled."""
        return self._spec

    @property
    def reference(self):
        """The :class:`ClusteringResult` that was scored."""
        return self._reference

    @property
    def records(self):
        """Tuple of :class:`ClusterStabilityRecord` in ascending label order."""
        return self._records

    @property
    def jaccard(self):
        """Read-only ``(K, R)`` float64 array of best Jaccard values."""
        return self._jaccard

    @property
    def mean_jaccard(self):
        """Unweighted mean of the evaluated records' means, NaN if none.

        Unweighted on purpose: a size-weighted figure is dominated by the
        largest cluster and hides the small clusters resampling exists to
        expose; the records carry the sizes for a weighted figure.
        """
        means = [record.mean_jaccard for record in self._records
                 if record.evaluated > 0]
        return float(np.mean(means)) if means else math.nan

    @property
    def agreement(self):
        """Tuple of one adjusted Rand index per resample, NaN when undefined."""
        return self._agreement

    @property
    def mean_agreement(self):
        """Mean of the defined agreement entries, NaN when none is defined."""
        defined = [value for value in self._agreement if not math.isnan(value)]
        return float(np.mean(defined)) if defined else math.nan

    @property
    def resamples(self):
        """Number of resamples drawn."""
        return self._resamples

    @property
    def fraction(self):
        """Fraction of items kept per resample, as a float."""
        return self._fraction

    @property
    def seed(self):
        """The seed as resolved: an int, or None for fresh entropy."""
        return self._seed

    @property
    def indices(self):
        """Tuple of sorted ``intp`` position arrays, or None if not kept."""
        return self._indices

    @property
    def labels(self):
        """Tuple of ``intp`` label arrays matching ``indices``, or None."""
        return self._labels

    @property
    def columns(self):
        """Column names of :meth:`to_table`."""
        return ("label", "size", "mean_jaccard", "dissolved", "recovered",
                "evaluated")

    def to_table(self):
        """Return a fresh list of tuples, one per record, matching
        :attr:`columns`."""
        return [tuple(record) for record in self._records]

    def __repr__(self):
        header = (f"ClusterStability(spec={self._spec!r}, "
                  f"resamples={self._resamples}, fraction={self._fraction}, "
                  f"mean_jaccard={self.mean_jaccard:.4f}, "
                  f"mean_agreement={self.mean_agreement:.4f})")
        names = self.columns
        cells = [(str(record.label), str(record.size),
                  f"{record.mean_jaccard:.4f}", f"{record.dissolved:.4f}",
                  f"{record.recovered:.4f}", str(record.evaluated))
                 for record in self._records]
        widths = [max([len(name)] + [len(cell[i]) for cell in cells])
                  for i, name in enumerate(names)]
        lines = [header,
                 "   " + "  ".join(f"{name:<{width}}"
                                   for name, width in zip(names, widths))]
        for cell in cells:
            line = "   " + "  ".join(f"{text:<{width}}"
                                     for text, width in zip(cell, widths))
            lines.append(line.rstrip())
        return "\n".join(lines)


def _item_count(items):
    package = _package()
    if isinstance(items, package.SymmetricDistanceMatrix):
        count = items.num_samples
    elif _is_fingerprint_batch(items):
        count = items.size
    else:
        raise TypeError(
            "cluster_stability() expects a SymmetricDistanceMatrix or an "
            f"oefp.OEFPBatch, not {type(items).__name__}; compute a matrix "
            "from a raw item sequence with pdist first")
    if count < 2:
        raise ValueError(
            f"cluster_stability() requires at least 2 items, got {count}")
    return count


def _resample_count(resamples):
    if isinstance(resamples, bool) or not isinstance(resamples, numbers.Integral):
        raise TypeError("resamples must be an int")
    if resamples < 1:
        raise ValueError("resamples must be at least 1")
    return int(resamples)


def _fraction_value(fraction):
    if isinstance(fraction, bool) or not isinstance(fraction, numbers.Real):
        raise TypeError("fraction must be a real number")
    fraction = float(fraction)
    if math.isnan(fraction) or not 0.0 < fraction <= 1.0:
        raise ValueError("fraction must be in (0, 1]")
    return fraction


def _seed_value(seed):
    if seed is None:
        return None
    if isinstance(seed, bool) or not isinstance(seed, numbers.Integral):
        raise TypeError("seed must be a non-negative int or None")
    if seed < 0:
        raise ValueError(
            "seed must be non-negative, as numpy.random.default_rng requires")
    return int(seed)


def _reference_rows(reference_labels):
    """Record labels in ascending order and every item's global row.

    Built once from the full reference, with noise mapped to -1: a resample
    that lacks a lower or a middle cluster must not shift the later clusters
    into other rows, which a per-resample encoding would do.
    """
    labels = np.asarray(reference_labels)
    unique, inverse = np.unique(labels, return_inverse=True)
    noise_count = int(np.count_nonzero(unique < 0))  # negatives sort first
    rows = np.asarray(inverse, dtype=np.intp) - noise_count
    rows[rows < 0] = -1
    return unique[noise_count:], rows


def _dense_codes(labels):
    """Codes ``0..K'-1`` for one resample's non-negative labels, -1 kept.

    Labels are arbitrary non-negative integers, so ``bincount`` on the raw
    values would allocate by the largest label rather than by the sample
    count, and :func:`partition_agreement` accepts only 32-bit labels.
    """
    labels = np.asarray(labels)
    codes = np.full(labels.shape, -1, dtype=np.intp)
    clustered = labels >= 0
    if clustered.any():
        codes[clustered] = np.unique(labels[clustered], return_inverse=True)[1]
    return codes


def _best_jaccard(ref_rows, lab_codes, num_rows):
    """Best Jaccard of every reference row against one resample's clusters.

    NaN for a row with no member in the resample; 0.0 for one whose members
    all became noise, or when the resample has no cluster at all. Only the
    observed contingency cells are formed (at most one per item), so the
    workspace is O(m + K + K') and no K x K' table ever exists.
    """
    column = np.full(num_rows, np.nan, dtype=np.float64)
    in_reference = ref_rows >= 0
    if not in_reference.any():
        return column
    reference_sizes = np.bincount(ref_rows[in_reference], minlength=num_rows)
    column[reference_sizes > 0] = 0.0
    in_both = in_reference & (lab_codes >= 0)
    if not in_both.any():
        return column
    cluster_sizes = np.bincount(lab_codes[lab_codes >= 0])
    num_codes = cluster_sizes.size
    cells, counts = np.unique(
        ref_rows[in_both] * num_codes + lab_codes[in_both], return_counts=True)
    rows = cells // num_codes
    codes = cells % num_codes
    union = reference_sizes[rows] + cluster_sizes[codes] - counts
    np.maximum.at(column, rows, counts / union)
    return column


def _records(record_labels, row_of_label, jaccard):
    sizes = np.bincount(row_of_label[row_of_label >= 0],
                        minlength=record_labels.size)
    records = []
    for index, label in enumerate(record_labels):
        scores = jaccard[index]
        scores = scores[~np.isnan(scores)]
        evaluated = int(scores.size)
        if evaluated:
            mean = float(scores.mean())
            dissolved = float(
                np.count_nonzero(scores < _DISSOLVED_BELOW) / evaluated)
            recovered = float(
                np.count_nonzero(scores > _RECOVERED_ABOVE) / evaluated)
        else:
            mean = dissolved = recovered = math.nan
        records.append(ClusterStabilityRecord(
            int(label), int(sizes[index]), mean, dissolved, recovered,
            evaluated))
    return tuple(records)


def cluster_stability(algorithm, items, *, resamples=100, fraction=0.5, seed=0,
                      reference=None, noise="singletons", keep_partitions=True,
                      num_threads=0):
    """Score how well each cluster survives resampling of the items.

    Draws ``resamples`` subsets of ``round(fraction * N)`` items (at least
    two) without replacement, reruns the spec on each, and matches every
    reference cluster to its best counterpart by Jaccard overlap. Per
    cluster: the mean best Jaccard over the resamples it appeared in, and
    the fractions of those in which it dissolved (best Jaccard below 0.5) or
    was recovered (above 0.75), Hennig's clusterboot statistics. Per
    resample: the adjusted Rand index between the reference restricted to
    the subset and the resample's partition.

    Noise enters in two fixed ways. The Jaccard matching never treats noise
    as a cluster: a reference cluster's members that a resample labels noise
    count against its overlap through the union. ``noise`` governs only the
    agreement call, exactly as it does in :func:`partition_agreement`.

    :param algorithm: A :class:`ClusteringSpec`, a roster name, or a callable;
        a name or callable is wrapped as ``ClusteringSpec(algorithm)``.
    :param items: A ``SymmetricDistanceMatrix`` or an ``oefp.OEFPBatch`` of
        at least two items.
    :param resamples: Number of resamples, at least 1.
    :param fraction: Fraction of items kept per resample, in ``(0, 1]``;
        ``1`` selects every item every time.
    :param seed: A non-negative int for a reproducible draw, or None for
        fresh entropy.
    :param reference: The :class:`ClusteringResult` to score, over all
        ``items``; None runs the spec once on the full input.
    :param noise: ``"singletons"``, ``"grouped"`` or ``"excluded"``,
        forwarded to :func:`partition_agreement`; checked there on the first
        resample.
    :param keep_partitions: Whether the result retains each resample's
        positions and labels (two ``intp`` arrays per resample).
    :param num_threads: Worker threads for the dense and memory-mapped
        matrix gather; the spec's algorithm threads through its own options.
    :returns: A :class:`ClusterStability`.
    :raises TypeError: For an argument of the wrong type, including ``bool``
        where an int or float is expected.
    :raises ValueError: For an argument out of range, a reference of another
        size, a callable whose result does not match the item count, or one
        that returns the reference object again or relabels it during
        resampling, or labels that are not one-dimensional.
    """
    spec = (algorithm if isinstance(algorithm, ClusteringSpec)
            else ClusteringSpec(algorithm))
    package = _package()
    num_items = _item_count(items)
    resamples = _resample_count(resamples)
    fraction = _fraction_value(fraction)
    seed = _seed_value(seed)
    if not isinstance(keep_partitions, (bool, np.bool_)):
        raise TypeError(
            "keep_partitions must be True or False, not "
            f"{type(keep_partitions).__name__}")
    keep_partitions = bool(keep_partitions)
    num_threads = _thread_count(num_threads)
    if reference is not None:
        if not isinstance(reference, package.ClusteringResult):
            raise TypeError(
                "reference must be a ClusteringResult or None, not "
                f"{type(reference).__name__}")
        if reference.num_samples != num_items:
            raise ValueError(
                f"reference has {reference.num_samples} labels but items has "
                f"{num_items}")
    else:
        reference = spec.run(items)
        if reference.num_samples != num_items:
            raise ValueError(
                f"{spec!r} returned {reference.num_samples} labels for the "
                f"{num_items} reference items")

    subset_size = max(2, round(fraction * num_items))
    # A copy: the object behind ``reference`` may be mutated by a callable.
    reference_labels = np.array(reference.labels, dtype=np.intp)
    if reference_labels.ndim != 1:
        raise ValueError(
            f"reference labels must be one-dimensional, not shape "
            f"{reference_labels.shape}")
    record_labels, row_of_label = _reference_rows(reference_labels)
    jaccard = np.full((record_labels.size, resamples), np.nan, dtype=np.float64)
    agreement = []
    kept_indices = []
    kept_labels = []
    rng = np.random.default_rng(seed)
    for index in range(resamples):
        chosen = np.sort(rng.choice(num_items, subset_size, replace=False))
        chosen = chosen.astype(np.intp, copy=False)
        result = spec.run(take(items, chosen, num_threads=num_threads))
        if result is reference:
            raise ValueError(
                f"resample {index}: {spec!r} returned the reference result "
                "object itself; a clustering callable must return a new "
                "ClusteringResult on every call")
        if result.num_samples != subset_size:
            raise ValueError(
                f"resample {index}: {spec!r} returned {result.num_samples} "
                f"labels for {subset_size} items")
        ref_rows = row_of_label[chosen]
        # A column-shaped label array passes the sample-count check yet
        # broadcasts the scoring masks to (m, m).
        labels = np.asarray(result.labels)
        if labels.ndim != 1:
            raise ValueError(
                f"resample {index}: {spec!r} returned labels of shape "
                f"{labels.shape}; labels must be one-dimensional")
        codes = _dense_codes(labels)
        jaccard[:, index] = _best_jaccard(ref_rows, codes, record_labels.size)
        agreement.append(package.partition_agreement(
            ref_rows, codes, noise=noise).adjusted_rand_index)
        if keep_partitions:
            kept_indices.append(chosen)
            # A copy, not a view: a callable that reuses one result object
            # and relabels it in place must not rewrite earlier partitions.
            kept_labels.append(np.array(labels, dtype=np.intp))
    # Two checks because they catch different aliasing: the identity check
    # names the offending resample as soon as a callable hands back the
    # object it already returned, and this comparison catches a callable that
    # relabels the reference through a retained alias without returning it.
    # Either would leave ``reference`` describing a partition the statistics
    # never scored.
    if not np.array_equal(reference.labels, reference_labels):
        raise ValueError(
            f"{spec!r} changed the reference's labels while resampling; the "
            "scored partition no longer matches stability.reference")
    return ClusterStability(
        spec, reference, _records(record_labels, row_of_label, jaccard),
        jaccard, agreement, resamples, fraction, seed,
        tuple(kept_indices) if keep_partitions else None,
        tuple(kept_labels) if keep_partitions else None)
