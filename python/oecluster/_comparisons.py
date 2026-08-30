"""
Comparison construction shared by :func:`oecluster.pdist` and
:func:`oecluster.cdist`.

The module owns two registries keyed on the lowercase comparison name: a
builder that turns items plus keyword options into a C++ comparison object,
and an optional normalizer that rewrites the item list before labels and set
sizes are fixed. Keeping both here lets a comparison add an item-rewriting
step (conformer expansion, complete-case filtering) without ``pdist`` and
``cdist`` growing a branch per comparison.
"""

from collections.abc import Callable
from typing import Any

from . import oecluster as _oecluster
from .oecluster import FingerprintComparison as _FingerprintComparison
from .oecluster import FingerprintOptions, ROCSOptions, SuperposeOptions
from .oecluster import ROCSComparison as _ROCSComparison
from .oecluster import SuperposeComparison as _SuperposeComparison

Builder = Callable[[list[Any], bool, dict[str, Any], bool], tuple[Any, str]]
Normalizer = Callable[[list[Any], dict[str, Any]],
                      tuple[list[Any], list[list[Any]]]]

_BUILDERS: dict[str, Builder] = {}
_NORMALIZERS: dict[str, Normalizer] = {}


def register_comparison(name, builder, normalizer=None):
    """
    Register a comparison builder, and optionally a normalizer, under a name.

    :param name: Lowercase comparison name accepted by :func:`oecluster.pdist`
                 and :func:`oecluster.cdist`.
    :param builder: Callable ``(items, similarity, kwargs, symmetric)``
                    returning ``(comparison_obj, comparison_name)``. It must
                    pop every keyword it consumes from ``kwargs`` and raise
                    :class:`TypeError` for anything left over.
    :param normalizer: Optional callable ``(items, kwargs)`` returning
                       ``(kept_items, excluded)``. It must read ``kwargs``
                       without popping, because the builder runs afterwards.
    """
    _BUILDERS[name] = builder
    if normalizer is not None:
        _NORMALIZERS[name] = normalizer


def supported_comparisons():
    """
    List the registered comparison names.

    :returns: Sorted list of names accepted by ``pdist`` and ``cdist``.
    """
    return sorted(_BUILDERS)


def _resolve_name(comparison):
    """
    Normalize and validate a comparison name.

    :param comparison: Comparison name in any case.
    :returns: The lowercase registered name.
    :raises ValueError: If no builder is registered under the name.
    """
    name = comparison.lower()
    if name not in _BUILDERS:
        valid = ", ".join(repr(known) for known in supported_comparisons())
        raise ValueError(
            f"Unknown comparison: {comparison!r}. Valid options: {valid}")
    return name


def extract_labels(items):
    """
    Extract labels from molecular items.

    :param items: List of molecules or design units.
    :returns: List of labels (molecule titles or indices).
    """
    labels = []

    try:
        # Try to import openeye to check types
        from openeye import oechem

        for idx, item in enumerate(items):
            if isinstance(item, oechem.OEMolBase):
                title = item.GetTitle()
                labels.append(title if title else f"mol_{idx}")
            else:
                labels.append(f"item_{idx}")
    except (ImportError, AttributeError):
        # If openeye not available or not a molecule, use indices
        labels = [f"item_{idx}" for idx in range(len(items))]

    return labels


def normalize_items(comparison, items, kwargs):
    """
    Apply a comparison's item-normalization step.

    Called by ``pdist`` and ``cdist`` before labels are extracted and before
    set sizes are fixed, so a comparison may both grow the item list
    (conformer expansion) and shrink it (complete-case filtering).

    :param comparison: Comparison name.
    :param items: Input items.
    :param kwargs: Comparison keyword options, read but never consumed.
    :returns: Tuple of ``(kept_items, excluded)``, where ``excluded`` is a
              list of ``[original_index, reason]`` pairs.
    :raises ValueError: If the comparison name is unknown.
    """
    name = _resolve_name(comparison)
    normalizer = _NORMALIZERS.get(name)
    if normalizer is None:
        return list(items), []
    return normalizer(list(items), kwargs)


def build_comparison(items, comparison, similarity, kwargs, *, symmetric):
    """
    Build a C++ comparison object from items and a comparison name.

    :param items: Normalized item list.
    :param comparison: Comparison name.
    :param similarity: Whether to compute similarities instead of distances.
    :param kwargs: Comparison-specific options, consumed by the builder.
    :param symmetric: True when the caller is ``pdist``. Builders use it to
                      reject options that are only valid for a rectangular
                      result.
    :returns: Tuple of ``(comparison_obj, comparison_name, params)``.
    :raises ValueError: If the comparison or an option value is unknown.
    :raises TypeError: If unknown kwargs remain after the builder runs.
    """
    name = _resolve_name(comparison)
    params = {'comparison_type': name, 'similarity': similarity}
    comparison_obj, comparison_name = _BUILDERS[name](
        items, similarity, kwargs, symmetric)
    # A correct builder pops all consumed options and raises TypeError for the
    # rest. If kwargs is non-empty here, the builder has a bug.
    if kwargs:
        raise RuntimeError(
            f"the {name!r} builder returned without consuming {list(kwargs)}. "
            f"A builder must pop every option it uses and raise TypeError for "
            f"the rest; reaching here is a bug in the builder, not bad input.")
    return comparison_obj, comparison_name, params


_FINGERPRINT_KEYS = (
    'fp_type',
    'storage',
    'numbits',
    'metric',
    'radius',
    'min_distance',
    'max_distance',
    'torsion_atom_count',
    'use_chirality',
    'p',
    'tversky_alpha',
    'tversky_beta',
)


_SPARSE_STORAGES = ('sparse', 'sparse_count')

# Which families read each per-family option, and which option replaces it for
# a family that does not. Keyed on every family so an unrecognized fp_type
# falls through to the C++ constructor's own error instead of being blamed here.
_FAMILY_ONLY_KEYS = {
    'radius': ('morgan',),
    'min_distance': ('atom_pair', 'topological_atom_pair'),
    'max_distance': ('atom_pair', 'topological_atom_pair'),
    'torsion_atom_count': ('topological_torsions',),
}

_FAMILY_REPLACEMENT = {
    'morgan': 'radius',
    'atom_pair': 'min_distance/max_distance',
    'topological_atom_pair': 'min_distance/max_distance',
    'topological_torsions': 'torsion_atom_count',
}

_METRIC_ONLY_KEYS = {
    'p': 'minkowski',
    'tversky_alpha': 'tversky',
    'tversky_beta': 'tversky',
}


def reject_inapplicable_fingerprint_kwargs(named, *, fp_type, storage, metric):
    """
    Reject an explicitly named option the rest of the configuration ignores.

    Only the Python surface can apply these rules. C++ receives a fully
    populated struct and cannot tell a deliberate value from a default, so it
    ignores unused fields instead -- a default-constructed ``FingerprintOptions``
    must stay a working Morgan configuration. Every message names both the
    option that would have been silently ignored and the one that replaces it.

    ``use_chirality`` applies to all four families and is never rejected.

    :param named: Option names the caller passed explicitly.
    :param fp_type: Selected fingerprint family.
    :param storage: Selected storage.
    :param metric: Selected metric name.
    :raises TypeError: If a named option does not apply to this configuration.
    """
    named = set(named)
    family = (fp_type or 'morgan').lower()
    store = (storage or 'binary').lower()
    metric_name = (metric or 'tanimoto').lower()

    if 'numbits' in named and store in _SPARSE_STORAGES:
        raise TypeError(
            f"numbits does not apply to storage={store!r}: a sparse "
            f"fingerprint keeps its family's own domain rather than folding to "
            f"a chosen width. Drop numbits, or use storage='binary' or "
            f"storage='count'.")

    if family in _FAMILY_REPLACEMENT:
        for key, families in _FAMILY_ONLY_KEYS.items():
            if key in named and family not in families:
                raise TypeError(
                    f"{key} does not apply to fp_type={family!r}; it belongs "
                    f"to {' and '.join(repr(f) for f in families)}. Use "
                    f"{_FAMILY_REPLACEMENT[family]} instead, or select one of "
                    f"those families.")

    for key, owner in _METRIC_ONLY_KEYS.items():
        if key in named and metric_name != owner:
            raise TypeError(
                f"{key} does not apply to metric={metric_name!r}; it belongs "
                f"to {owner!r}. Drop {key}, or select metric={owner!r}.")


def _build_fingerprint(items, similarity, kwargs, symmetric):
    """Build a :class:`FingerprintComparison` from keyword options."""
    # Captured before the pops: the explicitness rules turn on which options
    # the caller named, which is exactly what the pops destroy.
    named = set(kwargs)

    opts = FingerprintOptions()
    opts.similarity = similarity
    for key in _FINGERPRINT_KEYS:
        if key in kwargs:
            setattr(opts, key, kwargs.pop(key))
    if kwargs:
        raise TypeError(
            f"Unknown kwargs for fingerprint comparison: {list(kwargs)}")

    reject_inapplicable_fingerprint_kwargs(
        named, fp_type=opts.fp_type, storage=opts.storage, metric=opts.metric)

    # Tversky is symmetric only when alpha == beta, which the condensed pdist
    # form requires. Rejecting here names both parameters; the C++
    # ValidateForPDist message can only report the metric.
    if (symmetric and opts.metric.lower() == 'tversky'
            and opts.tversky_alpha != opts.tversky_beta):
        raise ValueError(
            "pdist requires a symmetric metric, but tversky with "
            f"tversky_alpha={opts.tversky_alpha} and "
            f"tversky_beta={opts.tversky_beta} is asymmetric. Use equal "
            "alpha and beta, or compute a rectangular result with cdist.")

    return _FingerprintComparison(items, opts), "fingerprint"


def _build_rocs(items, similarity, kwargs, symmetric):
    """Build a :class:`ROCSComparison` from keyword options."""
    opts = ROCSOptions()
    opts.similarity = similarity
    if 'score_type' in kwargs:
        st = kwargs.pop('score_type')
        score_map = {
            'combo_norm': _oecluster.ROCSScoreType_ComboNorm,
            'combo': _oecluster.ROCSScoreType_Combo,
            'shape': _oecluster.ROCSScoreType_Shape,
            'color': _oecluster.ROCSScoreType_Color,
        }
        if st not in score_map:
            raise ValueError(f"Unknown ROCS score type: {st}")
        opts.score_type = score_map[st]
    if 'color_ff_type' in kwargs:
        opts.color_ff_type = kwargs.pop('color_ff_type')
    if kwargs:
        raise TypeError(f"Unknown kwargs for rocs comparison: {list(kwargs)}")
    return _ROCSComparison(items, opts), "rocs"


def _build_superpose(items, similarity, kwargs, symmetric, *, default_method):
    """Build a :class:`SuperposeComparison` from keyword options."""
    opts = SuperposeOptions()
    opts.similarity = similarity
    method = kwargs.pop('method', None) or default_method
    if method is not None:
        method_map = {
            'global_carbon_alpha': _oecluster.SuperposeMethod_GlobalCarbonAlpha,
            'global': _oecluster.SuperposeMethod_Global,
            'ddm': _oecluster.SuperposeMethod_DDM,
            'weighted': _oecluster.SuperposeMethod_Weighted,
            'sse': _oecluster.SuperposeMethod_SSE,
            'sitehopper': _oecluster.SuperposeMethod_SiteHopper,
        }
        if method not in method_map:
            raise ValueError(f"Unknown superpose method: {method}")
        opts.method = method_map[method]
    if 'score_type' in kwargs:
        st = kwargs.pop('score_type')
        st_map = {
            'auto': _oecluster.SuperposeScoreType_Auto,
            'rmsd': _oecluster.SuperposeScoreType_RMSD,
            'tanimoto': _oecluster.SuperposeScoreType_Tanimoto,
            'patch_score': _oecluster.SuperposeScoreType_PatchScore,
        }
        if st not in st_map:
            raise ValueError(f"Unknown superpose score type: {st}")
        opts.score_type = st_map[st]
    if 'predicate' in kwargs:
        opts.predicate = kwargs.pop('predicate')
    if 'ref_predicate' in kwargs:
        opts.ref_predicate = kwargs.pop('ref_predicate')
    if 'fit_predicate' in kwargs:
        opts.fit_predicate = kwargs.pop('fit_predicate')
    if kwargs:
        raise TypeError(
            f"Unknown kwargs for superpose comparison: {list(kwargs)}")
    comparison_obj = _SuperposeComparison(items, opts)
    return comparison_obj, comparison_obj.ComparisonName()


register_comparison("fingerprint", _build_fingerprint)
register_comparison("rocs", _build_rocs)
register_comparison(
    "superpose",
    lambda items, similarity, kwargs, symmetric: _build_superpose(
        items, similarity, kwargs, symmetric, default_method=None))
register_comparison(
    "sitehopper",
    lambda items, similarity, kwargs, symmetric: _build_superpose(
        items, similarity, kwargs, symmetric, default_method="sitehopper"))
