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


def _set_if_given(opts, key, value):
    """Assign ``value`` to ``opts.key`` unless it is ``None``.

    ``None`` means "not specified". That is already the convention on the
    public comparison constructors -- ``FingerprintComparison.__new__`` and its
    siblings in ``__init__.py`` declare every option as ``None`` and skip the
    assignment rather than pushing ``None`` at a SWIG setter. The registry
    builders have to agree with them, or a config-driven caller who passes
    ``fp_type=None`` to mean "use the default" is met with ``invalid null
    reference in method 'FingerprintOptions_fp_type_set'``.
    """
    if value is not None:
        setattr(opts, key, value)


def _default_selector(value, default):
    """Fold a selector to lowercase, defaulting only when it is ``None``.

    ``None`` means unspecified and takes the default. Every other value -- an
    empty string included -- is something the caller chose, and folding it to
    the default would let this module give advice about a selection that was
    never made. See ``_set_if_given``: the real value reaches C++ either way,
    so the two have to agree on what "unspecified" means.
    """
    return default if value is None else value.lower()


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

# Mirrors normalize_family in src/comparisons/FingerprintComparison.cpp:124.
# The C++ constructor folds these spellings onto three generators; a rule
# keyed on the raw spelling would let an alias skip the family checks.
_FAMILY_ALIASES = {
    'morgan': 'morgan',
    'atom_pair': 'atom_pair',
    'atompair': 'atom_pair',
    'topological_atom_pair': 'atom_pair',
    'topological_torsions': 'topological_torsions',
    'topological_torsion': 'topological_torsions',
}

# Which families read each per-family option, and which option replaces it for
# a family that does not. Keyed on canonical families only.
_FAMILY_ONLY_KEYS = {
    'radius': ('morgan',),
    'min_distance': ('atom_pair',),
    'max_distance': ('atom_pair',),
    'torsion_atom_count': ('topological_torsions',),
}

_FAMILY_REPLACEMENT = {
    'morgan': 'radius',
    'atom_pair': 'min_distance/max_distance',
    'topological_torsions': 'torsion_atom_count',
}

_METRIC_ONLY_KEYS = {
    'p': 'minkowski',
    'tversky_alpha': 'tversky',
    'tversky_beta': 'tversky',
}

# Mirrors the on_fingerprint rows of METRIC_TABLE
# (src/comparisons/MetricTable.cpp:35). The descriptor-only rows --
# standardized_euclidean, seuclidean, mahalanobis -- are deliberately absent:
# C++ rejects them on this surface, so this module must treat them as
# unrecognized and stand aside.
_FINGERPRINT_METRICS = frozenset({
    'jaccard', 'tanimoto', 'dice', 'sokal_sneath', 'matching',
    'rogers_tanimoto', 'russell_rao', 'kulsinski', 'sokal_michener',
    'euclidean', 'manhattan', 'chebyshev', 'hamming', 'canberra',
    'bray_curtis', 'minkowski', 'tversky',
})

# Mirrors normalize_storage (src/comparisons/FingerprintComparison.cpp:153),
# which takes no aliases.
_FINGERPRINT_STORAGES = frozenset({'binary', 'count', 'sparse', 'sparse_count'})


def canonical_fingerprint_family(fp_type):
    """Return the canonical family for a spelling, or ``None`` if unrecognized.

    Mirrors ``normalize_family``
    (``src/comparisons/FingerprintComparison.cpp:124``). ``None`` means the C++
    constructor owns the error, and every Python rule keyed on the family has
    to stand aside so the caller sees it rather than advice about an unrelated
    option.
    """
    return _FAMILY_ALIASES.get(_default_selector(fp_type, 'morgan'))


def unrecognized_fingerprint_selector(fp_type, storage, metric):
    """Name the first selector the C++ constructor will reject, else ``None``.

    ``fp_type``, ``storage`` and ``metric`` are the three *authoritative*
    selectors: C++ owns whether each value exists, and this module only mirrors
    those tables. When any one of them is unrecognized, every advisory rule in
    this module has to stand aside, because the remedy an advisory message names
    -- a different metric, a different storage -- cannot make the call valid, and
    printing it hides the one thing the caller has to change.

    Standing aside on *any* unrecognized selector, rather than only on the axis
    a given rule reads, is deliberate. Two authoritative errors have no ordering
    between them, so the rule that is always right is to defer to C++ whenever
    C++ is going to speak.

    Mirror staleness is benign in one direction only. A value C++ has added but
    this module has not mirrored merely loses its advisory rules. A value C++
    has *removed* that this module still lists is not benign: it makes this
    function answer ``None`` for something C++ will reject, so an advisory rule
    fires ahead of the authoritative error -- the defect the mirrors exist to
    prevent. The mirror-drift tests in ``test_comparison_registry.py`` cover
    that direction.

    :param fp_type: Selected fingerprint family, or ``None`` for the default.
    :param storage: Selected storage, or ``None`` for the default.
    :param metric: Selected metric name, or ``None`` for the default.
    :returns: ``'fp_type'``, ``'storage'`` or ``'metric'`` for the first
        unrecognized selector, or ``None`` when all three are recognized.
    """
    if canonical_fingerprint_family(fp_type) is None:
        return 'fp_type'
    if _default_selector(storage, 'binary') not in _FINGERPRINT_STORAGES:
        return 'storage'
    if _default_selector(metric, 'tanimoto') not in _FINGERPRINT_METRICS:
        return 'metric'
    return None


def reject_inapplicable_fingerprint_kwargs(named, *, fp_type, storage, metric):
    """
    Reject an explicitly named option the rest of the configuration ignores.

    Only the Python surface can apply these rules. C++ receives a fully
    populated struct and cannot tell a deliberate value from a default, so it
    ignores unused fields instead -- a default-constructed ``FingerprintOptions``
    must stay a working Morgan configuration. Every message names both the
    option that would have been silently ignored and the one that replaces it.

    ``use_chirality`` applies to all four families and is never rejected.

    A value in any of ``fp_type``, ``storage`` or ``metric`` that this module
    does not recognize is left entirely to the C++ constructor: no rule here
    fires, so the caller sees the authoritative error instead of advice about an
    unrelated option. See ``unrecognized_fingerprint_selector``.

    :param named: Option names the caller passed explicitly.
    :param fp_type: Selected fingerprint family.
    :param storage: Selected storage.
    :param metric: Selected metric name.
    :raises TypeError: If a named option does not apply to this configuration.
    """
    named = set(named)
    if unrecognized_fingerprint_selector(fp_type, storage, metric) is not None:
        # An unrecognized selector is the C++ constructor's error to report.
        return
    spelling = _default_selector(fp_type, 'morgan')
    family = canonical_fingerprint_family(fp_type)
    store = _default_selector(storage, 'binary')
    metric_name = _default_selector(metric, 'tanimoto')

    if 'numbits' in named and store in _SPARSE_STORAGES:
        raise TypeError(
            f"numbits does not apply to storage={store!r}: a sparse "
            f"fingerprint keeps its family's own domain rather than folding to "
            f"a chosen width. Drop numbits, or use storage='binary' or "
            f"storage='count'.")

    for key, families in _FAMILY_ONLY_KEYS.items():
        if key in named and family not in families:
            raise TypeError(
                f"{key} does not apply to fp_type={spelling!r}; it belongs "
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
    # the caller named, which is exactly what the pops destroy. A ``None``
    # value is not a named option -- see _set_if_given -- so a rule keyed on an
    # option the caller declined to set must not fire.
    named = {key for key, value in kwargs.items() if value is not None}

    opts = FingerprintOptions()
    opts.similarity = similarity
    for key in _FINGERPRINT_KEYS:
        if key in kwargs:
            _set_if_given(opts, key, kwargs.pop(key))
    if kwargs:
        raise TypeError(
            f"Unknown kwargs for fingerprint comparison: {list(kwargs)}")

    unrecognized = unrecognized_fingerprint_selector(
        opts.fp_type, opts.storage, opts.metric)
    reject_inapplicable_fingerprint_kwargs(
        named, fp_type=opts.fp_type, storage=opts.storage, metric=opts.metric)

    # Tversky is symmetric only when alpha == beta, which the condensed pdist
    # form requires. Rejecting here names both parameters; the C++
    # ValidateForPDist message can only report the metric. Skipped whenever any
    # selector is unrecognized: equal weights or cdist cannot make an
    # unsupported fp_type, storage or metric valid, so that error goes first.
    # This guard runs after the rejector and needs its own precondition -- a
    # check that covers only the rejector leaves this path open, which is
    # exactly how the previous two rounds each missed half the defect.
    #
    # Weights outside [0, 1] are skipped for the same reason, and it is the
    # sharpest case: both remedies this message names -- equal weights, or
    # cdist -- still fail the range bound in validate_params
    # (src/comparisons/MetricTable.cpp:100), so firing here would name two
    # remedies that cannot work and bury the one that can. Written as a
    # positive range test so a NaN weight, which compares false against
    # everything, also falls through to C++.
    if (unrecognized is None and symmetric
            and opts.metric.lower() == 'tversky'
            and 0.0 <= opts.tversky_alpha <= 1.0
            and 0.0 <= opts.tversky_beta <= 1.0
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
        if st is not None:
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
        _set_if_given(opts, 'color_ff_type', kwargs.pop('color_ff_type'))
    if kwargs:
        raise TypeError(f"Unknown kwargs for rocs comparison: {list(kwargs)}")
    return _ROCSComparison(items, opts), "rocs"


def _build_superpose(items, similarity, kwargs, symmetric, *, default_method):
    """Build a :class:`SuperposeComparison` from keyword options."""
    opts = SuperposeOptions()
    opts.similarity = similarity
    method = kwargs.pop('method', None)
    if method is None:
        method = default_method
    elif not method and default_method is not None:
        # The sitehopper alias implies its own method, so an unusable value
        # falls back to the one the alias is named for -- the base chain's
        # behavior, preserved deliberately. On plain superpose there is no
        # implied method, so default_method is None and an unusable value
        # falls through to the ValueError below rather than being replaced by
        # a C++ default the caller never asked for.
        method = default_method
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
        if st is not None:
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
        _set_if_given(opts, 'predicate', kwargs.pop('predicate'))
    if 'ref_predicate' in kwargs:
        _set_if_given(opts, 'ref_predicate', kwargs.pop('ref_predicate'))
    if 'fit_predicate' in kwargs:
        _set_if_given(opts, 'fit_predicate', kwargs.pop('fit_predicate'))
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
