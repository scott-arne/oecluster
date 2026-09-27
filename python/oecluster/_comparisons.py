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

import numpy as np

from . import oecluster as _oecluster
from .oecluster import DescriptorComparison as _DescriptorComparison
from .oecluster import (
    DescriptorOptions,
    FingerprintOptions,
    MCSOptions,
    RMSDOptions,
    ROCSOptions,
    SuperposeOptions,
)
from .oecluster import FingerprintComparison as _FingerprintComparison
from .oecluster import MCSComparison as _MCSComparison
from .oecluster import RMSDComparison as _RMSDComparison
from .oecluster import ROCSComparison as _ROCSComparison
from .oecluster import SuperposeComparison as _SuperposeComparison

Builder = Callable[[list[Any], bool, dict[str, Any], bool], tuple[Any, str]]
Normalizer = Callable[[list[Any], dict[str, Any]],
                      tuple[list[Any], list[list[Any]]]]
Validator = Callable[[bool, dict[str, Any]], None]

_BUILDERS: dict[str, Builder] = {}
_NORMALIZERS: dict[str, Normalizer] = {}
_VALIDATORS: dict[str, Validator] = {}


def register_comparison(name, builder, normalizer=None, validator=None):
    """
    Register a comparison builder, and optional hooks, under a name.

    :param name: Lowercase comparison name accepted by :func:`oecluster.pdist`
                 and :func:`oecluster.cdist`.
    :param builder: Callable ``(items, similarity, kwargs, symmetric)``
                    returning ``(comparison_obj, comparison_name)``. It must
                    pop every keyword it consumes from ``kwargs`` and raise
                    :class:`TypeError` for anything left over.
    :param normalizer: Optional callable ``(items, kwargs)`` returning
                       ``(kept_items, excluded)``. It must read ``kwargs``
                       without popping, because the builder runs afterwards.
    :param validator: Optional callable ``(similarity, kwargs)`` raising on an
                      argument no input could make valid. It runs before the
                      normalizer and must also read ``kwargs`` without
                      popping. Register one whenever a later guard would
                      otherwise report a different problem first -- a
                      normalizer that empties the input, or the cutoff and
                      arrived-empty checks :func:`oecluster.cdist` makes after
                      validation. The builder remains the enforcing copy.
    """
    _BUILDERS[name] = builder
    if normalizer is not None:
        _NORMALIZERS[name] = normalizer
    if validator is not None:
        _VALIDATORS[name] = validator


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


def _default_selector(value, default, name):
    """Fold a selector to lowercase, defaulting only when it is ``None``.

    ``None`` means unspecified and takes the default. Every other value -- an
    empty string included -- is something the caller chose, and folding it to
    the default would let this module give advice about a selection that was
    never made. See ``_set_if_given``: the real value reaches C++ either way,
    so the two have to agree on what "unspecified" means.

    :param value: The selector as the caller gave it.
    :param default: The selector to use when ``value`` is ``None``.
    :param name: The argument name to blame in the failure message.
    :returns: ``default``, or ``value`` folded to lowercase.
    :raises TypeError: If ``value`` is neither ``None`` nor a string.
    """
    if value is None:
        return default
    # A non-string reaches here only from a direct call to one of this
    # module's rule functions: on the pdist and cdist paths the SWIG setters
    # have already refused anything but a string. Naming the argument still
    # beats the bare "'int' object has no attribute 'lower'" that the fold
    # would otherwise produce, which blames a method the caller never wrote.
    if not isinstance(value, str):
        raise TypeError(
            f"{name} must be a string or None, not "
            f"{type(value).__name__} ({value!r})")
    return value.lower()


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


def validate_request(comparison, similarity, kwargs):
    """
    Run a comparison's pre-normalization argument check, if it has one.

    Called by ``pdist`` and ``cdist`` ahead of :func:`normalize_items` so that
    an argument no input could make valid is reported before a normalizer gets
    the chance to fail on the item list instead. A normalizer that filters can
    empty the list, and the resulting complaint describes a filter the caller
    never asked for rather than the argument they must change.

    :param comparison: Comparison name.
    :param similarity: Whether the caller asked for similarities.
    :param kwargs: Comparison keyword options, read but never consumed.
    :raises ValueError: If the comparison name is unknown, or a validator
                        rejects an argument it owns outright.
    :raises TypeError: If a validator refuses a keyword. Among the reasons: a
                       name it does not know; a value whose type the
                       comparison's options object will not take; and a
                       Python-only option, one that never reaches an options
                       object, given a value of the wrong type.
    :raises RuntimeError: If a validator hands an option value to C++ and C++
                          refuses it.
    """
    validator = _VALIDATORS.get(_resolve_name(comparison))
    if validator is not None:
        validator(similarity, kwargs)


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


def canonical_fingerprint_family(fp_type):
    """Return the canonical family for a spelling, or ``None`` if unrecognized.

    Mirrors ``normalize_family``
    (``src/comparisons/FingerprintComparison.cpp:124``). Callers reach this only
    after the C++ constructor has accepted ``fp_type``, so ``None`` means the
    alias table above is behind C++ rather than that the spelling is invalid.
    The family-keyed rules stand aside in that case: a family this module cannot
    name has no advisory rules to apply.
    """
    return _FAMILY_ALIASES.get(_default_selector(fp_type, 'morgan', 'fp_type'))


def numbits_is_inapplicable(named, storage):
    """Report whether an explicitly named ``numbits`` has no meaning here.

    One predicate serves both call sites in ``_build_fingerprint`` -- the
    pre-construction reset and the post-construction rejection -- and they must
    agree exactly. If the reset could ever fire where the rejection does not, a
    caller-supplied ``numbits`` would be silently discarded on a call that
    succeeds, which is the one outcome the explicitness rules exist to prevent.
    Sharing the predicate makes that drift impossible rather than merely
    unlikely.

    :param named: Option names the caller passed explicitly.
    :param storage: Selected storage, or ``None`` for the default.
    :returns: ``True`` when ``numbits`` was named and the storage is sparse.
    """
    return ('numbits' in named
            and _default_selector(storage, 'binary', 'storage') in _SPARSE_STORAGES)


def reject_inapplicable_fingerprint_kwargs(named, *, fp_type, storage, metric):
    """
    Reject an explicitly named option the rest of the configuration ignores.

    Only the Python surface can apply these rules. C++ receives a fully
    populated struct and cannot tell a deliberate value from a default, so it
    ignores unused fields instead -- a default-constructed ``FingerprintOptions``
    must stay a working Morgan configuration. Every message names the option
    that would have been silently ignored together with the setting that made
    it inert, and offers a way out -- a replacement option to pass, or the
    setting to change.

    ``use_chirality`` applies to all three families and is never rejected.
    Six spellings fold onto those three; see ``_FAMILY_ALIASES``.

    Every rule here is advisory: it reports an option C++ accepts and silently
    ignores. Callers reach this function only after the C++ constructor has
    accepted the configuration, which is what keeps an advisory message from
    pre-empting an authoritative one. See ``_build_fingerprint``.

    :param named: Option names the caller passed explicitly.
    :param fp_type: Selected fingerprint family.
    :param storage: Selected storage.
    :param metric: Selected metric name.
    :raises TypeError: If a named option does not apply to this configuration.
    """
    named = set(named)
    spelling = _default_selector(fp_type, 'morgan', 'fp_type')
    family = canonical_fingerprint_family(fp_type)
    store = _default_selector(storage, 'binary', 'storage')
    metric_name = _default_selector(metric, 'tanimoto', 'metric')

    if numbits_is_inapplicable(named, storage):
        raise TypeError(
            f"numbits does not apply to storage={store!r}: a sparse "
            f"fingerprint keeps its family's own domain rather than folding to "
            f"a chosen width. Drop numbits, or use storage='binary'. "
            f"storage='count' folds to a width too, but only together with a "
            f"numeric metric such as 'manhattan': the default 'tanimoto' is a "
            f"bit-set metric and counted storage refuses it.")

    # A family C++ accepted but this module does not know yields family=None.
    # Standing aside then makes the alias table's staleness benign in both
    # directions: an unknown family loses its advisory rules and nothing else.
    # Firing instead would reject a key against a family we cannot name.
    for key, families in (() if family is None else _FAMILY_ONLY_KEYS.items()):
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

    # Morgan's sparse generators validate a num_bits they do not read: with
    # storage='sparse' and numbits=0 the constructor raises "Morgan num_bits
    # must be greater than zero", which names the one remedy that cannot work
    # -- a positive numbits is rejected below, and only dropping it succeeds.
    # atom_pair and topological_torsions accept numbits=0 under sparse storage,
    # which is what makes this an upstream OEFP inconsistency rather than a
    # rule to mirror here. Reset the field so C++ never validates a value it
    # will not read.
    #
    # No caller-supplied value is silently discarded: numbits_is_inapplicable
    # guards both this reset and the rejection below, so every call that
    # reaches this line goes on to raise by name -- from the constructor if
    # something else is wrong, from the rejector otherwise.
    if numbits_is_inapplicable(named, opts.storage):
        opts.numbits = FingerprintOptions().numbits

    # C++ speaks first, and that ordering is the whole design. The constructor
    # validates every selector and every metric parameter it owns; anything it
    # rejects raises here, so no advisory rule below can pre-empt an
    # authoritative error by naming a remedy that cannot make the call valid.
    # Seven fix rounds tried to get this right by having Python predict the
    # verdict, and each one missed a case the next one found.
    #
    # Construction is not free -- it computes the fingerprints, and the cost
    # scales with molecule size: roughly 10 ms for 400 molecules of 120 heavy
    # atoms, which is nearly all of what that pdist costs. Ordering it first is
    # still free on the success path, because pdist builds this same object
    # anyway and reuses it. What it costs is one wasted fingerprint pass before
    # an advisory raise, about 17 ms for 2000 small molecules. That is the
    # price of an authoritative message, and it is worth paying.
    comparison = _FingerprintComparison(items, opts)

    reject_inapplicable_fingerprint_kwargs(
        named, fp_type=opts.fp_type, storage=opts.storage, metric=opts.metric)

    # Tversky is symmetric only when alpha == beta, which the condensed pdist
    # form requires. Rejecting here names both parameters; the C++
    # ValidateForPDist message can only report the metric. No precondition is
    # needed any more: an unrecognized metric or an out-of-range weight has
    # already raised at construction above.
    if (symmetric and opts.metric.lower() == 'tversky'
            and opts.tversky_alpha != opts.tversky_beta):
        raise ValueError(
            "pdist requires a symmetric metric, but tversky with "
            f"tversky_alpha={opts.tversky_alpha} and "
            f"tversky_beta={opts.tversky_beta} is asymmetric. Use equal "
            "alpha and beta, or compute a rectangular result with cdist.")

    return comparison, "fingerprint"


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


def _string_vector(values):
    """Convert a Python sequence into a SWIG ``StringVector``."""
    vector = _oecluster.StringVector()
    for value in values:
        vector.push_back(str(value))
    return vector


def _double_vector(values):
    """Convert a Python sequence into a SWIG ``DoubleVector``."""
    vector = _oecluster.DoubleVector()
    for value in values:
        vector.push_back(float(value))
    return vector


def _empty_option_error(name):
    """
    Build the refusal for a descriptor option supplied as an empty sequence.

    Emptiness is the native layer's own encoding of "no override":
    ``DescriptorComparison.cpp`` keys ``has_override`` and the
    ``variances``/``inverse_covariance`` pairing rules off ``.empty()``, and a
    selection option left empty resolves to the same default schema as one
    never set. An empty sequence therefore arrives in C++ indistinguishable
    from an omitted option, so Python is the only place the two can still be
    told apart.

    Returned rather than raised so the caller chooses *when* it fires. This is
    an advisory rule of this layer's own and it ranks last among the checks
    that need no molecules, so a caller with an authoritative verdict to run
    first has to be able to hold it back.

    :param name: The keyword option name, used in the message.
    :returns: The ``ValueError`` to raise.
    """
    return ValueError(
        f"{name}= was given as an empty sequence, which the descriptor "
        f"layer cannot tell apart from never passing it. Pass the values "
        f"you want, or omit {name}= to take the default.")


_DESCRIPTOR_KEYS = (
    'sources',
    'columns',
    'groups',
    'metric',
    'variances',
    'inverse_covariance',
    'missing',
    'p',
)

# The subset carrying a sequence, which is the subset an empty value can be
# passed to. Each name is also the ``DescriptorOptions`` field name.
_DESCRIPTOR_SEQUENCE_KEYS = (
    'sources',
    'columns',
    'groups',
    'variances',
    'inverse_covariance',
)


def _descriptor_options_and_empties(kwargs):
    """
    Build a :class:`DescriptorOptions` and report its empty sequence options.

    Conversion is separated from the refusal because deciding a sequence is
    empty requires converting it, while the refusal itself must wait: it ranks
    last among the checks that read no molecules, and raising during the
    conversion would put it ahead of every one of them. Callers with a verdict
    to give precedence to run it against the returned options and raise
    afterwards; callers with none use :func:`descriptor_options`.

    An empty vector is assigned rather than skipped. Every sequence field of
    ``DescriptorOptions`` default-constructs empty, so the two leave the object
    in the same state, and the validator therefore rules on exactly the request
    the caller made.

    :param kwargs: Comparison keyword options, read but never consumed.
    :returns: A ``(options, empty)`` pair, where ``empty`` names the sequence
        options that were supplied as empty sequences, in signature order.
    """
    opts = DescriptorOptions()
    if kwargs.get('sources') is not None:
        opts.sources = _string_vector(kwargs['sources'])
    if kwargs.get('columns') is not None:
        opts.columns = _string_vector(kwargs['columns'])
    if kwargs.get('groups') is not None:
        opts.groups = _string_vector(kwargs['groups'])
    if kwargs.get('metric') is not None:
        opts.metric = str(kwargs['metric'])
    if kwargs.get('variances') is not None:
        opts.variances = _double_vector(kwargs['variances'])
    if kwargs.get('inverse_covariance') is not None:
        opts.inverse_covariance = _double_vector(
            np.asarray(kwargs['inverse_covariance'], dtype=np.float64).ravel())
    if kwargs.get('missing') is not None:
        opts.missing = str(kwargs['missing']).lower()
    if kwargs.get('p') is not None:
        opts.p = float(kwargs['p'])

    # Read back off the options rather than off kwargs, so that the 2-D
    # inverse_covariance form is judged by what the ravel produced.
    empty = [name for name in _DESCRIPTOR_SEQUENCE_KEYS
             if kwargs.get(name) is not None and len(getattr(opts, name)) == 0]
    return opts, empty


def descriptor_options(kwargs):
    """
    Build a :class:`DescriptorOptions` from keyword options.

    Reads ``kwargs`` without consuming it, because the normalizer runs before
    the builder and both need the same options.

    Refuses an empty sequence option on the spot, which is correct only for a
    caller that has no argument-level verdict of its own left to give. The
    callers that do -- ``_validate_descriptor``, and the prebuilt
    ``DescriptorComparison`` -- use
    :func:`_descriptor_options_and_empties` and raise after that verdict.

    :param kwargs: Comparison keyword options.
    :returns: A populated ``DescriptorOptions``.
    :raises ValueError: If any sequence-valued option was passed as an empty
        sequence, which C++ cannot distinguish from an omitted option.
    """
    opts, empty = _descriptor_options_and_empties(kwargs)
    if empty:
        raise _empty_option_error(empty[0])
    return opts


def _normalize_descriptor(items, kwargs):
    """Drop molecules missing a descriptor value under the complete-case policy."""
    # Only ``None`` is unspecified, as in ``_default_selector``. Folding any
    # other value onto the default would run the complete-case filter for a
    # policy the caller never chose, and that filter can empty the item list
    # and raise before C++ ever reports the value as invalid.
    missing = kwargs.get('missing')
    policy = 'complete_case' if missing is None else str(missing).lower()
    if policy != 'complete_case':
        return items, []

    excluded = sorted(
        int(index) for index in _oecluster.descriptor_excluded_indices(
            items, descriptor_options(kwargs)))
    if not excluded:
        return items, []

    dropped = set(excluded)
    kept = [item for idx, item in enumerate(items) if idx not in dropped]
    return kept, [[idx, "missing-descriptor"] for idx in excluded]


def _validate_descriptor(similarity, kwargs):
    """Reject descriptor arguments before any molecule is read.

    Three of the rules are this layer's own -- no similarity form, no unknown
    keyword, no empty sequence option. The option values go to
    ``validate_descriptor_options``, which does not answer for every mistake
    C++ can decide without molecules; the comment on that call says which it
    leaves downstream.

    The rules run in the order the caller has to fix them. The two names come
    first because neither depends on a value being readable at all. The
    emptiness refusal comes last, after the C++ verdicts: it is this layer's
    own advisory rule, and a caller who also misspelled the metric must be told
    about the metric, which no edit to the empty sequence would resolve. What
    that ordering does not reach is a source, column or group name on a request
    with no override, which the validator resolves no schema for and so leaves
    downstream; the guard therefore does precede *that* report.

    Registered as the comparison's validator so the whole band also runs before
    ``_normalize_descriptor``, whose complete-case filter can empty the item
    list and have ``pdist`` refuse the shape instead of the argument. Reads
    ``kwargs`` without popping: the normalizer and then the builder still need
    every key, and ``build_comparison`` treats leftovers as a builder bug.

    :param similarity: Whether the caller asked for similarities.
    :param kwargs: Comparison keyword options, read but never consumed.
    :raises ValueError: If similarities were requested, or if a sequence-valued
        option was passed as an empty sequence.
    :raises TypeError: If any keyword option is not a descriptor option.
    :raises RuntimeError: If C++ refuses an option value outright.
    """
    if similarity:
        raise ValueError(
            "the descriptor comparison has no similarity form: it measures "
            "distance in descriptor space. Use similarity=False.")

    unknown = [key for key in kwargs if key not in _DESCRIPTOR_KEYS]
    if unknown:
        raise TypeError(
            f"Unknown kwargs for descriptor comparison: {unknown}")

    # Only C++ can rule on a value: the metric table lives in a private header
    # and is not reachable from here, so a Python copy of it would be a second
    # source of truth that drifts. ``validate_descriptor_options`` decides with
    # no molecules in hand, which is what lets it run before the
    # filter -- including the rules matching ``variances`` and
    # ``inverse_covariance`` to the selected columns, since the schema comes
    # from ``sources`` alone. One of those is the semidefinite verdict on an
    # ``inverse_covariance`` override, which costs an eigendecomposition and is
    # paid again by ``descriptor_excluded_indices`` and by the constructor,
    # neither of which reuses this one. It stops short of the source, column
    # and group names on a request with no override: resolving those needs a
    # schema, and with no override nothing on that path asks for one to be
    # built, so C++ reports them further down instead. Everything else left
    # downstream reads the input; among those, the minimum input size, the
    # complete-case row check, the refusal when every selected column is
    # constant. None of it can be decided here.
    opts, empty = _descriptor_options_and_empties(kwargs)
    _oecluster.validate_descriptor_options(opts)

    # Last in the band, so none of the verdicts above is pre-empted. Building
    # the options is what discovers the emptiness, which is why the refusal is
    # carried this far rather than raised where it was found.
    #
    # Not a repeat of the builder's copy. ``cdist`` refuses an input set that
    # arrived empty between here and the builder, so deleting this raise makes
    # ``cdist([], b, "descriptor", columns=[])`` report the shape instead.
    if empty:
        raise _empty_option_error(empty[0])


def _build_descriptor(items, similarity, kwargs, symmetric):
    """Build a :class:`DescriptorComparison` from keyword options."""
    # ``build_comparison`` is a module-level entry point, reachable without the
    # ``validate_request`` call ``pdist`` and ``cdist`` make first, so the
    # builder stays the enforcing copy.
    _validate_descriptor(similarity, kwargs)

    opts = descriptor_options(kwargs)
    for key in _DESCRIPTOR_KEYS:
        kwargs.pop(key, None)

    return _DescriptorComparison(items, opts), "descriptor"


# expand_conformers is a Python-only normalizer key, never reaching RMSDOptions.
_RMSD_KEYS = ('overlay', 'automorph', 'heavy_only', 'expand_conformers')


def rmsd_options(kwargs):
    """
    Build an :class:`RMSDOptions` from keyword options.

    Reads ``kwargs`` without consuming it.

    :param kwargs: Comparison keyword options.
    :returns: A populated ``RMSDOptions``.
    """
    opts = RMSDOptions()
    _set_if_given(opts, 'overlay', kwargs.get('overlay'))
    _set_if_given(opts, 'automorph', kwargs.get('automorph'))
    _set_if_given(opts, 'heavy_only', kwargs.get('heavy_only'))
    return opts


def _expand_conformers_requested(kwargs):
    """
    Read the ``expand_conformers`` flag, refusing a value that is not boolean.

    ``numpy.bool_`` counts as boolean here, matching
    :func:`oecluster._gate.require_metric`. The three C++ options do not accept
    it, because their SWIG setters take a plain ``bool``.

    :param kwargs: Comparison keyword options, read but never consumed.
    :returns: Whether multi-conformer molecules should be expanded. Omitting
        the flag, or passing ``None``, means yes.
    :raises TypeError: If ``expand_conformers`` is given a value other than
        ``None``, a bool, or a ``numpy.bool_``.
    """
    expand = kwargs.get('expand_conformers')
    if expand is None:
        return True
    if not isinstance(expand, (bool, np.bool_)):
        raise TypeError(
            "expand_conformers must be True or False, "
            f"not {type(expand).__name__} ({expand!r}). "
            "Accepting anything else would let a value like 'false' or 0 "
            "settle the expansion by truthiness rather than by what the "
            "caller wrote.")
    return bool(expand)


def _normalize_rmsd(items, kwargs):
    """
    Expand multi-conformer molecules into one molecule per pose.

    ``OERMSD`` compares a single pose against a single pose, so each conformer
    becomes its own item, titled ``"<title>:conf<n>"``. Multi-conformer inputs
    are expanded into copies; single-conformer inputs pass through unchanged.
    Either way the caller's molecules are untouched, because the comparison
    takes its own snapshot.
    """
    if not _expand_conformers_requested(kwargs):
        return items, []

    from openeye import oechem

    expanded = []
    for idx, item in enumerate(items):
        num_confs = item.NumConfs() if hasattr(item, "NumConfs") else 1
        if num_confs <= 1:
            expanded.append(item)
            continue
        title = item.GetTitle() or f"mol_{idx}"
        for conf_index, conf in enumerate(item.GetConfs()):
            single = oechem.OEMol()
            oechem.OEAddMols(single, oechem.OEGraphMol(item))
            src = oechem.OEFloatArray(conf.GetMaxAtomIdx() * 3)
            conf.GetCoords(src)
            dst = oechem.OEFloatArray(single.GetMaxAtomIdx() * 3)
            # strict=True because a length mismatch would otherwise be silent:
            # OEFloatArray zero-initializes, so the unwritten atoms would sit at
            # the origin and produce a plausible-looking distance.
            for src_atom, dst_atom in zip(
                    conf.GetAtoms(), single.GetAtoms(), strict=True):
                for axis in range(3):
                    dst[dst_atom.GetIdx() * 3 + axis] = src[
                        src_atom.GetIdx() * 3 + axis]
            single.SetCoords(dst)
            single.SetTitle(f"{title}:conf{conf_index}")
            expanded.append(single)
    return expanded, []


def _validate_rmsd(similarity, kwargs):
    """Reject RMSD arguments before any molecule is read.

    Registered as the comparison's validator so it runs ahead of
    ``_normalize_rmsd`` and, in :func:`oecluster.cdist`, ahead of the cutoff
    and arrived-empty guards. Without it,
    ``cdist(a, b, "rmsd", similarity=True, cutoff=0.5)`` would report the
    cutoff and name a remedy that cannot work, instead of naming the argument
    the caller has to fix.

    :param similarity: Whether the caller asked for similarities.
    :param kwargs: Comparison keyword options, read but never consumed.
    :raises ValueError: If similarities were requested.
    :raises TypeError: If any keyword option is not an RMSD option, if
        ``expand_conformers`` is given a value other than ``None``, a bool, or
        a ``numpy.bool_``, or if an option value is of a type ``RMSDOptions``
        will not take.
    """
    if similarity:
        raise ValueError(
            "the rmsd comparison has no similarity form: RMSD is a distance "
            "in angstroms. Use similarity=False.")
    unknown = [key for key in kwargs if key not in _RMSD_KEYS]
    if unknown:
        raise TypeError(f"Unknown kwargs for rmsd comparison: {unknown}")
    _expand_conformers_requested(kwargs)
    # Building the options here is what puts the option-value refusals ahead of
    # cdist's guards too. Without it, cdist([], b, "rmsd", overlay="false")
    # reports the empty input set -- the same misdirection the similarity check
    # above exists to prevent.
    rmsd_options(kwargs)


def _build_rmsd(items, similarity, kwargs, symmetric):
    """Build an :class:`RMSDComparison` from keyword options."""
    # ``build_comparison`` is a module-level entry point, reachable without the
    # ``validate_request`` call ``pdist`` and ``cdist`` make first, so the
    # builder stays the enforcing copy.
    _validate_rmsd(similarity, kwargs)

    opts = rmsd_options(kwargs)
    for key in _RMSD_KEYS:
        kwargs.pop(key, None)

    return _RMSDComparison(items, opts), "rmsd"


# similarity is deliberately absent: it is not a keyword option but its own
# positional argument on the builder, and MCSOptions leaves it False.
_MCS_KEYS = ('search_mode', 'match_level', 'max_matches')

# The scoped C++ enums flatten to module-level names in the generated wrapper,
# the way SuperposeMethod does above.
_MCS_SEARCH_MODES = {
    'approximate': _oecluster.MCSSearchMode_Approximate,
    'exhaustive': _oecluster.MCSSearchMode_Exhaustive,
}

_MCS_MATCH_LEVELS = {
    'default': _oecluster.MCSMatchLevel_Default,
    'exact': _oecluster.MCSMatchLevel_Exact,
    'loose': _oecluster.MCSMatchLevel_Loose,
}


def mcs_options(kwargs):
    """
    Build an :class:`MCSOptions` from keyword options.

    Reads ``kwargs`` without consuming it, so the validator and the builder can
    both call it. ``similarity`` is never set here, because it is not an
    ``_MCS_KEYS`` member: each of the two construction sites assigns it from its
    own argument.

    :param kwargs: Comparison keyword options.
    :returns: A populated ``MCSOptions``.
    :raises TypeError: If ``search_mode`` or ``match_level`` is neither a string
        nor ``None``, if ``max_matches`` is a bool, or if ``max_matches`` is of
        a type ``MCSOptions`` will not take.
    :raises ValueError: If ``search_mode`` or ``match_level`` names a mode the
        comparison does not have. Names are resolved by lookup, never by
        truthiness.
    """
    opts = MCSOptions()
    search_mode = kwargs.get('search_mode')
    if search_mode is not None:
        key = _default_selector(search_mode, 'approximate', 'search_mode')
        if key not in _MCS_SEARCH_MODES:
            raise ValueError(f"Unknown MCS search mode: {search_mode}")
        opts.search_mode = _MCS_SEARCH_MODES[key]
    match_level = kwargs.get('match_level')
    if match_level is not None:
        key = _default_selector(match_level, 'default', 'match_level')
        if key not in _MCS_MATCH_LEVELS:
            raise ValueError(f"Unknown MCS match level: {match_level}")
        opts.match_level = _MCS_MATCH_LEVELS[key]
    # bool is an int subclass, so the SWIG ``unsigned int`` typemap -- the
    # type check every other option here leans on -- reads True as a budget
    # of 1. That is the most destructive value the option has: it moves
    # morphine against penicillin G from 0.717949 to 0.958333. Refuse it
    # explicitly rather than let the typemap wave it through.
    max_matches = kwargs.get('max_matches')
    if isinstance(max_matches, bool):
        raise TypeError(
            "max_matches must be an integer, not a bool: True would be read "
            "as a match budget of 1, which silently corrupts scores")
    _set_if_given(opts, 'max_matches', max_matches)
    return opts


def _validate_mcs(similarity, kwargs):
    """Reject MCS arguments before any molecule is read.

    Registered as the comparison's validator so that, in :func:`oecluster.cdist`,
    an unknown keyword or an unusable ``search_mode`` is reported ahead of the
    cutoff and arrived-empty guards, which would otherwise name a remedy that
    cannot fix it.

    :param similarity: Whether the caller asked for similarities. MCS supports
        both orientations -- bond Tanimoto is natively a similarity and the
        distance is the derived form -- so nothing is refused on it. The
        parameter is here because the validator contract passes it.
    :param kwargs: Comparison keyword options, read but never consumed.
    :raises TypeError: If any keyword option is not an MCS option, if
        ``max_matches`` is a bool, or if an option value is of a type
        ``MCSOptions`` will not take.
    :raises ValueError: If ``search_mode`` or ``match_level`` names a mode the
        comparison does not have.
    """
    unknown = [key for key in kwargs if key not in _MCS_KEYS]
    if unknown:
        raise TypeError(f"Unknown kwargs for mcs comparison: {unknown}")
    # Building the options here is what puts the enum-name and option-value
    # refusals ahead of cdist's guards too.
    mcs_options(kwargs)


def _build_mcs(items, similarity, kwargs, symmetric):
    """Build an :class:`MCSComparison` from keyword options."""
    # ``build_comparison`` is a module-level entry point, reachable without the
    # ``validate_request`` call ``pdist`` and ``cdist`` make first, so the
    # builder stays the enforcing copy.
    _validate_mcs(similarity, kwargs)

    opts = mcs_options(kwargs)
    # mcs_options never sets this, because similarity is not an _MCS_KEYS
    # member. The top-level MCSComparison wrapper carries its own copy of this
    # line for the same reason. Assigned raw, as every sibling builder does, so
    # the SWIG bool typemap refuses a string the caller meant as false.
    opts.similarity = similarity
    for key in _MCS_KEYS:
        kwargs.pop(key, None)

    return _MCSComparison(items, opts), "mcs"


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
register_comparison("descriptor", _build_descriptor, _normalize_descriptor,
                    _validate_descriptor)
register_comparison("rmsd", _build_rmsd, _normalize_rmsd, _validate_rmsd)
# No normalizer: MCS is topological, so multi-conformer input passes through
# untouched and 3D coordinates are irrelevant.
register_comparison("mcs", _build_mcs, None, _validate_mcs)
