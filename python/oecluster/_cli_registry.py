"""Algorithm registry for the command line, derived from the clustering roster.

The roster carries no type annotations, so every fact here comes from a
parameter's kind and its default value, plus a small declared table for what
introspection cannot reach. Deriving rather than hardcoding means a new roster
entry appears on the command line without touching this module.
"""
import difflib
import inspect
import math

from ._parameter_selection import _roster

#: First-parameter names that mean "this entry can consume a distance matrix".
_MATRIX_KINDS = ("distance_matrix", "items")

#: Routed by their own flags, so never settable through ``--set``.
_FLAG_OPTIONS = {"num_threads": "--threads", "allow_nonmetric": "--allow-nonmetric"}

#: A tuning knob for the parallel passes, not a clustering parameter; the
#: library's own default is what callers should use.
_HIDDEN = {"chunk_size"}

#: Required options have no default to infer a type from, so each declares
#: one. Leaving them to a permissive fallback let ``threshold=true`` become
#: True, which butina then silently read as 1.0.
_REQUIRED_TYPES = {
    "butina.threshold": float,
    "dbscan.eps": float,
    "sphere_exclusion.threshold": float,
    "jarvis_patrick.kmin": int,
}

#: Meaningless when the input is already a matrix: they describe how to build
#: one from raw items.
_EXCLUDED_FOR_ITEMS = {"comparison", "similarity"}

#: The library coerces these with a bare ``int()`` (``__init__.py:2884``), so a
#: float would be silently truncated rather than refused.
_INTEGER = {
    "agglomerative.n_clusters", "hdbscan.max_cluster_size",
    "hdbscan.min_cluster_size", "hdbscan.min_samples", "jarvis_patrick.k",
    "jarvis_patrick.kmin", "k_medoids.max_iterations", "k_medoids.n_clusters",
    "leiden.k", "leiden.n_iterations",
}

#: Default ``None`` in the signature, but required when the input is a matrix.
_CONDITIONALLY_REQUIRED = {"jarvis_patrick.k", "leiden.k"}

#: Not expressible in a ``key=value`` grammar.
_SEQUENCE = {"k_medoids.initial_medoids"}


class Entry:
    """One algorithm's command-line schema."""

    def __init__(self, name, kind, required, optional, sequence, nonmetric, summary):
        self.name = name
        self.kind = kind
        self.eligible = kind in _MATRIX_KINDS
        self.required = required
        self.optional = optional
        self.sequence = sequence
        self.accepts_nonmetric = nonmetric
        self.summary = summary

    def known(self):
        """:returns: Every option name this algorithm accepts from ``--set``."""
        return sorted([*self.required, *self.optional])


def _summary(fn):
    doc = inspect.getdoc(fn) or ""
    return doc.splitlines()[0].strip() if doc else ""


def build():
    """Introspect the roster into ``{name: Entry}``.

    :returns: Mapping of algorithm name to :class:`Entry`.
    """
    registry = {}
    for name, fn in _roster().items():
        params = list(inspect.signature(fn).parameters.values())
        kind = params[0].name
        required, optional, sequence, nonmetric = [], {}, [], False
        for param in params[1:]:
            # Variadics are dropped by KIND, not by name: jarvis_patrick,
            # leiden and sphere_exclusion all carry **kwargs, and under the
            # "no default means required" rule a missed one becomes a bogus
            # required option.
            if param.kind in (param.VAR_POSITIONAL, param.VAR_KEYWORD):
                continue
            if param.name == "allow_nonmetric":
                nonmetric = True
                continue
            if param.name in _FLAG_OPTIONS or param.name in _HIDDEN or (
                    kind == "items" and param.name in _EXCLUDED_FOR_ITEMS):
                continue
            key = f"{name}.{param.name}"
            if key in _SEQUENCE:
                sequence.append(param.name)
            elif param.default is param.empty or key in _CONDITIONALLY_REQUIRED:
                required.append(param.name)
            else:
                optional[param.name] = param.default
        registry[name] = Entry(name, kind, required, optional, sequence,
                               nonmetric, _summary(fn))
    return registry


def _as_int(option, raw):
    """:raises ValueError: If ``raw`` is not a finite exact integer."""
    try:
        number = float(raw)
    except ValueError:
        raise ValueError(f"{option} must be an integer, got {raw!r}") from None
    # int(float("inf")) raises OverflowError, and int(float("nan")) raises
    # ValueError; both would escape as something the CLI never promised.
    if not math.isfinite(number) or number != int(number):
        raise ValueError(f"{option} must be an integer, got {raw!r}")
    return int(number)


def _as_float(option, raw):
    """:raises ValueError: If ``raw`` is not a finite number."""
    try:
        number = float(raw)
    except ValueError:
        raise ValueError(f"{option} must be a number, got {raw!r}") from None
    # An infinite or NaN threshold, epsilon or alpha is not something the
    # native layer rejects with a message a user can act on, so every float
    # path refuses it here rather than only the declared required ones.
    if not math.isfinite(number):
        raise ValueError(f"{option} must be finite, got {raw!r}")
    return number


def option_type(algorithm, option, entry):
    """:returns: A short type name for the schema table."""
    key = f"{algorithm}.{option}"
    if key in _REQUIRED_TYPES:
        return _REQUIRED_TYPES[key].__name__
    if key in _INTEGER:
        return "int"
    default = entry.optional.get(option)
    if isinstance(default, bool):
        return "bool"
    if isinstance(default, str):
        return "str"
    if isinstance(default, int):
        return "int"
    if isinstance(default, float):
        return "float"
    return "number"


def coerce(algorithm, option, raw, entry):
    """Convert one ``--set`` value using the option's declared or default type.

    :raises ValueError: If the value does not fit the option's type.
    """
    key = f"{algorithm}.{option}"
    if key in _REQUIRED_TYPES:
        want = _REQUIRED_TYPES[key]
        if want is int:
            return _as_int(option, raw)
        return _as_float(option, raw)
    if key in _INTEGER:
        return _as_int(option, raw)
    default = entry.optional.get(option)
    if isinstance(default, bool):
        if raw.lower() in ("true", "false"):
            return raw.lower() == "true"
        raise ValueError(f"{option} must be true or false, got {raw!r}")
    if isinstance(default, str):
        return raw
    if isinstance(default, int):
        try:
            return int(raw)
        except ValueError:
            raise ValueError(f"{option} must be an integer, got {raw!r}") from None
    if isinstance(default, float):
        return _as_float(option, raw)
    # Only options whose default is None reach here, and their type is
    # genuinely unknown, so the permissive pass stays -- but it never sees a
    # required option, which is where guessing did real damage.
    for convert in (int, float):
        try:
            value = convert(raw)
        except ValueError:
            continue
        # Unknown type or not, a value that parses as a number is a number,
        # and the finiteness rule applies to it as much as to a declared one.
        if convert is float and not math.isfinite(value):
            raise ValueError(f"{option} must be finite, got {raw!r}")
        return value
    if raw.lower() in ("true", "false"):
        return raw.lower() == "true"
    return raw


def resolve(algorithm, assignments, registry, *, swept=None):
    """Validate ``--set`` assignments and return the options mapping.

    :param algorithm: Roster name, already known to exist.
    :param assignments: Raw ``KEY=VALUE`` strings.
    :param registry: Mapping from :func:`build`.
    :param swept: Option name supplied per value by a sweep, if any.
    :returns: Coerced option mapping for :class:`ClusteringSpec`.
    :raises ValueError: For an unknown, excluded, duplicated or malformed key,
        a bad value, or a missing required option.
    """
    entry = registry[algorithm]
    options = {}
    for item in assignments:
        if "=" not in item:
            raise ValueError(f"--set expects KEY=VALUE, got {item!r}")
        key, raw = item.split("=", 1)
        # Both sides are stripped: a shell-quoted `--set "linkage = ward"`
        # would otherwise carry the space into the value, where a string
        # option smuggles it through and a bool option refuses outright.
        key, raw = key.strip(), raw.strip()
        if key in _FLAG_OPTIONS:
            raise ValueError(f"set {key} with {_FLAG_OPTIONS[key]}, not --set")
        if key in entry.sequence:
            raise ValueError(
                f"{key} takes a sequence; use the Python API for it")
        if key not in entry.known():
            hint = difflib.get_close_matches(key, entry.known(), n=1)
            suffix = f"; did you mean {hint[0]!r}?" if hint else ""
            raise ValueError(f"{algorithm} has no option {key!r}{suffix}")
        if key in options:
            raise ValueError(f"--set {key} given twice")
        if swept is not None and key == swept:
            raise ValueError(f"{key} is the swept parameter; drop it from --set")
        options[key] = coerce(algorithm, key, raw, entry)
    for name in entry.required:
        if name not in options and name != swept:
            raise ValueError(f"{algorithm} requires --set {name}=…")
    return options
