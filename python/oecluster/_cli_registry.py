"""Algorithm registry for the command line, derived from the clustering roster.

The roster carries no type annotations, so every fact here comes from a
parameter's kind and its default value, plus a small declared table for what
introspection cannot reach. Deriving rather than hardcoding means a new roster
entry appears on the command line without touching this module.
"""
import decimal
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
#: one from raw items, or bound the graph built from a comparison.
#: ``distance_matrix`` is the keyword alias for the matrix itself.
_EXCLUDED_FOR_ITEMS = {"comparison", "distance_matrix", "max_graph_bytes",
                       "similarity"}

#: The library coerces these with a bare ``int()`` (dbscan's
#: ``options.min_samples = int(min_samples)``), so a float would be silently truncated rather than refused.
_INTEGER = {
    "agglomerative.n_clusters", "hdbscan.max_cluster_size",
    "hdbscan.min_cluster_size", "hdbscan.min_samples", "jarvis_patrick.k",
    "jarvis_patrick.kmin", "k_medoids.max_iterations", "k_medoids.n_clusters",
    "leiden.k", "leiden.n_iterations",
}

#: Float-valued but defaulting to ``None``, so there is no default to infer a
#: type from. Declared rather than guessed: guessing from the value shape read
#: ``distance_threshold=ward`` as a string and ``=true`` as a bool, which
#: agglomerative then takes as 1.0.
_FLOAT = {"agglomerative.distance_threshold"}

#: Default ``None``, or the library's "not given" sentinel, in the signature,
#: but required when the input is a matrix.
_CONDITIONALLY_REQUIRED = {"butina.threshold", "dbscan.eps", "jarvis_patrick.k",
                           "leiden.k"}

#: Not expressible in a ``key=value`` grammar.
_SEQUENCE = {"k_medoids.initial_medoids"}

#: Widest magnitude any native integer parameter can hold. The clustering
#: signatures use ``size_t`` for the counts, ``int64_t`` for
#: ``leiden.n_iterations`` and ``uint64_t`` for ``leiden.seed``, so no single
#: native type covers them all; ``uint64``'s maximum is the widest of the
#: three and is therefore the only bound that refuses nothing the library
#: accepts. The per-parameter limits are left to the library, which already
#: checks each one and names it ("n_iterations must be between -1 and ...").
#: The bound here exists to keep a mistyped exponent from being built into an
#: integer at all, not to second-guess those checks.
_MAGNITUDE_MAX = 2**64 - 1

#: ``Decimal.adjusted()`` of ``_MAGNITUDE_MAX``: anything larger needs more
#: than 20 digits and is refused before ``int()`` sees it.
_MAGNITUDE_DIGITS = 19


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
    # Decimal, not float: the value is decimal text, and a binary float
    # cannot hold it faithfully. Past 2**53 consecutive integers collapse
    # onto the same float, so 2**53+1 came back one short, and 1e-324
    # underflows to 0.0, which reads as the integer 0 rather than as the
    # fraction it is. Decimal keeps the written value exact.
    try:
        number = decimal.Decimal(raw)
    except decimal.InvalidOperation:
        raise ValueError(f"{option} must be an integer, got {raw!r}") from None
    # Decimal("inf") and Decimal("nan") parse, so finiteness is checked here
    # rather than left to int(), which raises OverflowError on the first and
    # an opaque ValueError on the second.
    if not number.is_finite() or number != number.to_integral_value():
        raise ValueError(f"{option} must be an integer, got {raw!r}")
    # The magnitude is bounded before int() is allowed to materialize the
    # digits. int(Decimal("1e999999999")) builds a billion-digit integer and
    # hangs, and anything past 4300 digits trips CPython's conversion limit,
    # whose message tells the user about sys.set_int_max_str_digits. Both
    # reach a user who merely mistyped a number. adjusted() is the base-10
    # exponent of the leading digit, so it bounds the size without building
    # anything; past 20 digits the value cannot fit any native parameter.
    if number.adjusted() > _MAGNITUDE_DIGITS or abs(int(number)) > _MAGNITUDE_MAX:
        raise ValueError(f"{option} is out of range, got {raw!r}")
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
    if key in _FLOAT:
        return "float"
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
    if key in _FLOAT:
        return _as_float(option, raw)
    default = entry.optional.get(option)
    if isinstance(default, bool):
        if raw.lower() in ("true", "false"):
            return raw.lower() == "true"
        raise ValueError(f"{option} must be true or false, got {raw!r}")
    if isinstance(default, str):
        return raw
    if isinstance(default, int):
        # Same path as a declared integer: an option is no less an integer
        # for having inferred its type from its default, and a bare int()
        # here would take neither the exactness nor the range guard.
        return _as_int(option, raw)
    if isinstance(default, float):
        return _as_float(option, raw)
    # Only an option defaulting to None reaches here, and every one the
    # roster has is numeric and declared in a table above, so this is a
    # registry gap rather than user error. It used to guess from the value
    # shape, which accepted `distance_threshold=ward` as a string and
    # `=true` as a bool. Refusing keeps the guess from coming back; the
    # companion test over None defaults turns a new undeclared option into a
    # suite failure rather than a message a user has to decipher.
    raise ValueError(f"{option} has no declared type; declare one in "
                     f"_cli_registry for {key}")


def resolve(algorithm, assignments, registry, *, swept=None, remedy="--set"):
    """Validate ``--set`` assignments and return the options mapping.

    :param algorithm: Roster name, already known to exist.
    :param assignments: Raw ``KEY=VALUE`` strings.
    :param registry: Mapping from :func:`build`.
    :param swept: Option name supplied per value by a sweep, if any.
    :param remedy: The syntax the caller actually accepts an option in,
        named in every message that tells the user what to do about one.
        ``consensus`` validates each ``--member``'s options through here and
        refuses ``--set`` in that mode, so a member missing a required
        option was answered "butina requires --set threshold=…" -- advice
        the same command rejects.
    :returns: Coerced option mapping for :class:`ClusteringSpec`.
    :raises ValueError: For an unknown, excluded, duplicated or malformed key,
        a bad value, or a missing required option.
    """
    entry = registry[algorithm]
    options = {}
    for item in assignments:
        if "=" not in item:
            raise ValueError(f"{remedy} expects KEY=VALUE, got {item!r}")
        key, raw = item.split("=", 1)
        # Both sides are stripped: a shell-quoted `--set "linkage = ward"`
        # would otherwise carry the space into the value, where a string
        # option smuggles it through and a bool option refuses outright.
        key, raw = key.strip(), raw.strip()
        if key in _FLAG_OPTIONS:
            raise ValueError(
                f"set {key} with {_FLAG_OPTIONS[key]}, not {remedy}")
        if key in entry.sequence:
            raise ValueError(
                f"{key} takes a sequence; use the Python API for it")
        if key not in entry.known():
            hint = difflib.get_close_matches(key, entry.known(), n=1)
            suffix = f"; did you mean {hint[0]!r}?" if hint else ""
            raise ValueError(f"{algorithm} has no option {key!r}{suffix}")
        if key in options:
            raise ValueError(f"{key} given twice in {remedy}")
        if swept is not None and key == swept:
            raise ValueError(
                f"{key} is the swept parameter; drop it from {remedy}")
        options[key] = coerce(algorithm, key, raw, entry)
    for name in entry.required:
        if name not in options and name != swept:
            raise ValueError(
                f"{algorithm} requires {name}=…; supply it with {remedy}")
    return options
