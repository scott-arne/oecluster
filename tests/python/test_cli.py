"""The oecluster command line: registry, inputs, commands and output."""
import time

import pytest
from click.testing import CliRunner
from oecluster import _cli_registry
from oecluster._cli_main import cli

# Later steps and tasks append test functions at the END of this file and
# grow the block above -- never beside the tests, which ruff rejects as E402
# and reorders under I001. Import only what the step actually uses: an
# import added ahead of its first use is rejected as F401. Each step below
# prints the block verbatim at that point. `oecluster` is third-party here
# (the package lives under `python/`), so it sorts with numpy and pytest.


def _registry():
    return _cli_registry.build()


def test_eligibility_is_derived_from_the_first_parameter():
    # Not a hardcoded list: `distance_matrix` and `items` entries can consume
    # a matrix, `fingerprints` and `mols` entries cannot, so a new roster
    # entry classifies itself.
    registry = _registry()
    eligible = sorted(n for n, e in registry.items() if e.eligible)
    assert eligible == ["agglomerative", "butina", "dbscan", "hdbscan",
                        "jarvis_patrick", "k_medoids", "leiden",
                        "sphere_exclusion"]
    refused = sorted(n for n, e in registry.items() if not e.eligible)
    assert refused == ["bitbirch", "bitbirch_recluster", "bitbirch_refine",
                       "murcko"]


def test_keyword_only_parameters_without_defaults_are_required():
    # jarvis_patrick.kmin is KEYWORD_ONLY with no default, so a requiredness
    # check that inspects only positional kinds misses it.
    registry = _registry()
    assert "kmin" in registry["jarvis_patrick"].required
    assert "threshold" in registry["butina"].required
    assert "eps" in registry["dbscan"].required


def test_conditionally_required_options_are_required_for_a_matrix():
    registry = _registry()
    assert "k" in registry["jarvis_patrick"].required
    assert "k" in registry["leiden"].required


def test_no_variadic_parameter_reaches_any_schema():
    # jarvis_patrick, leiden and sphere_exclusion all carry **kwargs, and
    # under "no default means required" a missed one becomes a bogus
    # required option.
    for entry in _registry().values():
        assert "kwargs" not in entry.known()
        assert "args" not in entry.known()


def test_flag_routed_options_are_not_settable():
    registry = _registry()
    for entry in registry.values():
        assert "num_threads" not in entry.known()
        assert "allow_nonmetric" not in entry.known()


def test_the_four_algorithms_that_accept_nonmetric_are_detected():
    registry = _registry()
    accepting = sorted(n for n, e in registry.items() if e.accepts_nonmetric)
    assert accepting == ["agglomerative", "butina", "dbscan", "hdbscan"]


def test_items_entries_hide_matrix_inapplicable_options():
    registry = _registry()
    assert "comparison" not in registry["leiden"].known()
    assert "similarity" not in registry["leiden"].known()


def test_resolve_coerces_by_the_defaults_type():
    registry = _registry()
    options = _cli_registry.resolve(
        "butina", ["threshold=0.3", "reordering=true"], registry)
    assert options == {"threshold": 0.3, "reordering": True}
    assert isinstance(options["reordering"], bool)
    linkage = _cli_registry.resolve(
        "agglomerative", ["linkage=ward", "n_clusters=4"], registry)
    assert linkage == {"linkage": "ward", "n_clusters": 4}


def test_whitespace_around_a_value_is_ignored():
    # A shell-quoted `--set "linkage = ward"` must not smuggle a leading
    # space into the value: the key is stripped, so the value has to be too,
    # or a bool assignment spelled that way is refused outright.
    registry = _registry()
    assert _cli_registry.resolve(
        "agglomerative", ["linkage = ward"], registry) == {"linkage": "ward"}
    assert _cli_registry.resolve(
        "butina", ["threshold=0.3", "reordering = true"],
        registry) == {"threshold": 0.3, "reordering": True}


@pytest.mark.parametrize("algorithm, assignment", [
    ("hdbscan", "min_samples=1.5"),
    # 1e-324 underflows to 0.0 in binary floating point, so a coercion that
    # round-trips through float reads it as the integer 0 and accepts it.
    ("jarvis_patrick", "kmin=1e-324"),
])
def test_an_integer_option_refuses_a_non_integral_value(algorithm, assignment):
    # The library coerces with a bare int(), so 1.5 would silently become 1.
    with pytest.raises(ValueError, match="must be an integer"):
        _cli_registry.resolve(algorithm, [assignment], _registry())


@pytest.mark.parametrize("assignment, expected", [
    ("n_iterations=12", 12),
    # Accepted deliberately: an integral float spelling is unambiguous.
    ("n_iterations=12.0", 12),
    # Past 2**53 a binary float can no longer represent consecutive
    # integers, so a float round-trip silently returns a different number
    # than the user typed.
    ("n_iterations=9007199254740993", 2**53 + 1),
    ("n_iterations=9223372036854775807", 2**63 - 1),
])
def test_an_integer_option_is_preserved_exactly(assignment, expected):
    options = _cli_registry.resolve(
        "leiden", [assignment, "k=5"], _registry())
    assert options["n_iterations"] == expected


@pytest.mark.parametrize("value", [
    # An exponent typo in --set used to be built out digit by digit before
    # anything looked at the magnitude: 1e100000 hit CPython's 4300-digit
    # conversion limit, whose message names sys.set_int_max_str_digits at a
    # user who only mistyped a number.
    "1e100000", "1e400",
    # 2**64, the first magnitude past the widest native parameter.
    "18446744073709551616", "-18446744073709551616",
])
def test_an_out_of_range_integer_option_is_refused(value):
    with pytest.raises(ValueError, match="out of range"):
        _cli_registry.resolve("leiden", [f"n_iterations={value}", "k=5"],
                              _registry())


def test_refusing_an_enormous_exponent_is_prompt():
    # The bound has to be applied before int() materializes the digits.
    # Building 10**999999999 consumes unbounded CPU and memory, so the
    # command hung before reading a single file and had to be killed.
    start = time.monotonic()
    with pytest.raises(ValueError, match="out of range"):
        _cli_registry.resolve("leiden", ["n_iterations=1e999999999", "k=5"],
                              _registry())
    assert time.monotonic() - start < 1.0


def test_the_largest_in_range_integer_is_accepted_exactly():
    # leiden's seed is a uint64 natively, so 2**64-1 is the widest value any
    # of these parameters takes and the bound must not refuse it.
    options = _cli_registry.resolve(
        "leiden", ["seed=18446744073709551615", "k=5"], _registry())
    assert options["seed"] == 2**64 - 1


def test_an_integer_option_coerces_the_same_way_whatever_its_default():
    # dbscan.min_samples has an int default and hdbscan.min_samples a None
    # default with a declared type. Both are integer options, so they must
    # not disagree about what an integer looks like or share none of the
    # range guard.
    registry = _registry()
    assert _cli_registry.resolve(
        "dbscan", ["eps=0.5", "min_samples=4.0"], registry)["min_samples"] == 4
    with pytest.raises(ValueError, match="out of range"):
        _cli_registry.resolve("dbscan", ["eps=0.5", "min_samples=1e400"],
                              registry)


@pytest.mark.parametrize("assignment", [
    "distance_threshold=true", "distance_threshold=false",
    "distance_threshold=ward",
])
def test_a_numeric_option_refuses_a_non_numeric_value(assignment):
    # distance_threshold defaults to None, so there is no default to infer a
    # type from. Guessing from the value shape accepted a bool, which
    # agglomerative then reads as 1.0, and a string, which only blows up
    # after the matrix has been loaded -- defeating the point of validating
    # before any file is read.
    with pytest.raises(ValueError, match="must be a number"):
        _cli_registry.resolve("agglomerative", [assignment], _registry())


def test_every_option_defaulting_to_none_declares_a_type():
    # The companion to the override tripwire: an option whose default is
    # None has no type to infer, so a new roster entry with one must declare
    # it here rather than reach a guess.
    declared = (_cli_registry._INTEGER | _cli_registry._FLOAT
                | _cli_registry._REQUIRED_TYPES.keys())
    for name, entry in _registry().items():
        for option, default in entry.optional.items():
            if default is None:
                assert f"{name}.{option}" in declared, f"{name}.{option}"


def test_the_chunk_size_tuning_knob_stays_hidden():
    # Emptying _HIDDEN would leak chunk_size into every schema with the rest
    # of the suite still green. Both halves are pinned, so the assertion
    # cannot go vacuous if the roster drops the parameter instead.
    import inspect

    from oecluster._parameter_selection import _roster
    carriers = [name for name, fn in _roster().items()
                if "chunk_size" in inspect.signature(fn).parameters]
    assert carriers, "no roster entry takes chunk_size; this test is vacuous"
    registry = _registry()
    for entry in registry.values():
        assert "chunk_size" not in entry.known()
    with pytest.raises(ValueError, match="no option"):
        _cli_registry.resolve("butina", ["threshold=0.3", "chunk_size=64"],
                              registry)


@pytest.mark.parametrize("assignments, needle", [
    (["threshold=0.3", "thresold=0.4"], "did you mean"),
    (["threshold=0.3", "num_threads=4"], "--threads"),
    (["reordering=true"], "requires"),
    (["threshold=0.3", "threshold=0.4"], "given twice"),
    (["bogus"], "KEY=VALUE"),
])
def test_resolve_refuses_bad_assignments(assignments, needle):
    with pytest.raises(ValueError, match=needle):
        _cli_registry.resolve("butina", assignments, _registry())


def test_a_sequence_option_points_at_the_python_api():
    with pytest.raises(ValueError, match="sequence"):
        _cli_registry.resolve("k_medoids", ["initial_medoids=1,2"], _registry())


def test_the_swept_parameter_satisfies_requiredness():
    # select-parameter supplies threshold per sweep value, so demanding it in
    # --set would reject the command's normal form.
    registry = _registry()
    assert _cli_registry.resolve("butina", [], registry, swept="threshold") == {}
    with pytest.raises(ValueError, match="swept parameter"):
        _cli_registry.resolve("butina", ["threshold=0.3"], registry,
                              swept="threshold")


@pytest.mark.parametrize("assignment", ["n_clusters=inf", "n_clusters=nan"])
def test_a_non_finite_integer_option_is_refused(assignment):
    # A bare int(float("inf")) raises OverflowError and int(float("nan"))
    # raises ValueError with an opaque message; neither is a usage error the
    # command layer can translate.
    with pytest.raises(ValueError, match="must be an integer"):
        _cli_registry.resolve("k_medoids", [assignment], _registry())


@pytest.mark.parametrize("algorithm, assignment", [
    # One case per coercion branch, because the guard has to sit on all
    # three: a declared required float, an option whose default is None and
    # so takes the permissive fallback, and one whose default is a float.
    ("butina", "threshold=inf"),
    ("agglomerative", "distance_threshold=inf"),
    ("hdbscan", "cluster_selection_epsilon=inf"),
    ("hdbscan", "alpha=nan"),
])
def test_a_non_finite_float_option_is_refused(algorithm, assignment):
    # float("inf") parses happily, so without the guard an infinite
    # threshold reaches the native layer.
    with pytest.raises(ValueError, match="must be finite"):
        _cli_registry.resolve(algorithm, [assignment], _registry())


def test_every_declared_override_names_a_real_parameter():
    # A roster rename must fail the suite rather than leave a stale schema.
    import inspect

    from oecluster._parameter_selection import _roster
    roster = _roster()
    declared = (_cli_registry._INTEGER | _cli_registry._CONDITIONALLY_REQUIRED
                | _cli_registry._SEQUENCE | _cli_registry._REQUIRED_TYPES.keys())
    for key in declared:
        algorithm, option = key.split(".", 1)
        assert algorithm in roster, key
        assert option in inspect.signature(roster[algorithm]).parameters, key


def _run(*args):
    return CliRunner().invoke(cli, list(args))


def test_algorithms_lists_every_roster_entry():
    result = _run("algorithms")
    assert result.exit_code == 0
    assert "butina" in result.output
    assert "bitbirch" in result.output


def test_algorithms_shows_one_algorithms_schema():
    result = _run("algorithms", "butina")
    assert result.exit_code == 0
    # Assert on the row, not on the words anywhere in the output: "required"
    # is also a column header, printed even for an algorithm whose required
    # list is empty, so the required/optional split could vanish entirely
    # and a bare substring check would still hold.
    rows = [line for line in result.output.splitlines() if "threshold" in line]
    assert len(rows) == 1
    assert "float" in rows[0]
    assert "required" in rows[0]


def test_algorithms_explains_an_ineligible_entry():
    result = _run("algorithms", "murcko")
    assert result.exit_code == 0
    assert "not a distance matrix" in result.output


def test_algorithms_refuses_an_unknown_name():
    result = _run("algorithms", "nosuch")
    assert result.exit_code == 2
