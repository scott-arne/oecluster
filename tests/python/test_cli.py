"""The oecluster command line: registry, inputs, commands and output."""
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


def test_an_integer_option_refuses_a_non_integral_value():
    # The library coerces with a bare int(), so 1.5 would silently become 1.
    registry = _registry()
    with pytest.raises(ValueError, match="must be an integer"):
        _cli_registry.resolve("hdbscan", ["min_samples=1.5"], registry)


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


def test_a_non_finite_float_option_is_refused():
    # float("inf") parses happily, so without the guard an infinite
    # threshold reaches the native layer.
    with pytest.raises(ValueError, match="must be finite"):
        _cli_registry.resolve("butina", ["threshold=inf"], _registry())


def test_every_declared_override_names_a_real_parameter():
    # A roster rename must fail the suite rather than leave a stale schema.
    import inspect

    from oecluster._parameter_selection import _roster
    roster = _roster()
    declared = (_cli_registry._INTEGER | _cli_registry._CONDITIONALLY_REQUIRED
                | _cli_registry._SEQUENCE)
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
    assert "threshold" in result.output
    assert "required" in result.output


def test_algorithms_explains_an_ineligible_entry():
    result = _run("algorithms", "murcko")
    assert result.exit_code == 0
    assert "not a distance matrix" in result.output


def test_algorithms_refuses_an_unknown_name():
    result = _run("algorithms", "nosuch")
    assert result.exit_code == 2
