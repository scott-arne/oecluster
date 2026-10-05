"""The oecluster command line: registry, inputs, commands and output."""
import json
import os
import time
import traceback
import warnings
from pathlib import Path

import click
import numpy as np
import oecluster
import pytest
from click.testing import CliRunner
from oecluster import _cli_input, _cli_main, _cli_registry, _cli_render
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


# Fixtures below build oepdist-shaped files: a condensed .npy or raw .bin
# plus the JSON sidecar oepdist writes beside it.
def _condensed(n=12):
    rng = np.random.default_rng(0)
    points = np.vstack([rng.normal(centre, 0.2, (n // 2, 2))
                        for centre in (0.0, 4.0)])
    return np.array([float(np.linalg.norm(points[i] - points[j]))
                     for i in range(n) for j in range(i + 1, n)])


def _write_sidecar(path, n, **overrides):
    body = {"mode": "pdist", "comparison": "tanimoto",
            "params": {"similarity": False}, "n_rows": n, "n_cols": n,
            "row_labels": [f"m{i}" for i in range(n)], "col_labels": []}
    body.update(overrides)
    with open(os.path.splitext(path)[0] + ".json", "w", encoding="utf-8") as f:
        json.dump(body, f)


def _npy(tmp_path, **overrides):
    values = _condensed()
    path = str(tmp_path / "d.npy")
    np.save(path, values)
    _write_sidecar(path, 12, **overrides)
    return path


def test_an_npy_and_its_sidecar_load_with_provenance(tmp_path):
    matrix = _cli_input.load(_npy(tmp_path), warn=lambda message: None)
    assert matrix.num_samples == 12
    assert list(matrix.labels)[:2] == ["m0", "m1"]
    assert matrix.comparison_name == "tanimoto"


def test_bin_and_npy_agree(tmp_path):
    values = _condensed()
    npy = str(tmp_path / "a.npy")
    np.save(npy, values)
    _write_sidecar(npy, 12)
    raw = str(tmp_path / "b.bin")
    values.tofile(raw)
    _write_sidecar(raw, 12)
    first = _cli_input.load(npy, warn=lambda message: None)
    second = _cli_input.load(raw, warn=lambda message: None)
    assert np.allclose(list(first.condensed), list(second.condensed))


def test_an_npz_round_trips(tmp_path):
    source = _cli_input.load(_npy(tmp_path), warn=lambda message: None)
    out = str(tmp_path / "m.npz")
    source.to_file(out)
    seen = []
    loaded = _cli_input.load(out, warn=seen.append)
    assert loaded.num_samples == source.num_samples
    # from_condensed stamps is_distance "unknown" unconditionally, so every
    # .npz written from a raw input warns on reload. Asserted here because
    # this is the only test that reaches the .npz warning at all: with the
    # callback discarded, deleting the branch left the suite green, and the
    # orientation warning is what stands between a similarity file and
    # silently inverted clusters.
    assert "orientation unproven" in seen[0]


def test_a_similarity_file_is_refused(tmp_path):
    # The most dangerous input there is: clustering it would invert every
    # result, and nothing in the numbers reveals it.
    path = _npy(tmp_path, params={"similarity": True})
    with pytest.raises(_cli_input.InputError, match="similarities"):
        _cli_input.load(path, warn=lambda message: None)


def test_an_absent_similarity_flag_warns_but_proceeds(tmp_path):
    # The rocs case. The library refuses only a proven similarity, so the CLI
    # matches it rather than being stricter, and warns instead.
    path = _npy(tmp_path, params={})
    seen = []
    matrix = _cli_input.load(path, warn=seen.append)
    assert matrix.num_samples == 12
    assert "orientation unproven" in seen[0]


@pytest.mark.parametrize("overrides, needle", [
    ({"mode": "cdist"}, "cross-distance"),
    ({"n_rows": 99, "n_cols": 99}, "sidecar says"),
    ({"n_rows": 12, "n_cols": 7}, "not square"),
])
def test_a_bad_sidecar_is_refused(tmp_path, overrides, needle):
    path = _npy(tmp_path, **overrides)
    with pytest.raises(_cli_input.InputError, match=needle):
        _cli_input.load(path, warn=lambda message: None)


def test_a_missing_sidecar_is_refused(tmp_path):
    path = _npy(tmp_path)
    os.remove(os.path.splitext(path)[0] + ".json")
    with pytest.raises(_cli_input.InputError, match="sidecar"):
        _cli_input.load(path, warn=lambda message: None)


def test_an_unparseable_sidecar_names_the_likely_cause(tmp_path):
    # oepdist's JsonEscape handles only quotes and backslashes, so a title
    # with a newline produces invalid JSON.
    path = _npy(tmp_path)
    with open(os.path.splitext(path)[0] + ".json", "w", encoding="utf-8") as f:
        f.write('{"mode": "pdist", "row_labels": ["a\nb"]}')
    with pytest.raises(_cli_input.InputError, match="newline"):
        _cli_input.load(path, warn=lambda message: None)


def test_a_non_triangular_count_is_refused(tmp_path):
    path = str(tmp_path / "odd.npy")
    np.save(path, np.zeros(5))
    _write_sidecar(path, 12)
    with pytest.raises(_cli_input.InputError, match="condensed"):
        _cli_input.load(path, warn=lambda message: None)


def test_csv_is_refused_with_a_remedy(tmp_path):
    path = tmp_path / "m.csv"
    path.write_text(",a,b\na,0,1\nb,1,0\n", encoding="utf-8")
    with pytest.raises(_cli_input.InputError, match="re-run oepdist"):
        _cli_input.load(str(path), warn=lambda message: None)


def _patched_npz(tmp_path, name, *, drop=(), **patch):
    """A .npz written by the API, rewritten with one field corrupted."""
    matrix = oecluster.SymmetricDistanceMatrix.from_condensed(
        np.array([0.1, 0.9, 0.2, 0.8, 0.3, 0.7]), labels=["a", "b", "c", "d"])
    good = str(tmp_path / "good.npz")
    matrix.to_file(good)
    with np.load(good, allow_pickle=False) as archive:
        fields = {key: archive[key] for key in archive.files}
    for key in drop:
        del fields[key]
    fields.update(patch)
    path = str(tmp_path / f"{name}.npz")
    np.savez(path, **fields)
    return path


def test_malformed_npz_metadata_is_an_input_error(tmp_path):
    # facts_json is decoded JSON of any shape. A scalar reaches dict.update
    # deep in the library and raises TypeError -- not one of the types the
    # loader was catching, so it escaped the InputError contract and
    # _translate would have reported an unreadable file as a usage error.
    path = _patched_npz(tmp_path, "facts", facts_json=np.array("42"))
    with pytest.raises(_cli_input.InputError, match="cannot read"):
        _cli_input.load(path, warn=lambda message: None)


@pytest.mark.parametrize("labels, count", [
    (["a"], 1),
    (["a", "b", "c", "d", "e", "f"], 6),
])
def test_an_npz_label_count_mismatch_is_refused(tmp_path, labels, count):
    # from_condensed enforces this, but the .npz path restores state
    # directly and does not. A short list silently truncated every per-item
    # output: a four-item matrix with one label wrote one CSV row.
    path = _patched_npz(tmp_path, "labels", labels=np.array(labels))
    with pytest.raises(_cli_input.InputError, match=f"{count} labels"):
        _cli_input.load(path, warn=lambda message: None)


def test_a_non_object_params_sidecar_is_refused(tmp_path):
    # `params` is raw JSON from oepdist, so it can be any type. A truthy
    # non-object reaches params.get(...) and raises AttributeError, which
    # the command layer never promised and would not translate.
    path = _npy(tmp_path, params="oops")
    with pytest.raises(_cli_input.InputError, match="non-object"):
        _cli_input.load(path, warn=lambda message: None)


def test_a_cross_distance_npz_is_refused(tmp_path):
    # Build one through the public API so the test pins the real shape.
    # CrossDistanceMatrix has no from_values: the real API is
    # __init__(matrix, comparison_name, ...) and from_file.
    cross = oecluster.CrossDistanceMatrix(np.zeros((2, 3)), "precomputed")
    out = str(tmp_path / "cross.npz")
    cross.to_file(out)
    with pytest.raises(_cli_input.InputError, match="symmetric"):
        _cli_input.load(out, warn=lambda message: None)


@pytest.mark.parametrize("dtype", [
    "complex128", "datetime64[s]", "bool", "<U4",
    [("a", "<f8"), ("b", "<i4")],
])
def test_a_non_real_numeric_npy_is_refused(tmp_path, dtype):
    # Every one of these used to reach float64 and be clustered: complex
    # kept its real part behind a ComplexWarning (defeating the refusal
    # from_condensed places at __init__.py:1350 for exactly this reason),
    # datetime64 became epoch seconds, and bool became 0/1 -- which for an
    # adjacency matrix is a similarity, the inversion this module exists to
    # refuse. The structured dtype raised a bare TypeError from the cast.
    path = str(tmp_path / "odd.npy")
    np.save(path, np.zeros(6, dtype=dtype))
    _write_sidecar(path, 4)
    with pytest.raises(_cli_input.InputError, match="real numbers"):
        _cli_input.load(path, warn=lambda message: None)


def test_a_complex_npy_is_refused_before_any_warning(tmp_path):
    # The cast emitted ComplexWarning, so under -W error the warning itself
    # escaped load() ahead of any refusal: not an InputError, and not even
    # an exception the command layer could name a file for.
    path = str(tmp_path / "c.npy")
    np.save(path, np.array([1 + 1j, 2 + 2j, 3 + 3j]))
    _write_sidecar(path, 3)
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        with pytest.raises(_cli_input.InputError, match="real numbers"):
            _cli_input.load(path, warn=lambda message: None)


@pytest.mark.parametrize("name, drop, patch", [
    # The first two escaped the enumerated catch: a 0-d condensed indexes
    # shape[0] inside from_file, and a negative num_samples reaches
    # DenseStorage's size_t. The rest already refused and are pinned so a
    # future narrowing of the net cannot quietly let them through.
    ("zero_d", (), {"condensed": np.array(1.0)}),
    ("negative", (), {"num_samples": np.array(-4), "condensed": np.zeros(10)}),
    ("rank", (), {"condensed": np.zeros((2, 3))}),
    ("text", (), {"condensed": np.array(["a"] * 6)}),
    ("complex", (), {"condensed": np.zeros(6, dtype="complex128")}),
    ("flat_labels", (), {"labels": np.zeros(())}),
    ("no_values", ("condensed",), {}),
])
def test_a_structurally_broken_npz_is_an_input_error(tmp_path, name, drop,
                                                     patch):
    path = _patched_npz(tmp_path, name, drop=drop, **patch)
    with pytest.raises(_cli_input.InputError, match="cannot read"):
        _cli_input.load(path, warn=lambda message: None)


def test_an_npz_archive_under_a_npy_name_is_refused(tmp_path):
    # np.load reads the bytes, not the name, and hands back a lazy NpzFile
    # for any zip: the dtype check then read .dtype on an archive object and
    # raised AttributeError. Found by fuzzing the fix, not by review.
    source = _cli_input.load(_npy(tmp_path), warn=lambda message: None)
    real = str(tmp_path / "real.npz")
    source.to_file(real)
    renamed = tmp_path / "renamed.npy"
    renamed.write_bytes(Path(real).read_bytes())
    _write_sidecar(str(renamed), 12)
    with pytest.raises(_cli_input.InputError, match="archive"):
        _cli_input.load(str(renamed), warn=lambda message: None)


@pytest.mark.parametrize("flag", [1, 0, "true", "false", [], {}])
def test_a_sidecar_similarity_that_is_not_a_boolean_is_refused(tmp_path, flag):
    # The orientation gate tests `is True` and `is None`, so any other value
    # fell between the two branches: 1 and "true" were clustered as proven
    # distances on the strength of a flag that says the opposite, and
    # without even the unproven-orientation warning. Found by fuzzing.
    path = _npy(tmp_path, params={"similarity": flag})
    with pytest.raises(_cli_input.InputError, match="neither true nor false"):
        _cli_input.load(path, warn=lambda message: None)


def test_a_memory_error_is_not_reported_as_a_bad_file(tmp_path, monkeypatch):
    # The inverted catch has to let this one through: the command layer maps
    # MemoryError to exit 1 with the --mmap hint, and a file too large to
    # hold is not a malformed one.
    def explode(*args, **kwargs):
        raise MemoryError("Unable to allocate 32.0 GiB")

    monkeypatch.setattr(_cli_input.np, "load", explode)
    with pytest.raises(MemoryError):
        _cli_input.load(_npy(tmp_path), warn=lambda message: None)


def test_a_library_bug_is_blamed_on_the_file_but_stays_diagnosable(
        tmp_path, monkeypatch):
    # The price of the inverted net, pinned rather than wished away:
    # load_distance_matrix is this project's own parser, so a bug inside it
    # is caught by the .npz guard and reported as an unusable file. What the
    # guard must not do is erase it -- the exception's own text is carried
    # into the message verbatim, which is what keeps the misattributed
    # report diagnosable.
    source = _cli_input.load(_npy(tmp_path), warn=lambda message: None)
    out = str(tmp_path / "m.npz")
    source.to_file(out)

    def explode(path):
        raise NameError("name 'maht' is not defined")

    monkeypatch.setattr(_cli_input, "load_distance_matrix", explode)
    with pytest.raises(_cli_input.InputError, match="'maht' is not defined"):
        _cli_input.load(out, warn=lambda message: None)


def test_the_guards_do_not_wrap_whole_functions(tmp_path, monkeypatch):
    # _items_from_pairs is called outside every guard, so a fault in it
    # escapes as itself. Read this narrowly: it pins only that the guards
    # stay around individual library calls, and would fail if one were ever
    # widened to the body of load(). It proves nothing about code reached
    # inside a guard -- the test above states what happens there.
    def explode(count):
        raise NameError("name 'maht' is not defined")

    monkeypatch.setattr(_cli_input, "_items_from_pairs", explode)
    with pytest.raises(NameError):
        _cli_input.load(_npy(tmp_path), warn=lambda message: None)


@pytest.mark.parametrize("wrap", [str, Path, os.fsencode])
def test_any_path_like_argument_is_accepted(tmp_path, wrap):
    # Tasks 3-5 declare the argument with click.Path(path_type=Path), and
    # every format check here is a string operation. Bytes cannot arrive
    # from click, but os.fspath passes them through unchanged, so decoding
    # is what keeps the str-only checks below honest.
    matrix = _cli_input.load(wrap(_npy(tmp_path)), warn=lambda message: None)
    assert matrix.num_samples == 12


def test_an_npz_whose_orientation_fact_is_not_a_boolean_warns(tmp_path):
    # The sidecar's similarity flag and the archive's is_distance fact are
    # the same gate read from two files, and both are arbitrary JSON. The
    # library refuses only `is_distance is False`, so a corrupted fact of 0
    # -- falsy, and plainly not a proven distance -- passed the gate, and an
    # equality test against "unknown" meant the CLI did not warn either.
    matrix = oecluster.SymmetricDistanceMatrix.from_condensed(
        np.array([0.1, 0.9, 0.2, 0.8, 0.3, 0.7]), labels=["a", "b", "c", "d"])
    good = str(tmp_path / "good.npz")
    matrix.to_file(good)
    with np.load(good, allow_pickle=False) as archive:
        fields = {key: archive[key] for key in archive.files}
    fields["facts_json"] = np.array(json.dumps({"is_distance": 0}))
    path = str(tmp_path / "fact.npz")
    np.savez(path, **fields)
    seen = []
    assert _cli_input.load(path, warn=seen.append).num_samples == 4
    assert "orientation unproven" in seen[0]


@pytest.mark.parametrize("labels", [
    "abcdefghijkl",
    {f"m{i}": i for i in range(12)},
    5,
])
def test_sidecar_row_labels_that_are_not_an_array_are_refused(tmp_path,
                                                              labels):
    # from_condensed runs list() over whatever it is given, so a string
    # spells itself out and a dict hands over its keys. Both of these have
    # twelve entries, so every length check passed and the per-item output
    # was labelled with identities the user never wrote.
    path = _npy(tmp_path, row_labels=labels)
    with pytest.raises(_cli_input.InputError, match="row_labels"):
        _cli_input.load(path, warn=lambda message: None)


@pytest.mark.parametrize("entry", [0, None, ["m3"], 1.5, True])
def test_a_sidecar_row_label_that_is_not_a_string_is_refused(tmp_path, entry):
    # Refused rather than coerced: str(None) is "None" and str(0) is "0",
    # which are invented identities exactly like the cases above. oepdist
    # writes titles as JSON strings, so a non-string is a broken producer.
    labels = [f"m{i}" for i in range(12)]
    labels[3] = entry
    path = _npy(tmp_path, row_labels=labels)
    with pytest.raises(_cli_input.InputError, match="row_labels"):
        _cli_input.load(path, warn=lambda message: None)


def test_an_npz_holding_a_negative_distance_is_refused(tmp_path):
    # Reachable through the public API, not only by editing an archive:
    # condensed is writeable and to_file keeps the facts it was built with,
    # so the file claims a proven distance matrix while holding a negative
    # value. The clustering gate does not cover it either -- it re-scans for
    # non-finite values, not for negatives -- so butina clustered this.
    matrix = oecluster.SymmetricDistanceMatrix.from_condensed(
        np.array([0.1, 0.9, 0.2, 0.8, 0.3, 0.7]), labels=["a", "b", "c", "d"])
    matrix.condensed[0] = -1.0
    out = str(tmp_path / "neg.npz")
    matrix.to_file(out)
    with pytest.raises(_cli_input.InputError, match="negative"):
        _cli_input.load(out, warn=lambda message: None)


def test_a_raw_file_holding_a_negative_distance_is_refused(tmp_path):
    # The twin of the case above, pinned so the two paths cannot drift: it
    # is the asymmetry between them that let the .npz through.
    values = _condensed()
    values[0] = -1.0
    path = str(tmp_path / "neg.npy")
    np.save(path, values)
    _write_sidecar(path, 12)
    with pytest.raises(_cli_input.InputError, match="negative"):
        _cli_input.load(path, warn=lambda message: None)


@pytest.mark.parametrize("params", [[], 0, "", False])
def test_a_falsy_non_object_params_sidecar_is_refused(tmp_path, params):
    # `sidecar.get("params") or {}` applied the default before the type
    # check, so every falsy non-object was read as "no params" instead of
    # as the producer fault it is -- and then drew the unproven-orientation
    # warning rather than a refusal.
    path = _npy(tmp_path, params=params)
    with pytest.raises(_cli_input.InputError, match="non-object"):
        _cli_input.load(path, warn=lambda message: None)


@pytest.mark.parametrize("suffix", [".NPY", ".BIN"])
def test_an_uppercase_extension_is_accepted(tmp_path, suffix):
    # oepdist lowercases the extension before dispatching on it
    # (tools/OutputWriter.cpp:19-24), so `-o out.NPY` writes a real .NPY
    # that this loader was refusing as an unsupported format.
    values = _condensed()
    if suffix == ".NPY":
        np.save(str(tmp_path / "source.npy"), values)
        data = (tmp_path / "source.npy").read_bytes()
    else:
        data = values.tobytes()
    path = tmp_path / f"d{suffix}"
    path.write_bytes(data)
    _write_sidecar(str(path), 12)
    matrix = _cli_input.load(str(path), warn=lambda message: None)
    assert matrix.num_samples == 12


def _matrix(tmp_path, n=40):
    rng = np.random.default_rng(0)
    points = np.vstack([rng.normal(centre, 0.25, (n // 4, 2))
                        for centre in (0.0, 4.0, 8.0, 12.0)])
    values = np.array([float(np.linalg.norm(points[i] - points[j]))
                       for i in range(n) for j in range(i + 1, n)])
    matrix = oecluster.SymmetricDistanceMatrix.from_condensed(
        values, labels=[f"m{i}" for i in range(n)])
    path = str(tmp_path / "m.npz")
    matrix.to_file(path)
    return path


def test_cluster_finds_the_four_blobs(tmp_path):
    result = _run("cluster", _matrix(tmp_path), "--algorithm", "butina",
                  "--set", "threshold=1.0")
    assert result.exit_code == 0
    assert "clusters=4" in result.stdout


def test_cluster_writes_json_matching_the_schema(tmp_path):
    out = str(tmp_path / "r.json")
    path = _matrix(tmp_path)
    result = _run("cluster", path, "--algorithm", "butina",
                  "--set", "threshold=1.0", "--output", out)
    assert result.exit_code == 0
    document = json.loads(Path(out).read_text(encoding="utf-8"))
    # The whole envelope, not a key count: a renamed or dropped field is
    # what a downstream reader breaks on, and a count survives a rename.
    assert sorted(document) == ["command", "input", "result",
                                "schema_version", "spec"]
    assert document["schema_version"] == 1
    assert document["command"] == "cluster"
    assert document["input"] == {"path": path, "num_items": 40,
                                 "orientation": "unknown"}
    assert document["spec"] == {
        "algorithm": "butina",
        "options": {"threshold": 1.0, "num_threads": 0}}
    assert sorted(document["result"]) == ["ids", "labels", "num_clusters",
                                          "num_noise"]
    assert document["result"]["ids"] == [f"m{i}" for i in range(40)]
    assert document["result"]["num_clusters"] == 4
    assert document["result"]["num_noise"] == 0
    assert len(document["result"]["labels"]) == 40


def test_cluster_writes_csv_with_the_documented_header(tmp_path):
    out = str(tmp_path / "r.csv")
    _run("cluster", _matrix(tmp_path), "--algorithm", "butina",
         "--set", "threshold=1.0", "--output", out)
    lines = Path(out).read_text(encoding="utf-8").splitlines()
    assert lines[0] == "id,label"
    # Every pair, not just the first: a truncated or misaligned id column
    # reads perfectly well from row 1 alone.
    assert len(lines) == 41
    rows = [line.split(",") for line in lines[1:]]
    assert [row[0] for row in rows] == [f"m{i}" for i in range(40)]
    assigned = [row[1] for row in rows]
    assert sorted(assigned.count(value)
                  for value in set(assigned)) == [10, 10, 10, 10]


def test_quiet_silences_the_table_but_not_the_warning(tmp_path):
    # click 8.5 has no mix_stderr; result.stdout and result.stderr are
    # separate, and the orientation warning must survive --quiet.
    out = str(tmp_path / "r.csv")
    result = _run("cluster", _matrix(tmp_path), "--algorithm", "butina",
                  "--set", "threshold=1.0", "--quiet", "--output", out)
    assert result.exit_code == 0
    assert "clusters=" not in result.stdout
    assert "orientation unproven" in result.stderr


def test_an_unsupported_output_extension_is_refused(tmp_path):
    result = _run("cluster", _matrix(tmp_path), "--algorithm", "butina",
                  "--set", "threshold=1.0", "--output",
                  str(tmp_path / "r.txt"))
    assert result.exit_code == 2
    assert ".csv or .json" in result.output


@pytest.mark.parametrize("extra, code, needle", [
    ([], 2, "requires"),
    (["--set", "threshold=x"], 2, "number"),
    (["--set", "threshold=1.0", "--threads", "-1"], 2, "non-negative"),
])
def test_cluster_validation_exits_two(tmp_path, extra, code, needle):
    result = _run("cluster", _matrix(tmp_path), "--algorithm", "butina",
                  *extra)
    assert result.exit_code == code
    assert needle in result.output


def test_an_ineligible_algorithm_exits_two(tmp_path):
    result = _run("cluster", _matrix(tmp_path), "--algorithm", "murcko")
    assert result.exit_code == 2
    assert "not a distance matrix" in result.output


def test_a_missing_input_exits_one(tmp_path):
    result = _run("cluster", str(tmp_path / "gone.npz"),
                  "--algorithm", "butina", "--set", "threshold=1.0")
    assert result.exit_code == 1


def test_nan_and_infinity_encode_distinctly():
    assert _cli_render.encode(float("nan")) is None
    assert _cli_render.encode(float("inf")) == "inf"
    assert _cli_render.encode(float("-inf")) == "-inf"
    assert _cli_render.encode(0.5) == 0.5


def test_an_output_over_the_input_matrix_is_refused(tmp_path):
    # Checked before anything is loaded: the run would otherwise truncate
    # the file it is about to read.
    path = _matrix(tmp_path)
    result = _run("cluster", path, "--algorithm", "butina",
                  "--set", "threshold=1.0", "--output", path)
    assert result.exit_code == 2
    assert "matrix" in result.output


def test_an_output_over_the_input_sidecar_is_refused(tmp_path):
    # With input d.npy the obvious --output d.json names exactly the sidecar
    # that carries the shape and provenance the raw file does not.
    path = _npy(tmp_path)
    result = _run("cluster", path, "--algorithm", "butina", "--set",
                  "threshold=1.0", "--output",
                  os.path.splitext(path)[0] + ".json")
    assert result.exit_code == 2
    assert "sidecar" in result.output


def test_traceback_keeps_the_underlying_failure_in_the_chain(tmp_path):
    # _cli_input raises `from None`, which only hides a context it has
    # already recorded. --traceback promises the full traceback, and the
    # frames worth having are the failing library call's, not the frame that
    # renamed the failure.
    path = _patched_npz(tmp_path, "facts", facts_json=np.array("42"))
    result = _run("cluster", path, "--algorithm", "butina",
                  "--set", "threshold=1.0", "--traceback")
    assert result.exit_code == 1
    assert isinstance(result.exception, _cli_input.InputError)
    rendered = "".join(traceback.format_exception(result.exception))
    assert "During handling of the above exception" in rendered


def test_without_traceback_the_same_failure_is_a_clean_exit(tmp_path):
    path = _patched_npz(tmp_path, "facts", facts_json=np.array("42"))
    result = _run("cluster", path, "--algorithm", "butina",
                  "--set", "threshold=1.0")
    assert result.exit_code == 1
    assert result.exception.__class__ is SystemExit
    assert "cannot read" in result.output


def test_a_bytes_label_is_decoded_rather_than_repred(tmp_path):
    # A .npz carries whatever labels the Python API was handed, and numpy
    # brings a bytes label back as np.bytes_, whose str() is "b'a'" -- an
    # identity the user never wrote. The sidecar path refuses non-strings;
    # this path accepts them by design, so the output layer has to render
    # them rather than assume str.
    path = _patched_npz(tmp_path, "bytes",
                        labels=np.array([b"a", b"b", b"c", b"d"]))
    csv_path = str(tmp_path / "r.csv")
    json_path = str(tmp_path / "r.json")
    for out in (csv_path, json_path):
        result = _run("cluster", path, "--algorithm", "butina",
                      "--set", "threshold=1.0", "--output", out)
        assert result.exit_code == 0, result.output
    rows = Path(csv_path).read_text(encoding="utf-8").splitlines()
    assert rows[1].startswith("a,")
    document = json.loads(Path(json_path).read_text(encoding="utf-8"))
    assert document["result"]["ids"] == ["a", "b", "c", "d"]


def test_an_unlabelled_matrix_falls_back_to_indices(tmp_path):
    # _labels_of has to tell "no labels" from a real empty list: an
    # unlabelled matrix reports [], and the id column becomes positions.
    matrix = oecluster.SymmetricDistanceMatrix.from_condensed(
        np.array([0.1, 0.9, 0.2, 0.8, 0.3, 0.7]))
    path = str(tmp_path / "plain.npz")
    matrix.to_file(path)
    csv_path = str(tmp_path / "r.csv")
    json_path = str(tmp_path / "r.json")
    for out in (csv_path, json_path):
        assert _run("cluster", path, "--algorithm", "butina", "--set",
                    "threshold=1.0", "--output", out).exit_code == 0
    rows = Path(csv_path).read_text(encoding="utf-8").splitlines()
    assert rows[1].startswith("0,")
    document = json.loads(Path(json_path).read_text(encoding="utf-8"))
    assert "ids" not in document["result"]


def _case_folding(tmp_path):
    """:returns: True if this filesystem treats two casings as one file."""
    probe = tmp_path / "CaseProbe"
    probe.write_text("x", encoding="utf-8")
    return (tmp_path / "caseprobe").exists()


def _uppercase_npy(tmp_path):
    """An oepdist-shaped input named U.NPY, with its U.json sidecar."""
    np.save(str(tmp_path / "source.npy"), _condensed())
    path = tmp_path / "U.NPY"
    path.write_bytes((tmp_path / "source.npy").read_bytes())
    _write_sidecar(str(path), 12)
    return str(path)


def test_an_uppercase_input_extension_still_protects_the_sidecar(tmp_path):
    # _cli_input dispatches on the lowercased suffix because `oepdist -o
    # out.NPY` writes a real .NPY. Compared case-sensitively here, U.NPY
    # looked like a format that has no sidecar, so the run overwrote the
    # one file that makes the input readable at all.
    path = _uppercase_npy(tmp_path)
    sidecar = str(tmp_path / "U.json")
    result = _run("cluster", path, "--algorithm", "butina",
                  "--set", "threshold=1.0", "--output", sidecar)
    assert result.exit_code == 2
    assert "sidecar" in result.output
    assert json.loads(Path(sidecar).read_text(encoding="utf-8"))["mode"] \
        == "pdist"


def test_a_recased_destination_cannot_clobber_the_sidecar(tmp_path):
    # realpath compares bytes, so it reports d.json and D.JSON as two files
    # even where the filesystem has only one. The sidecar has to survive
    # either way: refused where the names collide, written beside it where
    # they do not.
    path = _npy(tmp_path)
    sidecar = os.path.splitext(path)[0] + ".json"
    recased = os.path.join(os.path.dirname(sidecar),
                           os.path.basename(sidecar).upper())
    result = _run("cluster", path, "--algorithm", "butina",
                  "--set", "threshold=1.0", "--output", recased)
    assert json.loads(Path(sidecar).read_text(encoding="utf-8"))["mode"] \
        == "pdist"
    if _case_folding(tmp_path):
        assert result.exit_code == 2
        assert "sidecar" in result.output
    else:
        assert result.exit_code == 0


@pytest.mark.parametrize("error", [click.exceptions.Exit(0), click.Abort()])
def test_translate_passes_clicks_own_control_flow_through(error):
    # Exit and Abort both subclass RuntimeError, so the arm that maps a
    # RuntimeError to exit 1 turned ctx.exit(0) into exit 1 with the
    # message "0", and an aborted confirmation into a failed run.
    @_cli_main._translate
    def body():
        raise error

    with pytest.raises(type(error)):
        body()


def test_write_output_emits_strict_json_for_non_finite_values(tmp_path):
    # Bare NaN and Infinity are a Python extension that strict parsers
    # reject, and consensus and stability produce NaN routinely.
    out = str(tmp_path / "r.json")
    _cli_render.write_output(out, {"document": {
        "x": float("nan"), "y": float("inf"),
        "z": [float("-inf"), 0.5], "nested": {"w": float("nan")}}})
    body = Path(out).read_text(encoding="utf-8")
    assert "NaN" not in body
    assert "Infinity" not in body
    assert json.loads(body) == {"x": None, "y": "inf",
                                "z": ["-inf", 0.5], "nested": {"w": None}}


def test_a_bad_output_extension_is_refused_before_any_work(tmp_path):
    # The writer validated the extension at write time, so the whole run
    # happened first: for cluster that wastes seconds, for the commands
    # that reuse the writer it wastes minutes.
    result = _run("cluster", _matrix(tmp_path), "--algorithm", "butina",
                  "--set", "threshold=1.0", "--output",
                  str(tmp_path / "r.txt"))
    assert result.exit_code == 2
    assert "clusters=" not in result.output
    assert "cluster" not in result.stdout


@pytest.mark.parametrize("value, chained", [("", False), ("0", False),
                                            ("1", True)])
def test_the_traceback_variable_treats_empty_and_zero_as_off(
        tmp_path, monkeypatch, value, chained):
    # Any non-empty value was truthy, so OECLUSTER_CLI_TRACEBACK=0 turned
    # tracebacks on.
    monkeypatch.setenv("OECLUSTER_CLI_TRACEBACK", value)
    path = _patched_npz(tmp_path, "facts", facts_json=np.array("42"))
    result = _run("cluster", path, "--algorithm", "butina",
                  "--set", "threshold=1.0")
    assert result.exit_code == 1
    if chained:
        assert isinstance(result.exception, _cli_input.InputError)
    else:
        assert result.exception.__class__ is SystemExit


def test_an_npz_label_that_is_not_utf8_is_refused(tmp_path):
    # The .npz path takes the labels the Python API was handed, so a bytes
    # label need not be UTF-8. The raw path already refuses a non-UTF-8
    # title in its sidecar; refusing here keeps the two from disagreeing
    # about the same bad title.
    path = _patched_npz(tmp_path, "raw",
                        labels=np.array([b"\xff", b"b", b"c", b"d"]))
    with pytest.raises(_cli_input.InputError, match="UTF-8"):
        _cli_input.load(path, warn=lambda message: None)


def test_two_undecodable_byte_labels_are_refused_not_collapsed(tmp_path):
    # errors="replace" maps every undecodable byte onto the one replacement
    # character, so b"\xff" and b"\xfe" exported the same id and the result
    # could no longer be joined back to the user's items.
    path = _patched_npz(tmp_path, "raw",
                        labels=np.array([b"\xff", b"\xfe", b"c", b"d"]))
    out = str(tmp_path / "r.csv")
    result = _run("cluster", path, "--algorithm", "butina",
                  "--set", "threshold=1.0", "--output", out)
    assert result.exit_code == 1
    assert "UTF-8" in result.output
    assert not Path(out).exists()


def test_the_renderer_never_collapses_two_byte_labels():
    # The loader refuses undecodable bytes, so this is the second line of
    # defence rather than the first: whatever reaches the renderer, two
    # distinct labels must not leave it as one id.
    assert _cli_render.text(b"\xff") != _cli_render.text(b"\xfe")
    assert _cli_render.text(b"a") == "a"


@pytest.mark.parametrize("threads", ["18446744073709551616", str(2**100)])
def test_an_oversized_thread_count_is_refused_before_the_load(tmp_path,
                                                              threads):
    # _thread_count bounds below but not above, so the value reached a
    # native size_t setter after the matrix had been loaded and escaped as
    # an OverflowError traceback. The input here does not exist, so a
    # refusal that reads the file first would exit 1 instead of 2.
    result = _run("cluster", str(tmp_path / "gone.npz"), "--algorithm",
                  "butina", "--set", "threshold=1.0", "--threads", threads)
    assert result.exit_code == 2
    assert "size_t" in result.output


def test_the_largest_in_range_thread_count_is_still_accepted(tmp_path):
    # The bound must refuse nothing the native parameter can hold, or it is
    # second-guessing the library rather than guarding the conversion.
    registry = _cli_registry.build()
    spec = _cli_main._spec("butina", ["threshold=1.0"], registry,
                           oecluster._SIZE_T_MAX, False)
    assert spec.options["num_threads"] == oecluster._SIZE_T_MAX


def test_an_empty_output_path_is_refused(tmp_path):
    # `--output "$OUT"` with an unset variable discarded the whole run:
    # both the validation and the write tested truthiness, so "" read as
    # "no output requested" and the command exited 0 having written
    # nothing. Under --quiet that was completely silent.
    result = _run("cluster", _matrix(tmp_path), "--algorithm", "butina",
                  "--set", "threshold=1.0", "--quiet", "--output", "")
    assert result.exit_code == 2
    assert "empty" in result.output
    assert "clusters=" not in result.output


def test_omitting_output_entirely_is_still_fine(tmp_path):
    # The other half of the distinction: not supplied is not the same as
    # supplied empty, and only the second is an error.
    result = _run("cluster", _matrix(tmp_path), "--algorithm", "butina",
                  "--set", "threshold=1.0")
    assert result.exit_code == 0


@pytest.mark.parametrize("algorithm, assignment", [
    ("k_medoids", "n_clusters=18446744073709551615"),
    ("hdbscan", "min_cluster_size=18446744073709551615"),
    ("agglomerative", "n_clusters=18446744073709551615"),
])
def test_a_set_integer_at_the_native_maximum_exits_cleanly(tmp_path,
                                                           algorithm,
                                                           assignment):
    # The companion to the --threads bound, and the reason --set needs no
    # equivalent: the registry caps magnitude at the widest native integer
    # and the library range-checks each option by name, so the overflow
    # --threads produced cannot be reached this way. Pinned so that
    # widening either cap reopens it as a failure here rather than as a
    # traceback in front of a user.
    result = _run("cluster", _matrix(tmp_path), "--algorithm", algorithm,
                  "--quiet", "--set", assignment)
    assert result.exit_code == 2
    assert result.exception.__class__ is SystemExit


def test_distinct_array_labels_that_render_alike_are_refused(tmp_path):
    # The class, not the type: numpy prints array([1.000000001]) as "[1.]",
    # so three distinct labels exported one id. Guarding bytes did nothing
    # for this, and a third type would have found the same hole.
    path = _patched_npz(tmp_path, "arrays", labels=np.array(
        [[1.000000001], [1.000000002], [3.0], [4.0]]))
    out = str(tmp_path / "r.csv")
    result = _run("cluster", path, "--algorithm", "butina",
                  "--set", "threshold=1.0", "--output", out)
    assert result.exit_code == 1
    assert "render" in result.output
    assert not Path(out).exists()


def test_two_genuinely_duplicate_labels_still_export_two_rows(tmp_path):
    # The converse, and it is deliberate: two items really called mol1 share
    # an id on purpose. The rule is "distinct labels that render alike", not
    # "duplicate ids", and a fix that refused this would be wrong.
    path = _patched_npz(tmp_path, "dup",
                        labels=np.array(["mol1", "mol1", "c", "d"]))
    out = str(tmp_path / "r.csv")
    result = _run("cluster", path, "--algorithm", "butina",
                  "--set", "threshold=1.0", "--output", out)
    assert result.exit_code == 0
    lines = Path(out).read_text(encoding="utf-8").splitlines()
    assert [line.split(",")[0] for line in lines[1:]] == ["mol1", "mol1",
                                                          "c", "d"]


@pytest.mark.parametrize("labels, ids", [
    (["a", "b", "c", "d"], ["a", "b", "c", "d"]),
    ([b"a", b"b", b"c", b"d"], ["a", "b", "c", "d"]),
    # Scalars, at the precision that defeats the array rendering: numpy
    # scalars print in full, so these must not be caught by the guard.
    ([1.000000001, 1.000000002, 3.0, 4.0],
     ["1.000000001", "1.000000002", "3.0", "4.0"]),
    ([1, 2, 3, 4], ["1", "2", "3", "4"]),
])
def test_label_types_that_render_distinctly_are_accepted(tmp_path, labels,
                                                         ids):
    path = _patched_npz(tmp_path, "kinds", labels=np.array(labels))
    out = str(tmp_path / "r.csv")
    result = _run("cluster", path, "--algorithm", "butina",
                  "--set", "threshold=1.0", "--output", out)
    assert result.exit_code == 0, result.output
    lines = Path(out).read_text(encoding="utf-8").splitlines()
    assert [line.split(",")[0] for line in lines[1:]] == ids


def test_a_label_that_cannot_be_written_is_refused_before_the_run(tmp_path):
    # A lone surrogate is a str, renders as itself and stays distinct, so
    # every other check passes -- and then the UTF-8 encode inside the
    # writer fails, after the header row has already reached the disk.
    path = _patched_npz(tmp_path, "sur",
                        labels=np.array(["a", "\udcff", "c", "d"]))
    out = str(tmp_path / "r.csv")
    result = _run("cluster", path, "--algorithm", "butina",
                  "--set", "threshold=1.0", "--output", out)
    assert result.exit_code == 1
    assert "UTF-8" in result.output
    assert not Path(out).exists()


@pytest.mark.parametrize("name, payload, error", [
    ("r.csv", {"header": ["id", "label"], "rows": [("a", 0), ("\udcff", 1)]},
     UnicodeEncodeError),
    ("r.json", {"document": {"x": object()}}, TypeError),
])
def test_a_failed_write_leaves_no_partial_file(tmp_path, name, payload,
                                               error):
    # A header-only CSV on disk reads as a run that succeeded and found
    # nothing, which is worse than no file at all.
    out = str(tmp_path / name)
    with pytest.raises(error):
        _cli_render.write_output(out, payload)
    assert not Path(out).exists()


def test_check_ids_states_the_rule_it_enforces():
    # Read as the specification: equal labels may share an id, different
    # labels may not.
    assert _cli_render.check_ids(["a", "a", "b"]) == ["a", "a", "b"]
    assert _cli_render.check_ids([]) == []
    with pytest.raises(ValueError, match="render"):
        _cli_render.check_ids([np.array([1.000000001]),
                               np.array([1.000000002])])
    # Multi-element arrays that really are equal cannot be settled cheaply,
    # and the guard refuses what it cannot settle rather than exporting an
    # id two items might share.
    assert _cli_render.check_ids([np.array([1.0, 2.0])]) == ["[1. 2.]"]


def _nonmetric(tmp_path, n=12):
    """A matrix with a proven triangle inequality violation."""
    values = np.array([0.1 if i // (n // 2) == j // (n // 2) else 0.2
                       for i in range(n) for j in range(i + 1, n)])
    values[0] = 0.9
    matrix = oecluster.SymmetricDistanceMatrix.from_condensed(
        values, labels=[f"m{i}" for i in range(n)])
    path = str(tmp_path / "nm.npz")
    matrix.to_file(path)
    return path


def test_select_parameter_sweeps_and_names_a_winner(tmp_path):
    result = _run("select-parameter", _matrix(tmp_path), "--algorithm",
                  "butina", "--parameter", "threshold", "--values",
                  "0.8,1.0,1.2")
    assert result.exit_code == 0, result.output
    assert "winner threshold=" in result.stdout


def test_select_parameter_does_not_require_the_swept_option(tmp_path):
    # The normal form omits --set threshold=, because the sweep supplies it.
    result = _run("select-parameter", _matrix(tmp_path), "--algorithm",
                  "butina", "--parameter", "threshold", "--values", "0.8,1.0")
    assert result.exit_code == 0, result.output


def test_passing_the_swept_option_again_is_an_error(tmp_path):
    result = _run("select-parameter", _matrix(tmp_path), "--algorithm",
                  "butina", "--parameter", "threshold", "--values", "0.8,1.0",
                  "--set", "threshold=0.9")
    assert result.exit_code == 2
    assert "swept parameter" in result.output


def test_select_parameter_csv_marks_the_winner(tmp_path):
    out = str(tmp_path / "s.csv")
    _run("select-parameter", _matrix(tmp_path), "--algorithm", "butina",
         "--parameter", "threshold", "--values", "0.8,1.0,1.2",
         "--output", out)
    lines = Path(out).read_text(encoding="utf-8").splitlines()
    assert lines[0].startswith("winner,")
    assert sum(1 for line in lines[1:] if line.startswith("true,")) == 1


def test_select_parameter_json_carries_the_winner_index(tmp_path):
    out = str(tmp_path / "s.json")
    _run("select-parameter", _matrix(tmp_path), "--algorithm", "butina",
         "--parameter", "threshold", "--values", "0.8,1.0", "--output", out)
    document = json.loads(Path(out).read_text(encoding="utf-8"))
    assert document["command"] == "select-parameter"
    assert isinstance(document["result"]["winner_index"], int)


def test_an_unknown_swept_option_is_refused(tmp_path):
    result = _run("select-parameter", _matrix(tmp_path), "--algorithm",
                  "butina", "--parameter", "nosuch", "--values", "1")
    assert result.exit_code == 2
    assert "no option" in result.output


def test_a_misspelled_swept_option_is_named_before_requiredness(tmp_path):
    # Reversing the two checks reports "butina requires --set threshold=…"
    # for a misspelled --parameter, because the swept name excused the real
    # option from requiredness.
    result = _run("select-parameter", _matrix(tmp_path), "--algorithm",
                  "butina", "--parameter", "treshold", "--values", "1.0")
    assert result.exit_code == 2
    assert "no option" in result.output
    assert "requires" not in result.output


@pytest.mark.parametrize("algorithm, parameter, values, needle", [
    ("butina", "threshold", "0.8,abc", "number"),
    ("butina", "threshold", "inf", "finite"),
    ("butina", "threshold", "nan", "finite"),
    ("butina", "threshold", "1e999", "finite"),
    ("agglomerative", "n_clusters", "1e999", "range"),
    ("agglomerative", "n_clusters", "2.5", "integer"),
])
def test_a_bad_sweep_value_is_refused_before_the_load(tmp_path, algorithm,
                                                      parameter, values,
                                                      needle):
    # --values does not pass through _cli_registry.resolve(), so without a
    # coercion of its own the grid reached the library untyped: a string
    # threshold, an infinity, or a mistyped exponent that int() expands into
    # a billion digits. The input here does not exist, so a refusal that
    # loaded the matrix first would exit 1 instead of 2.
    result = _run("select-parameter", str(tmp_path / "gone.npz"),
                  "--algorithm", algorithm, "--parameter", parameter,
                  "--values", values)
    assert result.exit_code == 2
    assert needle in result.output


def test_select_parameter_routes_nonmetric_to_the_scorer_too(tmp_path):
    # cluster_report runs its own metric gate and select_parameter forwards
    # report_options to it unchanged, so a flag that reaches only the spec
    # clusters fine and then fails in scoring.
    result = _run("select-parameter", _nonmetric(tmp_path), "--algorithm",
                  "butina", "--parameter", "threshold", "--values", "0.15,0.5",
                  "--allow-nonmetric")
    assert result.exit_code == 0, result.output
    assert "winner threshold=" in result.stdout


def test_a_nonmetric_matrix_without_the_flag_is_refused(tmp_path):
    result = _run("select-parameter", _nonmetric(tmp_path), "--algorithm",
                  "butina", "--parameter", "threshold", "--values", "0.15,0.5")
    assert result.exit_code == 2
    assert "triangle inequality" in result.output


def test_stability_scores_every_cluster(tmp_path):
    result = _run("stability", _matrix(tmp_path), "--algorithm", "butina",
                  "--set", "threshold=1.0", "--resamples", "8")
    assert result.exit_code == 0, result.output
    assert "resamples=8" in result.stdout
    assert "mean_jaccard=" in result.stdout


def test_stability_json_exports_reference_labels(tmp_path):
    # Named explicitly: ClusterStability also has per-resample labels paired
    # with indices, and a bare "labels" would be ambiguous between them.
    out = str(tmp_path / "st.json")
    _run("stability", _matrix(tmp_path), "--algorithm", "butina",
         "--set", "threshold=1.0", "--resamples", "6", "--output", out)
    document = json.loads(Path(out).read_text(encoding="utf-8"))
    assert len(document["result"]["reference_labels"]) == 40
    assert "labels" not in document["result"]


def test_stability_csv_carries_the_record_columns(tmp_path):
    out = str(tmp_path / "st.csv")
    _run("stability", _matrix(tmp_path), "--algorithm", "butina",
         "--set", "threshold=1.0", "--resamples", "6", "--output", out)
    lines = Path(out).read_text(encoding="utf-8").splitlines()
    assert lines[0] == "label,size,mean_jaccard,dissolved,recovered,evaluated"
    assert len(lines) == 5


def test_a_bad_stability_output_extension_is_refused_before_the_run(tmp_path):
    # The run is minutes long, so the extension has to be refused in front
    # of it; --output is the only destination the contract applies to.
    result = _run("stability", _matrix(tmp_path), "--algorithm", "butina",
                  "--set", "threshold=1.0", "--resamples", "4", "--output",
                  str(tmp_path / "st.txt"))
    assert result.exit_code == 2
    assert ".csv or .json" in result.output


def _result(tmp_path, name, *args):
    """Run a command with a JSON --output and return the result document."""
    out = str(tmp_path / name)
    result = _run(*args, "--quiet", "--output", out)
    assert result.exit_code == 0, result.output
    return json.loads(Path(out).read_text(encoding="utf-8"))["result"]


def test_the_criterion_reaches_the_scorer(tmp_path):
    # The criterion names the second column, so dropping criterion= at the
    # call site silently reverts the whole ranking to silhouette.
    document = _result(tmp_path, "c.json", "select-parameter",
                       _matrix(tmp_path), "--algorithm", "butina",
                       "--parameter", "threshold", "--values", "0.8,1.0",
                       "--criterion", "davies_bouldin_medoid")
    assert document["columns"][1] == "davies_bouldin_medoid"


def test_max_noise_fraction_rejects_only_the_noisy_value(tmp_path):
    # dbscan at eps 0.2 leaves 33 of 40 items as noise and at 0.5 leaves
    # none, so the bound has to discriminate between the two rows rather
    # than reject or keep both.
    document = _result(tmp_path, "n.json", "select-parameter",
                       _matrix(tmp_path), "--algorithm", "dbscan",
                       "--parameter", "eps", "--values", "0.2,0.5",
                       "--set", "min_samples=5",
                       "--max-noise-fraction", "0.5")
    assert "max_noise_fraction" in document["rows"][0][-1]
    assert document["rows"][1][-1] is None
    assert document["winner_index"] == 1


def test_max_clusters_rejects_an_oversized_partition(tmp_path):
    # The rejection is read by name, not as a bare "no winner": swapping
    # min_clusters and max_clusters at the call site leaves every row
    # eligible here and every row rejected in the test below, so each one
    # catches the swap the other does not.
    document = _result(tmp_path, "x.json", "select-parameter",
                       _matrix(tmp_path), "--algorithm", "butina",
                       "--parameter", "threshold", "--values", "0.8,1.0",
                       "--max-clusters", "3")
    assert all("max_clusters" in row[-1] for row in document["rows"])
    assert document["winner_index"] is None


def test_no_eligible_value_is_a_clean_outcome(tmp_path):
    # winner is legitimately None when every row was ineligible. The run
    # succeeded and the table still shows why, so this is exit 0 with a
    # null winner_index, not a failure.
    out = str(tmp_path / "w.json")
    result = _run("select-parameter", _matrix(tmp_path), "--algorithm",
                  "butina", "--parameter", "threshold", "--values", "0.8,1.0",
                  "--min-clusters", "99", "--output", out)
    assert result.exit_code == 0, result.output
    assert "no winner" in result.stdout
    document = json.loads(Path(out).read_text(encoding="utf-8"))["result"]
    assert document["winner_index"] is None
    assert all("min_clusters" in row[-1] for row in document["rows"])


def _noisy_stability(tmp_path, *args):
    """Run stability over a reference partition that has real noise."""
    result = _run("stability", _matrix(tmp_path), "--algorithm", "dbscan",
                  "--set", "eps=0.3", "--set", "min_samples=5",
                  "--resamples", "6", *args)
    assert result.exit_code == 0, result.output
    return result.stdout


def _statistic(summary, name):
    # Matched with the "=" attached: "mean_jaccard" is also a column header
    # in the table above the summary, and that bare word matches first.
    # The summary wraps at the console width, but rich breaks on spaces, so
    # each key=value stays one token.
    field = next(part for part in summary.split()
                 if part.startswith(f"{name}="))
    return field.split("=", 1)[1]


def test_the_seed_changes_the_resampling_draw(tmp_path):
    first = _noisy_stability(tmp_path, "--seed", "1")
    second = _noisy_stability(tmp_path, "--seed", "2")
    assert _statistic(first, "mean_jaccard") != _statistic(second,
                                                           "mean_jaccard")


def test_the_fraction_changes_how_much_is_drawn(tmp_path):
    half = _noisy_stability(tmp_path, "--seed", "1", "--fraction", "0.5")
    most = _noisy_stability(tmp_path, "--seed", "1", "--fraction", "0.9")
    assert _statistic(half, "mean_jaccard") != _statistic(most,
                                                          "mean_jaccard")


def test_the_noise_mode_reaches_the_agreement(tmp_path):
    # --noise governs the agreement and nothing else: the Jaccard matching
    # never treats noise as a cluster, so mean_jaccard is identical across
    # all three modes and only mean_agreement moves.
    summaries = {mode: _noisy_stability(tmp_path, "--seed", "0", "--noise",
                                        mode)
                 for mode in ("singletons", "grouped", "excluded")}
    agreements = {_statistic(summary, "mean_agreement")
                  for summary in summaries.values()}
    assert len(agreements) == 3
    assert len({_statistic(summary, "mean_jaccard")
                for summary in summaries.values()}) == 1


def test_the_resample_count_reaches_the_library(tmp_path):
    # The summary echoes --resamples straight back, so it proves nothing on
    # its own; the evaluated column is the library's own count.
    out = str(tmp_path / "e.csv")
    _run("stability", _matrix(tmp_path), "--algorithm", "butina", "--set",
         "threshold=1.0", "--resamples", "7", "--fraction", "0.9",
         "--quiet", "--output", out)
    lines = Path(out).read_text(encoding="utf-8").splitlines()[1:]
    assert [line.split(",")[-1] for line in lines] == ["7"] * 4


@pytest.mark.parametrize("values, needle", [
    ("0.8,,1.0", "field 2 of 3"),
    (",0.8,1.0", "field 1 of 3"),
    ("0.8,1.0,", "field 3 of 3"),
    ("0.8, ,1.0", "field 2 of 3"),
    (",", "field 1 of 2"),
    ("", "unset shell variable"),
    ("   ", "unset shell variable"),
])
def test_an_empty_sweep_field_is_refused_by_position(tmp_path, values,
                                                     needle):
    # Dropping empty fields instead meant `--values 0.8,,1.0` evaluated two
    # points and named a winner from a grid the user never wrote; the
    # two-point run and the three-point run were indistinguishable in the
    # output. The input here does not exist, so exit 2 rather than 1 proves
    # the refusal precedes the load.
    result = _run("select-parameter", str(tmp_path / "gone.npz"),
                  "--algorithm", "butina", "--parameter", "threshold",
                  "--values", values)
    assert result.exit_code == 2
    assert needle in result.output


def test_every_sweep_field_written_is_a_point_evaluated(tmp_path):
    # The positive half: three fields in, three rows out.
    document = _result(tmp_path, "g.json", "select-parameter",
                       _matrix(tmp_path), "--algorithm", "butina",
                       "--parameter", "threshold", "--values", "0.8,1.0,1.2")
    assert len(document["rows"]) == 3


def test_an_eligible_but_unscorable_grid_is_named_apart(tmp_path):
    # "no eligible" and "no scorable" call for opposite responses -- loosen
    # the bounds, or pick a criterion the input supports -- so the one
    # message for both was telling half the users the wrong thing. dbscan
    # at eps 0.2 finds a single cluster, whose silhouette is undefined.
    result = _run("select-parameter", _matrix(tmp_path), "--algorithm",
                  "dbscan", "--parameter", "eps", "--values", "0.2",
                  "--set", "min_samples=5")
    assert result.exit_code == 0, result.output
    assert "no scorable" in result.stdout


def _stability_payload(tmp_path, suffix, mode):
    out = str(tmp_path / f"{mode}.{suffix}")
    result = _run("stability", _matrix(tmp_path), "--algorithm", "dbscan",
                  "--set", "eps=0.3", "--set", "min_samples=5",
                  "--resamples", "6", "--noise", mode, "--quiet",
                  "--output", out)
    assert result.exit_code == 0, result.output
    return Path(out).read_text(encoding="utf-8")


def test_the_noise_mode_reaches_the_json_output(tmp_path):
    # --quiet --output is the whole point of the command, and all three
    # modes wrote byte-identical JSON: the one statistic --noise governs
    # was absent from the machine-readable surface, so the option did
    # nothing at all in a script.
    modes = ("singletons", "grouped", "excluded")
    documents = {mode: json.loads(_stability_payload(tmp_path, "json", mode))
                 for mode in modes}
    assert len({document["result"]["mean_agreement"]
                for document in documents.values()}) == 3
    for mode, document in documents.items():
        # The mode travels with the statistic because it is its unit: the
        # same partitions score differently under each.
        assert document["result"]["noise"] == mode


def test_the_stability_csv_stays_the_per_cluster_table(tmp_path):
    # Deliberate, not an oversight: the scalars are JSON-only, as consensus
    # already keeps its own mean_agreement out of a one-row-per-item CSV.
    # Asserted so the omission is a decision on the record rather than
    # something a later edit discovers by accident.
    payloads = {_stability_payload(tmp_path, "csv", mode)
                for mode in ("singletons", "grouped", "excluded")}
    assert len(payloads) == 1
    assert payloads.pop().splitlines()[0] == (
        "label,size,mean_jaccard,dissolved,recovered,evaluated")


def test_an_ineligible_algorithm_reads_the_same_in_every_command(tmp_path):
    # select-parameter duplicated only the membership half of the guard and
    # ran it before the eligibility check, so --algorithm murcko was
    # answered "murcko has no option 'threshold' to sweep".
    path = _matrix(tmp_path)
    for args in (("cluster", path, "--algorithm", "murcko"),
                 ("select-parameter", path, "--algorithm", "murcko",
                  "--parameter", "threshold", "--values", "1"),
                 ("stability", path, "--algorithm", "murcko")):
        result = _run(*args)
        assert result.exit_code == 2
        assert "not a distance matrix" in result.output, args


def test_consensus_bootstrap_mode(tmp_path):
    result = _run("consensus", _matrix(tmp_path), "--algorithm", "butina",
                  "--set", "threshold=1.0", "--resamples", "8")
    assert result.exit_code == 0, result.output
    assert "partitions=8" in result.stdout


def test_consensus_cross_algorithm_mode(tmp_path):
    result = _run("consensus", _matrix(tmp_path),
                  "--member", "butina;threshold=1.0",
                  "--member", "dbscan;eps=1.5;min_samples=3")
    assert result.exit_code == 0, result.output
    assert "partitions=2" in result.stdout


@pytest.mark.parametrize("extra, needle", [
    (["--algorithm", "butina", "--set", "threshold=1.0", "--resamples", "4",
      "--member", "butina;threshold=1.0"], "cross-algorithm mode"),
    ([], "either --algorithm"),
])
def test_the_two_modes_are_mutually_exclusive(tmp_path, extra, needle):
    result = _run("consensus", _matrix(tmp_path), *extra)
    assert result.exit_code == 2
    assert needle in result.output


def test_a_malformed_member_is_refused(tmp_path):
    result = _run("consensus", _matrix(tmp_path), "--member", "butina;oops")
    assert result.exit_code == 2
    assert "key=value" in result.output


def test_a_member_typo_is_caught_before_any_work(tmp_path):
    result = _run("consensus", _matrix(tmp_path),
                  "--member", "butina;thresold=1.0")
    assert result.exit_code == 2
    assert "did you mean" in result.output


def test_consensus_json_carries_the_ensemble_spec(tmp_path):
    # A singular algorithm field cannot describe a cross-algorithm ensemble,
    # so spec takes a members shape in that mode.
    out = str(tmp_path / "c.json")
    _run("consensus", _matrix(tmp_path), "--member", "butina;threshold=1.0",
         "--member", "dbscan;eps=1.5;min_samples=3", "--output", out)
    document = json.loads(Path(out).read_text(encoding="utf-8"))
    assert [m["algorithm"] for m in document["spec"]["members"]] == [
        "butina", "dbscan"]
    assert document["result"]["num_partitions"] == 2
    assert isinstance(document["result"]["records"][0]["label"], int)


def test_a_singleton_records_consensus_encodes_as_null(tmp_path):
    out = str(tmp_path / "c.json")
    _run("consensus", _matrix(tmp_path), "--member", "butina;threshold=0.05",
         "--output", out)
    document = json.loads(Path(out).read_text(encoding="utf-8"))
    values = [r["cluster_consensus"] for r in document["result"]["records"]]
    assert None in values


def _consensus_payload(tmp_path, mode):
    # butina at two thresholds plus dbscan disagree, which is what makes the
    # three modes score differently; an ensemble that agrees scores 1.0
    # under all three and the comparison below goes vacuous.
    out = str(tmp_path / f"{mode}.json")
    result = _run("consensus", _matrix(tmp_path), "--member",
                  "butina;threshold=1.0", "--member", "butina;threshold=0.3",
                  "--member", "dbscan;eps=0.3", "--noise", mode, "--quiet",
                  "--output", out)
    assert result.exit_code == 0, result.output
    return json.loads(Path(out).read_text(encoding="utf-8"))


def test_the_consensus_noise_mode_reaches_the_json_output(tmp_path):
    # Three distinct values, not merely the echoed mode: --noise is copied
    # straight from the parameter into the document, so asserting only
    # result["noise"] passes with the kernel argument deleted. This is the
    # shape the stability test already uses, for the same reason.
    documents = {mode: _consensus_payload(tmp_path, mode)
                 for mode in ("singletons", "grouped", "excluded")}
    assert len({document["result"]["mean_agreement"]
                for document in documents.values()}) == 3
    for mode, document in documents.items():
        # The mode travels beside the statistic as its unit: the same
        # partitions score differently under each, so a mean_agreement
        # recorded without its convention cannot be compared across runs.
        assert document["result"]["noise"] == mode


def test_the_consensus_csv_stays_one_row_per_item(tmp_path):
    out = str(tmp_path / "c.csv")
    result = _run("consensus", _matrix(tmp_path), "--member",
                  "butina;threshold=1.0", "--quiet", "--output", out)
    assert result.exit_code == 0, result.output
    lines = Path(out).read_text(encoding="utf-8").splitlines()
    assert lines[0] == "id,label"
    assert len(lines) == 41
    assert [line.split(",")[0] for line in lines[1:]] == [
        f"m{i}" for i in range(40)]


def test_the_panel_tells_two_members_of_one_algorithm_apart(tmp_path):
    # Rendering member names alone printed "butina, butina" for an ensemble
    # whose entire content is that the two members differ.
    result = _run("consensus", _matrix(tmp_path), "--member",
                  "butina;threshold=1.0", "--member", "butina;threshold=0.5")
    assert result.exit_code == 0, result.output
    assert "butina threshold=1.0" in result.stdout
    assert "butina threshold=0.5" in result.stdout


@pytest.mark.parametrize("extra, needle", [
    (["--member", "butina;threshold=1.0", "--set", "threshold=1.0"],
     "bootstrap mode"),
    (["--algorithm", "butina", "--set", "threshold=1.0", "--resamples", "0"],
     "positive"),
    (["--member", "butina;threshold=1.0", "--threshold", "1.5"], "[0, 1]"),
    (["--member", "butina;"], "empty option"),
    (["--member", ""], "needs an algorithm name"),
])
def test_consensus_validation_exits_two(tmp_path, extra, needle):
    # The input does not exist, so exit 2 rather than 1 also proves each
    # refusal precedes the load.
    result = _run("consensus", str(tmp_path / "gone.npz"), *extra)
    assert result.exit_code == 2
    assert needle in result.output


def _capture(monkeypatch, name):
    """Record the arguments one library entry point is called with."""
    seen = {}
    original = getattr(oecluster, name)

    def recorder(*args, **kwargs):
        seen.update(kwargs)
        seen["_args"] = args
        return original(*args, **kwargs)

    monkeypatch.setattr(oecluster, name, recorder)
    return seen


def _specs_used(monkeypatch):
    """Record the options of every ClusteringSpec actually run."""
    used = []
    original = oecluster.ClusteringSpec.run

    def recorder(self, items, *args, **kwargs):
        used.append(dict(self.options))
        return original(self, items, *args, **kwargs)

    monkeypatch.setattr(oecluster.ClusteringSpec, "run", recorder)
    return used


def test_threads_reaches_the_spec_for_cluster(tmp_path, monkeypatch):
    # Asserting only exit_code would pass with the routing deleted.
    used = _specs_used(monkeypatch)
    result = _run("cluster", _matrix(tmp_path), "--algorithm", "butina",
                  "--set", "threshold=1.0", "--threads", "3")
    assert result.exit_code == 0, result.output
    assert used == [{"threshold": 1.0, "num_threads": 3}]


def test_threads_reaches_stability_as_well_as_the_spec(tmp_path, monkeypatch):
    # cluster_stability takes its own num_threads for matrix gathering while
    # the algorithm's threads live in the spec, so one flag has two targets.
    seen = _capture(monkeypatch, "cluster_stability")
    result = _run("stability", _matrix(tmp_path), "--algorithm", "butina",
                  "--set", "threshold=1.0", "--resamples", "4",
                  "--threads", "3")
    assert result.exit_code == 0, result.output
    assert seen["num_threads"] == 3
    assert seen["_args"][0].options["num_threads"] == 3


def test_threads_reaches_every_consensus_destination(tmp_path, monkeypatch):
    # Cross-algorithm mode: the kernel and every member spec.
    seen = _capture(monkeypatch, "consensus")
    used = _specs_used(monkeypatch)
    result = _run("consensus", _matrix(tmp_path), "--member",
                  "butina;threshold=1.0", "--member", "dbscan;eps=1.5",
                  "--threads", "2")
    assert result.exit_code == 0, result.output
    assert seen["num_threads"] == 2
    assert [options["num_threads"] for options in used] == [2, 2]


def test_threads_reaches_consensus_bootstrap_destinations(tmp_path,
                                                          monkeypatch):
    # Bootstrap mode: the kernel, cluster_stability, and the spec.
    agreed = _capture(monkeypatch, "consensus")
    resampled = _capture(monkeypatch, "cluster_stability")
    result = _run("consensus", _matrix(tmp_path), "--algorithm", "butina",
                  "--set", "threshold=1.0", "--resamples", "4",
                  "--threads", "2")
    assert result.exit_code == 0, result.output
    assert agreed["num_threads"] == 2
    assert resampled["num_threads"] == 2
    assert resampled["_args"][0].options["num_threads"] == 2


def test_select_parameter_does_not_pass_num_threads(tmp_path, monkeypatch):
    # select_parameter has no num_threads parameter; passing one would be a
    # TypeError, so the flag must reach it only through the spec.
    seen = _capture(monkeypatch, "select_parameter")
    result = _run("select-parameter", _matrix(tmp_path), "--algorithm",
                  "butina", "--parameter", "threshold", "--values",
                  "0.8,1.0", "--threads", "3")
    assert result.exit_code == 0, result.output
    assert "num_threads" not in seen
    assert seen["_args"][0].options["num_threads"] == 3


@pytest.mark.parametrize("algorithm, option", [
    ("butina", "threshold=0.5"),
    ("dbscan", "eps=0.5"),
    ("agglomerative", "n_clusters=2"),
    ("hdbscan", "min_cluster_size=2"),
])
def test_allow_nonmetric_is_required_and_sufficient(tmp_path, algorithm,
                                                    option):
    # On a matrix with proven violations, each of the four accepting
    # algorithms must refuse without the flag and succeed with it. A metric
    # fixture would pass either way and prove nothing.
    path = _nonmetric(tmp_path)
    refused = _run("cluster", path, "--algorithm", algorithm, "--set", option)
    assert refused.exit_code != 0
    allowed = _run("cluster", path, "--algorithm", algorithm, "--set", option,
                   "--allow-nonmetric")
    assert allowed.exit_code == 0, allowed.output


def test_k_medoids_accepts_the_flag_without_receiving_it(tmp_path,
                                                         monkeypatch):
    # k_medoids never appeals to the triangle inequality and has no such
    # option, so the flag must be accepted and simply not forwarded.
    used = _specs_used(monkeypatch)
    result = _run("cluster", _nonmetric(tmp_path), "--algorithm", "k_medoids",
                  "--set", "n_clusters=2", "--allow-nonmetric")
    assert result.exit_code == 0, result.output
    assert "allow_nonmetric" not in used[0]


def test_setting_allow_nonmetric_points_at_the_flag(tmp_path):
    result = _run("cluster", _matrix(tmp_path), "--algorithm", "butina",
                  "--set", "threshold=1.0", "--set", "allow_nonmetric=true")
    assert result.exit_code == 2
    assert "--allow-nonmetric" in result.output


def test_allow_nonmetric_reaches_the_selection_scorer(tmp_path, monkeypatch):
    # cluster_report runs its own metric gate and select_parameter forwards
    # report_options to it unchanged, so the flag needs both destinations.
    seen = _capture(monkeypatch, "select_parameter")
    result = _run("select-parameter", _nonmetric(tmp_path), "--algorithm",
                  "butina", "--parameter", "threshold", "--values",
                  "0.15,0.5", "--allow-nonmetric")
    assert result.exit_code == 0, result.output
    assert seen["report_options"] == {"allow_nonmetric": True}
    assert seen["_args"][0].options["allow_nonmetric"] is True


def test_selection_without_the_flag_is_refused_on_a_nonmetric_matrix(
        tmp_path):
    result = _run("select-parameter", _nonmetric(tmp_path), "--algorithm",
                  "butina", "--parameter", "threshold", "--values",
                  "0.15,0.5")
    assert result.exit_code != 0


def test_nothing_is_read_before_validation_finishes(tmp_path, monkeypatch):
    # The recording stub the spec asks for: if validation is ordered
    # correctly, a bad option is reported without the loader ever running.
    def refuse(*args, **kwargs):
        raise AssertionError("the matrix was loaded before validation")

    monkeypatch.setattr("oecluster._cli_main._cli_input.load", refuse)
    for args in (
        ["cluster", "x.npz", "--algorithm", "butina", "--set", "nope=1"],
        ["cluster", "x.npz", "--algorithm", "nosuch", "--set", "a=1"],
        ["consensus", "x.npz", "--member", "butina;thresold=1.0"],
        ["consensus", "x.npz", "--algorithm", "butina", "--set",
         "threshold=1.0", "--resamples", "4", "--member", "dbscan;eps=1"],
        ["select-parameter", "x.npz", "--algorithm", "butina",
         "--parameter", "nosuch", "--values", "1"],
    ):
        result = _run(*args)
        assert result.exit_code == 2, args


@pytest.mark.parametrize("command, extra", [
    ("cluster", ["--algorithm", "butina", "--set", "threshold=1.0"]),
    ("stability", ["--algorithm", "butina", "--set", "threshold=1.0"]),
    ("select-parameter", ["--algorithm", "butina", "--parameter",
                          "threshold", "--values", "0.8,1.0"]),
    ("consensus", ["--member", "butina;threshold=1.0"]),
])
def test_no_command_will_overwrite_the_input_sidecar(tmp_path, monkeypatch,
                                                     command, extra):
    # With input d.npy the obvious --output d.json truncates the sidecar the
    # run depends on. Every writing command must refuse, before any work.
    def refuse(*args, **kwargs):
        raise AssertionError("work started before the path check")

    monkeypatch.setattr("oecluster._cli_main._cli_input.load", refuse)
    result = _run(command, str(tmp_path / "d.npy"), *extra,
                  "--output", str(tmp_path / "d.json"))
    assert result.exit_code == 2
    assert "sidecar" in result.output


@pytest.mark.parametrize("alias", ["same-string", "dot-segment", "symlink"])
def test_consensus_refuses_output_and_mmap_on_one_file(tmp_path, monkeypatch,
                                                       alias):
    # Two identical strings never exercise realpath, so the aliases are the
    # point: delete the normalisation and the last two cases stop failing.
    # pathlib would collapse a "." component on construction, so the
    # dot-segment case is built with os.path.join.
    def refuse(*args, **kwargs):
        raise AssertionError("work started before the path check")

    monkeypatch.setattr("oecluster._cli_main._cli_input.load", refuse)
    target = tmp_path / "x.json"
    if alias == "same-string":
        other = str(target)
    elif alias == "dot-segment":
        other = os.path.join(str(tmp_path), ".", "x.json")
    else:
        target.write_text("", encoding="utf-8")
        link = tmp_path / "link.json"
        link.symlink_to(target)
        other = str(link)
    result = _run("consensus", str(tmp_path / "m.npz"), "--member",
                  "butina;threshold=1.0", "--output", str(target),
                  "--mmap", other)
    assert result.exit_code == 2
    assert "different files" in result.output


def test_ids_appear_only_when_the_matrix_is_labelled(tmp_path):
    # An unlabelled matrix reports [] rather than None; treating that as
    # real labels wrote a results file with no rows at all.
    values = np.array([0.1, 0.9, 0.2, 0.8, 0.3, 0.7])
    bare = oecluster.SymmetricDistanceMatrix.from_condensed(values)
    path = str(tmp_path / "bare.npz")
    bare.to_file(path)
    out = str(tmp_path / "bare.csv")
    result = _run("cluster", path, "--algorithm", "butina", "--set",
                  "threshold=0.5", "--output", out)
    assert result.exit_code == 0, result.output
    lines = Path(out).read_text(encoding="utf-8").splitlines()
    assert len(lines) == 5, lines          # header plus one row per item
    assert lines[1].startswith("0,")
    document_path = str(tmp_path / "bare.json")
    _run("cluster", path, "--algorithm", "butina", "--set", "threshold=0.5",
         "--output", document_path)
    document = json.loads(Path(document_path).read_text(encoding="utf-8"))
    assert "ids" not in document["result"]
    assert len(document["result"]["labels"]) == 4


def test_translation_hides_the_exception_by_default(tmp_path):
    # The control for the two tests below. Note what is NOT asserted:
    # `result.exception is not None` is true here too, because Click records
    # the SystemExit(1) that a translated error raises. Only the exception's
    # type discriminates.
    result = _run("cluster", str(tmp_path / "gone.npz"), "--algorithm",
                  "butina", "--set", "threshold=1.0")
    assert result.exit_code == 1
    assert not isinstance(result.exception, _cli_input.InputError)


def test_the_traceback_flag_re_raises_the_original_exception(tmp_path):
    result = _run("cluster", str(tmp_path / "gone.npz"), "--algorithm",
                  "butina", "--set", "threshold=1.0", "--traceback")
    assert isinstance(result.exception, _cli_input.InputError)


def test_the_traceback_environment_variable_does_too(tmp_path, monkeypatch):
    monkeypatch.setenv("OECLUSTER_CLI_TRACEBACK", "1")
    result = _run("cluster", str(tmp_path / "gone.npz"), "--algorithm",
                  "butina", "--set", "threshold=1.0")
    assert isinstance(result.exception, _cli_input.InputError)


def test_consensus_writes_json_matching_the_schema(tmp_path):
    # The whole envelope and the whole result key set, not a key count:
    # consensus has the longest schema in the design and a count survives a
    # rename. The ensemble deliberately disagrees, so every statistic below
    # is a real number rather than the 1.0 a unanimous ensemble produces --
    # a constant substituted for any of them fails here.
    path = _matrix(tmp_path)
    out = str(tmp_path / "s.json")
    result = _run("consensus", path, "--member", "butina;threshold=1.0",
                  "--member", "butina;threshold=0.3", "--member",
                  "dbscan;eps=0.3", "--quiet", "--output", out)
    assert result.exit_code == 0, result.output
    document = json.loads(Path(out).read_text(encoding="utf-8"))
    assert sorted(document) == ["command", "input", "result",
                                "schema_version", "spec"]
    assert document["schema_version"] == 1
    assert document["command"] == "consensus"
    assert document["input"] == {"path": path, "num_items": 40,
                                 "orientation": "unknown"}
    outcome = document["result"]
    assert sorted(outcome) == [
        "agreement", "ids", "item_consensus", "labels", "mean_agreement",
        "noise", "num_clusters", "num_noise", "num_partitions", "records",
        "threshold", "unobserved_pairs"]
    assert outcome["ids"] == [f"m{i}" for i in range(40)]
    assert len(outcome["labels"]) == 40
    assert outcome["num_clusters"] == 12
    assert outcome["num_noise"] == sum(1 for v in outcome["labels"] if v < 0)
    assert outcome["num_partitions"] == 3
    # Zero here and non-zero for a bootstrap ensemble below, so no one
    # constant satisfies both.
    assert outcome["unobserved_pairs"] == 0
    assert outcome["threshold"] == 0.5
    assert outcome["noise"] == "singletons"
    # One per member, all different from each other and from 1.0.
    assert len(outcome["agreement"]) == 3
    assert len(set(outcome["agreement"])) == 3
    assert all(0.0 < value < 1.0 for value in outcome["agreement"])
    assert 0.0 < outcome["mean_agreement"] < 1.0
    assert len(outcome["item_consensus"]) == 40
    assert len({str(value) for value in outcome["item_consensus"]}) > 1
    assert len(outcome["records"]) == 12
    assert sorted(outcome["records"][0]) == ["cluster_consensus", "label",
                                             "size"]
    # Cross-checked against the labels rather than counted: a sum of 40 and
    # a length of 40 are both satisfied by forty zeroes, so without this
    # the test accepts a labels array that describes a different partition
    # from the one the records describe.
    sizes = {label: outcome["labels"].count(label)
             for label in set(outcome["labels"])}
    assert sorted(sizes) == list(range(12))
    assert {record["label"]: record["size"]
            for record in outcome["records"]} == sizes


def test_a_bootstrap_ensemble_reports_the_pairs_no_resample_observed(
        tmp_path):
    # The other half of the unobserved_pairs pin: resamples are partial, so
    # some pairs are seen by no member at all, and a constant that satisfied
    # the full-member case above cannot satisfy this one.
    out = str(tmp_path / "u.json")
    result = _run("consensus", _matrix(tmp_path), "--algorithm", "butina",
                  "--set", "threshold=1.0", "--resamples", "6", "--quiet",
                  "--output", out)
    assert result.exit_code == 0, result.output
    document = json.loads(Path(out).read_text(encoding="utf-8"))
    # resamples sits in spec because it describes the ensemble that was
    # built rather than how any one number is read. Nothing else records
    # it: num_partitions would survive its removal unchanged.
    assert document["spec"] == {
        "algorithm": "butina",
        "options": {"threshold": 1.0, "num_threads": 0},
        "resamples": 6}
    outcome = document["result"]
    assert outcome["num_partitions"] == 6
    assert 0 < outcome["unobserved_pairs"] <= 40 * 39 // 2


@pytest.mark.parametrize("threshold, clusters", [("0.0", 1), ("0.6", 12),
                                                 ("1.0", 27)])
def test_the_extraction_threshold_changes_the_partition(tmp_path, threshold,
                                                        clusters):
    # Dropping threshold= at the call site reverts every run to the library
    # default of 0.5, which yields 12 -- the same as the 0.6 row. Only a
    # grid whose outcomes differ from the default catches that, so all
    # three are pinned by value.
    document = _result(tmp_path, f"t{threshold}.json", "consensus",
                       _matrix(tmp_path), "--member", "butina;threshold=1.0",
                       "--member", "butina;threshold=0.3", "--member",
                       "dbscan;eps=0.3", "--threshold", threshold)
    assert document["threshold"] == float(threshold)
    assert document["num_clusters"] == clusters


@pytest.mark.parametrize("extra", [
    ["--member", "butina;threshold=1.0", "--member", "dbscan;eps=1.5"],
    ["--algorithm", "butina", "--set", "threshold=1.0", "--resamples", "4"],
])
def test_the_mmap_file_holds_the_whole_co_association_matrix(tmp_path, extra):
    # Both modes, because each builds the matrix through a different
    # ensemble shape. The file is the condensed upper triangle as float64,
    # so its size is fixed by the item count alone -- a truncated or
    # unwritten mapping is the failure this catches. The values are
    # co-association distances, so a unanimous pair is 0.0 and a spot check
    # of the first bytes can look like an empty file.
    mapped = str(tmp_path / "coassoc.bin")
    result = _run("consensus", _matrix(tmp_path), *extra, "--mmap", mapped,
                  "--quiet")
    assert result.exit_code == 0, result.output
    assert Path(mapped).stat().st_size == 40 * 39 // 2 * 8
    values = np.fromfile(mapped, dtype=np.float64)
    assert len(values) == 40 * 39 // 2
    assert np.all(np.isfinite(values))
    assert values.min() >= 0.0
    assert values.max() <= 1.0


@pytest.mark.parametrize("destination", ["--output", "--mmap"])
def test_a_destination_in_a_missing_directory_is_refused_first(tmp_path,
                                                               monkeypatch,
                                                               destination):
    # open() only fails where the file is written, which for consensus is
    # after every member has clustered and the matrix has been built. The
    # same rule that refuses a bad extension before the run applies here.
    def refuse(*args, **kwargs):
        raise AssertionError("work started before the path check")

    monkeypatch.setattr("oecluster._cli_main._cli_input.load", refuse)
    result = _run("consensus", _matrix(tmp_path), "--member",
                  "butina;threshold=1.0", destination,
                  str(tmp_path / "missing" / "x.json"))
    assert result.exit_code == 2
    assert "directory" in result.output


def test_a_member_option_value_may_not_contain_an_equals_sign(tmp_path):
    # No roster option takes a value containing '=' or ';', so the grammar
    # refuses one with that reason rather than inventing an escape syntax
    # for a case that does not exist.
    result = _run("consensus", str(tmp_path / "gone.npz"), "--member",
                  "butina;threshold=1=0")
    assert result.exit_code == 2
    assert "may not contain" in result.output


@pytest.mark.parametrize("extra", [
    ["--algorithm", "butina", "--set", "threshold=1.0"],
    ["--resamples", "4"],
])
def test_bootstrap_mode_needs_both_halves(tmp_path, extra):
    # Either one alone selects bootstrap mode, so neither can be inferred
    # from the other: a lone --algorithm must not silently resample once,
    # and a lone --resamples has nothing to resample.
    result = _run("consensus", str(tmp_path / "gone.npz"), *extra)
    assert result.exit_code == 2
    assert "needs both" in result.output


def test_a_member_is_told_how_to_supply_its_own_options(tmp_path):
    # The registry's advice is written for --set, which this same command
    # refuses in cross-algorithm mode: `--member butina` was answered
    # "butina requires --set threshold=…", which cannot be followed.
    result = _run("consensus", str(tmp_path / "gone.npz"), "--member",
                  "butina")
    assert result.exit_code == 2
    assert "requires" in result.output
    assert "--member" in result.output
    assert "--set" not in result.output


def test_two_absent_destinations_that_are_one_file_are_refused(tmp_path):
    # The data-loss case, and the reason the identity probe exists.
    # `--output r.json --mmap R.JSON` with neither file present: realpath
    # compares bytes, so the two names looked like two files. The run built
    # the 6240-byte co-association matrix, mapped it, and then let the JSON
    # writer truncate it in place -- exiting 0 with the artifact the user
    # asked for destroyed and no message at all.
    #
    # Both directions are asserted, because the obvious fix is wrong on a
    # case-sensitive filesystem, where r.json and R.JSON really are two
    # files and the invocation is legitimate. CI runs on Linux, so a
    # blanket case fold would pass here and refuse a working command there.
    first = tmp_path / "r.json"
    second = tmp_path / "R.JSON"
    result = _run("consensus", _matrix(tmp_path), "--member",
                  "butina;threshold=1.0", "--quiet", "--output", str(first),
                  "--mmap", str(second))
    if _case_folding(tmp_path):
        assert result.exit_code == 2
        assert "different files" in result.output
        assert not first.exists()
    else:
        assert result.exit_code == 0, result.output
        assert json.loads(first.read_text(encoding="utf-8"))["command"] \
            == "consensus"
        assert second.stat().st_size == 40 * 39 // 2 * 8


def test_two_absent_destinations_with_distinct_names_are_accepted(tmp_path):
    # The converse the probe must not break: making both paths exist to ask
    # the filesystem must not make two genuinely different names collide.
    out = tmp_path / "r.json"
    mapped = tmp_path / "s.bin"
    result = _run("consensus", _matrix(tmp_path), "--member",
                  "butina;threshold=1.0", "--quiet", "--output", str(out),
                  "--mmap", str(mapped))
    assert result.exit_code == 0, result.output
    assert mapped.stat().st_size == 40 * 39 // 2 * 8
    assert json.loads(out.read_text(encoding="utf-8"))["command"] == "consensus"


def test_a_refused_destination_leaves_no_placeholder(tmp_path):
    # The probe creates each absent destination so the filesystem can
    # answer. Nothing has been written at that point, so a placeholder left
    # behind is an empty results file for a run that never happened --
    # which reads as a run that succeeded and found nothing.
    out = tmp_path / "r.txt"
    result = _run("cluster", _matrix(tmp_path), "--algorithm", "butina",
                  "--set", "threshold=1.0", "--output", str(out))
    assert result.exit_code == 2
    assert not out.exists()
    assert sorted(p.name for p in tmp_path.iterdir()) == ["m.npz"]


def test_the_probe_does_not_disturb_an_existing_destination(tmp_path):
    # O_EXCL leaves an existing file alone, so a destination the user means
    # to overwrite still holds its old bytes until the writer replaces them
    # -- and a refusal for some other reason leaves it untouched.
    out = tmp_path / "r.txt"
    out.write_text("keep me", encoding="utf-8")
    result = _run("cluster", _matrix(tmp_path), "--algorithm", "butina",
                  "--set", "threshold=1.0", "--output", str(out))
    assert result.exit_code == 2
    assert out.read_text(encoding="utf-8") == "keep me"


def test_a_dangling_symlink_into_a_missing_directory_is_refused(tmp_path,
                                                                monkeypatch):
    # The link's own directory exists, so checking dirname(abspath(path))
    # passed it and the open failed only at write time. The write goes
    # through the resolved target, so that is what has to be checked.
    def refuse(*args, **kwargs):
        raise AssertionError("work started before the path check")

    monkeypatch.setattr("oecluster._cli_main._cli_input.load", refuse)
    link = tmp_path / "out.json"
    link.symlink_to(tmp_path / "missing" / "x.json")
    result = _run("consensus", _matrix(tmp_path), "--member",
                  "butina;threshold=1.0", "--output", str(link))
    assert result.exit_code == 2
    assert "directory" in result.output


def test_a_consensus_member_honours_allow_nonmetric(tmp_path, monkeypatch):
    # The only routing hole left: every other nonmetric case goes through
    # `cluster`, so dropping the flag from the member path left the suite
    # green while each member spec silently lost it.
    used = _specs_used(monkeypatch)
    path = _nonmetric(tmp_path)
    members = ["--member", "butina;threshold=0.15",
               "--member", "agglomerative;n_clusters=2"]
    refused = _run("consensus", path, *members, "--quiet")
    assert refused.exit_code != 0
    used.clear()
    allowed = _run("consensus", path, *members, "--allow-nonmetric", "--quiet")
    assert allowed.exit_code == 0, allowed.output
    assert [options["allow_nonmetric"] for options in used] == [True, True]


def test_two_dangling_links_onto_one_target_are_refused(tmp_path):
    # The fourth route into the same data loss, and the reason the checks
    # resolve first. O_CREAT | O_EXCL fails on a link whatever it points
    # at, so reserving the typed paths reserved nothing: there were no
    # inodes to compare, the two realpath strings differed by case, and
    # the run wrote the JSON document over the mapped matrix.
    #
    # Both directions again: on a case-sensitive filesystem the two
    # targets really are two files and the invocation is legitimate.
    first = tmp_path / "a.json"
    first.symlink_to(tmp_path / "result.json")
    second = tmp_path / "b.json"
    second.symlink_to(tmp_path / "RESULT.JSON")
    result = _run("consensus", _matrix(tmp_path), "--member",
                  "butina;threshold=1.0", "--quiet", "--output", str(first),
                  "--mmap", str(second))
    if _case_folding(tmp_path):
        assert result.exit_code == 2
        assert "different files" in result.output
        # Neither the target nor a reservation placeholder is left behind,
        # and the links the user made are still links.
        assert not (tmp_path / "result.json").exists()
        assert first.is_symlink()
        assert second.is_symlink()
    else:
        assert result.exit_code == 0, result.output
        assert json.loads((tmp_path / "result.json").read_text(
            encoding="utf-8"))["command"] == "consensus"
        assert (tmp_path / "RESULT.JSON").stat().st_size == 40 * 39 // 2 * 8


def test_a_dangling_link_is_written_through_to_its_target(tmp_path):
    # The converse the resolution must not break: distinct targets are
    # accepted, and the write lands on the target rather than replacing
    # the link. This also proves the reservation released the target --
    # a placeholder left in place would be overwritten here and hide it,
    # but a target deleted along with the link would not reappear.
    out = tmp_path / "out.json"
    out.symlink_to(tmp_path / "written.json")
    mapped = tmp_path / "map.bin"
    mapped.symlink_to(tmp_path / "written.bin")
    result = _run("consensus", _matrix(tmp_path), "--member",
                  "butina;threshold=1.0", "--quiet", "--output", str(out),
                  "--mmap", str(mapped))
    assert result.exit_code == 0, result.output
    assert out.is_symlink()
    assert mapped.is_symlink()
    assert json.loads((tmp_path / "written.json").read_text(
        encoding="utf-8"))["command"] == "consensus"
    assert (tmp_path / "written.bin").stat().st_size == 40 * 39 // 2 * 8


def test_a_refused_link_destination_leaves_no_target_behind(tmp_path):
    # The reservation now creates the link's target, which is a different
    # path from the one the user typed. Cleaning up the typed path would
    # delete the link and leave the placeholder standing -- the exact
    # inverse of what the cleanup is for.
    link = tmp_path / "out.txt"
    link.symlink_to(tmp_path / "written.txt")
    result = _run("cluster", _matrix(tmp_path), "--algorithm", "butina",
                  "--set", "threshold=1.0", "--output", str(link))
    assert result.exit_code == 2
    assert ".csv or .json" in result.output
    assert link.is_symlink()
    assert not (tmp_path / "written.txt").exists()
