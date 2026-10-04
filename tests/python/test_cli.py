"""The oecluster command line: registry, inputs, commands and output."""
import json
import os
import time
import warnings
from pathlib import Path

import numpy as np
import oecluster
import pytest
from click.testing import CliRunner
from oecluster import _cli_input, _cli_registry
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
