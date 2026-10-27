"""Native sparse accessors, and the callers that read entries through them."""

import math
import tracemalloc

import numpy as np
import oecluster
import pytest
from oecluster import _cli_input, _sparse_columns, _sparse_dense, _stability


def _storage(n, cutoff, entries):
    storage = oecluster.SparseStorage(n, cutoff)
    for i, j, value in entries:
        storage.Set(i, j, value)
    storage.Finalize()
    return storage


def _old_densify(storage):
    """The Python loop ``.condensed`` ran before release D."""
    n = storage.NumSamples()
    condensed = np.zeros(storage.NumPairs(), dtype=np.float64)
    for i, j, val in storage._entries():
        if i < j:
            idx = n * i + j - ((i + 2) * (i + 1)) // 2
        else:
            idx = n * j + i - ((j + 2) * (j + 1)) // 2
        condensed[idx] = val
    return condensed


_CASES = {
    "empty": (4, 0.9, []),
    "finite": (4, 0.9, [(0, 1, 0.1), (2, 3, 0.2), (1, 3, 0.3)]),
    "non_finite": (4, math.inf,
                   [(0, 1, math.nan), (0, 2, math.inf), (1, 2, -math.inf)]),
    "superseded": (4, 0.9, [(0, 1, math.nan), (0, 1, 0.3), (3, 2, 0.4)]),
}


@pytest.mark.parametrize("case", sorted(_CASES))
def test_entry_columns_match_the_entry_list(case):
    storage = _storage(*_CASES[case])
    rows, cols, values = _sparse_columns(storage)
    entries = storage._entries()
    assert (rows.dtype, cols.dtype, values.dtype) == (
        np.int64, np.int64, np.float64)
    assert rows.tolist() == [entry[0] for entry in entries]
    assert cols.tolist() == [entry[1] for entry in entries]
    np.testing.assert_array_equal(
        values, np.array([entry[2] for entry in entries], dtype=np.float64))


@pytest.mark.parametrize("case", sorted(_CASES))
def test_densify_matches_the_old_python_loop(case):
    storage = _storage(*_CASES[case])
    np.testing.assert_array_equal(_sparse_dense(storage),
                                  _old_densify(storage))


def test_the_accessors_return_bytearrays_and_writable_arrays():
    storage = _storage(*_CASES["finite"])
    columns = storage._entry_columns()
    assert isinstance(columns, tuple) and len(columns) == 3
    assert all(isinstance(column, bytearray) for column in columns)
    assert isinstance(storage._densify(), bytearray)
    assert all(array.flags.writeable for array in _sparse_columns(storage))
    assert _sparse_dense(storage).flags.writeable


def test_entry_columns_carry_an_index_above_32_bits():
    n = 2**32
    storage = _storage(n, 0.9, [(0, n - 1, 0.5)])
    rows, cols, values = _sparse_columns(storage)
    assert (rows.tolist(), cols.tolist(), values.tolist()) == (
        [0], [n - 1], [0.5])
    assert storage._entries() == [(0, n - 1, 0.5)]


def test_densify_refuses_a_pair_count_no_bytearray_can_hold():
    storage = oecluster.SparseStorage(2**32, 0.9)
    storage.Finalize()
    with pytest.raises(MemoryError):
        storage._densify()


def _refuse_entry_list(monkeypatch):
    def refuse(self):
        raise AssertionError("a caller built the per-entry Python list")
    monkeypatch.setattr(oecluster.SparseStorage, "_entries", refuse)


def _matrix(storage):
    return oecluster.SymmetricDistanceMatrix(storage, "test")


@pytest.mark.parametrize("case", sorted(_CASES))
def test_condensed_densifies_without_the_entry_list(case, monkeypatch):
    storage = _storage(*_CASES[case])
    expected = _old_densify(storage)
    _refuse_entry_list(monkeypatch)
    condensed = _matrix(storage).condensed
    np.testing.assert_array_equal(condensed, expected)
    if condensed.size:
        condensed[0] = 1.0  # still writable, as it was


@pytest.mark.parametrize("case", sorted(_CASES))
def test_sparse_export_keeps_every_entry_without_the_entry_list(
        case, monkeypatch):
    pytest.importorskip("scipy")
    storage = _storage(*_CASES[case])
    entries = storage._entries()
    _refuse_entry_list(monkeypatch)
    coo = _matrix(storage).sparse
    assert coo.shape == (storage.NumSamples(), storage.NumSamples())
    assert coo.row.tolist() == [entry[0] for entry in entries]
    assert coo.col.tolist() == [entry[1] for entry in entries]
    np.testing.assert_array_equal(
        coo.data, np.array([entry[2] for entry in entries], dtype=np.float64))


def test_to_file_writes_the_entries_without_the_entry_list(
        tmp_path, monkeypatch):
    storage = _storage(*_CASES["superseded"])
    entries = storage._entries()
    path = tmp_path / "sparse.npz"
    _refuse_entry_list(monkeypatch)
    _matrix(storage).to_file(str(path))
    with np.load(str(path), allow_pickle=False) as data:
        assert data["sparse_i"].dtype == np.int64
        assert data["sparse_j"].dtype == np.int64
        assert data["sparse_v"].dtype == np.float64
        assert data["sparse_i"].tolist() == [entry[0] for entry in entries]
        assert data["sparse_j"].tolist() == [entry[1] for entry in entries]
        np.testing.assert_array_equal(
            data["sparse_v"],
            np.array([entry[2] for entry in entries], dtype=np.float64))


def test_the_cli_negative_check_reads_values_without_the_entry_list(
        tmp_path, monkeypatch):
    clean = _matrix(_storage(4, 0.9, [(0, 1, 0.1)]))
    negative = _matrix(_storage(4, 0.9, [(0, 1, 0.1), (2, 3, -0.5)]))
    _refuse_entry_list(monkeypatch)
    _cli_input._refuse_negative(clean, str(tmp_path / "clean.npz"))
    with pytest.raises(_cli_input.InputError, match="negative"):
        _cli_input._refuse_negative(negative, str(tmp_path / "negative.npz"))


def test_stability_scans_a_sparse_subset_without_the_entry_list(monkeypatch):
    clean = _matrix(_storage(4, 0.9, [(0, 1, 0.1)]))
    holed = _matrix(_storage(4, 0.9, [(0, 1, math.nan)]))
    _refuse_entry_list(monkeypatch)
    source = {"data_integrity": "nan_present", "metric_probe": "not_run"}
    assert _stability._data_facts(holed, source)["data_integrity"] == "nan_present"
    # A non-descriptor subset with no NaN left is complete.
    assert _stability._data_facts(clean, source)["data_integrity"] == "complete"


def test_densifying_a_sparse_matrix_allocates_only_its_output(monkeypatch):
    n = 600
    rng = np.random.default_rng(0)
    pairs = set()
    while len(pairs) < 131_072:
        i, j = sorted(rng.integers(0, n, size=2).tolist())
        if i != j:
            pairs.add((i, j))
    storage = _storage(n, 0.9, [(i, j, 0.5) for i, j in sorted(pairs)])
    output = storage.NumPairs() * 8
    bound = output + 256 * 1024

    def overhead(call):
        tracemalloc.start()
        try:
            before = tracemalloc.get_traced_memory()[0]
            tracemalloc.reset_peak()
            call()
            return tracemalloc.get_traced_memory()[1] - before
        finally:
            tracemalloc.stop()

    matrix = _matrix(storage)
    # Control: the old loop builds a Python tuple per entry.
    assert overhead(lambda: _old_densify(storage)) > bound
    assert overhead(lambda: matrix.condensed) <= bound
