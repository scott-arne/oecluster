"""Native sparse accessors, and the callers that read entries through them."""

import itertools
import math
import subprocess
import sys
import tracemalloc

import numpy as np
import oecluster
import pytest
from oecluster import _cli_input, _gate, _sparse_columns, _sparse_dense, _stability


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
    assert all(array.flags.writeable for array in _sparse_columns(storage))
    assert _sparse_dense(storage).flags.writeable


def test_entry_columns_carry_an_index_above_32_bits():
    n = 2**32
    storage = _storage(n, 0.9, [(0, n - 1, 0.5)])
    rows, cols, values = _sparse_columns(storage)
    assert (rows.tolist(), cols.tolist(), values.tolist()) == (
        [0], [n - 1], [0.5])
    assert storage._entries() == [(0, n - 1, 0.5)]


def test_densify_refuses_a_pair_count_no_array_can_hold():
    storage = oecluster.SparseStorage(2**32, 0.9)
    storage.Finalize()
    with pytest.raises(Exception) as expected:
        np.zeros(storage.NumPairs(), dtype=np.float64)
    with pytest.raises(type(expected.value)):
        _sparse_dense(storage)


@pytest.mark.parametrize("chunk", [2, 1024])
def test_densify_keeps_the_last_entry_for_a_pair(chunk, monkeypatch):
    # With a chunk of 2 the duplicates of each pair fall in different
    # chunks; with 1024 they share one.
    monkeypatch.setattr(_gate, "_FINITE_SCAN_CHUNK", chunk)
    storage = _storage(4, 0.9, [(0, 1, math.nan), (2, 3, 0.2), (0, 1, 0.3),
                                (1, 3, 0.4), (3, 2, 0.5)])
    dense = _sparse_dense(storage)
    np.testing.assert_array_equal(dense, _old_densify(storage))
    assert dense[0] == 0.3  # pair (0, 1)
    assert dense[5] == 0.5  # pair (2, 3)


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


def test_densifying_a_sparse_matrix_allocates_its_output_columns_and_a_chunk(
        monkeypatch):
    # The default chunk holds more entries than this test stores, so the
    # whole scatter would be one chunk and the bound could not tell a
    # chunked scatter from an unchunked one. With 1024 entries, as in
    # test_bounded_scans.py, a chunk's scratch of about 65 bytes an entry
    # stays inside the 256 KiB allowance.
    monkeypatch.setattr(_gate, "_FINITE_SCAN_CHUNK", 1024)
    n = 600
    rng = np.random.default_rng(0)
    pairs = set()
    while len(pairs) < 131_072:
        i, j = sorted(rng.integers(0, n, size=2).tolist())
        if i != j:
            pairs.add((i, j))
    storage = _storage(n, 0.9, [(i, j, 0.5) for i, j in sorted(pairs)])
    output = storage.NumPairs() * 8
    # The output, the int64 row and column and float64 value columns
    # _entry_columns() returns, and the scratch of one chunk.
    columns = len(pairs) * 24
    bound = output + columns + 256 * 1024

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


_RESIDENCY_SCRIPT = """
import resource
import sys

import numpy as np
import oecluster


def peak():
    value = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return value if sys.platform == "darwin" else value * 1024


storage = oecluster.SparseStorage({n}, 0.9)
for i, j in {pairs}:
    storage.Set(i, j, 0.5)
storage.Finalize()
matrix = oecluster.SymmetricDistanceMatrix(storage, "test")
m0 = peak()
condensed = matrix.condensed
m1 = peak()
control = np.empty(storage.NumPairs(), dtype=np.float64)
control.fill(0.0)
m2 = peak()
print(m0, m1, m2)
"""


def test_densifying_a_sparse_matrix_leaves_untouched_pages_unbacked():
    # tracemalloc counts allocated bytes, so a densify that writes every page
    # of its output passes the bound above; resident memory tells them apart.
    pytest.importorskip("resource")
    n = 10_000
    pairs = [(k * 97, k * 97 + 1) for k in range(100)]
    offsets = [n * i + j - ((i + 2) * (i + 1)) // 2 for i, j in pairs]
    # 64 KiB apart, so each entry lands on its own page at any page size up
    # to 64 KiB.
    assert min(b - a for a, b in itertools.pairwise(offsets)) * 8 >= 65_536
    script = _RESIDENCY_SCRIPT.format(n=n, pairs=pairs)
    result = subprocess.run([sys.executable, "-X", "faulthandler", "-c", script],
                            capture_output=True, text=True, timeout=120,
                            check=False)
    report = (result.returncode, result.stdout, result.stderr[-2000:])
    assert result.returncode == 0, report
    m0, m1, m2 = (int(value) for value in result.stdout.split())
    bound = 128 * 2**20
    # ru_maxrss is a high-water mark, so a delta can only understate growth;
    # the densify is measured first so the control's peak cannot hide it.
    # 100 touched pages are at most 6.25 MiB even with 64 KiB pages, while
    # the output is about 381 MiB.
    assert m1 - m0 <= bound, report
    # The control writes all 381 MiB. If that does not exceed the bound, this
    # platform's ru_maxrss cannot see residency, and the test must not pass.
    assert m2 - m1 > bound, report
