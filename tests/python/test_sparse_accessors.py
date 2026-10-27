"""Native sparse accessors, and the callers that read entries through them."""

import math

import numpy as np
import oecluster
import pytest
from oecluster import _sparse_columns, _sparse_dense


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
