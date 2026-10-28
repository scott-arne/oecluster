"""Validation scans stay within a fixed memory bound.

Bounds are measured with tracemalloc, which numpy reports its array
allocations to. A bound test for a path that used to allocate in proportion
to its input is paired with a control that runs the replaced code and shows
it exceeds the same bound; a bound the old code also met would prove nothing.
The C- and F-order ``_all_finite`` cases have no control: they guard layouts
the replaced scan already met.
"""

import re
import tracemalloc

import numpy as np
import oecluster
import pytest
from oecluster import (
    SymmetricDistanceMatrix,
    _gate,
    _require_kernel_distances,
    _stability,
    _upper_triangle,
)

_BOUND = 256 * 1024
_CHUNK = 1024


def _overhead(call):
    """Peak bytes allocated during ``call()``, above what was traced before."""
    tracemalloc.start()
    try:
        before = tracemalloc.get_traced_memory()[0]
        tracemalloc.reset_peak()
        call()
        return tracemalloc.get_traced_memory()[1] - before
    finally:
        tracemalloc.stop()


def _large(name, size=2_097_152, seed=0):
    rng = np.random.default_rng(seed)
    if name == "c_order":
        return rng.random(size)
    if name == "f_order":
        return np.asfortranarray(rng.random((2048, size // 2048)))
    return rng.random(2 * size)[::2]


def _small():
    base = np.arange(60, dtype=np.float64)
    return {
        "one_d": base.copy(),
        "c_order_2d": base.reshape(6, 10).copy(),
        "f_order_2d": np.asfortranarray(base.reshape(6, 10)),
        "transposed": base.reshape(10, 6).T,
        "strided": base[::3],
        "empty": np.empty(0),
    }


@pytest.mark.parametrize("name", sorted(_small()))
def test_chunks_cover_the_array_in_c_order(name, monkeypatch):
    monkeypatch.setattr(_gate, "_FINITE_SCAN_CHUNK", 7)
    array = _small()[name]
    parts, offsets = [], []
    for offset, chunk in _gate._chunks(array):
        assert chunk.ndim == 1
        assert 0 < chunk.size <= 7
        offsets.append(offset)
        parts.append(chunk.copy())  # a chunk is valid until the next one
    joined = np.concatenate(parts) if parts else np.empty(0)
    np.testing.assert_array_equal(joined, array.ravel(order="C"))
    assert offsets == [sum(part.size for part in parts[:k])
                       for k in range(len(parts))]


def _with_nan(name, position):
    """The named small layout, its layout kept, with a NaN at C-order flat
    index ``position`` (negative counts from the end)."""
    base = np.arange(60, dtype=np.float64)
    if name == "transposed":
        view = base.reshape(10, 6).copy().T
    elif name == "strided":
        view = base.copy()[::3]
    elif name == "f_order_2d":
        view = np.asfortranarray(base.reshape(6, 10))
    elif name == "c_order_2d":
        view = base.reshape(6, 10).copy()
    else:
        view = base.copy()
    view[np.unravel_index(position % view.size, view.shape)] = np.nan
    return view


def test_all_finite_is_true_for_an_empty_array():
    assert _gate._all_finite(np.empty(0)) is True


@pytest.mark.parametrize("dtype", [object, "U3"], ids=["object", "str"])
def test_all_finite_refuses_an_empty_array_of_a_refused_dtype(dtype):
    """``np.isfinite`` refuses by dtype, so an empty array is no exception."""
    array = np.empty(0, dtype=dtype)
    with pytest.raises(TypeError):
        np.isfinite(array).all()
    with pytest.raises(TypeError):
        _gate._all_finite(array)
    assert _gate._all_finite(np.empty(0)) is True


def test_all_finite_scans_an_f_ordered_matrix_as_its_transpose(monkeypatch):
    """Finiteness does not depend on order, so an F-ordered matrix is handed
    to ``_chunks`` as its C-contiguous transpose, which takes the flat path
    rather than the slower C-order ``nditer`` walk."""
    seen = []
    real = _gate._chunks

    def spy(values):
        seen.append(values.flags.c_contiguous)
        return real(values)

    monkeypatch.setattr(_gate, "_chunks", spy)
    array = _small()["f_order_2d"]
    assert not array.flags.c_contiguous
    assert _gate._all_finite(array) is True
    assert seen == [True]


@pytest.mark.parametrize("name", sorted(set(_small()) - {"empty"}))
@pytest.mark.parametrize("position", [0, 6, 7, -1])
def test_all_finite_finds_a_nan_wherever_it_sits(name, position, monkeypatch):
    """Positions 0, 6, 7 and the last sit at a chunk start, a chunk end, the
    next chunk's start and the last partial chunk of a C-order scan with a
    chunk of 7. The F-ordered layouts are scanned as their C-contiguous
    transposes, whose chunk boundaries the C-order cases cover; for them this
    shows only that the NaN is found."""
    monkeypatch.setattr(_gate, "_FINITE_SCAN_CHUNK", 7)
    assert _gate._all_finite(_small()[name]) is True
    assert _gate._all_finite(_with_nan(name, position)) is False


def test_all_finite_refuses_an_object_array_as_isfinite_does(monkeypatch):
    monkeypatch.setattr(_gate, "_FINITE_SCAN_CHUNK", 7)
    array = np.array([1.0, 2.0, 3.0, 4.0], dtype=object)[::2]
    with pytest.raises(TypeError):
        np.isfinite(array).all()
    with pytest.raises(TypeError):
        _gate._all_finite(array)


@pytest.mark.parametrize("name", ["c_order", "f_order", "strided"])
def test_all_finite_stays_within_the_bound(name, monkeypatch):
    monkeypatch.setattr(_gate, "_FINITE_SCAN_CHUNK", _CHUNK)
    array = _large(name)
    assert _overhead(lambda: _gate._all_finite(array)) <= _BOUND


def test_the_release_c_scan_copies_a_strided_input():
    """Control for the strided bound: ``ravel(order="K")`` copied it whole."""
    array = _large("strided")

    def old():
        flat = array.ravel(order="K")
        return all(np.isfinite(flat[start:start + _CHUNK]).all()
                   for start in range(0, flat.size, _CHUNK))

    assert _overhead(old) > _BOUND


def _spy_on_all_finite(monkeypatch):
    calls = []
    real = _gate._all_finite

    def spy(values):
        calls.append(values.shape)
        return real(values)

    monkeypatch.setattr(_gate, "_all_finite", spy)
    return calls


def test_from_condensed_scans_through_all_finite(monkeypatch):
    calls = _spy_on_all_finite(monkeypatch)
    oecluster.SymmetricDistanceMatrix.from_condensed(np.array([0.1, 0.2, 0.3]))
    assert calls


def test_stability_scans_a_dense_subset_through_all_finite(monkeypatch):
    dense = oecluster.SymmetricDistanceMatrix.from_condensed(
        np.array([0.1, 0.2, 0.3]))
    calls = _spy_on_all_finite(monkeypatch)
    facts = _stability._data_facts(
        dense, {"data_integrity": "nan_present", "metric_probe": "not_run"})
    assert calls
    assert facts["data_integrity"] == "complete"


_N = 1449  # 1449 * 1449 = 2,099,601 elements, above 2,097,152


def _metric_square(n=_N, order="C", seed=1):
    rng = np.random.default_rng(seed)
    positions = rng.random(n)
    square = np.abs(positions[:, None] - positions[None, :])
    return np.asfortranarray(square) if order == "F" else square


def _old_checks(square, condensed):
    """The release-C value checks, for the controls."""
    np.all(np.isfinite(square))
    np.allclose(square, square.T, rtol=1e-5, atol=1e-8)
    float(np.abs(condensed).max())
    np.any(condensed < 0.0)


@pytest.mark.parametrize("order", ["C", "F"])
def test_the_square_value_checks_stay_within_the_bound(order, monkeypatch):
    monkeypatch.setattr(_gate, "_FINITE_SCAN_CHUNK", _CHUNK)
    square = _metric_square(order=order)
    condensed = _upper_triangle(square)
    assert _overhead(lambda: _old_checks(square, condensed)) > _BOUND
    assert _overhead(lambda: SymmetricDistanceMatrix._check_condensed_values(
        square, square, condensed)) <= _BOUND


def test_the_condensed_value_checks_stay_within_the_bound(monkeypatch):
    monkeypatch.setattr(_gate, "_FINITE_SCAN_CHUNK", _CHUNK)
    condensed = np.random.default_rng(2).random(2_097_152)

    def old():
        np.all(np.isfinite(condensed))
        np.any(condensed < 0.0)

    assert _overhead(old) > _BOUND
    assert _overhead(lambda: SymmetricDistanceMatrix._check_condensed_values(
        condensed, None, condensed)) <= _BOUND


@pytest.mark.parametrize("order", ["C", "F"])
def test_the_upper_triangle_is_copied_without_index_arrays(order):
    square = _metric_square(order=order)
    output = _N * (_N - 1) // 2 * 8
    old = _overhead(lambda: square[np.triu_indices(_N, k=1)])
    new = _overhead(lambda: _upper_triangle(square))
    assert old > output + _BOUND
    assert new <= output + _BOUND
    np.testing.assert_array_equal(_upper_triangle(square),
                                  square[np.triu_indices(_N, k=1)])


def test_from_condensed_refuses_both_mirrors_of_a_one_sided_tolerance(
        monkeypatch):
    a, b = 1.0, 1.00001001005006
    assert np.isclose(a, b, rtol=1e-5, atol=1e-8)
    assert not np.isclose(b, a, rtol=1e-5, atol=1e-8)
    upper_a = np.array([[0.0, a], [b, 0.0]])
    for chunk in (1, 1024):
        monkeypatch.setattr(_gate, "_FINITE_SCAN_CHUNK", chunk)
        for square in (upper_a, upper_a.T.copy()):
            with pytest.raises(ValueError, match="must be symmetric"):
                SymmetricDistanceMatrix.from_condensed(square)


@pytest.mark.parametrize("chunk", [1, 3, 1024])
def test_from_condensed_reports_the_same_first_failure(chunk, monkeypatch):
    monkeypatch.setattr(_gate, "_FINITE_SCAN_CHUNK", chunk)
    asymmetric_and_negative = np.array(
        [[0.0, -0.5, 0.2], [-0.7, 0.0, 0.3], [0.2, 0.3, 0.0]])
    diagonal_and_negative = np.array(
        [[1.0, -0.5, -0.5], [-0.5, 1.0, -0.5], [-0.5, -0.5, 1.0]])
    negative_last = np.array([0.1, 0.2, 0.3, 0.4, 0.5, -0.6])
    with pytest.raises(ValueError, match="must be symmetric"):
        SymmetricDistanceMatrix.from_condensed(asymmetric_and_negative)
    with pytest.raises(ValueError, match="must have a zero diagonal"):
        SymmetricDistanceMatrix.from_condensed(diagonal_and_negative)
    with pytest.raises(ValueError, match="contains negative values"):
        SymmetricDistanceMatrix.from_condensed(negative_last)


def test_from_condensed_reads_a_strided_square_as_before():
    square = _metric_square(n=40)
    wide = np.zeros((40, 80))
    wide[:, ::2] = square
    strided = wide[:, ::2]
    expected = SymmetricDistanceMatrix.from_condensed(square.copy())
    got = SymmetricDistanceMatrix.from_condensed(strided)
    np.testing.assert_array_equal(got.condensed, expected.condensed)


_ENTRIES = 1_048_576


def _valid_entries(n=2048, count=_ENTRIES, seed=3):
    rng = np.random.default_rng(seed)
    rows = rng.integers(0, n - 1, size=count).astype(np.int64)
    cols = rows + 1
    values = rng.random(count) * 0.5
    return rows, cols, values, n


def test_sparse_entry_validation_stays_within_the_bound(monkeypatch):
    monkeypatch.setattr(_gate, "_FINITE_SCAN_CHUNK", _CHUNK)
    rows, cols, values, n = _valid_entries()

    def old():
        np.flatnonzero(values > 0.9)
        np.flatnonzero((rows < 0) | (cols < 0) | (rows >= n) | (cols >= n)
                       | (rows == cols))

    assert _overhead(old) > _BOUND
    assert _overhead(lambda: SymmetricDistanceMatrix._validate_sparse_entries(
        rows, cols, values, n, 0.9, "test")) <= _BOUND


@pytest.mark.parametrize("chunk", [1, 3, 1024])
def test_sparse_entry_validation_reports_the_same_first_offender(
        chunk, monkeypatch):
    monkeypatch.setattr(_gate, "_FINITE_SCAN_CHUNK", chunk)
    # Entries 3 and 6 are bad pairs; entries 4 and 7 are over the cutoff. The
    # cutoff pass runs over every entry first, so entry 4 is reported ahead of
    # the earlier bad pair. Each pass has a second offender: at chunks 1 and 3
    # it sits in a later chunk, and at 1024 in the same chunk as the first, so
    # neither may be the one reported.
    rows = np.array([0, 1, 2, 3, 0, 1, 7, 2, 3], dtype=np.int64)
    cols = np.array([1, 2, 3, 3, 2, 3, 0, 4, 5], dtype=np.int64)
    values = np.full(9, 0.1)
    values[4] = 0.95
    values[7] = 0.97
    with pytest.raises(ValueError, match="sparse entry 4 has value 0.95"):
        SymmetricDistanceMatrix._validate_sparse_entries(
            rows, cols, values, 6, 0.9, "test")
    values[[4, 7]] = 0.5
    with pytest.raises(ValueError, match=r"sparse entry 3 is \(3, 3\)"):
        SymmetricDistanceMatrix._validate_sparse_entries(
            rows, cols, values, 6, 0.9, "test")


@pytest.mark.parametrize("dtype", ["<U1", "S1", "M8[s]"])
def test_sparse_entry_validation_refuses_an_empty_non_numeric_value_array(
        dtype):
    """The cutoff comparison refuses by dtype, so an empty array is no
    exception, though the chunked scan yields no chunk for it."""
    empty_indices = np.empty(0, dtype=np.int64)
    values = np.array([], dtype=dtype)
    with pytest.raises(TypeError):
        np.greater(values, 0.9)
    with pytest.raises(TypeError):
        SymmetricDistanceMatrix._validate_sparse_entries(
            empty_indices, empty_indices, values, 5, 0.9, "test")
    SymmetricDistanceMatrix._validate_sparse_entries(
        empty_indices, empty_indices, np.array([], dtype=np.float64), 5, 0.9,
        "test")


def _kernel_matrix(n=2049, seed=4):
    values = np.random.default_rng(seed).random(n * (n - 1) // 2) * 0.9
    return SymmetricDistanceMatrix.from_condensed(values, check=False)


@pytest.mark.parametrize("kernel", ["complement", "laplacian"])
def test_the_kernel_range_check_stays_within_the_bound(kernel, monkeypatch):
    monkeypatch.setattr(_gate, "_FINITE_SCAN_CHUNK", _CHUNK)
    matrix = _kernel_matrix()
    condensed = matrix.condensed

    def old():
        mask = condensed < 0.0
        if kernel == "complement":
            mask |= condensed > 1.0
        mask.any()

    assert _overhead(old) > _BOUND
    assert _overhead(lambda: _require_kernel_distances(
        matrix, kernel, "vendi_score")) <= _BOUND


@pytest.mark.parametrize("chunk", [1, 3, 1024])
@pytest.mark.parametrize("position, value, kernel", [
    (0, -0.1, "laplacian"), (2, 1.5, "complement"), (5, -0.2, "complement"),
])
def test_the_kernel_range_check_names_the_same_first_pair(
        chunk, position, value, kernel, monkeypatch):
    monkeypatch.setattr(_gate, "_FINITE_SCAN_CHUNK", chunk)
    values = np.full(6, 0.5)
    values[position] = value
    if position != 5:
        values[5] = -0.3  # a later offender must not be the one reported
    matrix = SymmetricDistanceMatrix.from_condensed(values, check=False)
    pairs = {0: "d(0, 1)", 2: "d(0, 3)", 5: "d(2, 3)"}
    with pytest.raises(ValueError, match=re.escape(f"{pairs[position]} = {value!r}")):
        _require_kernel_distances(matrix, kernel, "vendi_score")


def _spy_on_chunks(monkeypatch):
    calls = []
    real = _gate._chunks

    def spy(values):
        calls.append(values.shape)
        return real(values)

    monkeypatch.setattr(_gate, "_chunks", spy)
    return calls


def test_sparse_entry_validation_scans_through_chunks(monkeypatch):
    rows, cols, values, n = _valid_entries(count=10)
    calls = _spy_on_chunks(monkeypatch)
    SymmetricDistanceMatrix._validate_sparse_entries(
        rows, cols, values, n, 0.9, "test")
    assert calls


def test_the_kernel_range_check_scans_through_chunks(monkeypatch):
    matrix = _kernel_matrix(n=10)
    calls = _spy_on_chunks(monkeypatch)
    _require_kernel_distances(matrix, "complement", "vendi_score")
    assert calls
