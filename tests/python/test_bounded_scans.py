"""Validation scans stay within a fixed memory bound.

Bounds are measured with tracemalloc, which numpy reports its array
allocations to. A bound test for a path that used to allocate in proportion
to its input is paired with a control that runs the replaced code and shows
it exceeds the same bound; a bound the old code also met would prove nothing.
The C- and F-order ``_all_finite`` cases have no control: they guard layouts
the replaced scan already met.
"""

import tracemalloc

import numpy as np
import oecluster
import pytest
from oecluster import _gate, _stability

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
