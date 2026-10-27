"""Counts in archived matrices are non-negative integer scalars."""

import errno
import json

import numpy as np
import oecluster
import pytest

_PREFIX = "must be a non-negative integer scalar"


def _rewrite(path, **changes):
    """Rewrite an archive with keys replaced, or removed where the value is None."""
    with np.load(str(path), allow_pickle=False) as data:
        arrays = {key: data[key] for key in data.files}
    for key, value in changes.items():
        if value is None:
            arrays.pop(key, None)
        else:
            arrays[key] = value
    np.savez_compressed(str(path), **arrays)


def _symmetric(tmp_path, n=3):
    path = tmp_path / "symmetric.npz"
    values = np.arange(1, n * (n - 1) // 2 + 1, dtype=np.float64) / 10.0
    oecluster.SymmetricDistanceMatrix.from_condensed(values).to_file(str(path))
    return path


def _cross(tmp_path):
    path = tmp_path / "cross.npz"
    matrix = np.array([[0.1, 0.2, 0.3], [0.4, 0.5, 0.6]])
    oecluster.CrossDistanceMatrix(matrix, "test").to_file(str(path))
    return path


def _bad(n):
    return {
        "float": np.array(float(n)),
        "bool": np.array(True),
        "string": np.array(str(n)),
        "object": np.array(n, dtype=object),
        "array": np.array([n]),
        "negative": np.array(-1),
    }


@pytest.mark.parametrize("case", sorted(_bad(3)))
def test_a_malformed_num_samples_is_refused(tmp_path, case):
    path = _symmetric(tmp_path)
    _rewrite(path, num_samples=_bad(3)[case])
    with pytest.raises(ValueError,
                       match=f"Malformed symmetric matrix: num_samples {_PREFIX}"):
        oecluster.SymmetricDistanceMatrix.from_file(str(path))


@pytest.mark.parametrize("case", sorted(_bad(3)))
def test_a_malformed_legacy_num_items_is_refused(tmp_path, case):
    path = _symmetric(tmp_path)
    _rewrite(path, num_samples=None, num_items=_bad(3)[case])
    with pytest.raises(ValueError,
                       match=f"Malformed symmetric matrix: num_items {_PREFIX}"):
        oecluster.SymmetricDistanceMatrix.from_file(str(path))


@pytest.mark.parametrize("key, n", [("n_a", 2), ("n_b", 3)])
@pytest.mark.parametrize("case", sorted(_bad(2)))
def test_a_malformed_cross_count_is_refused(tmp_path, key, n, case):
    path = _cross(tmp_path)
    _rewrite(path, **{key: _bad(n)[case]})
    with pytest.raises(ValueError,
                       match=f"Malformed cross matrix: {key} {_PREFIX}"):
        oecluster.CrossDistanceMatrix.from_file(str(path))


def test_a_negative_count_that_matches_the_length_is_refused(tmp_path):
    # -1 * -2 // 2 == 1, so one condensed value used to pass the length check
    # and reach DenseStorage(-1), which raised OverflowError.
    path = _symmetric(tmp_path, n=2)
    _rewrite(path, num_samples=np.array(-1))
    with pytest.raises(ValueError,
                       match=f"Malformed symmetric matrix: num_samples {_PREFIX}"):
        oecluster.SymmetricDistanceMatrix.from_file(str(path))


_DTYPES = [np.int32, np.int64, np.uint16, np.uint64]


@pytest.mark.parametrize("dtype", _DTYPES)
@pytest.mark.parametrize("key", ["num_samples", "num_items"])
def test_any_integer_scalar_symmetric_count_loads(tmp_path, key, dtype):
    path = _symmetric(tmp_path)
    changes = {key: np.array(3, dtype=dtype)}
    if key == "num_items":
        changes["num_samples"] = None  # a legacy archive has no num_samples
    _rewrite(path, **changes)
    assert oecluster.SymmetricDistanceMatrix.from_file(str(path)).num_samples == 3


@pytest.mark.parametrize("dtype", _DTYPES)
@pytest.mark.parametrize("key, n", [("n_a", 2), ("n_b", 3)])
def test_any_integer_scalar_cross_count_loads(tmp_path, key, n, dtype):
    path = _cross(tmp_path)
    _rewrite(path, **{key: np.array(n, dtype=dtype)})
    assert oecluster.CrossDistanceMatrix.from_file(str(path)).shape == (2, 3)


@pytest.mark.parametrize("key", ["num_samples", "num_items"])
def test_a_zero_item_symmetric_archive_loads(tmp_path, key):
    path = tmp_path / "empty.npz"
    oecluster.SymmetricDistanceMatrix.from_condensed(
        np.empty((0, 0))).to_file(str(path))
    changes = {key: np.array(0)}
    if key == "num_items":
        changes["num_samples"] = None
    _rewrite(path, **changes)
    assert oecluster.SymmetricDistanceMatrix.from_file(str(path)).num_samples == 0


@pytest.mark.parametrize("shape", [(0, 3), (2, 0)])
def test_a_zero_count_cross_archive_loads(tmp_path, shape):
    path = tmp_path / "empty_cross.npz"
    np.savez_compressed(
        str(path), matrix_kind=np.array("cross"), matrix=np.empty(shape),
        n_a=np.array(shape[0]), n_b=np.array(shape[1]),
        comparison_name=np.array("test"))
    assert oecluster.CrossDistanceMatrix.from_file(str(path)).shape == shape


def test_num_samples_wins_and_num_items_is_not_read(tmp_path):
    path = _symmetric(tmp_path)
    # An object array cannot be read under allow_pickle=False, so this loads
    # only if num_items is never read.
    _rewrite(path, num_items=np.array(99, dtype=object))
    assert oecluster.SymmetricDistanceMatrix.from_file(str(path)).num_samples == 3


def test_a_damaged_count_member_is_a_malformed_file(tmp_path):
    import zipfile

    path = _symmetric(tmp_path)
    with zipfile.ZipFile(str(path)) as archive:
        info = archive.getinfo("num_samples.npy")
    raw = bytearray(path.read_bytes())
    offset = info.header_offset
    name_length = int.from_bytes(raw[offset + 26:offset + 28], "little")
    extra_length = int.from_bytes(raw[offset + 28:offset + 30], "little")
    data_start = offset + 30 + name_length + extra_length
    raw[data_start + info.compress_size - 1] ^= 0xFF
    path.write_bytes(bytes(raw))
    with pytest.raises(ValueError,
                       match=f"Malformed symmetric matrix: num_samples {_PREFIX}") as error:
        oecluster.SymmetricDistanceMatrix.from_file(str(path))
    assert error.value.__cause__ is not None


def test_an_archive_with_no_count_is_refused(tmp_path):
    path = _symmetric(tmp_path)
    _rewrite(path, num_samples=None, num_items=None)
    with pytest.raises(ValueError,
                       match="missing required key 'num_samples'"):
        oecluster.SymmetricDistanceMatrix.from_file(str(path))


def test_an_oversized_sparse_count_is_a_value_error_from_the_native_refusal(
        tmp_path):
    path = tmp_path / "sparse.npz"
    np.savez_compressed(
        str(path),
        comparison_name=np.array("test"),
        params_json=np.array(json.dumps({})),
        labels=np.array([]),
        num_samples=np.array(2**32 + 1),
        storage_kind=np.array("sparse"),
        sparse_cutoff=np.array(0.5, dtype=np.float64),
        sparse_i=np.array([], dtype=np.int64),
        sparse_j=np.array([], dtype=np.int64),
        sparse_v=np.array([], dtype=np.float64),
    )
    with pytest.raises(ValueError) as info:
        oecluster.SymmetricDistanceMatrix.from_file(str(path))
    assert str(info.value) == (
        "Malformed symmetric matrix: SparseStorage refused num_samples "
        "4294967297")
    assert isinstance(info.value.__cause__, RuntimeError)


def test_a_memory_error_reading_a_count_is_not_a_malformed_file(
        tmp_path, monkeypatch):
    path = _symmetric(tmp_path)
    original = np.lib.npyio.NpzFile.__getitem__

    def getitem(self, key):
        if key == "num_samples":
            raise MemoryError("simulated")
        return original(self, key)

    monkeypatch.setattr(np.lib.npyio.NpzFile, "__getitem__", getitem)
    with pytest.raises(MemoryError, match="simulated"):
        oecluster.SymmetricDistanceMatrix.from_file(str(path))


def test_an_os_error_reading_a_count_is_not_a_malformed_file(
        tmp_path, monkeypatch):
    path = _symmetric(tmp_path)
    original = np.lib.npyio.NpzFile.__getitem__

    def getitem(self, key):
        if key == "num_samples":
            raise OSError(errno.EIO, "I/O error")
        return original(self, key)

    monkeypatch.setattr(np.lib.npyio.NpzFile, "__getitem__", getitem)
    with pytest.raises(OSError, match="I/O error"):
        oecluster.SymmetricDistanceMatrix.from_file(str(path))


def test_a_count_member_without_the_npy_magic_is_refused(tmp_path):
    import zipfile

    path = _symmetric(tmp_path)
    with zipfile.ZipFile(str(path)) as archive:
        members = {name: archive.read(name) for name in archive.namelist()}
    members["num_samples.npy"] = b"3"  # numpy hands back raw bytes for this
    with zipfile.ZipFile(str(path), "w") as archive:
        for name, content in members.items():
            archive.writestr(name, content)
    with pytest.raises(ValueError,
                       match=f"Malformed symmetric matrix: num_samples {_PREFIX}"):
        oecluster.SymmetricDistanceMatrix.from_file(str(path))
