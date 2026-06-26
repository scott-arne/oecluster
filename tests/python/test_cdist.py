import numpy as np
import pytest


def test_cross_distance_matrix_basics():
    """CrossDistanceMatrix wraps a rectangular array with metadata."""
    from oecluster import CrossDistanceMatrix, DistanceMatrix

    mat = np.array([[0.1, 0.2, 0.3], [0.4, 0.5, 0.6]], dtype=np.float64)
    cdm = CrossDistanceMatrix(mat, "fingerprint", ["a0", "a1"], ["b0", "b1", "b2"], {})

    assert isinstance(cdm, DistanceMatrix)
    assert cdm.shape == (2, 3)
    assert cdm.labels_a == ["a0", "a1"]
    assert cdm.labels_b == ["b0", "b1", "b2"]
    assert cdm.comparison_name == "fingerprint"
    assert len(cdm) == 2
    np.testing.assert_array_equal(np.asarray(cdm), mat)
    assert repr(cdm) == "CrossDistanceMatrix(comparison='fingerprint', shape=(2, 3))"


def test_cross_distance_matrix_roundtrip(tmp_path):
    """CrossDistanceMatrix serializes and reloads via .npz."""
    from oecluster import CrossDistanceMatrix

    mat = np.array([[0.1, 0.2], [0.3, 0.4], [0.5, 0.6]], dtype=np.float64)
    cdm = CrossDistanceMatrix(mat, "rocs", ["a", "b", "c"], ["x", "y"], {"k": 1})
    path = tmp_path / "cross.npz"
    cdm.to_file(str(path))

    loaded = CrossDistanceMatrix.from_file(str(path))
    np.testing.assert_array_equal(loaded.matrix, mat)
    assert loaded.shape == (3, 2)
    assert loaded.labels_a == ["a", "b", "c"]
    assert loaded.labels_b == ["x", "y"]
    assert loaded.comparison_name == "rocs"
    assert loaded.params == {"k": 1}


def test_cross_from_file_rejects_symmetric(tmp_path):
    """CrossDistanceMatrix.from_file rejects a symmetric (no matrix_kind) file."""
    from oecluster import CrossDistanceMatrix, DenseStorage, SymmetricDistanceMatrix

    storage = DenseStorage(3)
    storage.Set(0, 1, 0.5)
    storage.Set(0, 2, 0.25)
    storage.Set(1, 2, 0.75)
    sdm = SymmetricDistanceMatrix(storage, "test", ["a", "b", "c"], {})
    path = tmp_path / "sym.npz"
    sdm.to_file(str(path))

    with pytest.raises(ValueError, match="symmetric"):
        CrossDistanceMatrix.from_file(str(path))


def test_cross_from_file_rejects_malformed(tmp_path):
    """A cross file with a wrong-shape matrix raises a descriptive ValueError."""
    from oecluster import CrossDistanceMatrix

    path = tmp_path / "bad.npz"
    # matrix_kind says cross, but matrix shape disagrees with n_a/n_b.
    np.savez_compressed(
        path,
        matrix_kind=np.array("cross"),
        matrix=np.zeros((2, 2), dtype=np.float64),
        n_a=np.array(2),
        n_b=np.array(3),
        labels_a=np.array(["a", "b"]),
        labels_b=np.array(["x", "y", "z"]),
        comparison_name=np.array("fingerprint"),
        params_json=np.array("{}"),
    )
    with pytest.raises(ValueError, match="Malformed cross matrix"):
        CrossDistanceMatrix.from_file(str(path))
