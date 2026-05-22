"""scikit-learn parity tests for native DBSCAN."""

import numpy as np
import pytest


def _dense_distance_matrix(square):
    """Build a dense DistanceMatrix from a square distance matrix."""
    from oecluster import DenseStorage
    from oecluster import DistanceMatrix

    square = np.asarray(square, dtype=np.float64)
    storage = DenseStorage(square.shape[0])
    for i in range(square.shape[0]):
        for j in range(i + 1, square.shape[0]):
            storage.Set(i, j, float(square[i, j]))
    return DistanceMatrix(storage, "test", [f"item_{i}" for i in range(square.shape[0])], {})


def test_dbscan_matches_sklearn_precomputed_toy():
    import oecluster
    from sklearn.cluster import dbscan as sklearn_dbscan

    x = np.array([0, 2, 3, 4, 6, 8, 10], dtype=np.float64)
    square = np.abs(x[:, None] - x[None, :])
    dm = _dense_distance_matrix(square)

    expected_core, expected_labels = sklearn_dbscan(
        square, eps=1.0, min_samples=3, metric="precomputed"
    )
    observed = oecluster.dbscan(dm, eps=1.0, min_samples=3)

    assert observed.core_sample_indices == tuple(int(i) for i in expected_core)
    assert observed.labels.tolist() == expected_labels.tolist()
    assert observed.clusters == ((1, 2, 3),)


def test_dbscan_rejects_invalid_inputs():
    import oecluster

    dm = _dense_distance_matrix([[0.0, 1.0], [1.0, 0.0]])

    with pytest.raises(ValueError, match="eps"):
        oecluster.dbscan(dm, eps=-0.1)
    with pytest.raises(ValueError, match="min_samples"):
        oecluster.dbscan(dm, eps=0.5, min_samples=0)
    with pytest.raises(TypeError, match="DistanceMatrix"):
        oecluster.dbscan(object(), eps=0.5)
