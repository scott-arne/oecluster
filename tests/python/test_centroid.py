"""Tests for cluster representative selection."""

import pytest


def _dense_distance_matrix(values, n):
    """Build a dense DistanceMatrix from explicit upper-triangle values."""
    from oecluster import DenseStorage
    from oecluster import DistanceMatrix

    storage = DenseStorage(n)
    for (i, j), value in values.items():
        storage.Set(i, j, value)
    return DistanceMatrix(storage, "test", [f"item_{i}" for i in range(n)], {})


def test_centroid_methods_select_expected_representatives():
    """Choose first, medoid/mean, and minimax representatives."""
    import oecluster

    dm = _dense_distance_matrix(
        {
            (0, 1): 0.1,
            (0, 2): 0.1,
            (0, 3): 1.0,
            (1, 2): 0.6,
            (1, 3): 0.6,
            (2, 3): 0.6,
        },
        4,
    )
    cluster = (0, 1, 2, 3)

    assert oecluster.centroid(cluster, dm, method="first") == 0
    assert oecluster.centroid(cluster, dm, method="medoid") == 0
    assert oecluster.centroid(cluster, dm, method="mean") == 0
    assert oecluster.centroid(cluster, dm, method="minimax") == 1


def test_centroid_tie_breaks_by_cluster_order():
    """Tie-breaking preserves the original cluster member order."""
    import oecluster

    dm = _dense_distance_matrix(
        {
            (0, 1): 1.0,
            (0, 2): 1.0,
            (1, 2): 0.1,
        },
        3,
    )

    assert oecluster.centroid((0, 1, 2), dm, method="medoid") == 1
    assert oecluster.centroid((0, 2, 1), dm, method="medoid") == 2


def test_centroid_accepts_butina_cluster_result():
    """Use centroid() directly on clusters returned by butina()."""
    import oecluster

    dm = _dense_distance_matrix(
        {
            (0, 1): 0.1,
            (0, 2): 0.1,
            (1, 2): 0.8,
        },
        3,
    )
    cluster = oecluster.butina(dm, threshold=0.2)[0]

    assert oecluster.centroid(cluster, dm, method="first") == cluster[0]
    assert oecluster.centroid(cluster, dm, method="medoid") == 0


def test_centroid_rejects_invalid_inputs():
    """Raise Python-level errors for invalid wrapper inputs."""
    import oecluster

    dm = _dense_distance_matrix({(0, 1): 0.1}, 2)

    with pytest.raises(ValueError, match="empty"):
        oecluster.centroid((), dm)
    with pytest.raises(ValueError, match="method"):
        oecluster.centroid((0, 1), dm, method="bogus")
    with pytest.raises(TypeError, match="DistanceMatrix"):
        oecluster.centroid((0, 1), object())


def test_centroid_sparse_storage_only_supports_first_method():
    """Distance-based representatives need complete pairwise distances."""
    import oecluster

    storage = oecluster.SparseStorage(3, 0.2)
    storage.Set(0, 1, 0.1)
    storage.Set(0, 2, 0.3)
    storage.Set(1, 2, 0.3)
    storage.Finalize()
    dm = oecluster.DistanceMatrix(storage, "test", [], {})

    assert oecluster.centroid((1, 0, 2), dm, method="first") == 1
    with pytest.raises(RuntimeError, match="SparseStorage"):
        oecluster.centroid((1, 0, 2), dm, method="medoid")
