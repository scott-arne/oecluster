"""Shared Python clustering result hierarchy behavior."""

import numpy as np


def _toy_distance_matrix():
    """Build a small DistanceMatrix with two clear clusters of two items."""
    from oecluster import DenseStorage, DistanceMatrix

    storage = DenseStorage(4)
    storage.Set(0, 1, 0.1)
    storage.Set(0, 2, 0.9)
    storage.Set(0, 3, 0.9)
    storage.Set(1, 2, 0.9)
    storage.Set(1, 3, 0.9)
    storage.Set(2, 3, 0.1)
    return DistanceMatrix(storage, "test", ["a", "b", "c", "d"], {})


def test_base_normalizes_labels_and_clusters():
    import oecluster

    result = oecluster.ClusteringResult(
        labels=[0, -1, 1, 0],
        clusters=((0, 3), (2,)),
    )

    assert isinstance(result.labels, np.ndarray)
    assert result.labels.dtype == np.intp
    assert result.labels.tolist() == [0, -1, 1, 0]
    assert result.clusters == ((0, 3), (2,))


def test_dunders_reflect_clusters():
    import oecluster

    result = oecluster.ClusteringResult(
        labels=[0, 0, 1, 1],
        clusters=((0, 1), (2, 3)),
    )

    assert len(result) == 2
    assert result.num_clusters == 2
    assert result.num_items == 4
    assert list(iter(result)) == [(0, 1), (2, 3)]
    assert result[0] == (0, 1)
    assert result[1] == (2, 3)
    assert repr(result) == "ClusteringResult(num_clusters=2, num_items=4)"


def test_butina_returns_butina_result_subclass():
    import oecluster

    dm = _toy_distance_matrix()
    result = oecluster.butina(dm, threshold=0.2)

    assert isinstance(result, oecluster.ButinaResult)
    assert isinstance(result, oecluster.ClusteringResult)
    # Butina carries no algorithm-specific fields.
    assert not hasattr(result, "centroids")
    assert not hasattr(result, "probabilities")
    assert not hasattr(result, "core_sample_indices")
    assert repr(result).startswith("ButinaResult(")


def test_dbscan_returns_dbscan_result_subclass():
    import oecluster

    dm = _toy_distance_matrix()
    result = oecluster.dbscan(dm, eps=0.2, min_samples=2)

    assert isinstance(result, oecluster.DBSCANResult)
    assert isinstance(result, oecluster.ClusteringResult)
    assert isinstance(result.core_sample_indices, tuple)
    # DBSCAN does not carry centroids or probabilities.
    assert not hasattr(result, "centroids")
    assert not hasattr(result, "probabilities")
