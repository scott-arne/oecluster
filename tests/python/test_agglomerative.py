"""scikit-learn parity tests for native agglomerative clustering."""

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


def _cluster_members(labels):
    """Return label-independent cluster member tuples."""
    labels = np.asarray(labels)
    return tuple(
        sorted(
            tuple(int(index) for index in np.flatnonzero(labels == label))
            for label in np.unique(labels)
        )
    )


def _sklearn_agglomerative(square, *, linkage, n_clusters=2, distance_threshold=None):
    from sklearn.cluster import AgglomerativeClustering

    kwargs = {
        "linkage": linkage,
        "compute_full_tree": True,
    }
    if distance_threshold is None:
        kwargs["n_clusters"] = n_clusters
    else:
        kwargs["n_clusters"] = None
        kwargs["distance_threshold"] = distance_threshold

    try:
        model = AgglomerativeClustering(metric="precomputed", **kwargs)
    except TypeError:
        model = AgglomerativeClustering(affinity="precomputed", **kwargs)  # pyright: ignore[reportCallIssue]
    return model.fit(square)


@pytest.mark.parametrize("linkage", ["single", "complete", "average"])
def test_agglomerative_matches_sklearn_precomputed_linkages(linkage):
    import oecluster

    x = np.array([0.0, 0.1, 0.2, 5.0, 5.2, 5.4], dtype=np.float64)
    square = np.abs(x[:, None] - x[None, :])
    dm = _dense_distance_matrix(square)

    expected = _sklearn_agglomerative(square, linkage=linkage, n_clusters=2)
    observed = oecluster.agglomerative(dm, linkage=linkage, n_clusters=2)

    assert isinstance(observed, oecluster.AgglomerativeResult)
    assert isinstance(observed, oecluster.ClusteringResult)
    assert _cluster_members(observed.labels) == _cluster_members(expected.labels_)
    assert observed.clusters == ((0, 1, 2), (3, 4, 5))
    assert len(observed.children) == square.shape[0] - 1
    assert observed.distances is not None
    assert observed.distances.shape == (square.shape[0] - 1,)
    assert observed.cluster_sizes[-1] == square.shape[0]


def test_agglomerative_distance_threshold_matches_sklearn():
    import oecluster

    square = np.array(
        [
            [0.0, 0.1, 5.0, 5.2],
            [0.1, 0.0, 4.8, 5.0],
            [5.0, 4.8, 0.0, 0.2],
            [5.2, 5.0, 0.2, 0.0],
        ],
        dtype=np.float64,
    )
    dm = _dense_distance_matrix(square)

    expected = _sklearn_agglomerative(
        square,
        linkage="average",
        distance_threshold=0.15,
    )
    observed = oecluster.agglomerative(
        dm,
        linkage="average",
        distance_threshold=0.15,
    )

    assert _cluster_members(observed.labels) == _cluster_members(expected.labels_)
    assert observed.clusters == ((0, 1), (2,), (3,))


def test_agglomerative_weighted_linkage_returns_tree_metadata():
    import oecluster

    square = np.array(
        [
            [0.0, 0.1, 0.2, 10.0, 11.0],
            [0.1, 0.0, 0.3, 12.0, 13.0],
            [0.2, 0.3, 0.0, 20.0, 21.0],
            [10.0, 12.0, 20.0, 0.0, 0.4],
            [11.0, 13.0, 21.0, 0.4, 0.0],
        ],
        dtype=np.float64,
    )
    dm = _dense_distance_matrix(square)

    observed = oecluster.agglomerative(dm, linkage="weighted", n_clusters=2)

    assert observed.labels.tolist() == [0, 0, 0, 1, 1]
    assert observed.clusters == ((0, 1, 2), (3, 4))
    assert observed.children == ((0, 1), (2, 5), (3, 4), (6, 7))
    assert observed.distances is not None
    np.testing.assert_allclose(observed.distances, [0.1, 0.25, 0.4, 16.0])
    assert observed.cluster_sizes == (2, 3, 2, 5)


def test_agglomerative_can_stop_before_full_tree():
    import oecluster

    square = np.array(
        [
            [0.0, 0.1, 5.0, 5.2],
            [0.1, 0.0, 4.8, 5.0],
            [5.0, 4.8, 0.0, 0.2],
            [5.2, 5.0, 0.2, 0.0],
        ],
        dtype=np.float64,
    )
    dm = _dense_distance_matrix(square)

    observed = oecluster.agglomerative(
        dm,
        linkage="average",
        n_clusters=2,
        compute_full_tree=False,
    )

    assert observed.labels.tolist() == [0, 0, 1, 1]
    assert observed.children == ((0, 1), (2, 3))
    assert observed.distances is not None
    np.testing.assert_allclose(observed.distances, [0.1, 0.2])
    assert observed.cluster_sizes == (2, 2)


def test_agglomerative_rejects_invalid_inputs():
    import oecluster

    dm = _dense_distance_matrix([[0.0, 1.0], [1.0, 0.0]])

    with pytest.raises(ValueError, match="n_clusters"):
        oecluster.agglomerative(dm, n_clusters=0)
    with pytest.raises(ValueError, match="distance_threshold"):
        oecluster.agglomerative(dm, distance_threshold=-0.1)
    with pytest.raises(ValueError, match="linkage"):
        oecluster.agglomerative(dm, linkage="ward")
    with pytest.raises(TypeError, match="DistanceMatrix"):
        oecluster.agglomerative(object())
