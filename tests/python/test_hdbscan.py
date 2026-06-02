"""scikit-learn parity tests for native HDBSCAN."""

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


def test_hdbscan_matches_sklearn_precomputed_toy():
    import oecluster
    from sklearn.cluster import HDBSCAN as SklearnHDBSCAN

    x = np.array([0.0, 0.1, 0.2, 5.0, 5.1, 5.2], dtype=np.float64)
    square = np.abs(x[:, None] - x[None, :])
    dm = _dense_distance_matrix(square)

    expected = SklearnHDBSCAN(
        min_cluster_size=2,
        min_samples=2,
        metric="precomputed",
        cluster_selection_method="leaf",
        allow_single_cluster=False,
        copy=True,
    ).fit(square)
    observed = oecluster.hdbscan(
        dm,
        min_cluster_size=2,
        min_samples=2,
        cluster_selection_method="leaf",
    )

    assert isinstance(observed, oecluster.HDBSCANResult)
    assert isinstance(observed, oecluster.ClusteringResult)
    assert observed.labels.tolist() == expected.labels_.tolist()
    np.testing.assert_allclose(observed.probabilities, expected.probabilities_)
    assert observed.clusters == ((0, 1, 2), (3, 4, 5))


def test_hdbscan_matches_sklearn_with_noise_and_eom():
    import oecluster
    from sklearn.cluster import HDBSCAN as SklearnHDBSCAN

    x = np.array([0.0, 0.1, 0.2, 5.0, 5.1, 5.2, 20.0], dtype=np.float64)
    square = np.abs(x[:, None] - x[None, :])
    dm = _dense_distance_matrix(square)

    expected = SklearnHDBSCAN(
        min_cluster_size=2,
        min_samples=2,
        metric="precomputed",
        cluster_selection_method="eom",
        allow_single_cluster=False,
        copy=True,
    ).fit(square)
    observed = oecluster.hdbscan(
        dm,
        min_cluster_size=2,
        min_samples=2,
        cluster_selection_method="eom",
    )

    assert observed.labels.tolist() == expected.labels_.tolist()
    np.testing.assert_allclose(observed.probabilities, expected.probabilities_)
    assert observed.clusters == ((0, 1, 2), (3, 4, 5))


def test_hdbscan_matches_sklearn_labels_for_tied_mst_edges():
    import oecluster
    from sklearn.cluster import HDBSCAN as SklearnHDBSCAN

    rng = np.random.default_rng(7)
    n_samples = 120
    n_clusters = 10
    centers = np.linspace(0.0, 10.0 * (n_clusters - 1), n_clusters)
    assignments = np.arange(n_samples) % n_clusters
    x = np.sort(centers[assignments] + rng.normal(0.0, 0.03, n_samples))
    square = np.abs(x[:, None] - x[None, :])
    dm = _dense_distance_matrix(square)

    expected = SklearnHDBSCAN(
        min_cluster_size=5,
        min_samples=5,
        metric="precomputed",
        cluster_selection_method="eom",
        allow_single_cluster=False,
        copy=True,
    ).fit(square)
    observed = oecluster.hdbscan(
        dm,
        min_cluster_size=5,
        min_samples=5,
        cluster_selection_method="eom",
    )

    assert observed.labels.tolist() == expected.labels_.tolist()
    assert observed.probabilities.shape == expected.probabilities_.shape
    assert np.all(observed.probabilities >= 0.0)
    assert np.all(observed.probabilities <= 1.0)


def test_hdbscan_rejects_invalid_inputs():
    import oecluster

    dm = _dense_distance_matrix([[0.0, 1.0], [1.0, 0.0]])

    with pytest.raises(ValueError, match="min_cluster_size"):
        oecluster.hdbscan(dm, min_cluster_size=1)
    with pytest.raises(ValueError, match="min_samples"):
        oecluster.hdbscan(dm, min_cluster_size=2, min_samples=0)
    with pytest.raises(ValueError, match="cluster_selection_method"):
        oecluster.hdbscan(dm, min_cluster_size=2, cluster_selection_method="flat")
    with pytest.raises(TypeError, match="DistanceMatrix"):
        oecluster.hdbscan(object(), min_cluster_size=2)
