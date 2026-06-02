"""RDKit parity tests for Butina clustering."""

import pytest


def _dense_distance_matrix(values, n):
    """Build a dense DistanceMatrix from explicit upper-triangle values."""
    from oecluster import DenseStorage
    from oecluster import DistanceMatrix

    storage = DenseStorage(n)
    for (i, j), value in values.items():
        storage.Set(i, j, value)
    return DistanceMatrix(storage, "test", [f"item_{i}" for i in range(n)], {})


def _rdkit_butina(distance_matrix, threshold, *, reordering=False):
    """Return RDKit Butina clusters using square form to avoid order drift."""
    pytest.importorskip("rdkit")
    from rdkit.ML.Cluster import Butina

    return Butina.ClusterData(
        distance_matrix.squareform(),
        distance_matrix.num_items,
        threshold,
        isDistData=True,
        reordering=reordering,
    )


def _assert_labels_match_clusters(result, num_items):
    """Every member's label equals its cluster index, over the full item range."""
    assert len(result.labels) == num_items
    for cluster_index, cluster in enumerate(result.clusters):
        for member in cluster:
            assert result.labels[member] == cluster_index


def test_butina_matches_rdkit_without_reordering():
    """Match RDKit's greedy Butina clustering without count updates."""
    import oecluster

    dm = _dense_distance_matrix(
        {
            (0, 1): 0.10,
            (0, 2): 0.40,
            (0, 3): 0.90,
            (0, 4): 0.90,
            (1, 2): 0.20,
            (1, 3): 0.80,
            (1, 4): 0.90,
            (2, 3): 0.30,
            (2, 4): 0.90,
            (3, 4): 0.90,
        },
        5,
    )

    observed = oecluster.butina(dm, threshold=0.35, reordering=False)
    expected = _rdkit_butina(dm, 0.35, reordering=False)

    assert observed.clusters == expected
    _assert_labels_match_clusters(observed, dm.num_items)


def test_butina_matches_rdkit_with_reordering():
    """Match RDKit's optional post-cluster count reordering."""
    import oecluster

    dm = _dense_distance_matrix(
        {
            (0, 1): 0.10,
            (0, 2): 0.10,
            (0, 3): 0.90,
            (0, 4): 0.90,
            (0, 5): 0.90,
            (1, 2): 0.90,
            (1, 3): 0.10,
            (1, 4): 0.90,
            (1, 5): 0.90,
            (2, 3): 0.90,
            (2, 4): 0.10,
            (2, 5): 0.90,
            (3, 4): 0.90,
            (3, 5): 0.10,
            (4, 5): 0.10,
        },
        6,
    )

    observed = oecluster.butina(dm, threshold=0.2, reordering=True)
    expected = _rdkit_butina(dm, 0.2, reordering=True)

    assert observed.clusters == expected
    _assert_labels_match_clusters(observed, dm.num_items)


def test_butina_fingerprint_distance_matrix_matches_rdkit(
    aspirin_mol,
    ethanol_mol,
):
    """Cluster a real fingerprint distance matrix with RDKit parity."""
    import oecluster

    dm = oecluster.pdist([aspirin_mol, ethanol_mol], "fingerprint")

    observed = oecluster.butina(dm, threshold=1.0)
    expected = _rdkit_butina(dm, 1.0)

    assert observed.clusters == expected
    _assert_labels_match_clusters(observed, dm.num_items)


def test_butina_rejects_negative_threshold():
    """Reject invalid thresholds before calling the C++ algorithm."""
    import oecluster

    dm = _dense_distance_matrix({(0, 1): 0.1}, 2)

    with pytest.raises(ValueError, match="threshold"):
        oecluster.butina(dm, threshold=-0.1)


def test_butina_rejects_sparse_cutoff_below_threshold():
    """Reject sparse inputs that cannot contain all threshold neighbors."""
    import oecluster

    storage = oecluster.SparseStorage(3, 0.2)
    storage.Set(0, 1, 0.1)
    storage.Set(0, 2, 0.3)
    storage.Set(1, 2, 0.3)
    storage.Finalize()
    dm = oecluster.DistanceMatrix(storage, "test", [], {})

    with pytest.raises(RuntimeError, match="cutoff"):
        oecluster.butina(dm, threshold=0.3)
