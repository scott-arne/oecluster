"""Tests for cluster representative selection."""

import math

import pytest


def _dense_distance_matrix(values, n):
    """Build a dense DistanceMatrix from explicit upper-triangle values."""
    from oecluster import DenseStorage
    from oecluster import DistanceMatrix

    storage = DenseStorage(n)
    for (i, j), value in values.items():
        storage.Set(i, j, value)
    return DistanceMatrix(storage, "test", [f"item_{i}" for i in range(n)], {})


def test_centroid_api_is_not_exported():
    """Representative terminology replaces the old centroid helper."""
    import oecluster

    assert "centroid" not in oecluster.__all__
    assert not hasattr(oecluster, "centroid")


def test_representative_methods_select_expected_members():
    """Choose medoid, minimax, and highest-neighborhood representatives."""
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

    assert oecluster.representative(cluster, dm, method="medoid") == 0
    assert oecluster.representative(cluster, dm, method="minimax") == 1
    assert (
        oecluster.representative(
            cluster,
            dm,
            method="highest_neighborhood",
            threshold=0.2,
        )
        == 0
    )


def test_rank_representatives_reports_quality_metrics():
    """Return a complete score-ranked list with per-representative metrics."""
    import oecluster

    dm = _dense_distance_matrix(
        {
            (0, 1): 0.1,
            (0, 2): 0.1,
            (0, 3): 1.0,
            (0, 4): 0.9,
            (1, 2): 0.6,
            (1, 3): 0.6,
            (1, 4): 0.7,
            (2, 3): 0.6,
            (2, 4): 0.8,
            (3, 4): 0.2,
        },
        5,
    )

    ranked = oecluster.rank_representatives(
        (0, 1, 2, 3),
        dm,
        method="medoid",
        threshold=0.2,
        scaffold_labels=("A", "A", "B", "B", "external"),
    )

    assert [entry.member for entry in ranked] == [0, 1, 2, 3]
    assert ranked[0].score == pytest.approx(0.4)
    metrics = ranked[0].metrics
    assert metrics.mean_distance_to_cluster == pytest.approx(0.4)
    assert metrics.max_distance_to_cluster == pytest.approx(1.0)
    assert metrics.median_distance_to_cluster == pytest.approx(0.1)
    assert metrics.neighbor_fraction_at_threshold == pytest.approx(0.75)
    assert metrics.nearest_external_distance == pytest.approx(0.9)
    assert metrics.cluster_radius == pytest.approx(1.0)
    assert metrics.cluster_diameter == pytest.approx(1.0)
    assert metrics.silhouette_like_score == pytest.approx((0.9 - 0.4) / 0.9)
    assert metrics.scaffold_purity == pytest.approx(0.5)
    assert metrics.representative_rank == 1


def test_weighted_medoid_uses_penalty_and_priority_vectors():
    """Bias medoid ranking toward practical metadata when requested."""
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

    ranked = oecluster.rank_representatives(
        (0, 1, 2, 3),
        dm,
        method="weighted_medoid",
        alpha=1.0,
        beta=1.0,
        gamma=1.0,
        liability_penalties=(0.0, 0.0, 0.0, 0.0),
        priority_scores=(0.0, 0.4, 0.0, 0.0),
    )

    assert ranked[0].member == 1
    assert ranked[0].score == pytest.approx((0.1 + 0.6 + 0.6) / 3.0 - 0.4)


def test_select_representatives_supports_score_and_diversity_modes():
    """Allow users to choose top-score or coverage-oriented k representatives."""
    import oecluster

    dm = _dense_distance_matrix(
        {
            (0, 1): 0.1,
            (0, 2): 0.8,
            (0, 3): 0.8,
            (1, 2): 0.8,
            (1, 3): 0.8,
            (2, 3): 0.2,
        },
        4,
    )
    cluster = (0, 1, 2, 3)

    score_selected = oecluster.select_representatives(
        cluster,
        dm,
        k=2,
        method="medoid",
        selection="score",
    )
    diversity_selected = oecluster.select_representatives(
        cluster,
        dm,
        k=2,
        method="medoid",
        selection="diversity",
    )

    assert [entry.member for entry in score_selected] == [0, 1]
    assert [entry.member for entry in diversity_selected] == [0, 2]
    assert [entry.metrics.representative_rank for entry in diversity_selected] == [1, 3]


def test_representative_rejects_invalid_inputs():
    """Raise Python-level errors for invalid wrapper inputs."""
    import oecluster

    dm = _dense_distance_matrix({(0, 1): 0.1}, 2)

    with pytest.raises(ValueError, match="empty"):
        oecluster.representative((), dm)
    with pytest.raises(ValueError, match="method"):
        oecluster.representative((0, 1), dm, method="bogus")
    with pytest.raises(ValueError, match="threshold"):
        oecluster.representative((0, 1), dm, method="highest_neighborhood")
    with pytest.raises(ValueError, match="selection"):
        oecluster.select_representatives((0, 1), dm, k=1, selection="bogus")
    with pytest.raises(ValueError, match="same length"):
        oecluster.rank_representatives(
            (0, 1),
            dm,
            method="weighted_medoid",
            liability_penalties=(0.0,),
        )
    with pytest.raises(TypeError, match="DistanceMatrix"):
        oecluster.representative((0, 1), object())


def test_metrics_without_optional_inputs_use_nan_or_infinity():
    """Optional quality metrics have explicit missing-value behavior."""
    import oecluster

    dm = _dense_distance_matrix({(0, 1): 0.1}, 2)
    ranked = oecluster.rank_representatives((0, 1), dm, method="medoid")

    assert math.isnan(ranked[0].metrics.neighbor_fraction_at_threshold)
    assert math.isinf(ranked[0].metrics.nearest_external_distance)
    assert ranked[0].metrics.silhouette_like_score == pytest.approx(1.0)
    assert math.isnan(ranked[0].metrics.scaffold_purity)
