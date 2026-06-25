"""Tests for the clustering-quality report."""

import math

import pytest


def _two_cluster_dm():
    """Two clusters {0,1},{2,3}; intra 0.2, cross 0.8."""
    from oecluster import DenseStorage, DistanceMatrix

    s = DenseStorage(4)
    s.Set(0, 1, 0.2)
    s.Set(2, 3, 0.2)
    s.Set(0, 2, 0.8)
    s.Set(0, 3, 0.8)
    s.Set(1, 2, 0.8)
    s.Set(1, 3, 0.8)
    return DistanceMatrix(s, "test", ["a", "b", "c", "d"], {})


def test_report_basic_and_compactness():
    import oecluster

    dm = _two_cluster_dm()
    result = oecluster.butina(dm, threshold=0.5)
    report = oecluster.cluster_report(result, dm)

    assert isinstance(report, oecluster.ClusterReport)
    assert report.num_samples == 4
    assert report.num_clusters == 2
    assert report.num_noise == 0
    assert report.largest_cluster_fraction == 0.5
    assert math.isclose(report.silhouette, 0.75, rel_tol=1e-9)
    assert math.isclose(report.mean_intra_distance, 0.2, rel_tol=1e-9)


def test_coverage_alignment_and_overrides():
    import oecluster

    dm = _two_cluster_dm()
    result = oecluster.butina(dm, threshold=0.5)
    report = oecluster.cluster_report(
        result, dm, coverage_thresholds=[0.1, 0.3])

    assert report.coverage_thresholds == (0.1, 0.3)
    assert len(report.coverage_at) == 2
    # Every point is within 0.2 of its medoid -> coverage 1.0 at 0.3, and at
    # 0.1 the non-medoid member (0.2 away) is not covered -> 0.5.
    assert math.isclose(report.coverage_at[1], 1.0, rel_tol=1e-9)
    assert math.isclose(report.coverage_at[0], 0.5, rel_tol=1e-9)


def test_preset_changes_thresholds():
    import oecluster

    dm = _two_cluster_dm()
    result = oecluster.butina(dm, threshold=0.5)
    tight = oecluster.cluster_report(result, dm, preset="tight")
    diversity = oecluster.cluster_report(result, dm, preset="diversity")

    assert tight.coverage_thresholds == (0.20, 0.30, 0.40)
    assert diversity.coverage_thresholds == (0.40, 0.50, 0.60)


def test_invalid_preset_and_types():
    import oecluster

    dm = _two_cluster_dm()
    result = oecluster.butina(dm, threshold=0.5)
    with pytest.raises(ValueError, match="preset"):
        oecluster.cluster_report(result, dm, preset="nope")
    with pytest.raises(TypeError):
        oecluster.cluster_report(result, "not a matrix")


def test_report_is_read_only():
    import oecluster

    dm = _two_cluster_dm()
    report = oecluster.cluster_report(oecluster.butina(dm, threshold=0.5), dm)
    with pytest.raises(AttributeError):
        report.num_clusters = 99
    with pytest.raises(AttributeError):
        report.coverage_at = ()  # pyright: ignore[reportAttributeAccessIssue]


def test_result_method_names():
    import oecluster

    dm = _two_cluster_dm()
    assert oecluster.butina(dm, threshold=0.5).method == "butina"
    assert oecluster.dbscan(dm, eps=0.5, min_samples=1).method == "dbscan"
    assert oecluster.agglomerative(dm, n_clusters=2).method == "agglomerative"


def test_report_captures_method():
    import oecluster

    dm = _two_cluster_dm()
    rb = oecluster.cluster_report(oecluster.butina(dm, threshold=0.5), dm)
    assert rb.method == "butina"
    assert "method='butina'" in repr(rb)


def test_compare_reports_multi_and_labels():
    import oecluster

    dm = _two_cluster_dm()
    rb = oecluster.cluster_report(oecluster.butina(dm, threshold=0.5), dm)
    rd = oecluster.cluster_report(oecluster.dbscan(dm, eps=0.5, min_samples=1), dm)
    ra = oecluster.cluster_report(oecluster.agglomerative(dm, n_clusters=2), dm)

    cmp = oecluster.compare_reports(rb, rd, ra)
    assert isinstance(cmp, oecluster.ClusterReportComparison)
    assert cmp.reports == (rb, rd, ra)
    assert not hasattr(cmp, "a")

    rows = cmp.to_table()
    names = [r[0] for r in rows]
    assert "num_clusters" in names
    assert "num_noise" in names
    # one value per report -> length 4
    assert all(len(r) == 4 for r in rows)
    # method names appear as column headers in the repr
    header = repr(cmp).splitlines()[0]
    assert "butina" in header
    assert "dbscan" in header
    assert "agglomerative" in header


def test_compare_reports_validation():
    import oecluster

    dm = _two_cluster_dm()
    rb = oecluster.cluster_report(oecluster.butina(dm, threshold=0.5), dm)
    with pytest.raises(ValueError):
        oecluster.compare_reports(rb)            # fewer than two
    with pytest.raises(TypeError):
        oecluster.compare_reports(rb, "not a report")
