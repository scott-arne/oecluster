"""Tests for the clustering-quality report."""

import math

import pytest


def _two_cluster_dm():
    """Two clusters {0,1},{2,3}; intra 0.2, cross 0.8."""
    from oecluster import DenseStorage, SymmetricDistanceMatrix

    s = DenseStorage(4)
    s.Set(0, 1, 0.2)
    s.Set(2, 3, 0.2)
    s.Set(0, 2, 0.8)
    s.Set(0, 3, 0.8)
    s.Set(1, 2, 0.8)
    s.Set(1, 3, 0.8)
    return SymmetricDistanceMatrix(s, "test", ["a", "b", "c", "d"], {})


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


def test_scalar_fields_mirror_the_native_struct():
    """A field added in C++ and forgotten in _SCALAR_FIELDS is invisible with a
    green suite, which is the gap this closes."""
    import oecluster
    from oecluster import oecluster as _native

    native = _native.ClusterReport()
    scalar_names = {
        name
        for name in dir(native)
        if not name.startswith("_")
        and isinstance(getattr(native, name), (int, float))
        and not isinstance(getattr(native, name), bool)
    }
    assert set(oecluster.ClusterReport._SCALAR_FIELDS) == scalar_names


def test_new_types_are_public():
    import oecluster

    assert "ClusterRecord" in oecluster.__all__
    assert "ClusterReportRequested" in oecluster.__all__
    assert oecluster.ClusterRecord is not None
    assert oecluster.ClusterReportRequested is not None


def test_default_constructed_record_reads_undefined_not_zero():
    """Exporting ClusterRecord exports its default constructor, so a caller can
    build or preallocate one. The four undefined-valued fields must read NaN
    there, not 0.0 -- a zero nearest_cluster_distance would say another cluster
    sits at zero distance. This is the Python half of the C++ assertion in
    ClusterReportTest.NewSurfaceDefaultsToUnrequested.
    """
    import math

    import oecluster
    import oecluster.oecluster as _native

    # ClusterRecordVector comes off the extension module: Task 2 instantiates
    # the SWIG template but no task adds the vector type to the package's
    # __all__, so only ClusterRecord itself is re-exported at package level.
    for record in (
        oecluster.ClusterRecord(),
        _native.ClusterRecordVector(1)[0],
    ):
        assert math.isnan(record.mean_intra_distance)
        assert math.isnan(record.median_intra_distance)
        assert math.isnan(record.nearest_cluster_distance)
        assert math.isnan(record.silhouette)
        # -1 literal, not oecluster.NO_NEAREST_CLUSTER: the sentinel is a C++
        # constant and no task exports it to the Python package namespace.
        assert record.nearest_cluster == -1
        # 0.0 is the real singleton value for these two, not a placeholder.
        assert record.radius == 0.0
        assert record.diameter == 0.0


@pytest.mark.parametrize(
    "flag", ["compute_pair_rank_indices", "compute_per_cluster_records"]
)
def test_non_bool_flags_are_rejected_by_name(flag):
    """The match includes "must be a bool" deliberately. Before the signature
    gains the keyword, Python's own "unexpected keyword argument" TypeError
    also names the flag, so matching on the flag alone would pass green
    against an unimplemented feature."""
    import oecluster

    dm = _two_cluster_dm()
    result = oecluster.butina(dm, threshold=0.5)
    with pytest.raises(TypeError, match=f"{flag} must be a bool"):
        oecluster.cluster_report(result, dm, **{flag: "yes"})


def test_flag_type_error_outranks_the_matrix_pairing_check():
    """The local-argument block runs before the pairing check, so an
    authoritative complaint about what was typed is never pre-empted by an
    advisory one about the matrix. As above, the match pins the body's own
    message rather than the flag name, which an unexpected-keyword TypeError
    would also carry."""
    import oecluster
    from oecluster import DenseStorage, SymmetricDistanceMatrix

    dm = _two_cluster_dm()
    result = oecluster.butina(dm, threshold=0.5)
    bigger = SymmetricDistanceMatrix(
        DenseStorage(6), "test", ["a", "b", "c", "d", "e", "f"], {})
    with pytest.raises(
            TypeError, match="compute_pair_rank_indices must be a bool"):
        oecluster.cluster_report(result, bigger, compute_pair_rank_indices=1.5)


def test_sparse_storage_still_outranks_a_non_bool_allow_nonmetric():
    """The other half of the ordering rule, and the half that must not move.
    The sparse refusal already sits ahead of the metric gate, and the gate is
    where a malformed allow_nonmetric is caught; inserting the two new flag
    checks above both of them must leave this pair in the order it had."""
    import oecluster
    from oecluster import SparseStorage, SymmetricDistanceMatrix

    dm = _two_cluster_dm()
    result = oecluster.butina(dm, threshold=0.5)

    # Four samples, so the result/matrix pairing check -- which outranks the
    # storage refusal -- passes and leaves sparse as the first real problem.
    storage = SparseStorage(4, 0.5)
    storage.Set(0, 1, 0.2)
    storage.Set(2, 3, 0.2)
    storage.Finalize()
    sparse_dm = SymmetricDistanceMatrix(
        storage, "test", ["a", "b", "c", "d"], {})

    # ValueError, not TypeError: if the gate ran first this would be a
    # TypeError about allow_nonmetric and the test would fail on the type.
    with pytest.raises(ValueError, match="SparseStorage"):
        oecluster.cluster_report(result, sparse_dm, allow_nonmetric="yes")


def test_records_and_requested_round_trip():
    import oecluster

    dm = _two_cluster_dm()
    result = oecluster.butina(dm, threshold=0.5)

    plain = oecluster.cluster_report(result, dm)
    assert plain.records == ()
    assert plain.requested.per_cluster_records is False
    assert plain.requested.pair_rank_indices is False

    detailed = oecluster.cluster_report(
        result, dm, compute_per_cluster_records=True)
    assert len(detailed.records) == detailed.num_clusters
    assert detailed.requested.per_cluster_records is True
    for ordinal, record in enumerate(detailed.records):
        assert record.label == ordinal
    with pytest.raises(AttributeError):
        detailed.records[0].label = 7


def test_noise_coverage_parallels_coverage():
    import oecluster

    dm = _two_cluster_dm()
    report = oecluster.cluster_report(oecluster.butina(dm, threshold=0.5), dm)
    assert len(report.noise_coverage_at) == len(report.coverage_at)


def test_partition_error_surfaces_as_runtime_error():
    """Every std::exception maps to SWIG_RuntimeError; A1 adds no ValueError
    mapping. The message must survive so the caller learns which index is
    wrong. The result must still cover four samples, or the pairing check
    at __init__.py:3094 fires first and this tests nothing."""
    import oecluster

    dm = _two_cluster_dm()
    broken = oecluster.ClusteringResult([0, 0, 0, 0], [[0, 1]])
    with pytest.raises(RuntimeError, match="sample 2"):
        oecluster.cluster_report(broken, dm)


def test_matrix_gate_outranks_the_native_partition_refusal():
    """The companion to the test above, and the one that keeps the docstring
    honest. Native cluster_report ranks the partition error above the
    finiteness one, so C++ given this same pair of bad inputs names sample 2.
    Python never reaches that ordering: _gate.require_metric runs first, scans
    the whole matrix, and refuses with a ValueError that names no index. The
    docstring documents this stricter contract, and prose is all that says so
    unless something pins it."""
    import oecluster
    from oecluster import DenseStorage, SymmetricDistanceMatrix

    # Built inline rather than through _two_cluster_dm(), which cannot produce
    # a non-finite matrix. The NaN goes into the storage before construction;
    # SymmetricDistanceMatrix does not validate finiteness on the way in.
    storage = DenseStorage(4)
    storage.Set(0, 1, 0.2)
    storage.Set(2, 3, 0.2)
    storage.Set(0, 2, 0.8)
    storage.Set(0, 3, 0.8)
    storage.Set(1, 2, 0.8)
    storage.Set(1, 3, float("nan"))
    dm = SymmetricDistanceMatrix(storage, "test", ["a", "b", "c", "d"], {})

    broken = oecluster.ClusteringResult([0, 0, 0, 0], [[0, 1]])
    with pytest.raises(ValueError, match="non-finite entries"):
        oecluster.cluster_report(broken, dm)
