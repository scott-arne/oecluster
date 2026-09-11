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


def test_comparison_table_carries_the_request_flag_and_noise_coverage():
    """The row set is a contract, and nothing pinned it: a comparison checking
    only that two known names appear lets rows be added or dropped in silence.
    The request indicator earns its row because c_index and baker_hubert_gamma
    are the only rows whose NaN can mean "nobody asked" rather than "asked and
    undefined", and the comparison path is the one place ClusterReportRequested
    is otherwise unavailable."""
    import oecluster

    dm = _two_cluster_dm()
    result = oecluster.butina(dm, threshold=0.5)
    asked = oecluster.cluster_report(result, dm, compute_pair_rank_indices=True)
    plain = oecluster.cluster_report(result, dm)

    rows = oecluster.compare_reports(asked, plain).to_table()
    names = [row[0] for row in rows]
    values = {row[0]: row[1:] for row in rows}
    assert len(values) == len(rows)

    # 27 scalar metrics, the request indicator, then a coverage row and a
    # noise-coverage row for each of the default preset's three thresholds.
    assert len(rows) == 34
    # Directly beneath the two rows it explains: c_index and
    # baker_hubert_gamma are the last two entries of _SCALAR_FIELDS.
    assert names[names.index("baker_hubert_gamma") + 1] == (
        "requested_pair_rank_indices")
    # No matching row for the other flag, deliberately: per_cluster_records
    # governs records, which this table does not carry.
    assert "requested_per_cluster_records" not in values

    assert values["requested_pair_rank_indices"] == (True, False)
    assert values["c_index"][0] == 0.0
    assert math.isnan(values["c_index"][1])
    assert values["baker_hubert_gamma"][0] == 1.0
    assert math.isnan(values["baker_hubert_gamma"][1])

    for threshold in asked.coverage_thresholds:
        assert values[f"coverage_at[{threshold}]"] == (1.0, 1.0)
        # Compared as reprs rather than through all(), so a failure prints the
        # offending entry. This clustering has no noise, so the noise curve is
        # NaN everywhere while the ordinary one is complete -- which is what
        # makes a noise row sourced from coverage_at visible.
        assert [repr(value)
                for value in values[f"noise_coverage_at[{threshold}]"]] == (
                    ["nan", "nan"])


def test_comparison_table_pads_a_report_with_no_clusters():
    """An all-noise report keeps its threshold list but leaves both coverage
    vectors empty, so the row still has to be produced with a value in it. The
    lookup's fallback covers that, and the noise rows now depend on it too."""
    import oecluster

    dm = _two_cluster_dm()
    clustered = oecluster.cluster_report(
        oecluster.butina(dm, threshold=0.5), dm)
    all_noise = oecluster.cluster_report(
        oecluster.ClusteringResult([-1, -1, -1, -1], []), dm)
    assert all_noise.num_clusters == 0
    assert all_noise.coverage_at == ()
    assert all_noise.noise_coverage_at == ()

    values = {row[0]: row[1:]
              for row in oecluster.compare_reports(clustered, all_noise).to_table()}
    for threshold in clustered.coverage_thresholds:
        assert math.isnan(values[f"coverage_at[{threshold}]"][1])
        assert math.isnan(values[f"noise_coverage_at[{threshold}]"][1])


def test_scalar_fields_mirror_the_native_struct():
    """A field added in C++ and forgotten in _SCALAR_FIELDS is invisible with a
    green suite, which is the gap this closes."""
    import oecluster
    from oecluster import oecluster as _native

    native = _native.ClusterReport()
    # thisown is excluded by name rather than by dropping every bool: it is
    # SWIG's ownership flag and the only bool on the struct today, and a filter
    # on the type would silently swallow a genuinely bool-valued metric added
    # in C++ later.
    scalar_names = {
        name
        for name in dir(native)
        if not name.startswith("_")
        and name != "thisown"
        and isinstance(getattr(native, name), (int, float))
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

    # ClusterRecordVector comes off the extension module: the SWIG template is
    # instantiated there but the vector type is not added to the package's
    # __all__, so only ClusterRecord itself is re-exported at package level.
    py_record = oecluster.ClusterRecord()
    native_record = _native.ClusterRecordVector(1)[0]

    for record in (py_record, native_record):
        assert math.isnan(record.mean_intra_distance)
        assert math.isnan(record.median_intra_distance)
        assert math.isnan(record.nearest_cluster_distance)
        assert math.isnan(record.silhouette)
        # -1 literal, not oecluster.NO_NEAREST_CLUSTER: the sentinel is a C++
        # constant that is not exported to the Python package namespace.
        assert record.nearest_cluster == -1
        # 0.0 is the real singleton value for these two, not a placeholder.
        assert record.radius == 0.0
        assert record.diameter == 0.0

    # The literals above say what seven of the values should be, which a
    # mirror cannot; this says the two sides agree on all twelve and in the
    # same order, which the literals cannot. repr rather than ==, because a
    # NaN does not compare equal to itself.
    assert [repr(value) for value in py_record] == [
        repr(getattr(native_record, name))
        for name in oecluster.ClusterRecord._fields
    ]


@pytest.mark.parametrize(
    "flag", ["compute_pair_rank_indices", "compute_per_cluster_records"]
)
def test_non_bool_flags_are_rejected_by_name(flag):
    """The match includes "must be True or False" deliberately. Before the
    signature gains the keyword, Python's own "unexpected keyword argument"
    TypeError also names the flag, so matching on the flag alone would pass
    green against an unimplemented feature."""
    import oecluster

    dm = _two_cluster_dm()
    result = oecluster.butina(dm, threshold=0.5)
    with pytest.raises(
            TypeError, match=f"{flag} must be True or False") as excinfo:
        oecluster.cluster_report(result, dm, **{flag: "yes"})
    # The offending value is carried as well as its type, which is what keeps
    # the message legible for a type whose own __name__ is "bool" -- numpy's
    # is, and a bare type name would read as "must be a bool, got bool".
    assert "'yes'" in str(excinfo.value)


@pytest.mark.parametrize(
    ("flag", "attribute"),
    [
        ("compute_pair_rank_indices", "pair_rank_indices"),
        ("compute_per_cluster_records", "per_cluster_records"),
    ],
)
def test_numpy_bools_are_accepted_and_take_effect(flag, attribute):
    """numpy.bool_ is what ``arr.any()`` and every comparison of numpy scalars
    returns, and allow_nonmetric on this same call has always taken it, so
    refusing it here gave one call two admissibility rules for its bool
    keywords. Acceptance alone is not the whole contract: the value has to
    reach the native option too, which is what requested reports back."""
    import numpy as np
    import oecluster

    dm = _two_cluster_dm()
    result = oecluster.butina(dm, threshold=0.5)

    on = oecluster.cluster_report(result, dm, **{flag: np.True_})
    assert getattr(on.requested, attribute) is True
    # The false half matters as much: a check that coerced with truthiness
    # instead of admitting the type would pass the True case and lose this one.
    off = oecluster.cluster_report(result, dm, **{flag: np.False_})
    assert getattr(off.requested, attribute) is False


@pytest.mark.parametrize(
    "flag", ["compute_pair_rank_indices", "compute_per_cluster_records"]
)
def test_flag_type_error_outranks_the_matrix_pairing_check(flag):
    """The local-argument block runs before the pairing check, so an
    authoritative complaint about what was typed is never pre-empted by an
    advisory one about the matrix. As above, the match pins the body's own
    message rather than the flag name, which an unexpected-keyword TypeError
    would also carry. Both flags are covered because the ordering is a
    property of each check's position, not of the pair."""
    import oecluster
    from oecluster import DenseStorage, SymmetricDistanceMatrix

    dm = _two_cluster_dm()
    result = oecluster.butina(dm, threshold=0.5)
    bigger = SymmetricDistanceMatrix(
        DenseStorage(6), "test", ["a", "b", "c", "d", "e", "f"], {})
    with pytest.raises(TypeError, match=f"{flag} must be True or False"):
        oecluster.cluster_report(result, bigger, **{flag: 1.5})


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
    # Both clusters hold two members 0.2 apart, 0.8 from the other cluster, so
    # every row carries the same values and any field sourced from a different
    # native field lands on a number this fixture does not produce.
    for ordinal, record in enumerate(detailed.records):
        assert record.label == ordinal
        assert record.size == 2
        # Written against the ordinal rather than as the constant it would be
        # for either row on its own: that is what tells nearest_cluster apart
        # from label, the two fields most easily crossed.
        assert record.nearest_cluster == 1 - ordinal
        # Either member of a two-member cluster can be the medoid, so pin what
        # is invariant -- that the representative is one of the sample indices
        # a medoid search can return here -- rather than which one it picked.
        assert record.representative in (1, 3)
        assert math.isclose(record.mean_intra_distance, 0.2, rel_tol=1e-9)
        assert math.isclose(record.median_intra_distance, 0.2, rel_tol=1e-9)
        assert math.isclose(record.radius, 0.2, rel_tol=1e-9)
        assert math.isclose(record.diameter, 0.2, rel_tol=1e-9)
        assert math.isclose(
            record.mean_representative_distance, 0.2, rel_tol=1e-9)
        assert math.isclose(record.nearest_cluster_distance, 0.8, rel_tol=1e-9)
        assert math.isclose(record.silhouette, 0.75, rel_tol=1e-9)
        assert record.boundary_violations == 0
    with pytest.raises(AttributeError):
        detailed.records[0].label = 7


def test_pair_rank_indices_are_computed_only_on_request():
    """The two pair-rank metrics are the only opt-in scalars, so whether the
    flag reaches the native options struct is the whole difference between a
    value and a NaN. The negative half is asserted as well, which is what
    makes this a wiring test rather than a value test."""
    import oecluster

    dm = _two_cluster_dm()
    result = oecluster.butina(dm, threshold=0.5)

    unasked = oecluster.cluster_report(result, dm)
    assert math.isnan(unasked.c_index)
    assert math.isnan(unasked.baker_hubert_gamma)

    asked = oecluster.cluster_report(
        result, dm, compute_pair_rank_indices=True)
    assert asked.requested.pair_rank_indices is True
    # Both are at their extremes because every within-cluster pair is closer
    # than every between-cluster pair on this fixture: C reaches its 0.0 floor
    # and gamma its 1.0 ceiling. An extreme is a weaker witness than an
    # interior value would be, but it still separates every wiring fault --
    # a dropped or crossed flag leaves both fields NaN.
    assert asked.c_index == 0.0
    assert asked.baker_hubert_gamma == 1.0


def test_noise_coverage_parallels_coverage():
    """Equal lengths are guaranteed by the C++ struct, so that assertion is
    true of the right vector and of any wrong one. The values are what tell
    the two apart: this clustering has no noise at all, so every noise entry
    must read NaN while ordinary coverage is complete."""
    import oecluster

    dm = _two_cluster_dm()
    report = oecluster.cluster_report(oecluster.butina(dm, threshold=0.5), dm)
    assert len(report.noise_coverage_at) == len(report.coverage_at)
    assert report.num_noise == 0
    assert all(math.isnan(value) for value in report.noise_coverage_at)
    # Pinning the ordinary curve is what makes the NaN assertion above
    # load-bearing: the two vectors genuinely differ here, rather than both
    # happening to be unset.
    assert report.coverage_at == (1.0, 1.0, 1.0)


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
