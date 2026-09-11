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


def _noise_bearing_dm():
    """Six points: a mutually close trio {0,1,2}, then 3, 4 and 5 far from all.

    Under DBSCAN at eps 0.2 with min_samples 3, points 0, 1 and 2 lie within
    eps of one another and each counts three neighbours including itself, so
    they are core points and form the single cluster; 3, 4 and 5 reach nobody
    and become noise. Point 0 is the cluster's medoid, its distances to the
    other two members summing to 0.2 against 0.3 for each of them, so every
    sample's distance to the nearest representative is simply its distance to
    point 0: 0.0, 0.1, 0.1, 0.3, 0.4 and 0.5.
    """
    from oecluster import DenseStorage, SymmetricDistanceMatrix

    storage = DenseStorage(6)
    for (left, right), distance in {
        (0, 1): 0.1, (0, 2): 0.1, (1, 2): 0.2,
        (0, 3): 0.3, (1, 3): 0.35, (2, 3): 0.35,
        (0, 4): 0.4, (1, 4): 0.45, (2, 4): 0.45,
        (0, 5): 0.5, (1, 5): 0.55, (2, 5): 0.55,
        (3, 4): 0.6, (3, 5): 0.6, (4, 5): 0.6,
    }.items():
        storage.Set(left, right, distance)
    return SymmetricDistanceMatrix(
        storage, "test", ["a", "b", "c", "d", "e", "f"], {})


def _row(table, name):
    """The value cells of the row named ``name``, without the name itself."""
    for row in table:
        if row[0] == name:
            return row[1:]
    raise AssertionError(f"no row named {name!r} in {[r[0] for r in table]}")


def _repr_row(rendered, name):
    """The rendered cells of the row named ``name``, without the name column."""
    for line in rendered.splitlines():
        if line.startswith(name) and not line[len(name):len(name) + 1].strip():
            return line[len(name):]
    raise AssertionError(f"no rendered row named {name!r} in:\n{rendered}")


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


def test_repr_column_labels_follow_the_report_order():
    """Which column a label sits over, not merely that the label is present.

    The check above looks for three method names anywhere in the header, which
    any permutation of them satisfies -- including a reversal, the arrangement
    that puts every number under someone else's name while the table still
    reads as a plausible comparison. Three clusterings that disagree on
    ``num_clusters`` make the pairing observable: the header and that row are
    split into cells and compared position by position.
    """
    import oecluster

    dm = _noise_bearing_dm()
    reports = (
        oecluster.cluster_report(oecluster.butina(dm, threshold=0.25), dm),
        oecluster.cluster_report(
            oecluster.dbscan(dm, eps=0.2, min_samples=3), dm),
        oecluster.cluster_report(oecluster.agglomerative(dm, n_clusters=3), dm),
    )
    assert [report.method for report in reports] == [
        "butina", "dbscan", "agglomerative"]
    # Distinct, and deliberately not symmetric under reversal: a palindrome
    # here would leave the label order unpinned by the row below.
    assert [report.num_clusters for report in reports] == [4, 1, 3]

    rendered = repr(oecluster.compare_reports(*reports))
    assert rendered.splitlines()[0].split() == [
        "metric", "butina", "dbscan", "agglomerative"]
    assert _repr_row(rendered, "num_clusters").split() == ["4", "1", "3"]


def test_repr_renders_a_float_cell_at_its_own_value():
    """The formatter is the last thing between a metric and a reader.

    Nothing else asserts a rendered number: the repr checks elsewhere look for
    a token such as ``nan`` or ``--`` in a named row, which any arithmetic on
    the finite cells leaves alone. A cell that renders some other number is
    the quietest failure this class of code has -- the table is well formed,
    the columns line up, and every value is wrong.
    """
    import oecluster

    dm = _two_cluster_dm()
    result = oecluster.butina(dm, threshold=0.5)
    reports = [oecluster.cluster_report(result, dm) for _ in range(2)]
    rendered = repr(oecluster.compare_reports(*reports))

    # Both clusters hold two members 0.2 apart and lie 0.8 from each other.
    assert _repr_row(rendered, "mean_intra_distance").split() == ["0.2", "0.2"]
    assert _repr_row(rendered, "silhouette").split() == ["0.75", "0.75"]
    # An integer field alongside, because the float branch is the only one the
    # rendering above exercises and a shift applied to both would otherwise
    # look like a single formatter's convention rather than an error.
    assert _repr_row(rendered, "num_samples").split() == ["4", "4"]


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
    are the only rows that can read None for "nobody asked" rather than NaN for
    "asked and undefined", and the comparison path is the one place
    ClusterReportRequested is otherwise unavailable."""
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
    # None rather than NaN in the second column: plain did not ask, so the cell
    # holds no answer at all. NaN there would be the undefined-value reading.
    assert values["c_index"][0] == 0.0
    assert values["c_index"][1] is None
    assert values["baker_hubert_gamma"][0] == 1.0
    assert values["baker_hubert_gamma"][1] is None

    # The claim is the gate's scope, not just its effect on the two rows above,
    # so it has to range over every scalar: widening the gate to cover a metric
    # that is always computed refuses an answer the report holds, which is as
    # wrong as failing to refuse one it does not. Only plain's column can show
    # this. Every cell of asked's is non-None whatever the gate covers, because
    # asked requested the one flag any of them could be gated on.
    for name in oecluster.ClusterReport._SCALAR_FIELDS:
        expected_none = name in ("c_index", "baker_hubert_gamma")
        assert (values[name][1] is None) is expected_none, name

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
    # None, not NaN: the report carries the threshold but answered nothing at
    # it, which is the same "never asked" state as a threshold it does not
    # carry -- and distinct from the NaN of a curve that was computed and came
    # out undefined.
    for threshold in clustered.coverage_thresholds:
        assert values[f"coverage_at[{threshold}]"][1] is None
        assert values[f"noise_coverage_at[{threshold}]"][1] is None


def test_comparison_table_carries_populated_coverage_curves():
    """The table looks its coverage values up by threshold index, and the two
    fixtures the tests above use are palindromes -- (1.0, 1.0, 1.0) beside three
    NaNs, and a column padded with None throughout -- so neither can tell an
    intact curve from a rearranged one. A report with real noise supplies a pair
    of curves that rise strictly and differ from each other at every threshold,
    which makes a rearrangement of either row family visible; pairing it with a
    report whose coverage is flat and whose noise curve is all NaN keeps the two
    columns from being confused with one another.
    """
    import oecluster

    noise_dm = _noise_bearing_dm()
    noisy = oecluster.cluster_report(
        oecluster.dbscan(noise_dm, eps=0.2, min_samples=3), noise_dm)
    clean_dm = _two_cluster_dm()
    clean = oecluster.cluster_report(
        oecluster.butina(clean_dm, threshold=0.5), clean_dm)
    # Both reports use the default preset, so the union of thresholds is that
    # preset's own and the row names below are the ones the table emits.
    assert noisy.coverage_thresholds == (0.25, 0.35, 0.45)
    assert clean.coverage_thresholds == (0.25, 0.35, 0.45)

    values = {row[0]: row[1:]
              for row in oecluster.compare_reports(noisy, clean).to_table()}

    # Asserted threshold by threshold rather than as whole tuples so a failure
    # names the threshold that moved.
    assert math.isclose(values["coverage_at[0.25]"][0], 3 / 6, rel_tol=1e-9)
    assert math.isclose(values["coverage_at[0.35]"][0], 4 / 6, rel_tol=1e-9)
    assert math.isclose(values["coverage_at[0.45]"][0], 5 / 6, rel_tol=1e-9)
    assert math.isclose(values["coverage_at[0.25]"][1], 1.0, rel_tol=1e-9)
    assert math.isclose(values["coverage_at[0.35]"][1], 1.0, rel_tol=1e-9)
    assert math.isclose(values["coverage_at[0.45]"][1], 1.0, rel_tol=1e-9)

    assert math.isclose(values["noise_coverage_at[0.25]"][0], 0.0, rel_tol=1e-9)
    assert math.isclose(
        values["noise_coverage_at[0.35]"][0], 1 / 3, rel_tol=1e-9)
    assert math.isclose(
        values["noise_coverage_at[0.45]"][0], 2 / 3, rel_tol=1e-9)
    # The clean report has no noise, so its whole noise curve is NaN. Compared
    # as reprs rather than through math.isnan so a failure prints the value.
    assert [repr(v) for v in values["noise_coverage_at[0.25]"][1:]] == ["nan"]
    assert [repr(v) for v in values["noise_coverage_at[0.35]"][1:]] == ["nan"]
    assert [repr(v) for v in values["noise_coverage_at[0.45]"][1:]] == ["nan"]


def test_comparison_scalar_cells_come_from_their_own_report():
    """Every scalar row against the three reports' own attributes.

    The table's whole purpose is to put one report's number beside another's,
    and nothing else in this file checks that a cell came from its own column's
    report: a lookup that read one fixed report for every column would publish
    one clustering's metrics under all the labels and satisfy every other test
    here. All three reports request the pair-rank indices, so no cell is gated
    and the oracle is uniform across all 27 scalars.

    The three clusterings are chosen so that no scalar holds the same value in
    all three columns -- a noisy DBSCAN run, a Butina split of the same matrix,
    and a Butina run tight enough to leave every point a singleton. A scalar the
    fixtures agreed on everywhere would be a scalar whose cell could be read
    from any report and still match, so the per-field check below is what stops
    the loop going quietly vacuous on the metrics that happen to coincide. It
    has to be per field: an aggregate count of how many scalars differ is
    satisfied by the ones that already discriminate. Non-constancy is as far as
    real clusterings reach, so cell-by-cell provenance -- every column
    distinguishable from every other on every field -- is pinned separately by
    ``test_comparison_scalar_cells_track_their_report_cell_by_cell``.
    """
    import oecluster

    noise_dm = _noise_bearing_dm()
    noisy = oecluster.cluster_report(
        oecluster.dbscan(noise_dm, eps=0.2, min_samples=3), noise_dm,
        compute_pair_rank_indices=True)
    split = oecluster.cluster_report(
        oecluster.butina(noise_dm, threshold=0.25), noise_dm,
        compute_pair_rank_indices=True)
    clean_dm = _two_cluster_dm()
    atomised = oecluster.cluster_report(
        oecluster.butina(clean_dm, threshold=0.1), clean_dm,
        compute_pair_rank_indices=True)

    reports = (noisy, split, atomised)
    table = oecluster.compare_reports(*reports).to_table()
    for name in oecluster.ClusterReport._SCALAR_FIELDS:
        # repr rather than == so a NaN cell compares equal to its own NaN.
        expected = [repr(getattr(report, name)) for report in reports]
        assert len(set(expected)) > 1, f"{name} is equal in every fixture"
        assert [repr(cell) for cell in _row(table, name)] == expected, name


def test_comparison_table_pads_each_report_against_its_own_thresholds():
    """Compared reports need not carry the same thresholds -- the row set is
    their union -- so each column has to be padded against its own list. Every
    other fixture in this file gives both reports identical thresholds, which
    leaves reading one report's list for every column indistinguishable from
    reading each report's own; overlapping-but-unequal lists separate the two,
    and a column padded against the wrong list reports one report's value under
    another's threshold while dropping the value it does have.

    Both columns come from the same clustering, so the threshold list is the
    only thing that differs between them. The expectations are derived from the
    per-sample distances to the nearest representative that the fixture's own
    docstring works out: coverage counts all six samples at or under each
    threshold, noise coverage counts only points 3, 4 and 5, and a threshold a
    report does not carry pads to None in that report's column alone.
    """
    import oecluster

    dm = _noise_bearing_dm()
    result = oecluster.dbscan(dm, eps=0.2, min_samples=3)
    low = oecluster.cluster_report(
        result, dm, coverage_thresholds=[0.15, 0.35])
    high = oecluster.cluster_report(
        result, dm, coverage_thresholds=[0.35, 0.45])
    # Asserted before any row is read so an override that stopped taking effect
    # fails legibly here rather than as a KeyError on a hardcoded row name.
    assert low.coverage_thresholds == (0.15, 0.35)
    assert high.coverage_thresholds == (0.35, 0.45)

    values = {row[0]: row[1:]
              for row in oecluster.compare_reports(low, high).to_table()}

    # Asserted threshold by threshold, both columns together, so a failure
    # names the cell that moved. The pads are None rather than NaN: neither
    # report was ever asked about the other's threshold, and both curves here
    # are fully populated, so a NaN in any of these cells would be wrong twice
    # over.
    assert math.isclose(values["coverage_at[0.15]"][0], 3 / 6, rel_tol=1e-9)
    assert values["coverage_at[0.15]"][1] is None
    assert math.isclose(values["coverage_at[0.35]"][0], 4 / 6, rel_tol=1e-9)
    assert math.isclose(values["coverage_at[0.35]"][1], 4 / 6, rel_tol=1e-9)
    assert values["coverage_at[0.45]"][0] is None
    assert math.isclose(values["coverage_at[0.45]"][1], 5 / 6, rel_tol=1e-9)

    assert math.isclose(
        values["noise_coverage_at[0.15]"][0], 0.0, rel_tol=1e-9)
    assert values["noise_coverage_at[0.15]"][1] is None
    assert math.isclose(
        values["noise_coverage_at[0.35]"][0], 1 / 3, rel_tol=1e-9)
    assert math.isclose(
        values["noise_coverage_at[0.35]"][1], 1 / 3, rel_tol=1e-9)
    assert values["noise_coverage_at[0.45]"][0] is None
    assert math.isclose(
        values["noise_coverage_at[0.45]"][1], 2 / 3, rel_tol=1e-9)


def test_pair_rank_cells_distinguish_unasked_from_undefined():
    """Both states of an opt-in cell in one place, which is what makes either
    one mean anything: a report that did not ask reads None, and a report that
    asked and got no answer reads NaN. A table that collapsed the two would
    satisfy half of this and fail the other. The repr is checked as well
    because the two states are useless to a reader if they print alike.

    The rendering is read row by row rather than swept for a token: a search of
    the whole table finds a -- wherever it comes from, so it passes just as well
    when the formatter blanks some unrelated cell it should have printed. Naming
    the row ties each token to the cell that must carry it."""
    import oecluster

    dm = _two_cluster_dm()
    result = oecluster.butina(dm, threshold=0.5)
    cheap = oecluster.cluster_report(result, dm)
    rich = oecluster.cluster_report(result, dm, compute_pair_rank_indices=True)

    cells = _row(oecluster.compare_reports(cheap, rich).to_table(), "c_index")
    assert cells[0] is None
    assert isinstance(cells[1], float)

    # One cluster of four: every distance is a within-pair, so S_min == S_max
    # and no couple is concordant or discordant. Asked, and undefined.
    single = oecluster.ClusteringResult([0, 0, 0, 0], [[0, 1, 2, 3]])
    undefined = oecluster.cluster_report(
        single, dm, compute_pair_rank_indices=True)
    comparison = oecluster.compare_reports(cheap, undefined)
    undefined_cells = _row(comparison.to_table(), "c_index")
    assert undefined_cells[0] is None
    assert math.isnan(undefined_cells[1])

    rendered = repr(comparison)
    c_index_line = _repr_row(rendered, "c_index")
    assert "--" in c_index_line
    assert "nan" in c_index_line
    # The other two rows are what keep the assertions above from being
    # satisfied by any -- anywhere in the table. A formatter that tested
    # falsiness rather than identity would blank the False flag and the zero
    # count, claiming nobody asked about values the reports did answer.
    assert "False" in _repr_row(rendered, "requested_pair_rank_indices")
    assert "0" in _repr_row(rendered, "num_noise")
    assert "--" not in _repr_row(rendered, "num_noise")


def test_coverage_rows_use_none_for_thresholds_a_report_never_used():
    """Two reports with disjoint threshold lists, so every row in the union is
    one report's own and the other's blank. The third assertion is what keeps
    the first two honest: a lookup that returned None for everything would
    satisfy them both while reporting nothing at all."""
    import oecluster

    dm = _two_cluster_dm()
    result = oecluster.butina(dm, threshold=0.5)
    a = oecluster.cluster_report(result, dm, coverage_thresholds=[0.3])
    b = oecluster.cluster_report(result, dm, coverage_thresholds=[0.7])

    table = oecluster.compare_reports(a, b).to_table()
    assert _row(table, "coverage_at[0.3]")[1] is None
    assert _row(table, "coverage_at[0.7]")[0] is None
    assert _row(table, "coverage_at[0.3]")[0] is not None


def test_coverage_lookup_indexes_a_report_s_own_threshold_order():
    """A report keeps the threshold order it was given; the rows do not.

    Row names come from the sorted union across reports, but a value is found
    by walking one report's own ``coverage_thresholds`` and taking the curve
    entry at the matching position. Sorting that walk would keep the row names
    identical and still read the wrong entry, and every other fixture here
    passes its thresholds already ascending, which makes the two indistinct.

    Passed descending, the report's own order is (0.45, 0.25, 0.35), so the
    curve entries are in that order too. The expectations are the fixture's:
    of the six samples, three lie within 0.25 of the nearest representative,
    four within 0.35 and five within 0.45, and of the three noise points none,
    one and two do.
    """
    import oecluster

    dm = _noise_bearing_dm()
    result = oecluster.dbscan(dm, eps=0.2, min_samples=3)
    descending = oecluster.cluster_report(
        result, dm, coverage_thresholds=[0.45, 0.25, 0.35])
    # The premise: the native preserves the caller's order rather than sorting
    # it, so the lookup has a wrong answer available to give.
    assert descending.coverage_thresholds == (0.45, 0.25, 0.35)

    ascending = oecluster.cluster_report(
        result, dm, coverage_thresholds=[0.25, 0.35, 0.45])
    values = {row[0]: row[1:]
              for row in oecluster.compare_reports(
                  descending, ascending).to_table()}

    # Both columns describe the same clustering at the same thresholds, so the
    # two cells of each row must agree however the report stored them.
    for threshold, coverage, noise_coverage in (
            (0.25, 3 / 6, 0.0), (0.35, 4 / 6, 1 / 3), (0.45, 5 / 6, 2 / 3)):
        for cell in values[f"coverage_at[{threshold}]"]:
            assert math.isclose(cell, coverage, rel_tol=1e-9), threshold
        for cell in values[f"noise_coverage_at[{threshold}]"]:
            assert math.isclose(cell, noise_coverage, rel_tol=1e-9), threshold


def test_noise_coverage_rows_are_nan_for_a_noise_free_report():
    """The preserved half of the distinction, and the reason a blanket sweep of
    NaN to None would be wrong. This report carries the threshold and its noise
    curve is full length, so the cell was asked and came back undefined; only a
    cell nobody asked about may read None."""
    import oecluster

    dm = _two_cluster_dm()
    report = oecluster.cluster_report(
        oecluster.butina(dm, threshold=0.5), dm, coverage_thresholds=[0.3])
    table = oecluster.compare_reports(report, report).to_table()
    assert math.isnan(_row(table, "noise_coverage_at[0.3]")[0])


def test_all_noise_report_contributes_threshold_rows_that_are_all_none():
    """ClusterReport copies coverage_thresholds whatever K is, and the row-key
    union reads that field rather than the curve, so an all-noise report's own
    threshold becomes a row whose every cell is None -- the honest reading of
    "no report answered here"."""
    import oecluster

    dm = _two_cluster_dm()
    clustered = oecluster.cluster_report(
        oecluster.butina(dm, threshold=0.5), dm, coverage_thresholds=[0.3])
    all_noise = oecluster.cluster_report(
        oecluster.ClusteringResult([-1, -1, -1, -1], []), dm,
        coverage_thresholds=[0.55])

    table = oecluster.compare_reports(clustered, all_noise).to_table()
    assert _row(table, "coverage_at[0.55]") == (None, None)
    assert _row(table, "coverage_at[0.3]")[1] is None
    assert _row(table, "coverage_at[0.3]")[0] is not None


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

    # The declared order mirrors the native struct's, and it is the order of
    # tuple(record) and of the columns pandas.DataFrame(report.records) builds,
    # so it is public. Pinned against a literal because every other assertion
    # on a record reads a field by name: a swap of two same-typed fields --
    # label, size and representative are all int, and radius, diameter and
    # mean_representative_distance are all float -- is invisible to all of
    # them, and to the mirror below, whose two sides are both driven by this
    # same declaration and so move together.
    assert oecluster.ClusterRecord._fields == (
        "label", "size", "representative", "mean_intra_distance",
        "median_intra_distance", "radius", "diameter",
        "mean_representative_distance", "nearest_cluster",
        "nearest_cluster_distance", "silhouette", "boundary_violations")

    # The literals above say what seven of the values should be, which a
    # mirror cannot; this says the two sides agree on all twelve values, with
    # the literal order above holding the order they are read in. repr rather
    # than ==, because a NaN does not compare equal to itself.
    assert [repr(value) for value in py_record] == [
        repr(getattr(native_record, name))
        for name in oecluster.ClusterRecord._fields
    ]


@pytest.mark.parametrize(
    "flag",
    ["compute_pair_rank_indices", "compute_per_cluster_records",
     "treat_noise_as_singletons"],
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


def test_treat_noise_as_singletons_admits_numpy_bools_and_keeps_both_effects():
    """The gate must not become an over-refusal. ``requested`` does not record
    this keyword, so the effect is read where it lands: on a DBSCAN result with
    three noise points, folding noise into the singleton accounting moves
    singleton_fraction and leaving it out does not. Both halves are pinned,
    because a check that coerced with truthiness rather than admitting the type
    would still pass the True case."""
    import numpy as np
    import oecluster

    dm = _noise_bearing_dm()
    result = oecluster.dbscan(dm, eps=0.2, min_samples=3)
    assert oecluster.cluster_report(
        result, dm, treat_noise_as_singletons=np.True_).singleton_fraction == 0.75
    assert oecluster.cluster_report(
        result, dm, treat_noise_as_singletons=np.False_).singleton_fraction == 0.0
    # The Python bools behave identically, which is what makes the numpy pair
    # above an admission rather than a second behaviour.
    assert oecluster.cluster_report(
        result, dm, treat_noise_as_singletons=True).singleton_fraction == 0.75
    assert oecluster.cluster_report(
        result, dm, treat_noise_as_singletons=False).singleton_fraction == 0.0


def test_treat_noise_as_singletons_is_refused_in_signature_order():
    """The local-argument block reports its keywords in the order the signature
    declares them, so a call that is wrong in two ways names the one the caller
    wrote first. treat_noise_as_singletons sits between representative_method
    and num_threads, and both neighbours are asserted: placing the new check
    beside its two bool siblings at the end of the block would have inverted
    the num_threads pair."""
    import oecluster

    dm = _two_cluster_dm()
    result = oecluster.butina(dm, threshold=0.5)

    with pytest.raises(ValueError, match="Unknown representative method"):
        oecluster.cluster_report(
            result, dm, representative_method="bogus",
            treat_noise_as_singletons="no")
    with pytest.raises(
            TypeError, match="treat_noise_as_singletons must be True or False"):
        oecluster.cluster_report(
            result, dm, treat_noise_as_singletons="no", num_threads=-1)


@pytest.mark.parametrize(
    "flag",
    ["compute_pair_rank_indices", "compute_per_cluster_records",
     "treat_noise_as_singletons"],
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
    # The declared order is the public tuple's order, and every other
    # assertion here reads a field by name -- which a swap of the two
    # declarations leaves untouched while reversing positional unpacking and
    # anything that serialises the tuple. The order is not arbitrary: it
    # mirrors the native struct, which declares pair_rank_indices first.
    assert oecluster.ClusterReportRequested._fields == (
        "pair_rank_indices", "per_cluster_records")

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


def test_record_fields_separate_on_an_unevenly_spaced_cluster():
    """Two members 0.2 apart give a cluster whose mean, median, radius,
    diameter and mean representative distance are all the same number, so a
    record field sourced from any of the other four still reads correctly.
    Four members at unequal distances pull the five apart -- 0.45, 0.50, 0.55,
    0.60 and 0.35 here -- and pairwise distinct values are what make every
    crossing among them fail.

    The matrix is built inline because the shared two-cluster helper cannot
    produce this shape. Every expectation is derived from the definitions:
    cluster 0's six within-pair distances are 0.15, 0.35, 0.45, 0.55, 0.60 and
    0.60, which sum to 2.70 for a mean of 0.45, whose middle two sorted entries
    average to a median of 0.50, and whose largest is the diameter, 0.60. The
    members' distance sums are 1.20, 1.40, 1.05 and 1.75, so sample 2 is the
    medoid; its mean to the other three is 1.05 / 3 = 0.35 and its farthest
    member lies at 0.55, the radius. Every member's nearest other cluster sits
    at 1.00, and the mean of the members' own-cluster mean distances is the
    mean within-pair distance itself, so the cluster silhouette is
    1 - 0.45 / 1.00 = 0.55.
    """
    import oecluster
    from oecluster import DenseStorage, SymmetricDistanceMatrix

    storage = DenseStorage(6)
    within = {
        (0, 1): 0.45, (0, 2): 0.15, (0, 3): 0.60,
        (1, 2): 0.35, (1, 3): 0.60, (2, 3): 0.55,
        (4, 5): 0.10,
    }
    for (left, right), distance in within.items():
        storage.Set(left, right, distance)
    # One value for every cross pair keeps each member's nearest-other-cluster
    # mean at 1.00, which is what makes the silhouette derivable by hand, and
    # puts every between-cluster distance beyond the 0.30 boundary threshold so
    # no violation is counted.
    for left in range(4):
        for right in (4, 5):
            storage.Set(left, right, 1.00)
    dm = SymmetricDistanceMatrix(
        storage, "test", ["a", "b", "c", "d", "e", "f"], {})

    # Partitioned by hand rather than by an algorithm: the point is the record
    # arithmetic, and a threshold that happened to split these four differently
    # would change the values without changing the test.
    result = oecluster.ClusteringResult(
        [0, 0, 0, 0, 1, 1], [[0, 1, 2, 3], [4, 5]])
    report = oecluster.cluster_report(
        result, dm, compute_per_cluster_records=True)
    uneven, pair = report.records

    assert uneven.label == 0
    assert uneven.size == 4
    assert uneven.representative == 2
    assert math.isclose(uneven.mean_intra_distance, 0.45, rel_tol=1e-9)
    assert math.isclose(uneven.median_intra_distance, 0.50, rel_tol=1e-9)
    assert math.isclose(uneven.radius, 0.55, rel_tol=1e-9)
    assert math.isclose(uneven.diameter, 0.60, rel_tol=1e-9)
    assert math.isclose(
        uneven.mean_representative_distance, 0.35, rel_tol=1e-9)
    assert uneven.nearest_cluster == 1
    assert math.isclose(uneven.nearest_cluster_distance, 1.00, rel_tol=1e-9)
    assert math.isclose(uneven.silhouette, 0.55, rel_tol=1e-9)
    assert uneven.boundary_violations == 0

    # The second cluster exists so the first has a neighbour at all. Its own
    # five floats coincide at 0.10, which is the collapsed case again.
    assert pair.label == 1
    assert pair.size == 2
    # Either member is an equally good medoid of a two-member cluster.
    assert pair.representative in (4, 5)
    assert math.isclose(pair.mean_intra_distance, 0.10, rel_tol=1e-9)
    assert pair.nearest_cluster == 0
    assert math.isclose(pair.nearest_cluster_distance, 1.00, rel_tol=1e-9)
    assert math.isclose(pair.silhouette, 0.90, rel_tol=1e-9)
    assert pair.boundary_violations == 0


def test_record_boundary_violations_are_the_cluster_s_own_not_the_total():
    """The one record field whose name is also a report field's name.

    Every other fixture in this file counts no violations at all, so a record
    column filled from ``ClusterReport.boundary_violations`` instead of the
    record's own copies a 0 onto a 0 and nothing notices. Three equal clusters
    all beneath one boundary separate them: each of the twelve cross pairs is
    counted once in the report and in both of its clusters, so the report reads
    12 while each record reads 8, and the record column sums to twice the
    report -- the relation the record field documents.

    Equal clusters are deliberate. They make the three records identical, so
    the failure this catches cannot be mistaken for a misrouting between rows;
    what is being read is the field, not the row.
    """
    import oecluster
    from oecluster import DenseStorage, SymmetricDistanceMatrix

    storage = DenseStorage(6)
    for left in range(6):
        for right in range(left + 1, 6):
            # Paired by index: {0,1}, {2,3}, {4,5}.
            storage.Set(left, right,
                        0.2 if left // 2 == right // 2 else 0.8)
    dm = SymmetricDistanceMatrix(
        storage, "test", ["a", "b", "c", "d", "e", "f"], {})

    result = oecluster.butina(dm, threshold=0.5)
    report = oecluster.cluster_report(
        result, dm, boundary_threshold=0.9, compute_per_cluster_records=True)

    assert report.num_clusters == 3
    assert report.boundary_violations == 12
    per_cluster = [record.boundary_violations for record in report.records]
    assert per_cluster == [8, 8, 8]
    # Stated as the relation rather than left implicit in the numbers above, so
    # the two quantities cannot drift into agreement without this failing.
    assert sum(per_cluster) == 2 * report.boundary_violations


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
    # Read positionally as well, on the suite's only asymmetric request: the
    # field names can be declared correctly while the tuple they produce comes
    # out reversed, and only an unequal pair of flags can tell the two apart.
    assert tuple(asked.requested) == (True, False)
    # Both are at their extremes because every within-cluster pair is closer
    # than every between-cluster pair on this fixture: C reaches its 0.0 floor
    # and gamma its 1.0 ceiling. An extreme is a weaker witness than an
    # interior value would be, but it still separates every wiring fault --
    # a dropped or crossed flag leaves both fields NaN.
    assert asked.c_index == 0.0
    assert asked.baker_hubert_gamma == 1.0


def test_requested_holds_for_a_pair_rank_index_that_came_back_undefined():
    """requested records the ask, not the outcome, and this is the case that
    separates the two. A single cluster has no between-cluster pair for either
    index to rank, so both come back NaN however the caller asked -- and a flag
    derived from whether a value arrived would report False here, collapsing
    "asked and undefined" into "nobody asked". The NaN assertions are what make
    the True load-bearing; without them the report is indistinguishable from
    any other successful request.

    Built inline because the shared two-cluster helper has two clusters by
    construction and so cannot leave a pair-rank index undefined.
    """
    import oecluster
    from oecluster import DenseStorage, SymmetricDistanceMatrix

    storage = DenseStorage(3)
    storage.Set(0, 1, 0.2)
    storage.Set(0, 2, 0.3)
    storage.Set(1, 2, 0.25)
    dm = SymmetricDistanceMatrix(storage, "test", ["a", "b", "c"], {})

    # Every pair lies inside the threshold, so all three points join one
    # cluster and the between-pair array both indices need is empty.
    report = oecluster.cluster_report(
        oecluster.butina(dm, threshold=0.5), dm,
        compute_pair_rank_indices=True)

    assert report.num_clusters == 1
    assert math.isnan(report.c_index)
    assert math.isnan(report.baker_hubert_gamma)
    assert report.requested.pair_rank_indices is True


def test_requested_holds_when_a_records_request_yields_no_records():
    """The same rule on the other flag. An all-noise clustering has no cluster
    to describe, so records is empty whatever was asked for, and a flag sourced
    from the tuple's emptiness would report False for a request that was made.
    """
    import oecluster
    from oecluster import DenseStorage, SymmetricDistanceMatrix

    # Four mutually distant points: no neighbourhood at this eps holds anyone
    # but the point itself, so nothing reaches core status and all four are
    # noise.
    storage = DenseStorage(4)
    for left in range(4):
        for right in range(left + 1, 4):
            storage.Set(left, right, 0.9)
    dm = SymmetricDistanceMatrix(storage, "test", ["a", "b", "c", "d"], {})

    report = oecluster.cluster_report(
        oecluster.dbscan(dm, eps=0.1, min_samples=3), dm,
        compute_per_cluster_records=True)

    assert report.num_clusters == 0
    assert report.num_noise == 4
    assert report.records == ()
    assert report.requested.per_cluster_records is True


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
    # Compared as a materialised list of reprs rather than through all(), which
    # reports only that a generator was falsy and hides which entry was wrong.
    assert [repr(value) for value in report.noise_coverage_at] == (
        ["nan"] * len(report.coverage_at))
    # Pinning the ordinary curve is what makes the NaN assertion above
    # load-bearing: the two vectors genuinely differ here, rather than both
    # happening to be unset.
    assert report.coverage_at == (1.0, 1.0, 1.0)


def test_noise_coverage_climbs_with_distance_from_the_representative():
    """The NaN case above pins a vector empty of information: three NaNs are
    three NaNs however they are reordered, so no rearrangement of the curve is
    visible there, and the ordinary curve on that fixture is (1.0, 1.0, 1.0),
    its own reverse. A populated curve is needed to see one. Here both curves
    rise strictly and differ from each other at every entry, which makes a
    reversal of either visible and also makes sourcing either from the other
    visible.

    Derived from the definitions rather than read off a run. Given the
    per-sample distances to the nearest representative that the fixture's own
    docstring works out, coverage counts the samples at or under each threshold
    over all six -- three, then four, then five -- while noise coverage counts
    only points 3, 4 and 5, of which none, then one, then two are covered.
    """
    import oecluster

    dm = _noise_bearing_dm()

    report = oecluster.cluster_report(
        oecluster.dbscan(dm, eps=0.2, min_samples=3), dm)

    assert report.num_clusters == 1
    assert report.num_noise == 3
    assert report.coverage_thresholds == (0.25, 0.35, 0.45)
    # Asserted entry by entry rather than as a whole tuple so a failure names
    # the threshold that moved.
    assert math.isclose(report.coverage_at[0], 3 / 6, rel_tol=1e-9)
    assert math.isclose(report.coverage_at[1], 4 / 6, rel_tol=1e-9)
    assert math.isclose(report.coverage_at[2], 5 / 6, rel_tol=1e-9)
    assert math.isclose(report.noise_coverage_at[0], 0.0, rel_tol=1e-9)
    assert math.isclose(report.noise_coverage_at[1], 1 / 3, rel_tol=1e-9)
    assert math.isclose(report.noise_coverage_at[2], 2 / 3, rel_tol=1e-9)


def test_partition_error_surfaces_as_runtime_error():
    """Every std::exception maps to SWIG_RuntimeError; A1 adds no ValueError
    mapping. The message must survive so the caller learns which index is
    wrong. The result must still cover four samples, or the result/matrix
    pairing check fires first and this tests nothing."""
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


def test_nan_coverage_threshold_is_refused():
    """A NaN threshold reads as an ordinary number that simply never matches.

    It passes the non-negative check, because ``nan < 0.0`` is false, and the
    native then answers it: the report holds a real coverage value. But the
    comparison table finds a report's value by ``t == threshold``, which NaN
    never satisfies, so the cell publishes as None -- and None in that table
    means the question was never asked. The report would hold an answer while
    the table said nobody asked for one, so the call has to be refused here,
    at the last point where the caller can still be told.

    Matched on the whole message rather than the shared tail: two adjacent
    checks whose messages differ only in the subject are exactly the pair a
    reader has to tell apart, and a match on "must not be NaN" alone accepts
    either one -- so a caller who passed a NaN coverage threshold could be sent
    to inspect boundary_threshold with the suite still green.
    """
    import oecluster

    dm = _two_cluster_dm()
    result = oecluster.butina(dm, threshold=0.5)
    with pytest.raises(ValueError,
                       match="coverage thresholds must not be NaN"):
        oecluster.cluster_report(result, dm, coverage_thresholds=[math.nan])


def test_nan_boundary_threshold_is_refused():
    """The same hole on the boundary threshold, with a worse consequence.

    No distance compares true against NaN, so every point clears the boundary
    and the report comes back with ``boundary_violations = 0``: a clean bill of
    health for a question that was never answerable. Unlike the coverage case
    there is not even a None to hint at it, so nothing downstream can recover
    the fact that the threshold was meaningless.

    Matched on the subject as well as the rule, for the reason given on the
    coverage test above: the two messages are each other's nearest neighbour,
    and only the subject distinguishes them.
    """
    import oecluster

    dm = _two_cluster_dm()
    result = oecluster.butina(dm, threshold=0.5)
    with pytest.raises(ValueError,
                       match="boundary_threshold must not be NaN"):
        oecluster.cluster_report(result, dm, boundary_threshold=math.nan)


def test_negative_infinity_still_refused_as_non_negative():
    """Negative and NaN are disjoint, and the messages must stay disjoint too.

    ``-inf`` is refused for being negative, not for being non-finite, and a
    caller reading the message is being told which rule it broke. Matching only
    on ValueError would pass whichever check fired, so this pins the wording:
    if the NaN check ever widened to cover the negatives, or took their message
    over, this is what notices.
    """
    import oecluster

    dm = _two_cluster_dm()
    result = oecluster.butina(dm, threshold=0.5)
    with pytest.raises(ValueError, match="must be non-negative"):
        oecluster.cluster_report(result, dm,
                                 coverage_thresholds=[float("-inf")])


def test_infinite_coverage_threshold_is_accepted_and_reaches_the_table():
    """Over-refusal is the mirror defect, and ``inf`` is a coherent question.

    Everything is within an infinite distance, so coverage at ``inf`` is 1.0,
    and ``inf == inf`` holds, so the comparison table's threshold match finds
    it and the cell carries a real number. A NaN check written as a blanket
    non-finite ban would refuse this call, and one written to refuse further
    downstream would leave the cell as None -- rendered ``--``, meaning nobody
    asked. Asserting only that the call did not raise would miss the second.
    """
    import oecluster

    dm = _two_cluster_dm()
    result = oecluster.butina(dm, threshold=0.5)
    reports = [
        oecluster.cluster_report(result, dm, coverage_thresholds=[float("inf")])
        for _ in range(2)
    ]
    for report in reports:
        assert report.coverage_thresholds == (float("inf"),)

    table = oecluster.compare_reports(*reports).to_table()
    cells = _row(table, "coverage_at[inf]")
    assert None not in cells
    assert cells == (1.0, 1.0)


def test_infinite_boundary_threshold_is_accepted_and_counts_every_pair():
    """The same mirror defect on the other threshold, where it is easier to hit.

    ``inf`` is refused by ``not math.isfinite`` and admitted by ``math.isnan``,
    and the two read alike on every NaN, so only a finite non-finite value can
    tell the intended check from the over-broad one. The question is coherent:
    every between-cluster distance lies within an infinite boundary, so the
    answer is every cross pair. Asserting only that the call returned would
    pass just as well if the threshold had been dropped on the way through and
    the count came back 0, so the count itself is what this pins.

    Four points in two pairs give four cross pairs, each counted once in the
    report total, and each counted in both of its clusters -- so the two
    records read 4 apiece and their sum is twice the total.
    """
    import oecluster

    dm = _two_cluster_dm()
    result = oecluster.butina(dm, threshold=0.5)
    report = oecluster.cluster_report(
        result, dm, boundary_threshold=float("inf"),
        compute_per_cluster_records=True)

    assert report.boundary_violations == 4
    assert [record.boundary_violations for record in report.records] == [4, 4]

    # The counterpart under the default preset's finite boundary of 0.30: the
    # cross distances are 0.8, so nothing violates. Without it, a threshold
    # ignored entirely and a threshold of inf would be told apart only if the
    # ignored default happened to count something, which here it does not.
    default_boundary = oecluster.cluster_report(
        result, dm, compute_per_cluster_records=True)
    assert default_boundary.boundary_violations == 0


def _sentinel_report(index):
    """A report whose every scalar is a value unique to it and to that field.

    Built from a stand-in for the native report rather than from a clustering:
    no real clustering can be relied on to give 27 metrics that all differ
    across every pair of columns, and several are legitimately NaN in more than
    one column at once. The wrapper reads the native report by plain attribute
    access, so a namespace carrying the same names is enough to construct a
    genuine report whose cells are individually identifiable.
    """
    import types

    import oecluster

    native = types.SimpleNamespace(
        coverage_thresholds=(0.3,), coverage_at=(0.5,),
        noise_coverage_at=(0.25,), records=(),
        requested=types.SimpleNamespace(pair_rank_indices=True,
                                        per_cluster_records=False))
    for position, name in enumerate(oecluster.ClusterReport._SCALAR_FIELDS):
        setattr(native, name, 100.0 * index + position)
    return oecluster.ClusterReport(native, method=f"stub{index}")


def test_comparison_scalar_cells_track_their_report_cell_by_cell():
    """Each of the 81 scalar cells against the one value only its report holds.

    The companion test over real clusterings can only tell a cell from its
    neighbour where the two clusterings disagree on that metric, so a
    misrouting that happens to land on a column holding the same number stays
    invisible there. Here every report/field pair carries a different value, so
    a cell read from the wrong report is wrong whichever report it came from
    and whichever field it was.
    """
    import oecluster

    fields = oecluster.ClusterReport._SCALAR_FIELDS
    reports = tuple(_sentinel_report(index) for index in range(3))
    # The oracle is only as strong as the sentinels are distinct. Were any two
    # of the cells to hold the same value, a misroute between them would match
    # and this test would weaken without ever failing -- the failure mode it
    # exists to rule out, reappearing one level up in its own fixture.
    stored = [getattr(report, name) for report in reports for name in fields]
    assert len(set(stored)) == len(stored)

    table = oecluster.compare_reports(*reports).to_table()
    for position, name in enumerate(fields):
        expected = tuple(100.0 * index + position for index in range(3))
        assert _row(table, name) == expected, name
