"""Tests for the SAR coherence metrics."""

import concurrent.futures
import math

import oecluster
import pytest
from oecluster import (
    ClusteringResult,
    DenseStorage,
    SparseStorage,
    SymmetricDistanceMatrix,
)


def _line_dm(coordinates, facts=None):
    """A distance matrix over a 1-D embedding, with distances ``|xi - xj|``.

    Absolute differences along a line are a metric by construction, so these
    fixtures carry exact literal distances without depending on what a
    fingerprint happens to score for a given SMILES.
    """
    n = len(coordinates)
    storage = DenseStorage(n)
    for i in range(n):
        for j in range(i + 1, n):
            storage.Set(i, j, abs(coordinates[i] - coordinates[j]))
    labels = [f"m{i}" for i in range(n)]
    return SymmetricDistanceMatrix(storage, "test", labels, {}, facts)


# Three points 0.25 apart end to end: d(0,1) = d(1,2) = 0.25, d(0,2) = 0.5.
_SALI_COORDS = [0.0, 0.25, 0.5]
_SALI_ACTIVITY = [0.0, 1.0, 3.0]

# Four points with classes A, A, B, A. Each of 0 and 1 has the other as its
# nearest neighbour; 2 and 3 have each other and disagree.
_MODI_COORDS = [0.0, 0.1, 0.5, 0.7]
_MODI_CLASSES = ["A", "A", "B", "A"]

# Three clusters of two, activity rising by two between neighbours.
_COHERENCE_LABELS = [0, 0, 1, 1, 2, 2]
_COHERENCE_ACTIVITY = [1.0, 3.0, 5.0, 7.0, 9.0, 11.0]


def test_sar_coherence_decomposes_a_three_cluster_fixture():
    coherence = oecluster.sar_coherence(_COHERENCE_LABELS, _COHERENCE_ACTIVITY)

    assert coherence.num_samples == 6
    assert coherence.num_scored == 6
    assert coherence.num_clusters == 3
    assert coherence.eta_squared == pytest.approx(0.9142857142857143, abs=1e-12)
    assert coherence.omega_squared == pytest.approx(0.8333333333333334,
                                                    abs=1e-12)


def test_sar_coherence_reports_the_per_cluster_table():
    coherence = oecluster.sar_coherence(_COHERENCE_LABELS, _COHERENCE_ACTIVITY)

    assert isinstance(coherence.clusters, tuple)
    assert [row.label for row in coherence.clusters] == [0, 1, 2]
    assert [row.num_scored for row in coherence.clusters] == [2, 2, 2]
    assert [row.mean_activity for row in coherence.clusters] == [2.0, 6.0, 10.0]
    # Population standard deviation of {1, 3} and of each other pair.
    assert [row.stddev_activity for row in coherence.clusters] == [1.0, 1.0, 1.0]
    assert isinstance(coherence.clusters[0], oecluster.ClusterActivity)


def test_sar_coherence_accepts_a_clustering_result(monkeypatch):
    """A ClusteringResult takes the result overload, and the answers agree.

    Agreement alone does not pin the dispatch. ``_agreement_labels`` opens with
    ``getattr(value, "labels", value)``, so a ClusteringResult handed to the
    label branch is decomposed and scored perfectly well -- deleting the
    ``isinstance`` branch from ``sar_coherence`` leaves every assertion below
    the spy passing. The spy is the only thing here that says the native
    overload taking a result was the one called.

    That the two *native* overloads agree on a shared labeling is pinned in
    ``tests/python/test_native_bindings.py``; what this test pins is the
    Python-side choice between them.
    """
    calls = []
    real_conversion = oecluster._native_clustering_result

    def spy(result):
        calls.append(result)
        return real_conversion(result)

    monkeypatch.setattr(oecluster, "_native_clustering_result", spy)

    dm = _line_dm([0.0, 0.1, 0.2, 5.0, 5.1, 5.2])
    result = oecluster.dbscan(dm, 0.3, min_samples=2)
    activity = [1.0, 1.2, 0.9, 7.0, 7.4, 7.1]

    from_result = oecluster.sar_coherence(result, activity)
    from_labels = oecluster.sar_coherence(list(result.labels), activity)

    # Once, for the result call. The label call must not reach the conversion,
    # or the two branches are not the two branches.
    assert len(calls) == 1
    assert calls[0] is result

    # Two degenerate answers agree as readily as two correct ones: had dbscan
    # returned one cluster, or all noise, the three equalities below would hold
    # without either overload having decomposed anything. It splits the fixture
    # in two, and the split explains almost all of the variance.
    assert from_result.num_clusters == 2
    assert from_result.eta_squared == pytest.approx(0.9976426214049976,
                                                    abs=1e-12)

    assert from_result.num_clusters == from_labels.num_clusters
    assert from_result.eta_squared == from_labels.eta_squared
    assert from_result.omega_squared == from_labels.omega_squared


def test_sar_coherence_excludes_noise_by_default():
    coherence = oecluster.sar_coherence([-1, 0, 0, 1, 1],
                                        [9.0, 1.0, 3.0, 5.0, 7.0])

    assert coherence.num_samples == 5
    assert coherence.num_scored == 4
    assert coherence.num_clusters == 2


def test_sar_coherence_honours_every_noise_spelling():
    """The three spellings decompose three different partitions.

    Two noise samples, not one: with a single noise sample, grouped and
    singletons produce the identical partition, and every count agrees between
    them no matter what the implementation does with the spelling. The second
    noise sample is what separates "pool the noise into one cluster" from "give
    each noise sample its own", so an implementation that collapsed singletons
    onto grouped would have to change a number here.

    The row tables are asserted as well as the counts, because the counts alone
    do not show that grouped pools the two noise samples into a single row whose
    mean is their average.
    """
    labels = [-1, -1, 0, 0, 1, 1]
    activity = [9.0, 10.0, 1.0, 3.0, 5.0, 7.0]

    excluded = oecluster.sar_coherence(labels, activity, noise="excluded")
    grouped = oecluster.sar_coherence(labels, activity, noise="grouped")
    singletons = oecluster.sar_coherence(labels, activity, noise="singletons")

    assert excluded.num_clusters == 2
    assert grouped.num_clusters == 3
    assert singletons.num_clusters == 4
    assert excluded.num_scored == 4
    assert grouped.num_scored == 6
    assert singletons.num_scored == 6

    def rows(coherence):
        return [(row.label, row.num_scored, row.mean_activity)
                for row in coherence.clusters]

    assert rows(excluded) == [(0, 2, 2.0), (1, 2, 6.0)]
    # One noise row, holding both samples, at their mean.
    assert rows(grouped) == [(-1, 2, 9.5), (0, 2, 2.0), (1, 2, 6.0)]
    # Two noise rows, each its own sample, each at its own value -- and both
    # still labelled -1, which is why the rows are read by position.
    assert rows(singletons) == [(-1, 1, 9.0), (-1, 1, 10.0), (0, 2, 2.0),
                                (1, 2, 6.0)]

    # The fourth field, which the row tables above leave out: a spread over one
    # sample is not defined, so each singleton noise row reports NaN while every
    # two-member row reports a real number.
    assert [row.stddev_activity for row in excluded.clusters] == [1.0, 1.0]
    assert [row.stddev_activity for row in grouped.clusters] == [0.5, 1.0, 1.0]
    singleton_stddevs = [row.stddev_activity for row in singletons.clusters]
    assert [math.isnan(value) for value in singleton_stddevs] == [True, True,
                                                                  False, False]
    assert singleton_stddevs[2:] == [1.0, 1.0]


def test_sar_coherence_rejects_an_unknown_noise_spelling():
    with pytest.raises(ValueError, match="Unknown noise handling"):
        oecluster.sar_coherence(_COHERENCE_LABELS, _COHERENCE_ACTIVITY,
                                noise="drop")


def test_sar_coherence_treats_nan_activity_as_missing():
    coherence = oecluster.sar_coherence([0, 0, 1, 1],
                                        [1.0, float("nan"), 5.0, 7.0])

    assert coherence.num_samples == 4
    assert coherence.num_scored == 3
    assert coherence.clusters[0].num_scored == 1
    # One scored member, so no spread is defined for that cluster.
    assert math.isnan(coherence.clusters[0].stddev_activity)


def test_sar_coherence_rejects_a_length_mismatch():
    with pytest.raises(ValueError, match="4 samples"):
        oecluster.sar_coherence([0, 0, 1, 1], [1.0, 2.0])


def test_sar_coherence_reports_an_empty_labeling_as_a_length_mismatch():
    """An empty labeling is a length mismatch, and both spellings say so alike.

    By the time the overload split is reached the activity is already known to
    be non-empty, so there is no empty labeling that is not also a length
    mismatch. A separate emptiness refusal on the sequence branch -- which the
    ClusteringResult branch has no equivalent of -- would give the same input
    shape two different messages depending on how the caller spelled the empty
    clustering.
    """
    with pytest.raises(ValueError) as from_sequence:
        oecluster.sar_coherence([], [1.0])

    with pytest.raises(ValueError) as from_result:
        oecluster.sar_coherence(ClusteringResult([], []), [1.0])

    assert str(from_sequence.value) == (
        "activity has 1 entries but the clustering has 0 samples")
    assert str(from_result.value) == str(from_sequence.value)


def test_sar_coherence_rejects_an_empty_clustering_result():
    """The length check alone cannot catch this: 0 == 0 agrees.

    Without the emptiness check ahead of the overload split, the refusal
    happens in C++ and arrives as RuntimeError, which is not what the
    docstring promises for an empty activity.
    """
    with pytest.raises(ValueError, match="non-empty activity"):
        oecluster.sar_coherence(ClusteringResult([], []), [])


def test_sar_coherence_rejects_an_out_of_range_result_label():
    """Both overloads report an unrepresentable label the same way.

    ClusteringResult keeps its labels in an intp array and validates no
    range, so the failure only happens while filling the native int vector.
    """
    result = ClusteringResult([0, 2 ** 40], [[0], [1]])

    with pytest.raises(ValueError, match="32-bit signed int"):
        oecluster.sar_coherence(result, [1.0, 2.0])


def test_sar_coherence_rejects_a_mapping_activity():
    """A Mapping iterates its keys, so a dict would score the indices."""
    with pytest.raises(TypeError, match="not a mapping"):
        oecluster.sar_coherence([0, 0], {0: 1.0, 1: 2.0})


def test_sar_coherence_rejects_a_string_activity_value():
    with pytest.raises(TypeError, match="sequence of floats"):
        oecluster.sar_coherence([0, 0], ["1.0", "2.0"])


def test_the_activity_entry_points_reject_a_bytes_like_activity():
    """b"12" iterates as ints, so without the guard it scores as [49.0, 50.0].

    The bare-str guard beside it catches the text form of the same mistake;
    modelability already refuses bytes through its per-item str requirement.
    """
    for bad in (b"12", bytearray(b"12"), memoryview(b"12")):
        with pytest.raises(TypeError, match="not a bytes-like object"):
            oecluster.sar_coherence([0, 0], bad)
    with pytest.raises(TypeError, match="not a bytes-like object"):
        oecluster.activity_landscape(_line_dm(_SALI_COORDS), b"123")


def test_sar_coherence_accepts_any_iterable_of_activity_values():
    tuple_form = oecluster.sar_coherence(
        _COHERENCE_LABELS, tuple(_COHERENCE_ACTIVITY))
    generator_form = oecluster.sar_coherence(
        _COHERENCE_LABELS, (value for value in _COHERENCE_ACTIVITY))

    for coherence in (tuple_form, generator_form):
        assert coherence.eta_squared == pytest.approx(0.9142857142857143,
                                                      abs=1e-12)


def test_sar_coherence_rejects_a_misspelled_keyword():
    """A typo'd keyword must refuse rather than tune nothing.

    The native options proxy accepts any attribute name and ignores it, so a
    misspelling that reached it would be inert. Keeping the tuning in keyword
    arguments puts the refusal in Python, and the control below shows the
    refusal is about the spelling rather than the option being unsupported.
    """
    assert oecluster.sar_coherence(_COHERENCE_LABELS, _COHERENCE_ACTIVITY,
                                   noise="grouped").num_clusters == 3

    # The misspelling is the subject of the test, so pyright flagging it is the
    # static half of the same finding rather than a defect to fix.
    with pytest.raises(TypeError, match="noize"):
        oecluster.sar_coherence(
            _COHERENCE_LABELS, _COHERENCE_ACTIVITY,
            noize="grouped")  # pyright: ignore[reportCallIssue]


def test_sar_coherence_to_table_and_repr():
    coherence = oecluster.sar_coherence(_COHERENCE_LABELS, _COHERENCE_ACTIVITY)
    table = coherence.to_table()

    assert [name for name, _ in table] == [
        "num_samples", "num_scored", "num_clusters",
        "eta_squared", "omega_squared",
    ]
    # The names alone are not the check: a to_table() emitting the right five
    # labels against empty cells would satisfy the list above and both
    # substring assertions below.
    assert dict(table)["num_clusters"] == 3
    assert dict(table)["eta_squared"] == pytest.approx(0.9142857142857143,
                                                       abs=1e-12)

    rendered = repr(coherence)
    assert rendered.splitlines()[0].startswith("metric")
    assert "eta_squared" in rendered
    assert "0.9143" in rendered
    # The per-cluster table is a sequence, not a metric row: a table of tables
    # does not render. Read off the row names rather than searched for as a
    # substring of the whole rendering, because "num_clusters" ends in
    # "clusters" and no implementation could satisfy the substring form.
    assert [line.split()[0] for line in rendered.splitlines()[1:]] == [
        "num_samples", "num_scored", "num_clusters",
        "eta_squared", "omega_squared",
    ]


def test_the_row_tables_survive_the_native_result():
    """Rows read off a temporary must still carry their values.

    The native result owns its member vectors, so reading one off a temporary
    parent yields an empty vector: ``native.sar_coherence(...).clusters`` is
    length zero even though the held parent reports two rows. The Pythonic
    layer therefore copies every row into its own record while the native
    result is still bound to a local name. Written in the temporary idiom on
    purpose -- indexing the rows straight off the returned scorecard -- so a
    future refactor that starts handing back the native vector fails here with
    an IndexError instead of shipping silent truncation.
    """
    cluster_row = oecluster.sar_coherence(_COHERENCE_LABELS,
                                          _COHERENCE_ACTIVITY).clusters[1]
    assert cluster_row == oecluster.ClusterActivity(
        label=1, num_scored=2, mean_activity=6.0, stddev_activity=1.0)

    class_row = oecluster.modelability(_line_dm(_MODI_COORDS),
                                       _MODI_CLASSES).classes[1]
    assert class_row == oecluster.ClassConcordance(
        label="B", num_members=1, fraction_same_class=0.0)


def test_activity_landscape_matches_a_hand_fixture():
    """SALI is |da| / d: 1/0.25, 3/0.5 and 2/0.25 -- 4, 6 and 8.

    Two pairs sit within the 0.30 distance threshold with at least one log
    unit between them, so two of the three pairs are cliffs. Every pair's
    activity difference exceeds the 0.625-sigma band, so no molecule's nearest
    in-band neighbour is closer than its nearest out-of-band one and rmodi is
    zero.
    """
    landscape = oecluster.activity_landscape(_line_dm(_SALI_COORDS),
                                             _SALI_ACTIVITY)

    assert landscape.num_samples == 3
    assert landscape.num_scored == 3
    assert landscape.num_pairs_scored == 3
    assert landscape.num_zero_distance_pairs == 0
    assert landscape.num_cliffs == 2
    assert landscape.cliff_density == pytest.approx(2.0 / 3.0, abs=1e-12)
    assert landscape.max_sali == 8.0
    assert landscape.mean_sali == 6.0
    assert landscape.rmodi == 0.0
    assert landscape.activity_stddev == pytest.approx(1.2472191289246473,
                                                      abs=1e-12)


def test_activity_landscape_counts_zero_distance_pairs_apart_from_sali():
    """A zero-distance pair has no SALI, but it is still a cliff.

    The ratio is undefined there, so the pair is counted and left out of the
    SALI statistics rather than contributing an infinity that would poison
    both the maximum and the mean.
    """
    landscape = oecluster.activity_landscape(_line_dm([0.0, 0.0, 0.5]),
                                             [0.0, 1.0, 2.0])

    assert landscape.num_pairs_scored == 3
    assert landscape.num_zero_distance_pairs == 1
    assert landscape.num_cliffs == 1
    assert landscape.max_sali == 4.0
    assert landscape.mean_sali == 3.0


def test_activity_landscape_reports_rmodi_one_for_flat_activity():
    """Flat activity has a zero band, and every pair sits inside it."""
    landscape = oecluster.activity_landscape(_line_dm(_SALI_COORDS),
                                             [2.0, 2.0, 2.0])

    assert landscape.activity_stddev == 0.0
    assert landscape.rmodi == 1.0
    assert landscape.num_cliffs == 0
    assert landscape.max_sali == 0.0

    # These three values are not what the geometry always says. The same three
    # points with a rising activity score the opposite on every one of them,
    # so the assertions above are reading the activity and not a constant.
    rising = oecluster.activity_landscape(_line_dm(_SALI_COORDS),
                                          _SALI_ACTIVITY)
    assert rising.rmodi == 0.0
    assert rising.num_cliffs == 2
    assert rising.max_sali == 8.0


def test_activity_landscape_thresholds_move_the_cliff_count():
    """Both thresholds are read, and each moves the count off its default.

    The default calls two of the three pairs cliffs. Widening the distance
    threshold past 0.5 admits the third; raising the activity threshold above
    the largest in-range difference of 2.0 drops both. Asserting one setting
    on its own would also pass for an implementation that ignored the keyword
    and happened to agree with it.
    """
    default = oecluster.activity_landscape(_line_dm(_SALI_COORDS),
                                           _SALI_ACTIVITY)
    wider = oecluster.activity_landscape(
        _line_dm(_SALI_COORDS), _SALI_ACTIVITY, distance_threshold=0.6)
    stricter = oecluster.activity_landscape(
        _line_dm(_SALI_COORDS), _SALI_ACTIVITY, activity_threshold=2.5)

    assert default.num_cliffs == 2
    assert wider.num_cliffs == 3
    assert wider.cliff_density == 1.0
    assert stricter.num_cliffs == 0
    assert stricter.cliff_density == 0.0


def test_activity_landscape_rmodi_delta_moves_rmodi():
    """The third tuning keyword is read, and it moves RMODI off its default.

    The other two thresholds are covered by the test above; rmodi_delta is the
    one keyword no other test in this file reads. Setting ``rmodi_delta =
    0.625`` unconditionally inside the entry point -- ignoring whatever the
    caller passed -- leaves every other assertion in this file passing.

    Three settings, not one. Asserting a single value would also pass for an
    implementation that ignored the keyword and happened to agree with it, so
    the default case is deliberately spelled as a call passing no keyword at
    all, and the two bracketing cases have to disagree with it. The band is
    ``rmodi_delta`` standard deviations either side of a molecule's own
    activity: narrow it to 0.05 and almost nothing shares a band; widen it to
    3.0 and everything does.
    """
    dm = _line_dm([0.0, 0.1, 0.3, 0.6, 1.0, 1.5])
    activity = [1.0, 1.2, 3.0, 3.1, 8.0, 8.4]

    default = oecluster.activity_landscape(dm, activity)
    narrow = oecluster.activity_landscape(dm, activity, rmodi_delta=0.05)
    wide = oecluster.activity_landscape(dm, activity, rmodi_delta=3.0)

    assert default.rmodi == pytest.approx(0.8333333333333334, abs=1e-12)
    assert narrow.rmodi == pytest.approx(0.16666666666666666, abs=1e-12)
    assert wide.rmodi == pytest.approx(1.0, abs=1e-12)

    # The keyword moves the band, not the data: the standard deviation the band
    # is measured in is the same in all three, so the three RMODI values differ
    # because of the setting and nothing else.
    for landscape in (default, narrow, wide):
        assert landscape.activity_stddev == pytest.approx(2.998008598312479,
                                                          abs=1e-12)


def test_a_zero_distance_pair_is_a_cliff_only_above_the_activity_threshold():
    """The pair is counted as zero-distance either way; whether it is also a
    cliff still depends on the activity difference."""
    dm = _line_dm([0.0, 0.0, 1.0])

    equal = oecluster.activity_landscape(dm, [5.0, 5.0, 9.0])
    differing = oecluster.activity_landscape(dm, [5.0, 9.0, 9.0])

    assert equal.num_zero_distance_pairs == 1
    assert differing.num_zero_distance_pairs == 1
    assert equal.num_cliffs == 0
    assert differing.num_cliffs == 1


def test_activity_landscape_to_table_and_repr():
    landscape = oecluster.activity_landscape(_line_dm(_SALI_COORDS),
                                             _SALI_ACTIVITY)
    table = landscape.to_table()

    assert [name for name, _ in table] == [
        "num_samples", "num_scored", "num_pairs_scored", "num_cliffs",
        "cliff_density", "num_zero_distance_pairs", "max_sali", "mean_sali",
        "rmodi", "activity_stddev",
    ]
    # As in the coherence case: the right ten names against empty cells would
    # satisfy the list above and the substring check below.
    assert dict(table)["num_cliffs"] == 2
    assert dict(table)["cliff_density"] == pytest.approx(2.0 / 3.0, abs=1e-12)

    rendered = repr(landscape)
    assert "cliff_density" in rendered
    assert "0.6667" in rendered


def test_activity_landscape_rejects_a_non_matrix():
    with pytest.raises(TypeError, match="SymmetricDistanceMatrix"):
        oecluster.activity_landscape([[0.0, 0.25], [0.25, 0.0]], [1.0, 2.0])


def test_activity_landscape_rejects_sparse_storage():
    """ValueError, not TypeError: the argument's type is right, its storage
    is not."""
    storage = SparseStorage(3, 0.9)
    storage.Set(0, 1, 0.25)
    storage.Set(1, 2, 0.25)
    storage.Finalize()
    dm = SymmetricDistanceMatrix(storage, "test", ["a", "b", "c"], {})

    with pytest.raises(ValueError, match="SparseStorage"):
        oecluster.activity_landscape(dm, _SALI_ACTIVITY)


def test_activity_landscape_rejects_a_length_mismatch():
    with pytest.raises(ValueError, match="3 samples"):
        oecluster.activity_landscape(_line_dm(_SALI_COORDS), [1.0, 2.0])


def test_activity_landscape_rejects_an_empty_activity():
    """A zero-sample matrix is constructible, so 0 == 0 satisfies the length
    check and only the emptiness check stands between this call and a
    RuntimeError from C++."""
    with pytest.raises(ValueError, match="non-empty activity"):
        oecluster.activity_landscape(_line_dm([]), [])


def test_activity_landscape_rejects_a_negative_num_threads():
    with pytest.raises(ValueError, match="num_threads must be non-negative"):
        oecluster.activity_landscape(_line_dm(_SALI_COORDS), _SALI_ACTIVITY,
                                     num_threads=-1)


def test_activity_landscape_refuses_every_bad_threshold_as_value_error():
    """All three thresholds are refused in Python, not left to C++.

    SWIG maps every native exception to ``RuntimeError``, so a threshold that
    only C++ validates comes back as the wrong exception type for a condition
    the caller could see for itself. These mirror ``validate_landscape_options``
    in ``src/clustering/SARCoherence.cpp``: the same three names, finiteness
    tested before sign, and the same two messages.

    Infinity as well as NaN for each keyword. A check written as ``math.isnan``
    would satisfy the NaN rows and let every infinity through -- and an infinite
    ``distance_threshold`` does not even fail downstream; it silently calls
    every pair structurally near.
    """
    dm = _line_dm(_SALI_COORDS)
    # Each keyword is spelled out in its own call rather than unpacked from a
    # mapping, so that a typo here is a type error rather than a TypeError the
    # pytest.raises below would report as a missing ValueError.
    thresholds = (
        ("distance_threshold",
         lambda value: oecluster.activity_landscape(
             dm, _SALI_ACTIVITY, distance_threshold=value)),
        ("activity_threshold",
         lambda value: oecluster.activity_landscape(
             dm, _SALI_ACTIVITY, activity_threshold=value)),
        ("rmodi_delta",
         lambda value: oecluster.activity_landscape(
             dm, _SALI_ACTIVITY, rmodi_delta=value)),
    )

    for keyword, call in thresholds:
        # Negative infinity is the row that pins the ordering: it is the only
        # input both checks would refuse, so it is the only one whose message
        # differs depending on which runs first. Dropping it as redundant with
        # the other two infinities would leave the order unpinned again.
        for bad, expected in ((math.nan, "must be finite"),
                              (math.inf, "must be finite"),
                              (-math.inf, "must be finite"),
                              (-1.0, "must be non-negative")):
            with pytest.raises(ValueError, match=f"{keyword} {expected}"):
                call(bad)

    # An int beyond double range fails in the cast rather than at isfinite, and
    # would otherwise leave as OverflowError.
    for keyword, call in thresholds:
        with pytest.raises(ValueError, match=f"{keyword} must be finite"):
            call(10 ** 1000)

    # Zero is the boundary and has to be accepted. Without this, a validator
    # that refused all three keywords outright would satisfy the block above.
    for _, call in thresholds:
        assert call(0.0).num_samples == 3


def test_activity_landscape_rejects_a_misspelled_keyword():
    """A typo'd tuning keyword refuses; see the coherence case for why."""
    assert oecluster.activity_landscape(
        _line_dm(_SALI_COORDS), _SALI_ACTIVITY,
        distance_threshold=0.6).num_cliffs == 3

    with pytest.raises(TypeError, match="distance_treshold"):
        oecluster.activity_landscape(
            _line_dm(_SALI_COORDS), _SALI_ACTIVITY,
            distance_treshold=0.6)  # pyright: ignore[reportCallIssue]


def test_activity_landscape_refuses_subset_scored_distances():
    """The gate reaches the public entry point, not just _gate's own tests."""
    dm = _line_dm(_SALI_COORDS, facts={'data_integrity': "subset_scored"})

    with pytest.raises(ValueError) as excinfo:
        oecluster.activity_landscape(dm, _SALI_ACTIVITY)

    message = str(excinfo.value)
    assert "activity_landscape" in message
    assert "not mutually comparable" in message
    assert "missing='complete_case'" in message

    # The refusal is the recorded fact's doing and not something about this
    # fixture: the identical geometry without the stamp scores.
    assert oecluster.activity_landscape(_line_dm(_SALI_COORDS),
                                        _SALI_ACTIVITY).num_cliffs == 2


def test_modelability_matches_a_hand_fixture():
    """0 and 1 are each other's nearest neighbour and share class A; 2 and 3
    are each other's and do not. So A scores 2/3, B scores 0, MODI is 1/3."""
    report = oecluster.modelability(_line_dm(_MODI_COORDS), _MODI_CLASSES)

    assert report.num_samples == 4
    assert report.num_scored == 4
    assert report.num_classes == 2
    assert report.modi == pytest.approx(1.0 / 3.0, abs=1e-12)
    assert [row.label for row in report.classes] == ["A", "B"]
    assert [row.num_members for row in report.classes] == [3, 1]
    assert report.classes[0].fraction_same_class == pytest.approx(
        2.0 / 3.0, abs=1e-12)
    assert report.classes[1].fraction_same_class == 0.0
    assert isinstance(report.classes[0], oecluster.ClassConcordance)


def test_modelability_reports_nan_for_a_single_class():
    """With one class no molecule has a neighbour that could differ."""
    report = oecluster.modelability(_line_dm([0.0, 0.1, 0.5]),
                                    ["A", "A", "A"])

    assert report.num_classes == 1
    assert math.isnan(report.modi)
    assert math.isnan(report.classes[0].fraction_same_class)

    # NaN is not this geometry's standing answer. Relabelling the far point
    # gives the same three molecules two classes and a real score, so the two
    # NaN assertions above are reading the annotation.
    two_classes = oecluster.modelability(_line_dm([0.0, 0.1, 0.5]),
                                         ["A", "A", "B"])
    assert two_classes.modi == 0.5


def test_modelability_ignores_empty_class_strings():
    """An empty string is a missing annotation, not a category.

    Dropping sample 1 leaves 0, 2 and 3, whose nearest scored neighbours are
    2, 3 and 2 -- every one of them a class change.
    """
    report = oecluster.modelability(_line_dm(_MODI_COORDS),
                                    ["A", "", "B", "A"])

    assert report.num_samples == 4
    assert report.num_scored == 3
    assert report.num_classes == 2
    assert report.modi == 0.0

    # Zero is not what these four points always score: annotating sample 1 as
    # A rather than dropping it moves the index to 1/3.
    annotated = oecluster.modelability(_line_dm(_MODI_COORDS), _MODI_CLASSES)
    assert annotated.modi == pytest.approx(1.0 / 3.0, abs=1e-12)


def test_modelability_rejects_a_bare_str():
    """A str is iterable, so without the guard "AAB" is three annotations."""
    with pytest.raises(TypeError, match="not a single str"):
        oecluster.modelability(_line_dm([0.0, 0.1, 0.5]), "AAB")


def test_modelability_accepts_any_iterable_of_class_strings():
    """The other half of "reject str only": a tuple, a one-shot generator and a
    NumPy object array are all accepted, and all score the same as the list."""
    dm, _ = _thread_fixture()
    block = ["A"] * 30 + ["B"] * 30
    expected = 0.9666666666666667

    assert oecluster.modelability(dm, block).modi == pytest.approx(
        expected, abs=1e-12)
    assert oecluster.modelability(dm, tuple(block)).modi == pytest.approx(
        expected, abs=1e-12)
    assert oecluster.modelability(
        dm, (label for label in block)).modi == pytest.approx(
            expected, abs=1e-12)
    numpy = pytest.importorskip("numpy")
    assert oecluster.modelability(
        dm, numpy.array(block, dtype=object)).modi == pytest.approx(
            expected, abs=1e-12)


def test_modelability_refuses_a_non_matrix_and_a_negative_thread_count():
    """Both guards are documented on the signature and neither was pinned;
    without them the caller gets AttributeError and a SWIG OverflowError."""
    with pytest.raises(TypeError, match="expects a SymmetricDistanceMatrix"):
        oecluster.modelability([[0.0]], _MODI_CLASSES)
    with pytest.raises(ValueError, match="num_threads must be non-negative"):
        oecluster.modelability(_line_dm(_MODI_COORDS), _MODI_CLASSES,
                               num_threads=-1)


def test_modelability_rejects_a_non_string_class():
    with pytest.raises(TypeError, match="sequence of class strings"):
        oecluster.modelability(_line_dm([0.0, 0.1, 0.5]), [0, 0, 1])


def test_modelability_rejects_a_length_mismatch():
    with pytest.raises(ValueError, match="4 samples"):
        oecluster.modelability(_line_dm(_MODI_COORDS), ["A", "B"])


def test_modelability_rejects_empty_activity_classes():
    """Per entry point, like the SparseStorage refusal: the zero-sample case
    slips past a length check that compares 0 with 0."""
    with pytest.raises(ValueError, match="non-empty activity_classes"):
        oecluster.modelability(_line_dm([]), [])


def test_modelability_rejects_sparse_storage():
    """The refusal is per entry point, so testing it through
    activity_landscape alone would leave this one free to forget it."""
    storage = SparseStorage(4, 0.9)
    storage.Set(0, 1, 0.1)
    storage.Set(1, 2, 0.4)
    storage.Set(2, 3, 0.2)
    storage.Finalize()
    dm = SymmetricDistanceMatrix(storage, "test", ["a", "b", "c", "d"], {})

    with pytest.raises(ValueError, match="SparseStorage"):
        oecluster.modelability(dm, _MODI_CLASSES)


def test_modelability_refuses_subset_scored_distances():
    dm = _line_dm(_MODI_COORDS, facts={'data_integrity': "subset_scored"})

    with pytest.raises(ValueError) as excinfo:
        oecluster.modelability(dm, _MODI_CLASSES)

    # The caller name as well as the refusal: the gate interpolates whatever
    # string it is handed, so a message asserted on "not mutually comparable"
    # alone reads identically when this entry point passes the wrong name.
    assert "modelability" in str(excinfo.value)
    assert "not mutually comparable" in str(excinfo.value)

    # As in the landscape case: the same geometry without the stamp scores.
    assert oecluster.modelability(
        _line_dm(_MODI_COORDS), _MODI_CLASSES).modi == pytest.approx(
            1.0 / 3.0, abs=1e-12)


def test_modelability_rejects_a_misspelled_keyword():
    """A typo'd tuning keyword refuses; see the coherence case for why."""
    assert oecluster.modelability(_line_dm(_MODI_COORDS), _MODI_CLASSES,
                                  num_threads=1).num_classes == 2

    with pytest.raises(TypeError, match="num_theads"):
        oecluster.modelability(
            _line_dm(_MODI_COORDS), _MODI_CLASSES,
            num_theads=1)  # pyright: ignore[reportCallIssue]


def test_modelability_to_table_and_repr():
    report = oecluster.modelability(_line_dm(_MODI_COORDS), _MODI_CLASSES)
    table = report.to_table()

    assert [name for name, _ in table] == [
        "num_samples", "num_scored", "num_classes", "modi",
    ]
    # The four names against empty cells would satisfy the list above and the
    # substring check below.
    assert dict(table)["num_classes"] == 2
    assert dict(table)["modi"] == pytest.approx(1.0 / 3.0, abs=1e-12)

    rendered = repr(report)
    assert "modi" in rendered
    assert "0.3333" in rendered
    # The per-class table does not render, checked against the row names for
    # the reason given in test_sar_coherence_to_table_and_repr: "num_classes"
    # ends in "classes".
    assert [line.split()[0] for line in rendered.splitlines()[1:]] == [
        "num_samples", "num_scored", "num_classes", "modi",
    ]

    # A NaN field renders as a cell rather than blanking the row or raising in
    # the formatter. The fixture above is all-finite and cannot show it; a
    # single-class matrix gives a NaN modi.
    flat = repr(oecluster.modelability(_line_dm(_MODI_COORDS),
                                       ["A", "A", "A", "A"]))
    assert flat.splitlines()[-1].split() == ["modi", "nan"]


def _thread_fixture():
    """Sixty points on a line with a rising, slightly noisy activity."""
    coordinates = [i * 0.05 for i in range(60)]
    activity = [math.sin(i) + 0.1 * i for i in range(60)]
    return _line_dm(coordinates), activity


def test_the_matrix_metrics_agree_across_thread_counts():
    """Bit for bit, not to a tolerance: the sweep combines its per-row
    partials in ascending row order after the join precisely so the thread
    count cannot move the last digit."""
    dm, activity = _thread_fixture()

    one = oecluster.activity_landscape(dm, activity, num_threads=1)
    many = oecluster.activity_landscape(dm, activity, num_threads=4)

    # Two empty or all-zero tables would agree as readily as two correct ones.
    # This fixture finds cliffs and leaves every reported scalar finite, so the
    # equality below is an agreement about real values.
    assert one.num_cliffs == 131
    assert all(math.isfinite(value) for _, value in one.to_table())

    assert one.to_table() == many.to_table()

    # The block annotation, not the alternating one: the alternating annotation
    # scores exactly 0.0, which a threaded path that computed nothing would also
    # return, so an equality over it agrees about silence. Pinning the
    # single-thread value first makes the equality an agreement about a number.
    block = ["A"] * 30 + ["B"] * 30
    one_thread = oecluster.modelability(dm, block, num_threads=1).modi
    assert one_thread == pytest.approx(0.9666666666666667, abs=1e-12)
    assert one_thread == oecluster.modelability(dm, block, num_threads=4).modi

    # The degenerate annotation still agrees across thread counts, but it is no
    # longer the only thing the modelability half of this test rests on.
    assert (oecluster.modelability(dm, ["A", "B"] * 30, num_threads=1).modi ==
            oecluster.modelability(dm, ["A", "B"] * 30, num_threads=4).modi)


def test_the_matrix_metrics_are_safe_under_concurrent_calls():
    """Four Python threads in one process, scoring the same matrix.

    This is a thread-safety test, not a GIL-release test: the assertion is
    that nothing is shared across calls, so every thread gets the same answer
    whether or not the calls actually overlapped. The GIL release itself is
    asserted against the interface file, in
    ``tests/python/test_native_bindings.py``.
    """
    dm, activity = _thread_fixture()
    # The block annotation rather than an alternating one, for the reason given
    # in the thread-count test: an alternating annotation scores exactly 0.0, so
    # pooled calls that computed nothing would satisfy the final equality.
    classes = ["A"] * 30 + ["B"] * 30
    expected = oecluster.activity_landscape(dm, activity).to_table()
    expected_modi = oecluster.modelability(dm, classes).modi

    # Same control as the thread-count test: pin what the threads are agreeing
    # about, so that a constant or empty answer cannot satisfy the equalities.
    assert dict(expected)["num_cliffs"] == 131
    assert all(math.isfinite(value) for _, value in expected)
    assert expected_modi == pytest.approx(0.9666666666666667, abs=1e-12)

    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
        futures = [pool.submit(oecluster.activity_landscape, dm, activity)
                   for _ in range(8)]
        modi_futures = [pool.submit(oecluster.modelability, dm, classes)
                        for _ in range(8)]
        results = [future.result().to_table() for future in futures]
        modis = [future.result().modi for future in modi_futures]

    assert len(results) == 8
    assert len(modis) == 8
    assert all(result == expected for result in results)
    assert all(modi == expected_modi for modi in modis)


def test_native_refusals_surface_as_runtime_error():
    """SWIG maps every ``std::exception`` to ``RuntimeError``, and A3 adds no
    ``ValueError`` mapping -- the same contract
    ``test_partition_error_surfaces_as_runtime_error`` pins for A1.

    So the conditions Python cannot see for itself arrive as ``RuntimeError``,
    and the docstrings say so. A negative distance is the interesting one: the
    gate measures finiteness, which a negative number passes, so the refusal
    can only come from C++.
    """
    with pytest.raises(RuntimeError, match=r"activity\[1\] is infinite"):
        oecluster.sar_coherence([0, 0, 1], [1.0, math.inf, 3.0])

    negative_storage = DenseStorage(3)
    negative_storage.Set(0, 1, -0.25)
    negative_storage.Set(0, 2, 0.5)
    negative_storage.Set(1, 2, 0.25)
    negative = SymmetricDistanceMatrix(
        negative_storage, "test", ["m0", "m1", "m2"], {}, None)

    with pytest.raises(RuntimeError, match="finite and non-negative"):
        oecluster.modelability(negative, ["A", "A", "B"])


def test_the_sar_coherence_surface_is_exported():
    """Attribute access does not consult __all__, so every test above passes
    with the names missing from it; star-import and API discovery do not."""
    exported = ("sar_coherence", "activity_landscape", "modelability",
                "SARCoherence", "ClusterActivity", "ActivityLandscape",
                "Modelability", "ClassConcordance")

    missing = [name for name in exported if name not in oecluster.__all__]
    assert missing == []
    assert all(hasattr(oecluster, name) for name in exported)
