"""Python surface of ClusteringSpec and select_parameter."""
import math

import numpy as np
import oecluster
import oefp
import pytest
from openeye import oechem

FP_SMILES = ["CCO", "CCCO", "CCCCO", "c1ccccc1", "Cc1ccccc1", "CCc1ccccc1",
             "CC(=O)O", "CC(=O)OC", "CCN", "CCCN", "C1CCCCC1", "c1ccncc1"]

ROSTER = ("butina", "dbscan", "hdbscan", "agglomerative", "k_medoids",
          "bitbirch", "bitbirch_recluster", "bitbirch_refine",
          "sphere_exclusion", "jarvis_patrick", "leiden", "murcko")

# Butina thresholds that give ten singletons, the two blobs, and one cluster.
GRID = (0.05, 0.2, 0.95)


def _mols(smiles_list):
    mols = []
    for idx, smi in enumerate(smiles_list):
        mol = oechem.OEGraphMol()
        oechem.OESmilesToMol(mol, smi)
        mol.SetTitle(f"mol{idx}")
        mols.append(mol)
    return mols


def _blobs():
    """Ten points, two blobs of five: intra-blob 0.1, inter-blob 0.9."""
    values = [0.1 if (i < 5) == (j < 5) else 0.9
              for i in range(10) for j in range(i + 1, 10)]
    return oecluster.SymmetricDistanceMatrix.from_condensed(np.array(values))


def _batch(bits, rows):
    return oefp.OEFPBatch.from_fingerprints(
        [oefp.OEFP.from_on_bits(bits, list(on)) for on in rows])


def _batch_from_bits(bits):
    arr = np.asarray(bits, dtype=np.uint8)
    fingerprints = []
    for row in arr:
        on_bits = np.flatnonzero(row).astype(int).tolist()
        fingerprints.append(oefp.OEFP.from_on_bits(arr.shape[1], on_bits))
    return oefp.OEFPBatch.from_fingerprints(fingerprints)


def _fps():
    """Two groups of five sharing an eight-bit block, one private bit each."""
    rows = ([set(range(8)) | {16 + i} for i in range(5)]
            + [set(range(32, 40)) | {48 + i} for i in range(5)])
    return _batch(64, rows)


def _pairs():
    """p0 == p1 and p2 == p3 (distance 0); every other distance 0.5."""
    zero = {(0, 1), (2, 3)}
    values = [0.0 if (i, j) in zero else 0.5
              for i in range(4) for j in range(i + 1, 4)]
    return oecluster.SymmetricDistanceMatrix.from_condensed(np.array(values))


def split_pairs(items, split):
    """A foreign clusterer: 'across' pairs coincident medoids, 'together' does not."""
    if split == "across":
        return oecluster.ClusteringResult([0, 1, 0, 1], ((0, 2), (1, 3)))
    return oecluster.ClusteringResult([0, 0, 1, 1], ((0, 1), (2, 3)))


def _recording(calls):
    """A clusterer that records the options it was called with."""
    def record(items, **options):
        calls.append(options)
        return oecluster.ClusteringResult([0, 0, 0, 0], ((0, 1, 2, 3),))
    return record


# --- ClusteringSpec ---------------------------------------------------------

@pytest.mark.parametrize("name", ROSTER)
def test_spec_resolves_every_roster_name(name):
    spec = oecluster.ClusteringSpec(name)
    assert spec.name == name
    assert spec.algorithm is getattr(oecluster, name)


def test_spec_name_matching_is_case_insensitive():
    assert oecluster.ClusteringSpec("Butina").name == "butina"


def test_spec_maps_a_roster_function_back_to_its_name():
    spec = oecluster.ClusteringSpec(oecluster.sphere_exclusion)
    assert spec.name == "sphere_exclusion"
    assert spec.algorithm is oecluster.sphere_exclusion


def test_spec_keeps_a_foreign_callable_and_its_name():
    spec = oecluster.ClusteringSpec(split_pairs, split="across")
    assert spec.name == "split_pairs"
    assert spec.algorithm is split_pairs
    assert dict(spec.options) == {"split": "across"}


def test_spec_refuses_an_unknown_name_and_lists_the_roster():
    with pytest.raises(ValueError, match="butina"):
        oecluster.ClusteringSpec("kmeans")


def test_spec_refuses_a_non_callable():
    with pytest.raises(TypeError, match="callable"):
        oecluster.ClusteringSpec(42)


def test_spec_options_are_a_read_only_copy():
    options = {"num_threads": 2}
    spec = oecluster.ClusteringSpec("butina", **options)
    options["num_threads"] = 8
    assert dict(spec.options) == {"num_threads": 2}
    with pytest.raises(TypeError):
        spec.options["num_threads"] = 4


def test_run_passes_fixed_options_and_overrides_win():
    calls = []
    spec = oecluster.ClusteringSpec(_recording(calls), threshold=0.2,
                                    num_threads=1)
    spec.run(_blobs())
    spec.run(_blobs(), threshold=0.5)
    assert calls == [{"threshold": 0.2, "num_threads": 1},
                     {"threshold": 0.5, "num_threads": 1}]


def test_run_returns_the_algorithm_result():
    result = oecluster.ClusteringSpec("butina").run(_blobs(), threshold=0.2)
    assert isinstance(result, oecluster.ButinaResult)
    assert result.num_clusters == 2


def test_run_refuses_a_callable_that_returns_no_result():
    spec = oecluster.ClusteringSpec(lambda items, **options: None)
    with pytest.raises(TypeError, match="ClusteringResult"):
        spec.run(_blobs())


def test_run_refuses_a_graph_as_a_partition():
    # knn_graph is public and callable but returns a KNNGraph, not a result.
    with pytest.raises(TypeError, match="ClusteringResult"):
        oecluster.ClusteringSpec(oecluster.knn_graph).run(_blobs(), k=3)


def test_spec_equality_and_unhashability():
    spec = oecluster.ClusteringSpec("butina", num_threads=2)
    assert spec == oecluster.ClusteringSpec(oecluster.butina, num_threads=2)
    assert spec != oecluster.ClusteringSpec("butina", num_threads=3)
    assert spec != oecluster.ClusteringSpec("dbscan", num_threads=2)
    assert spec != "butina"
    with pytest.raises(TypeError, match="unhashable"):
        hash(spec)


def test_spec_equality_handles_array_valued_options():
    # numpy's == is elementwise, so plain dict equality would raise here.
    left = oecluster.ClusteringSpec("k_medoids", initial_medoids=np.array([0, 5]))
    right = oecluster.ClusteringSpec("k_medoids", initial_medoids=np.array([0, 5]))
    assert left == right
    assert left != oecluster.ClusteringSpec("k_medoids", initial_medoids=np.array([0, 6]))
    assert left != oecluster.ClusteringSpec("k_medoids", initial_medoids=np.array([0, 5, 9]))
    assert left != oecluster.ClusteringSpec("k_medoids", initial_medoids=5)


def test_spec_repr():
    assert (repr(oecluster.ClusteringSpec("butina", num_threads=4))
            == "ClusteringSpec('butina', num_threads=4)")
    assert repr(oecluster.ClusteringSpec("leiden")) == "ClusteringSpec('leiden')"


def test_package_exports_clustering_spec():
    assert "ClusteringSpec" in oecluster.__all__


# --- select_parameter over a matrix -----------------------------------------

def test_default_criterion_picks_the_best_silhouette():
    selection = oecluster.select_parameter("butina", _blobs(), "threshold", GRID)
    assert selection.criterion == "silhouette"
    assert selection.parameter == "threshold"
    assert selection.spec == oecluster.ClusteringSpec("butina")
    assert [row.value for row in selection.rows] == list(GRID)
    scores = [row.score for row in selection.rows]
    assert scores[0] == 0.0
    assert scores[1] == pytest.approx(0.8889, abs=1e-4)
    assert math.isnan(scores[2])
    assert selection.winner_index == 1
    assert selection.winner is selection.rows[1]
    assert selection.winner.value == 0.2
    assert selection.winner.result.num_clusters == 2
    assert isinstance(selection.winner.report, oecluster.ClusterReport)
    assert all(row.eligible and row.rejection is None for row in selection.rows)


def test_rows_carry_the_partition_and_report_of_each_value():
    selection = oecluster.select_parameter(oecluster.butina, _blobs(),
                                           "threshold", GRID)
    assert [row.result.num_clusters for row in selection.rows] == [10, 2, 1]
    assert [row.report.num_clusters for row in selection.rows] == [10, 2, 1]
    assert all(row.result.method == "butina" for row in selection.rows)
    assert all(isinstance(row, oecluster.SweepRow) for row in selection.rows)


def test_sweep_row_fields():
    assert oecluster.SweepRow._fields == (
        "value", "result", "report", "score", "eligible", "rejection")


def test_a_spec_is_run_with_its_fixed_options():
    spec = oecluster.ClusteringSpec("butina", num_threads=1)
    selection = oecluster.select_parameter(spec, _blobs(), "threshold", GRID)
    assert selection.spec is spec
    assert selection.winner.value == 0.2


def test_equal_scores_keep_grid_order():
    selection = oecluster.select_parameter("butina", _blobs(), "threshold",
                                           (0.2, 0.2))
    assert selection.winner_index == 0


def test_all_nan_scores_give_no_winner():
    selection = oecluster.select_parameter("butina", _blobs(), "threshold",
                                           (0.95,))
    assert selection.winner is None
    assert selection.winner_index is None
    assert len(selection.rows) == 1


def test_a_value_the_algorithm_refuses_propagates_its_own_error():
    with pytest.raises(ValueError, match="non-negative"):
        oecluster.select_parameter("butina", _blobs(), "threshold", (0.2, -1.0))


def test_values_may_be_a_generator():
    selection = oecluster.select_parameter("butina", _blobs(), "threshold",
                                           (t for t in GRID))
    assert [row.value for row in selection.rows] == list(GRID)


def test_report_options_are_forwarded_to_the_scorer():
    selection = oecluster.select_parameter(
        "butina", _blobs(), "threshold", (0.2,),
        report_options={"preset": "tight"})
    assert selection.rows[0].report.coverage_thresholds == (0.2, 0.3, 0.4)


# --- bounds and ranking -----------------------------------------------------

def test_davies_bouldin_without_bounds_rewards_the_singleton_partition():
    selection = oecluster.select_parameter(
        "butina", _blobs(), "threshold", GRID, criterion="davies_bouldin_medoid")
    assert selection.winner.value == 0.05
    assert selection.winner.score == 0.0


def test_davies_bouldin_under_a_cluster_cap_picks_the_two_blobs():
    selection = oecluster.select_parameter(
        "butina", _blobs(), "threshold", GRID, criterion="davies_bouldin_medoid",
        max_clusters=2)
    assert selection.winner.value == 0.2
    assert selection.winner.score == pytest.approx(0.1778, abs=1e-4)
    assert selection.rows[0].eligible is False
    assert selection.rows[0].rejection == "num_clusters 10 > max_clusters 2"


def test_min_clusters_rejects_rows_with_the_exact_text():
    selection = oecluster.select_parameter("butina", _blobs(), "threshold",
                                           GRID, min_clusters=3)
    assert [row.eligible for row in selection.rows] == [True, False, False]
    assert selection.rows[1].rejection == "num_clusters 2 < min_clusters 3"
    assert selection.rows[2].rejection == "num_clusters 1 < min_clusters 3"
    assert selection.winner_index == 0
    assert selection.winner.score == 0.0


def test_min_clusters_normalizes_a_numpy_integer_bound():
    selection = oecluster.select_parameter("butina", _blobs(), "threshold",
                                           GRID, min_clusters=np.int64(3))
    assert selection.rows[1].rejection == "num_clusters 2 < min_clusters 3"


def test_a_bound_leaving_only_nan_rows_gives_no_winner():
    selection = oecluster.select_parameter("butina", _blobs(), "threshold",
                                           GRID, max_clusters=1)
    assert [row.eligible for row in selection.rows] == [False, False, True]
    assert selection.winner is None


def test_max_noise_fraction_rejects_an_all_noise_row():
    spec = oecluster.ClusteringSpec("dbscan", min_samples=6)
    selection = oecluster.select_parameter(spec, _blobs(), "eps", (0.2,),
                                           max_noise_fraction=0.5)
    row, = selection.rows
    assert row.report.noise_fraction == 1.0
    assert row.eligible is False
    assert row.rejection == "noise_fraction 1.0 > max_noise_fraction 0.5"
    assert selection.winner is None


def test_the_first_violated_bound_is_the_one_reported():
    spec = oecluster.ClusteringSpec("dbscan", min_samples=6)
    selection = oecluster.select_parameter(spec, _blobs(), "eps", (0.2,),
                                           max_noise_fraction=0.5, min_clusters=2)
    assert selection.rows[0].rejection == (
        "noise_fraction 1.0 > max_noise_fraction 0.5")


def test_a_noise_cap_met_exactly_keeps_the_row_eligible():
    spec = oecluster.ClusteringSpec("dbscan", min_samples=3)
    selection = oecluster.select_parameter(spec, _blobs(), "eps", (0.2,),
                                           max_noise_fraction=0.0)
    assert selection.rows[0].eligible
    assert selection.winner.value == 0.2


def test_an_infinite_score_is_rankable_and_loses_to_a_finite_one():
    selection = oecluster.select_parameter(
        split_pairs, _pairs(), "split", ("across", "together"),
        criterion="davies_bouldin_medoid")
    assert selection.rows[0].score == math.inf
    assert selection.rows[0].eligible
    assert selection.rows[1].score == 0.0
    assert selection.winner_index == 1


def test_an_infinite_score_wins_when_it_is_the_only_row():
    selection = oecluster.select_parameter(
        split_pairs, _pairs(), "split", ("across",),
        criterion="davies_bouldin_medoid")
    assert selection.winner_index == 0
    assert selection.winner.score == math.inf


# --- optional stages --------------------------------------------------------

@pytest.mark.parametrize("criterion", ["c_index", "baker_hubert_gamma"])
def test_a_pair_rank_criterion_switches_the_stage_on(criterion):
    options = {"preset": "default"}
    selection = oecluster.select_parameter(
        "butina", _blobs(), "threshold", GRID, criterion=criterion,
        report_options=options)
    assert all(row.report.requested.pair_rank_indices for row in selection.rows)
    assert math.isnan(selection.rows[0].score)
    assert math.isnan(selection.rows[2].score)
    assert selection.winner.value == 0.2
    assert options == {"preset": "default"}


def test_c_index_is_minimized():
    selection = oecluster.select_parameter(
        "butina", _blobs(), "threshold", GRID, criterion="c_index")
    assert selection.winner.score == 0.0


@pytest.mark.parametrize("criterion", ["c_index", "baker_hubert_gamma"])
def test_a_stage_the_caller_switched_off_is_refused(criterion):
    calls = []
    with pytest.raises(ValueError, match="compute_pair_rank_indices"):
        oecluster.select_parameter(
            _recording(calls), _blobs(), "threshold", GRID, criterion=criterion,
            report_options={"compute_pair_rank_indices": False})
    assert calls == []


def test_a_stage_the_caller_switched_on_is_kept():
    selection = oecluster.select_parameter(
        "butina", _blobs(), "threshold", (0.2,),
        report_options={"compute_pair_rank_indices": True})
    assert selection.criterion == "silhouette"
    assert selection.rows[0].report.requested.pair_rank_indices is True


# --- validation before the first run ----------------------------------------

def test_parameter_must_be_a_non_empty_string():
    calls = []
    with pytest.raises(TypeError, match="parameter"):
        oecluster.select_parameter(_recording(calls), _blobs(), 3, GRID)
    with pytest.raises(ValueError, match="parameter"):
        oecluster.select_parameter(_recording(calls), _blobs(), "", GRID)
    assert calls == []


def test_values_must_be_a_non_empty_iterable_other_than_a_string():
    calls = []
    with pytest.raises(TypeError, match="values"):
        oecluster.select_parameter(_recording(calls), _blobs(), "threshold",
                                   "0.2")
    with pytest.raises(TypeError, match="values"):
        oecluster.select_parameter(_recording(calls), _blobs(), "threshold", 0.2)
    with pytest.raises(ValueError, match="values"):
        oecluster.select_parameter(_recording(calls), _blobs(), "threshold", ())
    assert calls == []


def test_report_options_must_be_a_mapping():
    calls = []
    with pytest.raises(TypeError, match="report_options"):
        oecluster.select_parameter(_recording(calls), _blobs(), "threshold",
                                   GRID, report_options=[("preset", "tight")])
    assert calls == []


@pytest.mark.parametrize("kwargs, error, match", [
    ({"max_noise_fraction": "0.5"}, TypeError, "max_noise_fraction"),
    ({"max_noise_fraction": True}, TypeError, "max_noise_fraction"),
    ({"max_noise_fraction": 1.5}, ValueError, "max_noise_fraction"),
    ({"max_noise_fraction": -0.1}, ValueError, "max_noise_fraction"),
    ({"max_noise_fraction": float("nan")}, ValueError, "max_noise_fraction"),
    ({"min_clusters": 2.0}, TypeError, "min_clusters"),
    ({"min_clusters": True}, TypeError, "min_clusters"),
    ({"min_clusters": 0}, ValueError, "min_clusters"),
    ({"max_clusters": 0}, ValueError, "max_clusters"),
    ({"min_clusters": 3, "max_clusters": 2}, ValueError, "min_clusters"),
])
def test_bounds_are_validated_before_clustering(kwargs, error, match):
    calls = []
    with pytest.raises(error, match=match):
        oecluster.select_parameter(_recording(calls), _blobs(), "threshold",
                                   GRID, **kwargs)
    assert calls == []


def test_criterion_must_be_a_validity_index():
    calls = []
    with pytest.raises(ValueError, match="num_clusters"):
        oecluster.select_parameter(_recording(calls), _blobs(), "threshold",
                                   GRID, criterion="num_clusters")
    # The refusal lists the allowed fields.
    with pytest.raises(ValueError, match="silhouette"):
        oecluster.select_parameter(_recording(calls), _blobs(), "threshold",
                                   GRID, criterion="mean_intra_distance")
    with pytest.raises(TypeError, match="criterion"):
        oecluster.select_parameter(_recording(calls), _blobs(), "threshold",
                                   GRID, criterion=3)
    assert calls == []


def test_unsupported_item_kinds_are_refused_before_clustering():
    calls = []
    blobs = _blobs()
    with pytest.raises(TypeError, match="prebuilt comparison"):
        oecluster.select_parameter(_recording(calls), [1, 2, 3], "threshold",
                                   GRID)
    with pytest.raises(TypeError):
        oecluster.select_parameter(_recording(calls), oecluster.knn_graph(blobs, 3),
                                   "threshold", GRID)
    with pytest.raises(TypeError, match="CrossDistanceMatrix"):
        oecluster.select_parameter(
            _recording(calls),
            oecluster.CrossDistanceMatrix(np.zeros((2, 3)), "test"),
            "threshold", GRID)
    assert calls == []


def test_a_sparse_matrix_is_refused_before_clustering():
    storage = oecluster.SparseStorage(4, 0.5)
    storage.Set(0, 1, 0.2)
    storage.Set(2, 3, 0.2)
    storage.Finalize()
    sparse = oecluster.SymmetricDistanceMatrix(storage, "test",
                                               ["a", "b", "c", "d"], {})
    calls = []
    with pytest.raises(ValueError, match="SparseStorage"):
        oecluster.select_parameter(_recording(calls), sparse, "threshold", GRID)
    assert calls == []


# --- the direction table ----------------------------------------------------

def test_every_criterion_is_a_scorecard_field_with_a_direction():
    from oecluster._parameter_selection import _CRITERIA

    cluster_fields = set(oecluster.ClusterReport._SCALAR_FIELDS)
    isim_fields = set(oecluster.ISimReport._SCALAR_FIELDS)
    assert set(_CRITERIA) == {
        "silhouette", "isim_silhouette", "dunn_index",
        "dunn_mean_separation_mean_diameter",
        "dunn_medoid_separation_medoid_spread", "calinski_harabasz_medoid",
        "point_biserial", "baker_hubert_gamma", "davies_bouldin_medoid",
        "c_index"}
    assert set(_CRITERIA.values()) <= {1, -1}
    assert {name for name, direction in _CRITERIA.items() if direction == -1} == {
        "davies_bouldin_medoid", "c_index"}
    assert set(_CRITERIA) - {"isim_silhouette"} <= cluster_fields
    assert {"isim_silhouette", "dunn_medoid_separation_medoid_spread",
            "calinski_harabasz_medoid", "davies_bouldin_medoid"} <= isim_fields
    for profile in ("num_clusters", "noise_fraction", "mean_intra_distance"):
        assert profile not in _CRITERIA


# --- the selection object ---------------------------------------------------

def test_selection_is_read_only():
    selection = oecluster.select_parameter("butina", _blobs(), "threshold", (0.2,))
    for name in ("spec", "parameter", "criterion", "rows", "winner",
                 "winner_index", "columns", "_rows", "_winner_index", "extra"):
        with pytest.raises(AttributeError, match="read-only"):
            setattr(selection, name, None)
    assert selection.winner_index == 0


def test_to_table_and_columns():
    selection = oecluster.select_parameter("butina", _blobs(), "threshold",
                                           GRID, min_clusters=2)
    assert selection.columns == ("threshold", "silhouette", "num_clusters",
                                 "noise_fraction", "eligible", "rejection")
    table = selection.to_table()
    assert isinstance(table, list)
    assert len(table) == 3
    value, score, clusters, noise, eligible, rejection = table[1]
    assert (value, clusters, noise, eligible, rejection) == (0.2, 2, 0.0, True, None)
    assert score == pytest.approx(0.8889, abs=1e-4)
    assert table[2][0] == 0.95
    assert math.isnan(table[2][1])
    assert table[2][4] is False
    assert table[2][5] == "num_clusters 1 < min_clusters 2"
    assert table is not selection.to_table()


def test_repr_marks_the_winner_and_rejections():
    selection = oecluster.select_parameter("butina", _blobs(), "threshold",
                                           GRID, min_clusters=2)
    lines = repr(selection).splitlines()
    assert lines[0] == ("ParameterSelection(spec=ClusteringSpec('butina'), "
                        "parameter='threshold', criterion='silhouette')")
    assert lines[1].split() == ["threshold", "silhouette", "num_clusters",
                                "noise_fraction"]
    assert lines[2].startswith("   0.05")
    assert lines[3].startswith(" * 0.2")
    assert "0.8889" in lines[3]
    assert lines[4].startswith("   0.95")
    assert lines[4].endswith("rejected: num_clusters 1 < min_clusters 2")
    assert "nan" in lines[4]


def test_repr_formats_a_numpy_grid_value_without_the_np_wrapper():
    grid = np.linspace(0.05, 0.95, 3)
    selection = oecluster.select_parameter(
        "butina", _blobs(), "threshold", grid, min_clusters=2)
    text = repr(selection)
    assert "np.float64" not in text
    for value in grid:
        assert repr(value.item()) in text
    assert "rejected: num_clusters 1 < min_clusters 2" in text


def test_package_exports_the_four_names():
    for name in ("ClusteringSpec", "select_parameter", "ParameterSelection",
                 "SweepRow"):
        assert name in oecluster.__all__
        assert getattr(oecluster, name) is not None


# --- fingerprints -----------------------------------------------------------

def test_fingerprints_score_through_isim_report():
    selection = oecluster.select_parameter("bitbirch", _fps(), "threshold",
                                           (0.2, 0.6, 0.95))
    assert selection.criterion == "isim_silhouette"
    assert all(isinstance(row.report, oecluster.ISimReport)
               for row in selection.rows)
    assert all(row.report.requested.centroid_indices for row in selection.rows)
    assert [row.report.num_clusters for row in selection.rows] == [1, 2, 10]
    assert math.isnan(selection.rows[0].score)
    assert selection.rows[1].score == pytest.approx(0.8)
    assert selection.rows[2].score == 0.0
    assert selection.winner.value == 0.6


def test_fingerprint_options_mapping_is_left_unchanged():
    options = {}
    oecluster.select_parameter("bitbirch", _fps(), "threshold", (0.6,),
                               report_options=options)
    assert options == {}


def test_a_linear_criterion_leaves_the_centroid_stage_off():
    selection = oecluster.select_parameter(
        "bitbirch", _fps(), "threshold", (0.6,),
        criterion="calinski_harabasz_medoid")
    assert selection.rows[0].report.requested.centroid_indices is False
    assert selection.rows[0].score == pytest.approx(125.0)


@pytest.mark.parametrize("criterion, expected", [
    ("isim_silhouette", 0.8),
    ("davies_bouldin_medoid", 0.32),
    ("dunn_medoid_separation_medoid_spread", 3.125),
])
def test_a_centroid_stage_criterion_switches_the_stage_on(criterion, expected):
    options = {}
    selection = oecluster.select_parameter(
        "bitbirch", _fps(), "threshold", (0.6,), criterion=criterion,
        report_options=options)
    assert selection.rows[0].report.requested.centroid_indices is True
    assert selection.rows[0].score == pytest.approx(expected)
    assert options == {}


def test_davies_bouldin_on_fingerprints_under_a_cluster_cap():
    selection = oecluster.select_parameter(
        "bitbirch", _fps(), "threshold", (0.6, 0.95),
        criterion="davies_bouldin_medoid", max_clusters=5)
    assert selection.rows[1].rejection == "num_clusters 10 > max_clusters 5"
    assert selection.winner.value == 0.6


def test_a_criterion_the_scorer_does_not_produce_is_refused():
    calls = []
    with pytest.raises(ValueError, match="isim_report"):
        oecluster.select_parameter(_recording(calls), _fps(), "threshold",
                                   (0.6,), criterion="dunn_index")
    assert calls == []


# --- prebuilt comparisons ---------------------------------------------------

def _same_rows(a, b):
    """to_table() equality that treats NaN as equal to NaN."""
    if len(a) != len(b):
        return False
    for row_a, row_b in zip(a, b):
        for cell_a, cell_b in zip(row_a, row_b):
            both_nan = (isinstance(cell_a, float) and isinstance(cell_b, float)
                        and math.isnan(cell_a) and math.isnan(cell_b))
            if not both_nan and cell_a != cell_b:
                return False
    return True


def test_a_prebuilt_comparison_matches_the_matrix_path():
    mols = _mols(FP_SMILES)
    grid = (0.6, 0.8, 0.9)
    by_matrix = oecluster.select_parameter(
        "sphere_exclusion", oecluster.pdist(mols, "fingerprint"), "threshold",
        grid)
    prebuilt = oecluster.select_parameter(
        "sphere_exclusion", oecluster.FingerprintComparison(mols), "threshold",
        grid)
    assert [row.report.num_clusters for row in by_matrix.rows] == [8, 4, 3]
    assert _same_rows(prebuilt.to_table(), by_matrix.to_table())
    assert prebuilt.winner_index == by_matrix.winner_index == 1
    assert ([row.result.labels.tolist() for row in prebuilt.rows]
            == [row.result.labels.tolist() for row in by_matrix.rows])
    assert all(isinstance(row.report, oecluster.ClusterReport)
               for row in prebuilt.rows)


@pytest.mark.parametrize("criterion", ["c_index", "baker_hubert_gamma"])
def test_pair_rank_criteria_need_a_matrix(criterion):
    calls = []
    prebuilt = oecluster.FingerprintComparison(_mols(FP_SMILES))
    with pytest.raises(ValueError, match="SymmetricDistanceMatrix"):
        oecluster.select_parameter(_recording(calls), prebuilt, "threshold",
                                   (0.5,), criterion=criterion)
    assert calls == []


# --- a partition the scorer refuses ----------------------------------------

def test_a_partition_the_scorer_refuses_propagates_its_error():
    # bitbirch_refine can hand back an emptied leaf subcluster as an empty
    # member list; both scorecards refuse an empty cluster.
    bits = [[1, 1, 1], [0, 1, 1], [1, 0, 0], [1, 1, 0]]
    spec = oecluster.ClusteringSpec(
        "bitbirch_refine", branching_factor=2, merge_criterion="diameter",
        singly=False, redistribute_largest_cluster=True)
    with pytest.raises(RuntimeError, match="at least one member"):
        oecluster.select_parameter(spec, _batch_from_bits(bits), "threshold",
                                   (0.8,))
