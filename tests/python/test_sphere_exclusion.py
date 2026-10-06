"""Python surface of sphere_exclusion: orders, assignment, dispatch, refusals."""

import math

import numpy as np
import oecluster
import pytest
from openeye import oechem

# Default Morgan/Tanimoto distances with exact 1.0 ties; at 0.75 the input
# order keeps the same five centers as the #Circles sequential packing.
FP_SMILES = ["CCO", "CCCO", "CCCCO", "c1ccccc1", "Cc1ccccc1", "CCc1ccccc1",
             "CC(=O)O", "CC(=O)OC", "CCN", "CCCN", "C1CCCCC1", "c1ccncc1"]

# Prefixing water makes the descriptor complete-case mask drop position 0.
DESCRIPTOR_SMILES = ["CCO", "CCC", "CCCC", "c1ccccc1", "CCN", "CCOC",
                     "CC(=O)O", "CCCCCC"]

# Five points on a line at these coordinates; the distance is the gap.
_LINE = (0.0, 1.0, 3.0, 7.0, 8.0)


def _mols(smiles_list):
    mols = []
    for idx, smi in enumerate(smiles_list):
        mol = oechem.OEGraphMol()
        oechem.OESmilesToMol(mol, smi)
        mol.SetTitle(f"mol{idx}")
        mols.append(mol)
    return mols


def _dense_distance_matrix(square):
    """Build a dense SymmetricDistanceMatrix from a square distance list."""
    storage = oecluster.DenseStorage(len(square))
    for i in range(len(square)):
        for j in range(i + 1, len(square)):
            storage.Set(i, j, float(square[i][j]))
    return oecluster.SymmetricDistanceMatrix(
        storage, "test", [f"item_{i}" for i in range(len(square))], {})


def _points_matrix(points):
    return _dense_distance_matrix(
        [[abs(a - b) for b in points] for a in points])


def _line_matrix():
    return _points_matrix(_LINE)


def _compare_matrix(comparison):
    """The matrix a lazy path must reproduce: every pair read via Compare."""
    n = comparison.Size()
    storage = oecluster.DenseStorage(n)
    for i in range(n):
        for j in range(i + 1, n):
            storage.Set(i, j, comparison.Compare(i, j))
    return oecluster.SymmetricDistanceMatrix(
        storage, "test", [f"item_{i}" for i in range(n)], {})


def test_the_input_order_is_leader_clustering():
    result = oecluster.sphere_exclusion(_line_matrix(), 1.5)
    assert result.clusters == ((0, 1), (2,), (3, 4))
    assert result.labels.tolist() == [0, 0, 1, 2, 2]
    assert result.centers == (0, 2, 3)
    assert result.method == "sphere_exclusion"
    assert result.excluded == []


def test_the_neighbors_order_is_butina():
    matrix = _line_matrix()
    result = oecluster.sphere_exclusion(matrix, 1.5, order="NEIGHBORS")
    assert result.clusters == ((4, 3), (1, 0), (2,))
    assert result.clusters == oecluster.butina(matrix, 1.5).clusters
    assert result.centers == (4, 1, 2)


def test_a_permutation_order_takes_centers_in_that_order():
    expected = ((4, 3), (2,), (1, 0))
    for order in ([4, 3, 2, 1, 0], (4, 3, 2, 1, 0), np.array([4, 3, 2, 1, 0]),
                  range(4, -1, -1)):
        result = oecluster.sphere_exclusion(_line_matrix(), 1.5, order=order)
        assert result.clusters == expected
        assert result.centers == (4, 2, 1)
        assert result.labels.tolist() == [2, 2, 1, 0, 0]


def test_nearest_assignment_moves_items_to_their_closest_center():
    """Item 1 is 1.5 from both centers and stays with the earlier one."""
    matrix = _points_matrix((0.0, 1.5, 3.0, 1.6))
    first = oecluster.sphere_exclusion(matrix, 2.0)
    nearest = oecluster.sphere_exclusion(matrix, 2.0, assignment="Nearest")
    assert first.clusters == ((0, 1, 3), (2,))
    assert nearest.clusters == ((0, 1), (2, 3))
    assert nearest.labels.tolist() == [0, 0, 1, 1]
    assert nearest.centers == first.centers == (0, 2)


def test_the_three_call_styles_agree():
    mols = _mols(FP_SMILES)
    by_matrix = oecluster.sphere_exclusion(
        oecluster.pdist(mols, "fingerprint"), 0.75)
    by_object = oecluster.sphere_exclusion(
        oecluster.FingerprintComparison(mols), 0.75)
    by_name = oecluster.sphere_exclusion(mols, 0.75, comparison="fingerprint")

    assert list(by_matrix.centers) == [0, 3, 5, 6, 10]
    packing = oecluster.circles(mols, comparison="fingerprint",
                                threshold=0.75, method="sequential")
    assert list(by_matrix.centers) == packing.members
    for other in (by_object, by_name):
        assert other.clusters == by_matrix.clusters
        assert other.labels.tolist() == by_matrix.labels.tolist()
        assert other.centers == by_matrix.centers
        assert other.excluded == []


@pytest.mark.parametrize("reordering", [False, True])
def test_the_neighbors_order_matches_butina_on_fingerprints(reordering):
    matrix = oecluster.pdist(_mols(FP_SMILES), "fingerprint")
    butina = oecluster.butina(matrix, 0.75, reordering=reordering)
    result = oecluster.sphere_exclusion(matrix, 0.75, order="neighbors",
                                        reordering=reordering)
    assert result.clusters == butina.clusters
    assert result.labels.tolist() == butina.labels.tolist()


def test_positions_refer_to_the_callers_items_after_normalization_drops_one():
    mols = _mols(["O"] + DESCRIPTOR_SMILES)
    matrix = oecluster.pdist(mols, "descriptor")
    by_name = oecluster.sphere_exclusion(mols, 1.0, comparison="descriptor")
    by_matrix = oecluster.sphere_exclusion(matrix, 1.0)

    assert by_name.excluded == [[0, "missing-descriptor"]]
    assert by_name.num_samples == 9
    assert by_name.labels[0] == -1
    assert by_name.labels[1:].tolist() == by_matrix.labels.tolist()
    assert by_name.clusters == tuple(
        tuple(i + 1 for i in cluster) for cluster in by_matrix.clusters)
    assert by_name.centers == tuple(i + 1 for i in by_matrix.centers)

    # A permutation names every caller position, dropped ones included;
    # the dropped position is skipped when it is mapped to native indices.
    reversed_name = oecluster.sphere_exclusion(
        mols, 1.0, comparison="descriptor", order=list(reversed(range(9))))
    reversed_matrix = oecluster.sphere_exclusion(
        matrix, 1.0, order=list(reversed(range(8))))
    assert reversed_name.clusters == tuple(
        tuple(i + 1 for i in cluster) for cluster in reversed_matrix.clusters)
    assert reversed_name.labels[0] == -1


def test_the_result_is_a_clustering_result():
    result = oecluster.sphere_exclusion(_line_matrix(), 1.5)
    assert isinstance(result, oecluster.ClusteringResult)
    assert isinstance(result, oecluster.SphereExclusionResult)
    assert len(result) == 3
    assert list(result) == [(0, 1), (2,), (3, 4)]
    assert result[2] == (3, 4)
    assert repr(result) == (
        "SphereExclusionResult(num_clusters=3, num_samples=5)")


def test_excluded_is_a_copy():
    mols = _mols(["O"] + DESCRIPTOR_SMILES)
    result = oecluster.sphere_exclusion(mols, 1.0, comparison="descriptor")
    result.excluded.append([99, "tampered"])
    result.excluded[0][1] = "tampered"
    assert result.excluded == [[0, "missing-descriptor"]]


def test_empty_matrices_and_comparisons_give_an_empty_result():
    for items in (oecluster.pdist([], "fingerprint"),
                  oecluster.FingerprintComparison([])):
        result = oecluster.sphere_exclusion(items, 0.5)
        assert result.num_samples == 0
        assert result.clusters == ()
        assert result.centers == ()


def test_threading_options_are_forwarded(monkeypatch):
    native = oecluster.oecluster
    real = native.sphere_exclusion
    seen = []

    def spy(target, options):
        seen.append((options.num_threads, options.chunk_size))
        return real(target, options)

    monkeypatch.setattr(native, "sphere_exclusion", spy)
    result = oecluster.sphere_exclusion(_mols(FP_SMILES), 0.75,
                                        comparison="fingerprint",
                                        num_threads=3, chunk_size=2)
    assert seen == [(3, 2)]
    assert list(result.centers) == [0, 3, 5, 6, 10]


@pytest.mark.parametrize(("kwargs", "match"), [
    ({"threshold": -1.0}, "threshold must be non-negative"),
    ({"threshold": math.nan}, "not NaN"),
    ({"threshold": math.inf}, "threshold must be finite"),
    ({"similarity": True}, "similarity=True"),
    ({"order": "bogus"}, "Unknown sphere_exclusion order"),
    ({"order": [0, 1, 2, 3]}, "4 entries for 5 items"),
    ({"order": [0, 1, 2, 3, 5]}, "outside the item range"),
    ({"order": [-1, 1, 2, 3, 4]}, "outside the item range"),
    ({"order": [0, 1, 1, 3, 4]}, "repeats position 1"),
    ({"reordering": True}, "reordering requires order='neighbors'"),
    ({"reordering": True, "order": [0, 1, 2, 3, 4]},
     "reordering requires order='neighbors'"),
    ({"assignment": "bogus"}, "Unknown sphere_exclusion assignment"),
    ({"assignment": None}, "Unknown sphere_exclusion assignment"),
    ({"num_threads": -1}, "num_threads must be at least 0"),
    ({"chunk_size": 0}, "chunk_size must be at least 1"),
    ({"chunk_size": 2.5}, "chunk_size must be an integer"),
])
def test_invalid_arguments_are_value_errors(kwargs, match):
    kwargs = {"threshold": 1.5, **kwargs}
    threshold = kwargs.pop("threshold")
    with pytest.raises(ValueError, match=match):
        oecluster.sphere_exclusion(_line_matrix(), threshold, **kwargs)


@pytest.mark.parametrize("kwargs", [
    {"order": [0, True, 2, 3, 4]},
    {"order": [0, 1.0, 2, 3, 4]},
    {"order": None},
    {"order": 3},
    {"reordering": "yes", "order": "neighbors"},
])
def test_ill_typed_arguments_are_type_errors(kwargs):
    with pytest.raises(TypeError):
        oecluster.sphere_exclusion(_line_matrix(), 1.5, **kwargs)


@pytest.mark.parametrize("assignment", ["first", "nearest"])
@pytest.mark.parametrize("reordering", [False, True])
def test_the_neighbors_order_runs_on_the_lazy_paths(reordering, assignment):
    mols = _mols(FP_SMILES)
    prebuilt = oecluster.FingerprintComparison(mols)
    expected = oecluster.sphere_exclusion(
        _compare_matrix(prebuilt), 0.75, order="neighbors",
        reordering=reordering, assignment=assignment)
    for items, extra in ((prebuilt, {}),
                         (mols, {"comparison": "fingerprint"})):
        result = oecluster.sphere_exclusion(
            items, 0.75, order="neighbors", reordering=reordering,
            assignment=assignment, num_threads=4, chunk_size=3, **extra)
        assert result.clusters == expected.clusters
        assert result.labels.tolist() == expected.labels.tolist()
        assert result.centers == expected.centers


def test_the_lazy_neighbors_order_maps_positions_around_a_dropped_item():
    mols = _mols(["O", *DESCRIPTOR_SMILES])
    by_name = oecluster.sphere_exclusion(mols, 1.0, comparison="descriptor",
                                         order="neighbors")
    by_matrix = oecluster.sphere_exclusion(oecluster.pdist(mols, "descriptor"),
                                           1.0, order="neighbors")
    assert by_name.excluded == [[0, "missing-descriptor"]]
    assert by_name.labels[0] == -1
    assert by_name.labels[1:].tolist() == by_matrix.labels.tolist()
    assert by_name.clusters == tuple(
        tuple(i + 1 for i in cluster) for cluster in by_matrix.clusters)


def test_a_budget_reaches_the_lazy_neighbors_graph(monkeypatch):
    native = oecluster.oecluster
    real = native.sphere_exclusion
    seen = []

    def spy(target, options):
        seen.append(options.max_graph_bytes)
        return real(target, options)

    monkeypatch.setattr(native, "sphere_exclusion", spy)
    mols = _mols(FP_SMILES)
    oecluster.sphere_exclusion(mols, 0.75, comparison="fingerprint",
                               order="neighbors", max_graph_bytes=1 << 30)
    oecluster.sphere_exclusion(mols, 0.75, comparison="fingerprint",
                               order="neighbors")
    assert seen == [1 << 30, 0]


def test_a_budget_needs_a_comparison_and_the_neighbors_order():
    mols = _mols(FP_SMILES)
    with pytest.raises(TypeError, match="no graph budget"):
        oecluster.sphere_exclusion(oecluster.pdist(mols, "fingerprint"), 0.75,
                                   order="neighbors", max_graph_bytes=1 << 30)
    for order in ("input", list(range(len(mols)))):
        with pytest.raises(TypeError,
                           match="only order='neighbors' builds one"):
            oecluster.sphere_exclusion(mols, 0.75, comparison="fingerprint",
                                       order=order, max_graph_bytes=1 << 30)
    with pytest.raises(TypeError, match="not a bool"):
        oecluster.sphere_exclusion(mols, 0.75, comparison="fingerprint",
                                   order="neighbors", max_graph_bytes=True)
    with pytest.raises(ValueError, match="must be positive, got 0"):
        oecluster.sphere_exclusion(mols, 0.75, comparison="fingerprint",
                                   order="neighbors", max_graph_bytes=0)


def test_the_lazy_neighbors_order_refuses_rocs():
    # The neighbor order's graph needs repeatable scores, which ROCS does not
    # give; refused by name before the comparison is built.
    with pytest.raises(ValueError,
                       match="cannot cluster a ROCS comparison without a "
                             "matrix"):
        oecluster.sphere_exclusion(_mols(FP_SMILES), 0.5, comparison="rocs",
                                   order="neighbors")


def test_an_oversized_lazy_neighbors_graph_is_a_memory_error():
    mols = _mols(FP_SMILES)
    with pytest.raises(MemoryError,
                       match=r"sphere_exclusion would build a threshold graph "
                             r"of \d+ bytes for 12 items"):
        oecluster.sphere_exclusion(mols, 0.9, comparison="fingerprint",
                                   order="neighbors", max_graph_bytes=64)


def test_a_non_finite_matrix_entry_is_refused_by_the_gate():
    square = [[0.0, 1.0, 2.0], [1.0, 0.0, math.nan], [2.0, math.nan, 0.0]]
    with pytest.raises(ValueError, match="non-finite"):
        oecluster.sphere_exclusion(_dense_distance_matrix(square), 0.5)


def test_sparse_storage_is_refused():
    storage = oecluster.SparseStorage(4, 0.5)
    matrix = oecluster.SymmetricDistanceMatrix(storage, "test", list("abcd"),
                                               {})
    with pytest.raises(ValueError, match="SparseStorage is not supported"):
        oecluster.sphere_exclusion(matrix, 0.5)


@pytest.mark.parametrize(("missing", "match"), [
    ("propagate", "non-finite"),
    ("ignore", "per-pair feature subsets"),
])
def test_descriptor_missingness_that_cannot_be_ranked_is_refused(missing,
                                                                 match):
    options = {"missing": missing}
    if missing == "ignore":
        options["metric"] = "euclidean"
    mols = _mols(DESCRIPTOR_SMILES)
    with pytest.raises(ValueError, match=match):
        oecluster.sphere_exclusion(mols, 1.0, comparison="descriptor",
                                   **options)
    with pytest.raises(ValueError, match=match):
        oecluster.sphere_exclusion(
            oecluster.DescriptorComparison(mols, **options), 1.0)


def test_a_similarity_comparison_is_refused():
    prebuilt = oecluster.FingerprintComparison(_mols(FP_SMILES),
                                               similarity=True)
    with pytest.raises(ValueError, match="requires distances"):
        oecluster.sphere_exclusion(prebuilt, 0.5)


@pytest.mark.parametrize("call", [
    lambda mols: oecluster.sphere_exclusion(
        oecluster.pdist(mols, "fingerprint"), 0.5, comparison="fingerprint"),
    lambda mols: oecluster.sphere_exclusion(
        oecluster.pdist(mols, "fingerprint"), 0.5, radius=1),
    lambda mols: oecluster.sphere_exclusion(mols, 0.5),
    lambda mols: oecluster.sphere_exclusion(mols, 0.5, comparison=3),
    lambda mols: oecluster.sphere_exclusion(
        oecluster.cdist(mols[:2], mols[2:4], "fingerprint"), 0.5),
])
def test_arguments_that_fit_no_path_are_type_errors(call):
    with pytest.raises(TypeError):
        call(_mols(FP_SMILES))


def test_an_empty_sequence_is_refused():
    with pytest.raises(ValueError, match="requires at least one item"):
        oecluster.sphere_exclusion([], 0.5, comparison="fingerprint")
    with pytest.raises(ValueError, match="requires at least one item"):
        oecluster.sphere_exclusion(_mols(["O"]), 0.5, comparison="descriptor")


def test_the_sphere_exclusion_surface_is_exported():
    exported = ("sphere_exclusion", "SphereExclusionResult")
    missing = [name for name in exported if name not in oecluster.__all__]
    assert missing == []
    assert all(hasattr(oecluster, name) for name in exported)
