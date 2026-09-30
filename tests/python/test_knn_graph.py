"""Python surface of knn_graph: paths, KNNGraph members, refusals."""

import math

import numpy as np
import oecluster
import pytest
from openeye import oechem

FP_SMILES = ["CCO", "CCCO", "CCCCO", "c1ccccc1", "Cc1ccccc1", "CCc1ccccc1",
             "CC(=O)O", "CC(=O)OC", "CCN", "CCCN", "C1CCCCC1", "c1ccncc1"]

# Prefixing water makes the descriptor complete-case mask drop position 0.
DESCRIPTOR_SMILES = ["CCO", "CCC", "CCCC", "c1ccccc1", "CCN", "CCOC",
                     "CC(=O)O", "CCCCCC"]

# Six points on a line; with k = 2 row 2 ties items 0 and 3 at distance 2.
_POINTS = (0.0, 1.0, 2.0, 4.0, 7.0, 8.0)


def _mols(smiles_list):
    mols = []
    for idx, smi in enumerate(smiles_list):
        mol = oechem.OEGraphMol()
        oechem.OESmilesToMol(mol, smi)
        mol.SetTitle(f"mol{idx}")
        mols.append(mol)
    return mols


def _square(points):
    return [[abs(a - b) for b in points] for a in points]


def _dense_distance_matrix(square):
    storage = oecluster.DenseStorage(len(square))
    for i in range(len(square)):
        for j in range(i + 1, len(square)):
            storage.Set(i, j, float(square[i][j]))
    return oecluster.SymmetricDistanceMatrix(
        storage, "test", [f"item_{i}" for i in range(len(square))], {})


def _sparse_distance_matrix(square, cutoff):
    storage = oecluster.SparseStorage(len(square), cutoff)
    for i in range(len(square)):
        for j in range(i + 1, len(square)):
            storage.Set(i, j, float(square[i][j]))
    storage.Finalize()
    return oecluster.SymmetricDistanceMatrix(
        storage, "test", [f"item_{i}" for i in range(len(square))], {})


def _oracle(square, k):
    """Every other item sorted by (distance, index); the first k kept."""
    indices, distances = [], []
    for i, row in enumerate(square):
        ranked = sorted((row[j], j) for j in range(len(row)) if j != i)[:k]
        indices.append([j for _, j in ranked])
        distances.append([d for d, _ in ranked])
    return (np.array(indices, dtype=np.int64).reshape(len(square), k),
            np.array(distances, dtype=np.float64).reshape(len(square), k))


def test_a_dense_matrix_matches_the_oracle():
    square = _square(_POINTS)
    graph = oecluster.knn_graph(_dense_distance_matrix(square), 2)
    indices, distances = _oracle(square, 2)
    np.testing.assert_array_equal(graph.indices, indices)
    np.testing.assert_array_equal(graph.distances, distances)
    assert graph.indices.tolist() == [[1, 2], [0, 2], [1, 0], [2, 1], [5, 3],
                                      [4, 3]]
    assert graph.k == 2


def test_sparse_storage_matches_the_oracle():
    square = _square(_POINTS)
    for cutoff in (100.0, 4.0):
        graph = oecluster.knn_graph(_sparse_distance_matrix(square, cutoff), 2)
        indices, distances = _oracle(square, 2)
        np.testing.assert_array_equal(graph.indices, indices)
        np.testing.assert_array_equal(graph.distances, distances)


def test_a_sparse_short_item_is_a_runtime_error():
    matrix = _sparse_distance_matrix(_square(_POINTS), 1.0)
    with pytest.raises(RuntimeError, match="raise the cutoff or lower k"):
        oecluster.knn_graph(matrix, 2)


def test_the_three_call_styles_agree():
    mols = _mols(FP_SMILES)
    by_matrix = oecluster.knn_graph(oecluster.pdist(mols, "fingerprint"), 3)
    by_name = oecluster.knn_graph(mols, 3, comparison="fingerprint")
    prebuilt = oecluster.knn_graph(oecluster.FingerprintComparison(mols), 3)
    for other in (by_name, prebuilt):
        np.testing.assert_array_equal(other.indices, by_matrix.indices)
        np.testing.assert_allclose(other.distances, by_matrix.distances,
                                   rtol=0, atol=1e-12)
        assert other.excluded == []


def test_threading_options_do_not_change_the_graph():
    mols = _mols(FP_SMILES)
    expected = oecluster.knn_graph(mols, 3, comparison="fingerprint")
    for num_threads, chunk_size in ((1, 1), (4, 7), (2, 2**64 - 1)):
        graph = oecluster.knn_graph(mols, 3, comparison="fingerprint",
                                    num_threads=num_threads,
                                    chunk_size=chunk_size)
        np.testing.assert_array_equal(graph.indices, expected.indices)


def test_rows_follow_normalization_and_indices_are_caller_positions():
    mols = _mols(["O", *DESCRIPTOR_SMILES])
    by_name = oecluster.knn_graph(mols, 2, comparison="descriptor")
    by_matrix = oecluster.knn_graph(oecluster.pdist(mols, "descriptor"), 2)
    assert by_name.excluded == [[0, "missing-descriptor"]]
    assert by_name.positions.tolist() == list(range(1, 9))
    assert len(by_name) == 8
    np.testing.assert_array_equal(by_name.indices, by_matrix.indices + 1)
    assert 0 not in by_name.indices
    assert by_matrix.positions.tolist() == list(range(8))


def test_arrays_have_rows_by_k_shape_and_dtype():
    graph = oecluster.knn_graph(_dense_distance_matrix(_square(_POINTS)), 3)
    assert graph.indices.shape == (6, 3)
    assert graph.indices.dtype == np.int64
    assert graph.distances.shape == (6, 3)
    assert graph.distances.dtype == np.float64
    assert graph.positions.shape == (6,)
    assert graph.positions.dtype == np.int64
    assert len(graph) == 6
    assert repr(graph) == "KNNGraph(num_rows=6, k=3)"


def test_arrays_and_excluded_are_copies():
    mols = _mols(["O", *DESCRIPTOR_SMILES])
    graph = oecluster.knn_graph(mols, 2, comparison="descriptor")
    graph.indices[0, 0] = -99
    graph.distances[0, 0] = -99.0
    graph.positions[0] = -99
    graph.excluded.append([99, "tampered"])
    assert graph.indices[0, 0] != -99
    assert graph.distances[0, 0] != -99.0
    assert graph.positions[0] == 1
    assert graph.excluded == [[0, "missing-descriptor"]]


def test_empty_matrices_and_comparisons_give_an_empty_graph():
    for items in (oecluster.pdist([], "fingerprint"),
                  oecluster.FingerprintComparison([])):
        graph = oecluster.knn_graph(items, 3)
        assert len(graph) == 0
        assert graph.k == 3
        assert graph.indices.shape == (0, 3)


def test_a_k_no_numpy_array_can_hold_is_refused_even_without_items():
    empty = oecluster.pdist([], "fingerprint")
    assert oecluster.knn_graph(empty, 10**6).indices.shape == (0, 10**6)
    with pytest.raises(ValueError, match="largest NumPy array dimension"):
        oecluster.knn_graph(empty, 2**63)


def test_the_constructor_is_private():
    with pytest.raises(TypeError, match="knn_graph"):
        oecluster.KNNGraph()


@pytest.mark.parametrize(("kwargs", "match"), [
    ({"k": 0}, "k must be between 1 and 5 for 6 items, got 0"),
    ({"k": 6}, "k must be between 1 and 5 for 6 items, got 6"),
    ({"k": 2.5}, "k must be an integer"),
    ({"k": True}, "k must be an integer"),
    ({"k": -1}, "k must be at least 0"),
    ({"k": 2, "num_threads": -1}, "num_threads must be at least 0"),
    ({"k": 2, "chunk_size": 0}, "chunk_size must be at least 1"),
    ({"k": 2, "similarity": True}, "similarity=True is not supported"),
])
def test_invalid_options_are_value_errors(kwargs, match):
    matrix = _dense_distance_matrix(_square(_POINTS))
    with pytest.raises(ValueError, match=match):
        oecluster.knn_graph(matrix, **kwargs)


def test_a_single_item_is_refused():
    with pytest.raises(ValueError, match="a single item has no neighbors"):
        oecluster.knn_graph(_dense_distance_matrix([[0.0]]), 1)


@pytest.mark.parametrize("call", [
    lambda mols: oecluster.knn_graph(
        oecluster.pdist(mols, "fingerprint"), 2, comparison="fingerprint"),
    lambda mols: oecluster.knn_graph(oecluster.pdist(mols, "fingerprint"), 2,
                                     radius=1),
    lambda mols: oecluster.knn_graph(mols, 2),
    lambda mols: oecluster.knn_graph(
        oecluster.cdist(mols[:2], mols[2:4], "fingerprint"), 1),
])
def test_arguments_that_fit_no_path_are_type_errors(call):
    with pytest.raises(TypeError):
        call(_mols(FP_SMILES))


def test_an_empty_sequence_is_refused():
    with pytest.raises(ValueError, match="requires at least one item"):
        oecluster.knn_graph([], 1, comparison="fingerprint")


def test_matrices_the_gate_refuses_are_value_errors():
    similarity = oecluster.pdist(_mols(FP_SMILES), "fingerprint",
                                 similarity=True)
    with pytest.raises(ValueError, match="holds a similarity"):
        oecluster.knn_graph(similarity, 2)
    stamped = _sparse_distance_matrix(_square(_POINTS), 100.0)
    stamped._facts['is_distance'] = False
    with pytest.raises(ValueError, match="holds a similarity"):
        oecluster.knn_graph(stamped, 2)
    square = _square(_POINTS)
    square[1][3] = square[3][1] = math.nan
    with pytest.raises(ValueError, match="non-finite"):
        oecluster.knn_graph(_dense_distance_matrix(square), 2)
    with pytest.raises(ValueError, match="non-finite"):
        oecluster.knn_graph(_sparse_distance_matrix(square, 100.0), 2)


def test_a_similarity_comparison_is_refused():
    prebuilt = oecluster.FingerprintComparison(_mols(FP_SMILES),
                                               similarity=True)
    with pytest.raises(ValueError, match="requires distances"):
        oecluster.knn_graph(prebuilt, 2)


def test_other_diversity_callers_still_refuse_sparse_storage():
    matrix = _sparse_distance_matrix(_square(_POINTS), 100.0)
    with pytest.raises(ValueError, match="SparseStorage is not supported"):
        oecluster.sphere_exclusion(matrix, 0.5)


def test_the_knn_graph_surface_is_exported():
    assert "knn_graph" in oecluster.__all__
    assert "KNNGraph" in oecluster.__all__
    assert oecluster.KNNGraph.__module__ == "oecluster"
