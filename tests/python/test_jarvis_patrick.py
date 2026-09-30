"""Python surface of jarvis_patrick: graph and raw input, results, refusals."""

import math

import numpy as np
import oecluster
import pytest
from openeye import oechem

FP_SMILES = ["CCO", "CCCO", "CCCCO", "c1ccccc1", "Cc1ccccc1", "CCc1ccccc1",
             "CC(=O)O", "CC(=O)OC", "CCN", "CCCN", "C1CCCCC1", "c1ccncc1"]

DESCRIPTOR_SMILES = ["CCO", "CCC", "CCCC", "c1ccccc1", "CCN", "CCOC",
                     "CC(=O)O", "CCCCCC"]

# With k = 2: 0, 1 and 2 are pairwise mutual and share one item; 4 and 5 are
# mutual and share 3; 3 names 2 and 1 but neither names it back.
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


def _oracle_labels(square, k, kmin):
    """Classic Jarvis-Patrick by brute force over a NumPy neighbor table."""
    n = len(square)
    rows = []
    for i in range(n):
        ranked = sorted((square[i][j], j) for j in range(n) if j != i)[:k]
        rows.append({j for _, j in ranked})
    component = list(range(n))
    for i in range(n):
        for j in range(i + 1, n):
            if j in rows[i] and i in rows[j] and len(rows[i] & rows[j]) >= kmin:
                old, new = max(component[i], component[j]), min(component[i],
                                                                component[j])
                component = [new if c == old else c for c in component]
    order = {root: label for label, root in enumerate(sorted(set(component)))}
    return [order[c] for c in component]


def test_a_hand_worked_example():
    result = oecluster.jarvis_patrick(
        _dense_distance_matrix(_square(_POINTS)), k=2, kmin=1)
    assert result.clusters == ((0, 1, 2), (3,), (4, 5))
    assert result.labels.tolist() == [0, 0, 0, 1, 2, 2]
    assert result.k == 2
    assert result.kmin == 1
    assert result.method == "jarvis_patrick"
    assert result.excluded == []


def test_every_path_matches_the_oracle():
    rng = np.random.default_rng(3)
    points = rng.integers(0, 12, size=14).astype(float).tolist()
    square = _square(points)
    for k, kmin in ((3, 1), (4, 2), (5, 3)):
        expected = _oracle_labels(square, k, kmin)
        dense = oecluster.jarvis_patrick(_dense_distance_matrix(square), k=k,
                                         kmin=kmin)
        sparse = oecluster.jarvis_patrick(
            _sparse_distance_matrix(square, 100.0), k=k, kmin=kmin)
        graph = oecluster.jarvis_patrick(
            oecluster.knn_graph(_dense_distance_matrix(square), k), kmin=kmin)
        for result in (dense, sparse, graph):
            assert result.labels.tolist() == expected


def test_the_three_call_styles_and_the_graph_agree():
    mols = _mols(FP_SMILES)
    by_matrix = oecluster.jarvis_patrick(oecluster.pdist(mols, "fingerprint"),
                                         k=3, kmin=1)
    by_name = oecluster.jarvis_patrick(mols, k=3, kmin=1,
                                       comparison="fingerprint")
    prebuilt = oecluster.jarvis_patrick(oecluster.FingerprintComparison(mols),
                                        k=3, kmin=1)
    graph = oecluster.jarvis_patrick(
        oecluster.knn_graph(mols, 3, comparison="fingerprint"), kmin=1)
    for other in (by_name, prebuilt, graph):
        assert other.clusters == by_matrix.clusters
        assert other.labels.tolist() == by_matrix.labels.tolist()


def test_a_graph_ignores_valid_threading_options():
    graph = oecluster.knn_graph(_dense_distance_matrix(_square(_POINTS)), 2)
    expected = oecluster.jarvis_patrick(graph, kmin=1)
    other = oecluster.jarvis_patrick(graph, kmin=1, k=2, num_threads=3,
                                     chunk_size=1)
    assert other.clusters == expected.clusters


def test_positions_refer_to_the_callers_items_after_normalization_drops_one():
    mols = _mols(["O", *DESCRIPTOR_SMILES])
    by_matrix = oecluster.jarvis_patrick(oecluster.pdist(mols, "descriptor"),
                                         k=2, kmin=0)
    for by_name in (
            oecluster.jarvis_patrick(mols, k=2, kmin=0,
                                     comparison="descriptor"),
            oecluster.jarvis_patrick(
                oecluster.knn_graph(mols, 2, comparison="descriptor"),
                kmin=0)):
        assert by_name.excluded == [[0, "missing-descriptor"]]
        assert by_name.num_samples == 9
        assert by_name.labels[0] == -1
        assert by_name.labels[1:].tolist() == by_matrix.labels.tolist()
        assert by_name.clusters == tuple(
            tuple(i + 1 for i in cluster) for cluster in by_matrix.clusters)


def test_empty_input_gives_an_empty_result():
    for items in (oecluster.pdist([], "fingerprint"),
                  oecluster.FingerprintComparison([])):
        result = oecluster.jarvis_patrick(items, k=3, kmin=7)
        assert result.num_samples == 0
        assert result.clusters == ()
        graph = oecluster.knn_graph(items, 3)
        assert oecluster.jarvis_patrick(graph, kmin=9).clusters == ()


@pytest.mark.parametrize(("kwargs", "match"), [
    ({"k": 2, "kmin": 2}, "kmin must be less than k = 2, got 2"),
    ({"k": 2, "kmin": 5}, "kmin must be less than k = 2, got 5"),
    ({"k": 6, "kmin": 1}, "k must be between 1 and 5 for 6 items, got 6"),
    ({"k": 0, "kmin": 0}, "k must be between 1 and 5 for 6 items, got 0"),
    ({"k": 2.0, "kmin": 1}, "k must be an integer"),
    ({"k": 2, "kmin": 1.5}, "kmin must be an integer"),
    ({"k": 2, "kmin": True}, "kmin must be an integer"),
    ({"k": 2, "kmin": -1}, "kmin must be at least 0"),
    ({"k": 2, "kmin": 1, "num_threads": -1}, "num_threads must be at least 0"),
    ({"k": 2, "kmin": 1, "chunk_size": 0}, "chunk_size must be at least 1"),
    ({"k": 2, "kmin": 1, "similarity": True},
     "similarity=True is not supported"),
])
def test_invalid_options_are_value_errors(kwargs, match):
    matrix = _dense_distance_matrix(_square(_POINTS))
    with pytest.raises(ValueError, match=match):
        oecluster.jarvis_patrick(matrix, **kwargs)


def test_graph_specific_refusals():
    graph = oecluster.knn_graph(_dense_distance_matrix(_square(_POINTS)), 2)
    with pytest.raises(ValueError, match="k=3 does not match the graph's k=2"):
        oecluster.jarvis_patrick(graph, k=3, kmin=1)
    with pytest.raises(ValueError, match="kmin must be less than k = 2"):
        oecluster.jarvis_patrick(graph, kmin=2)
    with pytest.raises(ValueError, match="chunk_size must be at least 1"):
        oecluster.jarvis_patrick(graph, kmin=1, chunk_size=0)
    with pytest.raises(TypeError, match="already fixes the neighbors"):
        oecluster.jarvis_patrick(graph, kmin=1, comparison="fingerprint")
    with pytest.raises(TypeError, match="already fixes the neighbors"):
        oecluster.jarvis_patrick(graph, kmin=1, radius=2)


def test_raw_input_requires_k():
    with pytest.raises(TypeError, match="requires k= unless items is a KNNGraph"):
        oecluster.jarvis_patrick(_dense_distance_matrix(_square(_POINTS)),
                                 kmin=1)


def test_k_and_kmin_are_keyword_only():
    with pytest.raises(TypeError):
        oecluster.jarvis_patrick(_dense_distance_matrix(_square(_POINTS)), 2, 1)


def test_a_single_item_is_refused():
    with pytest.raises(ValueError, match="a single item has no neighbors"):
        oecluster.jarvis_patrick(_dense_distance_matrix([[0.0]]), k=1, kmin=0)


def test_an_empty_sequence_is_refused():
    with pytest.raises(ValueError, match="requires at least one item"):
        oecluster.jarvis_patrick([], k=1, kmin=0, comparison="fingerprint")


def test_a_sparse_short_item_is_a_runtime_error():
    matrix = _sparse_distance_matrix(_square(_POINTS), 1.0)
    with pytest.raises(RuntimeError, match="raise the cutoff or lower k"):
        oecluster.jarvis_patrick(matrix, k=2, kmin=1)


def test_matrices_the_gate_refuses_are_value_errors():
    similarity = oecluster.pdist(_mols(FP_SMILES), "fingerprint",
                                 similarity=True)
    with pytest.raises(ValueError, match="holds a similarity"):
        oecluster.jarvis_patrick(similarity, k=2, kmin=1)
    stamped = _sparse_distance_matrix(_square(_POINTS), 100.0)
    stamped._facts['is_distance'] = False
    with pytest.raises(ValueError, match="holds a similarity"):
        oecluster.jarvis_patrick(stamped, k=2, kmin=1)
    square = _square(_POINTS)
    square[1][3] = square[3][1] = math.nan
    with pytest.raises(ValueError, match="non-finite"):
        oecluster.jarvis_patrick(_dense_distance_matrix(square), k=2, kmin=1)
    with pytest.raises(ValueError, match="non-finite"):
        oecluster.jarvis_patrick(_sparse_distance_matrix(square, 100.0), k=2,
                                 kmin=1)


def test_a_similarity_comparison_is_refused():
    prebuilt = oecluster.FingerprintComparison(_mols(FP_SMILES),
                                               similarity=True)
    with pytest.raises(ValueError, match="requires distances"):
        oecluster.jarvis_patrick(prebuilt, k=2, kmin=1)


@pytest.mark.parametrize("call", [
    lambda mols: oecluster.jarvis_patrick(
        oecluster.pdist(mols, "fingerprint"), k=2, kmin=1,
        comparison="fingerprint"),
    lambda mols: oecluster.jarvis_patrick(mols, k=2, kmin=1),
    lambda mols: oecluster.jarvis_patrick(
        oecluster.cdist(mols[:2], mols[2:4], "fingerprint"), k=1, kmin=0),
])
def test_arguments_that_fit_no_path_are_type_errors(call):
    with pytest.raises(TypeError):
        call(_mols(FP_SMILES))


def test_an_empty_similarity_matrix_returns_an_empty_result():
    matrix = oecluster.pdist([], "fingerprint", similarity=True)
    result = oecluster.jarvis_patrick(matrix, k=7, kmin=3)
    assert result.num_samples == 0
    assert result.clusters == ()
    assert result.k == 7
    assert result.kmin == 3


def test_an_empty_similarity_comparison_returns_an_empty_result():
    comparison = oecluster.FingerprintComparison([], similarity=True)
    result = oecluster.jarvis_patrick(comparison, k=7, kmin=3)
    assert result.num_samples == 0
    assert result.clusters == ()
    assert result.k == 7
    assert result.kmin == 3


def test_non_empty_similarity_input_is_still_refused():
    matrix = oecluster.pdist(_mols(FP_SMILES), "fingerprint", similarity=True)
    with pytest.raises(ValueError, match="holds a similarity"):
        oecluster.jarvis_patrick(matrix, k=2, kmin=1)


def test_the_jarvis_patrick_surface_is_exported():
    assert "jarvis_patrick" in oecluster.__all__
    assert "JarvisPatrickResult" in oecluster.__all__
    assert issubclass(oecluster.JarvisPatrickResult,
                      oecluster.ClusteringResult)
