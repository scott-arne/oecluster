"""Python surface of leiden: graph and raw input, results, refusals."""

import math

import numpy as np
import oecluster
import pytest
from openeye import oechem

FP_SMILES = ["CCO", "CCCO", "CCCCO", "c1ccccc1", "Cc1ccccc1", "CCc1ccccc1",
             "CC(=O)O", "CC(=O)OC", "CCN", "CCCN", "C1CCCCC1", "c1ccncc1"]

DESCRIPTOR_SMILES = ["CCO", "CCC", "CCCC", "c1ccccc1", "CCN", "CCOC",
                     "CC(=O)O", "CCCCCC"]

# With k = 3 each block of four is its own neighborhood: every row names the
# other three members, so every N+ is the block, each in-block edge weighs
# 4 / (8 - 4) = 1, and no edge crosses the gap. Each block has internal
# weight 6 and strength 12 of 2m = 24, so modularity is
# 2 * (6 / 12 - (12 / 24)**2) = 0.5, and CPM at resolution 0.5 is
# 2 * (6 - 0.5 * 4 * 3 / 2) = 6. Singletons, each of strength 3, have
# modularity -8 * (3 / 24)**2 = -0.125.
_POINTS = (0.0, 1.0, 2.0, 3.0, 10.0, 11.0, 12.0, 13.0)
_BLOCKS = [0, 0, 0, 0, 1, 1, 1, 1]


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


def test_a_hand_worked_example():
    result = oecluster.leiden(_dense_distance_matrix(_square(_POINTS)), k=3)
    assert result.labels.tolist() == _BLOCKS
    assert result.clusters == ((0, 1, 2, 3), (4, 5, 6, 7))
    assert result.quality == 0.5
    assert result.iterations == 2
    assert result.objective == "modularity"
    assert result.resolution == 1.0
    assert result.k == 3
    assert result.method == "leiden"
    assert result.excluded == []


def test_cpm_reports_its_own_quality():
    result = oecluster.leiden(_dense_distance_matrix(_square(_POINTS)), k=3,
                              objective="cpm", resolution=0.5)
    assert result.labels.tolist() == _BLOCKS
    assert result.quality == 6.0
    assert result.objective == "cpm"
    assert result.resolution == 0.5


def test_zero_iterations_return_singletons():
    result = oecluster.leiden(_dense_distance_matrix(_square(_POINTS)), k=3,
                              n_iterations=0)
    assert result.labels.tolist() == list(range(8))
    assert result.iterations == 0
    assert result.quality == -0.125


def test_every_input_path_agrees():
    square = _square(_POINTS)
    dense = oecluster.leiden(_dense_distance_matrix(square), k=3)
    sparse = oecluster.leiden(_sparse_distance_matrix(square, 100.0), k=3)
    graph = oecluster.leiden(
        oecluster.knn_graph(_dense_distance_matrix(square), 3))
    for result in (dense, sparse, graph):
        assert result.labels.tolist() == _BLOCKS
        assert result.quality == 0.5


def test_the_three_call_styles_and_the_graph_agree():
    mols = _mols(FP_SMILES)
    by_matrix = oecluster.leiden(oecluster.pdist(mols, "fingerprint"), k=3)
    by_name = oecluster.leiden(mols, k=3, comparison="fingerprint")
    prebuilt = oecluster.leiden(oecluster.FingerprintComparison(mols), k=3)
    graph = oecluster.leiden(
        oecluster.knn_graph(mols, 3, comparison="fingerprint"))
    for other in (by_name, prebuilt, graph):
        assert other.clusters == by_matrix.clusters
        assert other.labels.tolist() == by_matrix.labels.tolist()
        assert other.quality == by_matrix.quality
        assert other.iterations == by_matrix.iterations


def test_the_same_seed_gives_the_same_result_and_threads_do_not_matter():
    rng = np.random.default_rng(7)
    points = rng.normal(size=40).tolist()
    matrix = _dense_distance_matrix(_square(points))
    first = oecluster.leiden(matrix, k=6, seed=11, num_threads=1)
    for other in (oecluster.leiden(matrix, k=6, seed=11, num_threads=1),
                  oecluster.leiden(matrix, k=6, seed=11, num_threads=4,
                                   chunk_size=7)):
        assert other.labels.tolist() == first.labels.tolist()
        assert other.quality == first.quality
        assert other.iterations == first.iterations


def test_a_graph_ignores_valid_threading_options():
    graph = oecluster.knn_graph(_dense_distance_matrix(_square(_POINTS)), 3)
    expected = oecluster.leiden(graph)
    other = oecluster.leiden(graph, k=3, num_threads=3, chunk_size=1)
    assert other.clusters == expected.clusters
    assert other.k == 3


def test_positions_refer_to_the_callers_items_after_normalization_drops_one():
    mols = _mols(["O", *DESCRIPTOR_SMILES])
    by_matrix = oecluster.leiden(oecluster.pdist(mols, "descriptor"), k=2)
    for by_name in (
            oecluster.leiden(mols, k=2, comparison="descriptor"),
            oecluster.leiden(
                oecluster.knn_graph(mols, 2, comparison="descriptor"))):
        assert by_name.excluded == [[0, "missing-descriptor"]]
        assert by_name.num_samples == 9
        assert by_name.labels[0] == -1
        assert by_name.labels[1:].tolist() == by_matrix.labels.tolist()
        assert by_name.clusters == tuple(
            tuple(i + 1 for i in cluster) for cluster in by_matrix.clusters)


def test_empty_input_gives_an_empty_result():
    for items in (oecluster.pdist([], "fingerprint"),
                  oecluster.FingerprintComparison([])):
        result = oecluster.leiden(items, k=3, objective="cpm",
                                  resolution=0.25)
        assert result.num_samples == 0
        assert result.clusters == ()
        assert result.quality == 0.0
        assert result.iterations == 0
        assert result.objective == "cpm"
        assert result.resolution == 0.25
        assert result.k == 3
        graph = oecluster.knn_graph(items, 3)
        from_graph = oecluster.leiden(graph)
        assert from_graph.clusters == ()
        assert from_graph.k == 3


@pytest.mark.parametrize(("kwargs", "match"), [
    ({"k": 8}, "k must be between 1 and 7 for 8 items, got 8"),
    ({"k": 0}, "k must be between 1 and 7 for 8 items, got 0"),
    ({"k": 2.0}, "k must be an integer"),
    ({"k": 3, "objective": "CPM"},
     "objective must be 'modularity' or 'cpm', got 'CPM'"),
    ({"k": 3, "objective": None},
     "objective must be 'modularity' or 'cpm', got None"),
    ({"k": 3, "resolution": -0.5},
     r"resolution must be finite and non-negative, got -0\.5"),
    ({"k": 3, "resolution": math.inf},
     "resolution must be finite and non-negative"),
    ({"k": 3, "resolution": 10**400},
     "resolution must be finite and non-negative"),
    ({"k": 3, "prune": 1.0}, r"prune must be finite and in \[0, 1\), got 1\.0"),
    ({"k": 3, "prune": -0.1}, r"prune must be finite and in \[0, 1\)"),
    ({"k": 3, "prune": math.nan}, r"prune must be finite and in \[0, 1\)"),
    ({"k": 3, "theta": 0.0}, "theta must be finite and positive, got 0.0"),
    ({"k": 3, "theta": math.inf}, "theta must be finite and positive"),
    ({"k": 3, "resolution": True}, "resolution must be finite and non-negative"),
    ({"k": 3, "resolution": "0.5"}, "resolution must be finite and non-negative"),
    ({"k": 3, "resolution": None}, "resolution must be finite and non-negative"),
    ({"k": 3, "resolution": "abc"}, "resolution must be finite and non-negative"),
    ({"k": 3, "prune": True}, r"prune must be finite and in \[0, 1\)"),
    ({"k": 3, "prune": "0.5"}, r"prune must be finite and in \[0, 1\)"),
    ({"k": 3, "prune": None}, r"prune must be finite and in \[0, 1\)"),
    ({"k": 3, "prune": "abc"}, r"prune must be finite and in \[0, 1\)"),
    ({"k": 3, "theta": True}, "theta must be finite and positive"),
    ({"k": 3, "theta": "0.5"}, "theta must be finite and positive"),
    ({"k": 3, "theta": None}, "theta must be finite and positive"),
    ({"k": 3, "theta": "abc"}, "theta must be finite and positive"),
    ({"k": 3, "n_iterations": -2},
     "n_iterations must be between -1 and 9223372036854775807, got -2"),
    ({"k": 3, "n_iterations": 2**63},
     ("n_iterations must be between -1 and 9223372036854775807, got "
      "9223372036854775808")),
    ({"k": 3, "n_iterations": 1.5}, "n_iterations must be an integer"),
    ({"k": 3, "n_iterations": True}, "n_iterations must be an integer"),
    ({"k": 3, "seed": -1},
     "seed must be between 0 and 18446744073709551615, got -1"),
    ({"k": 3, "seed": 2**64},
     ("seed must be between 0 and 18446744073709551615, got "
      "18446744073709551616")),
    ({"k": 3, "seed": False}, "seed must be an integer"),
    ({"k": 3, "num_threads": -1}, "num_threads must be at least 0"),
    ({"k": 3, "chunk_size": 0}, "chunk_size must be at least 1"),
    ({"k": 3, "similarity": True}, "similarity=True is not supported"),
])
def test_invalid_options_are_value_errors(kwargs, match):
    matrix = _dense_distance_matrix(_square(_POINTS))
    with pytest.raises(ValueError, match=match):
        oecluster.leiden(matrix, **kwargs)


@pytest.mark.parametrize("kwargs", [
    {"n_iterations": -1},
    {"seed": 0},
    {"seed": 2**64 - 1},
    {"objective": "cpm", "resolution": 0.0},
    {"prune": 0.0},
])
def test_boundary_values_are_accepted(kwargs):
    result = oecluster.leiden(_dense_distance_matrix(_square(_POINTS)), k=3,
                              **kwargs)
    assert result.labels.tolist() == _BLOCKS
    assert result.iterations == 2


def test_the_largest_iteration_count_is_accepted_and_converts():
    # Running 2**63 - 1 passes never finishes, so the validator and the
    # native conversion are checked separately.
    largest = 2**63 - 1
    assert oecluster._leiden_int(largest, "n_iterations", -1, largest) == largest
    options = oecluster.oecluster.LeidenOptions()
    options.n_iterations = largest
    assert options.n_iterations == largest


def test_graph_specific_refusals():
    graph = oecluster.knn_graph(_dense_distance_matrix(_square(_POINTS)), 3)
    with pytest.raises(ValueError, match="k=2 does not match the graph's k=3"):
        oecluster.leiden(graph, k=2)
    with pytest.raises(ValueError, match="chunk_size must be at least 1"):
        oecluster.leiden(graph, chunk_size=0)
    with pytest.raises(ValueError, match="seed must be between"):
        oecluster.leiden(graph, seed=-1)
    with pytest.raises(TypeError, match="already fixes the neighbors"):
        oecluster.leiden(graph, comparison="fingerprint")
    with pytest.raises(TypeError, match="already fixes the neighbors"):
        oecluster.leiden(graph, radius=2)


def test_raw_input_requires_k():
    with pytest.raises(TypeError, match="requires k= unless items is a KNNGraph"):
        oecluster.leiden(_dense_distance_matrix(_square(_POINTS)))


def test_options_are_keyword_only():
    with pytest.raises(TypeError):
        oecluster.leiden(_dense_distance_matrix(_square(_POINTS)), 3)


def test_more_than_int_max_items_are_refused():
    oecluster._leiden_check_size(2147483647)
    with pytest.raises(ValueError,
                       match="supports at most 2147483647 items, got "
                             "2147483648"):
        oecluster._leiden_check_size(2147483648)


def test_a_single_item_is_refused():
    with pytest.raises(ValueError, match="a single item has no neighbors"):
        oecluster.leiden(_dense_distance_matrix([[0.0]]), k=1)


def test_an_empty_sequence_is_refused():
    with pytest.raises(ValueError, match="requires at least one item"):
        oecluster.leiden([], k=1, comparison="fingerprint")


def test_a_sparse_short_item_is_a_runtime_error():
    matrix = _sparse_distance_matrix(_square(_POINTS), 1.0)
    with pytest.raises(RuntimeError, match="raise the cutoff or lower k"):
        oecluster.leiden(matrix, k=3)


def test_matrices_the_gate_refuses_are_value_errors():
    similarity = oecluster.pdist(_mols(FP_SMILES), "fingerprint",
                                 similarity=True)
    with pytest.raises(ValueError, match="holds a similarity"):
        oecluster.leiden(similarity, k=2)
    stamped = _sparse_distance_matrix(_square(_POINTS), 100.0)
    stamped._facts['is_distance'] = False
    with pytest.raises(ValueError, match="holds a similarity"):
        oecluster.leiden(stamped, k=2)
    square = _square(_POINTS)
    square[1][3] = square[3][1] = math.nan
    with pytest.raises(ValueError, match="non-finite"):
        oecluster.leiden(_dense_distance_matrix(square), k=2)
    with pytest.raises(ValueError, match="non-finite"):
        oecluster.leiden(_sparse_distance_matrix(square, 100.0), k=2)


def test_a_similarity_comparison_is_refused():
    prebuilt = oecluster.FingerprintComparison(_mols(FP_SMILES),
                                               similarity=True)
    with pytest.raises(ValueError, match="requires distances"):
        oecluster.leiden(prebuilt, k=2)


@pytest.mark.parametrize("call", [
    lambda mols: oecluster.leiden(
        oecluster.pdist(mols, "fingerprint"), k=2, comparison="fingerprint"),
    lambda mols: oecluster.leiden(mols, k=2),
    lambda mols: oecluster.leiden(
        oecluster.cdist(mols[:2], mols[2:4], "fingerprint"), k=1),
])
def test_arguments_that_fit_no_path_are_type_errors(call):
    with pytest.raises(TypeError):
        call(_mols(FP_SMILES))


def test_an_empty_similarity_matrix_returns_an_empty_result():
    matrix = oecluster.pdist([], "fingerprint", similarity=True)
    result = oecluster.leiden(matrix, k=7)
    assert result.num_samples == 0
    assert result.clusters == ()
    assert result.k == 7


def test_an_empty_similarity_comparison_returns_an_empty_result():
    comparison = oecluster.FingerprintComparison([], similarity=True)
    result = oecluster.leiden(comparison, k=7)
    assert result.num_samples == 0
    assert result.clusters == ()
    assert result.k == 7


def test_the_result_type_defaults_and_stores_every_field():
    default = oecluster.LeidenResult([], [])
    assert (default.quality, default.iterations, default.objective,
            default.resolution, default.k, default.excluded) == (
        0.0, 0, "modularity", 0.0, 0, [])
    result = oecluster.LeidenResult(
        [0, -1], [[0]], quality=0.25, iterations=3, objective="cpm",
        resolution=0.5, k=4, excluded=[[1, "missing-descriptor"]])
    assert (result.quality, result.iterations, result.objective,
            result.resolution, result.k) == (0.25, 3, "cpm", 0.5, 4)
    assert result.excluded == [[1, "missing-descriptor"]]
    assert result.method == "leiden"


def test_the_leiden_surface_is_exported():
    assert "leiden" in oecluster.__all__
    assert "LeidenResult" in oecluster.__all__
    assert issubclass(oecluster.LeidenResult, oecluster.ClusteringResult)
