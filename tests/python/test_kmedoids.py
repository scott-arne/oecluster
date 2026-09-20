"""Python surface, validation mirror and degenerate cases for k_medoids."""

import ctypes
import inspect

import numpy as np
import pytest


def _dense_distance_matrix(square):
    """Build a dense DistanceMatrix from a square distance matrix."""
    from oecluster import DenseStorage, SymmetricDistanceMatrix

    square = np.asarray(square, dtype=np.float64)
    storage = DenseStorage(square.shape[0])
    for i in range(square.shape[0]):
        for j in range(i + 1, square.shape[0]):
            storage.Set(i, j, float(square[i, j]))
    return SymmetricDistanceMatrix(
        storage, "test", [f"item_{i}" for i in range(square.shape[0])], {})


def _two_triples():
    """Two tight triples of three items separated by a wide gap.

    Every position is a dyadic rational, so the optimal cost is exactly 1.0
    rather than approximately so.
    """
    positions = np.array([0.0, 0.25, 0.5, 10.0, 10.25, 10.5])
    return _dense_distance_matrix(np.abs(positions[:, None] - positions[None, :]))


def _asymmetric_pair():
    """Four points where BUILD and farthest-first initialization disagree.

    Two tight pairs (0, 0.5) and (10, 11) separated by a large gap. For
    k=2, BUILD picks the center of each pair → medoids (1, 2), while
    farthest-first picks the two farthest points → medoids (1, 3). Both are
    distance-symmetric, so the fixture is deliberately asymmetric to ensure
    the two strategies provably differ.
    """
    positions = np.array([0.0, 0.5, 10.0, 11.0])
    return _dense_distance_matrix(np.abs(positions[:, None] - positions[None, :]))


def _fractional_cost():
    """Four points whose optimal cost is not an integer.

    Two pairs (0, 0.3) and (10, 10.7) with fractional separations, ensuring
    the cost cannot pass through int(self._cost) unchanged.
    """
    positions = np.array([0.0, 0.3, 10.0, 10.7])
    return _dense_distance_matrix(np.abs(positions[:, None] - positions[None, :]))


def test_k_medoids_finds_the_medoid_of_each_tight_triple():
    import oecluster

    result = oecluster.k_medoids(_two_triples(), n_clusters=2)

    assert result.medoids == (1, 4)
    assert result.labels.tolist() == [0, 0, 0, 1, 1, 1]
    assert result.clusters == ((0, 1, 2), (3, 4, 5))
    assert result.cost == 1.0
    assert result.converged is True
    assert result.method == "k_medoids"


def test_the_result_is_a_clustering_result_subclass():
    import oecluster

    result = oecluster.k_medoids(_two_triples(), n_clusters=2)

    assert isinstance(result, oecluster.KMedoidsResult)
    assert isinstance(result, oecluster.ClusteringResult)
    assert len(result) == 2
    assert result[0] == (0, 1, 2)
    assert result.num_samples == 6
    assert repr(result) == "KMedoidsResult(num_clusters=2, num_samples=6)"
    # k-medoids carries no centroids: a result never holds a field it cannot fill.
    assert not hasattr(result, "centroids")


def test_result_scalars_have_exact_types():
    """Scalar properties must be exact types, not bool/int subclasses."""
    import oecluster

    result = oecluster.k_medoids(_fractional_cost(), n_clusters=2)

    # type(x) is T rejects subclasses; isinstance accepts them.
    assert type(result.cost) is float
    assert type(result.n_iterations) is int
    assert type(result.converged) is bool


def test_cost_is_recomputed_and_matches_distances():
    """Cost must be the sum of item-to-medoid distances, not a fabrication."""
    import oecluster

    dm = _fractional_cost()
    result = oecluster.k_medoids(dm, n_clusters=2)

    # Recompute the cost independently: sum of each item's distance to its medoid.
    expected_cost = 0.0
    for i, label in enumerate(result.labels):
        medoid = result.medoids[label]
        if i != medoid:
            expected_cost += dm.storage.Get(min(i, medoid), max(i, medoid))

    assert result.cost == expected_cost
    # The cost must not be an integer, so int(self._cost) cannot pass.
    assert result.cost != int(result.cost)


# Case, not punctuation: the mirror lowercases the name and looks it up, so
# "FARTHEST_FIRST" matches and a CamelCase "FarthestFirst" would not -- its
# lowercase form drops the underscore. The documented spellings are the two
# snake_case ones; no CamelCase alias is offered.
@pytest.mark.parametrize("init",
                         ["build", "BUILD", "farthest_first", "FARTHEST_FIRST"])
def test_init_accepts_both_spellings_case_insensitively(init):
    import oecluster

    result = oecluster.k_medoids(_two_triples(), n_clusters=2, init=init)

    assert result.medoids == (1, 4)


def test_explicit_initialization_reaches_the_same_optimum():
    import oecluster

    result = oecluster.k_medoids(
        _two_triples(), n_clusters=2, init="explicit", initial_medoids=[0, 1])

    assert result.medoids == (1, 4)
    assert result.converged is True
    assert result.n_iterations > 0


def test_build_and_farthest_first_initialization_strategies_differ():
    """BUILD and farthest-first must route to different enumerators.

    On a symmetric fixture they may accidentally agree; this asymmetric one
    ensures they provably diverge.
    """
    import oecluster

    result_build = oecluster.k_medoids(_asymmetric_pair(), n_clusters=2, init="build")
    result_ff = oecluster.k_medoids(_asymmetric_pair(), n_clusters=2,
                                    init="farthest_first")

    assert result_build.medoids == (1, 2)
    assert result_ff.medoids == (1, 3)
    assert result_build.medoids != result_ff.medoids


def test_every_item_receives_a_label():
    import oecluster

    result = oecluster.k_medoids(_two_triples(), n_clusters=3)

    assert -1 not in result.labels.tolist()
    assert result.num_clusters == 3
    assert all(len(cluster) > 0 for cluster in result.clusters)


def test_k_equals_n_returns_the_identity_partition():
    import oecluster

    result = oecluster.k_medoids(_two_triples(), n_clusters=6)

    assert result.medoids == (0, 1, 2, 3, 4, 5)
    assert result.cost == 0.0
    assert result.n_iterations == 0
    assert result.converged is True


def test_a_capped_run_reports_unconverged():
    import oecluster

    result = oecluster.k_medoids(
        _two_triples(), n_clusters=2, init="explicit",
        initial_medoids=[0, 1], max_iterations=1)

    assert result.converged is False
    assert result.n_iterations == 1


def test_the_matrix_argument_must_be_a_distance_matrix():
    import oecluster

    with pytest.raises(TypeError, match="SymmetricDistanceMatrix"):
        oecluster.k_medoids([[0.0, 1.0], [1.0, 0.0]], n_clusters=1)


def test_an_unknown_init_names_the_three_spellings():
    import oecluster

    with pytest.raises(ValueError, match="farthest_first"):
        oecluster.k_medoids(_two_triples(), n_clusters=2, init="kmeans++")


# The unknown-init row is LAST in the native table, so every other invalid
# option outranks it. These pin the order rather than merely the messages: a
# mirror that rejects a bad init eagerly passes the test above and fails here,
# and would report a different first failure than a direct C++ caller sees.
@pytest.mark.parametrize("kwargs,message", [
    ({"n_clusters": 2, "chunk_size": 0}, "chunk_size"),
    ({"n_clusters": 2, "max_iterations": 0}, "max_iterations"),
    ({"n_clusters": 0}, "at least one"),
    ({"n_clusters": 7}, "at most the item count"),
    ({"n_clusters": 2, "initial_medoids": [0, 3]},
     "requires an explicit initialization"),
])
def test_every_other_invalid_option_outranks_an_unknown_init(kwargs, message):
    import oecluster

    with pytest.raises(ValueError, match=message):
        oecluster.k_medoids(_two_triples(), init="kmeans++", **kwargs)


@pytest.mark.parametrize("kwargs", [
    {"n_clusters": 2.5},
    {"max_iterations": 1.5},
    {"num_threads": "two"},
    {"chunk_size": None},
])
def test_non_integer_arguments_raise_type_error(kwargs):
    import oecluster

    with pytest.raises(TypeError):
        oecluster.k_medoids(_two_triples(), **kwargs)


def test_non_integer_seeds_raise_type_error():
    import oecluster

    with pytest.raises(TypeError, match="initial_medoids"):
        oecluster.k_medoids(
            _two_triples(), n_clusters=2, init="explicit",
            initial_medoids=[0.0, 1.0])


@pytest.mark.parametrize("kwargs,message", [
    ({"n_clusters": -1}, "n_clusters must be non-negative"),
    ({"max_iterations": -1}, "max_iterations must be non-negative"),
    ({"num_threads": -1}, "num_threads must be non-negative"),
    ({"chunk_size": -1}, "chunk_size must be non-negative"),
    ({"n_clusters": 2, "init": "explicit", "initial_medoids": [-1, 0]},
     "initial_medoids must be non-negative"),
])
def test_negative_integer_arguments_raise_value_error(kwargs, message):
    import oecluster

    with pytest.raises(ValueError, match=message):
        oecluster.k_medoids(_two_triples(), **kwargs)


@pytest.mark.parametrize("kwargs,message", [
    ({"max_iterations": 1 << 100}, "max_iterations exceeds size_t maximum"),
    ({"num_threads": 1 << 100}, "num_threads exceeds size_t maximum"),
    ({"chunk_size": 1 << 100}, "chunk_size exceeds size_t maximum"),
])
def test_oversized_integer_arguments_raise_value_error(kwargs, message):
    import oecluster

    with pytest.raises(ValueError, match=message):
        oecluster.k_medoids(_two_triples(), **kwargs)


def test_oversized_n_clusters_still_reports_item_count_bound():
    import oecluster

    with pytest.raises(ValueError, match="at most the item count"):
        oecluster.k_medoids(_two_triples(), n_clusters=1 << 100)


def test_oversized_check_fires_before_gate():
    """An oversized integer is refused before the gate runs.

    If the gate ran first, an oversized chunk_size on a non-comparable matrix
    would report the gate's refusal instead of the field's bound.
    """
    import oecluster
    from openeye import oechem

    # Create a subset_scored matrix that the gate would refuse.
    mols = []
    for smi in ["C", "CC", "CCC", "CCCC"]:
        mol = oechem.OEGraphMol()
        oechem.OESmilesToMol(mol, smi)
        mols.append(mol)
    dm = oecluster.pdist(mols, "descriptor", metric="euclidean", missing="ignore")
    assert dm.data_integrity == "subset_scored"

    with pytest.raises(ValueError, match="chunk_size exceeds size_t maximum"):
        oecluster.k_medoids(dm, n_clusters=2, chunk_size=1 << 100)


def test_size_t_maximum_is_accepted():
    """The exact size_t maximum is accepted and reaches the setter."""
    import oecluster

    size_t_max = (1 << (8 * ctypes.sizeof(ctypes.c_size_t))) - 1

    # Reaching the setter without OverflowError is the test.
    # The call will fail at the gate (subset_scored) or elsewhere, but not
    # on the bound check or the SWIG setter.
    from openeye import oechem
    mols = []
    for smi in ["C", "CC", "CCC", "CCCC"]:
        mol = oechem.OEGraphMol()
        oechem.OESmilesToMol(mol, smi)
        mols.append(mol)
    dm = oecluster.pdist(mols, "descriptor", metric="euclidean", missing="ignore")

    with pytest.raises(ValueError, match="subset"):
        oecluster.k_medoids(dm, n_clusters=2, chunk_size=size_t_max)


def test_sparse_storage_is_refused():
    import oecluster

    storage = oecluster.SparseStorage(4, 0.5)
    dm = oecluster.SymmetricDistanceMatrix(storage, "test", list("abcd"), {})

    with pytest.raises(ValueError, match="SparseStorage is not supported"):
        oecluster.k_medoids(dm, n_clusters=2)


def test_a_zero_chunk_size_is_refused():
    import oecluster

    with pytest.raises(ValueError, match="chunk_size must be at least one"):
        oecluster.k_medoids(_two_triples(), n_clusters=2, chunk_size=0)


def test_zero_max_iterations_is_refused():
    import oecluster

    with pytest.raises(ValueError, match="max_iterations must be at least one"):
        oecluster.k_medoids(_two_triples(), n_clusters=2, max_iterations=0)


def test_zero_clusters_is_refused():
    import oecluster

    with pytest.raises(ValueError, match="n_clusters must be at least one"):
        oecluster.k_medoids(_two_triples(), n_clusters=0)


def test_more_clusters_than_items_is_refused():
    import oecluster

    with pytest.raises(ValueError, match="n_clusters must be at most"):
        oecluster.k_medoids(_two_triples(), n_clusters=7)


def test_seeds_without_explicit_initialization_are_refused():
    import oecluster

    with pytest.raises(ValueError, match="requires an explicit initialization"):
        oecluster.k_medoids(_two_triples(), n_clusters=2, initial_medoids=[0, 3])


def test_explicit_initialization_without_seeds_is_refused():
    import oecluster

    with pytest.raises(ValueError, match="exactly n_clusters indices"):
        oecluster.k_medoids(_two_triples(), n_clusters=2, init="explicit")


def test_the_wrong_number_of_seeds_is_refused():
    import oecluster

    with pytest.raises(ValueError, match="exactly n_clusters indices"):
        oecluster.k_medoids(
            _two_triples(), n_clusters=2, init="explicit", initial_medoids=[0])


# IndexError, not RuntimeError: the native layer raises std::out_of_range, which
# the GIL wrapper cannot preserve, so the mirror is what keeps the type right.
def test_an_out_of_range_seed_raises_index_error():
    import oecluster

    with pytest.raises(IndexError, match="outside the storage range"):
        oecluster.k_medoids(
            _two_triples(), n_clusters=2, init="explicit", initial_medoids=[0, 6])


def test_duplicate_seeds_are_refused():
    import oecluster

    with pytest.raises(ValueError, match="must be unique"):
        oecluster.k_medoids(
            _two_triples(), n_clusters=2, init="explicit", initial_medoids=[3, 3])


# Multiple simultaneous failures report the first one in native row order.
# A reordering that leaves single-failure tests green surfaces here. The two
# adjacent-pair cases below complete the set: every other adjacent pair is
# mutually exclusive and cannot be pinned by any input.
@pytest.mark.parametrize("kwargs,exc_type,message", [
    ({"n_clusters": 0, "max_iterations": 0, "chunk_size": 0},
     ValueError, "chunk_size must be at least one"),
    ({"n_clusters": 0, "max_iterations": 0},
     ValueError, "max_iterations must be at least one"),
    ({"n_clusters": 7, "chunk_size": 0},
     ValueError, "chunk_size must be at least one"),
    ({"n_clusters": 2, "init": "explicit", "initial_medoids": [6, 6]},
     IndexError, "outside the storage range"),
    ({"n_clusters": 7, "initial_medoids": [0]},
     ValueError, "at most the item count"),
    ({"n_clusters": 2, "init": "explicit", "initial_medoids": [0, 6, 1]},
     ValueError, "exactly n_clusters indices"),
])
def test_the_first_reported_failure_matches_the_native_row_order(kwargs, exc_type, message):
    import oecluster

    with pytest.raises(exc_type, match=message):
        oecluster.k_medoids(_two_triples(), **kwargs)


def test_sparse_storage_outranks_scalar_checks():
    import oecluster

    storage = oecluster.SparseStorage(4, 0.5)
    dm = oecluster.SymmetricDistanceMatrix(storage, "test", list("abcd"), {})

    with pytest.raises(ValueError, match="SparseStorage is not supported"):
        oecluster.k_medoids(dm, n_clusters=0, chunk_size=0)


# A missed mirror surfaces as RuntimeError. Nothing here may raise one.
@pytest.mark.parametrize("kwargs", [
    {"n_clusters": 0},
    {"n_clusters": 7},
    {"chunk_size": 0},
    {"max_iterations": 0},
    {"n_clusters": 2, "initial_medoids": [0, 3]},
    {"n_clusters": 2, "init": "explicit"},
    {"n_clusters": 2, "init": "explicit", "initial_medoids": [0, 6]},
    {"n_clusters": 2, "init": "explicit", "initial_medoids": [3, 3]},
])
def test_no_validation_failure_surfaces_as_runtime_error(kwargs):
    import oecluster

    with pytest.raises((ValueError, TypeError, IndexError)):
        oecluster.k_medoids(_two_triples(), **kwargs)


def test_the_surface_is_exported():
    import oecluster

    for name in ("k_medoids", "KMedoidsResult", "KMedoidsOptions"):
        assert hasattr(oecluster, name)
        assert name in oecluster.__all__


def test_the_signature_is_pinned():
    """Pin every parameter name, default, and the keyword-only shape."""
    import oecluster

    sig = inspect.signature(oecluster.k_medoids)
    params = sig.parameters

    # distance_matrix: positional-only would be an improvement but isn't
    # the current shape; for now it's positional-or-keyword, first position.
    assert list(params.keys()) == [
        "distance_matrix", "n_clusters", "init", "initial_medoids",
        "max_iterations", "num_threads", "chunk_size"
    ]

    # All except distance_matrix are keyword-only.
    assert params["distance_matrix"].kind == inspect.Parameter.POSITIONAL_OR_KEYWORD
    for name in ["n_clusters", "init", "initial_medoids", "max_iterations",
                 "num_threads", "chunk_size"]:
        assert params[name].kind == inspect.Parameter.KEYWORD_ONLY

    # Defaults
    assert params["n_clusters"].default == 2
    assert params["init"].default == "build"
    assert params["initial_medoids"].default is None
    assert params["max_iterations"].default == 100
    assert params["num_threads"].default == 0
    assert params["chunk_size"].default == 4096


def test_default_arguments_produce_expected_behavior():
    """Passing no keywords must use BUILD initialization and n_clusters=2."""
    import oecluster

    result = oecluster.k_medoids(_asymmetric_pair())

    # Default init="build" → BUILD strategy
    assert result.medoids == (1, 2)
    # Default n_clusters=2
    assert result.num_clusters == 2


def test_positive_num_threads_and_chunk_size_are_forwarded(monkeypatch):
    """Non-default positive values must reach the options unchanged."""
    import oecluster

    captured_options = []

    class MockResult:
        """Stub result that satisfies KMedoidsResult's unpacking."""
        def Labels(self):
            return np.array([0, 0])
        def Members(self):
            return ((0, 1),)
        def Medoids(self):
            return (0,)
        def Cost(self):
            return 0.0
        def NumIterations(self):
            return 0
        def Converged(self):
            return True

    def capture_and_return(storage, options):
        captured_options.append(options)
        return MockResult()

    monkeypatch.setattr(oecluster, "_k_medoids_cluster", capture_and_return)

    # Call with non-default positive values.
    dm = _two_triples()
    oecluster.k_medoids(dm, n_clusters=2, num_threads=4, chunk_size=1024)

    assert len(captured_options) == 1
    opts = captured_options[0]
    assert opts.num_threads == 4
    assert opts.chunk_size == 1024
