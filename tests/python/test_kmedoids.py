"""Python surface, validation mirror and degenerate cases for k_medoids."""

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
])
def test_negative_integer_arguments_raise_value_error(kwargs, message):
    import oecluster

    with pytest.raises(ValueError, match=message):
        oecluster.k_medoids(_two_triples(), **kwargs)


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
