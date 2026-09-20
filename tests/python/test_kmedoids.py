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

    Two pairs at unequal internal spacing: (0, 0.5) tight and (10, 11) loose.
    For k=2, BUILD picks the center of each pair → medoids (1, 2), while
    farthest-first picks the two farthest points → medoids (1, 3), so the
    strategies provably differ.
    """
    positions = np.array([0.0, 0.5, 10.0, 11.0])
    return _dense_distance_matrix(np.abs(positions[:, None] - positions[None, :]))


def _fractional_cost():
    """Four points whose optimal cost is not an integer.

    Two pairs (0, 0.25) and (10, 10.25) separated by a large gap. Every
    position is a dyadic rational, so the optimal cost is exactly 0.5 rather
    than approximately so, and int(0.5) cannot pass as the real value.
    """
    positions = np.array([0.0, 0.25, 10.0, 10.25])
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

    # End-to-end call on a comparable matrix succeeds with size_t maximum.
    result = oecluster.k_medoids(_two_triples(), n_clusters=2, chunk_size=size_t_max)

    assert result.medoids == (1, 4)
    assert result.cost == 1.0
    assert result.converged is True


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

    assert list(params.keys()) == [
        "distance_matrix", "n_clusters", "init", "initial_medoids",
        "max_iterations", "num_threads", "chunk_size"
    ]

    # distance_matrix is positional-or-keyword; all others are keyword-only.
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

    def capture_and_return(_storage, options):
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


def _mols():
    """Six molecules, copied from tests/python/test_metric_gate.py:30.

    Fingerprints need no coordinates, so these come straight from SMILES as
    OEGraphMol rather than through Omega.
    """
    from openeye import oechem

    smiles = ["CCO", "CCC", "CCCC", "c1ccccc1", "CCN", "CCOC"]
    mols = []
    for idx, smi in enumerate(smiles):
        mol = oechem.OEGraphMol()
        oechem.OESmilesToMol(mol, smi)
        mol.SetTitle(f"mol{idx}")
        mols.append(mol)
    return mols


# At a converged PAM solution every medoid is also the within-cluster medoid of
# its own cluster: otherwise swapping in the better member would lower the total
# and the verification pass would have taken that swap. This cross-checks two
# independently written implementations.
def test_every_medoid_is_its_own_cluster_s_representative():
    import oecluster

    dm = _two_triples()
    result = oecluster.k_medoids(dm, n_clusters=2)
    assert result.converged is True

    for label, cluster in enumerate(result.clusters):
        chosen = oecluster.representative(list(cluster), dm, method="medoid")
        assert chosen == result.medoids[label]


# The global optimum here is a tie: items 2 and 3 both sum to 30.0 across the
# whole set. The two implementations land on the same item by different rules --
# k_medoids breaks ties toward the smaller item index, while representative
# stable-sorts on score and so returns whichever tied item came first in the
# caller's list -- so they coincide only because the list passed here ascends.
# Each index is therefore asserted outright rather than merely against the
# other: a bare equality survives a fixture drifting out of the tie, since both
# implementations would then agree on whatever the new sole winner was.
def test_single_cluster_agrees_with_the_global_representative():
    import oecluster

    dm = _two_triples()
    result = oecluster.k_medoids(dm, n_clusters=1)

    assert result.medoids[0] == 2
    assert oecluster.representative(
        list(range(dm.num_samples)), dm, method="medoid") == 2
    # Reversing the input proves the rules really are distinct, so the comment
    # above cannot quietly rot into describing a shared index rule.
    assert oecluster.representative(
        list(reversed(range(dm.num_samples))), dm, method="medoid") == 3


# Dice is the repo's canonical non-metric fingerprint distance. This is the test
# that pins the require_comparable decision: if someone later "fixes" the gate to
# require_metric for consistency with the other clustering entry points, this is
# what tells them they changed the contract.
def test_a_non_metric_matrix_clusters_with_no_flag_at_all():
    import oecluster

    dm = oecluster.pdist(_mols(), "fingerprint", metric="dice")
    assert dm.metric_capabilities == {'zero_self': True, 'triangle': False}

    result = oecluster.k_medoids(dm, n_clusters=2)

    assert len(result.labels) == 6
    assert result.num_clusters == 2


# The parameter does not exist, and that is deliberate: there is no triangle
# assumption for it to override.
def test_allow_nonmetric_is_not_a_parameter():
    import oecluster

    dm = _two_triples()
    with pytest.raises(TypeError, match="allow_nonmetric"):
        oecluster.k_medoids(
            dm, n_clusters=2,
            allow_nonmetric=True)  # pyright: ignore[reportCallIssue]


def test_a_similarity_matrix_is_refused():
    import oecluster

    dm = oecluster.pdist(_mols(), "fingerprint", similarity=True)
    with pytest.raises(ValueError, match="similarity=False"):
        oecluster.k_medoids(dm, n_clusters=2)


# A tier-1 refusal require_comparable keeps and never waives, poked the way
# tests/python/test_metric_gate.py:24 pokes it.
def test_a_matrix_holding_non_finite_entries_is_refused():
    import math

    import oecluster

    dm = oecluster.pdist(_mols(), "fingerprint")
    dm.condensed[0] = math.nan

    with pytest.raises(ValueError, match="non-finite entries"):
        oecluster.k_medoids(dm, n_clusters=2)


# The other tier-1 refusal, and the one PAM depends on most directly: BUILD
# and every swap evaluation assume an item's distance to itself is the smallest
# it can be, so a nonzero diagonal would let a medoid lose its own cluster.
# Stamped directly, following tests/python/test_metric_gate.py:22.
def test_a_non_zero_self_distance_matrix_is_refused():
    import oecluster

    dm = oecluster.pdist(_mols(), "fingerprint")
    dm._facts['zero_self'] = False

    with pytest.raises(ValueError, match="zero self-distance"):
        oecluster.k_medoids(dm, n_clusters=2)


# The one refusal require_comparable keeps beyond the tier-1 checks, and the
# one that matters most here: PAM adds incomparable distances together rather
# than merely ranking them. Stamped directly, following
# tests/python/test_metric_gate.py:1413.
def test_a_subset_scored_matrix_is_refused():
    import oecluster

    dm = oecluster.pdist(_mols(), "fingerprint", metric="dice")
    dm._facts['data_integrity'] = "subset_scored"

    with pytest.raises(ValueError, match="subset"):
        oecluster.k_medoids(dm, n_clusters=2)


# The asymmetry from the design: k_medoids accepts a Dice matrix that
# cluster_report refuses, because the internal validity indices do lean on
# metric behavior where PAM does not. Tested so it stays documented behavior
# rather than an accident.
def test_cluster_report_still_requires_the_stronger_gate():
    import oecluster

    dm = oecluster.pdist(_mols(), "fingerprint", metric="dice")
    result = oecluster.k_medoids(dm, n_clusters=2)

    with pytest.raises(ValueError, match="triangle inequality"):
        oecluster.cluster_report(result, dm)
    oecluster.cluster_report(result, dm, allow_nonmetric=True)


def test_cluster_report_accepts_a_metric_k_medoids_result():
    import oecluster

    dm = oecluster.pdist(_mols(), "fingerprint")
    result = oecluster.k_medoids(dm, n_clusters=2)

    report = oecluster.cluster_report(result, dm)

    assert report.num_clusters == 2
