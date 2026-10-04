"""Python surface of consensus clustering."""
import math

import numpy as np
import oecluster
import pytest

# --- fixtures ----------------------------------------------------------------


def _blobs20():
    """Twenty points, two blobs of ten: intra-blob 0.1, inter-blob 0.9."""
    values = [0.1 if (i < 10) == (j < 10) else 0.9
              for i in range(20) for j in range(i + 1, 20)]
    return oecluster.SymmetricDistanceMatrix.from_condensed(np.array(values))


def _result(labels):
    """A ClusteringResult from labels alone; -1 is noise."""
    labels = np.asarray(labels)
    clusters = [tuple(int(p) for p in np.flatnonzero(labels == label))
                for label in sorted(set(labels.tolist())) if label >= 0]
    return oecluster.ClusteringResult(labels, clusters)


def _mols(smiles_list):
    from openeye import oechem

    mols = []
    for index, smiles in enumerate(smiles_list):
        mol = oechem.OEGraphMol()
        oechem.OESmilesToMol(mol, smiles)
        mol.SetTitle(f"mol{index}")
        mols.append(mol)
    return mols


def _members(ensemble, num_items):
    """Normalize an ensemble to ``(positions, labels)`` lists of ints."""
    normalized = []
    for element in ensemble:
        if isinstance(element, oecluster.ClusteringResult):
            normalized.append((list(range(num_items)),
                               [int(v) for v in element.labels]))
        else:
            positions, labels = element
            normalized.append(([int(p) for p in positions],
                               [int(v) for v in labels]))
    return normalized


def _slow_matrix(ensemble, num_items):
    """Co-association distances with dictionaries; shares nothing with the module."""
    members = _members(ensemble, num_items)
    distances = {}
    for i in range(num_items):
        for j in range(i + 1, num_items):
            co = 0
            obs = 0
            for positions, labels in members:
                if i not in positions or j not in positions:
                    continue
                obs += 1
                left = labels[positions.index(i)]
                right = labels[positions.index(j)]
                if left >= 0 and left == right:
                    co += 1
            distances[(i, j)] = 1.0 if obs == 0 else 1.0 - co / obs
    return distances


def _slow_components(distances, num_items, threshold):
    """Connected components of the pairs at or above `threshold`."""
    parent = list(range(num_items))

    def find(item):
        while parent[item] != item:
            parent[item] = parent[parent[item]]
            item = parent[item]
        return item

    for (i, j), distance in distances.items():
        if distance <= 1.0 - threshold:
            left, right = find(i), find(j)
            if left != right:
                parent[right] = left
    labels, seen = [0] * num_items, {}
    for item in range(num_items):
        root = find(item)
        if root not in seen:
            seen[root] = len(seen)
        labels[item] = seen[root]
    return labels


def _condensed(distances, num_items):
    return np.array([distances[(i, j)]
                     for i in range(num_items)
                     for j in range(i + 1, num_items)])


# --- native bindings ---------------------------------------------------------

def test_native_coassociation_distances_fills_a_destination():
    native = oecluster.oecluster
    destination = oecluster.DenseStorage(4)
    summary = native.coassociation_distances(
        4,
        native.SizeTVector([0, 4, 8]),
        native.SizeTVector([0, 1, 2, 3, 0, 1, 2, 3]),
        native.IntVector([0, 0, 1, 1, 0, 0, 0, 1]),
        destination)
    assert isinstance(summary, native.ConsensusMatrixSummary)
    assert summary.num_partitions == 2
    assert summary.unobserved_pairs == 0
    assert destination.Get(0, 1) == 0.0
    assert destination.Get(0, 2) == 0.5
    assert destination.Get(0, 3) == 1.0


def test_native_consensus_components_and_strength_round_trip():
    native = oecluster.oecluster
    destination = oecluster.DenseStorage(4)
    native.coassociation_distances(
        4,
        native.SizeTVector([0, 4]),
        native.SizeTVector([0, 1, 2, 3]),
        native.IntVector([0, 0, 1, 1]),
        destination)
    components = native.consensus_components(destination, 0.5)
    assert isinstance(components, native.IntVector)
    labels = list(components)
    assert labels == [0, 0, 1, 1]
    strength = native.consensus_strength(
        destination, native.IntVector(labels))
    assert isinstance(strength.item_consensus, native.DoubleVector)
    assert isinstance(strength.cluster_consensus, native.DoubleVector)
    assert list(strength.item_consensus) == [1.0, 1.0, 1.0, 1.0]
    assert list(strength.cluster_consensus) == [1.0, 1.0]


def test_native_refusal_is_a_runtime_error():
    native = oecluster.oecluster
    with pytest.raises(RuntimeError, match="more than once"):
        native.coassociation_distances(
            4,
            native.SizeTVector([0, 2]),
            native.SizeTVector([1, 1]),
            native.IntVector([0, 0]),
            oecluster.DenseStorage(4))


def test_native_strength_members_survive_a_temporary():
    native = oecluster.oecluster
    destination = oecluster.DenseStorage(4)
    native.coassociation_distances(
        4,
        native.SizeTVector([0, 4]),
        native.SizeTVector([0, 1, 2, 3]),
        native.IntVector([0, 0, 1, 1]),
        destination)
    labels = native.IntVector([0, 0, 1, 1])
    item = native.consensus_strength(destination, labels).item_consensus
    cluster = native.consensus_strength(destination, labels).cluster_consensus
    assert isinstance(item, native.DoubleVector)
    assert isinstance(cluster, native.DoubleVector)
    assert list(item) == [1.0, 1.0, 1.0, 1.0]
    assert list(cluster) == [1.0, 1.0]


# --- the consensus matrix ----------------------------------------------------

def test_matrix_matches_the_slow_reference():
    # Three members over four items; the third observes only 0, 1 and 3, and
    # item 2 is noise in the second.
    ensemble = [
        _result([0, 0, 1, 1]),
        oecluster.ClusteringResult([0, 0, -1, 1], ((0, 1), (3,))),
        (np.array([0, 1, 3]), np.array([0, 0, 1])),
    ]
    result = oecluster.consensus(ensemble, num_items=4)
    expected = _slow_matrix(ensemble, 4)
    np.testing.assert_allclose(result.matrix.condensed, _condensed(expected, 4))
    assert result.num_partitions == 3
    assert result.unobserved_pairs == 0


def test_a_pair_no_member_observed_together_is_distance_one_and_counted():
    ensemble = [(np.array([0, 1]), np.array([0, 0])),
                (np.array([2, 3]), np.array([0, 0]))]
    result = oecluster.consensus(ensemble, num_items=4)
    assert result.unobserved_pairs == 4
    assert result.matrix.squareform()[0][2] == 1.0
    assert result.matrix.squareform()[0][1] == 0.0


def test_noise_is_observed_but_never_co_clustered():
    result = oecluster.consensus([_result([0, 0, -1, -1])] * 3, num_items=4)
    square = result.matrix.squareform()
    assert square[0][1] == 0.0
    assert square[2][3] == 1.0
    assert result.unobserved_pairs == 0


def test_the_matrix_carries_consensus_facts_and_no_item_labels():
    result = oecluster.consensus([_result([0, 0, 1, 1])] * 2, num_items=4)
    matrix = result.matrix
    assert matrix.comparison_name == "consensus"
    assert matrix.params == {"num_partitions": 2}
    assert list(matrix.labels) == []
    facts = matrix.facts
    assert facts["is_distance"] is True
    assert facts["zero_self"] is True
    assert facts["triangle"] == "unknown"
    assert facts["data_integrity"] == "complete"
    assert facts["metric_probe"] in ("violations_found", "no_violations_found")


# --- extraction --------------------------------------------------------------

def test_the_default_extraction_matches_the_slow_components():
    ensemble = [
        _result([0, 0, 1, 1]),
        _result([0, 0, 0, 1]),
        (np.array([0, 1, 3]), np.array([0, 0, 1])),
    ]
    result = oecluster.consensus(ensemble, num_items=4)
    expected = _slow_components(_slow_matrix(ensemble, 4), 4, 0.5)
    assert result.labels.tolist() == expected
    assert result.threshold == 0.5
    assert result.spec is None


def test_threshold_zero_merges_everything_and_one_merges_only_unanimous_pairs():
    ensemble = [_result([0, 0, 1, 1]), _result([0, 0, 0, 1])]
    assert oecluster.consensus(ensemble, num_items=4,
                               threshold=0).num_clusters == 1
    unanimous = oecluster.consensus(ensemble, num_items=4, threshold=1)
    assert unanimous.labels.tolist() == [0, 0, 1, 2]


def test_a_method_spec_extracts_the_partition_instead():
    blobs = oecluster.butina(_blobs20(), threshold=0.2)
    spec = oecluster.ClusteringSpec("agglomerative", n_clusters=2,
                                    linkage="average")
    result = oecluster.consensus([blobs, blobs], method=spec)
    assert result.num_clusters == 2
    assert result.threshold is None
    assert result.spec == spec


def test_a_method_result_is_canonicalized():
    blobs = oecluster.butina(_blobs20(), threshold=0.2)

    def sparse_labels(items, **options):
        size = items.num_samples
        half = size // 2
        return oecluster.ClusteringResult(
            [2 ** 40] * half + [2 ** 40 + 7] * (size - half),
            (tuple(range(half)), tuple(range(half, size))))

    result = oecluster.consensus([blobs, blobs], method=sparse_labels)
    assert sorted(set(result.labels.tolist())) == [0, 1]
    assert result.clusters == (tuple(range(10)), tuple(range(10, 20)))


@pytest.mark.parametrize("name", ["bitbirch", "bitbirch_recluster",
                                  "bitbirch_refine", "murcko"])
def test_a_non_matrix_method_is_refused_by_name(name):
    with pytest.raises(ValueError, match="cannot read the consensus"):
        oecluster.consensus([_result([0, 0, 1, 1])] * 2, num_items=4,
                            method=name)


def test_threshold_and_method_are_mutually_exclusive():
    spec = oecluster.ClusteringSpec("agglomerative", n_clusters=2)
    ensemble = [_result([0, 0, 1, 1])] * 2
    with pytest.raises(ValueError, match="mutually exclusive"):
        oecluster.consensus(ensemble, num_items=4, threshold=0.5, method=spec)
    # Each alone is accepted: threshold=None is the unset sentinel.
    assert oecluster.consensus(ensemble, num_items=4, method=spec).spec == spec
    assert oecluster.consensus(ensemble, num_items=4,
                               threshold=0.5).threshold == 0.5


# --- ensemble sources --------------------------------------------------------

def test_a_cluster_stability_contributes_its_resamples():
    matrix = _blobs20()
    spec = oecluster.ClusteringSpec("butina", threshold=0.2)
    stability = oecluster.cluster_stability(spec, matrix, resamples=8, seed=0)
    result = oecluster.consensus(stability)
    assert result.num_partitions == 8
    assert result.num_clusters == 2


def test_a_cluster_stability_without_partitions_is_refused():
    matrix = _blobs20()
    spec = oecluster.ClusteringSpec("butina", threshold=0.2)
    stability = oecluster.cluster_stability(spec, matrix, resamples=4, seed=0,
                                            keep_partitions=False)
    with pytest.raises(ValueError, match="keep_partitions=True"):
        oecluster.consensus(stability)


def test_a_parameter_selection_contributes_every_row():
    selection = oecluster.select_parameter(
        "butina", _blobs20(), "threshold", [0.15, 0.2, 0.25])
    result = oecluster.consensus(selection)
    assert result.num_partitions == 3


def test_positions_may_arrive_in_any_order():
    ascending = oecluster.consensus(
        [(np.array([0, 1, 2, 3]), np.array([0, 0, 1, 1]))], num_items=4)
    descending = oecluster.consensus(
        [(np.array([3, 2, 1, 0]), np.array([1, 1, 0, 0]))], num_items=4)
    np.testing.assert_array_equal(ascending.matrix.condensed,
                                  descending.matrix.condensed)


# --- validation --------------------------------------------------------------

@pytest.mark.parametrize("kwargs, error, match", [
    ({"ensemble": "not an ensemble"}, TypeError, "ClusterStability"),
    ({"ensemble": []}, ValueError, "empty"),
    # num_items=None is the override, not an omission: the base arguments
    # below supply 4, which would make this call valid.
    ({"ensemble": [(np.array([0, 1]), np.array([0, 0]))], "num_items": None},
     ValueError, "num_items is required"),
    ({"num_items": 19}, ValueError, "19"),
    ({"num_items": True}, TypeError, "num_items"),
    ({"num_items": 1}, ValueError, "at least 2"),
    ({"threshold": 1.5}, ValueError, "between 0 and 1"),
    ({"threshold": math.nan}, ValueError, "between 0 and 1"),
    ({"threshold": True}, TypeError, "real number"),
    ({"threshold": "half"}, TypeError, "real number"),
    ({"noise": "cluster"}, ValueError, "noise"),
    ({"num_threads": 1.5}, TypeError, "num_threads"),
    ({"num_threads": -1}, ValueError, "non-negative"),
    ({"output": 7}, TypeError, "output"),
])
def test_validation_refuses_bad_arguments(kwargs, error, match):
    arguments = {"ensemble": [_result([0, 0, 1, 1])] * 2, "num_items": 4}
    arguments.update(kwargs)
    with pytest.raises(error, match=match):
        oecluster.consensus(**arguments)


def test_a_member_of_the_wrong_length_or_shape_is_refused():
    with pytest.raises(ValueError, match="one-dimensional"):
        oecluster.consensus(
            [oecluster.ClusteringResult([[0]] * 4, ((0,),))], num_items=4)
    with pytest.raises(ValueError, match="labels"):
        oecluster.consensus([(np.array([0, 1]), np.array([0, 0, 1]))],
                            num_items=4)


def test_a_method_result_of_the_wrong_size_or_shape_is_refused():
    ensemble = [_result([0, 0, 1, 1])] * 2

    def too_short(items, **options):
        return oecluster.ClusteringResult([0, 0], ((0, 1),))

    def column_shaped(items, **options):
        size = items.num_samples
        return oecluster.ClusteringResult([[0]] * size, (tuple(range(size)),))

    with pytest.raises(ValueError, match="consensus items"):
        oecluster.consensus(ensemble, num_items=4, method=too_short)
    with pytest.raises(ValueError, match="one-dimensional"):
        oecluster.consensus(ensemble, num_items=4, method=column_shaped)


@pytest.mark.parametrize("given", [19, 21])
def test_num_items_must_agree_with_a_stability_reference(given):
    spec = oecluster.ClusteringSpec("butina", threshold=0.2)
    stability = oecluster.cluster_stability(spec, _blobs20(), resamples=4,
                                            seed=0)
    with pytest.raises(ValueError, match="disagrees with the ensemble"):
        oecluster.consensus(stability, num_items=given)


@pytest.mark.parametrize("positions, error", [
    (np.array([0, 0]), ValueError),
    (np.array([0, 9]), ValueError),
    (np.array([-1, 1]), ValueError),
    (np.array([], dtype=np.intp), ValueError),
    (np.array([[0, 1]]), TypeError),
    (np.array([0.9, 1.1]), TypeError),
    (np.array([True, False]), TypeError),
])
def test_bad_member_positions_are_refused(positions, error):
    with pytest.raises(error):
        oecluster.consensus([(positions, np.array([0, 0]))], num_items=4)


def test_non_integer_member_labels_are_refused():
    with pytest.raises(TypeError, match="integers"):
        oecluster.consensus([(np.array([0, 1]), np.array([0.0, 1.0]))],
                            num_items=4)


def test_nothing_native_runs_before_validation_finishes(monkeypatch):
    monkeypatch.setattr(
        oecluster.oecluster, "coassociation_distances",
        lambda *args, **kwargs: pytest.fail("the matrix was built anyway"))
    ensemble = [_result([0, 0, 1, 1])] * 2
    for kwargs in ({"threshold": 1.5}, {"num_threads": -1}, {"output": 7},
                   {"noise": "cluster"}, {"method": "murcko"},
                   {"num_items": 19}):
        with pytest.raises((TypeError, ValueError)):
            oecluster.consensus(ensemble, num_items=4, **kwargs)


def test_num_items_is_validated_before_the_members_are_measured():
    # Two defects: a malformed num_items and a member of another size. The
    # order in section 5.1 puts the count first.
    mixed = [_result([0, 0, 1, 1]), _result([0, 0, 1])]
    with pytest.raises(TypeError, match="num_items"):
        oecluster.consensus(mixed, num_items=True)
    with pytest.raises(ValueError, match="labels"):
        oecluster.consensus(mixed)


def test_the_ensemble_kind_is_settled_before_anything_else():
    # Two defects at once: the order in section 5.1 decides which is
    # reported, and the ensemble comes first.
    with pytest.raises(TypeError, match="ClusterStability"):
        oecluster.consensus("not an ensemble", num_items=0, threshold=5.0)
    # Then the count, then the extraction arguments.
    with pytest.raises(ValueError, match="at least 2"):
        oecluster.consensus([(np.array([0]), np.array([0]))], num_items=1,
                            threshold=5.0)


def test_a_foreign_callable_is_not_refused_for_its_name():
    # The four refusals are roster entries, matched by identity: a caller's
    # own function called `murcko` is still a caller's function.
    def murcko(items, **options):
        size = items.num_samples
        return oecluster.ClusteringResult([0] * size, (tuple(range(size)),))

    result = oecluster.consensus([_result([0, 0, 1, 1])] * 2, num_items=4,
                                 method=murcko)
    assert result.num_clusters == 1


def test_an_ineligible_sweep_row_still_contributes():
    # min_clusters leaves the one-cluster row ineligible; consensus uses it
    # all the same, because an ineligible row is still a partition.
    selection = oecluster.select_parameter(
        "butina", _blobs20(), "threshold", [0.2, 0.95], min_clusters=2)
    assert any(not row.eligible for row in selection.rows)
    assert oecluster.consensus(selection).num_partitions == 2


def test_the_package_exports_what_this_task_adds():
    for name in ("consensus", "ConsensusResult"):
        assert name in oecluster.__all__
        assert hasattr(oecluster, name)


def test_a_caller_supplied_none_for_positions_is_refused_before_native_work(
        monkeypatch):
    # The module marks a full member with a private sentinel; a caller's
    # literal None must not borrow that meaning.
    monkeypatch.setattr(
        oecluster.oecluster, "coassociation_distances",
        lambda *args, **kwargs: pytest.fail("the matrix was built anyway"))
    with pytest.raises(TypeError, match="indices"):
        oecluster.consensus([(None, [0, 0, 1, 1])], num_items=4)
