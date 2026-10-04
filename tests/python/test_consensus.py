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


def _slow_strength(distances, labels, num_items):
    """Monti's item and cluster consensus, computed from the distance dict."""
    def co(i, j):
        return 1.0 - distances[(min(i, j), max(i, j))]

    item = []
    for i in range(num_items):
        peers = [j for j in range(num_items)
                 if j != i and labels[j] == labels[i] and labels[i] >= 0]
        item.append(sum(co(i, j) for j in peers) / len(peers)
                    if peers else math.nan)
    cluster = []
    for label in sorted({value for value in labels if value >= 0}):
        members = [i for i in range(num_items) if labels[i] == label]
        pairs = [(i, j) for i in members for j in members if i < j]
        cluster.append(sum(co(i, j) for i, j in pairs) / len(pairs)
                       if pairs else math.nan)
    return item, cluster


def _slow_agreement(ensemble, consensus_labels, num_items, noise):
    """The per-member adjusted Rand index over each member's own positions."""
    scores = []
    for positions, labels in _members(ensemble, num_items):
        observed = [int(consensus_labels[p]) for p in positions]
        scores.append(oecluster.partition_agreement(
            observed, labels, noise=noise).adjusted_rand_index)
    return scores


def _condensed(distances, num_items):
    return np.array([distances[(i, j)]
                     for i in range(num_items)
                     for j in range(i + 1, num_items)])


def _is_swig_vector(value, name):
    """True when ``value`` is the SWIG proxy for a vector, not a tuple.

    oecluster and oefp both register ``std::vector<double>`` in the SWIG
    runtime type table they share, so the proxy class for a returned vector
    belongs to whichever package registered it and import order decides --
    seven test modules here import oefp. Asserting the class identity would
    pin that accident. What the bindings promise is a vector proxy rather
    than the tuple SWIG returns by default, which is the defect these
    assertions exist to catch.
    """
    return not isinstance(value, (tuple, list)) and type(value).__name__ == name


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
    assert _is_swig_vector(components, "IntVector")
    labels = list(components)
    assert labels == [0, 0, 1, 1]
    strength = native.consensus_strength(
        destination, native.IntVector(labels))
    assert _is_swig_vector(strength.item_consensus, "DoubleVector")
    assert _is_swig_vector(strength.cluster_consensus, "DoubleVector")
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
    assert _is_swig_vector(item, "DoubleVector")
    assert _is_swig_vector(cluster, "DoubleVector")
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


# --- consensus strength and agreement ----------------------------------------

def test_strength_matches_the_slow_reference():
    ensemble = [
        _result([0, 0, 1, 1]),
        _result([0, 0, 0, 1]),
        (np.array([0, 1, 3]), np.array([0, 0, 1])),
    ]
    result = oecluster.consensus(ensemble, num_items=4)
    distances = _slow_matrix(ensemble, 4)
    item, cluster = _slow_strength(distances, result.labels.tolist(), 4)
    np.testing.assert_allclose(result.item_consensus, item, equal_nan=True)
    np.testing.assert_allclose(
        [record.cluster_consensus for record in result.records], cluster,
        equal_nan=True)


def test_agreement_matches_the_slow_reference_on_observed_positions():
    ensemble = [
        _result([0, 0, 1, 1]),
        _result([0, 0, 0, 1]),
        (np.array([0, 1, 3]), np.array([0, 0, 1])),
    ]
    result = oecluster.consensus(ensemble, num_items=4)
    expected = _slow_agreement(ensemble, result.labels, 4, "singletons")
    assert result.agreement == pytest.approx(expected, nan_ok=True)
    assert result.mean_agreement == pytest.approx(
        float(np.nanmean(result.agreement)))


def test_mean_agreement_averages_the_defined_entries_only():
    # Under "excluded" a member that calls every item noise leaves the index
    # undefined, and the mean must skip it rather than poison the figure.
    ensemble = [_result([0, 0, 1, 1]), _result([-1, -1, -1, -1])]
    result = oecluster.consensus(ensemble, num_items=4, noise="excluded")
    assert math.isnan(result.agreement[1])
    assert result.mean_agreement == pytest.approx(result.agreement[0])


def test_a_singleton_cluster_has_no_consensus_value():
    # Item 3 agrees with nobody, so it is its own cluster.
    ensemble = [_result([0, 0, 1, 2]), _result([0, 0, 1, 2])]
    result = oecluster.consensus(ensemble, num_items=4)
    singleton = result.labels[3]
    record = next(r for r in result.records if r.label == singleton)
    assert record.size == 1
    assert math.isnan(record.cluster_consensus)
    assert math.isnan(result.item_consensus[3])


def test_records_columns_and_to_table():
    result = oecluster.consensus([_result([0, 0, 1, 1])] * 2, num_items=4)
    assert result.columns == ("label", "size", "cluster_consensus")
    assert [record.label for record in result.records] == [0, 1]
    assert [record.size for record in result.records] == [2, 2]
    table = result.to_table()
    assert table == [tuple(record) for record in result.records]
    assert table is not result.to_table()


def test_repr_heads_the_table_with_the_extraction():
    result = oecluster.consensus([_result([0, 0, 1, 1])] * 2, num_items=4)
    text = repr(result)
    lines = text.splitlines()
    assert lines[0].startswith(
        "ConsensusResult(num_partitions=2, threshold=0.5, num_clusters=2, "
        "mean_agreement=")
    assert "label  size  cluster_consensus" in lines[1]
    assert len(lines) == 2 + len(result.records)

    spec = oecluster.ClusteringSpec("agglomerative", n_clusters=2)
    by_spec = oecluster.consensus([_result([0, 0, 1, 1])] * 2, num_items=4,
                                  method=spec)
    assert "spec=ClusteringSpec('agglomerative'" in repr(by_spec)


def test_item_consensus_refuses_in_place_edits():
    result = oecluster.consensus([_result([0, 0, 1, 1])] * 2, num_items=4)
    with pytest.raises(ValueError):
        result.item_consensus[0] = 0.0


def test_the_package_exports_the_record_type():
    assert "ConsensusRecord" in oecluster.__all__
    assert oecluster.ConsensusRecord is not None


# --- semantics over real ensembles -------------------------------------------

def test_an_identical_ensemble_reproduces_its_partition():
    blobs = oecluster.butina(_blobs20(), threshold=0.2)
    result = oecluster.consensus([blobs, blobs, blobs])
    assert result.num_clusters == blobs.num_clusters
    assert result.agreement == (1.0, 1.0, 1.0)
    assert result.mean_agreement == 1.0
    assert [record.cluster_consensus for record in result.records] == [1.0, 1.0]
    assert float(np.min(result.item_consensus)) == 1.0


def test_a_mutually_disagreeing_ensemble_yields_singletons():
    # Three partitions that share no pair: every co-association is 0.
    rotations = [
        _result([0, 0, 1, 1, 2, 2]),
        _result([0, 1, 2, 0, 1, 2]),
        _result([0, 1, 1, 2, 2, 0]),
    ]
    result = oecluster.consensus(rotations, num_items=6)
    assert result.num_clusters == 6
    assert all(math.isnan(value) for value in result.item_consensus)
    assert all(math.isnan(record.cluster_consensus)
               for record in result.records)


def test_a_single_member_reproduces_itself():
    member = _result([0, 0, 1, 1, 2])
    result = oecluster.consensus([member], num_items=5)
    assert result.labels.tolist() == [0, 0, 1, 1, 2]
    assert result.agreement == (1.0,)


def test_a_member_that_sees_only_noise_contributes_denominators_only():
    ensemble = [_result([0, 0, 1, 1]), _result([-1, -1, -1, -1])]
    result = oecluster.consensus(ensemble, num_items=4)
    # Support for (0,1) halves from 1/1 to 1/2, which still clears 0.5.
    assert result.matrix.squareform()[0][1] == 0.5
    assert result.unobserved_pairs == 0
    assert result.num_clusters == 2


def test_an_item_only_one_member_observed_is_normalized_by_that_member():
    ensemble = [
        (np.array([0, 1, 2]), np.array([0, 0, 1])),
        (np.array([0, 1]), np.array([0, 0])),
        (np.array([0, 1]), np.array([0, 1])),
    ]
    result = oecluster.consensus(ensemble, num_items=3)
    square = result.matrix.squareform()
    assert square[0][1] == pytest.approx(1.0 - 2.0 / 3.0)
    assert square[0][2] == 1.0  # one observer, never together


# --- the threshold boundary --------------------------------------------------

def _support(total, together):
    """`total` full members of two items, `together` of which co-cluster them."""
    return [oecluster.ClusteringResult(
                [0, 0] if index < together else [0, 1],
                ((0, 1),) if index < together else ((0,), (1,)))
            for index in range(total)]


@pytest.mark.parametrize("total, together, threshold, merged", [
    (10, 1, 0.1, True),
    (3, 2, 2 / 3, True),
    (2, 1, 0.5, True),
    (10, 4, 0.5, False),
    (10, 5, 0.5, True),
])
def test_a_support_at_the_threshold_merges(total, together, threshold, merged):
    result = oecluster.consensus(_support(total, together), num_items=2,
                                 threshold=threshold)
    assert (result.num_clusters == 1) is merged


def test_different_denominators_are_resolved_against_one_threshold():
    # Pair (0,1) is seen by 7 members and held together by 3 -> 3/7.
    # Pair (2,3) is seen by 9 members and held together by 4 -> 4/9.
    # A threshold between them must split the first and keep the second.
    members = []
    for index in range(9):
        if index < 7:
            positions = np.array([0, 1, 2, 3])
            labels = np.array([0, 0 if index < 3 else 1,
                               2, 2 if index < 4 else 3])
        else:
            positions = np.array([2, 3])
            labels = np.array([0, 0 if index < 4 else 1])
        members.append((positions, labels))
    result = oecluster.consensus(members, num_items=4, threshold=0.435)
    square = result.matrix.squareform()
    assert square[0][1] == pytest.approx(1.0 - 3.0 / 7.0)
    assert square[2][3] == pytest.approx(1.0 - 4.0 / 9.0)
    assert result.labels[0] != result.labels[1]
    assert result.labels[2] == result.labels[3]


# --- storage, determinism and noise ------------------------------------------

def test_a_memory_mapped_matrix_equals_the_dense_one(tmp_path):
    ensemble = [_result([0, 0, 1, 1]), _result([0, 0, 0, 1])]
    dense = oecluster.consensus(ensemble, num_items=4)
    path = tmp_path / "consensus.bin"
    mapped = oecluster.consensus(ensemble, num_items=4, output=str(path))
    assert isinstance(mapped.matrix.storage, oecluster.MMapStorage)
    np.testing.assert_array_equal(mapped.matrix.condensed,
                                  dense.matrix.condensed)
    np.testing.assert_array_equal(mapped.labels, dense.labels)
    assert path.exists()


def test_a_reused_output_path_does_not_accumulate(tmp_path):
    # MMapStorage keeps a file of the right size, so the kernel zeroes the
    # destination; without that the second run would double every count.
    ensemble = [_result([0, 0, 1, 1]), _result([0, 0, 0, 1])]
    dense = oecluster.consensus(ensemble, num_items=4)
    path = tmp_path / "consensus.bin"
    first = oecluster.consensus(ensemble, num_items=4, output=str(path))
    first_condensed = np.array(first.matrix.condensed)
    # The first mapping is released before the path is opened again: Windows
    # shares a mapped file for reading only, so a second writable open while
    # the first is alive would be denied.
    del first
    second = oecluster.consensus(ensemble, num_items=4, output=str(path))
    np.testing.assert_array_equal(first_condensed, dense.matrix.condensed)
    np.testing.assert_array_equal(second.matrix.condensed,
                                  dense.matrix.condensed)
    del second


@pytest.mark.parametrize("num_threads", [0, 1, 2])
def test_the_result_does_not_depend_on_the_thread_count(num_threads):
    blobs = oecluster.butina(_blobs20(), threshold=0.2)
    singletons = _result(list(range(20)))
    reference = oecluster.consensus([blobs, singletons, blobs])
    result = oecluster.consensus([blobs, singletons, blobs],
                                 num_threads=num_threads)
    np.testing.assert_array_equal(result.matrix.condensed,
                                  reference.matrix.condensed)
    np.testing.assert_array_equal(result.labels, reference.labels)
    assert result.agreement == reference.agreement


def test_noise_mode_changes_agreement_but_never_the_matrix():
    ensemble = [
        oecluster.ClusteringResult([0, 0, 0, 0, 1, 1, -1, -1],
                                   ((0, 1, 2, 3), (4, 5))),
        oecluster.ClusteringResult([0, 0, 0, -1, 1, 1, -1, 2],
                                   ((0, 1, 2), (4, 5), (7,))),
    ]
    results = {mode: oecluster.consensus(ensemble, num_items=8, noise=mode)
               for mode in ("singletons", "grouped", "excluded")}
    baseline = results["singletons"].matrix.condensed
    for mode, result in results.items():
        np.testing.assert_array_equal(result.matrix.condensed, baseline)
        expected = _slow_agreement(ensemble, result.labels, 8, mode)
        assert result.agreement == pytest.approx(expected, nan_ok=True)
    assert len({result.agreement for result in results.values()}) == 3


# --- interoperability --------------------------------------------------------

def test_a_member_with_huge_labels_is_scored_correctly():
    big = 2 ** 40
    ensemble = [
        oecluster.ClusteringResult([big, big, big + 1, big + 1],
                                   ((0, 1), (2, 3))),
        _result([0, 0, 1, 1]),
    ]
    result = oecluster.consensus(ensemble, num_items=4)
    np.testing.assert_allclose(result.matrix.condensed,
                               _condensed(_slow_matrix(ensemble, 4), 4))
    assert result.labels.tolist() == [0, 0, 1, 1]


def test_the_result_is_accepted_by_the_packages_own_consumers():
    # The partition is canonical precisely so these two never refuse it:
    # cluster_report rejects a label that is not its cluster's ordinal, and
    # both reject labels outside the signed 32-bit range.
    blobs = oecluster.butina(_blobs20(), threshold=0.2)

    def sparse_labels(items, **options):
        size = items.num_samples
        half = size // 2
        return oecluster.ClusteringResult(
            [2 ** 40] * half + [5] * (size - half),
            (tuple(range(half)), tuple(range(half, size))))

    result = oecluster.consensus([blobs, blobs], method=sparse_labels)
    assert isinstance(result, oecluster.ClusteringResult)
    assert oecluster.cluster_report(result, _blobs20()).num_clusters == 2
    assert oecluster.partition_agreement(
        result, blobs).adjusted_rand_index == pytest.approx(1.0)


def test_a_consumer_of_result_matrix_can_need_allow_nonmetric():
    # Passing result.matrix back in is not the same as passing the source
    # matrix the test above uses. Only members that observe different items
    # can break the triangle inequality: within one full member co-clustering
    # is transitive and every pair divides by the same member count, so an
    # ensemble of full partitions is always a metric. A ClusterStability's
    # members are partial by construction, and it is the documented example.
    matrix = _blobs20()
    spec = oecluster.ClusteringSpec("butina", threshold=0.2)
    stability = oecluster.cluster_stability(spec, matrix, resamples=8, seed=0)
    result = oecluster.consensus(stability)
    assert result.matrix.facts["metric_probe"] == "violations_found"
    with pytest.raises(ValueError, match="allow_nonmetric"):
        oecluster.cluster_report(result, result.matrix)
    report = oecluster.cluster_report(result, result.matrix,
                                      allow_nonmetric=True)
    assert report.num_clusters == result.num_clusters

    full = oecluster.consensus([oecluster.butina(matrix, threshold=0.2),
                                oecluster.butina(matrix, threshold=0.5)])
    assert full.matrix.facts["metric_probe"] == "no_violations_found"
    assert oecluster.cluster_report(
        full, full.matrix).num_clusters == full.num_clusters


def test_a_sparse_triangle_violation_hides_from_the_probe():
    """Known limitation, pinned deliberately: ``no_violations_found`` is not a
    certification, so a provably non-metric matrix can pass the metric gate.

    ``probe_triangle`` samples a bounded number of triples rather than
    enumerating all ``O(N^3)``, and it can only disprove. The test above is
    the case where a violation is dense enough to be drawn; this is the other
    side of the same coin, and the assertions below describe what the package
    really does rather than what a reader might hope it does. Making the probe
    exhaustive is not the fix -- the sampling trade-off is deliberate and
    predates consensus clustering.
    """
    num_items = 1000
    # The only co-clustered pairs are (0, 1) and (1, 2). Every other pair,
    # (0, 2) included, is unobserved and so takes distance 1.0.
    result = oecluster.consensus([([0, 1], [0, 0]), ([1, 2], [0, 0])],
                                 num_items=num_items)
    storage = result.matrix.storage

    # The violation is exact, not a rounding artefact: 1.0 > 0.0 + 0.0.
    assert storage.Get(0, 1) == 0.0
    assert storage.Get(1, 2) == 0.0
    assert storage.Get(0, 2) == 1.0
    assert storage.Get(0, 2) > storage.Get(0, 1) + storage.Get(1, 2)

    # The probe nevertheless reports nothing. One violating triple out of the
    # ~5e8 distinct inequalities a 1000-item matrix admits is essentially
    # never drawn by a 100,000-triple sample.
    facts = result.matrix.facts
    assert facts["metric_probe"] == "no_violations_found"
    assert facts["probe_violations"] == 0
    assert 0 < facts["probe_sampled"] <= 100000

    # So the gate accepts it: the refusal keys off "violations_found", and the
    # matrix never earned that stamp. No allow_nonmetric is needed, and none
    # would help -- there is nothing to override.
    report = oecluster.cluster_report(result, result.matrix)
    assert report.num_clusters == result.num_clusters
