"""Python surface of take and cluster_stability."""
import math
import tracemalloc

import numpy as np
import oecluster
import oefp
import pytest
from openeye import oechem

# --- fixtures ----------------------------------------------------------------

FP_ROWS = ([set(range(8)) | {16 + i} for i in range(5)]
           + [set(range(32, 40)) | {48 + i} for i in range(5)])


def _batch(rows):
    return oefp.OEFPBatch.from_fingerprints(
        [oefp.OEFP.from_on_bits(64, sorted(on)) for on in rows])


def _mols(smiles_list):
    mols = []
    for idx, smi in enumerate(smiles_list):
        mol = oechem.OEGraphMol()
        oechem.OESmilesToMol(mol, smi)
        mol.SetTitle(f"mol{idx}")
        mols.append(mol)
    return mols


def _fps():
    """Two groups of five sharing an eight-bit block, one private bit each."""
    return _batch(FP_ROWS)


def _blobs20():
    """Twenty points, two blobs of ten: intra-blob 0.1, inter-blob 0.9."""
    values = [0.1 if (i < 10) == (j < 10) else 0.9
              for i in range(20) for j in range(i + 1, 20)]
    return oecluster.SymmetricDistanceMatrix.from_condensed(np.array(values))


def _mmap_copy(matrix, path):
    """The same distances behind MMapStorage; pdist(output=) needs molecules."""
    n = matrix.num_samples
    storage = oecluster.MMapStorage(str(path), n)
    condensed = matrix.condensed
    k = 0
    for i in range(n):
        for j in range(i + 1, n):
            storage.Set(i, j, float(condensed[k]))
            k += 1
    return oecluster.SymmetricDistanceMatrix(
        storage, matrix.comparison_name, None, {}, matrix.facts)


def _sparse():
    """Six items, cutoff 0.5: (0,1), (2,3), (4,5) stored at 0.2, (1,4) at 0.4."""
    storage = oecluster.SparseStorage(6, 0.5)
    for i, j, value in ((0, 1, 0.2), (2, 3, 0.2), (4, 5, 0.2), (1, 4, 0.4)):
        storage.Set(i, j, value)
    storage.Finalize()
    return oecluster.SymmetricDistanceMatrix(
        storage, "test", ["a", "b", "c", "d", "e", "f"], {"cutoff": 0.5})


def _violating():
    """d(0,2) = 0.9 > d(0,1) + d(1,2) = 0.2; item 3 sits at 0.5 from everything."""
    return oecluster.SymmetricDistanceMatrix.from_condensed(
        np.array([0.1, 0.9, 0.5, 0.1, 0.5, 0.5]))


def _planted(n, pairs, facts, cutoff=None):
    """A matrix with hand-set pairs and planted facts, built on raw storage."""
    storage = (oecluster.SparseStorage(n, cutoff) if cutoff is not None
               else oecluster.DenseStorage(n))
    for i, j, value in pairs:
        storage.Set(i, j, value)
    storage.Finalize()
    return oecluster.SymmetricDistanceMatrix(storage, "test", None, {}, facts)


def _result(labels):
    """A ClusteringResult from labels alone; -1 is noise."""
    labels = np.asarray(labels)
    clusters = [tuple(int(p) for p in np.flatnonzero(labels == label))
                for label in sorted(set(labels.tolist())) if label >= 0]
    return oecluster.ClusteringResult(labels, clusters)


def _fixed(labels_for_size, calls=None):
    """A foreign clusterer whose partition depends only on the item count."""
    def clusterer(items, **options):
        if calls is not None:
            calls.append(options)
        size = (items.num_samples
                if isinstance(items, oecluster.SymmetricDistanceMatrix)
                else items.size)
        return _result(labels_for_size(size))
    return clusterer


def _thirds(size):
    """Three interleaved clusters with every seventh item noise."""
    return [-1 if p % 7 == 0 else p % 3 for p in range(size)]


# A hand-built reference over twenty items: two clusters, a singleton, noise.
HAND_REFERENCE = [0] * 8 + [1] * 7 + [-1] * 3 + [2] + [-1]


# --- slow reference implementation ------------------------------------------

def _slow_stability(reference_labels, indices, labels):
    """Hennig's statistics with Python sets; shares nothing with the module."""
    reference_labels = [int(label) for label in reference_labels]
    record_labels = sorted({label for label in reference_labels if label >= 0})
    jaccard = np.full((len(record_labels), len(indices)), math.nan)
    for r, (positions, resample) in enumerate(zip(indices, labels)):
        ref = [reference_labels[int(p)] for p in positions]
        resample = [int(label) for label in resample]
        for i, k in enumerate(record_labels):
            members = {p for p, label in enumerate(ref) if label == k}
            if not members:
                continue
            best = 0.0
            for d in {label for label in resample if label >= 0}:
                cluster = {p for p, label in enumerate(resample) if label == d}
                best = max(best, len(members & cluster) / len(members | cluster))
            jaccard[i, r] = best
    records = []
    for i, k in enumerate(record_labels):
        present = [value for value in jaccard[i] if not math.isnan(value)]
        count = len(present)
        records.append((
            k, reference_labels.count(k),
            sum(present) / count if count else math.nan,
            sum(value < 0.5 for value in present) / count if count else math.nan,
            sum(value > 0.75 for value in present) / count if count else math.nan,
            count))
    return jaccard, records


def _slow_agreement(reference_labels, indices, labels, noise):
    return [oecluster.partition_agreement(
                [int(reference_labels[int(p)]) for p in positions],
                [int(label) for label in resample],
                noise=noise).adjusted_rand_index
            for positions, resample in zip(indices, labels)]


def _assert_matches_slow(stability, reference_labels, noise="singletons"):
    jaccard, records = _slow_stability(
        reference_labels, stability.indices, stability.labels)
    np.testing.assert_allclose(stability.jaccard, jaccard, equal_nan=True)
    assert len(stability.records) == len(records)
    for got, want in zip(stability.records, records):
        assert (got.label, got.size, got.evaluated) == (want[0], want[1], want[5])
        for value, expected in zip(got[2:5], want[2:5]):
            assert ((math.isnan(value) and math.isnan(expected))
                    or value == pytest.approx(expected))
    assert stability.agreement == pytest.approx(
        _slow_agreement(reference_labels, stability.indices, stability.labels,
                        noise), nan_ok=True)


# --- native bindings ---------------------------------------------------------

def test_native_take_pairs_gathers_in_index_order():
    source = oecluster.DenseStorage(4)
    for i in range(4):
        for j in range(i + 1, 4):
            source.Set(i, j, 10 * i + j)
    destination = oecluster.DenseStorage(3)
    oecluster.oecluster.take_pairs(
        source, oecluster.oecluster.SizeTVector([3, 0, 2]), destination)
    assert destination.Get(0, 1) == source.Get(3, 0)
    assert destination.Get(0, 2) == source.Get(3, 2)
    assert destination.Get(1, 2) == source.Get(0, 2)


def test_native_take_pairs_refusal_is_a_runtime_error():
    source = oecluster.DenseStorage(3)
    with pytest.raises(RuntimeError, match="more than once"):
        oecluster.oecluster.take_pairs(
            source, oecluster.oecluster.SizeTVector([1, 1]),
            oecluster.DenseStorage(2))


def test_native_take_fingerprints_returns_an_oefp_native_batch():
    batch = _fps()
    native = oecluster.oecluster.take_fingerprints(
        batch, oecluster.oecluster.SizeTVector([7, 0]))
    assert type(native).__module__ == "oefp._native"
    subset = oefp.OEFPBatch._from_native(native)
    np.testing.assert_array_equal(subset.words, batch.words[[7, 0]])
    assert subset.spec == batch.spec


# --- take: matrices ----------------------------------------------------------

@pytest.mark.parametrize("kind", ["dense", "mmap"])
def test_take_follows_the_given_order(kind, tmp_path):
    matrix = _blobs20()
    if kind == "mmap":
        matrix = _mmap_copy(matrix, tmp_path / "blobs.bin")
    order = [15, 2, 7, 11]
    subset = oecluster.take(matrix, order)
    assert subset.num_samples == 4
    assert isinstance(subset.storage, oecluster.DenseStorage)
    np.testing.assert_array_equal(
        subset.squareform(), matrix.squareform()[np.ix_(order, order)])


def test_take_sparse_keeps_the_cutoff_and_stored_pairs_only():
    subset = oecluster.take(_sparse(), [4, 1, 0, 3])
    assert isinstance(subset.storage, oecluster.SparseStorage)
    assert subset.storage.Cutoff() == 0.5
    assert subset.storage.Get(1, 2) == 0.2    # source (1, 0)
    assert subset.storage.Get(0, 1) == 0.4    # source (4, 1)
    assert subset.storage.Get(0, 3) == 0.0    # source (4, 3): never stored
    assert subset.storage.Get(2, 3) == 0.0    # source (0, 3): never stored
    assert subset.labels == ["e", "b", "a", "d"]
    assert subset.params == {"cutoff": 0.5}


def test_take_preserves_name_params_and_slices_labels():
    matrix = oecluster.SymmetricDistanceMatrix.from_condensed(
        np.array([0.1, 0.2, 0.3]), labels=["x", "y", "z"],
        comparison_name="tanimoto", params={"bits": 1024})
    subset = oecluster.take(matrix, [2, 0])
    assert subset.comparison_name == "tanimoto"
    assert subset.params == {"bits": 1024}
    assert subset.labels == ["z", "x"]
    unlabelled = oecluster.take(_blobs20(), [1, 2])
    assert list(unlabelled.labels) == []


def test_take_thread_count_does_not_change_the_subset():
    matrix = _blobs20()
    order = list(range(19, -1, -1))[:13]
    one = oecluster.take(matrix, order, num_threads=1)
    auto = oecluster.take(matrix, order, num_threads=0)
    np.testing.assert_array_equal(one.condensed, auto.condensed)


# --- take: facts -------------------------------------------------------------

def test_take_inherits_the_comparison_facts():
    matrix = _planted(4, [(0, 1, 0.1), (0, 2, 0.2), (0, 3, 0.3), (1, 2, 0.1),
                          (1, 3, 0.2), (2, 3, 0.1)],
                      {"is_distance": True, "zero_self": True, "triangle": True,
                       "data_integrity": "complete"})
    facts = oecluster.take(matrix, [1, 3]).facts
    assert (facts["is_distance"], facts["zero_self"], facts["triangle"]) == (
        True, True, True)
    assert facts["data_integrity"] == "complete"
    assert facts["metric_probe"] == "not_run"


def test_take_of_every_item_inherits_every_fact_including_a_violation():
    bad = _violating()
    assert bad.facts["metric_probe"] == "violations_found"
    permuted = oecluster.take(bad, [3, 2, 1, 0])
    assert permuted.facts == bad.facts
    with pytest.raises(ValueError, match="triangle inequality"):
        oecluster.butina(permuted, threshold=0.3)


def test_take_reprobes_a_proper_subset_of_a_probed_source():
    bad = _violating()
    cleared = oecluster.take(bad, [0, 1, 3])
    assert cleared.facts["metric_probe"] == "no_violations_found"
    assert cleared.facts["probe_violations"] == 0
    assert oecluster.butina(cleared, threshold=0.3).num_samples == 3
    kept = oecluster.take(bad, [2, 0, 1])
    assert kept.facts["metric_probe"] == "violations_found"
    assert kept.facts["probe_violations"] == 1


def test_take_leaves_an_unprobed_source_unprobed():
    matrix = oecluster.SymmetricDistanceMatrix.from_condensed(
        _violating().condensed, probe_triples=0)
    assert matrix.facts["metric_probe"] == "not_run"
    assert oecluster.take(matrix, [0, 1, 2]).facts["metric_probe"] == "not_run"


def test_take_sparse_subset_inherits_the_probe_fields():
    planted = {"metric_probe": "violations_found", "probe_violations": 3,
               "probe_sampled": 9}
    matrix = _planted(5, [(0, 1, 0.2), (3, 4, 0.1)], planted, cutoff=0.5)
    facts = oecluster.take(matrix, [4, 3, 0]).facts
    assert {key: facts[key] for key in planted} == planted


def test_take_remeasures_nan_presence_on_a_dense_subset():
    matrix = _planted(4, [(0, 1, math.nan), (0, 2, 0.1), (0, 3, 0.1),
                          (1, 2, 0.1), (1, 3, 0.1), (2, 3, 0.1)],
                      {"data_integrity": "nan_present"})
    assert oecluster.take(matrix, [2, 3]).facts["data_integrity"] == "complete"
    assert oecluster.take(matrix, [0, 1, 2]).facts["data_integrity"] == "nan_present"


def test_take_remeasures_nan_presence_on_a_sparse_subset_without_densifying(
        monkeypatch):
    matrix = _planted(4, [(0, 1, math.nan), (2, 3, 0.2)],
                      {"data_integrity": "nan_present"}, cutoff=0.5)
    monkeypatch.setattr(
        oecluster.SymmetricDistanceMatrix, "condensed",
        property(lambda self: pytest.fail("a sparse subset was densified")))
    assert oecluster.take(matrix, [2, 3]).facts["data_integrity"] == "complete"
    assert oecluster.take(matrix, [0, 1, 3]).facts["data_integrity"] == "nan_present"


def test_take_keeps_subset_scored_integrity():
    matrix = _planted(3, [(0, 1, 0.1), (0, 2, 0.2), (1, 2, 0.3)],
                      {"data_integrity": "subset_scored"})
    assert oecluster.take(matrix, [0, 2]).facts["data_integrity"] == "subset_scored"


# --- take: refusals ----------------------------------------------------------

@pytest.mark.parametrize("indices, error", [
    ([0, 0], ValueError),
    ([0, 99], ValueError),
    ([-1], ValueError),
    ([], ValueError),
    ([2 ** 63], ValueError),
    ([[0, 1]], TypeError),
    ([0.9, 1.1], TypeError),
    ([True, False], TypeError),
    (["0", "1"], TypeError),
])
def test_take_refuses_bad_indices_before_any_native_call(indices, error,
                                                          monkeypatch):
    monkeypatch.setattr(oecluster.oecluster, "take_pairs",
                        lambda *args: pytest.fail("native gather ran"))
    with pytest.raises(error):
        oecluster.take(_blobs20(), indices)


def test_take_float_positions_are_refused_not_truncated():
    with pytest.raises(TypeError, match="coerced"):
        oecluster.take(_blobs20(), [0.9, 1.1])


@pytest.mark.parametrize("num_threads, error", [
    (1.5, TypeError), (True, TypeError), (-1, ValueError)])
def test_take_checks_num_threads(num_threads, error):
    with pytest.raises(error):
        oecluster.take(_blobs20(), [0, 1], num_threads=num_threads)


def test_take_refuses_other_kinds():
    with pytest.raises(TypeError, match="SymmetricDistanceMatrix or an oefp.OEFPBatch"):
        oecluster.take([[0.0, 0.1], [0.1, 0.0]], [0])


# --- take: fingerprint batches ----------------------------------------------

def test_take_batch_copies_rows_in_order_and_keeps_the_spec():
    batch = _fps()
    order = [7, 0, 3]
    subset = oecluster.take(batch, order)
    assert subset.size == 3
    np.testing.assert_array_equal(subset.words, batch.words[order])
    np.testing.assert_array_equal(subset.popcounts, batch.popcounts[order])
    assert subset.spec == batch.spec


def test_take_batch_clusters_like_a_batch_built_from_the_same_fingerprints():
    order = [9, 8, 1, 0, 5]
    subset = oecluster.take(_fps(), order)
    rebuilt = _batch([FP_ROWS[i] for i in order])
    np.testing.assert_array_equal(
        oecluster.bitbirch(subset, threshold=0.6).labels,
        oecluster.bitbirch(rebuilt, threshold=0.6).labels)


def test_take_batch_refuses_bad_indices_before_any_native_call(monkeypatch):
    monkeypatch.setattr(oecluster.oecluster, "take_fingerprints",
                        lambda *args: pytest.fail("native gather ran"))
    with pytest.raises(ValueError):
        oecluster.take(_fps(), [0, 10])
    with pytest.raises(TypeError):
        oecluster.take(_fps(), [0.0])


# --- cluster_stability: against the slow reference --------------------------

def test_stability_matches_the_slow_reference_for_butina_on_blobs():
    spec = oecluster.ClusteringSpec("butina", threshold=0.2)
    stability = oecluster.cluster_stability(spec, _blobs20(), resamples=25)
    _assert_matches_slow(stability, stability.reference.labels)


def test_stability_matches_the_slow_reference_for_bitbirch_on_fingerprints():
    spec = oecluster.ClusteringSpec("bitbirch", threshold=0.6)
    stability = oecluster.cluster_stability(spec, _fps(), resamples=25)
    _assert_matches_slow(stability, stability.reference.labels)


def test_stability_matches_the_slow_reference_on_a_hand_built_reference():
    stability = oecluster.cluster_stability(
        _fixed(_thirds), _blobs20(), resamples=30,
        reference=_result(HAND_REFERENCE))
    assert [record.label for record in stability.records] == [0, 1, 2]
    assert [record.size for record in stability.records] == [8, 7, 1]
    _assert_matches_slow(stability, HAND_REFERENCE)


def test_stability_keeps_every_score_in_the_row_of_its_own_label():
    reference = [3, 3, 3, 3, 7, 7, 12, 12, 12, 12, -1, -1]
    matrix = oecluster.SymmetricDistanceMatrix.from_condensed(
        np.full(12 * 11 // 2, 0.5))
    stability = oecluster.cluster_stability(
        _fixed(lambda size: [p % 2 for p in range(size)]), matrix,
        resamples=20, seed=0, reference=_result(reference))
    assert [record.label for record in stability.records] == [3, 7, 12]
    # Under this seed some resample holds neither item 4 nor 5, so the middle
    # row is NaN there while the later cluster still scores in its own row.
    missing = [r for r in range(20)
               if not ({4, 5} & set(stability.indices[r].tolist()))]
    assert missing
    for r in missing:
        assert math.isnan(stability.jaccard[1, r])
        assert not math.isnan(stability.jaccard[2, r])
    _assert_matches_slow(stability, reference)


# --- cluster_stability: validation ------------------------------------------

@pytest.mark.parametrize("kwargs, error", [
    ({"items": [[0.0, 0.1], [0.1, 0.0]]}, TypeError),
    ({"items": oecluster.SymmetricDistanceMatrix.from_condensed(np.array([]))},
     ValueError),
    ({"resamples": True}, TypeError),
    ({"resamples": 1.5}, TypeError),
    ({"resamples": 0}, ValueError),
    ({"fraction": True}, TypeError),
    ({"fraction": "0.5"}, TypeError),
    ({"fraction": 0}, ValueError),
    ({"fraction": 1.5}, ValueError),
    ({"fraction": math.nan}, ValueError),
    ({"seed": True}, TypeError),
    ({"seed": 1.5}, TypeError),
    ({"seed": -1}, ValueError),
    ({"keep_partitions": 1}, TypeError),
    ({"keep_partitions": "yes"}, TypeError),
    ({"num_threads": 1.5}, TypeError),
    ({"num_threads": -1}, ValueError),
    ({"reference": [0] * 20}, TypeError),
    ({"reference": _result([0] * 19)}, ValueError),
])
def test_validation_refuses_before_anything_runs(kwargs, error):
    calls = []
    arguments = {"items": _blobs20(), "resamples": 3}
    arguments.update(kwargs)
    with pytest.raises(error):
        oecluster.cluster_stability(_fixed(_thirds, calls), **arguments)
    assert calls == []


def test_items_error_points_a_raw_sequence_at_pdist():
    with pytest.raises(TypeError, match="pdist"):
        oecluster.cluster_stability("butina", [[0.0, 0.1], [0.1, 0.0]])


def test_a_bad_noise_mode_is_refused_by_partition_agreement():
    with pytest.raises(ValueError):
        oecluster.cluster_stability(_fixed(_thirds), _blobs20(), resamples=1,
                                    reference=_result(HAND_REFERENCE),
                                    noise="cluster")


def test_a_roster_name_and_a_callable_are_wrapped_as_specs():
    by_name = oecluster.cluster_stability("hdbscan", _blobs20(), resamples=2)
    assert by_name.spec == oecluster.ClusteringSpec("hdbscan")
    by_callable = oecluster.cluster_stability(
        _fixed(_thirds), _blobs20(), resamples=2)
    assert by_callable.spec.name == "clusterer"


# --- the result object -------------------------------------------------------

def test_columns_table_and_repr():
    stability = oecluster.cluster_stability(
        _fixed(_thirds), _blobs20(), resamples=4,
        reference=_result(HAND_REFERENCE))
    assert stability.columns == ("label", "size", "mean_jaccard", "dissolved",
                                 "recovered", "evaluated")
    table = stability.to_table()
    assert table == [tuple(record) for record in stability.records]
    assert table is not stability.to_table()
    text = repr(stability)
    assert text.startswith("ClusterStability(spec=ClusteringSpec('clusterer'), "
                           "resamples=4, fraction=0.5, mean_jaccard=")
    assert "label  size  mean_jaccard  dissolved  recovered  evaluated" in text
    assert len(text.splitlines()) == 2 + len(stability.records)


def test_repr_prints_nan_for_an_unevaluated_record():
    stability = oecluster.cluster_stability(
        _fixed(lambda size: [-1] * size), _blobs20(), resamples=1,
        reference=_result([0] * 19 + [1]), seed=0, fraction=0.5)
    # Seed 0 leaves item 19 out of the one half-size draw, so the singleton
    # cluster 1 is never evaluated.
    assert stability.records[1].evaluated == 0
    lines = repr(stability).splitlines()
    assert lines[-1].split()[2:5] == ["nan", "nan", "nan"]
    assert "mean_agreement=" in lines[0]


def test_result_is_read_only():
    stability = oecluster.cluster_stability(
        _fixed(_thirds), _blobs20(), resamples=2,
        reference=_result(HAND_REFERENCE))
    with pytest.raises(AttributeError):
        stability.resamples = 3
    with pytest.raises(AttributeError):
        stability._jaccard = None
    with pytest.raises(ValueError):
        stability.jaccard[0, 0] = 0.0


def test_package_exports_the_four_names():
    for name in ("take", "cluster_stability", "ClusterStability",
                 "ClusterStabilityRecord"):
        assert name in oecluster.__all__
        assert hasattr(oecluster, name)


# --- cluster_stability: a callable that aliases the reference ---------------

def test_a_callable_that_returns_its_reference_object_again_is_refused():
    shared = _result([0] * 20)

    def reusing(items, **options):
        # One object handed back on every call, relabelled in place: the
        # reference would silently become the last resample.
        shared.labels[:] = 1
        return shared

    with pytest.raises(ValueError, match="reference result object"):
        oecluster.cluster_stability(reusing, _blobs20(), resamples=2,
                                    fraction=1)


def test_a_reference_relabelled_during_resampling_is_refused():
    reference = _result([0] * 10 + [1] * 10)

    def relabelling(items, **options):
        # Mutates the caller's reference through a retained alias without
        # returning it.
        reference.labels[:] = 0
        return _result([p % 2 for p in range(items.num_samples)])

    with pytest.raises(ValueError, match="changed the reference"):
        oecluster.cluster_stability(relabelling, _blobs20(), resamples=2,
                                    reference=reference)


def test_column_shaped_resample_labels_are_refused():
    def columns(items, **options):
        size = items.num_samples
        return oecluster.ClusteringResult([[p % 2] for p in range(size)],
                                          ((0,), (1,)))

    with pytest.raises(ValueError, match="one-dimensional"):
        oecluster.cluster_stability(columns, _blobs20(), resamples=1,
                                    reference=_result(HAND_REFERENCE))


def test_column_shaped_reference_labels_are_refused():
    calls = []
    reference = oecluster.ClusteringResult([[p % 2] for p in range(20)],
                                           ((0,), (1,)))
    with pytest.raises(ValueError, match="one-dimensional"):
        oecluster.cluster_stability(_fixed(_thirds, calls), _blobs20(),
                                    resamples=1, reference=reference)
    assert calls == []


# --- cluster_stability: semantics -------------------------------------------

def test_a_perfectly_stable_clustering_scores_one_everywhere():
    spec = oecluster.ClusteringSpec("butina", threshold=0.2)
    stability = oecluster.cluster_stability(spec, _blobs20(), resamples=30)
    assert np.nanmin(stability.jaccard) == 1.0
    for record in stability.records:
        assert (record.mean_jaccard, record.dissolved, record.recovered) == (
            1.0, 0.0, 1.0)
    assert stability.mean_jaccard == 1.0
    for positions, agreement in zip(stability.indices, stability.agreement):
        if (positions < 10).any() and (positions >= 10).any():
            assert agreement == 1.0


def test_equal_seeds_reproduce_and_different_seeds_differ():
    spec = oecluster.ClusteringSpec("butina", threshold=0.2)
    first = oecluster.cluster_stability(spec, _blobs20(), resamples=5, seed=3)
    again = oecluster.cluster_stability(spec, _blobs20(), resamples=5, seed=3)
    other = oecluster.cluster_stability(spec, _blobs20(), resamples=5, seed=4)
    for a, b in zip(first.indices, again.indices):
        np.testing.assert_array_equal(a, b)
    for a, b in zip(first.labels, again.labels):
        np.testing.assert_array_equal(a, b)
    np.testing.assert_array_equal(first.jaccard, again.jaccard)
    assert first.agreement == again.agreement
    assert any(not np.array_equal(a, b)
               for a, b in zip(first.indices, other.indices))


def test_seed_none_draws_fresh_entropy_and_records_none(monkeypatch):
    seeds = []
    real_default_rng = np.random.default_rng

    def spying_default_rng(seed=None):
        seeds.append(seed)
        return real_default_rng(seed)

    monkeypatch.setattr(oecluster._stability.np.random, "default_rng",
                        spying_default_rng)
    # The patch is process-wide, so unrelated callers may also record seeds;
    # exactly one None means the resampling generator got fresh entropy.
    stability = oecluster.cluster_stability(
        _fixed(_thirds), _blobs20(), resamples=4, seed=None,
        reference=_result(HAND_REFERENCE))
    assert seeds.count(None) == 1
    assert stability.seed is None
    assert len(stability.indices) == 4
    assert stability.jaccard.shape == (len(set(HAND_REFERENCE) - {-1}), 4)


def test_full_fraction_reproduces_a_deterministic_reference():
    spec = oecluster.ClusteringSpec("butina", threshold=0.2)
    stability = oecluster.cluster_stability(spec, _blobs20(), resamples=3,
                                            fraction=1)
    assert stability.fraction == 1.0
    for positions, labels in zip(stability.indices, stability.labels):
        np.testing.assert_array_equal(positions, np.arange(20))
        np.testing.assert_array_equal(labels, stability.reference.labels)
    assert stability.agreement == (1.0, 1.0, 1.0)
    assert np.all(stability.jaccard == 1.0)


def test_full_fraction_measures_the_distance_to_a_supplied_reference():
    spec = oecluster.ClusteringSpec("butina", threshold=0.2)
    one_cluster = oecluster.butina(_blobs20(), threshold=0.95)
    assert one_cluster.num_clusters == 1
    stability = oecluster.cluster_stability(spec, _blobs20(), resamples=2,
                                            fraction=1, reference=one_cluster)
    assert np.all(stability.jaccard == 0.5)
    assert all(value < 1.0 for value in stability.agreement)


@pytest.mark.parametrize("fraction", (0.5, 0.1))
def test_subset_size_is_at_least_two(fraction):
    matrix = oecluster.SymmetricDistanceMatrix.from_condensed(
        np.array([0.1, 0.2, 0.3]))
    stability = oecluster.cluster_stability(
        _fixed(lambda size: [0] * size), matrix, resamples=3, fraction=fraction)
    assert all(positions.size == 2 for positions in stability.indices)


def test_noise_mode_changes_agreement_but_never_jaccard():
    reference = _result([0, 0, 0, 0, 1, 1, -1, -1])
    matrix = oecluster.SymmetricDistanceMatrix.from_condensed(np.full(28, 0.5))
    clusterer = _fixed(lambda size: [0, 0, 0, -1, 1, 1, -1, 2])
    results = {noise: oecluster.cluster_stability(
                   clusterer, matrix, resamples=2, fraction=1,
                   reference=reference, noise=noise)
               for noise in ("singletons", "grouped", "excluded")}
    for noise, stability in results.items():
        np.testing.assert_array_equal(stability.jaccard,
                                      results["singletons"].jaccard)
        _assert_matches_slow(stability, reference.labels, noise=noise)
    assert len({stability.agreement for stability in results.values()}) == 3


def test_excluded_noise_leaves_agreement_undefined_when_nothing_survives():
    reference = _result(HAND_REFERENCE)
    all_noise = _fixed(lambda size: [-1] * size)
    stability = oecluster.cluster_stability(
        all_noise, _blobs20(), resamples=4, reference=reference,
        noise="excluded")
    assert all(math.isnan(value) for value in stability.agreement)
    assert math.isnan(stability.mean_agreement)
    assert np.all(stability.jaccard[~np.isnan(stability.jaccard)] == 0.0)


def test_mean_agreement_averages_the_defined_entries_only():
    calls = []
    def alternating(items, **options):
        calls.append(options)
        size = items.num_samples
        labels = [-1] * size if len(calls) % 2 else [p % 2 for p in range(size)]
        return _result(labels)
    stability = oecluster.cluster_stability(
        alternating, _blobs20(), resamples=4, reference=_result(HAND_REFERENCE),
        noise="excluded")
    defined = [value for value in stability.agreement if not math.isnan(value)]
    assert len(defined) == 2
    assert stability.mean_agreement == pytest.approx(sum(defined) / 2)


def test_huge_labels_are_encoded_not_allocated():
    big = 2 ** 40
    clusterer = _fixed(lambda size: [big + (p % 2) for p in range(size)])
    stability = oecluster.cluster_stability(
        clusterer, _blobs20(), resamples=3, reference=_result(HAND_REFERENCE))
    jaccard, _ = _slow_stability(HAND_REFERENCE, stability.indices,
                                 stability.labels)
    np.testing.assert_allclose(stability.jaccard, jaccard, equal_nan=True)
    small = [np.asarray(labels) - big for labels in stability.labels]
    assert stability.agreement == pytest.approx(
        _slow_agreement(HAND_REFERENCE, stability.indices, small, "singletons"))


def test_singleton_heavy_partitions_need_no_contingency_table():
    size = 20000
    storage = oecluster.SparseStorage(size, 0.5)
    storage.Finalize()
    matrix = oecluster.SymmetricDistanceMatrix(storage, "test")
    reference = _result(np.arange(size))
    clusterer = _fixed(lambda count: list(range(count)))
    # The tracer may belong to the runner, so measure a delta and stop only a
    # tracer this test started.
    was_tracing = tracemalloc.is_tracing()
    if not was_tracing:
        tracemalloc.start()
    try:
        tracemalloc.reset_peak()
        baseline, _ = tracemalloc.get_traced_memory()
        stability = oecluster.cluster_stability(
            clusterer, matrix, resamples=1, reference=reference)
        _, peak = tracemalloc.get_traced_memory()
    finally:
        if not was_tracing:
            tracemalloc.stop()
    assert peak - baseline < 32 * 1024 * 1024
    assert stability.jaccard.shape == (size, 1)
    assert np.nansum(stability.jaccard) == size // 2


def test_a_resample_of_the_wrong_size_is_refused_by_index():
    clusterer = _fixed(lambda size: [0, 0, 1])
    with pytest.raises(ValueError, match="resample 0"):
        oecluster.cluster_stability(clusterer, _blobs20(), resamples=2,
                                    reference=_result(HAND_REFERENCE))


def test_an_internal_reference_of_the_wrong_size_is_refused_before_resampling():
    calls = []
    clusterer = _fixed(lambda size: [0, 0, 1], calls)
    with pytest.raises(ValueError, match="20 reference items"):
        oecluster.cluster_stability(clusterer, _blobs20(), resamples=2)
    assert len(calls) == 1


def test_keep_partitions_false_drops_only_the_partitions():
    spec = oecluster.ClusteringSpec("butina", threshold=0.2)
    kept = oecluster.cluster_stability(spec, _blobs20(), resamples=4, seed=1)
    dropped = oecluster.cluster_stability(spec, _blobs20(), resamples=4, seed=1,
                                          keep_partitions=False)
    assert dropped.indices is None and dropped.labels is None
    np.testing.assert_array_equal(dropped.jaccard, kept.jaccard)
    assert dropped.agreement == kept.agreement
    assert dropped.records == kept.records
    assert dropped.mean_jaccard == kept.mean_jaccard
    assert dropped.mean_agreement == kept.mean_agreement
    np.testing.assert_array_equal(dropped.reference.labels,
                                  kept.reference.labels)
    assert dropped.spec == kept.spec
    assert (dropped.resamples, dropped.fraction, dropped.seed) == (
        kept.resamples, kept.fraction, kept.seed)


def test_an_all_noise_reference_has_no_records():
    stability = oecluster.cluster_stability(
        _fixed(lambda size: [p % 2 for p in range(size)]), _blobs20(),
        resamples=3, reference=_result([-1] * 20))
    assert stability.records == ()
    assert stability.jaccard.shape == (0, 3)
    assert math.isnan(stability.mean_jaccard)
    assert len(stability.agreement) == 3
    assert stability.agreement == pytest.approx(
        _slow_agreement([-1] * 20, stability.indices, stability.labels,
                        "singletons"), nan_ok=True)
    assert stability.mean_agreement == pytest.approx(
        float(np.nanmean(stability.agreement)))
    assert "ClusterStability(" in repr(stability)


def test_a_supplied_reference_is_not_rerun():
    calls = []
    clusterer = _fixed(_thirds, calls)
    oecluster.cluster_stability(clusterer, _blobs20(), resamples=5,
                                reference=_result(HAND_REFERENCE))
    assert len(calls) == 5
    calls.clear()
    oecluster.cluster_stability(clusterer, _blobs20(), resamples=5)
    assert len(calls) == 6


def test_a_raising_resample_propagates():
    calls = []
    def flaky(items, **options):
        calls.append(options)
        if len(calls) == 3:
            raise RuntimeError("boom")
        return _result([0] * items.num_samples)
    with pytest.raises(RuntimeError, match="boom"):
        oecluster.cluster_stability(flaky, _blobs20(), resamples=4)


def test_retained_labels_are_snapshots_of_each_resample():
    calls = []
    shared = _result([0] * 10)

    def reusing(items, **options):
        # One result object relabelled in place between calls; the retained
        # partitions must not change after the fact.
        calls.append(options)
        shared.labels[:] = len(calls) % 2
        return shared

    stability = oecluster.cluster_stability(
        reusing, _blobs20(), resamples=2, reference=_result(HAND_REFERENCE))
    np.testing.assert_array_equal(stability.labels[0],
                                  np.ones(10, dtype=np.intp))
    np.testing.assert_array_equal(stability.labels[1],
                                  np.zeros(10, dtype=np.intp))


DESCRIPTOR_COLUMNS = ["FractionCsp3", "MolecularWeight"]


def _descriptor_matrix(missing):
    """Water has no carbon, so FractionCsp3 is present-and-NaN for it."""
    mols = _mols(["O", "CCO", "CCCO", "c1ccccc1", "CC(=O)O"])
    return oecluster.pdist(mols, "descriptor", missing=missing,
                           metric="euclidean", columns=DESCRIPTOR_COLUMNS)


def test_take_keeps_an_ignore_scored_descriptor_subset_subset_scored():
    source = _descriptor_matrix("ignore")
    assert source.params["missing"] == "ignore"
    assert source.facts["data_integrity"] == "nan_present"
    subset = oecluster.take(source, [1, 2, 3, 4])
    assert subset.facts["data_integrity"] == "subset_scored"
    with pytest.raises(ValueError):
        oecluster.k_medoids(subset, n_clusters=2)


def test_take_completes_a_propagated_descriptor_subset():
    source = _descriptor_matrix("propagate")
    assert source.params["missing"] == "propagate"
    assert source.facts["data_integrity"] == "nan_present"
    subset = oecluster.take(source, [1, 2, 3, 4])
    assert subset.facts["data_integrity"] == "complete"
    assert oecluster.k_medoids(subset, n_clusters=2).num_samples == 4


def test_take_treats_an_unrecorded_descriptor_policy_as_subset_scored():
    storage = oecluster.DenseStorage(3)
    storage.Set(0, 1, math.nan)
    storage.Set(0, 2, 0.4)
    storage.Set(1, 2, 0.5)
    source = oecluster.SymmetricDistanceMatrix(
        storage, "descriptor", None, {"comparison_type": "descriptor"},
        {"data_integrity": "nan_present"})
    assert oecluster.take(source, [1, 2]).facts["data_integrity"] == "subset_scored"
