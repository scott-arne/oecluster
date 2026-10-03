"""Python surface of take and cluster_stability."""
import math

import numpy as np
import oecluster
import oefp
import pytest

# --- fixtures ----------------------------------------------------------------

FP_ROWS = ([set(range(8)) | {16 + i} for i in range(5)]
           + [set(range(32, 40)) | {48 + i} for i in range(5)])


def _batch(rows):
    return oefp.OEFPBatch.from_fingerprints(
        [oefp.OEFP.from_on_bits(64, sorted(on)) for on in rows])


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
