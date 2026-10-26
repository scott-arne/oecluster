"""Parity tests for native BitBirch clustering."""

from __future__ import annotations

import importlib
import math
import sys
from pathlib import Path

import numpy as np
import oecluster
import oefp
import pytest

BITBIRCH_REPO = Path("/Users/johnss51/Development/python/bitbirch")
ISIM_REPO = Path("/Users/johnss51/Development/python/iSIM")

QUALITY_N = 3000
DETERMINISM_N = 5000

# Calibrated fast-vs-strict quality tolerances (see plan Task 6).
FAST_QUALITY_TOL = {
    "num_clusters_rel": 0.15,        # observed worst-case ~0.7%
    "intra_distance_abs": 0.05,      # observed: fast ~0.008–0.009 better (lower)
    "silhouette_abs": 0.05,          # observed worst-case ~0.012 lower
    "coverage_abs": 0.05,            # observed: no meaningful delta
}


def _load_reference_bitbirch():
    """Import the upstream BitBirch reference implementation.

    The reference is published as a source checkout rather than a distribution,
    so these parity tests can only run where the sibling repositories exist.
    Anywhere else, including CI and fresh clones, they skip.

    :returns: The reloaded ``bitbirch.bitbirch`` reference module.
    """
    for repo in (str(ISIM_REPO), str(BITBIRCH_REPO)):
        if repo not in sys.path:
            sys.path.insert(0, repo)
    try:
        module = importlib.import_module("bitbirch.bitbirch")
    except ImportError:
        pytest.skip(f"BitBirch reference checkout not found at {BITBIRCH_REPO}")
    return importlib.reload(module)


def _reference_result(
    bits,
    *,
    threshold=0.65,
    branching_factor=50,
    merge_criterion="diameter",
    tolerance=0.05,
    singly=True,
):
    bb = _load_reference_bitbirch()
    bb.set_merge(merge_criterion, tolerance=tolerance)
    model = bb.BitBirch(threshold=threshold, branching_factor=branching_factor)
    model.fit(np.asarray(bits, dtype=np.uint8), singly=singly)
    clusters = tuple(
        tuple(int(member) for member in cluster)
        for cluster in model.get_cluster_mol_ids()
    )
    labels = model.get_assignments(len(bits)) - 1
    centroid_map = model.get_centroids_mol_ids()
    centroid_by_cluster = {
        tuple(int(member) for member in members): np.asarray(centroid, dtype=np.uint8)
        for centroid, members in zip(
            centroid_map["centroids"],
            centroid_map["mol_ids"],
            strict=True,
        )
    }
    centroids = np.asarray([centroid_by_cluster[cluster] for cluster in clusters])
    return labels.astype(np.intp), clusters, centroids


def _random_bits(rows, cols, seed=13, density=0.125):
    rng = np.random.default_rng(seed)
    bits = (rng.random((rows, cols)) < density).astype(np.uint8)
    empty = np.flatnonzero(bits.sum(axis=1) == 0)
    for row in empty:
        bits[row, row % cols] = 1
    return bits


def _batch_from_bits(bits):
    arr = np.asarray(bits, dtype=np.uint8)
    fingerprints = []
    for row in arr:
        on_bits = np.flatnonzero(row).astype(int).tolist()
        fingerprints.append(oefp.OEFP.from_on_bits(arr.shape[1], on_bits))
    return oefp.OEFPBatch.from_fingerprints(fingerprints)


def _tanimoto_distance_matrix(bits):
    """Build a complete SymmetricDistanceMatrix of Tanimoto distances (1 - sim).

    Uses a vectorized intersection matmul; only the upper-triangle Set() calls
    iterate in Python. Intended for modest n (thousands of rows, not millions).
    """
    from oecluster import DenseStorage, SymmetricDistanceMatrix

    arr = np.asarray(bits, dtype=np.float64)
    n = arr.shape[0]
    popcounts = arr.sum(axis=1)
    inter = arr @ arr.T                                  # (n, n) intersection
    union = popcounts[:, None] + popcounts[None, :] - inter
    with np.errstate(invalid="ignore", divide="ignore"):
        sim = np.where(union == 0.0, 1.0, inter / union)
    dist = 1.0 - sim
    storage = DenseStorage(n)
    for i in range(n):
        for j in range(i + 1, n):
            storage.Set(i, j, float(dist[i, j]))
    labels = [str(k) for k in range(n)]
    return SymmetricDistanceMatrix(storage, "tanimoto", labels, {})


def _assert_non_degenerate(report):
    """Quality fixtures must yield defined metrics (multi-cluster, non-singleton)."""
    assert report.num_clusters >= 2
    assert math.isfinite(report.silhouette)
    assert math.isfinite(report.mean_intra_distance)


def _centroid_bits(result):
    words = np.asarray(result.centroids.words, dtype=np.uint64)
    out = np.zeros((words.shape[0], result.centroids.num_bits), dtype=np.uint8)
    for row_idx, row_words in enumerate(words):
        for bit in range(result.centroids.num_bits):
            if int(row_words[bit // 64]) & (1 << (bit % 64)):
                out[row_idx, bit] = 1
    return out


@pytest.mark.parametrize("merge_criterion", ["radius", "diameter", "tolerance"])
def test_bitbirch_matches_reference_duplicate_blocks(merge_criterion):
    bits = np.array(
        [
            [1, 1, 0, 0],
            [1, 1, 0, 0],
            [0, 0, 1, 1],
            [0, 0, 1, 1],
        ],
        dtype=np.uint8,
    )
    expected_labels, expected_clusters, expected_centroids = _reference_result(
        bits,
        threshold=0.75,
        branching_factor=2,
        merge_criterion=merge_criterion,
    )

    observed = oecluster.bitbirch(
        _batch_from_bits(bits),
        threshold=0.75,
        branching_factor=2,
        merge_criterion=merge_criterion,
    )

    assert isinstance(observed, oecluster.BitBirchResult)
    assert isinstance(observed, oecluster.ClusteringResult)
    assert observed.labels.tolist() == expected_labels.tolist()
    assert observed.clusters == expected_clusters
    assert observed.cluster_sizes == tuple(len(cluster) for cluster in expected_clusters)
    np.testing.assert_array_equal(_centroid_bits(observed), expected_centroids)


def test_bitbirch_matches_reference_split_tie_ordering():
    bits = np.array(
        [
            [1, 1, 0, 0],
            [1, 1, 0, 0],
            [0, 0, 1, 1],
            [0, 0, 1, 1],
            [1, 0, 1, 0],
            [1, 0, 1, 0],
        ],
        dtype=np.uint8,
    )
    expected_labels, expected_clusters, expected_centroids = _reference_result(
        bits,
        threshold=0.75,
        branching_factor=2,
        merge_criterion="diameter",
    )

    observed = oecluster.bitbirch(
        _batch_from_bits(bits),
        threshold=0.75,
        branching_factor=2,
        merge_criterion="diameter",
    )

    assert observed.labels.tolist() == expected_labels.tolist()
    assert observed.clusters == expected_clusters
    np.testing.assert_array_equal(_centroid_bits(observed), expected_centroids)


def test_bitbirch_rejects_invalid_inputs():
    batch = _batch_from_bits(np.array([[1, 0]], dtype=np.uint8))

    with pytest.raises(ValueError, match="threshold"):
        oecluster.bitbirch(batch, threshold=-0.1)
    with pytest.raises(ValueError, match="branching_factor"):
        oecluster.bitbirch(batch, branching_factor=0)
    with pytest.raises(ValueError, match="merge_criterion"):
        oecluster.bitbirch(batch, merge_criterion="ward")
    with pytest.raises(ValueError, match="mode"):
        oecluster.bitbirch(batch, mode="approximate")
    with pytest.raises(TypeError, match="OEFPBatch"):
        oecluster.bitbirch(object())


def _reference_recluster(
    bits,
    *,
    initial_threshold=0.65,
    second_threshold=0.7,
    second_tolerance=0.0,
):
    bb = _load_reference_bitbirch()
    arr = np.asarray(bits, dtype=np.uint8)
    bb.set_merge("diameter")
    first = bb.BitBirch(branching_factor=50, threshold=initial_threshold)
    first.fit(arr, singly=True)
    data, bigs = first.prepare_data_BFs(arr)
    bb.set_merge("tolerance", tolerance=second_tolerance)
    second = bb.BitBirch(branching_factor=50, threshold=second_threshold)
    second.fit_BFs(data)
    second.fit_BFs(bigs)
    clusters = tuple(
        tuple(int(member) for member in cluster)
        for cluster in second.get_cluster_mol_ids()
    )
    labels = second.get_assignments(len(arr)) - 1
    return labels.astype(np.intp), clusters


def test_bitbirch_recluster_matches_reference():
    bits = np.array(
        [
            [1, 1, 0, 0],
            [1, 1, 0, 0],
            [0, 0, 1, 1],
            [0, 0, 1, 1],
            [1, 0, 1, 0],
            [1, 0, 1, 0],
        ],
        dtype=np.uint8,
    )
    expected_labels, expected_clusters = _reference_recluster(
        bits,
        initial_threshold=0.65,
        second_threshold=0.75,
        second_tolerance=0.0,
    )

    observed = oecluster.bitbirch_recluster(
        _batch_from_bits(bits),
        initial_threshold=0.65,
        second_threshold=0.75,
        second_tolerance=0.0,
    )

    assert observed.labels.tolist() == expected_labels.tolist()
    assert observed.clusters == expected_clusters


def test_bitbirch_fast_mode_matches_strict_for_single_partition():
    bits = np.array(
        [
            [1, 1, 0, 0],
            [1, 1, 0, 0],
            [0, 0, 1, 1],
            [1, 0, 1, 0],
        ],
        dtype=np.uint8,
    )
    batch = _batch_from_bits(bits)

    strict = oecluster.bitbirch(batch, threshold=0.65, mode="strict_parity")
    fast = oecluster.bitbirch(batch, threshold=0.65, mode="fast")

    assert fast.labels.tolist() == strict.labels.tolist()
    assert fast.clusters == strict.clusters


def test_bitbirch_recluster_fast_mode_matches_strict_for_single_partition():
    bits = np.array(
        [
            [1, 1, 0, 0],
            [1, 1, 0, 0],
            [0, 0, 1, 1],
            [0, 0, 1, 1],
            [1, 0, 1, 0],
            [1, 0, 1, 0],
        ],
        dtype=np.uint8,
    )
    batch = _batch_from_bits(bits)

    strict = oecluster.bitbirch_recluster(
        batch,
        initial_threshold=0.65,
        second_threshold=0.75,
        second_tolerance=0.0,
        mode="strict_parity",
    )
    fast = oecluster.bitbirch_recluster(
        batch,
        initial_threshold=0.65,
        second_threshold=0.75,
        second_tolerance=0.0,
        mode="fast",
    )

    assert fast.labels.tolist() == strict.labels.tolist()
    assert fast.clusters == strict.clusters


def test_bitbirch_refine_fast_mode_matches_strict():
    bits = _prune_reassign_fixture_bits()
    batch = _batch_from_bits(bits)

    strict = oecluster.bitbirch_refine(
        batch,
        threshold=0.3,
        branching_factor=2,
        merge_criterion="diameter",
        singly=False,
        redistribute_largest_cluster=True,
        reassign_top_clusters=2,
        mode="strict_parity",
    )
    fast = oecluster.bitbirch_refine(
        batch,
        threshold=0.3,
        branching_factor=2,
        merge_criterion="diameter",
        singly=False,
        redistribute_largest_cluster=True,
        reassign_top_clusters=2,
        mode="fast",
    )

    assert fast.labels.tolist() == strict.labels.tolist()
    assert fast.clusters == strict.clusters


def test_bitbirch_fast_mode_docstrings_mark_reserved_alias():
    # bitbirch / bitbirch_recluster now document real fast behavior.
    fast_text = "partitions the input, fits chunk trees in parallel"
    assert fast_text in " ".join((oecluster.bitbirch.__doc__ or "").split())
    assert fast_text in " ".join((oecluster.bitbirch_recluster.__doc__ or "").split())
    # refine still runs the strict path; its docstring keeps the reserved note.
    refine_text = "refinement always runs the strict-parity path"
    assert refine_text in " ".join((oecluster.bitbirch_refine.__doc__ or "").split())


def _reference_refine_reassign(bits, *, top=2):
    bb = _load_reference_bitbirch()
    arr = np.asarray(bits, dtype=np.uint8)
    bb.set_merge("diameter")
    model = bb.BitBirch(branching_factor=3, threshold=0.65)
    model.fit(arr, singly=False)
    model.reassign(arr, top=top)
    clusters = tuple(
        tuple(int(member) for member in cluster)
        for cluster in model.get_cluster_mol_ids()
    )
    labels = model.get_assignments(len(arr)) - 1
    return labels.astype(np.intp), clusters


def test_bitbirch_refine_reassign_top_clusters_matches_reference():
    bits = np.array(
        [
            [0, 1, 1, 0],
            [1, 1, 0, 0],
            [1, 1, 1, 0],
            [0, 0, 0, 1],
            [1, 0, 1, 0],
            [1, 0, 1, 1],
        ],
        dtype=np.uint8,
    )
    expected_labels, expected_clusters = _reference_refine_reassign(bits, top=2)

    observed = oecluster.bitbirch_refine(
        _batch_from_bits(bits),
        threshold=0.65,
        branching_factor=3,
        merge_criterion="diameter",
        singly=False,
        reassign_top_clusters=2,
    )

    assert observed.labels.tolist() == expected_labels.tolist()
    assert observed.clusters == expected_clusters


def _reference_refine_prune(
    bits,
    *,
    threshold,
    branching_factor,
    reassign_top_clusters=0,
):
    bb = _load_reference_bitbirch()
    arr = np.asarray(bits, dtype=np.uint8)
    bb.set_merge("diameter")
    model = bb.BitBirch(branching_factor=branching_factor, threshold=threshold)
    model.fit(arr, singly=False)
    model.prune(arr)
    if reassign_top_clusters:
        model.reassign(arr, top=reassign_top_clusters)
    clusters = tuple(
        tuple(int(member) for member in cluster)
        for cluster in model.get_cluster_mol_ids()
    )
    labels = model.get_assignments(len(arr)) - 1
    return labels.astype(np.intp), clusters


def _prune_redistribution_fixture_bits():
    # This fixture makes prune visibly split member 0 out of the largest cluster.
    return np.array(
        [
            [0, 0, 1, 0, 0, 1],
            [1, 0, 0, 0, 0, 0],
            [0, 1, 1, 0, 0, 1],
            [0, 1, 0, 0, 0, 1],
            [0, 0, 0, 0, 0, 1],
            [1, 0, 0, 0, 0, 0],
            [0, 1, 0, 0, 0, 0],
            [0, 1, 0, 0, 0, 1],
        ],
        dtype=np.uint8,
    )


def _prune_reassign_fixture_bits():
    # This fixture makes reassign(top=2) visibly move member 0 after prune.
    return np.array(
        [
            [1, 1, 0, 0, 0, 0],
            [0, 1, 0, 0, 0, 0],
            [1, 0, 0, 0, 0, 0],
            [1, 0, 1, 1, 0, 0],
            [0, 0, 0, 0, 0, 1],
            [1, 0, 1, 0, 0, 0],
            [1, 1, 1, 0, 0, 0],
            [0, 1, 0, 0, 0, 0],
        ],
        dtype=np.uint8,
    )


def _duplicate_block_bits(*, rows=320, cols=128, density=0.2, prototypes=4):
    rng = np.random.default_rng(20260520)
    prototype_bits = (rng.random((prototypes, cols)) < density).astype(np.uint8)
    empty_rows = np.flatnonzero(prototype_bits.sum(axis=1) == 0)
    for row in empty_rows:
        prototype_bits[row, row % cols] = 1
    return prototype_bits[np.arange(rows) % prototypes].copy()


def test_bitbirch_refine_reassign_threads_match_serial_on_duplicate_blocks():
    bits = _duplicate_block_bits()
    batch = _batch_from_bits(bits)
    kwargs = {
        "threshold": 0.65,
        "branching_factor": 2,
        "merge_criterion": "diameter",
        "singly": False,
        "reassign_top_clusters": 4,
    }

    serial = oecluster.bitbirch_refine(batch, num_threads=1, **kwargs)
    threaded = oecluster.bitbirch_refine(batch, num_threads=2, **kwargs)

    assert threaded.labels.tolist() == serial.labels.tolist()
    assert threaded.clusters == serial.clusters
    assert threaded.cluster_sizes == serial.cluster_sizes


class _IndexOnly:
    """An integer through operator.index, and an error through int()."""

    def __index__(self):
        return 1

    def __int__(self):
        raise AssertionError("int() reached a value meant for operator.index")


_BITBIRCH_CALLS = [
    lambda batch, **kw: oecluster.bitbirch(batch, **kw),
    lambda batch, **kw: oecluster.bitbirch_recluster(
        batch, initial_threshold=0.5, second_threshold=0.5, **kw),
    lambda batch, **kw: oecluster.bitbirch_refine(batch, **kw),
]


@pytest.mark.parametrize("call", _BITBIRCH_CALLS)
@pytest.mark.parametrize("value", [2.5, "8", None], ids=["float", "str", "none"])
def test_bitbirch_refuses_a_non_integer_num_threads(call, value):
    batch = _batch_from_bits(_random_bits(30, 64, seed=5))
    with pytest.raises(TypeError):
        call(batch, num_threads=value)


@pytest.mark.parametrize("call", _BITBIRCH_CALLS)
def test_bitbirch_accepts_a_numpy_integer_or_index_only_num_threads(call):
    """_IndexOnly reaches every coercion, so bitbirch_refine's second one,
    which a refused value never reaches, cannot keep int() unnoticed."""
    batch = _batch_from_bits(_random_bits(30, 64, seed=5))
    expected = call(batch, num_threads=2).labels
    assert list(call(batch, num_threads=np.int64(2)).labels) == list(expected)
    call(batch, num_threads=_IndexOnly())


def _duplicate_prototype_bits(seed, *, prototypes=4, copies=3, cols=24, density=0.2):
    # Repeated prototypes give reassign equal centroids to select, and a tie
    # resolves to the lower centroid index, so the later of two equal centroids
    # can win no members at all. That is one route to a starved centroid, not
    # the only one: the fixture below starves one with no duplicate anywhere.
    rng = np.random.default_rng(seed)
    prototype_bits = (rng.random((prototypes, cols)) < density).astype(np.uint8)
    for row in np.flatnonzero(prototype_bits.sum(axis=1) == 0):
        prototype_bits[row, row % cols] = 1
    return np.repeat(prototype_bits, copies, axis=0)


def test_bitbirch_refine_reassign_drops_a_starved_cluster_instead_of_duplicating():
    bits = _duplicate_prototype_bits(0)

    observed = oecluster.bitbirch_refine(
        _batch_from_bits(bits),
        threshold=0.3,
        branching_factor=3,
        merge_criterion="diameter",
        singly=False,
        reassign_top_clusters=5,
    )

    # Five centroids are selected here and only four win members. Before the
    # fix the fifth subcluster kept the members it had already handed to
    # another, and the result carried sample 9 in two clusters.
    assert observed.clusters == ((3, 4, 5), (10, 11, 9), (0, 1, 2), (6, 7, 8))
    assert observed.labels.tolist() == [2, 2, 2, 0, 0, 0, 3, 3, 3, 1, 1, 1]
    assert observed.cluster_sizes == (3, 3, 3, 3)


def test_bitbirch_refine_reassign_starves_a_centroid_with_distinct_fingerprints():
    # A starved centroid needs no duplicate fingerprint and no pair of
    # near-identical centroids. All 24 rows here are distinct, the greatest
    # Tanimoto between any two of the five selected centroids is 0.0909, and
    # the ten-member cluster's majority centroid is all-zero, so it scores
    # zero against every row and is the nearest centroid of none of them.
    rng = np.random.default_rng(2)
    bits = (rng.random((24, 32)) < 0.18).astype(np.uint8)
    assert len(np.unique(bits, axis=0)) == len(bits)
    kwargs = {
        "threshold": 0.1,
        "branching_factor": 2,
        "merge_criterion": "diameter",
        "singly": False,
    }

    fitted = oecluster.bitbirch(_batch_from_bits(bits), **kwargs)
    observed = oecluster.bitbirch_refine(
        _batch_from_bits(bits), reassign_top_clusters=5, **kwargs)

    # Reassignment collapses the six fitted clusters to five. The starved
    # centroid is the largest selected subcluster's, but before the fix the
    # stale slot was the smallest: the old writer filled the first four
    # selected slots and left the fifth holding (4, 12), which the
    # reassignment had already placed in another cluster, so the result
    # reported six.
    assert len(fitted.clusters) == 6
    assert observed.clusters == (
        (14, 15, 19, 21, 22, 4, 12),
        (9, 20, 2, 3, 5, 6, 8),
        (7, 13, 10, 11, 17),
        (18, 0, 1, 16),
        (23,),
    )
    assert sorted(
        member for cluster in observed.clusters for member in cluster
    ) == list(range(len(bits)))


@pytest.mark.parametrize("reassign_top_clusters", [3, 4, 5])
@pytest.mark.parametrize("branching_factor", [2, 3])
@pytest.mark.parametrize("threshold", [0.3, 0.5, 0.7])
def test_bitbirch_refine_reassign_partitions_duplicate_bearing_inputs(
        threshold, branching_factor, reassign_top_clusters):
    """Every refined result assigns each input index to exactly one cluster."""
    for seed in range(25):
        bits = _duplicate_prototype_bits(seed)
        observed = oecluster.bitbirch_refine(
            _batch_from_bits(bits),
            threshold=threshold,
            branching_factor=branching_factor,
            merge_criterion="diameter",
            singly=False,
            reassign_top_clusters=reassign_top_clusters,
        )

        assigned = [member for cluster in observed.clusters for member in cluster]
        assert sorted(assigned) == list(range(len(bits))), (
            f"seed={seed} threshold={threshold} "
            f"branching_factor={branching_factor} "
            f"reassign_top_clusters={reassign_top_clusters} "
            f"clusters={observed.clusters}")
        assert observed.cluster_sizes == tuple(
            len(cluster) for cluster in observed.clusters)


def test_bitbirch_refine_redistribute_largest_cluster_requires_parent_pointers():
    bits = _prune_redistribution_fixture_bits()

    with pytest.raises(ValueError, match="singly=False"):
        oecluster.bitbirch_refine(
            _batch_from_bits(bits),
            threshold=0.4,
            branching_factor=2,
            merge_criterion="diameter",
            singly=True,
            redistribute_largest_cluster=True,
        )


def test_bitbirch_refine_redistribute_largest_cluster_matches_reference():
    bits = _prune_redistribution_fixture_bits()
    expected_labels, expected_clusters = _reference_refine_prune(
        bits,
        threshold=0.4,
        branching_factor=2,
    )

    observed = oecluster.bitbirch_refine(
        _batch_from_bits(bits),
        threshold=0.4,
        branching_factor=2,
        merge_criterion="diameter",
        singly=False,
        redistribute_largest_cluster=True,
    )

    assert observed.labels.tolist() == expected_labels.tolist()
    assert observed.clusters == expected_clusters
    assert observed.cluster_sizes == tuple(len(cluster) for cluster in expected_clusters)


def _emptied_subcluster_fixture_bits():
    # This fixture makes prune leave one leaf subcluster with no members: the
    # reference reports it as an empty cluster, and the native result drops it.
    return np.array(
        [
            [1, 1, 1],
            [0, 1, 1],
            [1, 0, 0],
            [1, 1, 0],
        ],
        dtype=np.uint8,
    )


def _refine_emptied_subcluster(bits):
    return oecluster.bitbirch_refine(
        _batch_from_bits(bits),
        threshold=0.8,
        branching_factor=2,
        merge_criterion="diameter",
        singly=False,
        redistribute_largest_cluster=True,
    )


def test_bitbirch_refine_reassign_leaves_an_already_empty_selected_slot_alone():
    # The boundary of the duplication fix. Prune empties a leaf subcluster,
    # which sorts last among the five selected, and reassignment starves one
    # centroid. The old writer left only the last selected slot unwritten,
    # which was that empty one, and the new writer empties the starved
    # centroid's own slot, which here is the same one. Neither leaves a slot
    # holding members, so this case never duplicated a member and does not
    # move now.
    bits = _emptied_subcluster_fixture_bits()
    pruned = _refine_emptied_subcluster(bits)

    reassigned = oecluster.bitbirch_refine(
        _batch_from_bits(bits),
        threshold=0.8,
        branching_factor=2,
        merge_criterion="diameter",
        singly=False,
        redistribute_largest_cluster=True,
        reassign_top_clusters=5,
    )

    assert reassigned.clusters == pruned.clusters
    assert reassigned.labels.tolist() == pruned.labels.tolist()
    assert reassigned.cluster_sizes == pruned.cluster_sizes


def test_bitbirch_refine_prune_matches_the_reference_but_drops_the_emptied_subcluster():
    bits = _emptied_subcluster_fixture_bits()
    expected_labels, expected_clusters = _reference_refine_prune(
        bits,
        threshold=0.8,
        branching_factor=2,
    )
    # A deliberate divergence from the reference: it reports the subcluster
    # prune emptied as an empty cluster, and the native result drops it. The
    # sort is size-descending, so the empties are trailing and no label moves
    # -- which is why the labels below are still compared as-is.
    assert expected_clusters[-1] == ()
    expected_clusters = tuple(
        cluster for cluster in expected_clusters if cluster)

    observed = _refine_emptied_subcluster(bits)

    assert observed.labels.tolist() == expected_labels.tolist()
    assert observed.clusters == expected_clusters


def test_cluster_report_accepts_a_refinement_that_emptied_a_subcluster():
    bits = _emptied_subcluster_fixture_bits()
    result = _refine_emptied_subcluster(bits)

    assert all(len(cluster) > 0 for cluster in result.clusters)
    report = oecluster.cluster_report(
        result, _tanimoto_distance_matrix(bits), allow_nonmetric=True)
    assert report.num_clusters == len(result.clusters)


def test_isim_report_accepts_a_refinement_that_emptied_a_subcluster():
    bits = _emptied_subcluster_fixture_bits()
    result = _refine_emptied_subcluster(bits)

    assert all(len(cluster) > 0 for cluster in result.clusters)
    report = oecluster.isim_report(result, _batch_from_bits(bits))
    assert report.num_clusters == len(result.clusters)


def test_bitbirch_refine_prune_then_reassign_matches_reference():
    bits = _prune_reassign_fixture_bits()
    expected_labels, expected_clusters = _reference_refine_prune(
        bits,
        threshold=0.3,
        branching_factor=2,
        reassign_top_clusters=2,
    )

    observed = oecluster.bitbirch_refine(
        _batch_from_bits(bits),
        threshold=0.3,
        branching_factor=2,
        merge_criterion="diameter",
        singly=False,
        redistribute_largest_cluster=True,
        reassign_top_clusters=2,
    )

    assert observed.labels.tolist() == expected_labels.tolist()
    assert observed.clusters == expected_clusters


def test_bitbirch_fast_small_n_matches_strict():
    bits = np.array(
        [[1, 1, 0, 0], [1, 1, 0, 0], [0, 0, 1, 1], [0, 0, 1, 1]],
        dtype=np.uint8,
    )
    batch = _batch_from_bits(bits)
    strict = oecluster.bitbirch(batch, threshold=0.75, branching_factor=2,
                                mode="strict_parity")
    fast = oecluster.bitbirch(batch, threshold=0.75, branching_factor=2,
                              mode="fast")
    assert fast.labels.tolist() == strict.labels.tolist()
    assert fast.clusters == strict.clusters
    assert fast.cluster_sizes == strict.cluster_sizes


@pytest.mark.parametrize("num_threads", [1, 2, 4, 0])
def test_bitbirch_fast_is_deterministic_across_threads(num_threads):
    bits = _random_bits(DETERMINISM_N, 64)
    batch = _batch_from_bits(bits)
    reference = oecluster.bitbirch(batch, threshold=0.5, branching_factor=50,
                                   mode="fast", num_threads=1)
    observed = oecluster.bitbirch(batch, threshold=0.5, branching_factor=50,
                                  mode="fast", num_threads=num_threads)
    assert observed.labels.tolist() == reference.labels.tolist()
    assert observed.clusters == reference.clusters


@pytest.mark.parametrize("merge_criterion", ["diameter", "tolerance"])
@pytest.mark.parametrize("singly", [True, False])
def test_bitbirch_fast_quality_equivalent_to_strict(merge_criterion, singly):
    bits = _random_bits(QUALITY_N, 64, seed=7)
    batch = _batch_from_bits(bits)
    dm = _tanimoto_distance_matrix(bits)

    common = {"threshold": 0.5, "branching_factor": 50,
              "merge_criterion": merge_criterion, "singly": singly}
    strict = oecluster.bitbirch(batch, mode="strict_parity", **common)
    fast = oecluster.bitbirch(batch, mode="fast", **common)

    strict_report = oecluster.cluster_report(strict, dm)
    fast_report = oecluster.cluster_report(fast, dm)
    _assert_non_degenerate(strict_report)
    _assert_non_degenerate(fast_report)

    assert abs(fast_report.num_clusters - strict_report.num_clusters) <= \
        FAST_QUALITY_TOL["num_clusters_rel"] * strict_report.num_clusters
    assert fast_report.mean_intra_distance <= \
        strict_report.mean_intra_distance + FAST_QUALITY_TOL["intra_distance_abs"]
    assert fast_report.median_intra_distance <= \
        strict_report.median_intra_distance + FAST_QUALITY_TOL["intra_distance_abs"]
    assert fast_report.silhouette >= \
        strict_report.silhouette - FAST_QUALITY_TOL["silhouette_abs"]
    for i in range(len(strict_report.coverage_at)):
        assert fast_report.coverage_at[i] >= \
            strict_report.coverage_at[i] - FAST_QUALITY_TOL["coverage_abs"]


@pytest.mark.parametrize("num_threads", [1, 2, 4, 0])
def test_bitbirch_recluster_fast_deterministic_across_threads(num_threads):
    bits = _random_bits(DETERMINISM_N, 64, seed=21)
    batch = _batch_from_bits(bits)
    ref = oecluster.bitbirch_recluster(batch, initial_threshold=0.5,
                                       second_threshold=0.5, branching_factor=50,
                                       mode="fast", num_threads=1)
    obs = oecluster.bitbirch_recluster(batch, initial_threshold=0.5,
                                       second_threshold=0.5, branching_factor=50,
                                       mode="fast", num_threads=num_threads)
    assert obs.labels.tolist() == ref.labels.tolist()
    assert obs.clusters == ref.clusters


def test_bitbirch_recluster_fast_quality_equivalent_to_strict():
    bits = _random_bits(QUALITY_N, 64, seed=21)
    batch = _batch_from_bits(bits)
    dm = _tanimoto_distance_matrix(bits)
    common = {"initial_threshold": 0.5, "second_threshold": 0.5, "branching_factor": 50}
    strict = oecluster.bitbirch_recluster(batch, mode="strict_parity", **common)
    fast = oecluster.bitbirch_recluster(batch, mode="fast", **common)
    strict_report = oecluster.cluster_report(strict, dm)
    fast_report = oecluster.cluster_report(fast, dm)
    _assert_non_degenerate(strict_report)
    _assert_non_degenerate(fast_report)
    assert abs(fast_report.num_clusters - strict_report.num_clusters) <= \
        FAST_QUALITY_TOL["num_clusters_rel"] * strict_report.num_clusters
    assert fast_report.mean_intra_distance <= \
        strict_report.mean_intra_distance + FAST_QUALITY_TOL["intra_distance_abs"]
    assert fast_report.silhouette >= \
        strict_report.silhouette - FAST_QUALITY_TOL["silhouette_abs"]


def test_bitbirch_refine_fast_equals_strict():
    bits = _random_bits(5000, 64, seed=33)
    batch = _batch_from_bits(bits)
    common = {"threshold": 0.5, "branching_factor": 50, "singly": False,
              "reassign_top_clusters": 2}
    strict = oecluster.bitbirch_refine(batch, mode="strict_parity", **common)
    fast = oecluster.bitbirch_refine(batch, mode="fast", **common)
    assert fast.labels.tolist() == strict.labels.tolist()
    assert fast.clusters == strict.clusters
