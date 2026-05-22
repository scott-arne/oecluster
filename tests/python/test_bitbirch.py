"""Parity tests for native BitBirch clustering."""

from __future__ import annotations

import importlib
import sys
from pathlib import Path

import numpy as np
import pytest

import oecluster
import oefp


BITBIRCH_REPO = Path("/Users/johnss51/Development/python/bitbirch")
ISIM_REPO = Path("/Users/johnss51/Development/python/iSIM")


def _load_reference_bitbirch():
    for repo in (str(ISIM_REPO), str(BITBIRCH_REPO)):
        if repo not in sys.path:
            sys.path.insert(0, repo)
    module = importlib.import_module("bitbirch.bitbirch")
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


def _batch_from_bits(bits):
    arr = np.asarray(bits, dtype=np.uint8)
    fingerprints = []
    for row in arr:
        on_bits = np.flatnonzero(row).astype(int).tolist()
        fingerprints.append(oefp.OEFP.from_on_bits(arr.shape[1], on_bits))
    return oefp.OEFPBatch.from_fingerprints(fingerprints)


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


def test_bitbirch_fast_mode_matches_strict_until_optimized():
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


def test_bitbirch_recluster_fast_mode_matches_strict_until_optimized():
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


def test_bitbirch_refine_fast_mode_matches_strict_until_optimized():
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


def test_bitbirch_refine_prune_matches_zero_sample_reference_centroid():
    bits = np.array(
        [
            [1, 1, 1],
            [0, 1, 1],
            [1, 0, 0],
            [1, 1, 0],
        ],
        dtype=np.uint8,
    )
    expected_labels, expected_clusters = _reference_refine_prune(
        bits,
        threshold=0.8,
        branching_factor=2,
    )

    observed = oecluster.bitbirch_refine(
        _batch_from_bits(bits),
        threshold=0.8,
        branching_factor=2,
        merge_criterion="diameter",
        singly=False,
        redistribute_largest_cluster=True,
    )

    assert observed.labels.tolist() == expected_labels.tolist()
    assert observed.clusters == expected_clusters


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
