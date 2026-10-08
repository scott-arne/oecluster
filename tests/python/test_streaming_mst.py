"""HDBSCAN and single-linkage agglomerative clustering from a comparison.

Exactness is defined against a matrix filled through the same comparison's
Compare(min, max), so every parity test builds that matrix rather than
calling pdist, whose batched kernels may differ in the last bit.
"""

from typing import Any

import oecluster
import pytest
from openeye import oechem

FP_SMILES = ["CCO", "CCCO", "CCCCO", "c1ccccc1", "Cc1ccccc1", "CCc1ccccc1",
             "CC(=O)O", "CC(=O)OC", "CCN", "CCCN", "C1CCCCC1", "c1ccncc1",
             "CCCCCO", "c1ccc2ccccc2c1", "CC(C)O", "OCCO", "CCOC(=O)C",
             "c1ccc(O)cc1", "c1ccc(N)cc1", "Clc1ccccc1"]

DESCRIPTOR_SMILES = ["CCO", "CCC", "CCCC", "c1ccccc1", "CCN", "CCOC",
                     "CC(=O)O", "CCCCCC"]


def _mols(smiles_list):
    mols = []
    for idx, smi in enumerate(smiles_list):
        mol = oechem.OEGraphMol()
        oechem.OESmilesToMol(mol, smi)
        mol.SetTitle(f"mol{idx}")
        mols.append(mol)
    return mols


def _compare_matrix(comparison):
    """The matrix a comparison form must reproduce: every pair via Compare."""
    n = comparison.Size()
    storage = oecluster.DenseStorage(n)
    for i in range(n):
        for j in range(i + 1, n):
            storage.Set(i, j, comparison.Compare(i, j))
    return oecluster.SymmetricDistanceMatrix(
        storage, "test", [f"item_{i}" for i in range(n)], {})


def _hdbscan_fields(result):
    return (result.labels.tolist(), result.clusters,
            result.probabilities.tolist())


def _agglomerative_fields(result):
    return (result.labels.tolist(), result.clusters, list(result.children),
            list(result.distances), list(result.cluster_sizes))


@pytest.mark.parametrize("chunk_size", [0, 1, 64])
@pytest.mark.parametrize("num_threads", [1, 4])
def test_hdbscan_three_inputs_match_a_compare_filled_matrix(num_threads,
                                                            chunk_size):
    mols = _mols(FP_SMILES)
    prebuilt = oecluster.FingerprintComparison(mols)
    matrix = _compare_matrix(prebuilt)
    for min_cluster_size, min_samples in ((2, 1), (2, None), (3, 4)):
        options: dict[str, Any] = {"min_cluster_size": min_cluster_size,
                                   "min_samples": min_samples}
        expected = _hdbscan_fields(oecluster.hdbscan(matrix, **options))
        for items, extra in ((prebuilt, {}),
                             (mols, {"comparison": "fingerprint"})):
            result = oecluster.hdbscan(items, num_threads=num_threads,
                                       chunk_size=chunk_size, **options,
                                       **extra)
            assert _hdbscan_fields(result) == expected


@pytest.mark.parametrize("num_threads", [1, 4])
def test_single_linkage_three_inputs_match_a_compare_filled_matrix(num_threads):
    mols = _mols(FP_SMILES)
    prebuilt = oecluster.FingerprintComparison(mols)
    matrix = _compare_matrix(prebuilt)
    cases: list[dict[str, Any]] = [
        {"n_clusters": 4}, {"n_clusters": 1},
        {"distance_threshold": 0.5},
        {"n_clusters": 3, "compute_full_tree": False}]
    for options in cases:
        expected = _agglomerative_fields(
            oecluster.agglomerative(matrix, linkage="single", **options))
        for items, extra in ((prebuilt, {}),
                             (mols, {"comparison": "fingerprint"})):
            result = oecluster.agglomerative(
                items, linkage="single", num_threads=num_threads, **options,
                **extra)
            assert _agglomerative_fields(result) == expected


def test_descriptors_match_a_compare_filled_matrix():
    mols = _mols(DESCRIPTOR_SMILES)
    prebuilt = oecluster.DescriptorComparison(mols, metric="euclidean")
    matrix = _compare_matrix(prebuilt)
    assert (_hdbscan_fields(oecluster.hdbscan(prebuilt, min_cluster_size=2))
            == _hdbscan_fields(oecluster.hdbscan(matrix, min_cluster_size=2)))
    assert (_agglomerative_fields(oecluster.agglomerative(
                prebuilt, linkage="single", n_clusters=3))
            == _agglomerative_fields(oecluster.agglomerative(
                matrix, linkage="single", n_clusters=3)))


def test_the_alias_still_takes_a_matrix():
    matrix = oecluster.pdist(_mols(FP_SMILES), "fingerprint")
    assert (_hdbscan_fields(oecluster.hdbscan(distance_matrix=matrix,
                                              min_cluster_size=2))
            == _hdbscan_fields(oecluster.hdbscan(matrix, min_cluster_size=2)))
    assert (_agglomerative_fields(oecluster.agglomerative(
                distance_matrix=matrix, n_clusters=3))
            == _agglomerative_fields(oecluster.agglomerative(
                matrix, n_clusters=3)))


@pytest.mark.parametrize("entry", ["hdbscan", "agglomerative"])
def test_the_input_is_required_and_the_alias_takes_only_a_matrix(entry):
    function = getattr(oecluster, entry)
    mols = _mols(FP_SMILES[:4])
    matrix = oecluster.pdist(mols, "fingerprint")
    with pytest.raises(TypeError, match="got both items and distance_matrix"):
        function(matrix, distance_matrix=matrix)
    with pytest.raises(TypeError, match="missing required argument: 'items'"):
        function()
    with pytest.raises(TypeError, match="expects a SymmetricDistanceMatrix"):
        function(distance_matrix=mols)
    with pytest.raises(TypeError, match="a sequence of items with comparison="):
        function(mols)


def test_there_is_one_sentinel_shared_by_every_entry_point():
    import inspect

    sentinel = oecluster._MISSING
    for function in (oecluster.hdbscan, oecluster.agglomerative):
        parameters = inspect.signature(function).parameters
        for name in ("items", "distance_matrix"):
            assert parameters[name].default is sentinel, (function, name)


def test_hdbscan_chunk_size_defaults_to_64():
    import inspect

    parameters = inspect.signature(oecluster.hdbscan).parameters
    assert parameters["chunk_size"].default == 64
    assert oecluster.HDBSCANOptions().chunk_size == 64


@pytest.mark.parametrize("linkage", ["average", "complete", "weighted"])
def test_a_comparison_needs_single_linkage(linkage):
    mols = _mols(FP_SMILES[:6])
    for items, extra in ((oecluster.FingerprintComparison(mols), {}),
                         (mols, {"comparison": "fingerprint"})):
        with pytest.raises(ValueError,
                           match="clusters a comparison only with "
                                 "linkage='single'"):
            oecluster.agglomerative(items, linkage=linkage, n_clusters=2,
                                    **extra)
    # The default linkage is average, so a comparison must name single.
    with pytest.raises(ValueError, match="linkage='single'"):
        oecluster.agglomerative(mols, comparison="fingerprint", n_clusters=2)


def test_the_linkage_is_refused_before_the_comparison_is_built():
    # "nonexistent" names no comparison, so had the build run first the error
    # would be about the name, not about the linkage.
    with pytest.raises(ValueError, match="linkage='single'"):
        oecluster.agglomerative(_mols(FP_SMILES[:4]), comparison="nonexistent",
                                n_clusters=2)


@pytest.mark.parametrize("min_samples", [1, 3])
def test_a_nan_alpha_is_refused_on_every_path(min_samples):
    # "nonexistent" names no comparison, so the refusal provably comes before
    # any build, let alone any pair.
    mols = _mols(FP_SMILES[:6])
    matrix = oecluster.pdist(mols, "fingerprint")
    for items, extra in ((matrix, {}), (oecluster.FingerprintComparison(mols), {}),
                         (mols, {"comparison": "nonexistent"})):
        with pytest.raises(ValueError, match="HDBSCAN alpha must be positive"):
            oecluster.hdbscan(items, min_cluster_size=2, min_samples=min_samples,
                              alpha=float("nan"), **extra)


def test_the_item_count_bounds_hold_for_comparisons():
    mols = _mols(FP_SMILES[:6])
    with pytest.raises(ValueError,
                       match=r"min_samples must be at most the item count \(6\)"):
        oecluster.hdbscan(mols, comparison="fingerprint", min_samples=7)
    with pytest.raises(ValueError,
                       match=r"min_samples must be at most the item count \(6\)"):
        oecluster.hdbscan(oecluster.FingerprintComparison(mols),
                          min_cluster_size=7)
    with pytest.raises(ValueError,
                       match=r"n_clusters must be at most the item count \(6\)"):
        oecluster.agglomerative(mols, comparison="fingerprint",
                                linkage="single", n_clusters=7)
    # A threshold cut ignores n_clusters.
    result = oecluster.agglomerative(mols, comparison="fingerprint",
                                     linkage="single", n_clusters=7,
                                     distance_threshold=0.5)
    assert len(result.labels) == 6


def test_a_dropped_item_is_refused():
    mols = _mols(["O", *DESCRIPTOR_SMILES])
    with pytest.raises(ValueError,
                       match=r"dropped item 0 \(missing-descriptor\)"):
        oecluster.hdbscan(mols, comparison="descriptor", metric="euclidean",
                          min_cluster_size=2)
    with pytest.raises(ValueError,
                       match=r"dropped item 0 \(missing-descriptor\)"):
        oecluster.agglomerative(mols, comparison="descriptor",
                                metric="euclidean", linkage="single")


def test_a_nonmetric_comparison_needs_allow_nonmetric():
    mols = _mols(FP_SMILES)
    with pytest.raises(ValueError, match="violates the triangle inequality"):
        oecluster.hdbscan(mols, comparison="fingerprint", metric="dice",
                          min_cluster_size=2)
    with pytest.raises(ValueError, match="violates the triangle inequality"):
        oecluster.agglomerative(mols, comparison="fingerprint", metric="dice",
                                linkage="single")
    assert len(oecluster.hdbscan(mols, comparison="fingerprint", metric="dice",
                                 min_cluster_size=2,
                                 allow_nonmetric=True).labels) == len(mols)
    assert len(oecluster.agglomerative(
        mols, comparison="fingerprint", metric="dice", linkage="single",
        allow_nonmetric=True).labels) == len(mols)


def test_a_truthy_allow_nonmetric_is_refused_before_dispatch():
    mols = _mols(FP_SMILES[:4])
    with pytest.raises(TypeError, match="allow_nonmetric must be True or"):
        oecluster.hdbscan(mols, comparison="nonexistent", min_cluster_size=2,
                          allow_nonmetric="False")
    with pytest.raises(TypeError, match="allow_nonmetric must be True or"):
        oecluster.agglomerative(mols, comparison="nonexistent",
                                linkage="single", allow_nonmetric="False")


def test_a_local_option_is_coerced_before_the_comparison_is_built():
    # "nonexistent" names no comparison, so had the build run first the error
    # would name the comparison, not the option the caller got wrong.
    mols = _mols(FP_SMILES[:4])
    with pytest.raises(TypeError, match="allow_single_cluster must be a bool"):
        oecluster.hdbscan(mols, comparison="nonexistent", min_cluster_size=2,
                          allow_single_cluster="yes")
    with pytest.raises(TypeError, match="compute_full_tree must be a bool"):
        oecluster.agglomerative(mols, comparison="nonexistent",
                                linkage="single", compute_full_tree="yes")


def test_rocs_is_refused_before_it_is_built():
    pytest.importorskip("openeye.oeomega")
    from openeye import oeomega

    omega = oeomega.OEOmega()
    omega.SetMaxConfs(1)
    omega.SetStrictStereo(False)
    mols = []
    for smi in ("c1ccccc1", "Cc1ccccc1", "c1ccc(O)cc1"):
        mol = oechem.OEMol()
        oechem.OESmilesToMol(mol, smi)
        assert omega(mol)
        mols.append(mol)
    prebuilt = oecluster.ROCSComparison(mols)
    for items, extra in ((mols, {"comparison": "ROCS"}), (prebuilt, {})):
        with pytest.raises(ValueError,
                           match="cannot cluster a ROCS comparison without a "
                                 "matrix.*two passes that must agree"):
            oecluster.hdbscan(items, min_cluster_size=2, **extra)
        with pytest.raises(ValueError,
                           match="cannot cluster a ROCS comparison without a "
                                 "matrix.*order the pairs are scored"):
            oecluster.agglomerative(items, linkage="single", **extra)
    matrix = oecluster.pdist(mols, "rocs")
    assert len(oecluster.agglomerative(matrix, linkage="single",
                                       allow_nonmetric=True).labels) == 3


def test_butina_keeps_its_rocs_message():
    # The threshold-graph text is unchanged by the shared helper.
    with pytest.raises(ValueError,
                       match="and the threshold graph compares every pair "
                             "twice. Compute the matrix"):
        oecluster.butina(_mols(FP_SMILES[:3]), 0.5, comparison="rocs")


def test_a_negative_matrix_distance_is_refused_by_hdbscan():
    storage = oecluster.DenseStorage(5)
    for i in range(5):
        for j in range(i + 1, 5):
            storage.Set(i, j, float(j - i))
    storage.Set(1, 3, -0.5)
    matrix = oecluster.SymmetricDistanceMatrix(
        storage, "test", [f"i{k}" for k in range(5)], {})
    for min_samples in (1, 2):
        with pytest.raises(RuntimeError,
                           match="read a negative distance between items 1 "
                                 "and 3"):
            oecluster.hdbscan(matrix, min_cluster_size=2,
                              min_samples=min_samples, allow_nonmetric=True)


def test_single_linkage_ties_follow_the_spanning_tree():
    # d(0,1) = d(1,2) = d(2,3) = 1, every other pair 2: the merge order at a
    # tied height comes from the tree, as documented for 5.20.0.
    storage = oecluster.DenseStorage(4)
    for i in range(4):
        for j in range(i + 1, 4):
            storage.Set(i, j, 1.0 if j == i + 1 else 2.0)
    matrix = oecluster.SymmetricDistanceMatrix(
        storage, "test", [f"i{k}" for k in range(4)], {})
    result = oecluster.agglomerative(matrix, linkage="single", n_clusters=2,
                                     allow_nonmetric=True)
    assert list(result.children) == [(0, 1), (2, 4), (3, 5)]
    assert result.clusters == ((0, 1, 2), (3,))
    cut = oecluster.agglomerative(matrix, linkage="single",
                                  distance_threshold=1.0,
                                  allow_nonmetric=True)
    assert cut.clusters == ((0, 1, 2, 3),)
