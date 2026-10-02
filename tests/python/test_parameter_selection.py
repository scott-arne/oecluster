"""Python surface of ClusteringSpec and select_parameter."""
import numpy as np
import oecluster
import oefp
import pytest
from openeye import oechem

FP_SMILES = ["CCO", "CCCO", "CCCCO", "c1ccccc1", "Cc1ccccc1", "CCc1ccccc1",
             "CC(=O)O", "CC(=O)OC", "CCN", "CCCN", "C1CCCCC1", "c1ccncc1"]

ROSTER = ("butina", "dbscan", "hdbscan", "agglomerative", "k_medoids",
          "bitbirch", "bitbirch_recluster", "bitbirch_refine",
          "sphere_exclusion", "jarvis_patrick", "leiden", "murcko")

# Butina thresholds that give ten singletons, the two blobs, and one cluster.
GRID = (0.05, 0.2, 0.95)


def _mols(smiles_list):
    mols = []
    for idx, smi in enumerate(smiles_list):
        mol = oechem.OEGraphMol()
        oechem.OESmilesToMol(mol, smi)
        mol.SetTitle(f"mol{idx}")
        mols.append(mol)
    return mols


def _blobs():
    """Ten points, two blobs of five: intra-blob 0.1, inter-blob 0.9."""
    values = [0.1 if (i < 5) == (j < 5) else 0.9
              for i in range(10) for j in range(i + 1, 10)]
    return oecluster.SymmetricDistanceMatrix.from_condensed(np.array(values))


def _batch(bits, rows):
    return oefp.OEFPBatch.from_fingerprints(
        [oefp.OEFP.from_on_bits(bits, list(on)) for on in rows])


def _batch_from_bits(bits):
    arr = np.asarray(bits, dtype=np.uint8)
    fingerprints = []
    for row in arr:
        on_bits = np.flatnonzero(row).astype(int).tolist()
        fingerprints.append(oefp.OEFP.from_on_bits(arr.shape[1], on_bits))
    return oefp.OEFPBatch.from_fingerprints(fingerprints)


def _fps():
    """Two groups of five sharing an eight-bit block, one private bit each."""
    rows = ([set(range(8)) | {16 + i} for i in range(5)]
            + [set(range(32, 40)) | {48 + i} for i in range(5)])
    return _batch(64, rows)


def _pairs():
    """p0 == p1 and p2 == p3 (distance 0); every other distance 0.5."""
    zero = {(0, 1), (2, 3)}
    values = [0.0 if (i, j) in zero else 0.5
              for i in range(4) for j in range(i + 1, 4)]
    return oecluster.SymmetricDistanceMatrix.from_condensed(np.array(values))


def split_pairs(items, split):
    """A foreign clusterer: 'across' pairs coincident medoids, 'together' does not."""
    if split == "across":
        return oecluster.ClusteringResult([0, 1, 0, 1], ((0, 2), (1, 3)))
    return oecluster.ClusteringResult([0, 0, 1, 1], ((0, 1), (2, 3)))


def _recording(calls):
    """A clusterer that records the options it was called with."""
    def record(items, **options):
        calls.append(options)
        return oecluster.ClusteringResult([0, 0, 0, 0], ((0, 1, 2, 3),))
    return record


# --- ClusteringSpec ---------------------------------------------------------

@pytest.mark.parametrize("name", ROSTER)
def test_spec_resolves_every_roster_name(name):
    spec = oecluster.ClusteringSpec(name)
    assert spec.name == name
    assert spec.algorithm is getattr(oecluster, name)


def test_spec_name_matching_is_case_insensitive():
    assert oecluster.ClusteringSpec("Butina").name == "butina"


def test_spec_maps_a_roster_function_back_to_its_name():
    spec = oecluster.ClusteringSpec(oecluster.sphere_exclusion)
    assert spec.name == "sphere_exclusion"
    assert spec.algorithm is oecluster.sphere_exclusion


def test_spec_keeps_a_foreign_callable_and_its_name():
    spec = oecluster.ClusteringSpec(split_pairs, split="across")
    assert spec.name == "split_pairs"
    assert spec.algorithm is split_pairs
    assert dict(spec.options) == {"split": "across"}


def test_spec_refuses_an_unknown_name_and_lists_the_roster():
    with pytest.raises(ValueError, match="butina"):
        oecluster.ClusteringSpec("kmeans")


def test_spec_refuses_a_non_callable():
    with pytest.raises(TypeError, match="callable"):
        oecluster.ClusteringSpec(42)


def test_spec_options_are_a_read_only_copy():
    options = {"num_threads": 2}
    spec = oecluster.ClusteringSpec("butina", **options)
    options["num_threads"] = 8
    assert dict(spec.options) == {"num_threads": 2}
    with pytest.raises(TypeError):
        spec.options["num_threads"] = 4


def test_run_passes_fixed_options_and_overrides_win():
    calls = []
    spec = oecluster.ClusteringSpec(_recording(calls), threshold=0.2,
                                    num_threads=1)
    spec.run(_blobs())
    spec.run(_blobs(), threshold=0.5)
    assert calls == [{"threshold": 0.2, "num_threads": 1},
                     {"threshold": 0.5, "num_threads": 1}]


def test_run_returns_the_algorithm_result():
    result = oecluster.ClusteringSpec("butina").run(_blobs(), threshold=0.2)
    assert isinstance(result, oecluster.ButinaResult)
    assert result.num_clusters == 2


def test_run_refuses_a_callable_that_returns_no_result():
    spec = oecluster.ClusteringSpec(lambda items, **options: None)
    with pytest.raises(TypeError, match="ClusteringResult"):
        spec.run(_blobs())


def test_run_refuses_a_graph_as_a_partition():
    # knn_graph is public and callable but returns a KNNGraph, not a result.
    with pytest.raises(TypeError, match="ClusteringResult"):
        oecluster.ClusteringSpec(oecluster.knn_graph).run(_blobs(), k=3)


def test_spec_equality_and_unhashability():
    spec = oecluster.ClusteringSpec("butina", num_threads=2)
    assert spec == oecluster.ClusteringSpec(oecluster.butina, num_threads=2)
    assert spec != oecluster.ClusteringSpec("butina", num_threads=3)
    assert spec != oecluster.ClusteringSpec("dbscan", num_threads=2)
    assert spec != "butina"
    with pytest.raises(TypeError, match="unhashable"):
        hash(spec)


def test_spec_equality_handles_array_valued_options():
    # numpy's == is elementwise, so plain dict equality would raise here.
    left = oecluster.ClusteringSpec("k_medoids", initial_medoids=np.array([0, 5]))
    right = oecluster.ClusteringSpec("k_medoids", initial_medoids=np.array([0, 5]))
    assert left == right
    assert left != oecluster.ClusteringSpec("k_medoids", initial_medoids=np.array([0, 6]))
    assert left != oecluster.ClusteringSpec("k_medoids", initial_medoids=np.array([0, 5, 9]))
    assert left != oecluster.ClusteringSpec("k_medoids", initial_medoids=5)


def test_spec_repr():
    assert (repr(oecluster.ClusteringSpec("butina", num_threads=4))
            == "ClusteringSpec('butina', num_threads=4)")
    assert repr(oecluster.ClusteringSpec("leiden")) == "ClusteringSpec('leiden')"


def test_package_exports_clustering_spec():
    assert "ClusteringSpec" in oecluster.__all__
