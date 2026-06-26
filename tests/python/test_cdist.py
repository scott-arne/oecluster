import numpy as np
import pytest


def test_cross_distance_matrix_basics():
    """CrossDistanceMatrix wraps a rectangular array with metadata."""
    from oecluster import CrossDistanceMatrix, DistanceMatrix

    mat = np.array([[0.1, 0.2, 0.3], [0.4, 0.5, 0.6]], dtype=np.float64)
    cdm = CrossDistanceMatrix(mat, "fingerprint", ["a0", "a1"], ["b0", "b1", "b2"], {})

    assert isinstance(cdm, DistanceMatrix)
    assert cdm.shape == (2, 3)
    assert cdm.labels_a == ["a0", "a1"]
    assert cdm.labels_b == ["b0", "b1", "b2"]
    assert cdm.comparison_name == "fingerprint"
    assert len(cdm) == 2
    np.testing.assert_array_equal(np.asarray(cdm), mat)
    assert repr(cdm) == "CrossDistanceMatrix(comparison='fingerprint', shape=(2, 3))"


def test_cross_distance_matrix_roundtrip(tmp_path):
    """CrossDistanceMatrix serializes and reloads via .npz."""
    from oecluster import CrossDistanceMatrix

    mat = np.array([[0.1, 0.2], [0.3, 0.4], [0.5, 0.6]], dtype=np.float64)
    cdm = CrossDistanceMatrix(mat, "rocs", ["a", "b", "c"], ["x", "y"], {"k": 1})
    path = tmp_path / "cross.npz"
    cdm.to_file(str(path))

    loaded = CrossDistanceMatrix.from_file(str(path))
    np.testing.assert_array_equal(loaded.matrix, mat)
    assert loaded.shape == (3, 2)
    assert loaded.labels_a == ["a", "b", "c"]
    assert loaded.labels_b == ["x", "y"]
    assert loaded.comparison_name == "rocs"
    assert loaded.params == {"k": 1}


def test_cross_from_file_rejects_symmetric(tmp_path):
    """CrossDistanceMatrix.from_file rejects a symmetric (no matrix_kind) file."""
    from oecluster import CrossDistanceMatrix, DenseStorage, SymmetricDistanceMatrix

    storage = DenseStorage(3)
    storage.Set(0, 1, 0.5)
    storage.Set(0, 2, 0.25)
    storage.Set(1, 2, 0.75)
    sdm = SymmetricDistanceMatrix(storage, "test", ["a", "b", "c"], {})
    path = tmp_path / "sym.npz"
    sdm.to_file(str(path))

    with pytest.raises(ValueError, match="symmetric"):
        CrossDistanceMatrix.from_file(str(path))


def test_cross_from_file_rejects_malformed(tmp_path):
    """A cross file with a wrong-shape matrix raises a descriptive ValueError."""
    from oecluster import CrossDistanceMatrix

    path = tmp_path / "bad.npz"
    # matrix_kind says cross, but matrix shape disagrees with n_a/n_b.
    np.savez_compressed(
        path,
        matrix_kind=np.array("cross"),
        matrix=np.zeros((2, 2), dtype=np.float64),
        n_a=np.array(2),
        n_b=np.array(3),
        labels_a=np.array(["a", "b"]),
        labels_b=np.array(["x", "y", "z"]),
        comparison_name=np.array("fingerprint"),
        params_json=np.array("{}"),
    )
    with pytest.raises(ValueError, match="Malformed cross matrix"):
        CrossDistanceMatrix.from_file(str(path))


def test_cdist_fingerprint_matches_scipy():
    """cdist cross-distances match a scipy reference for fingerprints."""
    from openeye import oechem
    import oecluster
    from oecluster import CrossDistanceMatrix

    smiles_a = ["c1ccccc1", "CCCCCCCC"]
    smiles_b = ["c1ccc(O)cc1", "c1ccncc1", "CCO"]

    def build(smis):
        mols = []
        for smi in smis:
            mol = oechem.OEGraphMol()
            oechem.OESmilesToMol(mol, smi)
            mols.append(mol)
        return mols

    mols_a = build(smiles_a)
    mols_b = build(smiles_b)

    result = oecluster.cdist(mols_a, mols_b, "fingerprint")
    assert isinstance(result, CrossDistanceMatrix)
    assert result.shape == (2, 3)
    assert result.labels_a == ["mol_0", "mol_1"]
    assert result.labels_b == ["mol_0", "mol_1", "mol_2"]

    # Reference: full symmetric pdist over the concatenated set, then slice the
    # rectangular A-vs-B block out of the squareform.
    all_mols = mols_a + mols_b
    full = oecluster.pdist(all_mols, "fingerprint").squareform()
    expected = full[:2, 2:]
    np.testing.assert_allclose(result.matrix, expected, atol=1e-9)


def test_cdist_orientation_is_a_rows_b_cols():
    """Entry [i, j] pairs items_a[i] with items_b[j]."""
    from openeye import oechem
    import oecluster

    def build(smis):
        mols = []
        for smi in smis:
            mol = oechem.OEGraphMol()
            oechem.OESmilesToMol(mol, smi)
            mols.append(mol)
        return mols

    a = build(["c1ccccc1", "CCCCCCCC", "CCO"])
    b = build(["c1ccccc1"])  # identical to a[0]
    result = oecluster.cdist(a, b, "fingerprint")
    assert result.shape == (3, 1)
    # a[0] vs b[0] are the same molecule -> distance ~0; others larger.
    assert result.matrix[0, 0] == pytest.approx(0.0, abs=1e-9)
    assert result.matrix[1, 0] > result.matrix[0, 0]


def test_cdist_empty_input_raises():
    """Empty set A or B raises ValueError before any C++ call."""
    import oecluster
    with pytest.raises(ValueError, match="empty"):
        oecluster.cdist([], [1], "fingerprint")
    with pytest.raises(ValueError, match="empty"):
        oecluster.cdist([1], [], "fingerprint")


def test_cdist_similarity_with_cutoff_raises():
    """cutoff > 0 with similarity=True is rejected."""
    import oecluster
    with pytest.raises(ValueError, match="cutoff"):
        oecluster.cdist([1], [1], "fingerprint", similarity=True, cutoff=0.5)


def test_cdist_rejects_prebuilt_comparison_object():
    """A non-string comparison is rejected by cdist."""
    from openeye import oechem
    import oecluster
    from oecluster import FingerprintComparison

    mol = oechem.OEGraphMol()
    oechem.OESmilesToMol(mol, "c1ccccc1")
    comp = FingerprintComparison([mol])
    with pytest.raises(TypeError, match="comparison"):
        oecluster.cdist([mol], [mol], comp)


def test_load_distance_matrix_dispatches(tmp_path):
    """load_distance_matrix returns the correct subclass by file kind."""
    import oecluster
    from oecluster import (CrossDistanceMatrix, DenseStorage,
                           SymmetricDistanceMatrix)

    storage = DenseStorage(3)
    storage.Set(0, 1, 0.5)
    storage.Set(0, 2, 0.25)
    storage.Set(1, 2, 0.75)
    sym = SymmetricDistanceMatrix(storage, "test", ["a", "b", "c"], {})
    sym_path = tmp_path / "sym.npz"
    sym.to_file(str(sym_path))

    mat = np.array([[0.1, 0.2], [0.3, 0.4]], dtype=np.float64)
    cross = CrossDistanceMatrix(mat, "rocs", ["a", "b"], ["x", "y"], {})
    cross_path = tmp_path / "cross.npz"
    cross.to_file(str(cross_path))

    loaded_sym = oecluster.load_distance_matrix(str(sym_path))
    loaded_cross = oecluster.load_distance_matrix(str(cross_path))

    assert isinstance(loaded_sym, SymmetricDistanceMatrix)
    assert isinstance(loaded_cross, CrossDistanceMatrix)
    assert loaded_sym.num_samples == 3
    assert loaded_cross.shape == (2, 2)


def test_cross_matrix_rejected_by_clustering():
    """Passing a CrossDistanceMatrix to any clustering/report entry point raises TypeError."""
    import oecluster
    from oecluster import CrossDistanceMatrix, DenseStorage, SymmetricDistanceMatrix

    mat = np.array([[0.1, 0.2], [0.3, 0.4]], dtype=np.float64)
    cross = CrossDistanceMatrix(mat, "fingerprint", ["a", "b"], ["x", "y"], {})

    # Guard sites whose first inspected argument is distance_matrix.
    with pytest.raises(TypeError):
        oecluster.butina(cross, 0.5)
    with pytest.raises(TypeError):
        oecluster.dbscan(cross, 0.5)
    with pytest.raises(TypeError):
        oecluster.hdbscan(cross)
    with pytest.raises(TypeError):
        oecluster.agglomerative(cross, n_clusters=2)
    # representative(cluster, distance_matrix): the distance_matrix guard runs
    # before any cluster work.
    with pytest.raises(TypeError):
        oecluster.representative((0, 1), cross)

    # cluster_report(result, distance_matrix): the result-type guard runs FIRST,
    # so a valid ClusteringResult is required to reach the distance_matrix guard.
    # Build a real result from a small symmetric matrix.
    storage = DenseStorage(4)
    for i in range(4):
        for j in range(i + 1, 4):
            storage.Set(i, j, 0.5)
    sym = SymmetricDistanceMatrix(storage, "test", ["a", "b", "c", "d"], {})
    result = oecluster.butina(sym, 0.6)
    with pytest.raises(TypeError):
        oecluster.cluster_report(result, cross)
