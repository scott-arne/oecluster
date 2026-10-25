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
    import oecluster
    from oecluster import CrossDistanceMatrix
    from openeye import oechem

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
    import oecluster
    from openeye import oechem

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
    import oecluster
    from oecluster import FingerprintComparison
    from openeye import oechem

    mol = oechem.OEGraphMol()
    oechem.OESmilesToMol(mol, "c1ccccc1")
    comp = FingerprintComparison([mol])
    with pytest.raises(TypeError, match="comparison"):
        oecluster.cdist([mol], [mol], comp)


def test_load_distance_matrix_dispatches(tmp_path):
    """load_distance_matrix returns the correct subclass by file kind."""
    import oecluster
    from oecluster import CrossDistanceMatrix, DenseStorage, SymmetricDistanceMatrix

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


def test_cdist_rocs_end_to_end():
    """The rectangular rocs route through the same molecule typemap."""
    pytest.importorskip("openeye.oeomega")
    import oecluster
    from openeye import oechem, oeomega

    omega = oeomega.OEOmega()
    omega.SetMaxConfs(1)
    omega.SetStrictStereo(False)

    def conformer(smiles, title):
        mol = oechem.OEMol()
        oechem.OESmilesToMol(mol, smiles)
        assert omega(mol)
        mol.SetTitle(title)
        return mol

    a = [conformer("c1ccccc1", "benzene")]
    b = [conformer("c1ccc(O)cc1", "phenol"), conformer("CCCCCCCC", "octane")]

    cross = oecluster.cdist(a, b, "rocs", score_type="shape")
    assert cross.comparison_name == "rocs"
    assert cross.shape == (1, 2)
    assert cross.labels_a == ["benzene"]
    assert cross.labels_b == ["phenol", "octane"]

    values = np.asarray(cross)
    assert np.all(np.isfinite(values))
    # Same ordering as the pdist test, reached through the rectangular path.
    assert values[0][0] < 0.1
    assert values[0][1] > 0.4


def _flexible_multiconformer():
    """A molecule whose conformers differ enough to be told apart by shape.

    The same construction as the fixture of the same name in
    ``test_native_bindings.py``, and it self-checks for the same reason: the
    test below tells a preserved ensemble from a collapsed one by the score
    alone, so an ensemble that Omega happened to generate tightly would make it
    agree for the wrong reason rather than fail.

    The check is on the one conformer the test actually reaches, the one
    ``_non_active_reference`` picks, measured against the active conformer that
    survives the collapse. Checking the widest gap in the ensemble instead
    would let a conformer the test never touches satisfy the guard.
    """
    pytest.importorskip("openeye.oeomega")
    import oecluster
    from openeye import oechem, oeomega

    omega = oeomega.OEOmega()
    omega.SetMaxConfs(3)
    omega.SetStrictStereo(False)
    multi = oechem.OEMol()
    oechem.OESmilesToMol(multi, "c1ccccc1CCCCc1ccccc1")
    assert omega(multi)
    assert multi.NumConfs() > 1

    active = oechem.OEMol(multi.GetConf(oechem.OEHasConfIdx(multi.GetActive().GetIdx())))
    gap = oecluster.cdist([_non_active_reference(multi)], [active], "rocs",
                          score_type="shape").matrix[0][0]
    assert gap > 0.1, (
        "the conformer the test using this fixture reaches is too similar to "
        "the active one to distinguish a dropped ensemble from a kept one: "
        f"gap={gap}")
    return multi


def _non_active_reference(multi):
    """An ``OEMol`` of a conformer that is *not* the active one.

    The active conformer is the only pose an ``OEMolBase`` view keeps, so a
    reference drawn from any other one is reachable through the ensemble and
    unreachable through the collapsed copy. That asymmetry is what makes the
    conformer loss visible as a number.
    """
    from openeye import oechem

    active_idx = multi.GetActive().GetIdx()
    others = [conf.GetIdx() for conf in multi.GetConfs() if conf.GetIdx() != active_idx]
    assert others, "fixture produced no conformer other than the active one"
    return oechem.OEMol(multi.GetConf(oechem.OEHasConfIdx(others[-1])))


def test_cdist_rocs_reads_set_b_by_its_own_type_whatever_set_a_holds():
    """cdist builds one comparison over ``items_a + items_b``, and that list's
    first element used to decide how set B was read.

    An ``OEGraphMol`` in set A put the concatenation on a constructor that
    scored set B's ``OEMol`` on its active conformer, and an ``OEMol`` in set A
    put it on one that refused set B's ``OEGraphMol``. The conversion now copies
    each element by its own type, so set B keeps its ensemble whatever set A
    holds, and the reverse arrangement is accepted.

    Asserted as a property rather than as literals. Pinning the numbers would
    let the property break while the figures were updated to match, which is
    the whole failure mode: every call returns a plausible score.
    """
    pytest.importorskip("openeye.oeomega")
    import oecluster
    from openeye import oechem

    multi = _flexible_multiconformer()
    reference = _non_active_reference(multi)
    graph_reference = oechem.OEGraphMol(reference)

    def shape(items_a, items_b):
        return oecluster.cdist(items_a, items_b, "rocs",
                               score_type="shape").matrix[0][0]

    through_graph = shape([graph_reference], [multi])
    through_oemol = shape([reference], [multi])
    assert through_graph == pytest.approx(through_oemol, abs=1e-6)

    # The ensemble is what makes them agree: collapsing set B to its active
    # conformer deliberately gives a different, worse number.
    collapsed = shape([reference], [oechem.OEMol(oechem.OEGraphMol(multi))])
    assert through_oemol < collapsed - 0.1

    # The reverse arrangement used to be refused. ``multi`` is now the
    # reference, and both forms of the fit agree.
    assert shape([multi], [graph_reference]) == pytest.approx(
        shape([multi], [reference]), abs=1e-6)


def test_cdist_mcs_accepts_mixed_sets_in_either_order():
    """The same cross-set case on ``mcs``, which needs no license.

    MCS never reads coordinates, so the order never changed a score; what
    it changed was whether the call was accepted at all. Kept here so the
    cross-set acceptance stays covered on a machine with no Omega license.
    """
    import oecluster
    from openeye import oechem

    benzene = oechem.OEGraphMol()
    oechem.OESmilesToMol(benzene, "c1ccccc1")
    toluene = oechem.OEMol()
    oechem.OESmilesToMol(toluene, "Cc1ccccc1")

    # Benzene's six bonds all match toluene's seven, so the distance is
    # 1 - 6/7 whichever set holds which.
    forward = oecluster.cdist([benzene], [toluene], "mcs").matrix[0][0]
    reverse = oecluster.cdist([toluene], [benzene], "mcs").matrix[0][0]
    assert forward == pytest.approx(1.0 / 7.0, abs=1e-6)
    assert reverse == pytest.approx(1.0 / 7.0, abs=1e-6)


def _two_mols():
    from openeye import oechem

    mols = []
    for smi in ("c1ccccc1", "Cc1ccccc1"):
        mol = oechem.OEGraphMol()
        oechem.OESmilesToMol(mol, smi)
        mols.append(mol)
    return mols


def test_cdist_refuses_a_negative_num_threads():
    """A negative count is a ValueError, as in butina, not the binding's
    OverflowError from a size_t field."""
    import oecluster

    mols = _two_mols()
    with pytest.raises(ValueError, match="num_threads must be non-negative"):
        oecluster.cdist(mols, mols, "fingerprint", num_threads=-1)


@pytest.mark.parametrize("value", [2.5, "8", None], ids=["float", "str", "none"])
def test_cdist_refuses_a_non_integer_num_threads(value):
    """A float, a string or None is refused rather than truncated."""
    import oecluster

    mols = _two_mols()
    with pytest.raises(TypeError):
        oecluster.cdist(mols, mols, "fingerprint", num_threads=value)


def test_cdist_accepts_a_numpy_integer_num_threads():
    import oecluster

    mols = _two_mols()
    expected = oecluster.cdist(mols, mols, "fingerprint", num_threads=2)
    got = oecluster.cdist(mols, mols, "fingerprint", num_threads=np.int64(2))  # pyright: ignore[reportArgumentType]
    np.testing.assert_array_equal(got.matrix, expected.matrix)
