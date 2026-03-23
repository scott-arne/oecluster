import pytest
import numpy as np

def test_import():
    import oecluster
    assert hasattr(oecluster, '__version__')
    assert hasattr(oecluster, 'pdist')

def test_dense_storage_roundtrip():
    from oecluster._oecluster import DenseStorage
    s = DenseStorage(4)
    s.Set(0, 1, 0.5)
    assert s.Get(0, 1) == pytest.approx(0.5)
    assert s.NumItems() == 4
    assert s.NumPairs() == 6

def test_pdist_fingerprint():
    """Test pdist with fingerprint metric on simple molecules."""
    from openeye import oechem
    import oecluster

    smiles = ["c1ccccc1", "c1ccc(O)cc1", "CCCCCCCC", "c1ccncc1"]
    mols = []
    for smi in smiles:
        mol = oechem.OEGraphMol()
        oechem.OESmilesToMol(mol, smi)
        mols.append(mol)

    dist = oecluster.pdist(mols, "fingerprint")
    assert dist.num_items == 4
    assert len(dist) == 6  # 4*3/2
    assert dist.metric_name == "fingerprint"

    # Check numpy array protocol
    arr = np.asarray(dist)
    assert arr.shape == (6,)
    assert arr.dtype == np.float64

    # Check squareform
    sq = dist.squareform()
    assert sq.shape == (4, 4)
    assert sq[0, 0] == 0.0
    assert sq[0, 1] == sq[1, 0]  # symmetric

def test_pdist_with_cutoff():
    """Test pdist with cutoff produces sparse result."""
    from openeye import oechem
    import oecluster

    smiles = ["c1ccccc1", "c1ccc(O)cc1", "CCCCCCCC"]
    mols = []
    for smi in smiles:
        mol = oechem.OEGraphMol()
        oechem.OESmilesToMol(mol, smi)
        mols.append(mol)

    dist = oecluster.pdist(mols, "fingerprint", cutoff=0.5)
    assert dist.num_items == 3

def test_pdist_progress():
    """Test progress callback is invoked."""
    from openeye import oechem
    import oecluster

    smiles = ["c1ccccc1", "c1ccc(O)cc1", "CCCCCCCC"]
    mols = []
    for smi in smiles:
        mol = oechem.OEGraphMol()
        oechem.OESmilesToMol(mol, smi)
        mols.append(mol)

    calls = []
    dist = oecluster.pdist(mols, "fingerprint",
                           progress=lambda d, t: calls.append((d, t)))
    assert len(calls) > 0

def test_distance_matrix_serialization(tmp_path):
    """Test DistanceMatrix save/load roundtrip."""
    from openeye import oechem
    import oecluster

    smiles = ["c1ccccc1", "c1ccc(O)cc1", "CCCCCCCC"]
    mols = []
    for smi in smiles:
        mol = oechem.OEGraphMol()
        oechem.OESmilesToMol(mol, smi)
        mols.append(mol)

    dist = oecluster.pdist(mols, "fingerprint")
    path = str(tmp_path / "test.npz")
    dist.to_file(path)

    loaded = oecluster.DistanceMatrix.from_file(path)
    np.testing.assert_array_almost_equal(
        np.asarray(dist), np.asarray(loaded)
    )


def test_pdist_fingerprint_similarity():
    """Test pdist with fingerprint metric in similarity mode."""
    from openeye import oechem
    import oecluster

    smiles = ["c1ccccc1", "c1ccc(O)cc1", "CCCCCCCC"]
    mols = []
    for smi in smiles:
        mol = oechem.OEGraphMol()
        oechem.OESmilesToMol(mol, smi)
        mols.append(mol)

    dist = oecluster.pdist(mols, "fingerprint", similarity=True)
    arr = np.asarray(dist)
    # Tanimoto similarity should be in [0, 1]
    assert np.all(arr >= 0.0)
    assert np.all(arr <= 1.0)

    # Compare with distance mode: sim + dist should equal 1.0
    dist_d = oecluster.pdist(mols, "fingerprint", similarity=False)
    arr_d = np.asarray(dist_d)
    np.testing.assert_array_almost_equal(arr + arr_d, np.ones_like(arr))


def test_pdist_superpose_sitehopper_alias():
    """Test that 'sitehopper' routes to superpose with method=sitehopper."""
    from openeye import oechem
    import oecluster

    files = [
        "tests/assets/spruce_5FQD_1_5FQD_1-ALIGNED_BC__DU__LVY_B-1438.oedu",
        "tests/assets/spruce_8G66_1_8G66_1-ALIGNED_BC__DU__YOT_B-502.oedu",
    ]
    dus = []
    for f in files:
        du = oechem.OEDesignUnit()
        if not oechem.OEReadDesignUnit(f, du):
            pytest.skip(f"Cannot read {f}")
        dus.append(du)

    dm = oecluster.pdist(dus, "sitehopper")
    assert "sitehopper" in dm.metric_name


def test_pdist_superpose_kwargs():
    """Test superpose with method kwarg."""
    from openeye import oechem
    import oecluster

    files = [
        "tests/assets/spruce_5FQD_1_5FQD_1-ALIGNED_BC__DU__LVY_B-1438.oedu",
        "tests/assets/spruce_8G66_1_8G66_1-ALIGNED_BC__DU__YOT_B-502.oedu",
    ]
    dus = []
    for f in files:
        du = oechem.OEDesignUnit()
        if not oechem.OEReadDesignUnit(f, du):
            pytest.skip(f"Cannot read {f}")
        dus.append(du)

    dm = oecluster.pdist(dus, "superpose", method="ddm")
    assert dm.metric_name == "superpose:ddm"


def test_pdist_superpose_predicate():
    """Test superpose with predicate kwarg."""
    from openeye import oechem
    import oecluster

    files = [
        "tests/assets/spruce_5FQD_1_5FQD_1-ALIGNED_BC__DU__LVY_B-1438.oedu",
        "tests/assets/spruce_8G66_1_8G66_1-ALIGNED_BC__DU__YOT_B-502.oedu",
    ]
    dus = []
    for f in files:
        du = oechem.OEDesignUnit()
        if not oechem.OEReadDesignUnit(f, du):
            pytest.skip(f"Cannot read {f}")
        dus.append(du)

    dm = oecluster.pdist(dus, "superpose", predicate="backbone")
    arr = np.asarray(dm)
    assert arr.shape == (1,)
    assert np.isfinite(arr[0])


def test_pdist_superpose_similarity():
    """Test superpose with similarity=True for SSE method."""
    from openeye import oechem
    import oecluster

    files = [
        "tests/assets/spruce_5FQD_1_5FQD_1-ALIGNED_BC__DU__LVY_B-1438.oedu",
        "tests/assets/spruce_8G66_1_8G66_1-ALIGNED_BC__DU__YOT_B-502.oedu",
    ]
    dus = []
    for f in files:
        du = oechem.OEDesignUnit()
        if not oechem.OEReadDesignUnit(f, du):
            pytest.skip(f"Cannot read {f}")
        dus.append(du)

    dm = oecluster.pdist(dus, "superpose", method="sse", similarity=True)
    arr = np.asarray(dm)
    # SSE Tanimoto similarity in [0, 1]
    assert np.all(arr >= 0.0)
    assert np.all(arr <= 1.0)


def test_pdist_unknown_kwargs_raises():
    """Test that unknown kwargs raise TypeError."""
    from openeye import oechem
    import oecluster

    mols = [oechem.OEGraphMol(), oechem.OEGraphMol()]
    oechem.OESmilesToMol(mols[0], "C")
    oechem.OESmilesToMol(mols[1], "CC")

    with pytest.raises(TypeError, match="Unknown kwargs"):
        oecluster.pdist(mols, "fingerprint", bogus_option=42)
