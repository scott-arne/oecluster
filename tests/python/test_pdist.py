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
