import numpy as np
import oecluster
import pytest
from openeye import oechem


def _pose(smiles, shift=0.0, title="mol"):
    """Build a single-conformer OEMol translated `shift` angstroms along x."""
    mol = oechem.OEMol()
    oechem.OESmilesToMol(mol, smiles)
    oechem.OEGenerate2DCoordinates(mol)
    coords = oechem.OEFloatArray(3 * mol.NumAtoms())
    mol.GetCoords(coords)
    if shift:
        for idx in range(0, len(coords), 3):
            coords[idx] += shift
        mol.SetCoords(coords)
    mol.SetTitle(title)
    return mol


def _multiconformer(smiles, shifts, title="mol"):
    """Build one OEMol carrying a conformer per shift."""
    mol = oechem.OEMol()
    oechem.OESmilesToMol(mol, smiles)
    oechem.OEGenerate2DCoordinates(mol)
    base = oechem.OEFloatArray(3 * mol.NumAtoms())
    mol.GetCoords(base)
    for shift in shifts[1:]:
        moved = oechem.OEFloatArray(list(base))
        for idx in range(0, len(moved), 3):
            moved[idx] += shift
        mol.NewConf(moved)
    mol.SetTitle(title)
    return mol


def test_pdist_rmsd_computes_distances():
    mols = [_pose("CCCO", shift, f"p{i}")
            for i, shift in enumerate((0.0, 1.0, 3.0))]
    dist = oecluster.pdist(mols, "rmsd")
    assert dist.num_samples == 3
    assert dist.comparison_name == "rmsd"
    assert np.all(dist.condensed >= 0.0)


def test_in_frame_rmsd_equals_the_translation():
    mols = [_pose("CCCO", 0.0, "a"), _pose("CCCO", 2.0, "b")]
    dist = oecluster.pdist(mols, "rmsd", automorph=False)
    np.testing.assert_allclose(dist.condensed[0], 2.0, atol=1e-4)


def test_overlay_removes_the_translation():
    mols = [_pose("CCCO", 0.0, "a"), _pose("CCCO", 2.0, "b")]
    dist = oecluster.pdist(mols, "rmsd", overlay=True)
    assert dist.condensed[0] < 1e-3


def test_rmsd_is_stamped_zero_self_with_an_unknown_triangle():
    """Geometric comparisons claim zero self-distance and nothing more.

    ``triangle`` is ``"unknown"`` rather than ``True``: OERMSD is not
    documented to attain the exact minimum over the automorphism group, so the
    quotient-metric argument that would prove ``True`` does not close. Unknown
    is permissive, so clustering still runs.
    """
    mols = [_pose("CCCO", shift, f"p{i}")
            for i, shift in enumerate((0.0, 1.0, 3.0))]
    dist = oecluster.pdist(mols, "rmsd")
    assert dist.is_distance is True
    assert dist.metric_capabilities == {'zero_self': True,
                                        'triangle': "unknown"}
    assert dist.data_integrity == "complete"
    oecluster.butina(dist, 1.5)


def test_mismatched_topology_is_rejected():
    mols = [_pose("CCCO", 0.0, "a"), _pose("c1ccccc1", 0.0, "b")]
    with pytest.raises(RuntimeError, match="rocs"):
        oecluster.pdist(mols, "rmsd")


def test_conformers_are_expanded_and_labeled():
    mols = [_multiconformer("CCCO", (0.0, 1.0, 2.0), "ligA"),
            _multiconformer("CCCO", (0.0, 3.0), "ligB")]
    dist = oecluster.pdist(mols, "rmsd")
    assert dist.num_samples == 5
    assert dist.labels == ["ligA:conf0", "ligA:conf1", "ligA:conf2",
                           "ligB:conf0", "ligB:conf1"]


def test_expansion_can_be_disabled():
    mols = [_multiconformer("CCCO", (0.0, 1.0, 2.0), "ligA"),
            _multiconformer("CCCO", (0.0, 3.0), "ligB")]
    dist = oecluster.pdist(mols, "rmsd", expand_conformers=False)
    assert dist.num_samples == 2
    assert dist.labels == ["ligA", "ligB"]


def test_expansion_does_not_modify_the_input():
    mol = _multiconformer("CCCO", (0.0, 1.0, 2.0), "ligA")
    oecluster.pdist([mol], "rmsd")
    assert mol.NumConfs() == 3
    assert mol.GetTitle() == "ligA"


def test_single_conformer_molecules_keep_their_titles():
    mols = [_pose("CCCO", 0.0, "a"), _pose("CCCO", 1.0, "b")]
    dist = oecluster.pdist(mols, "rmsd")
    assert dist.labels == ["a", "b"]


def test_cdist_expands_each_side():
    a = [_multiconformer("CCCO", (0.0, 1.0), "ligA")]
    b = [_pose("CCCO", 2.0, "ligB")]
    cross = oecluster.cdist(a, b, "rmsd")
    assert cross.shape == (2, 1)
    assert cross.labels_a == ["ligA:conf0", "ligA:conf1"]
    assert cross.labels_b == ["ligB"]


def test_cdist_expands_the_b_side():
    a = [_pose("CCCO", 0.0, "ligA")]
    b = [_multiconformer("CCCO", (1.0, 2.0, 3.0), "ligB")]
    cross = oecluster.cdist(a, b, "rmsd")
    assert cross.shape == (1, 3)
    assert cross.labels_a == ["ligA"]
    assert cross.labels_b == ["ligB:conf0", "ligB:conf1", "ligB:conf2"]


def test_cdist_expands_both_sides():
    a = [_multiconformer("CCCO", (0.0, 1.0), "ligA")]
    b = [_multiconformer("CCCO", (2.0, 3.0, 4.0), "ligB")]
    cross = oecluster.cdist(a, b, "rmsd")
    assert cross.shape == (2, 3)
    assert np.isfinite(np.asarray(cross)).all()


def test_a_topology_mismatch_across_the_cdist_boundary_is_rejected():
    """Each side is internally consistent; only A against B disagrees.

    The comparison is built from ``a + b``, so the upfront topology check sees
    the mismatch even though neither side alone would trip it.
    """
    a = [_pose("CCCO", 0.0, "a0"), _pose("CCCO", 1.0, "a1")]
    b = [_pose("c1ccccc1", 0.0, "b0")]
    with pytest.raises(RuntimeError, match="rocs"):
        oecluster.cdist(a, b, "rmsd")


def test_similarity_is_rejected_for_rmsd():
    mols = [_pose("CCCO", 0.0, "a"), _pose("CCCO", 1.0, "b")]
    with pytest.raises(ValueError, match="no similarity form"):
        oecluster.pdist(mols, "rmsd", similarity=True)


def test_unknown_rmsd_kwargs_are_rejected():
    mols = [_pose("CCCO", 0.0, "a"), _pose("CCCO", 1.0, "b")]
    with pytest.raises(TypeError, match="Unknown kwargs for rmsd"):
        oecluster.pdist(mols, "rmsd", bogus=1)


def test_the_factory_class_builds_a_usable_comparison():
    mols = [_pose("CCCO", 0.0, "a"), _pose("CCCO", 1.0, "b")]
    comparison = oecluster.RMSDComparison(mols, automorph=False)
    dist = oecluster.pdist(mols, comparison)
    np.testing.assert_allclose(dist.condensed[0], 1.0, atol=1e-4)


def test_thread_count_does_not_change_the_result():
    mols = [_pose("CCCO", float(i), f"p{i}") for i in range(8)]
    one = oecluster.pdist(mols, "rmsd", num_threads=1).condensed.copy()
    many = oecluster.pdist(mols, "rmsd", num_threads=8).condensed.copy()
    np.testing.assert_array_equal(one, many)
