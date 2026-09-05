import numpy as np
import oecluster
import pytest
from oecluster import _comparisons
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


def test_each_expanded_pose_carries_its_own_coordinates():
    """Counting the poses is not enough; the coordinates have to move too.

    An expansion that copied the source molecule's active conformer into every
    pose would satisfy the count and the labels above and still produce an
    all-zeros matrix. The shifts are chosen so each pairwise distance is the
    translation between two poses, which detects that each pose carries a
    distinct conformer's coordinates rather than a shared one.
    """
    mols = [_multiconformer("CCCO", (0.0, 1.0, 3.0), "ligA")]
    dist = oecluster.pdist(mols, "rmsd", automorph=False)
    np.testing.assert_allclose(dist.condensed, [1.0, 3.0, 2.0], atol=1e-4)


def test_expanded_poses_have_exact_source_conformer_coordinates():
    """Per-atom coordinate equality, not just a derived distance.

    Verifies that expansion copies coordinates correctly, atom by atom, even
    when atom indices have gaps. A per-atom check detects misassignments that
    could still yield plausible overall distances.
    """
    mol = oechem.OEMol()
    # Interleaved hydrogens create gaps after suppression
    oechem.OESmilesToMol(mol, "C([H])([H])([H])C([H])([H])O[H]")
    oechem.OEGenerate2DCoordinates(mol)

    # Add two more conformers with known shifts
    base = oechem.OEFloatArray(mol.GetMaxAtomIdx() * 3)
    mol.GetCoords(base)
    for shift in [0.5, 1.2]:
        moved = oechem.OEFloatArray(list(base))
        for idx in range(0, len(moved), 3):
            moved[idx] += shift
        mol.NewConf(moved)

    # Suppress hydrogens to create index gaps
    oechem.OESuppressHydrogens(mol)
    assert mol.GetMaxAtomIdx() > mol.NumAtoms()

    # Extract source conformer coordinates for comparison
    source_coords = []
    for conf in mol.GetConfs():
        per_atom = []
        for atom in conf.GetAtoms():
            coords = oechem.OEFloatArray(3)
            conf.GetCoords(atom, coords)
            per_atom.append((coords[0], coords[1], coords[2]))
        source_coords.append(per_atom)

    # Expand and verify coordinates match
    expanded, _ = _comparisons._normalize_rmsd([mol], {})
    assert len(expanded) == 3

    for i, pose in enumerate(expanded):
        for j, atom in enumerate(pose.GetAtoms()):
            coords = oechem.OEFloatArray(3)
            pose.GetCoords(atom, coords)
            actual = (coords[0], coords[1], coords[2])
            expected = source_coords[i][j]
            np.testing.assert_allclose(actual, expected, atol=1e-8,
                err_msg=f"Pose {i}, atom {j}: coordinates mismatch")


def test_expansion_with_atom_index_gaps():
    """Molecule with GetMaxAtomIdx() > NumAtoms() must still expand correctly.

    Atom deletions (hydrogen suppression, salt stripping) leave index gaps.
    Coordinate mapping must use iteration order, not raw indices, or the
    assignment silently scrambles atoms.

    The conformers here differ by a per-atom displacement rather than a rigid
    shift, and that is load-bearing rather than decorative. This construction
    builds the conformers before suppression, so a rigid shift moves the
    soon-to-be-dead slots along with the live ones; the raw-index misread then
    comes out translated by that same vector and every pairwise distance
    survives. Measured on this construction: with the raw-index bug reinstated,
    a rigid shift passes and the per-atom displacement fails.

    A consistent relabelling of the atoms is a different mutation, and this
    test does not catch it (measured). Others do, among them
    ``test_expanded_poses_have_exact_source_conformer_coordinates`` and
    ``test_mixing_single_and_multiconformer_molecules``.
    """
    mol = oechem.OEMol()
    # Interleaved hydrogens force GetMaxAtomIdx=9 when NumAtoms=3 after suppression
    oechem.OESmilesToMol(mol, "C([H])([H])([H])C([H])([H])O[H]")
    oechem.OEGenerate2DCoordinates(mol)

    # Build three conformers, each atom displaced by its own amount
    base = oechem.OEFloatArray(mol.GetMaxAtomIdx() * 3)
    mol.GetCoords(base)
    for scale in [1.0, 2.5]:
        moved = oechem.OEFloatArray(list(base))
        for slot in range(mol.GetMaxAtomIdx()):
            moved[slot * 3] += scale * (slot + 1)
            moved[slot * 3 + 1] += scale * 0.3 * slot
        mol.NewConf(moved)

    oechem.OESuppressHydrogens(mol)
    assert mol.NumAtoms() == 3
    assert mol.GetMaxAtomIdx() > 3  # gaps present

    # Ground truth: direct OERMSD on the source conformers
    expected_dists = []
    confs = list(mol.GetConfs())
    for i in range(len(confs)):
        for j in range(i + 1, len(confs)):
            d = oechem.OERMSD(confs[i], confs[j], True, True, False)
            expected_dists.append(d)

    dist = oecluster.pdist([mol], "rmsd")
    np.testing.assert_allclose(dist.condensed, expected_dists, atol=1e-4)


def test_mixing_single_and_multiconformer_molecules():
    """Single-conformer inputs pass through; multi-conformer are rebuilt.

    Mixing the two paths is what gives this test its reach: the single-conformer
    molecule bypasses the expansion, so the distances measured against it can
    register a mistake made inside it. Measured: reversing the atom pairing in
    the expansion fails this test, and the raw-index bug that
    ``test_expansion_with_atom_index_gaps`` covers does not.
    """
    single = _pose("CCCO", 0.5, "single")
    multi = _multiconformer("CCCO", (0.0, 2.0), "multi")

    dist = oecluster.pdist([single, multi], "rmsd", automorph=False)
    assert dist.num_samples == 3
    assert dist.labels == ["single", "multi:conf0", "multi:conf1"]
    # single at 0.5, multi:conf0 at 0.0 (base 2D), multi:conf1 at 0.0+2.0=2.0
    # Condensed order: (0,1), (0,2), (1,2)
    # = (single,multi:conf0), (single,multi:conf1), (multi:conf0,multi:conf1)
    # = dist(0.5,0.0), dist(0.5,2.0), dist(0.0,2.0)
    # = 0.5, 1.5, 2.0
    np.testing.assert_allclose(dist.condensed, [0.5, 1.5, 2.0], atol=1e-4)


def test_expansion_can_be_disabled():
    mols = [_multiconformer("CCCO", (0.0, 1.0, 2.0), "ligA"),
            _multiconformer("CCCO", (0.0, 3.0), "ligB")]
    dist = oecluster.pdist(mols, "rmsd", expand_conformers=False)
    assert dist.num_samples == 2
    assert dist.labels == ["ligA", "ligB"]


@pytest.mark.parametrize("value", ["false", "", 0, 1, 2.5])
def test_non_bool_expand_conformers_is_refused(value):
    """A non-bool must not settle the expansion by truthiness.

    Falsy strings and ints are refused alongside truthy ones: ``"false"``
    would enable the expansion and ``0`` would disable it, in both cases by
    coincidence rather than by what the caller wrote.
    """
    mols = [_multiconformer("CCCO", (0.0, 1.0), "ligA")]
    with pytest.raises(TypeError, match="expand_conformers must be True or"):
        oecluster.pdist(mols, "rmsd", expand_conformers=value)


def test_build_comparison_also_refuses_a_non_bool_expand_conformers():
    """The validator, not the normalizer, owns this refusal.

    ``build_comparison`` never calls the normalizer, so a check living only
    there would let the direct builder path accept the value and silently
    discard it.
    """
    mols = [_pose("CCCO", 0.0, "a"), _pose("CCCO", 1.0, "b")]
    with pytest.raises(TypeError, match="expand_conformers must be True or"):
        _comparisons.build_comparison(
            mols, "rmsd", False, {"expand_conformers": "false"},
            symmetric=True)


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
    cross = oecluster.cdist(a, b, "rmsd", automorph=False)
    assert cross.shape == (2, 1)
    assert cross.labels_a == ["ligA:conf0", "ligA:conf1"]
    assert cross.labels_b == ["ligB"]
    np.testing.assert_allclose(cross.matrix, [[2.0], [1.0]], atol=1e-4)


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
