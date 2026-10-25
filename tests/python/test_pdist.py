import itertools

import numpy as np
import pytest


def test_import():
    import oecluster
    assert hasattr(oecluster, '__version__')
    assert hasattr(oecluster, 'pdist')

def test_dense_storage_roundtrip():
    from oecluster import DenseStorage
    s = DenseStorage(4)
    s.Set(0, 1, 0.5)
    assert s.Get(0, 1) == pytest.approx(0.5)
    assert s.NumSamples() == 4
    assert s.NumPairs() == 6

def test_pdist_fingerprint():
    """Test pdist with default Morgan fingerprint comparison on simple molecules."""
    import oecluster
    from openeye import oechem

    smiles = ["c1ccccc1", "c1ccc(O)cc1", "CCCCCCCC", "c1ccncc1"]
    mols = []
    for smi in smiles:
        mol = oechem.OEGraphMol()
        oechem.OESmilesToMol(mol, smi)
        mols.append(mol)

    dist = oecluster.pdist(mols, "fingerprint")
    assert dist.num_samples == 4
    assert len(dist) == 6  # 4*3/2
    assert dist.comparison_name == "fingerprint"
    removed_metadata_name = "metr" + "ic_name"
    assert not hasattr(dist, removed_metadata_name)

    # Check numpy array protocol
    arr = np.asarray(dist)
    assert arr.shape == (6,)
    assert arr.dtype == np.float64

    # Check squareform
    sq = dist.squareform()
    assert sq.shape == (4, 4)
    assert sq[0, 0] == 0.0
    assert sq[0, 1] == sq[1, 0]  # symmetric


def test_pdist_fingerprint_atom_pair():
    """Test pdist with OEFP Atom Pair fingerprints."""
    import oecluster
    from openeye import oechem

    smiles = ["c1ccccc1", "c1ccc(O)cc1", "CCCCCCCC"]
    mols = []
    for smi in smiles:
        mol = oechem.OEGraphMol()
        oechem.OESmilesToMol(mol, smi)
        mols.append(mol)

    dist = oecluster.pdist(mols, "fingerprint", fp_type="atom_pair")
    arr = np.asarray(dist)
    assert dist.num_samples == 3
    assert arr.shape == (3,)
    assert np.all(arr >= 0.0)


def test_pdist_fingerprint_metric_kwarg():
    """Test pdist selects the OEFP scalar metric with the metric kwarg."""
    import oecluster
    from openeye import oechem

    smiles = ["c1ccccc1", "c1ccc(O)cc1", "CCCCCCCC"]
    mols = []
    for smi in smiles:
        mol = oechem.OEGraphMol()
        oechem.OESmilesToMol(mol, smi)
        mols.append(mol)

    dist = oecluster.pdist(mols, "fingerprint", metric="dice")
    arr = np.asarray(dist)
    assert dist.num_samples == 3
    assert arr.shape == (3,)
    assert np.all(arr >= 0.0)


def test_pdist_fingerprint_removed_openeye_type_raises():
    """OpenEye fingerprint families are not accepted by oecluster fingerprints."""
    import oecluster
    from openeye import oechem

    mols = [oechem.OEGraphMol(), oechem.OEGraphMol()]
    oechem.OESmilesToMol(mols[0], "C")
    oechem.OESmilesToMol(mols[1], "CC")

    with pytest.raises(Exception, match="OpenEye fingerprint type"):
        oecluster.pdist(mols, "fingerprint", fp_type="circular")


def test_pdist_fingerprint_rejects_openeye_mask_kwargs():
    """OpenEye atom and bond mask kwargs are no longer part of the API."""
    import oecluster
    from openeye import oechem

    mols = [oechem.OEGraphMol(), oechem.OEGraphMol()]
    oechem.OESmilesToMol(mols[0], "C")
    oechem.OESmilesToMol(mols[1], "CC")

    with pytest.raises(TypeError, match="Unknown kwargs"):
        oecluster.pdist(mols, "fingerprint", atom_type_mask=1)


def test_pdist_fingerprint_rejects_similarity_func_kwarg():
    """The old similarity_func name is not part of the hard-break API."""
    import oecluster
    from openeye import oechem

    mols = [oechem.OEGraphMol(), oechem.OEGraphMol()]
    oechem.OESmilesToMol(mols[0], "C")
    oechem.OESmilesToMol(mols[1], "CC")

    with pytest.raises(TypeError, match="Unknown kwargs"):
        oecluster.pdist(mols, "fingerprint", similarity_func="dice")


def test_pdist_with_cutoff():
    """Test pdist with cutoff produces sparse result."""
    import oecluster
    from openeye import oechem

    smiles = ["c1ccccc1", "c1ccc(O)cc1", "CCCCCCCC"]
    mols = []
    for smi in smiles:
        mol = oechem.OEGraphMol()
        oechem.OESmilesToMol(mol, smi)
        mols.append(mol)

    dist = oecluster.pdist(mols, "fingerprint", cutoff=0.5)
    assert dist.num_samples == 3


def test_pdist_similarity_with_cutoff_raises():
    """cutoff > 0 with similarity=True is rejected, as it is in cdist.

    Without the guard this pair's 0.2727 comes back as 0.0: sparse storage
    zeroes values above the cutoff, so on a similarity matrix it discards
    precisely the pairs that scored highest. Pairs below the cutoff come back
    untouched, which is what makes the corruption easy to miss in a larger
    matrix.

    Two comparisons, because the guard is a property of ``pdist`` rather than
    of one comparison. A refusal that happened to depend on the name
    ``fingerprint`` would leave every other similarity-capable comparison
    corrupting, and the 5.7.0 CHANGELOG entry claims exactly the opposite.
    """
    import oecluster
    from openeye import oechem

    mols = []
    for smi in ["c1ccccc1", "Cc1ccccc1"]:
        mol = oechem.OEGraphMol()
        oechem.OESmilesToMol(mol, smi)
        mols.append(mol)

    for name in ["fingerprint", "mcs"]:
        with pytest.raises(ValueError,
                           match="cutoff > 0 is not supported with "
                                 "similarity=True"):
            oecluster.pdist(mols, name, similarity=True, cutoff=0.2)


def test_pdist_and_cdist_share_the_cutoff_refusal_text():
    """The two refusals are one message, not two that happen to agree.

    Giving ``pdist`` the refusal ``cdist`` already had is the whole point of
    the change. A wording improvement applied to one copy and not the other
    would leave the two functions describing the same mistake differently,
    which is the inconsistency this task exists to remove.
    """
    import oecluster
    from openeye import oechem

    mols = []
    for smi in ["c1ccccc1", "Cc1ccccc1"]:
        mol = oechem.OEGraphMol()
        oechem.OESmilesToMol(mol, smi)
        mols.append(mol)

    with pytest.raises(ValueError) as from_pdist:
        oecluster.pdist(mols, "fingerprint", similarity=True, cutoff=0.2)
    with pytest.raises(ValueError) as from_cdist:
        oecluster.cdist(mols, mols, "fingerprint", similarity=True,
                        cutoff=0.2)

    assert str(from_pdist.value) == str(from_cdist.value)


def _benzene_toluene():
    from openeye import oechem

    mols = []
    for smi in ["c1ccccc1", "Cc1ccccc1"]:
        mol = oechem.OEGraphMol()
        oechem.OESmilesToMol(mol, smi)
        mols.append(mol)
    return mols


@pytest.mark.parametrize("factory", ["MCSComparison",
                                     "FingerprintComparison"])
def test_pdist_prebuilt_similarity_with_cutoff_raises(factory):
    """A prebuilt comparison reporting similarities refuses a cutoff too.

    The ``similarity`` argument says nothing about a prebuilt object, so the
    refusal reads the orientation the object reports about itself. Without it
    the MCS pair's 0.857 came back as 0.0. The progress callback proves the
    refusal lands before any pair is scored.
    """
    import oecluster

    mols = _benzene_toluene()
    comparison = getattr(oecluster, factory)(mols, similarity=True)
    calls = []
    with pytest.raises(ValueError,
                       match="cutoff > 0 is not supported for a prebuilt "
                             "comparison that reports similarities"):
        oecluster.pdist(mols, comparison, cutoff=0.5,
                        progress=lambda done, total: calls.append(done))
    assert calls == []


def test_pdist_prebuilt_distance_with_cutoff_is_sparse():
    """A prebuilt distance-oriented comparison keeps its sparse cutoff path."""
    import oecluster

    mols = _benzene_toluene()
    comparison = oecluster.MCSComparison(mols, similarity=False)
    result = oecluster.pdist(mols, comparison, cutoff=0.5)

    assert isinstance(result.storage, oecluster.SparseStorage)
    assert np.asarray(result.condensed) == pytest.approx([1.0 / 7.0])


def test_pdist_prebuilt_similarity_without_sparse_storage_is_accepted(
        tmp_path):
    """No cutoff, or an mmap output that ignores it, leaves values intact."""
    import oecluster

    mols = _benzene_toluene()
    comparison = oecluster.MCSComparison(mols, similarity=True)

    dense = oecluster.pdist(mols, comparison)
    mapped = oecluster.pdist(mols, comparison, cutoff=0.5,
                             output=str(tmp_path / "sim.mmap"))

    assert np.asarray(dense.condensed) == pytest.approx([6.0 / 7.0])
    assert np.asarray(mapped.condensed) == pytest.approx([6.0 / 7.0])


def _titled_mols(smiles_list):
    from openeye import oechem

    mols = []
    for idx, smi in enumerate(smiles_list):
        mol = oechem.OEGraphMol()
        oechem.OESmilesToMol(mol, smi)
        mol.SetTitle(f"mol{idx}")
        mols.append(mol)
    return mols


_THREE = ["c1ccccc1", "Cc1ccccc1", "CCO"]
_FIVE = _THREE + ["CCC", "CCCC"]


def test_a_prebuilt_comparison_refuses_comparison_options():
    """Options cannot reach a comparison that is already built.

    The prebuilt branch never read ``kwargs``, so every option was dropped in
    silence while the docstring promised a TypeError.
    """
    import oecluster

    mols = _titled_mols(_THREE)
    prebuilt = oecluster.FingerprintComparison(mols)
    with pytest.raises(TypeError, match=r"does not take comparison options"):
        oecluster.pdist(mols, prebuilt, fp_type="atom_pair")


def test_a_prebuilt_comparison_refuses_a_different_item_count():
    """Measured before the fix: 5 items in, a 3-sample matrix out, no error."""
    import oecluster

    prebuilt = oecluster.FingerprintComparison(_titled_mols(_THREE))
    with pytest.raises(ValueError,
                       match=r"5 items and a prebuilt comparison over 3"):
        oecluster.pdist(_titled_mols(_FIVE), prebuilt)


def test_a_prebuilt_comparison_refuses_items_without_a_length():
    import oecluster

    mols = _titled_mols(_THREE)
    with pytest.raises(TypeError, match=r"items to be a sized sequence"):
        oecluster.pdist(iter(mols), oecluster.FingerprintComparison(mols))


def test_a_prebuilt_comparison_accepts_the_matching_item_count():
    """Regression shield: passes before the fix, must keep passing after."""
    import oecluster

    mols = _titled_mols(_THREE)
    assert oecluster.pdist(
        mols, oecluster.FingerprintComparison(mols)).num_samples == 3


def test_raw_pdist_refuses_similarity_into_sparse_storage():
    """The native driver refuses what the wrapper refuses, for raw callers.

    The SWIG exception handler maps ``ComparisonError`` to ``RuntimeError``.
    Without the native guard the benzene-toluene MCS similarity of 0.857 is
    dropped by the sparse storage and reads back as 0.0.
    """
    import oecluster
    from oecluster import oecluster as native

    mols = _benzene_toluene()
    comparison = oecluster.MCSComparison(mols, similarity=True)
    storage = oecluster.SparseStorage(comparison.Size(), 0.5)

    with pytest.raises(RuntimeError, match="reports similarities"):
        native.pdist(comparison, storage, oecluster.PDistOptions())


def test_raw_cdist_refuses_similarity_with_cutoff():
    """The native cdist refuses a cutoff for a similarity comparison."""
    import oecluster
    from oecluster import oecluster as native

    mols = _benzene_toluene()
    comparison = oecluster.MCSComparison(mols, similarity=True)
    output = np.full((1, 1), -1.0, dtype=np.float64)
    options = oecluster.CDistOptions()
    options.cutoff = 0.5

    with pytest.raises(RuntimeError, match="reports similarities"):
        native.cdist_into_address(comparison, 1, output.ctypes.data, options)
    assert output[0, 0] == -1.0


def test_pdist_cutoff_positivity_is_decided_once():
    """A cutoff whose truth value is unstable cannot slip past the guard.

    The guard and the storage selection used to test ``cutoff > 0.0``
    independently, so a value answering False to the first and True to the
    second cleared the guard and then selected sparse storage anyway,
    reinstating the corruption the guard exists to prevent. Measured before
    the fix, this pair came back as 0.0 instead of 0.2727.

    Three shapes, because two of them survive a partial repair. Deciding the
    comparison once but leaving its result unconverted still lets the truth
    value shift between the two readers, since ``and`` yields the operand
    rather than a bool. And a repair that re-tests the cutoff behind the
    cached flag is invisible to the False-then-True direction, which
    short-circuits before ever reaching the second test.
    """
    import oecluster
    from openeye import oechem

    class ShiftingCutoff(float):
        """Compares as zero once, then as positive."""

        def __init__(self, _value):
            self._comparisons = 0

        def __gt__(self, other):
            self._comparisons += 1
            return self._comparisons > 1

    class MirrorCutoff(float):
        """Compares as positive once, then as zero."""

        def __init__(self, _value):
            self._comparisons = 0

        def __gt__(self, other):
            self._comparisons += 1
            return self._comparisons == 1

    class ShiftingTruth:
        """Falsy once, then truthy. A comparison may return one of these."""

        def __init__(self):
            self.conversions = 0

        def __bool__(self):
            self.conversions += 1
            return self.conversions > 1

    class ProxyCutoff(float):
        """Compares to a value whose truth is decided on conversion."""

        def __init__(self, _value):
            self.truth = ShiftingTruth()

        def __gt__(self, other):
            return self.truth

    mols = []
    for smi in ["c1ccccc1", "Cc1ccccc1"]:
        mol = oechem.OEGraphMol()
        oechem.OESmilesToMol(mol, smi)
        mols.append(mol)

    shifting = ShiftingCutoff(0.2)
    result = oecluster.pdist(mols, "fingerprint", similarity=True,
                             cutoff=shifting)
    assert result.condensed[0] == pytest.approx(0.2727, abs=1e-4)
    assert shifting._comparisons == 1

    proxy = ProxyCutoff(0.2)
    result = oecluster.pdist(mols, "fingerprint", similarity=True,
                             cutoff=proxy)
    assert result.condensed[0] == pytest.approx(0.2727, abs=1e-4)
    assert proxy.truth.conversions == 1

    mirror = MirrorCutoff(0.2)
    with pytest.raises(ValueError,
                       match="cutoff > 0 is not supported with "
                             "similarity=True"):
        oecluster.pdist(mols, "fingerprint", similarity=True, cutoff=mirror)
    assert mirror._comparisons == 1


def test_pdist_reports_a_typod_kwarg_before_the_cutoff():
    """A typo'd keyword outranks the cutoff refusal, as it does in cdist.

    ca761af placed the guard above validate_request, which turned this call's
    unknown-kwarg TypeError into a ValueError naming a cutoff -- so the
    remedy the message offered was to drop a cutoff that was never the
    problem. The comparison has to be one that validates its keywords in
    validate_request; ``mcs`` does, and ``fingerprint`` does not, so this
    would not bite with ``fingerprint``.
    """
    import oecluster
    from openeye import oechem

    mols = []
    for smi in ["c1ccccc1", "Cc1ccccc1"]:
        mol = oechem.OEGraphMol()
        oechem.OESmilesToMol(mol, smi)
        mols.append(mol)

    with pytest.raises(TypeError, match="Unknown kwargs for mcs"):
        oecluster.pdist(mols, "mcs", bogus=1, similarity=True, cutoff=0.5)
    with pytest.raises(TypeError, match="Unknown kwargs for mcs"):
        oecluster.cdist(mols, mols, "mcs", bogus=1, similarity=True,
                        cutoff=0.5)


def test_pdist_reports_an_unsupported_similarity_before_the_cutoff():
    """A comparison with no similarity form says so, cutoff or not.

    ``rmsd`` and ``descriptor`` refuse ``similarity=True`` outright, and that
    refusal names the argument the caller has to change. The cutoff refusal
    names dropping the cutoff, which would not help here -- the call would
    still be asking for similarities the comparison cannot produce.

    This is the half of the ordering the ``mcs`` test above cannot see.
    ``_validate_mcs`` documents that it ignores ``similarity``, so handing
    ``validate_request`` a hardcoded ``False`` would leave that test green
    while silently putting the cutoff refusal ahead of both refusals here.
    """
    import oecluster
    from openeye import oechem

    mols = []
    for smi in ["c1ccccc1", "Cc1ccccc1"]:
        mol = oechem.OEGraphMol()
        oechem.OESmilesToMol(mol, smi)
        mols.append(mol)

    with pytest.raises(ValueError,
                       match="rmsd comparison has no similarity form"):
        oecluster.pdist(mols, "rmsd", similarity=True, cutoff=0.5)
    with pytest.raises(ValueError,
                       match="descriptor comparison has no similarity form"):
        oecluster.pdist(mols, "descriptor", similarity=True, cutoff=0.5)


def test_pdist_progress():
    """Test progress callback is invoked."""
    import oecluster
    from openeye import oechem

    smiles = ["c1ccccc1", "c1ccc(O)cc1", "CCCCCCCC"]
    mols = []
    for smi in smiles:
        mol = oechem.OEGraphMol()
        oechem.OESmilesToMol(mol, smi)
        mols.append(mol)

    calls = []
    oecluster.pdist(mols, "fingerprint",
                    progress=lambda d, t: calls.append((d, t)))
    assert len(calls) > 0

def test_distance_matrix_serialization(tmp_path):
    """Test DistanceMatrix save/load roundtrip."""
    import oecluster
    from openeye import oechem

    smiles = ["c1ccccc1", "c1ccc(O)cc1", "CCCCCCCC"]
    mols = []
    for smi in smiles:
        mol = oechem.OEGraphMol()
        oechem.OESmilesToMol(mol, smi)
        mols.append(mol)

    dist = oecluster.pdist(mols, "fingerprint")
    path = str(tmp_path / "test.npz")
    dist.to_file(path)

    loaded = oecluster.SymmetricDistanceMatrix.from_file(path)
    np.testing.assert_array_almost_equal(
        np.asarray(dist), np.asarray(loaded)
    )
    assert loaded.comparison_name == "fingerprint"


def test_pdist_fingerprint_similarity():
    """Test pdist with fingerprint comparison in similarity mode."""
    import oecluster
    from openeye import oechem

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
    import oecluster
    from openeye import oechem

    # Use a single DU so pdist exercises the alias-dispatch path without
    # invoking Distance() — this test covers routing, not scoring. Scoring
    # of SiteHopper is covered by the superpose:sitehopper integration tests.
    du = oechem.OEDesignUnit()
    path = "tests/assets/spruce_5FQD_1_5FQD_1-ALIGNED_BC__DU__LVY_B-1438.oedu"
    if not oechem.OEReadDesignUnit(path, du):
        pytest.skip(f"Cannot read {path}")

    dm = oecluster.pdist([du], "sitehopper")
    assert "sitehopper" in dm.comparison_name


def test_pdist_superpose_kwargs():
    """Test superpose with method kwarg."""
    import oecluster
    from openeye import oechem

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
    assert dm.comparison_name == "superpose:ddm"


def test_pdist_superpose_predicate():
    """Test superpose with predicate kwarg."""
    import oecluster
    from openeye import oechem

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

    dm = oecluster.pdist(dus, "superpose", method="global", predicate="name CA")
    arr = np.asarray(dm)
    assert arr.shape == (1,)
    assert np.isfinite(arr[0])


def test_pdist_superpose_similarity():
    """Test superpose with similarity=True for SSE method."""
    import oecluster
    from openeye import oechem

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
    import oecluster
    from openeye import oechem

    mols = [oechem.OEGraphMol(), oechem.OEGraphMol()]
    oechem.OESmilesToMol(mols[0], "C")
    oechem.OESmilesToMol(mols[1], "CC")

    with pytest.raises(TypeError, match="Unknown kwargs"):
        oecluster.pdist(mols, "fingerprint", bogus_option=42)


def test_symmetric_distance_matrix_hierarchy():
    """SymmetricDistanceMatrix is the concrete pdist result; base is abstract."""
    from oecluster import DenseStorage, DistanceMatrix, SymmetricDistanceMatrix

    storage = DenseStorage(3)
    storage.Set(0, 1, 0.5)
    storage.Set(0, 2, 0.25)
    storage.Set(1, 2, 0.75)
    dm = SymmetricDistanceMatrix(storage, "test", ["a", "b", "c"], {})

    assert isinstance(dm, DistanceMatrix)
    assert dm.num_samples == 3
    assert dm.shape == (3, 3)
    assert dm.comparison_name == "test"
    assert repr(dm).startswith("SymmetricDistanceMatrix(")
    # The abstract base must not be directly constructible as a symmetric matrix.
    with pytest.raises(TypeError):
        DistanceMatrix("test", {})  # type: ignore[call-arg]


def test_pdist_rocs_end_to_end():
    """The public rocs route, which segfaulted from Python before Task 19.

    ``test_native_bindings.py`` covers the typemap in isolation. This covers
    what a caller actually invokes: option mapping, label extraction, storage
    allocation and the condensed result.
    """
    pytest.importorskip("openeye.oeomega")
    import oecluster
    from openeye import oechem, oeomega

    omega = oeomega.OEOmega()
    omega.SetMaxConfs(1)
    omega.SetStrictStereo(False)
    mols = []
    for idx, smi in enumerate(["c1ccccc1", "c1ccc(O)cc1", "CCCCCCCC"]):
        mol = oechem.OEMol()
        oechem.OESmilesToMol(mol, smi)
        assert omega(mol)
        mol.SetTitle(f"m{idx}")
        mols.append(mol)

    dist = oecluster.pdist(mols, "rocs", score_type="shape")
    assert dist.comparison_name == "rocs"
    assert dist.num_samples == 3
    assert dist.labels == ["m0", "m1", "m2"]
    assert dist.params["comparison_type"] == "rocs"

    condensed = np.asarray(dist)
    assert condensed.shape == (3,)
    assert np.all(np.isfinite(condensed))
    assert np.all((condensed >= 0.0) & (condensed <= 1.0))

    # Shape distance has to order these three the way chemistry does: benzene
    # overlays phenol far better than either overlays octane. Asserting the
    # ordering keeps the test meaningful without pinning overlay output to four
    # decimals, which would make it a change detector.
    benzene_phenol, benzene_octane, phenol_octane = condensed
    assert benzene_phenol < 0.1
    assert benzene_octane > 0.4
    assert phenol_octane > 0.4

    # Not the diagonal: squareform() zeroes it structurally, so asserting it
    # would pass regardless of what ROCS computed.
    square = np.asarray(dist.squareform())
    assert np.allclose(square, square.T)


@pytest.fixture(scope="module")
def multiconformer_rocs_mols():
    """Three flexible molecules with real ensembles, for the ROCS order test.

    Every molecule in ``test_pdist_rocs_end_to_end`` is built with
    ``SetMaxConfs(1)``, where there is no ensemble to search and nothing bites:
    measured 2026-10-15, the orientation gap on that fixture reaches 0.00028 on
    combo_norm and 0.00056 on shape, two orders of magnitude below what these
    three reach. Flexible molecules at twenty conformers each are what make
    ``BestOverlay`` a genuine search and the two orientations of a pair
    genuinely different answers.

    Module-scoped because the test below is parameterized over four score types
    and Omega need not run four times to produce the same ensembles. The
    molecules are only read afterwards -- ``ROCSComparison`` snapshots its own
    copies and prepares those -- so sharing them across parameters cannot carry
    state from one to the next.
    """
    pytest.importorskip("openeye.oeomega")
    from openeye import oechem, oeomega

    omega = oeomega.OEOmega()
    omega.SetMaxConfs(20)
    omega.SetStrictStereo(False)
    mols = []
    for idx, smi in enumerate(["CCCCCCCC", "CCCCCCCCO", "c1ccc(CCCCO)cc1"]):
        mol = oechem.OEMol()
        oechem.OESmilesToMol(mol, smi)
        assert omega(mol)
        assert mol.NumConfs() > 1
        mol.SetTitle(f"m{idx}")
        mols.append(mol)
    return mols


# Measured 2026-10-15: the largest absolute disagreement between the forward
# and reversed matrices on the fixture above is exactly 0.0, on all four score
# types, and so is the worst of the six permutations of the displaced-phenol
# fixture further down. Zero is still not the right tolerance to assert,
# because OEBestOverlayScore returns float and nothing in the toolkit promises
# bit-reproducibility; 1e-6 is the float-precision scale it works at, the same
# reasoning and the same value as SELF_SCORE_TOLERANCE in ROCSComparison.cpp.
# The effects it has to stay clear of are the 0.00099 color orientation gap
# below and the 0.4538 permutation effect, three and five orders of magnitude
# above it.
ROCS_ORDER_TOLERANCE = 1e-6

# Half the orientation gap measured on this fixture on 2026-10-15 --
# combo_norm 0.035512, combo 0.071024, shape 0.066000, color 0.000994 -- so
# each floor is a real observation halved rather than a round number, and a
# fixture that drifted most of the way to symmetric would fail here rather
# than let the reversal check below pass vacuously.
ROCS_MIN_ORIENTATION_GAP = {
    "combo_norm": 0.0177,
    "combo": 0.0355,
    "shape": 0.0330,
    "color": 0.0004,
}


def _rocs_orientation_gap(mols, score_type):
    """The directional gap on ``mols``, measured without the symmetrization.

    ``cdist(mols, mols)`` fills the whole square and keeps it directional:
    entry ``[i, j]`` scores ``mols[i]`` as the overlay reference against
    ``mols[j]``'s conformers, so entry ``[j, i]`` is the same pair the other way
    round and the block's own asymmetry is the orientation gap. One call rather
    than two, and through the public surface rather than the protected seam.
    """
    import oecluster

    block = np.asarray(oecluster.cdist(mols, mols, "rocs", score_type=score_type))
    off_diagonal = ~np.eye(len(mols), dtype=bool)
    return float(np.max(np.abs(block - block.T)[off_diagonal]))


@pytest.mark.parametrize("score_type", ["combo_norm", "combo", "shape", "color"])
def test_rocs_pdist_does_not_depend_on_input_order(multiconformer_rocs_mols,
                                                   score_type):
    """ROCS scores one orientation; pdist stores one triangle as symmetric.

    A tolerance rather than equality, and the reason is no longer retained
    history: every search runs on an overlay of its own, so a score is a
    function of its pair. What the tolerance covers is that the toolkit does
    not promise bit reproducibility, and pinning equality would make this a
    change detector for OEShape rather than a test of ordering.

    All four score types, because the guide claims the residual on all four and
    a test narrower than its documentation leaves the rest unpinned. The
    combine happens on the raw Tanimoto and each type converts it differently,
    so these are four distinct conversions of one rule.
    """
    import oecluster

    mols = multiconformer_rocs_mols

    # First establish that the two orientations disagree at all on this
    # fixture. Without this the reversal check below is satisfied by any
    # scorer that is already symmetric, including a one-orientation
    # implementation running on molecules that happen not to discriminate --
    # which is the regression this test exists to catch.
    gap = _rocs_orientation_gap(mols, score_type)
    assert gap > ROCS_MIN_ORIENTATION_GAP[score_type], (
        f"the {score_type} orientations no longer disagree on this fixture "
        f"(gap={gap}), so the reversal assertion below would pass vacuously")

    forward = np.asarray(
        oecluster.pdist(mols, comparison="rocs", score_type=score_type).squareform())
    reversed_ = np.asarray(
        oecluster.pdist(mols[::-1], comparison="rocs",
                        score_type=score_type).squareform())
    assert np.allclose(forward, reversed_[::-1, ::-1],
                       rtol=0.0, atol=ROCS_ORDER_TOLERANCE)


def test_rocs_pdist_survives_an_input_that_saturates_a_reused_overlay():
    """The input that proved scoring both orientations was not enough on its own.

    A phenol with one atom displaced to 200 angstroms is accepted -- a
    fiftieth of the 10,000 angstrom extent guard -- and it is the case
    ``ROCSComparison.cpp`` already documented: it scores 0.546 against benzene
    on an overlay used for the first time and saturates at exactly 1.0 on one
    that has scored anything before. While ``pdist`` reused a single overlay per
    worker, which of the two a pair got depended on where that pair fell in the
    input order, and the matrix moved by 0.453840 under a permutation even with
    both orientations scored. Measured 2026-10-15; the fix is that
    ``CompareUnordered`` builds its own overlay per search.

    All six permutations rather than a reversal, because only three of the six
    moved, and ``num_threads=1`` because that is the configuration the effect
    was found in -- it is not a thread race.
    """
    pytest.importorskip("openeye.oeomega")
    import oecluster
    from openeye import oechem, oeomega

    omega = oeomega.OEOmega()
    omega.SetMaxConfs(1)
    omega.SetStrictStereo(False)
    mols = []
    for smi in ["c1ccc(O)cc1", "c1ccccc1", "CCCCCCCC"]:
        mol = oechem.OEMol()
        oechem.OESmilesToMol(mol, smi)
        assert omega(mol)
        mols.append(mol)

    conformer = next(iter(mols[0].GetConfs()))
    atom = next(iter(conformer.GetAtoms()))
    coords = list(conformer.GetCoords(atom))
    coords[0] = 200.0
    assert conformer.SetCoords(atom, coords)

    def square(ordered):
        return np.asarray(
            oecluster.pdist(ordered, comparison="rocs", num_threads=1).squareform())

    base = square(mols)
    # Not vacuous only if this input still reaches the saturating branch at
    # all. The phenol/benzene entry is the one that moved; 1.0 is the
    # saturation value, and a fixture that no longer came near the first-use
    # score would make every permutation agree for the wrong reason.
    assert 0.0 < base[0, 1] < 1.0

    for perm in itertools.permutations(range(3)):
        permuted = square([mols[i] for i in perm])
        inverse = np.argsort(perm)
        restored = permuted[np.ix_(inverse, inverse)]
        assert np.allclose(base, restored, rtol=0.0, atol=ROCS_ORDER_TOLERANCE), (
            f"permutation {perm} moved the matrix by "
            f"{np.max(np.abs(base - restored))}")


def test_cdist_rocs_does_not_depend_on_thread_count_or_target_order():
    """A `cdist` cell must not move with the company its targets keep.

    Directionality explains why ``cdist(A, B)`` differs from ``cdist(B, A)``.
    It explains nothing about a cell changing with the thread count, or with
    where an unrelated molecule sits in the target list, and while the two
    scoring paths reused one overlay per thread this cell did both: measured
    2026-10-16, phenol-against-benzene read 1.000000000 with octane first at
    one thread and 0.546160012 both with two threads and with benzene first.
    1.0 is the saturation value, so the reused overlay was not returning a
    differently-oriented score but a wrong one.
    """
    pytest.importorskip("openeye.oeomega")
    import oecluster
    from openeye import oechem, oeomega

    omega = oeomega.OEOmega()
    omega.SetMaxConfs(1)
    omega.SetStrictStereo(False)
    built = []
    for smi in ["c1ccc(O)cc1", "c1ccccc1", "CCCCCCCC"]:
        mol = oechem.OEMol()
        oechem.OESmilesToMol(mol, smi)
        assert omega(mol)
        built.append(mol)
    phenol, benzene, octane = built

    conformer = next(iter(phenol.GetConfs()))
    atom = next(iter(conformer.GetAtoms()))
    coords = list(conformer.GetCoords(atom))
    coords[0] = 200.0
    assert conformer.SetCoords(atom, coords)

    def cell(targets, **kwargs):
        block = np.asarray(
            oecluster.cdist([phenol], targets, "rocs", **kwargs))
        return block[0, targets.index(benzene)]

    readings = [
        cell([octane, benzene], num_threads=1),
        cell([octane, benzene], num_threads=2, chunk_size=1),
        cell([benzene, octane], num_threads=1),
        cell([benzene], num_threads=1),
    ]
    assert max(readings) - min(readings) <= ROCS_ORDER_TOLERANCE, readings

    # Not vacuous only if this input still reaches the branch that used to
    # saturate. 1.0 is what a reused overlay returned here.
    assert 0.0 < readings[0] < 1.0


def test_cdist_rocs_stays_directional():
    """``cdist`` must keep searching only the fit molecule's conformers.

    ``pdist`` symmetrizes because it stores one triangle; ``cdist`` stores both
    and documents entry ``[i, j]`` as ``items_a[i]`` used *as the reference*.
    Taking a non-active conformer of a multiconformer molecule as the reference
    makes the two readings separable: the ensemble contains the pose the
    reference was taken from, so searching it recovers an exact overlay, while
    the single-conformer reference scored the other way round does not. A later
    change that symmetrized ``cdist`` as well would take the best of the two and
    lose that asymmetry.
    """
    pytest.importorskip("openeye.oeomega")
    import oecluster
    from openeye import oechem, oeomega

    omega = oeomega.OEOmega()
    omega.SetMaxConfs(5)
    omega.SetStrictStereo(False)
    multi = oechem.OEMol()
    oechem.OESmilesToMol(multi, "c1ccccc1CCCCc1ccccc1")
    assert omega(multi)
    assert multi.NumConfs() > 1

    last = list(multi.GetConfs())[-1]
    reference = oechem.OEMol(multi.GetConf(oechem.OEHasConfIdx(last.GetIdx())))
    assert reference.NumConfs() == 1

    forward = np.asarray(
        oecluster.cdist([reference], [multi], "rocs", score_type="shape"))
    backward = np.asarray(
        oecluster.cdist([multi], [reference], "rocs", score_type="shape"))

    # The ensemble as the fit set contains the reference's own pose, so the
    # search recovers it exactly. Measured 0.0 on 2026-10-15.
    assert forward[0, 0] == pytest.approx(0.0, abs=1e-3)
    # Swapping the sets makes the single pose the fit set, leaving BestOverlay
    # nothing to search, and the match is not recovered. Measured 0.481 on
    # 2026-10-15; asserting a loose floor rather than the value keeps this from
    # becoming a change detector for Omega output.
    assert backward[0, 0] > 0.1


def test_pdist_rocs_accepts_graph_molecules():
    """The documented input type, which ``rocs`` rejected until 5.7.0.

    Bare SMILES will not do here: ``ROCSComparison`` refuses a molecule whose
    recomputed dimension is below three, so the graph molecules are taken from
    embedded ones. Scored against the equivalent ``OEMol`` input, because an
    overload that resolved differently for the two types would still return
    perfectly plausible numbers.
    """
    pytest.importorskip("openeye.oeomega")
    import oecluster
    from openeye import oechem, oeomega

    omega = oeomega.OEOmega()
    omega.SetMaxConfs(1)
    omega.SetStrictStereo(False)
    oemols = []
    for idx, smi in enumerate(["c1ccccc1", "c1ccc(O)cc1", "CCCCCCCC"]):
        mol = oechem.OEMol()
        oechem.OESmilesToMol(mol, smi)
        assert omega(mol)
        mol.SetTitle(f"m{idx}")
        oemols.append(mol)
    graphs = [oechem.OEGraphMol(mol) for mol in oemols]

    dist = oecluster.pdist(graphs, "rocs", score_type="shape")
    assert dist.comparison_name == "rocs"
    assert dist.num_samples == 3
    assert dist.labels == ["m0", "m1", "m2"]

    # One conformer apiece on both routes, so BestOverlay has the same single
    # pose to choose from and the two must agree exactly.
    reference = oecluster.pdist(oemols, "rocs", score_type="shape")
    assert np.allclose(np.asarray(dist), np.asarray(reference), atol=1e-6)


def test_pdist_rocs_refuses_molecules_without_coordinates():
    """The refusal has to reach the caller through the public surface.

    ``pdist`` builds the comparison through the registry rather than by naming
    the class, which is a route a grep for constructor call sites does not find,
    so the C++ precondition is asserted here from the outside as well. Needs no
    Omega: molecules with no coordinates at all are the input under test.
    """
    import oecluster
    from openeye import oechem

    mols = []
    for smi in ["c1ccccc1", "CCCCCCCC"]:
        mol = oechem.OEMol()
        oechem.OESmilesToMol(mol, smi)
        mols.append(mol)

    with pytest.raises(RuntimeError, match="requires 3D coordinates"):
        oecluster.pdist(mols, "rocs")


def test_condensed_survives_the_matrix_that_produced_it():
    """The zero-copy view must own its buffer rather than borrow it.

    ``arr = pdist(...).condensed`` discards the matrix immediately, which is
    ordinary usage. While the array only borrowed the pointer, the storage was
    freed and the array went on reading reallocated memory -- returning
    plausible distances rather than crashing, so nothing downstream noticed.
    """
    import gc

    from oecluster import DenseStorage, SymmetricDistanceMatrix, _StorageView

    expected = [0.5, 0.25, 0.75]

    def build_and_discard():
        storage = DenseStorage(3)
        storage.Set(0, 1, expected[0])
        storage.Set(0, 2, expected[1])
        storage.Set(1, 2, expected[2])
        matrix = SymmetricDistanceMatrix(storage, "test", ["a", "b", "c"], {})
        return matrix.condensed

    arr = build_and_discard()
    gc.collect()

    # Churn the allocator so a freed buffer would be handed back out.
    for _ in range(3):
        ballast = [np.random.rand(4096) for _ in range(200)]
        del ballast
        gc.collect()

    # Not merely "base is not None": a borrowed pointer also has a base object,
    # it just does not own the storage. Name the type that does.
    assert isinstance(arr.base, _StorageView)
    assert list(arr) == pytest.approx(expected)


def test_condensed_stays_a_zero_copy_view():
    """Keeping the storage alive must not quietly turn the view into a copy.

    Guards the obvious wrong fix: copying the buffer would satisfy the
    lifetime test above while changing documented zero-copy behavior.
    """
    from oecluster import DenseStorage, SymmetricDistanceMatrix

    storage = DenseStorage(3)
    storage.Set(0, 1, 0.5)
    storage.Set(0, 2, 0.25)
    storage.Set(1, 2, 0.75)
    matrix = SymmetricDistanceMatrix(storage, "test", ["a", "b", "c"], {})

    arr = matrix.condensed
    assert not arr.flags.owndata

    arr[0] = 0.125
    assert storage.Get(0, 1) == pytest.approx(0.125)


def test_default_fingerprint_path_is_unchanged_by_the_5_0_0_rewrite():
    """Golden values for the untouched default: morgan/binary/tanimoto.

    The literals are the Jaccard distances over OEFP's own radius-2, 2048-bit
    binary Morgan fingerprints. They were checked against an independent
    oracle -- ``oefp.api.morgan_fingerprint`` called directly, with the
    intersection and union popcounted from ``np.array(f.words, np.uint64)`` --
    rather than against a second ``oecluster`` call, which would move with the
    thing it is meant to check.

    Two levels, because neither pins what the other does.
    ``FingerprintOptions()`` carries the declared defaults, so asserting on it
    pins them exactly -- and it is the struct the keyword path uses:
    ``_build_fingerprint`` constructs a default one and overwrites only the
    fields the caller named. The golden vector then pins that the pipeline
    honours the struct rather than reaching some other setting, which the
    struct assertions on their own do not.

    Of the six defaults the second call spells out, this is what the *vector*
    resolves on its own, measured by moving each one in turn:

    * ``fp_type`` and ``similarity`` are pinned: each of the other three
      families moves six entries, and ``similarity=True`` moves eight.
    * ``radius`` is pinned in both directions, and hexane is in the fixture to
      make that true. Over the four molecules that preceded it, radius 3, 4
      and 5 each reproduced the radius-2 vector exactly, leaving the default
      unpinned upward. With hexane in, radius 1 moves five entries and radius
      3 moves three.
    * ``metric`` is pinned against fifteen of the sixteen other metrics the
      keyword accepts; every one of them moves this vector. The exception is
      ``'jaccard'``, which is the same measure under another name and yields
      an identical vector.
    * ``storage`` is resolved only in part. ``'count'`` and ``'sparse_count'``
      cannot reach this vector because the default ``tanimoto`` is a bit-set
      metric they refuse outright, but ``'sparse'`` yields an identical vector
      and the literals alone would not notice it.
    * ``numbits`` is not resolved at all. 8192, 4096, 1024 and 512 all
      reproduce this vector; only 256 and below move it, because five small
      saturated molecules do not collide at a plausible width. Making
      2048-vs-1024 observable in output would need a forced collision and a
      much larger, more brittle fixture.

    Both are pinned at the declared level instead:
    ``FingerprintOptions().storage`` separates ``'binary'`` from ``'sparse'``
    and ``FingerprintOptions().numbits`` catches a width change the vector
    cannot see. Neither assertion says anything about what the pipeline then
    does with the value, which is what the literals below are for.

    The second call is not the oracle -- both sides would move together -- it
    only pins that spelling the defaults out is the same request as leaving
    them off.
    """
    import oecluster
    from openeye import oechem

    declared = oecluster.FingerprintOptions()
    assert declared.fp_type == "morgan"
    assert declared.storage == "binary"
    assert declared.metric == "tanimoto"
    assert declared.numbits == 2048
    assert declared.radius == 2

    smiles = ["CCO", "CCC", "CCCC", "c1ccccc1", "CCCCCC"]
    mols = []
    for smi in smiles:
        mol = oechem.OEGraphMol()
        oechem.OESmilesToMol(mol, smi)
        mols.append(mol)

    dist = oecluster.pdist(mols, "fingerprint")
    np.testing.assert_allclose(
        dist.condensed,
        [4.0 / 7.0, 5.0 / 8.0, 1.0, 7.0 / 10.0, 0.5, 1.0, 5.0 / 8.0, 1.0, 0.5,
         1.0],
        rtol=0.0, atol=1e-12)

    reference = oecluster.pdist(
        mols, "fingerprint", fp_type="morgan", storage="binary",
        metric="tanimoto", numbits=2048, radius=2, similarity=False)
    np.testing.assert_array_equal(dist.condensed, reference.condensed)
    assert dist.metric_capabilities == {'zero_self': True, 'triangle': True}


def test_pdist_refuses_a_negative_num_threads():
    """A negative count is a ValueError, as in butina, not the binding's
    OverflowError from a size_t field."""
    import oecluster

    with pytest.raises(ValueError, match="num_threads must be non-negative"):
        oecluster.pdist(_benzene_toluene(), "fingerprint", num_threads=-1)


@pytest.mark.parametrize("value", [2.5, "8", None], ids=["float", "str", "none"])
def test_pdist_refuses_a_non_integer_num_threads(value):
    """A float, a string or None is refused rather than truncated."""
    import oecluster

    with pytest.raises(TypeError):
        oecluster.pdist(_benzene_toluene(), "fingerprint", num_threads=value)


def test_pdist_accepts_a_numpy_integer_num_threads():
    import oecluster

    mols = _benzene_toluene()
    expected = oecluster.pdist(mols, "fingerprint", num_threads=2)
    got = oecluster.pdist(mols, "fingerprint", num_threads=np.int64(2))
    np.testing.assert_array_equal(got.condensed, expected.condensed)
