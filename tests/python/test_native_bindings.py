"""The SWIG surface the Python layer is built on.

These assertions are deliberately about names and shapes rather than behavior:
they fail loudly when an interface-file edit silently drops a symbol, which is
otherwise only visible as an AttributeError deep inside a builder.
"""

import math

import pytest


@pytest.fixture
def native():
    from oecluster import oecluster as _oecluster

    return _oecluster


def test_capability_enum_is_exposed(native):
    assert native.Capability_Unknown != native.Capability_No
    assert native.Capability_No != native.Capability_Yes
    assert native.Capability_Yes != native.Capability_Unknown


def test_data_integrity_enum_is_exposed(native):
    values = {
        native.DataIntegrity_Complete,
        native.DataIntegrity_NaNPresent,
        native.DataIntegrity_SubsetScored,
    }
    assert len(values) == 3


def test_gate_facts_defaults(native):
    facts = native.GateFacts()
    assert facts.is_distance == native.Capability_Unknown
    assert facts.zero_self == native.Capability_Unknown
    assert facts.triangle == native.Capability_Unknown
    assert facts.data_integrity == native.DataIntegrity_Complete


def _benzene_series():
    from openeye import oechem

    mols = []
    for smi in ("c1ccccc1", "c1ccc(O)cc1", "CCCCCCCC", "c1ccncc1"):
        mol = oechem.OEGraphMol()
        oechem.OESmilesToMol(mol, smi)
        mols.append(mol)
    return mols


def test_fingerprint_comparison_reports_facts(native):
    comparison = native.FingerprintComparison(_benzene_series(), native.FingerprintOptions())
    facts = comparison.Facts()
    assert facts.is_distance == native.Capability_Yes
    assert facts.zero_self == native.Capability_Yes
    assert facts.triangle == native.Capability_Yes
    assert facts.data_integrity == native.DataIntegrity_Complete


def test_descriptor_comparison_is_exposed(native):
    comparison = native.DescriptorComparison(_benzene_series(), native.DescriptorOptions())
    assert comparison.ComparisonName() == "descriptor"
    assert len(list(comparison.Columns())) > 0
    assert len(list(comparison.Variances())) == len(list(comparison.Columns()))
    assert comparison.Facts().data_integrity == native.DataIntegrity_Complete


def test_descriptor_excluded_indices_is_exposed(native):
    excluded = native.descriptor_excluded_indices(_benzene_series(), native.DescriptorOptions())
    assert list(excluded) == []


def test_descriptor_statistics_is_exposed(native):
    stats = native.descriptor_statistics(_benzene_series(), native.DescriptorStatisticsOptions())
    assert stats.num_rows == 4
    assert len(list(stats.columns)) == len(list(stats.variance))
    assert len(list(stats.dropped_columns)) == len(list(stats.dropped_reasons))
    assert list(stats.inverse_covariance) == []


def test_descriptor_statistics_options_accept_lists(native):
    options = native.DescriptorStatisticsOptions()
    options.sources = native.StringVector(["openeye"])
    options.inverse_covariance = True
    stats = native.descriptor_statistics(_benzene_series(), options)
    k = len(list(stats.columns))
    assert len(list(stats.inverse_covariance)) == k * k


def test_rmsd_options_are_exposed(native):
    assert hasattr(native, 'RMSDComparison')
    assert hasattr(native, 'RMSDOptions')
    options = native.RMSDOptions()
    assert options.overlay is False
    assert options.automorph is True
    assert options.heavy_only is True
    options.overlay = True
    options.automorph = False
    options.heavy_only = False
    assert options.overlay is True
    assert options.automorph is False
    assert options.heavy_only is False


def test_rmsd_comparison_from_molecules(native):
    from openeye import oechem

    mols = []
    for shift in range(3):
        mol = oechem.OEMol()
        oechem.OESmilesToMol(mol, "c1ccccc1")
        oechem.OEGenerate2DCoordinates(mol)
        for atom in mol.GetAtoms():
            x, y, z = mol.GetCoords(atom)
            mol.SetCoords(atom, (x + shift, y, z))
        mols.append(mol)

    comparison = native.RMSDComparison(mols, native.RMSDOptions())
    assert comparison.ComparisonName() == "rmsd"
    assert comparison.Size() == 3
    assert comparison.Compare(0, 0) == pytest.approx(0.0, abs=1e-6)
    assert comparison.Facts().zero_self == native.Capability_Yes


def _conformer_series():
    """Three benzenes with 3D coordinates, rigidly translated apart."""
    pytest.importorskip("openeye.oeomega")
    from openeye import oechem, oeomega

    omega = oeomega.OEOmega()
    omega.SetMaxConfs(1)
    omega.SetStrictStereo(False)

    mols = []
    for shift in range(3):
        mol = oechem.OEMol()
        oechem.OESmilesToMol(mol, "c1ccccc1")
        assert omega(mol)
        for atom in mol.GetAtoms():
            x, y, z = mol.GetCoords(atom)
            mol.SetCoords(atom, (x + shift, y, z))
        mols.append(mol)
    return mols


def test_rocs_comparison_accepts_molecules(native):
    """The typemap's other consumer. This path segfaulted before Task 19."""
    mols = _conformer_series()
    comparison = native.ROCSComparison(mols, native.ROCSOptions())
    assert comparison.ComparisonName() == "rocs"
    assert comparison.Size() == 3
    # Only that a finite score comes back. This is a typemap test, not a scoring
    # test: what it exists to prove is that molecules survive the crossing into
    # C++. The scores themselves are pinned in the C++ suite, and asserting a
    # value here would couple a binding test to the overlay numerics.
    value = comparison.Compare(0, 1)
    assert math.isfinite(value)
    assert 0.0 <= value <= 2.0


def test_rocs_shape_self_distance_is_zero(native):
    """Shape-only distance has a zero diagonal, asserted apart from combo.

    Shape never depended on the color preparation that combo needs, so keeping
    it as its own case means a regression in that preparation cannot make this
    assertion fail, and a failure here points at the overlay itself.
    """
    options = native.ROCSOptions()
    options.score_type = native.ROCSScoreType_Shape
    comparison = native.ROCSComparison(_conformer_series(), options)
    assert comparison.Compare(0, 0) == pytest.approx(0.0, abs=1e-5)


def test_typemap_preserves_conformers(native):
    """The shared_ptr typemap's copy must keep every conformer, not just one.

    ``ROCSComparison::Compare`` overlays the fit molecule with
    ``OEOverlay::BestOverlay``, which searches all of its conformers. A
    reference taken from a later conformer therefore scores a perfect overlay
    only if the copy kept them; against a single-conformer fit the same
    reference scores measurably worse. The second loop asserts that gap, so a
    fixture that stopped discriminating would fail rather than pass vacuously.
    """
    pytest.importorskip("openeye.oeomega")
    from openeye import oechem, oeomega

    omega = oeomega.OEOmega()
    omega.SetMaxConfs(3)
    omega.SetStrictStereo(False)
    multi = oechem.OEMol()
    oechem.OESmilesToMol(multi, "c1ccccc1CCCCc1ccccc1")
    assert omega(multi)
    assert multi.NumConfs() > 1

    def conformer_as_mol(conf):
        return oechem.OEMol(multi.GetConf(oechem.OEHasConfIdx(conf.GetIdx())))

    options = native.ROCSOptions()
    options.score_type = native.ROCSScoreType_Shape
    references = [conformer_as_mol(conf) for conf in multi.GetConfs()]

    for index, reference in enumerate(references):
        comparison = native.ROCSComparison([reference, multi], options)
        assert comparison.Compare(0, 1) == pytest.approx(0.0, abs=1e-3), (
            f"conformer {index} was not reachable in the copied molecule")

    single = references[0]
    gaps = [native.ROCSComparison([reference, single], options).Compare(0, 1)
            for reference in references[1:]]
    assert max(gaps) > 0.1, (
        "the conformers generated here are too similar for this test to "
        f"distinguish a dropped conformer from a kept one: gaps={gaps}")


def test_sar_coherence_is_bound(native):
    """The labels overload takes a plain list and returns the per-cluster table."""
    coherence = native.sar_coherence([0, 0, 1, 1], [1.0, 1.0, 5.0, 5.0])

    assert coherence.num_samples == 4
    assert coherence.num_scored == 4
    assert coherence.num_clusters == 2
    assert coherence.eta_squared == 1.0
    assert len(coherence.clusters) == 2
    assert coherence.clusters[0].label == 0
    assert coherence.clusters[0].mean_activity == 1.0


def test_sar_coherence_options_carry_the_excluded_default(native):
    options = native.SARCoherenceOptions()

    assert options.noise_handling == native.NoiseHandling_Excluded


def test_activity_landscape_is_bound(native):
    storage = native.DenseStorage(3)
    storage.Set(0, 1, 0.5)
    storage.Set(0, 2, 0.25)
    storage.Set(1, 2, 0.125)

    landscape = native.activity_landscape(storage, [0.0, 1.0, 3.0])

    assert landscape.num_samples == 3
    assert landscape.num_pairs_scored == 3
    assert landscape.max_sali == 16.0
    assert landscape.mean_sali == 10.0
    assert landscape.num_cliffs == 2


def test_activity_landscape_options_carry_the_published_defaults(native):
    options = native.ActivityLandscapeOptions()

    assert options.distance_threshold == 0.30
    assert options.activity_threshold == 1.0
    assert options.rmodi_delta == 0.625
    assert options.num_threads == 0


def test_modelability_is_bound(native):
    storage = native.DenseStorage(4)
    storage.Set(0, 1, 0.5)
    storage.Set(0, 2, 0.1)
    storage.Set(0, 3, 0.6)
    storage.Set(1, 2, 0.7)
    storage.Set(1, 3, 0.8)
    storage.Set(2, 3, 0.2)

    model = native.modelability(storage, ["A", "A", "B", "B"])

    assert model.num_classes == 2
    assert model.modi == 0.5
    assert len(model.classes) == 2
    assert model.classes[0].label == "A"
    assert model.classes[0].num_members == 2
    assert model.classes[0].fraction_same_class == 0.5


def test_the_row_vectors_are_wrapped_rather_than_opaque(native):
    """Reads a field off the first row of each member vector.

    A member vector wrapped as an opaque pointer still len()s and indexes from
    Python but hands back a SwigPyObject with no fields, so the existence of
    the two template names proves nothing on its own. Touching a row's field is
    what distinguishes a wrapped row from an opaque one, and it is the shape
    the Pythonic layer reads.
    """
    assert hasattr(native, "ClusterActivityVector")
    assert hasattr(native, "ClassConcordanceVector")

    coherence = native.sar_coherence([0, 0, 1, 1], [1.0, 1.0, 5.0, 5.0])
    assert coherence.clusters[0].num_scored == 2

    storage = native.DenseStorage(2)
    storage.Set(0, 1, 0.5)
    model = native.modelability(storage, ["A", "B"])
    assert model.classes[0].label == "A"


def test_the_new_entry_points_release_the_gil():
    """All three sweeps are O(n^2) or O(N) over native data and must not hold
    the interpreter while they run.

    Asserted against the interface file rather than at runtime. A timing test
    cannot separate "released the GIL" from "finished quickly" without a
    fixture large enough to hold the lock for a measurable stretch, and the
    only way to build one here is an O(n^2) Python loop that costs more suite
    time than the assertion buys. The interface file is the sole input SWIG
    reads for this, so a missing invocation here is a missing release in the
    generated wrapper.
    """
    import pathlib

    interface = pathlib.Path(__file__).resolve().parents[2] / "swig" / "oecluster.i"
    text = interface.read_text(encoding="utf-8")

    for name in ("sar_coherence", "activity_landscape", "modelability"):
        assert f"OECLUSTER_GIL_EXCEPTION(OECluster::{name}, {name})" in text
