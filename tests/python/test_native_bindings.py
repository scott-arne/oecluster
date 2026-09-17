"""The SWIG surface the Python layer is built on.

These assertions are deliberately about names and shapes rather than behavior:
they fail loudly when an interface-file edit silently drops a symbol, which is
otherwise only visible as an AttributeError deep inside a builder. A few of them
do pin exact values, where reading the number back is the only way to show that
a typemap carried real data across rather than a zero-initialized struct; those
values are chosen to be exact in binary so the bare == is deliberate.
"""

import math
import pathlib

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


def test_modelability_options_carry_the_hardware_concurrency_default(native):
    options = native.ModelabilityOptions()

    assert options.num_threads == 0


def test_sar_coherence_accepts_a_clustering_result(native):
    """The other sar_coherence overload: a derived result, bound as a base ref.

    butina_cluster returns ButinaResult; the overload takes
    ``const ClusteringResult&``. Nothing else in this file crosses that
    inheritance edge, and it is the form the Pythonic layer calls.

    The two-cluster split is forced by the fixture -- 0.1 within each pair,
    0.9 across -- so the effect size is the degenerate 1.0 rather than
    something the numerics could drift. The row order is the interesting part:
    Butina labels this input [1, 1, 0, 0], so a table ordered by first
    appearance among the scored samples starts at label 1, not at label 0.
    """
    storage = native.DenseStorage(4)
    storage.Set(0, 1, 0.1)
    storage.Set(0, 2, 0.9)
    storage.Set(0, 3, 0.9)
    storage.Set(1, 2, 0.9)
    storage.Set(1, 3, 0.9)
    storage.Set(2, 3, 0.1)

    butina = native.ButinaOptions()
    butina.distance_threshold = 0.5
    result = native.butina_cluster(storage, butina)
    assert list(result.Labels()) == [1, 1, 0, 0]

    coherence = native.sar_coherence(result, [1.0, 1.0, 5.0, 5.0])

    assert coherence.num_samples == 4
    assert coherence.num_scored == 4
    assert coherence.num_clusters == 2
    assert coherence.eta_squared == 1.0
    assert len(coherence.clusters) == 2
    assert coherence.clusters[0].label == 1
    assert coherence.clusters[0].mean_activity == 1.0
    assert coherence.clusters[1].label == 0
    assert coherence.clusters[1].mean_activity == 5.0


def test_the_options_object_changes_the_answer_it_governs(native):
    """The three-argument wrappers, with an option set to move the result.

    Every other call in this file takes the defaulted two-argument form, so the
    options wrapper is only reached here. Passing a *default-valued* options
    object would not test much: if the object failed to cross, the native
    defaults would apply and produce the same numbers, so the call could only
    fail by raising. Each option below is therefore set to a non-default value
    that changes a specific count, and the two forms are asserted to differ in
    exactly that count.

    sar_coherence: the fixture carries a -1, so noise_handling has something to
    do. The default drops it; Singletons promotes it, raising num_scored and
    num_clusters. A failure means either the options object did not cross or
    noise handling stopped being applied. Note this fixture has one noise
    sample, so it separates Excluded from the other two readings but not
    Singletons from Grouped.

    activity_landscape: raising distance_threshold past the 0.5 pair makes a
    third pair near, so num_cliffs rises. A failure means the threshold did not
    cross or is no longer consulted.
    """
    noisy = [0, 0, 1, -1]
    activity = [1.0, 1.0, 5.0, 9.0]

    dropped = native.sar_coherence(noisy, activity)
    assert dropped.num_samples == 4
    assert dropped.num_scored == 3
    assert dropped.num_clusters == 2
    assert [row.label for row in dropped.clusters] == [0, 1]

    singleton_options = native.SARCoherenceOptions()
    singleton_options.noise_handling = native.NoiseHandling_Singletons
    promoted = native.sar_coherence(noisy, activity, singleton_options)

    assert promoted.num_samples == 4
    assert promoted.num_scored == 4
    assert promoted.num_clusters == 3
    assert [row.label for row in promoted.clusters] == [0, 1, -1]

    def landscape_storage():
        storage = native.DenseStorage(3)
        storage.Set(0, 1, 0.5)
        storage.Set(0, 2, 0.25)
        storage.Set(1, 2, 0.125)
        return storage

    near = [0.0, 1.0, 3.0]
    narrow = native.activity_landscape(landscape_storage(), near)
    assert narrow.num_cliffs == 2

    wide_options = native.ActivityLandscapeOptions()
    wide_options.distance_threshold = 0.6
    wide = native.activity_landscape(landscape_storage(), near, wide_options)

    assert wide.num_cliffs == 3
    assert wide.cliff_density == 1.0


def test_the_modelability_options_overload_is_reachable(native):
    """Reachability only, because ModelabilityOptions has nothing that can move
    a result.

    Its single field is num_threads (SARCoherence.h:228), documented at :227
    and :240-243 as leaving the result identical bit for bit. There is no value
    that would make the two forms differ, so equality is the correct
    expectation and the only defect this can catch is the three-argument
    wrapper raising -- a bad typemap on the options parameter, or the overload
    not being generated at all. It cannot show that the object's contents were
    read.
    """
    storage = native.DenseStorage(4)
    storage.Set(0, 1, 0.5)
    storage.Set(0, 2, 0.1)
    storage.Set(0, 3, 0.6)
    storage.Set(1, 2, 0.7)
    storage.Set(1, 3, 0.8)
    storage.Set(2, 3, 0.2)
    classes = ["A", "A", "B", "B"]

    options = native.ModelabilityOptions()
    options.num_threads = 1
    threaded = native.modelability(storage, classes, options)

    assert threaded.modi == native.modelability(storage, classes).modi
    assert threaded.modi == 0.5


def test_the_row_vectors_are_wrapped_rather_than_opaque(native):
    """The member vectors resolve to the wrapped vector types, not to a pointer.

    Drop the %template for one of these and its member stops being a sequence
    altogether: SARCoherence.clusters comes back as a bare SwigPyObject that
    implements neither __len__ nor __getitem__, so len() and [0] raise
    TypeError before any row exists. That failure surfaces first in
    test_sar_coherence_is_bound, as an opaque TypeError on the len() line.

    The two module-level names below prove the vector types were instantiated
    somewhere, which is what names that cause -- but on their own they say
    nothing about how the members are typed. The isinstance checks are the part
    that ties the member to the template, and they are asserted nowhere else in
    the suite.
    """
    assert hasattr(native, "ClusterActivityVector")
    assert hasattr(native, "ClassConcordanceVector")

    coherence = native.sar_coherence([0, 0, 1, 1], [1.0, 1.0, 5.0, 5.0])
    assert isinstance(coherence.clusters, native.ClusterActivityVector)

    storage = native.DenseStorage(2)
    storage.Set(0, 1, 0.5)
    model = native.modelability(storage, ["A", "B"])
    assert isinstance(model.classes, native.ClassConcordanceVector)


def test_the_new_entry_points_release_the_gil():
    """All three sweeps are O(n^2) or O(N) over native data and must not hold
    the interpreter while they run.

    Asserted against the interface file rather than at runtime. A timing test
    cannot separate "released the GIL" from "finished quickly" without a
    fixture large enough to hold the lock for a measurable stretch, and the
    only way to build one here is an O(n^2) Python loop that costs more suite
    time than the assertion buys.

    Position is asserted as well as presence, and the position is the half that
    matters. ``%exception`` binds only to declarations SWIG parses after the
    invocation, so an invocation sitting below the ``%include`` that declares
    the function reaches nothing: the release disappears from every generated
    wrapper, SWIG says nothing about it, and a presence-only check stays green.

    What this cannot see: it reads the interface source, not the built
    extension, so it passes against a stale ``_oecluster.so`` whose wrappers
    predate the invocations. Only a rebuild rules that out.
    """
    interface = pathlib.Path(__file__).resolve().parents[2] / "swig" / "oecluster.i"
    text = interface.read_text(encoding="utf-8")

    include = '%include "oecluster/clustering/SARCoherence.h"'
    assert text.count(include) == 1, "the position check needs an unambiguous anchor"
    include_at = text.index(include)

    for name in ("sar_coherence", "activity_landscape", "modelability"):
        invocation = f"OECLUSTER_GIL_EXCEPTION(OECluster::{name}, {name})"
        assert invocation in text
        assert text.index(invocation) < include_at, (
            f"{invocation} must precede the %include that declares {name}, or "
            "SWIG applies the exception handler to nothing")
