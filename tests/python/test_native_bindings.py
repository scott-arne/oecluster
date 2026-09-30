"""The SWIG surface the Python layer is built on.

These assertions are deliberately about names and shapes rather than behavior:
they fail loudly when an interface-file edit silently drops a symbol, which is
otherwise only visible as an AttributeError deep inside a builder. A few of them
do pin exact values, where reading the number back is the only way to show that
a typemap carried real data across rather than a zero-initialized struct; those
values are either exact in binary or, as with ``distance_threshold``, the
nearest double to the same decimal literal the header writes, so the bare == is
deliberate.
"""

import math
import pathlib
import re

import pytest

# Both comment forms SWIG honours. Stripped before the interface file is
# searched, so a directive that is only present inside a comment does not
# satisfy a check that it is present.
_SWIG_COMMENT = re.compile(r"/\*.*?\*/|//[^\n]*", re.DOTALL)


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


def _graph_conformer_series():
    """``_conformer_series`` viewed as OEGraphMol, so each carries one pose."""
    from openeye import oechem

    return [oechem.OEGraphMol(mol) for mol in _conformer_series()]


def test_rocs_comparison_accepts_graph_molecules(native):
    """The OEMolBase overload, added because the docs all build OEGraphMol.

    Before it existed this list raised a bare SWIG "Wrong number or type of
    arguments" TypeError, naming no argument and suggesting no remedy.
    """
    comparison = native.ROCSComparison(_graph_conformer_series(), native.ROCSOptions())
    assert comparison.ComparisonName() == "rocs"
    assert comparison.Size() == 3
    # As in test_rocs_comparison_accepts_molecules: this is an overload test,
    # not a scoring test, so it asks only that the molecules survived the
    # crossing into C++.
    value = comparison.Compare(0, 1)
    assert math.isfinite(value)
    assert 0.0 <= value <= 2.0


def _flexible_multiconformer(native):
    """A molecule whose conformers differ enough to be told apart by shape.

    The ensemble is what the tests below use to tell a preserved
    multiconformer copy from a collapsed one, so "differ enough" is asserted
    rather than assumed: if Omega ever generates a tighter ensemble for this
    input, those tests would agree on a number for the wrong reason instead of
    failing. The same self-check, against the same threshold, guards the
    inlined fixture in ``test_typemap_preserves_conformers``.
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

    options = native.ROCSOptions()
    options.score_type = native.ROCSScoreType_Shape
    references = [oechem.OEMol(multi.GetConf(oechem.OEHasConfIdx(conf.GetIdx())))
                  for conf in multi.GetConfs()]
    gaps = [native.ROCSComparison([reference, references[0]], options).Compare(0, 1)
            for reference in references[1:]]
    assert max(gaps) > 0.1, (
        "the conformers generated here are too similar for the tests using "
        f"this fixture to distinguish a dropped ensemble from a kept one: "
        f"gaps={gaps}")
    return multi


def _non_active_reference(multi):
    """An ``OEMol`` of a conformer that is *not* the active one.

    The active conformer is the only pose an ``OEMolBase`` view of ``multi``
    keeps, so a reference drawn from any other one is reachable through the
    ensemble and unreachable through the collapsed copy. That asymmetry is what
    makes the conformer tests below able to tell the two apart.
    """
    from openeye import oechem

    active_idx = multi.GetActive().GetIdx()
    others = [conf.GetIdx() for conf in multi.GetConfs() if conf.GetIdx() != active_idx]
    assert others, "fixture produced no conformer other than the active one"
    return oechem.OEMol(multi.GetConf(oechem.OEHasConfIdx(others[-1])))


def test_oemol_input_still_binds_the_conformer_preserving_overload(native):
    """The hazard the OEMolBase overload introduces, pinned as a property.

    ``ROCSComparison`` has two constructors an ``OEMol`` list satisfies both
    of: ``vector<shared_ptr<OEMol>>`` and ``vector<OEMolBase*>``, the latter
    matching because an ``OEMol`` is an ``OEMolBase``. Only the first preserves
    conformers; the second copies through an ``OEMolBase&`` view, which for a
    multiconformer ``OEMol`` is its active conformer alone. Which one binds is
    decided by typemap precedence in ``swig/oecluster.i`` -- explicitly, since
    an equal precedence leaves it to a tie-break that was measured to rank them
    one way for one-argument calls and the opposite way for two-argument ones.

    Two arguments here on purpose. That is the form ``_comparisons.py`` builds,
    so it is the form ``pdist(..., "rocs")`` reaches, and it is the form that
    was broken while the one-argument form stayed correct.

    The assertion is a property, not a literal: the ensemble and the same
    molecule collapsed to its active conformer must score *differently*. A
    literal would pin the number while the property was what broke -- when the
    ensemble collapsed, these two scores became equal to the last digit.
    """
    from openeye import oechem

    multi = _flexible_multiconformer(native)
    options = native.ROCSOptions()
    options.score_type = native.ROCSScoreType_Shape
    reference = _non_active_reference(multi)

    # OEGraphMol is the active-conformer view, so round-tripping through it
    # produces exactly the loss the wrong overload would cause. Rebuilt as an
    # OEMol so that this list binds the same constructor the kept case does,
    # leaving the ensemble as the only difference between the two calls.
    collapsed_mol = oechem.OEMol(oechem.OEGraphMol(multi))
    assert collapsed_mol.NumConfs() == 1
    assert multi.NumConfs() > collapsed_mol.NumConfs()

    kept = native.ROCSComparison([reference, multi], options).Compare(0, 1)
    collapsed = native.ROCSComparison([reference, collapsed_mol], options).Compare(0, 1)

    # The property. Dropping the ensemble makes the first call compute the
    # second, and this margin is what goes to zero when it does.
    assert collapsed - kept > 0.1, (
        "the multiconformer OEMol scored as its active conformer alone, so the "
        f"ensemble was dropped crossing into C++: kept={kept}, "
        f"collapsed={collapsed}")

    # Corroboration, not the point: a reachable conformer overlays its own
    # reference exactly, which says the ensemble arrived intact rather than
    # merely arriving different.
    assert kept == pytest.approx(0.0, abs=1e-3)


def test_mixed_molecule_lists_resolve_asymmetrically_by_their_first_element(native):
    """Both orderings of a list mixing OEGraphMol with OEMol.

    Both typechecks sample element 0 only, so the first element alone decides
    which overload the whole list binds -- and the two orderings are therefore
    not symmetric.

    ``[OEGraphMol, OEMol]`` fails the strict check at element 0 and binds the
    OEMolBase overload, which accepts every element. This is the one ordering
    where a caller can still lose an ensemble: an OEMol later in the list is
    converted through an OEMolBase view. It is accepted rather than refused
    because refusing it would mean rejecting the OEGraphMol lists this overload
    exists to admit.

    ``[OEMol, OEGraphMol]`` passes the strict check at element 0, and that
    typemap validates *every* element rather than sampling one, so the
    OEGraphMol at index 1 is refused by name. Pinned on both the one- and
    two-argument forms: precedence used to order those two groups oppositely,
    and this refusal survived on one form while vanishing on the other.
    """
    from openeye import oechem

    multi = _flexible_multiconformer(native)
    options = native.ROCSOptions()
    options.score_type = native.ROCSScoreType_Shape
    reference = _non_active_reference(multi)
    graph_reference = oechem.OEGraphMol(reference)

    # Accepted, and the ensemble is lost -- element 0 put the whole list on the
    # OEMolBase overload. Asserted against the collapsed score rather than a
    # literal, which is what makes the conformer loss the claim.
    collapsed_mol = oechem.OEMol(oechem.OEGraphMol(multi))
    mixed = native.ROCSComparison([graph_reference, multi], options).Compare(0, 1)
    collapsed = native.ROCSComparison([reference, collapsed_mol], options).Compare(0, 1)
    assert mixed == pytest.approx(collapsed, abs=1e-6)

    # The reverse ordering is refused, on both argument-count forms.
    reverse = [reference, oechem.OEGraphMol(multi)]
    with pytest.raises(TypeError, match="List item is not an OEMol object"):
        native.ROCSComparison(reverse, options)
    with pytest.raises(TypeError, match="List item is not an OEMol object"):
        native.ROCSComparison(reverse)


def test_the_molbase_typecheck_is_ranked_behind_the_strict_one():
    """The precedence the three tests above rest on, asserted from the source.

    Those three all sit behind ``importorskip("openeye.oeomega")``, so on a
    machine without that license the entire guard against a precedence
    regression disappears. This one has no license gate.

    ``SWIG_TYPECHECK_POINTER`` is 0, and the strict
    ``vector<shared_ptr<OEMol>>`` typecheck carries it. Giving the permissive
    ``vector<OEMolBase*>`` typecheck the same value makes the two an exact tie,
    which SWIG then resolves by a rule it does not specify -- measured to rank
    them one way for one-argument calls and the opposite way for two-argument
    ones, which is how every ``pdist(..., "rocs")`` call came to drop
    conformers while the one-argument form stayed correct. Any value above
    zero breaks the tie in the strict typemap's favour, since lower precedence
    is examined first. Comments are stripped before the search so that the long
    note explaining this, which quotes both spellings, cannot satisfy it.

    Same limit as ``test_the_new_entry_points_release_the_gil``: this describes
    the interface, not the built artifact, so a stale ``_oecluster.so`` passes
    it. What it covers is the window before a rebuild -- where a tidy-up
    substituting the symbolic constant back would land. Where an Omega license
    exists, the behavioural tests above remain the real guard.
    """
    interface = pathlib.Path(__file__).resolve().parents[2] / "swig" / "oecluster.i"
    text = _SWIG_COMMENT.sub("", interface.read_text(encoding="utf-8"))

    def precedence_of(cpp_type):
        found = re.findall(
            r"%typemap\(typecheck,\s*precedence=([^)]+)\)\s*" + re.escape(cpp_type),
            text)
        assert len(found) == 1, (
            f"expected exactly one typecheck typemap for {cpp_type}, found {found}")
        return found[0].strip()

    strict = precedence_of("const std::vector<std::shared_ptr<OEChem::OEMol>>&")
    permissive = precedence_of("const std::vector<OEChem::OEMolBase*>&")

    assert strict == "SWIG_TYPECHECK_POINTER", (
        "this test compares the permissive typecheck against the strict one's "
        f"precedence; the strict one now reads {strict!r}")
    assert permissive != strict, (
        "equal precedence restores the tie whose unspecified resolution dropped "
        "conformers from every pdist(..., 'rocs') call")
    assert permissive.isdigit() and int(permissive) > 0, (
        "the permissive typecheck needs a numeric precedence above zero so the "
        f"strict one is examined first; found {permissive!r}")


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

    Asserted against the interface file rather than at runtime.

    Two ways a directive can be present but inert are both caught. Position:
    ``%exception`` binds only to declarations SWIG parses after the invocation,
    so an invocation sitting below the ``%include`` that declares the function
    reaches nothing. Commenting out: a ``//`` or a ``/* */`` around the
    invocation has the same effect. Either one removes the release from every
    generated wrapper and SWIG says nothing about it, so comments are stripped
    before the search and both anchors are measured on the stripped text.

    The limit of reading the source is that it describes the interface, not the
    artifact: a stale ``_oecluster.so`` whose wrappers predate the directives
    passes, as does a correctly placed directive naming a function that does
    not exist, which SWIG also accepts silently. Timing a call from a second
    thread would cover the first, but cannot separate "released the GIL" from
    "finished quickly" without a fixture large enough to hold the lock for a
    measurable stretch; the source is the better trade.
    """
    interface = pathlib.Path(__file__).resolve().parents[2] / "swig" / "oecluster.i"
    text = _SWIG_COMMENT.sub("", interface.read_text(encoding="utf-8"))

    include = '%include "oecluster/clustering/SARCoherence.h"'
    assert text.count(include) == 1, "the position check needs an unambiguous anchor"
    include_at = text.index(include)

    for name in ("sar_coherence", "activity_landscape", "modelability"):
        invocation = f"OECLUSTER_GIL_EXCEPTION(OECluster::{name}, {name})"
        assert invocation in text
        assert text.index(invocation) < include_at, (
            f"{invocation} must precede the %include that declares {name}, or "
            "SWIG applies the exception handler to nothing")


# Five points on a line at these coordinates; the distance is the gap. Farthest
# first from item 0 picks 4 (8 away), then 2 (3 from its nearest pick).
_LINE = (0.0, 1.0, 3.0, 7.0, 8.0)


def _line_storage(native, positions=_LINE):
    storage = native.DenseStorage(len(positions))
    for i in range(len(positions)):
        for j in range(i + 1, len(positions)):
            storage.Set(i, j, abs(positions[i] - positions[j]))
    return storage


def test_maxmin_select_is_exposed(native):
    options = native.MaxMinOptions()
    assert options.count == 0
    assert math.isnan(options.threshold)
    assert options.seed_mode == native.MaxMinSeed_Index
    assert options.seed == 0
    assert options.num_threads == 0
    assert options.chunk_size == 256

    options.count = 3
    result = native.maxmin_select(_line_storage(native), options)
    assert list(result.indices) == [0, 4, 2]
    assert math.isnan(result.pick_distances[0])
    assert list(result.pick_distances)[1:] == [8.0, 3.0]
    assert result.stop == native.MaxMinStop_Count


def test_maxmin_select_takes_an_initial_selection(native):
    options = native.MaxMinOptions()
    options.count = 3
    initial = native.SizeTVector()
    initial.push_back(1)
    initial.push_back(3)
    options.initial = initial
    result = native.maxmin_select(_line_storage(native), options)
    assert list(result.indices) == [1, 3, 2]
    assert list(result.pick_distances)[2] == 2.0


def test_maxmin_select_runs_on_a_comparison(native):
    options = native.MaxMinOptions()
    options.count = 2
    options.seed_mode = native.MaxMinSeed_Farthest
    comparison = native.FingerprintComparison(_benzene_series(),
                                              native.FingerprintOptions())
    result = native.maxmin_select(comparison, options)
    assert len(result.indices) == 2
    assert result.stop == native.MaxMinStop_Count


def test_the_stop_and_seed_enums_are_exposed(native):
    assert {native.MaxMinStop_Count, native.MaxMinStop_Threshold,
            native.MaxMinStop_Exhausted} == {0, 1, 2}
    assert {native.MaxMinSeed_Index, native.MaxMinSeed_Medoid,
            native.MaxMinSeed_Farthest} == {0, 1, 2}
    assert {native.CirclesMethod_MaxMin,
            native.CirclesMethod_Sequential} == {0, 1}


def test_circles_is_exposed(native):
    options = native.CirclesOptions()
    assert options.method == native.CirclesMethod_MaxMin
    assert options.num_threads == 0
    assert options.chunk_size == 256

    packing = native.circles(_line_storage(native), 2.5, options)
    assert packing.count == 3
    assert list(packing.members) == [0, 4, 2]
    assert packing.threshold == 2.5
    assert packing.method == native.CirclesMethod_MaxMin

    options.method = native.CirclesMethod_Sequential
    packing = native.circles(_line_storage(native), 2.5, options)
    assert list(packing.members) == [0, 2, 3]
    assert packing.method == native.CirclesMethod_Sequential


def test_maxmin_select_without_a_stop_condition_is_a_runtime_error(native):
    with pytest.raises(RuntimeError,
                       match="requires a count, a threshold, or both"):
        native.maxmin_select(_line_storage(native), native.MaxMinOptions())


def test_a_similarity_comparison_is_a_runtime_error(native):
    fingerprint_options = native.FingerprintOptions()
    fingerprint_options.similarity = True
    comparison = native.FingerprintComparison(_benzene_series(),
                                              fingerprint_options)
    options = native.MaxMinOptions()
    options.count = 2
    with pytest.raises(RuntimeError, match="requires distances"):
        native.maxmin_select(comparison, options)
    with pytest.raises(RuntimeError, match="requires distances"):
        native.circles(comparison, 0.5, native.CirclesOptions())


def test_the_medoid_seed_on_a_comparison_is_a_runtime_error(native):
    comparison = native.FingerprintComparison(_benzene_series(),
                                              native.FingerprintOptions())
    options = native.MaxMinOptions()
    options.count = 2
    options.seed_mode = native.MaxMinSeed_Medoid
    with pytest.raises(RuntimeError,
                       match="Medoid requires a distance matrix"):
        native.maxmin_select(comparison, options)


def test_a_nan_circles_threshold_is_a_runtime_error(native):
    with pytest.raises(RuntimeError, match="finite and non-negative"):
        native.circles(_line_storage(native), math.nan,
                       native.CirclesOptions())


def test_a_non_finite_distance_read_is_a_runtime_error(native):
    """The Python matrix path refuses this up front through the gate, so the
    native refusal is only reachable through the raw binding."""
    storage = _line_storage(native)
    storage.Set(0, 4, math.nan)
    options = native.MaxMinOptions()
    options.count = 2
    with pytest.raises(RuntimeError,
                       match="non-finite distance between items 0 and 4"):
        native.maxmin_select(storage, options)


def test_the_diversity_entry_points_release_the_gil():
    """Both entry points fold O(N) rows per pick over native data, and the
    comparison overloads run for as long as the comparisons do.

    Asserted against the interface file for the reasons given in
    test_the_new_entry_points_release_the_gil: a directive below the
    ``%include`` that declares its function, or one inside a comment, is
    silently inert.
    """
    interface = pathlib.Path(__file__).resolve().parents[2] / "swig" / "oecluster.i"
    text = _SWIG_COMMENT.sub("", interface.read_text(encoding="utf-8"))

    include = '%include "oecluster/clustering/DiversitySelection.h"'
    assert text.count(include) == 1, "the position check needs an unambiguous anchor"
    include_at = text.index(include)

    for name in ("maxmin_select", "circles"):
        invocation = f"OECLUSTER_GIL_EXCEPTION(OECluster::{name}, {name})"
        assert invocation in text
        assert text.index(invocation) < include_at, (
            f"{invocation} must precede the %include that declares {name}, or "
            "SWIG applies the exception handler to nothing")


def test_sphere_exclusion_releases_the_gil():
    """sphere_exclusion reads O(N*k) native rows or builds the O(N^2)
    threshold graph, and its comparison overload runs for as long as the
    comparisons do.

    Asserted against the interface file for the reasons given in
    test_the_new_entry_points_release_the_gil.
    """
    interface = pathlib.Path(__file__).resolve().parents[2] / "swig" / "oecluster.i"
    text = _SWIG_COMMENT.sub("", interface.read_text(encoding="utf-8"))

    include = '%include "oecluster/clustering/SphereExclusion.h"'
    assert text.count(include) == 1, "the position check needs an unambiguous anchor"
    invocation = "OECLUSTER_GIL_EXCEPTION(OECluster::sphere_exclusion, sphere_exclusion)"
    assert invocation in text
    assert text.index(invocation) < text.index(include), (
        f"{invocation} must precede the %include that declares sphere_exclusion, "
        "or SWIG applies the exception handler to nothing")
