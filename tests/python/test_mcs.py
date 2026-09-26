"""MCS comparison: registration, options and gate facts."""

import oecluster
import pytest
from oecluster import oecluster as native
from openeye import oechem

BENZENE = "c1ccccc1"
TOLUENE = "Cc1ccccc1"
CYCLOHEXANE = "C1CCCCC1"
MORPHINE = "CN1CC[C@]23c4c5ccc(O)c4O[C@H]2[C@@H](O)C=C[C@H]3[C@H]1C5"
PENICILLIN_G = "CC1(C)S[C@@H]2[C@H](NC(=O)Cc3ccccc3)C(=O)N2[C@H]1C(=O)O"
SUCROSE = ("OC[C@H]1O[C@@](CO)(O[C@H]2O[C@H](CO)[C@@H](O)[C@H](O)[C@H]2O)"
           "[C@@H](O)[C@@H]1O")
MACROLIDE = ("CC[C@H]1OC(=O)[C@H](C)[C@@H](O)[C@H](C)[C@@H](O)[C@](C)(O)C"
             "[C@@H](C)C(=O)[C@H](C)[C@@H](O)[C@]1(C)O")
DIBORANE = "[BH2]1[H][BH2][H]1"
SINGLE_BRIDGE = "[BH2][H][BH2]"
HEXANE = "CCCCCC"
ETHANE = "CC"
ETHENE = "C=C"
NEUTRAL_HYDRIDE = "[H][Li]"
CHARGED_HYDRIDE = "[H-][Li+]"
DIHYDROGEN = "[H][H]"


def _mol(smiles, title="mol"):
    """Parse a SMILES. MCS is topological, so no embedding step is needed."""
    mol = oechem.OEMol()
    oechem.OESmilesToMol(mol, smiles)
    mol.SetTitle(title)
    return mol


def _pair(first, second):
    return [_mol(first, "first"), _mol(second, "second")]


def _graph_mol(smiles, title="mol"):
    """Parse a SMILES into an ``OEGraphMol``, the type every doc example builds."""
    mol = oechem.OEGraphMol()
    oechem.OESmilesToMol(mol, smiles)
    mol.SetTitle(title)
    return mol


def _graph_pair(first, second):
    return [_graph_mol(first, "first"), _graph_mol(second, "second")]


def _multiconformer(smiles, title="mol"):
    """Build a three-conformer molecule without an embedding step.

    MCS never reads coordinates, so what the conformers contain does not
    matter -- only that the molecule carries more than one.
    """
    mol = _mol(smiles, title)
    oechem.OEAddExplicitHydrogens(mol)
    coords = oechem.OEFloatArray(3 * mol.GetMaxAtomIdx())
    mol.NewConf(coords)
    mol.NewConf(coords)
    return mol


def test_pdist_resolves_the_mcs_name():
    dm = oecluster.pdist(_pair(BENZENE, TOLUENE), "mcs")
    assert dm.num_samples == 2
    assert dm.condensed[0] == pytest.approx(0.142857, abs=1e-6)


def test_cdist_resolves_the_mcs_name():
    # cdist returns a CrossDistanceMatrix: .shape plus a numpy .matrix.
    result = oecluster.cdist([_mol(BENZENE)], [_mol(TOLUENE)], "mcs")
    assert result.shape == (1, 1)
    assert result.matrix[0, 0] == pytest.approx(0.142857, abs=1e-6)


def test_similarity_is_one_minus_distance():
    distances = oecluster.pdist(_pair(BENZENE, TOLUENE), "mcs")
    similarities = oecluster.pdist(_pair(BENZENE, TOLUENE), "mcs", similarity=True)
    assert similarities.condensed[0] == pytest.approx(
        1.0 - distances.condensed[0], abs=1e-9)


def test_similarity_mode_flips_both_orientation_facts():
    # The Python end of the similarity flag path. It fails if _build_mcs never
    # assigns opts.similarity, because the C++ facts are read off that flag.
    dm = oecluster.pdist(_pair(BENZENE, TOLUENE), "mcs", similarity=True)
    assert dm.is_distance is False
    assert dm.metric_capabilities['zero_self'] is False


def test_facts_report_unknown_triangle_and_complete_integrity():
    dm = oecluster.pdist(_pair(BENZENE, TOLUENE), "mcs")
    assert dm.is_distance is True
    assert dm.metric_capabilities == {'zero_self': True, 'triangle': "unknown"}
    assert dm.data_integrity == "complete"


def test_unknown_kwarg_names_the_mcs_comparison():
    with pytest.raises(TypeError, match="mcs comparison"):
        oecluster.pdist(_pair(BENZENE, TOLUENE), "mcs", overlay=True)


def test_bad_search_mode_is_refused_by_name():
    with pytest.raises(ValueError, match="Unknown MCS search mode"):
        oecluster.pdist(_pair(BENZENE, TOLUENE), "mcs", search_mode="thorough")


def test_bad_match_level_is_refused_by_name():
    with pytest.raises(ValueError, match="Unknown MCS match level"):
        oecluster.pdist(_pair(BENZENE, TOLUENE), "mcs", match_level="strict")


def test_max_matches_zero_is_refused():
    with pytest.raises(RuntimeError, match="max_matches"):
        oecluster.pdist(_pair(BENZENE, TOLUENE), "mcs", max_matches=0)


def test_pdist_and_cdist_agree_on_shared_pairs():
    mols = [_mol(BENZENE, "b"), _mol(TOLUENE, "t"), _mol(MORPHINE, "m")]
    # SymmetricDistanceMatrix has no per-pair accessor, so compare squareforms.
    square = oecluster.pdist(mols, "mcs").squareform()
    cross = oecluster.cdist(mols, mols, "mcs")
    for i in range(len(mols)):
        for j in range(len(mols)):
            assert cross.matrix[i, j] == pytest.approx(square[i, j], abs=1e-9)


def test_match_level_is_observable():
    # Benzene against cyclohexane: loose leaves bonds unconstrained so the whole
    # ring matches, while the default refuses aromatic against single.
    loose = oecluster.pdist(_pair(BENZENE, CYCLOHEXANE), "mcs", match_level="loose")
    assert loose.condensed[0] == pytest.approx(0.0, abs=1e-6)
    default = oecluster.pdist(_pair(BENZENE, CYCLOHEXANE), "mcs")
    assert default.condensed[0] == pytest.approx(1.0, abs=1e-6)

    # Benzene against toluene: exact adds hydrogen count and degree, dropping
    # the match from 6 bonds to 4.
    exact = oecluster.pdist(_pair(BENZENE, TOLUENE), "mcs", match_level="exact")
    assert exact.condensed[0] == pytest.approx(0.555556, abs=1e-6)


def test_search_mode_is_observable():
    # Sucrose against a macrolide fragment: exhaustive finds 17 matched bonds
    # where approximate finds 15. At 6.4 ms it is the cheapest measured pair
    # that discriminates the two modes.
    approximate = oecluster.pdist(_pair(SUCROSE, MACROLIDE), "mcs")
    assert approximate.condensed[0] == pytest.approx(0.605263, abs=1e-6)
    exhaustive = oecluster.pdist(_pair(SUCROSE, MACROLIDE), "mcs",
                                 search_mode="exhaustive")
    assert exhaustive.condensed[0] == pytest.approx(0.527778, abs=1e-6)


def test_max_matches_is_observable():
    # Morphine against penicillin G: a budget of one match cuts the count from
    # 11 bonds to 2.
    bounded = oecluster.pdist(_pair(MORPHINE, PENICILLIN_G), "mcs", max_matches=1)
    assert bounded.condensed[0] == pytest.approx(0.958333, abs=1e-6)
    unbounded = oecluster.pdist(_pair(MORPHINE, PENICILLIN_G), "mcs")
    assert unbounded.condensed[0] == pytest.approx(0.717949, abs=1e-6)


def test_mcs_comparison_is_exported():
    assert "MCSComparison" in oecluster.__all__


def test_top_level_wrapper_constructs_with_a_match_level():
    comparison = oecluster.MCSComparison(_pair(BENZENE, TOLUENE), match_level="exact")
    assert comparison.ComparisonName() == "mcs"
    assert comparison.Compare(0, 1) == pytest.approx(0.555556, abs=1e-6)


def test_top_level_wrapper_honours_similarity():
    # The wrapper does not route through _build_mcs, so it carries its own
    # opts.similarity assignment. Forgetting it is silent: the comparison would
    # return distances and claim to be one. The pdist tests above exercise
    # _build_mcs and would pass regardless.
    comparison = oecluster.MCSComparison(_pair(BENZENE, TOLUENE), similarity=True)
    facts = comparison.Facts()
    assert facts.is_distance == native.Capability_No
    assert facts.zero_self == native.Capability_No
    assert comparison.Compare(0, 1) == pytest.approx(6.0 / 7.0, abs=1e-6)


def _gate_mols():
    return [_mol(BENZENE, "b"), _mol(TOLUENE, "t"),
            _mol(CYCLOHEXANE, "c"), _mol(MORPHINE, "m")]


@pytest.mark.parametrize("call", [
    lambda dm: oecluster.butina(dm, 0.9),
    lambda dm: oecluster.dbscan(dm, 0.9),
    lambda dm: oecluster.hdbscan(dm, min_cluster_size=2),
    lambda dm: oecluster.agglomerative(dm, n_clusters=2),
    lambda dm: oecluster.k_medoids(dm, n_clusters=2),
    lambda dm: oecluster.activity_landscape(dm, [0.0, 1.0, 2.0, 3.0]),
    lambda dm: oecluster.modelability(dm, ["A", "A", "B", "B"]),
    lambda dm: oecluster.cluster_report(oecluster.butina(dm, 0.9), dm),
])
def test_every_entry_point_accepts_an_mcs_matrix(call):
    # triangle = "unknown" is permissive at the gate and data_integrity is
    # complete, so no entry point needs allow_nonmetric=True.
    dm = oecluster.pdist(_gate_mols(), "mcs")
    assert call(dm) is not None


def test_similarity_string_is_refused_by_the_wrapper():
    # bool("false") is True, so coercing here would silently return the
    # complement of the orientation the caller named. Assigning raw lets the
    # SWIG bool typemap refuse it, which is what every sibling comparison does.
    with pytest.raises(TypeError, match="similarity"):
        oecluster.MCSComparison(
            _pair(BENZENE, TOLUENE),
            similarity="false")  # pyright: ignore[reportArgumentType]


def test_similarity_string_is_refused_by_pdist():
    # The wrapper and _build_mcs are separate assignment sites, so the refusal
    # has to be pinned on both paths.
    with pytest.raises(TypeError, match="similarity"):
        oecluster.pdist(
            _pair(BENZENE, TOLUENE), "mcs",
            similarity="false")  # pyright: ignore[reportArgumentType]


def test_bridging_hydrogens_survive_suppression():
    # OESuppressHydrogens makes a hydrogen implicit on the heavy atom it hangs
    # off, so a hydrogen bonded to two atoms has nowhere to go. Diborane keeps
    # both bridging hydrogens: it has zero heavy-atom bonds, yet it clears the
    # zero-bond guard and scores on B-H bonds alone. A nonzero score is the
    # observable proof that the denominator is the snapshot's bond count.
    dm = oecluster.pdist(_pair(DIBORANE, SINGLE_BRIDGE), "mcs", similarity=True)
    assert dm.condensed[0] == pytest.approx(0.5)


def test_exact_match_level_enforces_ring_membership():
    # Ring membership is the exact preset's largest practical effect and was
    # missing from its documentation until it was measured: a ring never
    # matches a chain under "exact", however identical the bonds themselves.
    default = oecluster.pdist(_pair(CYCLOHEXANE, HEXANE), "mcs",
                              similarity=True, match_level="default")
    exact = oecluster.pdist(_pair(CYCLOHEXANE, HEXANE), "mcs",
                            similarity=True, match_level="exact")
    assert default.condensed[0] == pytest.approx(0.833333, abs=1e-6)
    assert exact.condensed[0] == pytest.approx(0.0)


def test_default_match_level_already_constrains_bond_order():
    # "exact" was documented as adding bond order, but DefaultBonds already
    # carries it and only "loose" ignores it. Pinning both ends keeps the
    # corrected wording honest.
    loose = oecluster.pdist(_pair(ETHANE, ETHENE), "mcs", similarity=True,
                            match_level="loose")
    default = oecluster.pdist(_pair(ETHANE, ETHENE), "mcs", similarity=True,
                              match_level="default")
    assert loose.condensed[0] == pytest.approx(1.0)
    assert default.condensed[0] == pytest.approx(0.0)


def test_multiconformer_molecule_stays_one_item():
    # The docs promise a multi-conformer OEMol is scored once rather than once
    # per pose. Every other fixture here has a single conformer, so a
    # regression routing mcs through the conformer-expanding normalizer would
    # leave them all green while changing the matrix's cardinality.
    first = _multiconformer(BENZENE, "first")
    second = _multiconformer(TOLUENE, "second")
    assert first.NumConfs() == 3
    assert second.NumConfs() == 3
    dm = oecluster.pdist([first, second], "mcs")
    assert dm.num_samples == 2
    assert len(dm.condensed) == 1
    # The single-conformer score, so the extra poses changed nothing about
    # the topology that was compared.
    assert dm.condensed[0] == pytest.approx(0.142857, abs=1e-6)
    assert oecluster.cdist([first], [second], "mcs").shape == (1, 1)


def test_pdist_accepts_graph_molecules():
    # Every example in README.md and docs/python-api.md builds an OEGraphMol,
    # and this comparison rejected that list outright until the OEMolBase
    # overload was added. Pinned against the OEMol score as well as against the
    # literal: an overload that resolved differently for the two input types
    # would still produce a number, so only comparing them catches it.
    graph = oecluster.pdist(_graph_pair(BENZENE, TOLUENE), "mcs")
    oemol = oecluster.pdist(_pair(BENZENE, TOLUENE), "mcs")
    assert graph.num_samples == 2
    assert graph.condensed[0] == pytest.approx(0.142857, abs=1e-6)
    assert graph.condensed[0] == pytest.approx(oemol.condensed[0], abs=1e-12)


def test_direct_constructor_accepts_graph_molecules():
    # The top-level wrapper hands its item list straight to the SWIG class, so
    # the new overload has to resolve on this route too and not only through
    # the registry that pdist goes by.
    comparison = oecluster.MCSComparison(_graph_pair(BENZENE, TOLUENE))
    assert comparison.Size() == 2
    assert comparison.Compare(0, 1) == pytest.approx(0.142857, abs=1e-6)
    native_comparison = native.MCSComparison(_graph_pair(BENZENE, TOLUENE),
                                             native.MCSOptions())
    assert native_comparison.Compare(0, 1) == pytest.approx(0.142857, abs=1e-6)

    # Both calls above pass options, and a defaulted argument expands into its
    # own SWIG dispatch case with its own overload ranking -- the ranking that
    # was measured to come out one way for one-argument calls and the opposite
    # way for two-argument ones. So the no-options form is pinned separately
    # rather than assumed to follow.
    defaulted = native.MCSComparison(_graph_pair(BENZENE, TOLUENE))
    assert defaulted.Compare(0, 1) == pytest.approx(0.142857, abs=1e-6)


def test_charged_hydrogen_survives_suppression():
    # Suppression folds a hydrogen into a heavy atom's implicit hydrogen
    # count, and a formal charge on the hydrogen has nowhere to go in such a
    # count, so the charged atom survives where the neutral one does not.
    # Both halves are pinned because the asymmetry is the whole claim: the
    # neutral hydride suppresses to nothing and is refused, while the charged
    # one keeps its single bond and is accepted.
    with pytest.raises(RuntimeError, match="no bonds after hydrogen"):
        oecluster.pdist(_pair(NEUTRAL_HYDRIDE, NEUTRAL_HYDRIDE), "mcs")
    dm = oecluster.pdist(_pair(CHARGED_HYDRIDE, CHARGED_HYDRIDE), "mcs")
    assert dm.num_samples == 2
    assert dm.condensed[0] == pytest.approx(0.0, abs=1e-6)


def test_molecular_hydrogen_is_the_degenerate_fold():
    # With no heavy atom in the component the fold has to stop somewhere, so
    # one hydrogen absorbs the other and survives as its owner -- neutral and
    # unbridged, which is why the rule is phrased around where the fold can
    # go rather than around the surviving atom's own properties. Nothing is
    # left to score, so the zero-bond guard refuses it.
    with pytest.raises(RuntimeError, match="no bonds after hydrogen"):
        oecluster.pdist(_pair(DIHYDROGEN, DIHYDROGEN), "mcs")

    # The same stray hydrogen alongside a real fragment is bondless, so it
    # cannot reach the bond Tanimoto: benzene scores against toluene exactly
    # as it does without it.
    stray = oecluster.pdist([_mol(f"{DIHYDROGEN}.{BENZENE}", "first"),
                             _mol(TOLUENE, "second")], "mcs")
    assert stray.condensed[0] == pytest.approx(0.142857, abs=1e-6)


def test_direct_constructor_honours_search_mode():
    # The wrapper builds its own kwargs dict rather than routing through
    # _build_mcs, so each option needs coverage here too: deleting
    # search_mode from that dict leaves every pdist test in this file green.
    default = oecluster.MCSComparison(_pair(SUCROSE, MACROLIDE))
    assert default.Compare(0, 1) == pytest.approx(0.605263, abs=1e-6)
    exhaustive = oecluster.MCSComparison(_pair(SUCROSE, MACROLIDE),
                                         search_mode="exhaustive")
    assert exhaustive.Compare(0, 1) == pytest.approx(0.527778, abs=1e-6)


def test_direct_constructor_honours_max_matches():
    # As above: deleting max_matches from the wrapper's kwargs dict is silent
    # against every other test in this file.
    default = oecluster.MCSComparison(_pair(MORPHINE, PENICILLIN_G))
    assert default.Compare(0, 1) == pytest.approx(0.717949, abs=1e-6)
    bounded = oecluster.MCSComparison(_pair(MORPHINE, PENICILLIN_G),
                                      max_matches=1)
    assert bounded.Compare(0, 1) == pytest.approx(0.958333, abs=1e-6)
