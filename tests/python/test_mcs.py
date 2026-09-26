"""MCS comparison: registration, options and gate facts."""

import oecluster
import pytest
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


def _mol(smiles, title="mol"):
    """Parse a SMILES. MCS is topological, so no embedding step is needed."""
    mol = oechem.OEMol()
    oechem.OESmilesToMol(mol, smiles)
    mol.SetTitle(title)
    return mol


def _pair(first, second):
    return [_mol(first, "first"), _mol(second, "second")]


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
