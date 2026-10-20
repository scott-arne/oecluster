"""ROCS and MCS read each molecule by its own type, wherever it sits in a list.

Both classes take a molecule list through one conversion that copies an
``OEMol`` with every conformer and anything else as a single conformer. Before
it, two competing constructors split the work and the first element picked
between them, so the same molecules gave the right answer, an active-conformer
answer, or a ``TypeError``, depending only on their order.
"""

import numpy as np
import oecluster
import pytest
from oecluster import oecluster as native

oechem = pytest.importorskip("openeye.oechem")


def _conformers(smiles, n):
    oeomega = pytest.importorskip("openeye.oeomega")
    mol = oechem.OEMol()
    oechem.OESmilesToMol(mol, smiles)
    omega = oeomega.OEOmega()
    omega.SetMaxConfs(n)
    omega.SetStrictStereo(False)
    assert omega(mol)
    return mol


@pytest.fixture(scope="module")
def mixed_pair():
    """A single-conformer query and a 30-conformer OEMol whose ACTIVE
    conformer is its worst overlay against that query. Only such a fixture
    can tell the all-conformer reading from the active-conformer one."""
    query = _conformers("CCCCCC(=O)NCCc1ccccc1", 1)
    query_graph = oechem.OEGraphMol(query)
    flex = _conformers("CCCCCC(=O)NCCc1ccccc1", 30)
    worst = max(
        (oecluster.ROCSComparison([query_graph, oechem.OEGraphMol(c)]).Compare(0, 1),
         c.GetIdx())
        for c in flex.GetConfs())
    for conf in flex.GetConfs():
        if conf.GetIdx() == worst[1]:
            flex.SetActive(conf)
    return query, query_graph, flex


def _pdist_cell(mols):
    return np.asarray(oecluster.pdist(mols, comparison="rocs",
                                      num_threads=1).squareform())[0, 1]


def test_rocs_reads_every_conformer_whatever_the_list_order(mixed_pair):
    query, query_graph, flex = mixed_pair
    reference = _pdist_cell([query, flex])
    active_only = _pdist_cell([query_graph, oechem.OEGraphMol(flex)])
    assert active_only - reference > 0.1, "fixture cannot discriminate"
    # The matrix is symmetric, so cell [0, 1] is the query-flex pair in both
    # orderings, including the one that puts flex at index 0.
    for mols in ([query_graph, flex], [flex, query_graph]):
        assert _pdist_cell(mols) == pytest.approx(reference, abs=1e-6)


def test_cdist_reads_set_b_by_its_own_type_whatever_set_a_holds(mixed_pair):
    """cdist builds one comparison over ``items_a + items_b``, so set A's
    first element used to decide how set B was read."""
    query, query_graph, flex = mixed_pair

    def cell(items_a, items_b):
        return oecluster.cdist(items_a, items_b, "rocs",
                               num_threads=1).matrix[0][0]

    reference = cell([query], [flex])
    active_only = cell([query_graph], [oechem.OEGraphMol(flex)])
    assert active_only - reference > 0.1, "fixture cannot discriminate"
    assert cell([query_graph], [flex]) == pytest.approx(reference, abs=1e-6)
    # The reverse sets used to be refused outright.
    assert cell([flex], [query_graph]) == pytest.approx(
        cell([flex], [query]), abs=1e-6)


def test_both_constructor_forms_read_every_conformer(mixed_pair):
    """SWIG emits a separate dispatch case for the defaulted options argument,
    so the one- and two-argument forms are each pinned. ``_comparisons.py``
    builds the two-argument one."""
    query, query_graph, flex = mixed_pair
    options = native.ROCSOptions()
    reference = native.ROCSComparison([query, flex], options).Compare(0, 1)
    active_only = native.ROCSComparison(
        [query_graph, oechem.OEGraphMol(flex)], options).Compare(0, 1)
    assert active_only - reference > 0.1, "fixture cannot discriminate"
    # Compare(i, j) overlays j onto i, so the reverse list is read as (1, 0).
    for mols, i, j in (([query_graph, flex], 0, 1), ([flex, query_graph], 1, 0)):
        assert native.ROCSComparison(mols, options).Compare(i, j) == pytest.approx(
            reference, abs=1e-6)
        assert native.ROCSComparison(mols).Compare(i, j) == pytest.approx(
            reference, abs=1e-6)


def _parsed(kind, smiles):
    mol = kind()
    oechem.OESmilesToMol(mol, smiles)
    return mol


def test_mcs_accepts_a_mixed_list_in_either_order():
    """MCS matches bonds, so its scores carry no conformer consequence; what
    used to depend on the order was whether the list was accepted at all.
    Benzene's six bonds all match toluene's seven, so the distance is 1 - 6/7.
    """
    benzene = _parsed(oechem.OEGraphMol, "c1ccccc1")
    toluene = _parsed(oechem.OEMol, "Cc1ccccc1")
    expected = 1.0 / 7.0
    for mols in ([benzene, toluene], [toluene, benzene]):
        assert oecluster.pdist(mols, "mcs").condensed[0] == pytest.approx(
            expected, abs=1e-6)
        assert oecluster.MCSComparison(mols).Compare(0, 1) == pytest.approx(
            expected, abs=1e-6)
        assert native.MCSComparison(mols).Compare(0, 1) == pytest.approx(
            expected, abs=1e-6)
        assert native.MCSComparison(mols, native.MCSOptions()).Compare(
            0, 1) == pytest.approx(expected, abs=1e-6)
    for items_a, items_b in (([benzene], [toluene]), ([toluene], [benzene])):
        assert oecluster.cdist(items_a, items_b, "mcs").matrix[0][0] == (
            pytest.approx(expected, abs=1e-6))


_CLASSES = {"rocs": "ROCSComparison", "mcs": "MCSComparison"}


def _call(comparison, route, mols, good):
    if route == "class":
        getattr(oecluster, _CLASSES[comparison])(mols)
    elif route == "pdist":
        oecluster.pdist(mols, comparison)
    elif route == "cdist":
        # The bad element sits in set A, where its index in the concatenated
        # list the comparison is built over is its index in set A.
        oecluster.cdist(mols, [good], comparison)
    else:
        native_class = getattr(native, _CLASSES[comparison])
        options = getattr(native, _CLASSES[comparison].replace(
            "Comparison", "Options"))()
        with pytest.raises(TypeError) as one_argument:
            native_class(mols)
        assert str(one_argument.value) == _expected_message(mols)
        native_class(mols, options)


def _expected_message(mols):
    index, bad = next((i, m) for i, m in enumerate(mols)
                      if not isinstance(m, oechem.OEMolBase))
    return f"molecule list item {index} is a {type(bad).__name__}, not an OEMolBase"


@pytest.mark.parametrize("comparison", ["rocs", "mcs"])
@pytest.mark.parametrize("route", ["class", "pdist", "cdist", "native"])
@pytest.mark.parametrize("position", ["first", "later"])
@pytest.mark.parametrize("kind", [oechem.OEGraphMol, oechem.OEMol],
                         ids=["graph", "oemol"])
def test_a_non_molecule_is_reported_by_position(comparison, route, position, kind):
    """Index, offending type and expected type, on every route, and never the
    dispatcher's ``new_<Class>`` message.

    A different non-molecule type at each position, so the type name has to
    come from the element at the reported index.
    """
    good = _parsed(kind, "c1ccccc1")
    other = _parsed(kind, "Cc1ccccc1")
    if position == "first":
        mols = ["x", good, other]
    else:
        mols = [good, other, 7]
    with pytest.raises(TypeError) as caught:
        _call(comparison, route, mols, good)
    message = str(caught.value)
    assert message == _expected_message(mols)
    assert "new_" not in message


@pytest.mark.parametrize("comparison", ["rocs", "mcs"])
@pytest.mark.parametrize("kind", [oechem.OEGraphMol, oechem.OEMol],
                         ids=["graph", "oemol"])
def test_a_molecule_that_wraps_nothing_says_so(comparison, kind):
    """An instance made by ``__new__`` passes the OEMolBase check but has no
    molecule behind it, so it is not called a non-OEMolBase."""
    empty = kind.__new__(kind)
    mols = [_parsed(kind, "c1ccccc1"), empty]
    native_class = getattr(native, _CLASSES[comparison])
    with pytest.raises(TypeError) as caught:
        native_class(mols)
    assert str(caught.value) == (
        f"molecule list item 1 is a {kind.__name__} that wraps no molecule")


def test_rmsd_still_refuses_a_graph_molecule():
    """The ROCS and MCS conversion is scoped to their constructors. RMSD
    measures conformers and keeps its strict ``OEMol`` conversion, on both
    constructor forms."""
    mols = [_parsed(oechem.OEMol, "c1ccccc1"), _parsed(oechem.OEGraphMol, "c1ccccc1")]
    with pytest.raises(TypeError, match=r"not an OEMol\b"):
        native.RMSDComparison(mols)
    with pytest.raises(TypeError, match=r"not an OEMol\b"):
        native.RMSDComparison(mols, native.RMSDOptions())
    with pytest.raises(TypeError, match=r"not an OEMol\b"):
        oecluster.RMSDComparison(mols)
