"""A bad element in a comparison's list is reported by its position.

Fingerprint, Descriptor and RMSD each have one constructor whose defaulted
options argument makes SWIG emit one- and two-argument forms, and Superpose
has a constructor for each of two list types, so all four sit behind SWIG's
overload dispatcher. A dispatcher that rejected a bad element raised its own
"new_<Class>" message, which names no element, and whether it did depended on
where the bad element sat. ROCS and MCS are covered the same way in
``test_rocs_conformer_order.py``.

``cdist`` builds one comparison over ``items_a + items_b``, so it names the
set and the index within that set, for every comparison.
"""

import oecluster
import pytest
from oecluster import _comparisons
from oecluster import oecluster as native

oechem = pytest.importorskip("openeye.oechem")


def _parsed(kind, smiles):
    mol = kind()
    oechem.OESmilesToMol(mol, smiles)
    return mol


# Class name, options name, the type a list element must have, and that
# type's name as the message gives it.
_CLASSES = {
    "fingerprint": ("FingerprintComparison", "FingerprintOptions",
                    oechem.OEMolBase, "OEMolBase"),
    "descriptor": ("DescriptorComparison", "DescriptorOptions",
                   oechem.OEMolBase, "OEMolBase"),
    "rmsd": ("RMSDComparison", "RMSDOptions", oechem.OEMol, "OEMol"),
    "superpose": ("SuperposeComparison", "SuperposeOptions",
                  oechem.OEMolBase, "OEMolBase"),
}

_ROUTES = ["class", "pdist", "cdist", "native"]


def _calls(comparison, route, items, good):
    """Each call one route makes with ``items``, and the name it reports the
    list by, or ``None`` where that is the conversion's own.

    The raw constructor is called in both forms, since SWIG emits a separate
    dispatch case for the defaulted options argument. On ``cdist`` the items
    are set A, so the reported index is the one in ``items``.
    """
    class_name, options_name, _, _ = _CLASSES[comparison]
    if route == "class":
        return [lambda: getattr(oecluster, class_name)(items)], None
    if route == "pdist":
        return [lambda: oecluster.pdist(items, comparison)], None
    if route == "cdist":
        return [lambda: oecluster.cdist(items, [good], comparison)], "items_a"
    native_class = getattr(native, class_name)
    options = getattr(native, options_name)()
    return [lambda: native_class(items),
            lambda: native_class(items, options)], None


def _expected(items, accepted, expected_name, list_name):
    index, bad = next((i, item) for i, item in enumerate(items)
                      if not isinstance(item, accepted))
    prefix = list_name or "molecule list"
    return (f"{prefix} item {index} is a {type(bad).__name__}, "
            f"not an {expected_name}")


def _assert_reported(calls, expected):
    for call in calls:
        with pytest.raises(TypeError) as caught:
            call()
        message = str(caught.value)
        assert message == expected
        assert "new_" not in message


@pytest.mark.parametrize("comparison", sorted(_CLASSES))
@pytest.mark.parametrize("route", _ROUTES)
@pytest.mark.parametrize("position", ["first", "later"])
def test_a_non_molecule_is_reported_by_position(comparison, route, position):
    """Index, offending type and expected type, on every route, and never the
    dispatcher's ``new_<Class>`` message.

    A different non-molecule type at each position, so the type name has to
    come from the element at the reported index.
    """
    _, _, accepted, expected_name = _CLASSES[comparison]
    kind = oechem.OEMol if comparison == "rmsd" else oechem.OEGraphMol
    good = _parsed(kind, "c1ccccc1")
    other = _parsed(kind, "Cc1ccccc1")
    items = ["x", good, other] if position == "first" else [good, other, 7]
    calls, list_name = _calls(comparison, route, items, good)
    _assert_reported(calls, _expected(items, accepted, expected_name, list_name))


def test_superpose_still_takes_both_list_types():
    """The molecule typecheck accepts any list, so only the order in which
    the dispatcher examines the two typechecks keeps a design-unit list on
    the design-unit constructor. Had the molecule one come first, this list
    would be refused as "not an OEMolBase"."""
    mols = [_parsed(oechem.OEGraphMol, "c1ccccc1"),
            _parsed(oechem.OEGraphMol, "Cc1ccccc1")]
    units = [oechem.OEDesignUnit(), oechem.OEDesignUnit()]
    for items in (mols, units):
        assert oecluster.SuperposeComparison(items).Size() == 2
        assert native.SuperposeComparison(items).Size() == 2
        assert native.SuperposeComparison(
            items, native.SuperposeOptions()).Size() == 2


def test_superpose_keeps_its_empty_list_message():
    """An empty list is accepted by both typechecks, so it reaches the
    design-unit constructor, which is examined first, and is refused as
    before."""
    for call in (lambda: native.SuperposeComparison([]),
                 lambda: native.SuperposeComparison([], native.SuperposeOptions())):
        with pytest.raises(RuntimeError,
                           match=r"^SuperposeComparison: empty structure list$"):
            call()


@pytest.mark.parametrize("route", _ROUTES)
def test_superpose_reports_a_mixed_list_at_index_one(route):
    """A list's first element picks the constructor, and the conversion then
    names the first element that does not belong."""
    mol = _parsed(oechem.OEGraphMol, "c1ccccc1")
    unit = oechem.OEDesignUnit()

    calls, list_name = _calls("superpose", route, [mol, unit], mol)
    _assert_reported(
        calls, f"{list_name or 'molecule list'} item 1 is a OEDesignUnit, "
               "not an OEMolBase")

    calls, list_name = _calls("superpose", route, [unit, "x"], unit)
    _assert_reported(
        calls, f"{list_name or 'design unit list'} item 1 is a str, "
               "not an OEDesignUnit")


class _PoseLike:
    """Answers ``NumConfs`` like a two-conformer molecule but is not one."""

    def NumConfs(self):
        return 2


@pytest.mark.parametrize("index", [0, 1])
def test_rmsd_checks_each_item_before_expanding_conformers(index):
    """Conformer expansion calls ``NumConfs`` and ``GetTitle`` on each item,
    so an object answering the first used to fail on the second, with an
    ``AttributeError`` naming a method rather than an element."""
    mol = _parsed(oechem.OEMol, "c1ccccc1")
    items = [mol, mol]
    items[index] = _PoseLike()
    expected = f"item {index} is a _PoseLike, not an OEMol"
    _assert_reported([lambda: oecluster.pdist(items, "rmsd")],
                     f"molecule list {expected}")
    _assert_reported([lambda: oecluster.cdist([mol], items, "rmsd")],
                     f"items_b {expected}")


@pytest.mark.parametrize("route", _ROUTES)
def test_rmsd_still_refuses_a_graph_molecule_by_position(route):
    items = [_parsed(oechem.OEMol, "c1ccccc1"),
             _parsed(oechem.OEGraphMol, "c1ccccc1")]
    calls, list_name = _calls("rmsd", route, items, items[0])
    _assert_reported(
        calls, f"{list_name or 'molecule list'} item 1 is a OEGraphMol, "
               "not an OEMol")


@pytest.mark.parametrize("route", _ROUTES)
def test_rmsd_reports_a_molecule_that_wraps_nothing(route):
    """An ``OEMol`` made by ``__new__`` passes the type check but has no
    molecule behind it. Expansion used to call its ``NumConfs``, which raised
    OEChem's own error naming a method."""
    empty = oechem.OEMol.__new__(oechem.OEMol)
    items = [_parsed(oechem.OEMol, "c1ccccc1"), empty]
    calls, list_name = _calls("rmsd", route, items, items[0])
    _assert_reported(
        calls,
        f"{list_name or 'molecule list'} item 1 is a OEMol that wraps no molecule")


# Each comparison, through both steps of cdist that can raise: the per-set
# normalization, where descriptor's complete-case filter and RMSD's item check
# report an index within the set, and the combined build, which reports one
# within ``items_a + items_b``. The last column is that raw index.
@pytest.mark.parametrize("comparison, kwargs, expected_name, raw_index", [
    ("fingerprint", {}, "OEMolBase", 5),
    ("rocs", {}, "OEMolBase", 5),
    ("superpose", {}, "OEMolBase", 5),
    ("sitehopper", {}, "OEMolBase", 5),
    ("descriptor", {}, "OEMolBase", 2),
    ("descriptor", {"missing": "propagate"}, "OEMolBase", 5),
    ("rmsd", {}, "OEMol", 2),
    ("rmsd", {"expand_conformers": False}, "OEMol", 5),
    ("mcs", {}, "OEMolBase", 5),
], ids=["fingerprint", "rocs", "superpose", "sitehopper", "descriptor",
        "descriptor-propagate", "rmsd", "rmsd-unexpanded", "mcs"])
def test_cdist_names_the_set_and_the_index_within_it(
        comparison, kwargs, expected_name, raw_index):
    """The bad element is set B's item 2. Set A also has an item 2, so an
    index with no set named is ambiguous, and its index in ``items_a +
    items_b``, 5, is a valid element of set B."""
    mol = _parsed(oechem.OEMol, "c1ccccc1")
    items_a = [mol, mol, mol]
    items_b = [mol, mol, 1, mol, mol, mol]
    with pytest.raises(TypeError) as caught:
        oecluster.cdist(items_a, items_b, comparison, **kwargs)
    assert str(caught.value) == f"items_b item 2 is a int, not an {expected_name}"
    # The report it replaced stays attached as the cause.
    assert str(caught.value.__cause__) == (
        f"molecule list item {raw_index} is a int, not an {expected_name}")


def test_cdist_names_a_design_unit_set():
    unit = oechem.OEDesignUnit()
    with pytest.raises(TypeError) as caught:
        oecluster.cdist([unit, unit, unit], [unit, unit, 1, unit, unit, unit],
                        "superpose")
    assert str(caught.value) == "items_b item 2 is a int, not an OEDesignUnit"


def test_cdist_passes_any_other_type_error_through_untouched():
    mol = _parsed(oechem.OEMol, "c1ccccc1")
    with pytest.raises(TypeError) as caught:
        oecluster.cdist([mol], [mol], "rocs", colour=1)
    assert str(caught.value) == "Unknown kwargs for rocs comparison: ['colour']"
    assert caught.value.__cause__ is None


class _TypeErrorSubclass(TypeError):
    pass


@pytest.mark.parametrize("error", [
    TypeError("see molecule list item 1 is a int, not an OEMolBase"),
    TypeError("molecule list item one is a int, not an OEMolBase"),
    TypeError("molecule list item 9 is a int, not an OEMolBase"),
    _TypeErrorSubclass("molecule list item 1 is a int, not an OEMolBase"),
], ids=["unanchored", "no-index", "past-both-sets", "subclass"])
def test_only_the_conversion_message_is_renamed(error):
    """Anything but the conversions' message raised as a plain ``TypeError``
    comes out as the same exception object."""
    with (pytest.raises(TypeError) as caught,
          _comparisons.items_named_by_set(("items_a", 2), ("items_b", 3))):
        raise error
    assert caught.value is error


def test_the_renamed_message_chains_the_conversion_message():
    original = TypeError("design unit list item 3 is a int, not an OEDesignUnit")
    with (pytest.raises(TypeError) as caught,
          _comparisons.items_named_by_set(("items_a", 2), ("items_b", 3))):
        raise original
    assert str(caught.value) == "items_b item 1 is a int, not an OEDesignUnit"
    assert caught.value.__cause__ is original
