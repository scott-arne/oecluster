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
# within ``items_a + items_b``. The last column says which step raises.
_CDIST_COMPARISONS = pytest.mark.parametrize(
    "comparison, kwargs, expected_name, while_normalizing", [
        ("fingerprint", {}, "OEMolBase", False),
        ("rocs", {}, "OEMolBase", False),
        ("superpose", {}, "OEMolBase", False),
        ("sitehopper", {}, "OEMolBase", False),
        ("descriptor", {}, "OEMolBase", True),
        ("descriptor", {"missing": "propagate"}, "OEMolBase", False),
        ("rmsd", {}, "OEMol", True),
        ("rmsd", {"expand_conformers": False}, "OEMol", False),
        ("mcs", {}, "OEMolBase", False),
    ], ids=["fingerprint", "rocs", "superpose", "sitehopper", "descriptor",
            "descriptor-propagate", "rmsd", "rmsd-unexpanded", "mcs"])

# Set sizes, where the bad element sits, and its index in ``items_a +
# items_b``. "b-middle" is set B's item 2 where set A also has an item 2, so
# an index with no set named is ambiguous, and its concatenated index, 5, is a
# valid element of set B. "b-first" and "a-last" sit either side of the
# boundary between the sets.
_PLACEMENTS = {
    "b-middle": (3, 6, "items_b", 2, 5),
    "b-first": (2, 2, "items_b", 0, 2),
    "a-last": (2, 2, "items_a", 1, 1),
}


@_CDIST_COMPARISONS
@pytest.mark.parametrize("placement", sorted(_PLACEMENTS))
def test_cdist_names_the_set_and_the_index_within_it(
        comparison, kwargs, expected_name, while_normalizing, placement):
    size_a, size_b, set_name, index, concatenated = _PLACEMENTS[placement]
    mol = _parsed(oechem.OEMol, "c1ccccc1")
    items_a = [mol] * size_a
    items_b = [mol] * size_b
    (items_a if set_name == "items_a" else items_b)[index] = 1
    with pytest.raises(TypeError) as caught:
        oecluster.cdist(items_a, items_b, comparison, **kwargs)
    assert str(caught.value) == (
        f"{set_name} item {index} is a int, not an {expected_name}")
    # The report it replaced stays attached as the cause, with the index the
    # raising step saw.
    raw_index = index if while_normalizing else concatenated
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


class _ForgedRefusal:
    """An element whose type check raises, on every read, a ``TypeError``
    worded exactly like a conversion's refusal of some other element."""

    def __init__(self):
        self.raised = []

    @property
    def __class__(self):
        error = TypeError("molecule list item 1 is a sentinel, not an OEMolBase")
        self.raised.append(error)
        raise error


@pytest.mark.parametrize("comparison", ["rocs", "descriptor", "rmsd"])
@pytest.mark.parametrize("route", ["pdist", "cdist"])
def test_an_element_s_own_type_error_is_never_renamed(comparison, route):
    """A check that raised has not answered, so its exception surfaces as it
    was, whatever it says. ROCS raises from the combined build, descriptor
    from its native filter and RMSD from its Python item check."""
    mol = _parsed(oechem.OEMol, "c1ccccc1")
    bad = _ForgedRefusal()
    with pytest.raises(TypeError) as caught:
        if route == "pdist":
            oecluster.pdist([mol, bad], comparison)
        else:
            oecluster.cdist([mol], [mol, bad], comparison)
    assert caught.value is bad.raised[-1]
    assert caught.value.__cause__ is None


class _NamelessMeta(type):
    @property
    def __name__(cls):
        raise LookupError("no name")


class _MisnamedMeta(type):
    @property
    def __name__(cls):
        return 5


@pytest.mark.parametrize("meta", [_NamelessMeta, _MisnamedMeta],
                         ids=["raising", "non-string"])
@pytest.mark.parametrize("kwargs", [{}, {"expand_conformers": False}],
                         ids=["expanded", "unexpanded"])
def test_rmsd_names_an_unnamed_type_as_the_conversion_does(meta, kwargs):
    """A type's ``__name__`` can raise or be no string. The native message
    then says ``<unnamed type>``, and RMSD's Python item check, which runs
    only when expanding, must say the same."""
    mol = _parsed(oechem.OEMol, "c1ccccc1")
    odd = meta("Odd", (), {})()
    expected = "item 1 is a <unnamed type>, not an OEMol"
    _assert_reported([lambda: oecluster.pdist([mol, odd], "rmsd", **kwargs)],
                     f"molecule list {expected}")
    _assert_reported(
        [lambda: oecluster.cdist([mol], [mol, odd], "rmsd", **kwargs)],
        f"items_b {expected}")


def _marked(message, index, cls=TypeError):
    error = cls(message)
    setattr(error, _comparisons._ITEM_INDEX_ATTRIBUTE, index)
    return error


class _TypeErrorSubclass(TypeError):
    pass


@pytest.mark.parametrize("make_error", [
    lambda: TypeError("molecule list item 1 is a int, not an OEMolBase"),
    lambda: _marked("molecule list item 1 is a int, not an OEMolBase", 1,
                    _TypeErrorSubclass),
    lambda: _marked("molecule list item 9 is a int, not an OEMolBase", 9),
    lambda: _marked("molecule list item 1 is a int, not an OEMolBase", "1"),
    lambda: _marked("see molecule list item 1 is a int, not an OEMolBase", 1),
], ids=["unmarked", "subclass", "past-both-sets", "non-int-mark",
        "unparsed-message"])
def test_only_a_conversion_s_own_error_is_renamed(make_error):
    """Anything but a plain ``TypeError`` carrying the conversions' mark and
    their message comes out as the same exception object, however it is
    worded."""
    error = make_error()
    with (pytest.raises(TypeError) as caught,
          _comparisons.items_named_by_set(("items_a", 2), ("items_b", 3))):
        raise error
    assert caught.value is error


def test_the_renamed_message_chains_the_conversion_message():
    original = _marked(
        "design unit list item 3 is a int, not an OEDesignUnit", 3)
    with (pytest.raises(TypeError) as caught,
          _comparisons.items_named_by_set(("items_a", 2), ("items_b", 3))):
        raise original
    assert str(caught.value) == "items_b item 1 is a int, not an OEDesignUnit"
    assert caught.value.__cause__ is original
