"""A cross-runtime object's C++ type comes from its SWIG descriptor.

oecluster's cross-runtime type check, which ``_swig_cast_verdict`` exposes,
finds the descriptor a ``this`` carries in the SWIG runtime that created it,
and converts through that runtime's cast graph. The verdict is 1 when the
object converts to the target, 0 when it does not, 2 when it wraps nothing,
3 when its descriptor is in no runtime this module can read.

The molecule and design-unit list conversions decide each element's type
with this check, after the ``isinstance`` call each still makes. Before
they used it, some inputs reached undefined behaviour in them. The tests of
those inputs run in a child interpreter under faulthandler, so a crash fails
the test instead of ending the session.
"""

import subprocess
import sys
import textwrap

import oecluster
import pytest
from oecluster import oecluster as native

oechem = pytest.importorskip("openeye.oechem")

_MOLBASE = "_p_OEChem__OEMolBase"
_MCMOLBASE = "_p_OEChem__OEMCMolBase"
_DESIGN_UNIT = "_p_OEBio__OEDesignUnit"
_BATCH = "_p_OEFP__OEFPBatch"


@pytest.mark.parametrize("make, target, verdict", [
    (oechem.OEGraphMol, _MCMOLBASE, 0),
    (oechem.OEMol, _MCMOLBASE, 1),
    (oechem.OEGraphMol, _MOLBASE, 1),
    (oechem.OEMol, _MOLBASE, 1),
    (oechem.OEQMol, _MOLBASE, 1),
    (oechem.OEDesignUnit, _MOLBASE, 0),
    (oechem.OEDesignUnit, _DESIGN_UNIT, 1),
], ids=["graph-mcmol", "mol-mcmol", "graph-molbase", "mol-molbase",
        "qmol-molbase", "du-molbase", "du-du"])
def test_the_verdict_follows_openeye_s_cast_graph(make, target, verdict):
    """Pins the OpenEye cast-graph facts these verdicts depend on. If OpenEye
    changes its type hierarchy, these fail instead of a verdict changing
    silently."""
    assert native._swig_cast_verdict(make(), target) == verdict


def test_an_oefp_batch_converts_through_this_module_s_runtime():
    oefp = pytest.importorskip("oefp")
    batch = oefp.OEFPBatch.from_fingerprints(
        [oefp.OEFP.from_on_bits(8, [1]), oefp.OEFP.from_on_bits(8, [1, 2])])
    assert native._swig_cast_verdict(batch._native, _BATCH) == 1


def test_an_object_made_by_new_wraps_nothing():
    assert native._swig_cast_verdict(
        oechem.OEGraphMol.__new__(oechem.OEGraphMol), _MOLBASE) == 2


def test_an_object_without_this_wraps_nothing():
    assert native._swig_cast_verdict(object(), _MOLBASE) == 2


def test_an_unknown_target_name_is_refused():
    with pytest.raises(ValueError, match="no descriptor"):
        native._swig_cast_verdict(oechem.OEGraphMol(), "_p_NoSuchType")


class _RaisingThis(oechem.OEGraphMol):
    @property
    def this(self):
        raise LookupError("no this")

    @this.setter
    def this(self, value):
        self.__dict__["_own_this"] = value


def test_an_exception_reading_this_propagates():
    with pytest.raises(LookupError, match="no this"):
        native._swig_cast_verdict(_RaisingThis(), _MOLBASE)


class _Holder:
    """Carries whatever ``this`` it is given. The helper reads only ``this``,
    so a holder needs no OpenEye base class to reach every step."""

    def __init__(self, this):
        self.this = this


def test_a_this_that_is_not_a_swig_object_is_not_the_type():
    """Pins step 2: a ``this`` whose type is not really named
    ``SwigPyObject`` is not the target type, and nothing in it is read."""
    assert native._swig_cast_verdict(_Holder(5), _MOLBASE) == 0


class _MisreportingMeta(type):
    """Makes every class it builds report itself, through ordinary attribute
    lookup, as a ``SwigPyObject`` far larger than any real one."""

    @property
    def __name__(cls):
        return "SwigPyObject"

    @property
    def __basicsize__(cls):
        return 1_000_000_000


@pytest.mark.parametrize("real_name, slots", [
    ("SwigPyObject", ()),
    ("_Small", ("ptr", "ty")),
], ids=["real-size-too-small", "real-name-wrong"])
def test_a_type_that_misreports_its_name_and_size_is_not_the_type(
        real_name, slots):
    """Pins step 2's reading of the real name and size through ``type``'s
    own descriptors, which a metaclass cannot override.

    Each case passes one real check and fails the other. The first is really
    named ``SwigPyObject`` but holds only an object header, too small for
    ``ptr`` and ``ty``. The second is large enough, with both slots set, but
    is really named ``_Small``. A check that believed the reported values
    would read past the first object's end, and would read the second's
    ``ty`` and answer 3.
    """
    cls = _MisreportingMeta(real_name, (), {"__slots__": slots})
    assert (cls.__name__, cls.__basicsize__) == ("SwigPyObject", 1_000_000_000)
    fake = cls()
    for slot in slots:
        setattr(fake, slot, object())
    assert native._swig_cast_verdict(_Holder(fake), _MOLBASE) == 0


def _unreadable_swig_object():
    """An object that passes step 2 genuinely, since its class is really
    named ``SwigPyObject`` and is large enough for ``ptr`` and ``ty``. Both
    slots hold Python objects, so ``ty`` is non-null and is no descriptor."""
    swig_object = type("SwigPyObject", (), {"__slots__": ("ptr", "ty")})
    fake = swig_object()
    fake.ptr = object()
    fake.ty = object()
    return fake


def test_a_descriptor_in_no_readable_runtime_is_unreadable():
    """Pins steps 3 and 4: a ``ty`` that neither runtime 5 nor runtime 4
    holds answers 3, whatever ``ptr`` holds, with both runtimes loaded."""
    assert {"swig_runtime_data4", "swig_runtime_data5"} <= set(sys.modules)
    assert native._swig_cast_verdict(
        _Holder(_unreadable_swig_object()), _MOLBASE) == 3


_RETRY_RUNTIME_4 = textwrap.dedent("""
    import sys

    import oecluster
    from oecluster import oecluster as native

    loaded = sorted(name for name in sys.modules
                    if name.startswith(("openeye", "swig_runtime_data4")))
    assert not loaded, loaded

    swig_object = type("SwigPyObject", (), {"__slots__": ("ptr", "ty")})
    fake = swig_object()
    fake.ptr = object()
    fake.ty = object()
    holder = type("Holder", (), {})()
    holder.this = fake
    print(native._swig_cast_verdict(holder, "_p_OEChem__OEMolBase"))

    from openeye import oechem

    mol = oechem.OEGraphMol()
    oechem.OESmilesToMol(mol, "c1ccccc1")
    print(native._swig_cast_verdict(mol, "_p_OEChem__OEMolBase"))
    oecluster.FingerprintComparison([mol, mol])
    print("converted")
""")


def test_a_failed_runtime_4_import_is_retried():
    """Pins step 3's rule that a failed lookup of runtime 4 is never cached.

    In a fresh interpreter that has not imported OpenEye, a ``ty`` outside
    runtime 5 makes the helper look for runtime 4, which is absent, so the
    verdict is 3. Once ``openeye.oechem`` is imported, the next lookup must
    find runtime 4: a molecule converts, and a comparison accepts molecules.
    A helper that kept the first failure would answer 3 for the molecule.
    """
    result = subprocess.run(
        [sys.executable, "-X", "faulthandler", "-c", _RETRY_RUNTIME_4],
        check=False, capture_output=True, text=True, timeout=300)
    assert result.returncode == 0, (result.returncode, result.stderr[-2000:])
    assert result.stdout.split() == ["3", "1", "converted"]


class _FreshThis:
    """Each read of ``this`` builds a benzene and returns its SwigPyObject,
    which owns the molecule and which nothing else refers to."""

    @property
    def this(self):
        fresh = oechem.OEGraphMol()
        oechem.OESmilesToMol(fresh, "c1ccccc1")
        return fresh.this


def test_a_this_produced_on_each_read_is_held():
    """Pins that the hook's path through the helper nets zero references.

    A ``this`` made on each read, which only the helper's reference keeps
    alive, gives 1 on every one of many reads. A molecule's ``this`` keeps
    its reference count across many calls, so nothing leaks and nothing is
    released twice, and the molecule is still usable afterwards. The
    hand-over in ``*held`` itself is not visible here: the hook discards the
    converted pointer and releases ``*held`` at once, so a helper that
    dropped ``this`` early would pass too.
    ``test_a_this_produced_on_each_read_is_computed_on`` is where the
    hand-over is visible, through the list conversions.
    """
    fresh = _FreshThis()
    for _ in range(200):
        assert native._swig_cast_verdict(fresh, _MOLBASE) == 1

    mol = oechem.OEGraphMol()
    oechem.OESmilesToMol(mol, "c1ccccc1")
    this = mol.this
    before = sys.getrefcount(this)
    for _ in range(1000):
        assert native._swig_cast_verdict(mol, _MOLBASE) == 1
    assert sys.getrefcount(this) == before
    assert oechem.OEMolToSmiles(mol) == "c1ccccc1"


# ---- The molecule and design-unit list conversions ----

# Run ahead of each child's own body. Each class is an object whose own
# account of itself a conversion must not trust: ``Spoof``'s ``__class__``
# names a type it is not, a ``Borrowing`` item's ``this`` is another
# object's, and ``Producing``'s ``this`` is a new object on each read, which
# nothing but the reader keeps alive. OpenEye's own ``__init__`` reads
# ``this`` while it constructs, and a ``this`` that answered another object
# then would be linked to the new one, so ``Borrowing`` and ``Producing``
# raise AttributeError for ``this``, as an instance with none does, until
# that ``__init__`` has stored its own.
_CHILD_PRELUDE = textwrap.dedent(r'''
    import sys

    import oecluster
    from oecluster import oecluster as native
    from openeye import oechem


    def report(label, call):
        """Prints the label, a tab, then what the call returned or raised."""
        try:
            value = call()
        except Exception as error:
            outcome = f"raised {type(error).__name__}: {error}"
        else:
            outcome = f"returned {value!r}"
        print(label, outcome, sep="\t", flush=True)


    def parsed(kind, smiles):
        mol = kind()
        oechem.OESmilesToMol(mol, smiles)
        return mol


    def read(path):
        stream = oechem.oemolistream(path)
        mols = [oechem.OEGraphMol(mol) for mol in stream.GetOEGraphMols()]
        stream.close()
        return mols


    class Spoof(oechem.OEGraphMol):
        """``isinstance`` calls it an OEMol, since its ``__class__`` answers
        OEMol, while its ``this`` wraps the OEGraphMol it really is."""

        __class__ = property(lambda self: oechem.OEMol)


    class Borrowing:
        """Once constructed, ``this`` answers the object given instead."""

        def __init__(self, borrowed):
            super().__init__()
            self.__dict__["_borrowed"] = borrowed

        @property
        def this(self):
            try:
                return self.__dict__["_borrowed"]
            except KeyError:
                raise AttributeError("this") from None

        @this.setter
        def this(self, value):
            self.__dict__["_own_this"] = value


    class BorrowingGraphMol(Borrowing, oechem.OEGraphMol):
        pass


    class BorrowingMol(Borrowing, oechem.OEMol):
        pass


    class BorrowingUnit(Borrowing, oechem.OEDesignUnit):
        pass


    class Producing(oechem.OEGraphMol):
        """Once constructed, each read of ``this`` copies source into a new
        OEGraphMol and answers the copy's ``this``, which owns the copy and
        which nothing else refers to."""

        def __init__(self, source):
            super().__init__()
            self.__dict__["_source"] = source

        @property
        def this(self):
            try:
                source = self.__dict__["_source"]
            except KeyError:
                raise AttributeError("this") from None
            return oechem.OEGraphMol(source).this

        @this.setter
        def this(self, value):
            self.__dict__["_own_this"] = value
''')


def _run_child(body, *args):
    """Runs the prelude and ``body`` in a fresh interpreter under
    faulthandler, with ``args`` as ``sys.argv[1:]``, and returns what each
    label reported, as ``{label: "returned ..." or "raised ..."}``."""
    script = _CHILD_PRELUDE + textwrap.dedent(body)
    result = subprocess.run(
        [sys.executable, "-X", "faulthandler", "-c", script,
         *(str(arg) for arg in args)],
        check=False, capture_output=True, text=True, timeout=300)
    assert result.returncode == 0, (
        result.returncode, result.stdout[-2000:], result.stderr[-2000:])
    return dict(line.split("\t", 1) for line in result.stdout.splitlines())


@pytest.fixture(scope="module")
def shaped_pair(tmp_path_factory):
    """The path of an OEB file holding two single-conformer molecules with
    3D coordinates, which ROCS requires, so every child reads the same two.
    The two have different Murcko frameworks."""
    oeomega = pytest.importorskip("openeye.oeomega")
    path = tmp_path_factory.mktemp("shaped") / "pair.oeb"
    stream = oechem.oemolostream(str(path))
    for smiles in ("c1ccc2ccccc2c1CCN", "Oc1ccc(cc1)C1CCCC1"):
        mol = oechem.OEMol()
        oechem.OESmilesToMol(mol, smiles)
        omega = oeomega.OEOmega()
        omega.SetMaxConfs(1)
        omega.SetStrictStereo(False)
        assert omega(mol)
        oechem.OEWriteMolecule(stream, mol)
    stream.close()
    return path


# Every way a list reaches a comparison: the class, pdist, cdist with the
# list's first element as set A, and the raw constructor in both forms, since
# SWIG emits a separate dispatch case for the defaulted options argument.
_ROUTES = ("class", "pdist", "cdist", "native", "native with options")

_SPOOF_SCORES = r'''
    comparison, name, path = sys.argv[1:4]
    native_class = getattr(native, comparison)
    options = getattr(native, comparison.replace("Comparison", "Options"))
    first, second = read(path)
    for label, items in (("plain", [first, second]),
                         ("spoof", [Spoof(first), second])):
        report(label + " class",
               lambda: getattr(oecluster, comparison)(items).Compare(0, 1))
        report(label + " pdist",
               lambda: oecluster.pdist(items, name).condensed[0])
        report(label + " cdist",
               lambda: oecluster.cdist(items[:1], items[1:], name).matrix[0][0])
        report(label + " native", lambda: native_class(items).Compare(0, 1))
        report(label + " native with options",
               lambda: native_class(items, options()).Compare(0, 1))
'''


@pytest.mark.parametrize("comparison, name", [
    ("ROCSComparison", "rocs"), ("MCSComparison", "mcs")], ids=["rocs", "mcs"])
def test_rocs_and_mcs_read_a_spoof_as_the_molecule_it_wraps(
        comparison, name, shaped_pair):
    """A spoof passes both of the conversion's ``isinstance`` checks, as an
    OEMolBase and as an OEMCMolBase, but wraps an OEGraphMol, which does not
    convert to OEMCMolBase. It is copied as the one-conformer molecule it is,
    so every route scores it exactly as the plain OEGraphMol."""
    outcomes = _run_child(_SPOOF_SCORES, comparison, name, shaped_pair)
    for route in _ROUTES:
        assert outcomes[f"plain {route}"].startswith("returned "), outcomes
        assert outcomes[f"spoof {route}"] == outcomes[f"plain {route}"], route


_RMSD_SPOOF = r'''
    first = read(sys.argv[1])[0]
    items = [oechem.OEMol(first), Spoof(first)]
    report("class", lambda: oecluster.RMSDComparison(items))
    report("pdist", lambda: oecluster.pdist(items, "rmsd"))
    report("cdist", lambda: oecluster.cdist(items, [items[0]], "rmsd"))
    report("native", lambda: native.RMSDComparison(items))
    report("native with options",
           lambda: native.RMSDComparison(items, native.RMSDOptions()))
'''


def test_rmsd_refuses_a_spoof_by_position(shaped_pair):
    """RMSD's ``isinstance(OEMol)`` accepts the spoof, which wraps an
    OEGraphMol, so the conversion refuses it as not an OEMol. ``cdist``
    names the set the spoof is in. The pose is the same molecule as the
    OEMol beside it, so a conversion that accepted the spoof would build a
    comparison."""
    outcomes = _run_child(_RMSD_SPOOF, shaped_pair)
    refused = "raised TypeError: {} item 1 is a Spoof, not an OEMol"
    assert outcomes == {
        route: refused.format("items_a" if route == "cdist" else "molecule list")
        for route in _ROUTES}


_REFUSAL = r'''
    comparison, holder, good, donor = sys.argv[1:5]
    native_class = getattr(native, comparison)
    options = getattr(native, comparison.replace("Comparison", "Options"))


    def made(kind):
        if kind == "OEDesignUnit":
            return oechem.OEDesignUnit()
        return parsed(getattr(oechem, kind), "c1ccccc1")


    borrowed = 5 if donor == "int" else made(donor).this
    items = [made(good), globals()[holder](borrowed)]
    report("native", lambda: native_class(items))
    report("native with options", lambda: native_class(items, options()))
'''

# Per comparison: the holder class, whose OpenEye base passes the
# conversion's isinstance check; the kind of a good item; the list's name;
# and the type the conversion requires, as its message names it.
# Fingerprint uses the OEMolBase list conversion, ROCS and MCS their own,
# RMSD the strict OEMol one and Superpose, given design units, the
# design-unit one.
_HOLDERS = {
    "FingerprintComparison": (
        "BorrowingGraphMol", "OEGraphMol", "molecule list", "OEMolBase"),
    "ROCSComparison": (
        "BorrowingGraphMol", "OEGraphMol", "molecule list", "OEMolBase"),
    "MCSComparison": (
        "BorrowingGraphMol", "OEGraphMol", "molecule list", "OEMolBase"),
    "RMSDComparison": ("BorrowingMol", "OEMol", "molecule list", "OEMol"),
    "SuperposeComparison": (
        "BorrowingUnit", "OEDesignUnit", "design unit list", "OEDesignUnit"),
}


def _assert_refused_by_position(comparison, donor):
    holder, good, list_name, required = _HOLDERS[comparison]
    outcomes = _run_child(_REFUSAL, comparison, holder, good, donor)
    refused = (f"raised TypeError: {list_name} item 1 is a {holder}, "
               f"not an {required}")
    assert outcomes == {"native": refused, "native with options": refused}


@pytest.mark.parametrize("comparison, donor", [
    ("FingerprintComparison", "OEDesignUnit"),
    ("ROCSComparison", "OEDesignUnit"),
    ("MCSComparison", "OEDesignUnit"),
    ("RMSDComparison", "OEDesignUnit"),
    ("SuperposeComparison", "OEGraphMol"),
], ids=["fingerprint", "rocs", "mcs", "rmsd", "superpose"])
def test_an_item_holding_another_object_s_this_is_refused(comparison, donor):
    """An item whose class passes the conversion's ``isinstance`` check, but
    whose ``this`` is another object's genuine SwigPyObject, is refused as
    the wrong type: a molecule holding a design unit, and a design unit
    holding a molecule."""
    _assert_refused_by_position(comparison, donor)


@pytest.mark.parametrize("comparison", sorted(_HOLDERS))
def test_an_item_whose_this_is_a_plain_object_is_refused(comparison):
    """An int is no SwigPyObject, so nothing in it is read as one."""
    _assert_refused_by_position(comparison, "int")


_PRODUCED = r'''
    call, path = sys.argv[1:3]
    calls = {
        "FingerprintComparison":
            lambda items: oecluster.FingerprintComparison(items).Compare(0, 1),
        "murcko_scaffolds":
            lambda items: oecluster.murcko_scaffolds(items, num_threads=2),
        "ROCSComparison":
            lambda items: oecluster.ROCSComparison(items).Compare(0, 1),
        "MCSComparison":
            lambda items: oecluster.MCSComparison(items).Compare(0, 1),
    }
    first, second = read(path)
    report("plain", lambda: calls[call]([first, second] * 40))
    report("produced",
           lambda: calls[call]([Producing(first), Producing(second)] * 40))
'''


@pytest.mark.parametrize("call", [
    "FingerprintComparison", "murcko_scaffolds", "ROCSComparison",
    "MCSComparison"])
def test_a_this_produced_on_each_read_is_computed_on(call, shaped_pair):
    """Each read of a ``Producing`` item's ``this`` makes a molecule that
    only that ``this`` keeps alive. The conversion holds the ``this`` it
    verified until it has no further use for the pointer: until the call
    returns for the OEMolBase list, which Fingerprint and Murcko take, and
    until the element is copied for ROCS and MCS. So each computes on the
    molecules the reads produced, as on the plain molecules they copy.

    The list repeats two items 40 times each. The OEMolBase list reads each
    element's ``this`` once, so 80 molecules are produced; ROCS and MCS read
    it twice, once per cast, so 160 are.

    What this pins for Fingerprint and Murcko is that the verified ``this``
    objects outlive the call, which reads the molecules with the GIL
    released. Eighty items make two of Murcko's 64-item chunks, and with
    ``num_threads=2`` it may start two workers, which it does when
    OpenEye's memory pool mode allows threads, as the default does. Workers
    claim chunks as they become free, so which worker reads which chunk is
    not fixed, and one may read both. The test therefore allows two workers
    to read at once but does not claim that they do."""
    outcomes = _run_child(_PRODUCED, call, shaped_pair)
    assert outcomes["plain"].startswith("returned "), outcomes
    assert outcomes["produced"] == outcomes["plain"]


def test_converting_one_list_repeatedly_gives_one_answer(shaped_pair):
    """Pins unchanged behaviour: the same molecules, through Fingerprint's
    OEMolBase list conversion and through ROCS's own, score the same every
    time. It pins stability, not caching: no verdict is cached."""
    stream = oechem.oemolistream(str(shaped_pair))
    mols = [oechem.OEGraphMol(mol) for mol in stream.GetOEGraphMols()]
    stream.close()
    for comparison in (oecluster.FingerprintComparison,
                       oecluster.ROCSComparison):
        scores = {comparison(mols).Compare(0, 1) for _ in range(50)}
        assert len(scores) == 1, (comparison.__name__, scores)


_UNREADABLE = r'''
    swig_object = type("SwigPyObject", (), {"__slots__": ("ptr", "ty")})
    fake = swig_object()
    fake.ptr = object()
    fake.ty = object()
    items = [parsed(oechem.OEGraphMol, "c1ccccc1"), BorrowingGraphMol(fake)]
    report("class", lambda: oecluster.FingerprintComparison(items))
    report("native", lambda: native.FingerprintComparison(items))
'''


def test_an_item_whose_descriptor_no_runtime_holds_is_refused():
    """An item whose ``this`` passes as a SwigPyObject but carries a ``ty``
    that neither runtime holds raises the list RuntimeError, which names
    the runtimes loaded and no index."""
    outcomes = _run_child(_UNREADABLE)
    refused = (
        "raised RuntimeError: molecule list item is a BorrowingGraphMol whose "
        "SWIG type is in no SWIG runtime oecluster can read (loaded: "
        "swig_runtime_data4, swig_runtime_data5)")
    assert outcomes == {"class": refused, "native": refused}


def _good_item(kind):
    if kind == "OEDesignUnit":
        return oechem.OEDesignUnit()
    mol = getattr(oechem, kind)()
    oechem.OESmilesToMol(mol, "c1ccccc1")
    return mol


@pytest.mark.parametrize("comparison", sorted(_HOLDERS))
def test_a_refused_list_releases_every_this_it_verified(comparison):
    """Guards the new reference handling, so it is not red before it: the
    conversions used to keep no reference to ``this`` at all.

    Each good item's ``this`` is verified and held, for the call by the
    OEMolBase list and until the item is copied by the others. The last
    item is refused, after every good one has been verified, so each of
    those holds must be released on the way out. A leak would add 200
    references to the good item's ``this`` per call."""
    _, good_kind, list_name, required = _HOLDERS[comparison]
    good = _good_item(good_kind)
    this = good.this
    before = sys.getrefcount(this)
    refused = rf"^{list_name} item 200 is a int, not an {required}$"
    for _ in range(50):
        with pytest.raises(TypeError, match=refused):
            getattr(native, comparison)([good] * 200 + [5])
    assert sys.getrefcount(this) == before


def test_a_converted_list_releases_every_this_it_verified():
    """The same guard on the path where the OEMolBase list converts: the
    tuple of verified ``this`` objects is released when the call returns."""
    good = _good_item("OEGraphMol")
    this = good.this
    before = sys.getrefcount(this)
    for _ in range(50):
        native.FingerprintComparison([good] * 200)
    assert sys.getrefcount(this) == before
