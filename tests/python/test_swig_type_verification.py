"""A cross-runtime object's C++ type comes from its SWIG descriptor.

oecluster's cross-runtime type check, which ``_swig_cast_verdict`` exposes,
finds the descriptor a ``this`` carries in the SWIG runtime that created it,
and converts through that runtime's cast graph. The verdict is 1 when the
object converts to the target, 0 when it does not, 2 when it wraps nothing,
3 when its descriptor is in no runtime this module can read.
"""

import subprocess
import sys
import textwrap

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
