"""The SWIG layer's input conversions: what they keep alive, what they report."""

import subprocess
import sys
import textwrap

import oecluster
import pytest

oechem = pytest.importorskip("openeye.oechem")

_CLEAR_MID_CALL = textwrap.dedent("""
    import threading, time, oecluster
    from openeye import oechem
    smiles = ["c1ccc2ccccc2c1CCN", "C1CCC(CC1)c1ccncc1",
              "O=C(Nc1ccccc1)c1ccco1", "c1ccc(cc1)C1CCNCC1"]
    mols = []
    for i in range(20000):
        m = oechem.OEGraphMol(); oechem.OESmilesToMol(m, smiles[i % 4])
        mols.append(m)
    stamps = {{}}
    def clear_later():
        time.sleep(0.2)
        mols.clear()
        stamps["cleared"] = time.perf_counter()
    worker = threading.Thread(target=clear_later)
    stamps["start"] = time.perf_counter()
    worker.start()
    out = {call}
    stamps["end"] = time.perf_counter()
    worker.join()
    print(stamps["start"], stamps["cleared"], stamps["end"], len(out))
""")


@pytest.mark.parametrize("call", [
    "oecluster.murcko_scaffolds(mols, num_threads=1)",
    "oecluster.murcko(mols, num_threads=1)",
])
def test_clearing_the_list_mid_call_does_not_crash(call):
    """Another thread emptying the list must not free molecules the call reads.

    Run in a child process: before the fix this crashed the interpreter, which
    would take the whole test session down with it. The clear is stamped after
    it finishes and the end by the statement after the call, so the ordering
    shows the list was empty before the main thread got past the call; without
    it a fast machine could finish first and pass without testing anything.
    """
    result = subprocess.run(
        [sys.executable, "-X", "faulthandler", "-c",
         _CLEAR_MID_CALL.format(call=call)],
        check=False, capture_output=True, text=True, timeout=600)
    assert result.returncode == 0, (result.returncode, result.stderr[-2000:])
    start, cleared, end, count = result.stdout.split()
    assert float(start) < float(cleared) < float(end), (
        "the list was not emptied during the call; enlarge the input")
    assert int(count) > 0


def test_a_refused_list_releases_its_molecules():
    """The error path must drop the references the conversion took.

    A leaked reference is otherwise invisible: nothing crashes, the molecules
    simply never die. The bad element sits last, so the refusal comes only
    after every good element has been converted. The native class is called
    directly so no Python-layer check can refuse the list before the typemap
    sees it; the anchored message proves the typemap is what refused it.
    """
    from oecluster import oecluster as native

    good = oechem.OEGraphMol()
    oechem.OESmilesToMol(good, "c1ccccc1")
    mols = [good] * 199 + [1]
    before = sys.getrefcount(good)
    for _ in range(2000):
        with pytest.raises(
                TypeError,
                match=r"^molecule list item 199 is a int, not an OEMolBase$"):
            native.FingerprintComparison(mols)
    assert sys.getrefcount(good) == before


class _RaisingClass:
    """An object whose ``isinstance`` check raises instead of answering.

    ``isinstance`` falls back to reading ``__class__``, which this refuses.
    """

    @property
    def __class__(self):
        raise LookupError("no class")


@pytest.mark.parametrize("comparison", [
    oecluster.ROCSComparison, oecluster.MCSComparison], ids=["rocs", "mcs"])
def test_a_raising_element_in_a_strict_oemol_list_propagates(comparison):
    """A check that raised has not answered, so its exception must surface.

    The OEMol at index 0 selects the strict ``OEMol`` list conversion; the
    element after it makes the isinstance check raise, which used to be
    replaced by a generic TypeError about the element's type.
    """
    mol = oechem.OEMol()
    oechem.OESmilesToMol(mol, "c1ccccc1")
    with pytest.raises(LookupError, match="^no class$"):
        comparison([mol, _RaisingClass()])


def test_a_raising_element_in_a_design_unit_list_propagates():
    """The same rule on the design-unit list conversion."""
    with pytest.raises(LookupError, match="^no class$"):
        oecluster.SuperposeComparison(
            [oechem.OEDesignUnit(), _RaisingClass()])


def test_a_raising_fingerprint_batch_check_propagates():
    """The same rule on the fingerprint-batch conversion.

    ``take_fingerprints`` has a single wrapper, so the object reaches the in
    typemap itself rather than an overload dispatcher's typecheck, which has
    no way to report an error.
    """
    pytest.importorskip("oefp")
    from oecluster import oecluster as native

    with pytest.raises(LookupError, match="^no class$"):
        native.take_fingerprints(_RaisingClass(), [0])


def _benzene(kind):
    mol = kind()
    oechem.OESmilesToMol(mol, "c1ccccc1")
    return mol


@pytest.mark.parametrize(
    "route", ["fingerprint", "rocs", "mcs", "superpose"])
def test_a_raising_first_element_propagates(route):
    """The overload dispatcher must not turn a raising check into "no".

    It reads element 0 alone, through a typecheck that cannot report an
    error. Answering "no" there left every overload ruled out and replaced
    the exception with the dispatcher's generic TypeError, before any
    conversion ran.
    """
    comparison, rest = {
        "fingerprint": (
            oecluster.FingerprintComparison, [_benzene(oechem.OEGraphMol)]),
        "rocs": (oecluster.ROCSComparison, [_benzene(oechem.OEMol)]),
        "mcs": (oecluster.MCSComparison, [_benzene(oechem.OEMol)]),
        "superpose": (oecluster.SuperposeComparison, [oechem.OEDesignUnit()]),
    }[route]
    with pytest.raises(LookupError, match="^no class$"):
        comparison([_RaisingClass(), *rest])


def test_a_raising_batch_check_on_an_overloaded_function_propagates():
    """``isim`` has two overloads, so its typecheck runs before the typemap."""
    pytest.importorskip("oefp")
    from oecluster import oecluster as native

    with pytest.raises(LookupError, match="^no class$"):
        native.isim(_RaisingClass())


def _raising_pointer(base, attribute):
    """An instance of a ``base`` subclass whose ``attribute`` raises.

    Built with ``__new__`` so ``__init__`` never tries to set the attribute.
    The isinstance check passes, so the failure comes from reading the
    wrapped pointer.
    """
    def refuse(self):
        raise LookupError("no pointer")

    subclass = type("Raising" + base.__name__, (base,),
                    {attribute: property(refuse)})
    return subclass.__new__(subclass)


@pytest.mark.parametrize("route", ["molbase", "oemol", "design_unit", "batch"])
def test_a_raising_pointer_read_propagates(route):
    """Reading the wrapped pointer can raise, too; that is not "no pointer".

    Only a missing attribute means the object wraps nothing. Any other
    exception used to be cleared and reported as a type or null-reference
    error naming the wrong cause.
    """
    if route == "molbase":
        def call():
            oecluster.murcko_scaffolds(
                [_benzene(oechem.OEGraphMol),
                 _raising_pointer(oechem.OEGraphMol, "this")])
    elif route == "oemol":
        def call():
            oecluster.ROCSComparison(
                [_benzene(oechem.OEMol),
                 _raising_pointer(oechem.OEMol, "this")])
    elif route == "design_unit":
        def call():
            oecluster.SuperposeComparison(
                [oechem.OEDesignUnit(),
                 _raising_pointer(oechem.OEDesignUnit, "this")])
    else:
        oefp = pytest.importorskip("oefp")
        from oecluster import oecluster as native

        def call():
            native.take_fingerprints(
                _raising_pointer(oefp.OEFPBatch, "_native"), [0])
    with pytest.raises(LookupError, match="^no pointer$"):
        call()


_DEGRADED_OEFP = textwrap.dedent("""
    import oefp
    from oecluster import oecluster as native
    batch = oefp.OEFPBatch.from_fingerprints(
        [oefp.OEFP.from_on_bits(8, [1]), oefp.OEFP.from_on_bits(8, [2])])
    del oefp.OEFPBatch
    try:
        native.take_fingerprints(batch, [0])
    except Exception as exc:
        print(type(exc).__name__, exc)
    else:
        print("accepted")
""")


def test_a_degraded_oefp_is_reported_as_itself():
    """A broken ``oefp`` must not read as "Expected an oefp.OEFPBatch".

    The batch-type lookup failing is the cause, and the caller is holding a
    real batch. Run in a child process: deleting ``oefp.OEFPBatch`` here
    would break every later test that builds a batch.
    """
    pytest.importorskip("oefp")
    result = subprocess.run(
        [sys.executable, "-c", _DEGRADED_OEFP],
        check=False, capture_output=True, text=True, timeout=600)
    assert result.returncode == 0, result.stderr[-2000:]
    assert result.stdout.startswith("AttributeError "), result.stdout
    assert "OEFPBatch" in result.stdout
