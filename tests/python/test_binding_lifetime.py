"""Lifetime of Python objects the SWIG layer hands to native code."""

import subprocess
import sys
import textwrap

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
        stamps["cleared"] = time.perf_counter()
        mols.clear()
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
    would take the whole test session down with it. The timestamps prove the
    list was emptied before the call returned; without them a fast machine
    could finish first and pass without testing anything.
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
