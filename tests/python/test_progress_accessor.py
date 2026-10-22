"""Reading PDistOptions.progress / CDistOptions.progress gives back the callable.

Before 5.26.1 the getter returned an owned std::function pointer SWIG had no
destructor for: an opaque object Python could not use, leaked on every read.
"""

import subprocess
import sys

import pytest
from oecluster import oecluster as native

OPTIONS = [native.PDistOptions, native.CDistOptions]


@pytest.mark.parametrize("options_type", OPTIONS)
def test_reading_progress_returns_the_callable_set(options_type):
    options = options_type()
    def callback(completed, total):
        pass
    options.progress = callback
    assert options.progress is callback


@pytest.mark.parametrize("options_type", OPTIONS)
def test_progress_reads_none_when_unset_or_cleared(options_type):
    options = options_type()
    assert options.progress is None
    options.progress = lambda completed, total: None
    options.progress = None
    assert options.progress is None


def test_a_cdist_progress_callable_is_called():
    """cdist builds its options through the same conversion, so a progress
    callable passed to it must still be called. pdist's is covered in
    test_pdist.py."""
    import oecluster
    oechem = pytest.importorskip("openeye.oechem")
    mols = []
    for smiles in ("CCO", "CCN", "CCC", "c1ccccc1"):
        mol = oechem.OEGraphMol()
        oechem.OESmilesToMol(mol, smiles)
        mols.append(mol)
    calls = []
    oecluster.cdist(mols[:2], mols[2:], "fingerprint",
                    progress=lambda completed, total: calls.append((completed, total)))
    assert calls, "the progress callable was never called"
    assert calls[-1][0] == calls[-1][1]


def test_reading_progress_leaks_nothing():
    script = (
        "from oecluster import oecluster as n\n"
        "def run():\n"
        "    for t in (n.PDistOptions, n.CDistOptions):\n"
        "        o = t()\n"
        "        o.progress = lambda a, b: None\n"
        "        for _ in range(3):\n"
        "            o.progress\n"
        "run()\n")
    result = subprocess.run([sys.executable, "-X", "faulthandler", "-c", script],
                            capture_output=True, text=True, timeout=120,
                            check=False)
    assert result.returncode == 0, result.stderr[-2000:]
    assert "memory leak of type 'std::function" not in result.stdout + result.stderr
