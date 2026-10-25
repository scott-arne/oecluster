"""SWIG accepts None as a null ``self``. Every wrapped accessor, method and
destructor refuses it with TypeError: getters and methods would otherwise
read through a null pointer, and setters and destructors would return
silently. Run in a child process, since the failure this guards against is a
crash."""

import subprocess
import sys

import pytest

_CALLS = {
    "getter": "n.PDistOptions.num_threads.fget(None)",
    "setter": "n.PDistOptions.num_threads.fset(None, 3)",
    "method": "n.DenseStorage.Finalize(None)",
    "const method": "n.DenseStorage.Get(None, 0, 1)",
    "overloaded method": "n.DoubleVector.resize(None, 3)",
    "destructor": "n._oecluster.delete_DenseStorage(None)",
}


@pytest.mark.parametrize("label", sorted(_CALLS))
def test_a_none_self_is_refused(label):
    script = (
        "from oecluster import oecluster as n\n"
        "try:\n"
        f"    {_CALLS[label]}\n"
        "except Exception as e:\n"
        "    print('raised', type(e).__name__, e)\n"
        "else:\n"
        "    print('returned')\n")
    result = subprocess.run([sys.executable, "-X", "faulthandler", "-c", script],
                            capture_output=True, text=True, timeout=120,
                            check=False)
    report = (result.returncode, result.stdout, result.stderr[-2000:])
    assert result.returncode == 0, report
    assert result.stdout.startswith("raised TypeError"), report
    assert "invalid null reference" in result.stdout, report
