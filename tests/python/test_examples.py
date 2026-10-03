"""Smoke tests for documented examples."""

import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def _run_example(script_name):
    env = os.environ.copy()
    # The example must import the same oecluster this test process did. In a
    # development checkout that is the source tree, which the subprocess only
    # sees through PYTHONPATH; on a wheel-testing CI runner it is the installed
    # package, and putting the source tree first would shadow it with a copy
    # that has no loadable extension.
    import oecluster

    source_package = ROOT / "python" / "oecluster"
    if Path(oecluster.__file__).resolve().parent == source_package.resolve():
        pythonpath = str(ROOT / "python")
        if env.get("PYTHONPATH"):
            pythonpath = pythonpath + os.pathsep + env["PYTHONPATH"]
        env["PYTHONPATH"] = pythonpath
    result = subprocess.run(
        [sys.executable, str(ROOT / "examples" / script_name)],
        check=False,
        capture_output=True,
        env=env,
        text=True,
    )
    # check=True would hide the script's own error; the traceback is the
    # evidence worth having when an example fails on a CI runner.
    assert result.returncode == 0, (
        f"{script_name} exited {result.returncode}\n"
        f"--- stdout ---\n{result.stdout}\n--- stderr ---\n{result.stderr}")
    return result


def test_quickstart_smiles_example_runs():
    """The README quickstart should run from inline SMILES data."""
    result = _run_example("quickstart_smiles.py")

    assert "clusters:" in result.stdout
    assert "representatives:" in result.stdout


def test_rank_representatives_example_runs():
    """The representative-ranking example should print top-k selections."""
    result = _run_example("rank_representatives.py")

    assert "score selection:" in result.stdout
    assert "diversity selection:" in result.stdout
