"""Smoke tests for documented examples."""

import os
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]


def _run_example(script_name):
    env = os.environ.copy()
    pythonpath = str(ROOT / "python")
    if env.get("PYTHONPATH"):
        pythonpath = pythonpath + os.pathsep + env["PYTHONPATH"]
    env["PYTHONPATH"] = pythonpath
    return subprocess.run(
        [sys.executable, str(ROOT / "examples" / script_name)],
        check=True,
        capture_output=True,
        env=env,
        text=True,
    )


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
