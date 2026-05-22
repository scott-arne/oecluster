"""Smoke tests for the BitBirch benchmark harness."""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]


def test_bitbirch_benchmark_supports_refinement_workflows():
    env = os.environ.copy()
    env["PYTHONPATH"] = (
        str(ROOT / "python")
        if not env.get("PYTHONPATH")
        else f"{ROOT / 'python'}{os.pathsep}{env['PYTHONPATH']}"
    )

    result = subprocess.run(
        [
            sys.executable,
            "benchmarks/bitbirch.py",
            "--sizes",
            "8",
            "--bits",
            "16",
            "--density",
            "0.2",
            "--threshold",
            "0.4",
            "--branching-factor",
            "2",
            "--reassign-top-clusters",
            "2",
            "--repeats",
            "1",
            "--warmups",
            "0",
            "--workflows",
            "cluster",
            "recluster",
            "reassign",
            "prune",
            "prune_reassign",
        ],
        cwd=ROOT,
        env=env,
        check=False,
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, result.stderr
    assert "| workflow |" in result.stdout
    assert "| cluster |" in result.stdout
    assert "| recluster |" in result.stdout
    assert "| reassign |" in result.stdout
    assert "| prune |" in result.stdout
    assert "| prune_reassign |" in result.stdout
