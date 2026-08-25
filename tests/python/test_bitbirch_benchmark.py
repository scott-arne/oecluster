"""Smoke tests for the BitBirch benchmark harness."""

from __future__ import annotations

import importlib.util
import os
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]

# Keep in step with benchmarks/bitbirch.py. Unless --native-only is passed the
# harness compares against the upstream BitBirch reference, which ships as a
# source checkout rather than a distribution, so those runs are skipped wherever
# the sibling repository is absent.
BITBIRCH_REPO = Path("/Users/johnss51/Development/python/bitbirch")

requires_bitbirch_reference = pytest.mark.skipif(
    not BITBIRCH_REPO.is_dir(),
    reason=f"BitBirch reference checkout not found at {BITBIRCH_REPO}",
)


def _load_benchmark_module():
    spec = importlib.util.spec_from_file_location(
        "oecluster_bitbirch_benchmark",
        ROOT / "benchmarks" / "bitbirch.py",
    )
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


@requires_bitbirch_reference
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
            "--regime",
            "duplicate_blocks",
            "--prototype-count",
            "4",
            "--threshold",
            "0.4",
            "--branching-factor",
            "2",
            "--reassign-top-clusters",
            "2",
            "--num-threads",
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


@requires_bitbirch_reference
def test_bitbirch_benchmark_handles_duplicate_block_reassign():
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
            "64",
            "--bits",
            "128",
            "--density",
            "0.2",
            "--regime",
            "duplicate_blocks",
            "--prototype-count",
            "4",
            "--threshold",
            "0.65",
            "--branching-factor",
            "2",
            "--reassign-top-clusters",
            "4",
            "--num-threads",
            "2",
            "--repeats",
            "1",
            "--warmups",
            "0",
            "--workflows",
            "reassign",
        ],
        cwd=ROOT,
        env=env,
        check=False,
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, result.stderr
    assert "| reassign |" in result.stdout


def test_bitbirch_benchmark_native_only_skips_reference_import(monkeypatch, capsys):
    benchmark = _load_benchmark_module()

    def fail_reference_import():
        raise AssertionError("native-only benchmarks must not import Python BitBirch")

    monkeypatch.setattr(benchmark, "load_reference_bitbirch", fail_reference_import)
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "bitbirch.py",
            "--native-only",
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
            "--repeats",
            "1",
            "--warmups",
            "0",
            "--workflows",
            "cluster",
        ],
    )

    benchmark.main()

    output = capsys.readouterr().out
    assert "| cluster |" in output
    assert "| n/a |" in output
