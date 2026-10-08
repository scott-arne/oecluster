"""Smoke tests for the streaming spanning-tree benchmark harness."""

from __future__ import annotations

import importlib.util
import json
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "benchmarks" / "streaming_mst.py"


def _load_benchmark_module():
    spec = importlib.util.spec_from_file_location(
        "oecluster_streaming_mst_benchmark", SCRIPT)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_the_module_imports_and_shares_the_threshold_library():
    module = _load_benchmark_module()
    assert module.library_smiles(50, seed=3) == module.library_smiles(50, seed=3)
    load = module.load_average()
    assert load is None or len(load) == 3


@pytest.mark.skipif(sys.platform == "win32",
                    reason="peak RSS is read through the POSIX resource module")
def test_a_tiny_run_reports_both_paths_and_skips_the_matrix_where_asked():
    completed = subprocess.run(
        [sys.executable, str(SCRIPT), "--sizes", "60",
         "--streaming-only-sizes", "80", "--threads", "0", "1",
         "--n-clusters", "4", "--json"],
        capture_output=True, text=True, check=True)
    rows = json.loads(completed.stdout)
    # Two sizes, two algorithms, two thread counts, two paths.
    assert len(rows) == 16
    by_key = {(row["algorithm"], row["n_items"], row["num_threads"],
               row["path"]): row for row in rows}
    for algorithm in ("hdbscan", "single"):
        for threads in (0, 1):
            streaming = by_key[(algorithm, 60, threads, "streaming")]
            matrix = by_key[(algorithm, 60, threads, "matrix")]
            assert streaming["attempted"]
            assert matrix["attempted"]
            assert streaming["rss_peak_bytes"] >= streaming["rss_inputs_bytes"] > 0
            assert streaming["pdist_seconds"] < matrix["pdist_seconds"] + 1.0
            assert len(streaming["load_average"]) == 3
            skipped = by_key[(algorithm, 80, threads, "matrix")]
            assert skipped["attempted"] is False
            assert skipped["matrix_bytes"] == 8 * (80 * 79 // 2)
            assert by_key[(algorithm, 80, threads, "streaming")]["attempted"]
    # A fixed n_clusters cut returns that many clusters on both paths.
    assert by_key[("single", 60, 0, "streaming")]["n_clusters_found"] == 4
    assert by_key[("single", 60, 0, "matrix")]["n_clusters_found"] == 4


@pytest.mark.skipif(sys.platform == "win32",
                    reason="peak RSS is read through the POSIX resource module")
def test_the_paths_and_algorithms_can_be_restricted():
    completed = subprocess.run(
        [sys.executable, str(SCRIPT), "--sizes", "40",
         "--streaming-only-sizes", "--algorithms", "single", "--n-clusters", "4",
         "--paths", "matrix", "--json"],
        capture_output=True, text=True, check=True)
    rows = json.loads(completed.stdout)
    assert [(row["algorithm"], row["path"]) for row in rows] == [
        ("single", "matrix")]
