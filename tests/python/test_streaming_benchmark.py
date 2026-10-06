"""Smoke tests for the streaming threshold benchmark harness."""

from __future__ import annotations

import importlib.util
import json
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "benchmarks" / "streaming_threshold.py"


def _load_benchmark_module():
    spec = importlib.util.spec_from_file_location(
        "oecluster_streaming_benchmark", SCRIPT)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_the_library_is_deterministic_and_every_smiles_parses():
    from openeye import oechem

    module = _load_benchmark_module()
    library = module.library_smiles(200, seed=3)
    assert library == module.library_smiles(200, seed=3)
    assert library != module.library_smiles(200, seed=4)
    assert len(set(library)) == 200
    mol = oechem.OEGraphMol()
    for smi in library:
        mol.Clear()
        assert oechem.OESmilesToMol(mol, smi), smi


def test_the_library_refuses_more_items_than_combinations():
    module = _load_benchmark_module()
    combinations = len(module.SCAFFOLDS) * len(module.SUBSTITUENTS) ** 3
    with pytest.raises(ValueError, match="combinations"):
        module.library_smiles(combinations + 1, seed=0)


@pytest.mark.skipif(sys.platform == "win32",
                    reason="peak RSS is read through the POSIX resource module")
def test_a_tiny_run_reports_both_paths_and_skips_the_matrix_where_asked():
    completed = subprocess.run(
        [sys.executable, str(SCRIPT), "--sizes", "60",
         "--streaming-only-sizes", "80", "--thresholds", "0.4519", "--json"],
        capture_output=True, text=True, check=True)
    rows = json.loads(completed.stdout)
    # Two sizes, two algorithms, one threshold, two paths.
    assert len(rows) == 8
    by_key = {(row["algorithm"], row["n_items"], row["path"]): row
              for row in rows}
    for algorithm in ("butina", "dbscan"):
        streaming = by_key[(algorithm, 60, "streaming")]
        matrix = by_key[(algorithm, 60, "matrix")]
        assert streaming["attempted"]
        assert matrix["attempted"]
        # 0.4519 sits far from every Tanimoto ratio, so pdist and Compare
        # cannot disagree about any pair and the paths cluster identically.
        assert streaming["n_clusters"] == matrix["n_clusters"]
        assert streaming["rss_peak_bytes"] >= streaming["rss_inputs_bytes"] > 0
        skipped = by_key[(algorithm, 80, "matrix")]
        assert skipped["attempted"] is False
        assert skipped["matrix_bytes"] == 8 * (80 * 79 // 2)
        assert by_key[(algorithm, 80, "streaming")]["attempted"]
