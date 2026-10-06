"""Benchmark HDBSCAN and single linkage from a comparison against the matrix path.

Both paths start from the same molecules and the same fingerprint comparison:
the streaming path passes the molecules with ``comparison="fingerprint"`` and
holds no matrix, the matrix path materializes ``pdist()`` first. Each
measurement runs in its own subprocess, because peak RSS is a process
high-water mark that whichever path ran first in a shared process would set
for both. At a size where the matrix would not fit, the matrix path is not
attempted; its would-be size is reported instead. Every record carries the
1-minute, 5-minute and 15-minute load averages taken as it started, because
the spanning-tree pass waits for its slowest thread at every step and a busy
machine slows it more than the load suggests.

Example::

    python benchmarks/streaming_mst.py --sizes 20000 \\
        --streaming-only-sizes 100000 --threads 0 1 --json

``--algorithms single --paths matrix`` restricts a run to what 5.19.0 can
also execute, for the before-and-after comparison of matrix single linkage.
"""

from __future__ import annotations

import argparse
import json
import os
import statistics
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

# The library generator and the RSS reader are shared with the threshold
# benchmark, which lives beside this file; a direct script run already has this
# directory on sys.path, an importlib load does not.
sys.path.insert(0, str(Path(__file__).resolve().parent))
from streaming_threshold import (
    library_smiles,
    matrix_bytes,
    peak_rss_bytes,
)

ALGORITHMS = ("hdbscan", "single")
PATHS = ("streaming", "matrix")


def load_average() -> list[float] | None:
    """Return the 1, 5 and 15-minute load averages, or None where unsupported.

    ``os.getloadavg`` does not exist on Windows.
    """
    getloadavg = getattr(os, "getloadavg", None)
    return list(getloadavg()) if getloadavg is not None else None


def measure(spec: dict[str, Any]) -> dict[str, Any]:
    """Run one measurement in this process and return its record.

    :param spec: ``algorithm``, ``path``, ``n_items``, ``seed``,
        ``num_threads``, ``min_cluster_size`` and ``n_clusters``.
    :returns: The spec plus ``seconds`` (everything after the molecules are
        parsed), ``pdist_seconds``, ``cluster_seconds``,
        ``n_clusters_found``, ``rss_inputs_bytes``, ``rss_peak_bytes`` and
        ``load_average``.
    """
    import oecluster
    from openeye import oechem

    mols = []
    for smi in library_smiles(spec["n_items"], spec["seed"]):
        mol = oechem.OEGraphMol()
        if not oechem.OESmilesToMol(mol, smi):
            raise ValueError(f"unparseable library SMILES: {smi}")
        mols.append(mol)
    rss_inputs = peak_rss_bytes()
    load = load_average()

    started = time.perf_counter()
    if spec["path"] == "streaming":
        source: Any = mols
        extra: dict[str, Any] = {"comparison": "fingerprint"}
    else:
        source = oecluster.pdist(mols, "fingerprint",
                                 num_threads=spec["num_threads"])
        extra = {}
    cluster_started = time.perf_counter()
    if spec["algorithm"] == "hdbscan":
        result = oecluster.hdbscan(
            source, min_cluster_size=spec["min_cluster_size"],
            num_threads=spec["num_threads"], **extra)
    else:
        result = oecluster.agglomerative(
            source, linkage="single", n_clusters=spec["n_clusters"],
            num_threads=spec["num_threads"], **extra)
    finished = time.perf_counter()
    return {**spec, "seconds": finished - started,
            "pdist_seconds": cluster_started - started,
            "cluster_seconds": finished - cluster_started,
            "n_clusters_found": result.num_clusters,
            "rss_inputs_bytes": rss_inputs, "rss_peak_bytes": peak_rss_bytes(),
            "load_average": load}


def measure_in_subprocess(spec: dict[str, Any]) -> dict[str, Any]:
    """Run :func:`measure` in a fresh interpreter and return its record."""
    completed = subprocess.run(
        [sys.executable, __file__, "--child", json.dumps(spec)],
        capture_output=True, text=True, check=False)
    if completed.returncode != 0:
        raise RuntimeError(
            f"benchmark child failed with exit status {completed.returncode} "
            f"for spec {spec}:\n{completed.stderr}")
    return json.loads(completed.stdout.strip().splitlines()[-1])


def benchmark(args: argparse.Namespace) -> list[dict[str, Any]]:
    """Run every requested measurement and return one record per row."""
    rows = []
    sizes = [(n, True) for n in args.sizes]
    sizes += [(n, False) for n in args.streaming_only_sizes]
    for n_items, matrix_fits in sizes:
        for algorithm in args.algorithms:
            for num_threads in args.threads:
                spec = {"algorithm": algorithm, "n_items": n_items,
                        "seed": args.seed, "num_threads": num_threads,
                        "min_cluster_size": args.min_cluster_size,
                        "n_clusters": args.n_clusters}
                for path in args.paths:
                    if path == "matrix" and not matrix_fits:
                        rows.append({**spec, "path": path, "attempted": False,
                                     "matrix_bytes": matrix_bytes(n_items)})
                        continue
                    runs = [measure_in_subprocess({**spec, "path": path})
                            for _ in range(args.repeats)]
                    row = dict(runs[0])
                    for key in ("seconds", "pdist_seconds", "cluster_seconds"):
                        row[key] = statistics.median(r[key] for r in runs)
                    row["rss_peak_bytes"] = max(r["rss_peak_bytes"] for r in runs)
                    row["attempted"] = True
                    row["matrix_bytes"] = matrix_bytes(n_items)
                    rows.append(row)
    return rows


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Benchmark HDBSCAN and single linkage from a comparison "
                    "against the matrix path.")
    parser.add_argument("--sizes", type=int, nargs="*", default=[20000],
                        help="sizes at which both paths run")
    parser.add_argument("--streaming-only-sizes", type=int, nargs="*",
                        default=[100000],
                        help="sizes at which the matrix path is not attempted")
    parser.add_argument("--algorithms", nargs="+", choices=ALGORITHMS,
                        default=list(ALGORITHMS))
    parser.add_argument("--paths", nargs="+", choices=PATHS,
                        default=list(PATHS))
    parser.add_argument("--threads", type=int, nargs="+", default=[0],
                        help="num_threads values; 0 is the library default")
    parser.add_argument("--min-cluster-size", type=int, default=5)
    parser.add_argument("--n-clusters", type=int, default=50)
    parser.add_argument("--repeats", type=int, default=1)
    parser.add_argument("--seed", type=int, default=7)
    parser.add_argument("--json", action="store_true", dest="as_json")
    parser.add_argument("--child", help=argparse.SUPPRESS)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    """Run the benchmark and print a table, or JSON with ``--json``."""
    args = parse_args(argv)
    if args.child is not None:
        print(json.dumps(measure(json.loads(args.child))))
        return
    rows = benchmark(args)
    if args.as_json:
        print(json.dumps(rows, indent=2))
        return
    print("| algorithm | n | threads | path | seconds | cluster s | clusters | "
          "peak RSS MB | load 1m | matrix MB |")
    print("| --- | ---: | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: |")
    for row in rows:
        matrix_mb = row["matrix_bytes"] / 1e6
        if not row["attempted"]:
            print(f"| {row['algorithm']} | {row['n_items']} | "
                  f"{row['num_threads']} | {row['path']} | not attempted | | | | | "
                  f"{matrix_mb:.0f} |")
            continue
        load = row["load_average"]
        load_text = f"{load[0]:.1f}" if load is not None else "n/a"
        print(f"| {row['algorithm']} | {row['n_items']} | {row['num_threads']} | "
              f"{row['path']} | {row['seconds']:.2f} | "
              f"{row['cluster_seconds']:.2f} | {row['n_clusters_found']} | "
              f"{row['rss_peak_bytes'] / 1e6:.0f} | {load_text} | "
              f"{matrix_mb:.0f} |")


if __name__ == "__main__":
    main()
