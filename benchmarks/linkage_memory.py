"""Benchmark the memory and time of complete, average and weighted linkage.

Each measurement runs in its own subprocess, because peak RSS is a process
high-water mark that whichever measurement ran first in a shared process would
set for all the others. The condensed distance array is built directly, one row
at a time, and never as a square matrix: a square ``N x N`` float64 array is
twice the condensed size and would count against the peak before clustering
starts. The matrix is built inside a helper that returns only the
``SymmetricDistanceMatrix``, so the NumPy source array is gone, and collected,
before the timed region; otherwise the source array, the native storage and the
algorithm's workspace would be three quadratic buffers at once. Every fixture is
a metric, because ``from_condensed`` probes the triangle inequality and
``agglomerative`` refuses a non-metric matrix. Every record carries the
1-minute, 5-minute and 15-minute load averages taken as it started.

``--cooldown SECONDS`` sleeps in the parent before each child is spawned. The
1-minute load average has roughly a 60-second time constant, so an all-core
row's contribution decays to about a fifth after 90 seconds; without that pause
the harness grades its own interference as ambient load.

Example::

    python benchmarks/linkage_memory.py --sizes 2000 5000 \\
        --linkages complete average weighted --fixtures random --json
"""

from __future__ import annotations

import argparse
import gc
import json
import os
import statistics
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

# The RSS reader is shared with the threshold benchmark, which lives beside this
# file; a direct script run already has this directory on sys.path, an
# importlib load does not.
sys.path.insert(0, str(Path(__file__).resolve().parent))
from streaming_threshold import matrix_bytes, peak_rss_bytes

LINKAGES = ("complete", "average", "weighted")
FIXTURES = ("random", "all-equal", "blocked", "duplicates")
BLOCKS = 20
DUPLICATE_GROUP = 8


def load_average() -> list[float] | None:
    """Return the 1, 5 and 15-minute load averages, or None where unsupported.

    ``os.getloadavg`` does not exist on Windows.
    """
    getloadavg = getattr(os, "getloadavg", None)
    return list(getloadavg()) if getloadavg is not None else None


def build_condensed(fixture: str, n_items: int, seed: int) -> Any:
    """Fill a condensed distance array for ``fixture`` without a square matrix.

    Rows are filled in scipy's upper-triangle row-major order. Point-based
    fixtures compute one row of Euclidean distances at a time, so the working
    set beyond the output is O(N).

    :param fixture: One of :data:`FIXTURES`.
    :param n_items: Number of items.
    :param seed: Seed for the point generator.
    :returns: A 1-D float64 array of ``n_items * (n_items - 1) / 2`` distances.
    """
    import numpy as np

    out = np.empty(n_items * (n_items - 1) // 2, dtype=np.float64)
    rng = np.random.default_rng(seed)
    if fixture == "random":
        points = rng.random((n_items, 2))
    elif fixture == "duplicates":
        groups = -(-n_items // DUPLICATE_GROUP)
        points = np.repeat(rng.random((groups, 2)), DUPLICATE_GROUP,
                           axis=0)[:n_items]
    elif fixture in ("all-equal", "blocked"):
        points = None
    else:
        raise ValueError(f"unknown fixture: {fixture}")
    block = (np.arange(n_items) * BLOCKS) // n_items

    offset = 0
    for i in range(n_items - 1):
        width = n_items - i - 1
        if points is not None:
            delta = points[i + 1:] - points[i]
            out[offset:offset + width] = np.sqrt((delta * delta).sum(axis=1))
        elif fixture == "all-equal":
            out[offset:offset + width] = 1.0
        else:
            out[offset:offset + width] = (block[i + 1:] != block[i])
        offset += width
    return out


def build_matrix(fixture: str, n_items: int, seed: int) -> Any:
    """Return a ``SymmetricDistanceMatrix`` and nothing else.

    Keeping the source array local to this function is what lets it be freed
    before the timed region.
    """
    import oecluster

    condensed = build_condensed(fixture, n_items, seed)
    return oecluster.SymmetricDistanceMatrix.from_condensed(condensed)


def measure(spec: dict[str, Any]) -> dict[str, Any]:
    """Run one measurement in this process and return its record.

    :param spec: ``linkage``, ``n_items``, ``fixture``, ``seed``,
        ``num_threads`` and ``n_clusters``.
    :returns: The spec plus ``seconds`` (the clustering call only),
        ``rss_peak_bytes``, ``matrix_bytes``, ``n_clusters_found`` and
        ``load_average``.
    """
    import oecluster

    dm = build_matrix(spec["fixture"], spec["n_items"], spec["seed"])
    gc.collect()
    load = load_average()

    started = time.perf_counter()
    result = oecluster.agglomerative(
        dm, linkage=spec["linkage"], n_clusters=spec["n_clusters"],
        num_threads=spec["num_threads"])
    seconds = time.perf_counter() - started
    return {**spec, "seconds": seconds,
            "rss_peak_bytes": peak_rss_bytes(),
            "matrix_bytes": matrix_bytes(spec["n_items"]),
            "n_clusters_found": result.num_clusters,
            "load_average": load, "attempted": True}


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
    for n_items in args.sizes:
        for fixture in args.fixtures:
            for linkage in args.linkages:
                for num_threads in args.threads:
                    spec = {"linkage": linkage, "n_items": n_items,
                            "fixture": fixture, "seed": args.seed,
                            "num_threads": num_threads,
                            "n_clusters": args.n_clusters}
                    runs = []
                    for _ in range(args.repeats):
                        # In the parent, so the sleep is outside the measured
                        # process and lets the previous row's load decay.
                        time.sleep(args.cooldown)
                        runs.append(measure_in_subprocess(spec))
                    row = dict(runs[0])
                    row["seconds"] = statistics.median(
                        r["seconds"] for r in runs)
                    row["rss_peak_bytes"] = max(
                        r["rss_peak_bytes"] for r in runs)
                    rows.append(row)
    return rows


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Benchmark the memory and time of complete, average and "
                    "weighted linkage.")
    parser.add_argument("--sizes", type=int, nargs="+",
                        default=[500, 2000, 5000])
    parser.add_argument("--linkages", nargs="+", choices=LINKAGES,
                        default=list(LINKAGES))
    parser.add_argument("--fixtures", nargs="+", choices=FIXTURES,
                        default=["random"])
    parser.add_argument("--threads", type=int, nargs="+", default=[0],
                        help="num_threads values; 0 is the library default")
    parser.add_argument("--n-clusters", type=int, default=50)
    parser.add_argument("--repeats", type=int, default=1)
    parser.add_argument("--seed", type=int, default=7)
    parser.add_argument("--cooldown", type=float, default=0.0,
                        metavar="SECONDS",
                        help="sleep this long in the parent before each "
                             "measurement; the 1-minute load average has a "
                             "~60 s time constant, so an all-core row's load "
                             "decays to about a fifth after 90 s, and without "
                             "it the harness grades its own interference as "
                             "ambient load")
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
    print("| linkage | n | fixture | threads | seconds | clusters | "
          "peak RSS MB | load 1m | matrix MB |")
    print("| --- | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: |")
    for row in rows:
        load = row["load_average"]
        load_text = f"{load[0]:.1f}" if load is not None else "n/a"
        print(f"| {row['linkage']} | {row['n_items']} | {row['fixture']} | "
              f"{row['num_threads']} | {row['seconds']:.2f} | "
              f"{row['n_clusters_found']} | "
              f"{row['rss_peak_bytes'] / 1e6:.0f} | {load_text} | "
              f"{row['matrix_bytes'] / 1e6:.0f} |")


if __name__ == "__main__":
    main()
