"""Benchmark streamed k-medoids against the matrix path, and price a comparison.

Both paths start from the same synthetic library and the same fingerprint
distances: the streaming path clusters a :class:`FingerprintComparison` and
holds only an ``O(n_clusters * item count)`` column cache, the matrix path
clusters a ``SymmetricDistanceMatrix`` that :func:`pdist` materialized first.
Each measurement runs in its own subprocess, because peak RSS is a process
high-water mark that whichever path ran first in a shared process would set for
both.

Input construction is outside the timed region on both paths, mirroring
``linkage_memory.py``: a helper returns only the clustering input, so the
molecule list is gone, and collected, before the call. The timed region is the
clustering call alone, which is what makes the two paths comparable -- the
streamed path's cost is passes over the pairs, not the one-off fingerprinting
that both paths pay. ``rss_inputs_bytes``, read immediately before the call,
carries what building the input peaked at, so the kernel's own growth is
``rss_peak_bytes - rss_inputs_bytes`` rather than a figure the matrix build
dominates.

Every record carries its own provenance -- host, OS release, core count,
requested *and* resolved thread count, oecluster version, chunk size, fixture,
seed, ``n_clusters``, init mode and the load averages taken as the measurement
started -- because a table assembled from rows that cannot say which machine
produced them cannot be re-checked, only re-measured.

**Peak RSS does not depend on machine load; wall time does.** A row whose
1-minute load average is at or above half the core count is still a valid
memory figure and not a valid timing. Each record decides that for itself and
reports it as ``timing_valid``, rather than leaving a reader to pair a load
average with a core count the row does not carry; the printed table marks a
failing row ``TAINTED``.

``--cooldown SECONDS`` sleeps in the parent before each child is spawned. The
1-minute load average has roughly a 60-second time constant, so an all-core
row's contribution decays to about a fifth after 90 seconds; without that pause
the harness grades its own interference as ambient load.

``--probe`` prices one comparison instead of running the grid, which is what
the default comparison budget has to be stated in terms of. See
:func:`probe` for why it is measured two ways.

Example::

    python benchmarks/kmedoids_streaming.py --sizes 2000 5000 --cooldown 90 --json
    python benchmarks/kmedoids_streaming.py --probe --probe-threads 1 --json
"""

from __future__ import annotations

import argparse
import gc
import json
import os
import platform
import statistics
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

# The RSS reader and the synthetic library are shared with the threshold
# benchmark, which lives beside this file; a direct script run already has this
# directory on sys.path, an importlib load does not.
sys.path.insert(0, str(Path(__file__).resolve().parent))
from streaming_threshold import library_smiles, matrix_bytes, peak_rss_bytes

PATHS = ("streaming", "matrix")
INITS = ("build", "farthest_first")

#: Budget handed to every streamed cell, far above anything a grid can spend.
#: Set explicitly rather than left to the native default, which a large enough
#: cell can exhaust mid-run: the swap loop would then stop between passes and
#: the row would time a truncated run while reporting the configuration asked
#: for. Measuring the algorithm means taking the guard off, not measuring the
#: guard. The default is 3e11 comparisons, which the published sizes clear by
#: two orders of magnitude -- but a benchmark that only happens to fit under a
#: limit it never states is one size change away from measuring the limit.
UNCONSTRAINED_COMPARISONS = 10**15


def load_average() -> list[float] | None:
    """Return the 1, 5 and 15-minute load averages, or None where unsupported.

    ``os.getloadavg`` does not exist on Windows.
    """
    getloadavg = getattr(os, "getloadavg", None)
    return list(getloadavg()) if getloadavg is not None else None


def resolved_workers(num_threads: int, n_items: int) -> int | None:
    """Return the worker count a run will actually get, or None if unknowable.

    Recorded alongside the requested ``num_threads`` because the requested
    value is usually 0, meaning "auto-detect", and a row carrying only the 0
    cannot be graded: the timing-validity rule is a load average against half
    the core count, and half of what is exactly the thing 0 withholds. Mirrors
    the native cap ``min(num_threads, n)`` and ``hardware_concurrency()`` for
    the auto case, for which ``os.cpu_count()`` is the usual answer.
    """
    cores = os.cpu_count()
    if num_threads == 0:
        return None if cores is None else min(cores, n_items)
    return min(num_threads, n_items)


def provenance(n_items: int, num_threads: int) -> dict[str, Any]:
    """Return the machine and build facts a row needs to be re-checkable.

    Emitted into every record rather than stated in the prose that quotes it.
    Retrofitting provenance cost a whole re-measurement once already, and a
    row merged into a table from another machine is indistinguishable without
    these: a 1-minute load of 20 fails the validity rule on 14 cores and
    passes it on 64.
    """
    import oecluster

    cores = os.cpu_count()
    load = load_average()
    workers = resolved_workers(num_threads, n_items)
    # None rather than False where either side is unknown: an unknowable gate
    # is not a passed one, and False would read as a measured failure.
    valid: bool | None = None
    if load is not None and cores is not None:
        valid = load[0] < cores / 2
    return {"host": platform.node(), "platform": platform.platform(),
            "cpu_count": cores, "resolved_workers": workers,
            "oecluster_version": oecluster.__version__,
            "load_average": load, "timing_valid": valid}


def sweep_comparisons(n_items: int) -> int:
    """Return the comparisons in one pass over every pair.

    Mirrors ``sweep_comparisons()`` in ``src/clustering/KMedoidsBudget.h``: a
    candidate's scan reads its whole row and nothing is retained between
    candidates, so each unordered pair is read twice and the diagonal is free.
    """
    return 0 if n_items < 2 else n_items * (n_items - 1)


def unavoidable_comparisons(n_items: int, n_clusters: int, init: str) -> int:
    """Return the comparisons a run must make whatever the swap loop does.

    Mirrors ``unavoidable_comparisons()`` in ``KMedoidsBudget.h``, which is
    what the native layer prices a streamed call against before it scores a
    single pair. The harness recomputes it so a cell that the budget would
    refuse is caught in the parent, by name, rather than surfacing as a failed
    child two sizes into a run.

    :param n_items: Item count.
    :param n_clusters: Medoid count.
    :param init: One of :data:`INITS`, or ``"explicit"``.
    :returns: ``(init sweeps + n_clusters) * n * (n - 1) + n * n_clusters``.
    :raises ValueError: If ``init`` is not a known strategy.
    """
    init_sweeps = {"build": n_clusters, "farthest_first": 1, "explicit": 0}
    if init not in init_sweeps:
        raise ValueError(f"unknown init: {init}")
    sweeps = init_sweeps[init] + n_clusters
    return sweeps * sweep_comparisons(n_items) + n_items * n_clusters


def build_molecules(n_items: int, seed: int) -> list[Any]:
    """Return ``n_items`` molecules from the shared synthetic library."""
    from openeye import oechem

    mols = []
    for smiles in library_smiles(n_items, seed):
        mol = oechem.OEGraphMol()
        if not oechem.OESmilesToMol(mol, smiles):
            raise ValueError(f"unparseable library SMILES: {smiles}")
        mols.append(mol)
    return mols


def build_input(path: str, n_items: int, seed: int, num_threads: int) -> Any:
    """Return the clustering input for ``path`` and nothing else.

    Keeping the molecule list local to this function is what lets it be freed
    before the timed region, so the fingerprints the input already holds are
    not counted twice against the peak.

    :param path: ``"streaming"`` for a comparison, ``"matrix"`` for a matrix.
    :param n_items: Item count.
    :param seed: Seed for the library draw.
    :param num_threads: Threads for the matrix fill; ignored when streaming.
    :returns: A ``FingerprintComparison`` or a ``SymmetricDistanceMatrix``.
    """
    import oecluster

    mols = build_molecules(n_items, seed)
    if path == "streaming":
        # Default options on both sides: FingerprintComparison(mols) and
        # pdist(mols, "fingerprint") produce the same distances, so the two
        # paths cluster the same problem rather than two similar ones.
        return oecluster.FingerprintComparison(mols)
    return oecluster.pdist(mols, "fingerprint", num_threads=num_threads)


def measure(spec: dict[str, Any]) -> dict[str, Any]:
    """Run one measurement in this process and return its record.

    :param spec: ``path``, ``n_items``, ``init``, ``n_clusters``, ``fixture``,
        ``seed``, ``num_threads``, ``chunk_size``, ``max_iterations`` and
        ``max_comparisons``.
    :returns: The spec plus ``seconds`` (the clustering call only),
        ``rss_inputs_bytes`` (peak before the call), ``rss_peak_bytes``,
        ``matrix_bytes``, ``n_clusters_found``, ``n_iterations``,
        ``converged``, ``cost`` and everything :func:`provenance` returns.
    """
    import oecluster

    target = build_input(spec["path"], spec["n_items"], spec["seed"],
                         spec["num_threads"])
    gc.collect()
    # Read before the timed call: peak RSS is a high-water mark, and the
    # matrix path already peaks at the pdist fill, so the end-to-end peak
    # alone cannot show what the kernel itself adds.
    rss_inputs = peak_rss_bytes()
    facts = provenance(spec["n_items"], spec["num_threads"])

    options: dict[str, Any] = {}
    if spec["path"] == "streaming":
        # Refused on a matrix, where the parameter does not apply.
        options["max_comparisons"] = spec["max_comparisons"]

    started = time.perf_counter()
    result = oecluster.k_medoids(
        target, n_clusters=spec["n_clusters"], init=spec["init"],
        max_iterations=spec["max_iterations"],
        num_threads=spec["num_threads"], chunk_size=spec["chunk_size"],
        **options)
    seconds = time.perf_counter() - started
    return {**spec, "seconds": seconds,
            "rss_inputs_bytes": rss_inputs,
            "rss_peak_bytes": peak_rss_bytes(),
            "matrix_bytes": matrix_bytes(spec["n_items"]),
            "n_clusters_found": result.num_clusters,
            "n_iterations": result.n_iterations,
            "converged": result.converged,
            "cost": result.cost, **facts, "attempted": True}


def probe(spec: dict[str, Any]) -> dict[str, Any]:
    """Price one comparison, two ways, and return the record.

    The budget counts ``Compare`` calls made inside the kernel, so that is
    what has to be priced, and a Python loop over ``Compare`` prices the SWIG
    round trip instead -- an upper bound worth recording precisely because the
    gap between the two is large enough to change what the default budget
    means.

    The in-kernel figure comes from a difference rather than a whole run. Two
    ``init="explicit"`` runs are timed, one given a budget of exactly the
    unavoidable comparisons and one given a single sweep more. The first
    cannot afford a pass and stops at the top of the swap loop; the second
    affords exactly one predicted-delta scan, which is one sweep over every
    pair. What the subtraction removes is the shared setup -- filling the
    medoid columns and assembling the result -- and at the default probe size
    that is 0.4 ms against a 141 ms sweep, so it is a tidiness rather than a
    protection. It removes nothing that happens *inside* the one sweep: the
    per-pass thread-pool construction and the first touch of the column cache
    are charged to the sweep, which is the conservative direction. So is the
    ``O(item count)`` column update of an accepted swap, two parts in ten
    thousand here.

    The divisor is the **priced** sweep, ``n * (n - 1)``, which is what
    ``max_comparisons`` counts. ``best_predicted_swap`` skips the
    ``n_clusters`` items that are already medoids, so the kernel makes
    slightly fewer calls than it is charged for and the figure is low by
    ``n_clusters / n`` -- 0.1% at the default probe ``n_clusters`` of 2. Per
    budget-counted comparison is the right unit for a number that exists to
    price a budget, but it is not quite per ``Compare`` call.

    :param spec: ``n_items``, ``n_clusters``, ``fixture``, ``seed``,
        ``num_threads``, ``chunk_size``, ``pairs`` and ``repeats``.
    :returns: The spec plus ``python_compare_microseconds``,
        ``kernel_compare_microseconds``, the medians and the full spreads
        behind it, ``sweep``, and everything :func:`provenance` returns.
    """
    import oecluster

    n_items = spec["n_items"]
    comparison = build_input("streaming", n_items, spec["seed"],
                             spec["num_threads"])
    gc.collect()
    facts = provenance(n_items, spec["num_threads"])

    # Strided so consecutive calls do not re-read one row, which would price a
    # warm cache the kernel's row scan does not get either.
    started = time.perf_counter()
    for index in range(spec["pairs"]):
        comparison.Compare(index % n_items, (index * 7 + 1) % n_items)
    python_seconds = time.perf_counter() - started

    k_probe = spec["n_clusters"]
    floor = unavoidable_comparisons(n_items, k_probe, "explicit")
    sweep = sweep_comparisons(n_items)

    def timed(budget: int) -> float:
        started = time.perf_counter()
        oecluster.k_medoids(
            comparison, n_clusters=k_probe, init="explicit",
            initial_medoids=list(range(k_probe)),
            num_threads=spec["num_threads"], chunk_size=spec["chunk_size"],
            max_comparisons=budget)
        return time.perf_counter() - started

    # Interleaved, not one batch then the other: a drifting machine would
    # otherwise land entirely on whichever batch ran during the busy spell and
    # move the difference by the whole drift.
    no_sweep = []
    one_sweep = []
    for _ in range(spec["repeats"]):
        no_sweep.append(timed(floor))
        one_sweep.append(timed(floor + sweep))
    delta = statistics.median(one_sweep) - statistics.median(no_sweep)

    return {**spec, "sweep": sweep,
            "python_compare_microseconds": python_seconds / spec["pairs"] * 1e6,
            "kernel_compare_microseconds": delta / sweep * 1e6,
            "no_sweep_seconds": statistics.median(no_sweep),
            "one_sweep_seconds": statistics.median(one_sweep),
            # Every repeat, not just the medians. This is the keystone number
            # of the measurement and a median alone cannot tell a stable
            # difference from a lucky draw out of a drifting machine.
            "no_sweep_seconds_all": no_sweep,
            "one_sweep_seconds_all": one_sweep,
            "kernel_compare_microseconds_range": [
                (min(one_sweep) - max(no_sweep)) / sweep * 1e6,
                (max(one_sweep) - min(no_sweep)) / sweep * 1e6],
            **facts}


def run_child(mode: str, payload: dict[str, Any]) -> dict[str, Any]:
    """Run ``mode`` in a fresh interpreter and return its record."""
    completed = subprocess.run(
        [sys.executable, __file__, f"--child-{mode}", json.dumps(payload)],
        capture_output=True, text=True, check=False)
    if completed.returncode != 0:
        raise RuntimeError(
            f"benchmark child failed with exit status {completed.returncode} "
            f"for {payload}:\n{completed.stderr}")
    return json.loads(completed.stdout.strip().splitlines()[-1])


def benchmark(args: argparse.Namespace) -> list[dict[str, Any]]:
    """Run every requested measurement and return one record per row.

    :raises ValueError: If a streamed cell's unavoidable comparisons exceed
        ``--max-comparisons``. Shrinking the cell to fit would publish a
        figure for a configuration nobody asked about, so the budget is the
        thing that has to move.
    """
    rows = []
    # Only when a streamed cell will actually run: the matrix path never
    # consults a comparison budget, so refusing --paths matrix over one would
    # block a measurement on a limit it does not obey.
    budgeted = "streaming" in args.paths
    for n_items in args.sizes:
        for init in args.inits:
            floor = unavoidable_comparisons(n_items, args.n_clusters, init)
            if budgeted and floor > args.max_comparisons:
                raise ValueError(
                    f"n_items={n_items} init={init} "
                    f"n_clusters={args.n_clusters} cannot run: its "
                    f"unavoidable {floor} comparisons exceed the "
                    f"--max-comparisons {args.max_comparisons}")
            spec = {"n_items": n_items, "init": init,
                    "n_clusters": args.n_clusters, "fixture": "library",
                    "seed": args.seed, "num_threads": args.num_threads,
                    "chunk_size": args.chunk_size,
                    "max_iterations": args.max_iterations,
                    "max_comparisons": args.max_comparisons}
            for path in args.paths:
                runs = []
                for _ in range(args.repeats):
                    # In the parent, so the sleep is outside the measured
                    # process and lets the previous row's load decay.
                    time.sleep(args.cooldown)
                    runs.append(run_child("measure", {**spec, "path": path}))
                rows.append(summarize(runs))
    return rows


def summarize(runs: list[dict[str, Any]]) -> dict[str, Any]:
    """Collapse the repeats of one cell into a single record."""
    row = dict(runs[0])
    row["seconds"] = statistics.median(r["seconds"] for r in runs)
    # Both readings come from the same run, the one whose growth past the
    # input peak is largest. Maxing them separately would pair a peak from one
    # run with a larger inputs reading from another and report a difference no
    # run measured -- always an understatement, in the number this benchmark
    # exists to establish.
    worst = max(runs, key=lambda r: r["rss_peak_bytes"] - r["rss_inputs_bytes"])
    row["rss_peak_bytes"] = worst["rss_peak_bytes"]
    row["rss_inputs_bytes"] = worst["rss_inputs_bytes"]
    # The heaviest run's load, not the first run's. Consumers gate timing
    # validity on this, and a median drawn partly from a busy run must not be
    # stamped with a quiet one. load_average() is None where os.getloadavg is
    # absent, so keep None rather than indexing it.
    loads = [r["load_average"] for r in runs if r["load_average"] is not None]
    row["load_average"] = max(loads, key=lambda la: la[0]) if loads else None
    # Re-derived from the load that was just substituted in, not inherited
    # from runs[0]: a row stamped with the busiest run's load and the
    # quietest run's verdict would say the opposite of its own evidence.
    cores = row["cpu_count"]
    row["timing_valid"] = (None if row["load_average"] is None or cores is None
                           else row["load_average"][0] < cores / 2)
    return row


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Benchmark streamed k-medoids against the matrix path.")
    parser.add_argument("--sizes", type=int, nargs="+",
                        default=[2000, 5000, 10000])
    parser.add_argument("--paths", nargs="+", choices=PATHS,
                        default=list(PATHS))
    parser.add_argument("--inits", nargs="+", choices=INITS,
                        default=list(INITS))
    parser.add_argument("--n-clusters", type=int, default=20)
    parser.add_argument("--max-iterations", type=int, default=100)
    parser.add_argument("--num-threads", type=int, default=0,
                        help="0 auto-detects hardware concurrency")
    parser.add_argument("--chunk-size", type=int, default=4096,
                        help="items per work unit; the library default is "
                             "4096, and a chunk at or above the item count "
                             "leaves the pass single-threaded")
    parser.add_argument("--max-comparisons", type=int,
                        default=UNCONSTRAINED_COMPARISONS,
                        help="budget for the streamed cells; the default is "
                             "far above the grid so no row times a run the "
                             "guard truncated")
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
    parser.add_argument("--probe", action="store_true",
                        help="price one comparison instead of running the grid")
    parser.add_argument("--probe-size", type=int, default=2000)
    parser.add_argument("--probe-clusters", type=int, default=2)
    parser.add_argument("--probe-pairs", type=int, default=200000)
    parser.add_argument("--probe-repeats", type=int, default=5)
    parser.add_argument("--json", action="store_true", dest="as_json")
    parser.add_argument("--child-measure", help=argparse.SUPPRESS)
    parser.add_argument("--child-probe", help=argparse.SUPPRESS)
    return parser.parse_args(argv)


def load_text(record: dict[str, Any]) -> str:
    """Return the 1-minute load with its verdict, for a table cell.

    The verdict travels with the number because the number alone does not
    carry one: the rule is a load against half the core count, so the same
    20 fails on 14 cores and passes on 64.
    """
    load = record["load_average"]
    if load is None:
        return "n/a"
    mark = {True: "ok", False: "TAINTED", None: "?"}[record["timing_valid"]]
    return f"{load[0]:.1f} {mark}"


def print_probe(record: dict[str, Any]) -> None:
    """Print the probe record as text."""
    low, high = record["kernel_compare_microseconds_range"]
    print(f"host {record['host']}  n {record['n_items']}  "
          f"fixture {record['fixture']}  "
          f"threads {record['num_threads']} "
          f"(resolved {record['resolved_workers']} of "
          f"{record['cpu_count']})  chunk {record['chunk_size']}  "
          f"load 1m {load_text(record)}")
    print(f"python Compare: {record['python_compare_microseconds']:.4f} us "
          f"per call over {record['pairs']} calls")
    print(f"in-kernel:      {record['kernel_compare_microseconds']:.4f} us "
          f"per budget-counted comparison over a {record['sweep']}-comparison "
          f"sweep")
    print(f"                worst-case spread over {record['repeats']} "
          f"repeats: {low:.4f} to {high:.4f} us")


def print_rows(rows: list[dict[str, Any]]) -> None:
    """Print the grid as a table."""
    print("| path | n | init | k | seconds | iters | clusters | peak RSS MB | "
          "inputs RSS MB | matrix MB | load 1m |")
    print("| --- | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | "
          "---: | --- |")
    for row in rows:
        print(f"| {row['path']} | {row['n_items']} | {row['init']} | "
              f"{row['n_clusters']} | {row['seconds']:.2f} | "
              f"{row['n_iterations']} | {row['n_clusters_found']} | "
              f"{row['rss_peak_bytes'] / 1e6:.0f} | "
              f"{row['rss_inputs_bytes'] / 1e6:.0f} | "
              f"{row['matrix_bytes'] / 1e6:.0f} | {load_text(row)} |")
    if any(row["timing_valid"] is not True for row in rows):
        print()
        print("TAINTED: the 1-minute load average was at or above half the "
              "core count when this row ran. Its peak RSS stands; its "
              "seconds are not a wall-clock claim.")


def main(argv: list[str] | None = None) -> None:
    """Run the benchmark or the probe, printing a table or ``--json``."""
    args = parse_args(argv)
    if args.child_measure is not None:
        print(json.dumps(measure(json.loads(args.child_measure))))
        return
    if args.child_probe is not None:
        print(json.dumps(probe(json.loads(args.child_probe))))
        return

    if args.probe:
        time.sleep(args.cooldown)
        record = run_child("probe", {
            "n_items": args.probe_size, "n_clusters": args.probe_clusters,
            "fixture": "library", "seed": args.seed,
            "num_threads": args.num_threads,
            "chunk_size": args.chunk_size, "pairs": args.probe_pairs,
            "repeats": args.probe_repeats})
        if args.as_json:
            print(json.dumps(record, indent=2))
        else:
            print_probe(record)
        return

    rows = benchmark(args)
    if args.as_json:
        print(json.dumps(rows, indent=2))
        return
    print_rows(rows)


if __name__ == "__main__":
    main()
