"""Benchmark Butina and DBSCAN from a comparison against the matrix path.

Both paths start from the same molecules and the same fingerprint comparison:
the streaming path passes the molecules with ``comparison="fingerprint"`` and
holds only the threshold graph, the matrix path materializes ``pdist()`` first.
Each measurement runs in its own subprocess, because peak RSS is a process
high-water mark that whichever path ran first in a shared process would set
for both. At a size where the matrix would not fit, the matrix path is not
attempted; its would-be size is reported instead.

Example::

    python benchmarks/streaming_threshold.py --sizes 20000 \\
        --streaming-only-sizes 100000 --json
"""

from __future__ import annotations

import argparse
import json
import random
import statistics
import subprocess
import sys
import time
from typing import Any

#: Three attachment points each. Substituents use ring-closure digit 2 only
#: and the scaffolds 1 and 3, so a substituent can never close a scaffold
#: ring that is still open.
SCAFFOLDS = (
    "c1cc({a})cc({b})c1{c}",
    "c1nc({a})cc({b})c1{c}",
    "C1CC({a})CC({b})C1{c}",
    "c1cc({a})c3cc({b})ccc3c1{c}",
    "O=C(c1cc({a})ccc1{b})N{c}",
    "c1cc({a})c3nc({b})[nH]c3c1{c}",
    "C1CN({a})CC({b})N1C(=O){c}",
    "c1sc({a})cc1C(=O)N({b}){c}",
    "C1N({a})CC({b})C1{c}",
    "c1cc({a})oc1C({b}){c}",
)

SUBSTITUENTS = (
    "C", "CC", "CCC", "C(C)C", "O", "OC", "OCC", "N", "NC", "N(C)C",
    "F", "Cl", "Br", "C#N", "C(=O)O", "C(=O)OC", "C(=O)N", "C(F)(F)F",
    "S(=O)(=O)C", "c2ccccc2", "c2ccncc2", "C2CC2", "C2CCCC2", "OC(=O)C",
    "NC(=O)C", "CO", "CCO", "CN", "SC", "[N+](=O)[O-]",
)

#: Thresholds of the form k / 10000 with k coprime to 10: a Tanimoto distance
#: on a 2048-bit fingerprint is a ratio with a denominator of at most 2048,
#: so it never sits within the 1e-12 by which pdist and Compare can disagree,
#: and the two paths cluster identically.
DEFAULT_THRESHOLDS = (0.3109, 0.4519)


def library_smiles(n_items: int, seed: int) -> list[str]:
    """Return a deterministic synthetic library of ``n_items`` SMILES.

    Molecules sharing a scaffold fall close together in fingerprint space, so
    the library is clustered by construction. Items are drawn without
    replacement from every scaffold and substituent combination.

    :param n_items: Number of molecules.
    :param seed: Seed for the draw.
    :returns: SMILES strings in ascending combination order.
    :raises ValueError: If more molecules are asked for than combinations exist.
    """
    width = len(SUBSTITUENTS)
    combinations = len(SCAFFOLDS) * width**3
    if n_items > combinations:
        raise ValueError(
            f"the library has {combinations} combinations; asked for {n_items}")
    smiles = []
    for index in sorted(random.Random(seed).sample(range(combinations), n_items)):
        rest, c = divmod(index, width)
        rest, b = divmod(rest, width)
        scaffold, a = divmod(rest, width)
        smiles.append(SCAFFOLDS[scaffold].format(
            a=SUBSTITUENTS[a], b=SUBSTITUENTS[b], c=SUBSTITUENTS[c]))
    return smiles


def peak_rss_bytes() -> int:
    """Return this process's peak resident set size in bytes.

    ``ru_maxrss`` is in bytes on macOS and in KiB on Linux. Imported here
    because ``resource`` is POSIX-only, and the library generator above must
    stay importable everywhere.
    """
    import resource

    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return int(peak) if sys.platform == "darwin" else int(peak) * 1024


def matrix_bytes(n_items: int) -> int:
    """Return the size of the condensed float64 matrix over ``n_items``."""
    return 8 * (n_items * (n_items - 1) // 2)


def measure(spec: dict[str, Any]) -> dict[str, Any]:
    """Run one measurement in this process and return its record.

    :param spec: ``algorithm``, ``path``, ``n_items``, ``threshold``,
        ``seed``, ``num_threads`` and ``min_samples``.
    :returns: The spec plus ``seconds``, ``n_clusters``, ``rss_inputs_bytes``
        and ``rss_peak_bytes``.
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

    cluster = oecluster.butina if spec["algorithm"] == "butina" else oecluster.dbscan
    options: dict[str, Any] = {"num_threads": spec["num_threads"]}
    if spec["algorithm"] == "dbscan":
        options["min_samples"] = spec["min_samples"]

    started = time.perf_counter()
    if spec["path"] == "streaming":
        result = cluster(mols, spec["threshold"], comparison="fingerprint",
                         **options)
    else:
        matrix = oecluster.pdist(mols, "fingerprint",
                                 num_threads=spec["num_threads"])
        result = cluster(matrix, spec["threshold"], **options)
    seconds = time.perf_counter() - started
    return {**spec, "seconds": seconds, "n_clusters": result.num_clusters,
            "rss_inputs_bytes": rss_inputs, "rss_peak_bytes": peak_rss_bytes()}


def measure_in_subprocess(spec: dict[str, Any]) -> dict[str, Any]:
    """Run :func:`measure` in a fresh interpreter and return its record."""
    completed = subprocess.run(
        [sys.executable, __file__, "--child", json.dumps(spec)],
        capture_output=True, text=True, check=True)
    return json.loads(completed.stdout.strip().splitlines()[-1])


def benchmark(args: argparse.Namespace) -> list[dict[str, Any]]:
    """Run every requested measurement and return one record per row."""
    rows = []
    sizes = [(n, True) for n in args.sizes]
    sizes += [(n, False) for n in args.streaming_only_sizes]
    for n_items, matrix_fits in sizes:
        for algorithm in ("butina", "dbscan"):
            for threshold in args.thresholds:
                spec = {"algorithm": algorithm, "n_items": n_items,
                        "threshold": threshold, "seed": args.seed,
                        "num_threads": args.num_threads,
                        "min_samples": args.min_samples}
                for path in ("streaming", "matrix"):
                    if path == "matrix" and not matrix_fits:
                        rows.append({**spec, "path": path, "attempted": False,
                                     "matrix_bytes": matrix_bytes(n_items)})
                        continue
                    runs = [measure_in_subprocess({**spec, "path": path})
                            for _ in range(args.repeats)]
                    row = dict(runs[0])
                    row["seconds"] = statistics.median(r["seconds"] for r in runs)
                    row["rss_peak_bytes"] = max(r["rss_peak_bytes"] for r in runs)
                    row["attempted"] = True
                    row["matrix_bytes"] = matrix_bytes(n_items)
                    rows.append(row)
    return rows


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Benchmark Butina and DBSCAN from a comparison against "
                    "the matrix path.")
    parser.add_argument("--sizes", type=int, nargs="*", default=[20000],
                        help="sizes at which both paths run")
    parser.add_argument("--streaming-only-sizes", type=int, nargs="*",
                        default=[100000],
                        help="sizes at which the matrix path is not attempted")
    parser.add_argument("--thresholds", type=float, nargs="+",
                        default=list(DEFAULT_THRESHOLDS))
    parser.add_argument("--min-samples", type=int, default=5)
    parser.add_argument("--repeats", type=int, default=1)
    parser.add_argument("--seed", type=int, default=7)
    parser.add_argument("--num-threads", type=int, default=0)
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
    print("| algorithm | n | threshold | path | seconds | clusters | "
          "peak RSS MB | inputs RSS MB | matrix MB |")
    print("| --- | ---: | ---: | --- | ---: | ---: | ---: | ---: | ---: |")
    for row in rows:
        matrix_mb = row["matrix_bytes"] / 1e6
        if not row["attempted"]:
            print(f"| {row['algorithm']} | {row['n_items']} | "
                  f"{row['threshold']} | {row['path']} | not attempted | | | | "
                  f"{matrix_mb:.0f} |")
            continue
        print(f"| {row['algorithm']} | {row['n_items']} | {row['threshold']} | "
              f"{row['path']} | {row['seconds']:.2f} | {row['n_clusters']} | "
              f"{row['rss_peak_bytes'] / 1e6:.0f} | "
              f"{row['rss_inputs_bytes'] / 1e6:.0f} | {matrix_mb:.0f} |")


if __name__ == "__main__":
    main()
