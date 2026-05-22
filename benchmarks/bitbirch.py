"""Benchmark native oecluster BitBirch against local Python BitBirch."""

from __future__ import annotations

import argparse
import importlib
import json
import statistics
import sys
import time
from collections.abc import Callable
from dataclasses import asdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import oecluster
import oefp


BITBIRCH_REPO = Path("/Users/johnss51/Development/python/bitbirch")
ISIM_REPO = Path("/Users/johnss51/Development/python/iSIM")


@dataclass(frozen=True)
class BenchmarkResult:
    """Single benchmark timing summary."""

    workflow: str
    mode: str
    n_samples: int
    n_bits: int
    density: float
    native_seconds: float
    reference_seconds: float

    @property
    def speedup(self) -> float:
        """Return Python BitBirch time divided by native time."""
        return self.reference_seconds / self.native_seconds


def load_reference_bitbirch():
    """Import the local Python BitBirch reference implementation."""
    for repo in (str(ISIM_REPO), str(BITBIRCH_REPO)):
        if repo not in sys.path:
            sys.path.insert(0, repo)
    module = importlib.import_module("bitbirch.bitbirch")
    return importlib.reload(module)


def make_bits(
    seed: int,
    rows: int,
    cols: int,
    density: float,
    regime: str = "random",
    prototype_count: int = 8,
) -> np.ndarray:
    """Create deterministic dense binary fingerprints."""
    rng = np.random.default_rng(seed)
    if regime == "random":
        bits = (rng.random((rows, cols)) < density).astype(np.uint8)
    elif regime == "duplicate_blocks":
        if prototype_count < 1:
            raise ValueError("prototype_count must be at least one")
        prototypes = (rng.random((prototype_count, cols)) < density).astype(np.uint8)
        prototype_rows = np.arange(rows, dtype=np.intp) % prototype_count
        bits = prototypes[prototype_rows].copy()
    else:
        raise ValueError(f"Unknown BitBirch benchmark regime: {regime!r}")

    empty_rows = np.flatnonzero(bits.sum(axis=1) == 0)
    for row in empty_rows:
        bits[row, row % cols] = 1
    return bits


def batch_from_bits(bits: np.ndarray):
    """Create an OEFP batch from a dense binary matrix."""
    fingerprints = [
        oefp.OEFP.from_on_bits(bits.shape[1], np.flatnonzero(row).astype(int).tolist())
        for row in bits
    ]
    return oefp.OEFPBatch.from_fingerprints(fingerprints)


def reference_result(model, n_items: int) -> tuple[np.ndarray, tuple[tuple[int, ...], ...]]:
    """Return reference labels and clusters in oecluster conventions."""
    labels = model.get_assignments(n_items) - 1
    clusters = tuple(
        tuple(int(member) for member in cluster)
        for cluster in model.get_cluster_mol_ids()
    )
    return labels.astype(np.intp), clusters


def native_result(result) -> tuple[np.ndarray, tuple[tuple[int, ...], ...]]:
    """Return native labels and clusters using a common shape."""
    return np.asarray(result.labels, dtype=np.intp), result.clusters


def assert_same_result(
    observed: tuple[np.ndarray, tuple[tuple[int, ...], ...]],
    expected: tuple[np.ndarray, tuple[tuple[int, ...], ...]],
    workflow: str,
) -> None:
    """Raise if native and reference outputs differ."""
    observed_labels, observed_clusters = observed
    expected_labels, expected_clusters = expected
    if observed_labels.tolist() != expected_labels.tolist():
        raise AssertionError(f"native BitBirch labels differ from reference for {workflow}")
    if observed_clusters != expected_clusters:
        raise AssertionError(f"native BitBirch clusters differ from reference for {workflow}")


def time_call(func: Callable[[], Any], repeats: int, warmups: int) -> float:
    """Return median execution time in seconds."""
    for _ in range(warmups):
        func()

    timings = []
    for _ in range(repeats):
        started = time.perf_counter()
        func()
        timings.append(time.perf_counter() - started)
    return statistics.median(timings)


def reference_cluster(
    bb,
    bits: np.ndarray,
    threshold: float,
    branching_factor: int,
) -> tuple[np.ndarray, tuple[tuple[int, ...], ...]]:
    """Return one-pass reference BitBirch clustering."""
    bb.set_merge("diameter")
    model = bb.BitBirch(threshold=threshold, branching_factor=branching_factor)
    model.fit(bits, singly=True)
    return reference_result(model, bits.shape[0])


def native_cluster(
    batch,
    threshold: float,
    branching_factor: int,
    mode: str,
    num_threads: int,
) -> tuple[np.ndarray, tuple[tuple[int, ...], ...]]:
    """Return one-pass native BitBirch clustering."""
    return native_result(
        oecluster.bitbirch(
            batch,
            threshold=threshold,
            branching_factor=branching_factor,
            merge_criterion="diameter",
            mode=mode,
            num_threads=num_threads,
        )
    )


def reference_recluster(
    bb,
    bits: np.ndarray,
    threshold: float,
    second_threshold: float,
    second_tolerance: float,
    branching_factor: int,
) -> tuple[np.ndarray, tuple[tuple[int, ...], ...]]:
    """Return reference two-stage BitBirch reclustering."""
    bb.set_merge("diameter")
    first = bb.BitBirch(threshold=threshold, branching_factor=branching_factor)
    first.fit(bits, singly=True)
    data, bigs = first.prepare_data_BFs(bits)
    bb.set_merge("tolerance", tolerance=second_tolerance)
    second = bb.BitBirch(threshold=second_threshold, branching_factor=branching_factor)
    second.fit_BFs(data)
    second.fit_BFs(bigs)
    return reference_result(second, bits.shape[0])


def native_recluster(
    batch,
    threshold: float,
    second_threshold: float,
    second_tolerance: float,
    branching_factor: int,
    mode: str,
    num_threads: int,
) -> tuple[np.ndarray, tuple[tuple[int, ...], ...]]:
    """Return native two-stage BitBirch reclustering."""
    return native_result(
        oecluster.bitbirch_recluster(
            batch,
            initial_threshold=threshold,
            second_threshold=second_threshold,
            second_tolerance=second_tolerance,
            branching_factor=branching_factor,
            mode=mode,
            num_threads=num_threads,
        )
    )


def reference_refine(
    bb,
    bits: np.ndarray,
    threshold: float,
    branching_factor: int,
    *,
    reassign_top_clusters: int = 0,
    redistribute_largest_cluster: bool = False,
) -> tuple[np.ndarray, tuple[tuple[int, ...], ...]]:
    """Return reference BitBirch refinement output."""
    bb.set_merge("diameter")
    model = bb.BitBirch(threshold=threshold, branching_factor=branching_factor)
    model.fit(bits, singly=False)
    if redistribute_largest_cluster:
        model.prune(bits)
    if reassign_top_clusters:
        model.reassign(bits, top=reassign_top_clusters)
    return reference_result(model, bits.shape[0])


def native_refine(
    batch,
    threshold: float,
    branching_factor: int,
    mode: str,
    num_threads: int,
    *,
    reassign_top_clusters: int = 0,
    redistribute_largest_cluster: bool = False,
) -> tuple[np.ndarray, tuple[tuple[int, ...], ...]]:
    """Return native BitBirch refinement output."""
    return native_result(
        oecluster.bitbirch_refine(
            batch,
            threshold=threshold,
            branching_factor=branching_factor,
            merge_criterion="diameter",
            singly=False,
            mode=mode,
            num_threads=num_threads,
            reassign_top_clusters=reassign_top_clusters,
            redistribute_largest_cluster=redistribute_largest_cluster,
        )
    )


def workflow_calls(
    workflow: str,
    bb,
    bits: np.ndarray,
    batch,
    args: argparse.Namespace,
) -> tuple[Callable[[], Any], Callable[[], Any]]:
    """Return reference and native callables for a benchmark workflow."""
    if workflow == "cluster":
        return (
            lambda: reference_cluster(bb, bits, args.threshold, args.branching_factor),
            lambda: native_cluster(
                batch,
                args.threshold,
                args.branching_factor,
                args.mode,
                args.num_threads,
            ),
        )
    if workflow == "recluster":
        return (
            lambda: reference_recluster(
                bb,
                bits,
                args.threshold,
                args.second_threshold,
                args.second_tolerance,
                args.branching_factor,
            ),
            lambda: native_recluster(
                batch,
                args.threshold,
                args.second_threshold,
                args.second_tolerance,
                args.branching_factor,
                args.mode,
                args.num_threads,
            ),
        )
    if workflow == "reassign":
        return (
            lambda: reference_refine(
                bb,
                bits,
                args.threshold,
                args.branching_factor,
                reassign_top_clusters=args.reassign_top_clusters,
            ),
            lambda: native_refine(
                batch,
                args.threshold,
                args.branching_factor,
                args.mode,
                args.num_threads,
                reassign_top_clusters=args.reassign_top_clusters,
            ),
        )
    if workflow == "prune":
        return (
            lambda: reference_refine(
                bb,
                bits,
                args.threshold,
                args.branching_factor,
                redistribute_largest_cluster=True,
            ),
            lambda: native_refine(
                batch,
                args.threshold,
                args.branching_factor,
                args.mode,
                args.num_threads,
                redistribute_largest_cluster=True,
            ),
        )
    if workflow == "prune_reassign":
        return (
            lambda: reference_refine(
                bb,
                bits,
                args.threshold,
                args.branching_factor,
                reassign_top_clusters=args.reassign_top_clusters,
                redistribute_largest_cluster=True,
            ),
            lambda: native_refine(
                batch,
                args.threshold,
                args.branching_factor,
                args.mode,
                args.num_threads,
                reassign_top_clusters=args.reassign_top_clusters,
                redistribute_largest_cluster=True,
            ),
        )
    raise ValueError(f"Unknown BitBirch benchmark workflow: {workflow!r}")


def benchmark(
    workflow: str,
    bb,
    bits: np.ndarray,
    batch,
    repeats: int,
    warmups: int,
    args: argparse.Namespace,
) -> BenchmarkResult:
    """Benchmark one data set."""
    # The reference stores running sums in the input dtype. Use signed integer
    # bits here so larger duplicate-block benchmarks measure algorithm parity
    # instead of NumPy unsigned overflow.
    reference_bits = bits.astype(np.int64, copy=False)
    reference_call, native_call = workflow_calls(workflow, bb, reference_bits, batch, args)
    expected = reference_call()
    observed = native_call()
    assert_same_result(observed, expected, workflow)

    reference_seconds = time_call(
        reference_call,
        repeats,
        warmups,
    )
    native_seconds = time_call(
        native_call,
        repeats,
        warmups,
    )
    return BenchmarkResult(
        workflow=workflow,
        mode=args.mode,
        n_samples=bits.shape[0],
        n_bits=bits.shape[1],
        density=float(bits.mean()),
        native_seconds=native_seconds,
        reference_seconds=reference_seconds,
    )


def parse_args() -> argparse.Namespace:
    """Parse command-line options."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--workflows",
        nargs="+",
        choices=["cluster", "recluster", "reassign", "prune", "prune_reassign"],
        default=["cluster"],
    )
    parser.add_argument("--mode", choices=["strict_parity", "fast"], default="strict_parity")
    parser.add_argument("--sizes", type=int, nargs="+", default=[250, 500, 1000])
    parser.add_argument("--bits", type=int, default=2048)
    parser.add_argument("--density", type=float, default=0.05)
    parser.add_argument("--regime", choices=["random", "duplicate_blocks"], default="random")
    parser.add_argument("--prototype-count", type=int, default=8)
    parser.add_argument("--threshold", type=float, default=0.65)
    parser.add_argument("--second-threshold", type=float, default=0.7)
    parser.add_argument("--second-tolerance", type=float, default=0.0)
    parser.add_argument("--branching-factor", type=int, default=50)
    parser.add_argument("--reassign-top-clusters", type=int, default=2)
    parser.add_argument("--num-threads", type=int, default=0)
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--warmups", type=int, default=1)
    parser.add_argument("--seed", type=int, default=13)
    parser.add_argument("--json", action="store_true", dest="as_json")
    return parser.parse_args()


def main() -> None:
    """Run benchmarks and print timing results."""
    args = parse_args()
    results = []
    bb = load_reference_bitbirch()
    for n_samples in args.sizes:
        bits = make_bits(
            args.seed,
            n_samples,
            args.bits,
            args.density,
            args.regime,
            args.prototype_count,
        )
        batch = batch_from_bits(bits)
        for workflow in args.workflows:
            results.append(
                benchmark(
                    workflow,
                    bb,
                    bits,
                    batch,
                    args.repeats,
                    args.warmups,
                    args,
                )
            )

    if args.as_json:
        payload = [asdict(result) | {"speedup": result.speedup} for result in results]
        print(json.dumps(payload, indent=2))
        return

    print("| workflow | mode | n | bits | density | native_s | python_s | speedup |")
    print("| :--- | :--- | ---: | ---: | ---: | ---: | ---: | ---: |")
    for result in results:
        print(
            f"| {result.workflow} | {result.mode} | "
            f"{result.n_samples} | {result.n_bits} | {result.density:.4f} | "
            f"{result.native_seconds:.6f} | {result.reference_seconds:.6f} | "
            f"{result.speedup:.2f}x |"
        )


if __name__ == "__main__":
    main()
