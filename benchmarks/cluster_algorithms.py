"""Benchmark native clustering against scikit-learn precomputed clustering.

The benchmark isolates clustering time: the synthetic distance matrix and
oecluster storage are built once, then each algorithm is timed repeatedly.
"""

from __future__ import annotations

import argparse
import json
import statistics
import time
from collections.abc import Callable
from dataclasses import asdict
from dataclasses import dataclass
from typing import Any

import numpy as np


@dataclass(frozen=True)
class BenchmarkResult:
    """Single benchmark timing summary."""

    algorithm: str
    n_samples: int
    native_seconds: float
    sklearn_seconds: float

    @property
    def speedup(self) -> float:
        """Return scikit-learn time divided by native time."""
        return self.sklearn_seconds / self.native_seconds


def make_clustered_points(n_samples: int, n_clusters: int, seed: int) -> np.ndarray:
    """Create deterministic two-dimensional clusters for precomputed distances."""
    rng = np.random.default_rng(seed)
    angles = np.linspace(0.0, 2.0 * np.pi, n_clusters, endpoint=False)
    centers = np.column_stack((np.cos(angles), np.sin(angles))) * 10.0
    assignments = np.arange(n_samples) % n_clusters
    jitter = rng.normal(loc=0.0, scale=0.2, size=(n_samples, 2))
    points = centers[assignments] + jitter
    order = np.lexsort((points[:, 1], points[:, 0]))
    return points[order].astype(np.float64)


def square_distances(points: np.ndarray) -> np.ndarray:
    """Return a dense square distance matrix."""
    diff = points[:, None, :] - points[None, :, :]
    return np.sqrt(np.sum(diff * diff, axis=2))


def dense_distance_matrix(square: np.ndarray):
    """Build an oecluster DistanceMatrix from a square matrix."""
    import oecluster

    storage = oecluster.DenseStorage(square.shape[0])
    for i in range(square.shape[0]):
        for j in range(i + 1, square.shape[0]):
            storage.Set(i, j, float(square[i, j]))
    return oecluster.DistanceMatrix(
        storage,
        "synthetic",
        [f"item_{i}" for i in range(square.shape[0])],
        {},
    )


def cluster_members(labels: np.ndarray) -> tuple[tuple[int, ...], ...]:
    """Normalize labels into sorted member tuples, ignoring noise."""
    clusters: list[tuple[int, ...]] = []
    for label in sorted(label for label in set(labels.tolist()) if label >= 0):
        members = tuple(int(index) for index in np.flatnonzero(labels == label))
        clusters.append(members)
    return tuple(sorted(clusters))


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


def check_dbscan_parity(square: np.ndarray, distance_matrix, eps: float, min_samples: int) -> None:
    """Verify native DBSCAN cluster membership against scikit-learn."""
    import oecluster
    from sklearn.cluster import dbscan as sklearn_dbscan

    expected_core, expected_labels = sklearn_dbscan(
        square,
        eps=eps,
        min_samples=min_samples,
        metric="precomputed",
    )
    observed = oecluster.dbscan(
        distance_matrix,
        eps=eps,
        min_samples=min_samples,
    )
    if tuple(int(i) for i in expected_core) != observed.core_sample_indices:
        raise AssertionError("DBSCAN core sample indices differ from scikit-learn")
    if cluster_members(observed.labels) != cluster_members(expected_labels):
        raise AssertionError("DBSCAN cluster memberships differ from scikit-learn")


def check_hdbscan_parity(
    square: np.ndarray,
    distance_matrix,
    min_cluster_size: int,
    min_samples: int,
) -> None:
    """Verify native HDBSCAN cluster membership against scikit-learn."""
    import oecluster
    from sklearn.cluster import HDBSCAN as SklearnHDBSCAN

    expected = SklearnHDBSCAN(
        min_cluster_size=min_cluster_size,
        min_samples=min_samples,
        metric="precomputed",
        cluster_selection_method="eom",
        copy=True,
    ).fit(square)
    observed = oecluster.hdbscan(
        distance_matrix,
        min_cluster_size=min_cluster_size,
        min_samples=min_samples,
        cluster_selection_method="eom",
    )
    if cluster_members(observed.labels) != cluster_members(expected.labels_):
        raise AssertionError("HDBSCAN cluster memberships differ from scikit-learn")
    # Larger synthetic sets can contain equal mutual-reachability ties. Cluster
    # membership is the stable parity contract; probabilities must remain valid.
    if observed.probabilities.shape != expected.probabilities_.shape:
        raise AssertionError("HDBSCAN probability shape differs from scikit-learn")
    if not np.all((observed.probabilities >= 0.0) & (observed.probabilities <= 1.0)):
        raise AssertionError("HDBSCAN probabilities are outside [0, 1]")


def benchmark_dbscan(
    square: np.ndarray,
    distance_matrix,
    repeats: int,
    warmups: int,
    num_threads: int,
) -> BenchmarkResult:
    """Benchmark native DBSCAN and scikit-learn DBSCAN."""
    import oecluster
    from sklearn.cluster import dbscan as sklearn_dbscan

    eps = 0.2
    min_samples = 5
    native_seconds = time_call(
        lambda: oecluster.dbscan(
            distance_matrix,
            eps=eps,
            min_samples=min_samples,
            num_threads=num_threads,
        ),
        repeats,
        warmups,
    )
    sklearn_seconds = time_call(
        lambda: sklearn_dbscan(
            square,
            eps=eps,
            min_samples=min_samples,
            metric="precomputed",
            n_jobs=None if num_threads == 0 else num_threads,
        ),
        repeats,
        warmups,
    )
    return BenchmarkResult("DBSCAN", square.shape[0], native_seconds, sklearn_seconds)


def benchmark_hdbscan(
    square: np.ndarray,
    distance_matrix,
    repeats: int,
    warmups: int,
    num_threads: int,
) -> BenchmarkResult:
    """Benchmark native HDBSCAN and scikit-learn HDBSCAN."""
    import oecluster
    from sklearn.cluster import HDBSCAN as SklearnHDBSCAN

    min_cluster_size = 5
    min_samples = 5
    native_seconds = time_call(
        lambda: oecluster.hdbscan(
            distance_matrix,
            min_cluster_size=min_cluster_size,
            min_samples=min_samples,
            cluster_selection_method="eom",
            num_threads=num_threads,
        ),
        repeats,
        warmups,
    )
    sklearn_seconds = time_call(
        lambda: SklearnHDBSCAN(
            min_cluster_size=min_cluster_size,
            min_samples=min_samples,
            metric="precomputed",
            cluster_selection_method="eom",
            n_jobs=None if num_threads == 0 else num_threads,
            copy=True,
        ).fit(square),
        repeats,
        warmups,
    )
    return BenchmarkResult("HDBSCAN", square.shape[0], native_seconds, sklearn_seconds)


def parse_args() -> argparse.Namespace:
    """Parse command-line options."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sizes", type=int, nargs="+", default=[250, 500, 1000])
    parser.add_argument("--clusters", type=int, default=10)
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--warmups", type=int, default=1)
    parser.add_argument("--seed", type=int, default=7)
    parser.add_argument("--num-threads", type=int, default=0)
    parser.add_argument("--skip-parity", action="store_true")
    parser.add_argument("--json", action="store_true", dest="as_json")
    return parser.parse_args()


def main() -> None:
    """Run benchmarks and print timing results."""
    args = parse_args()
    results = []

    for n_samples in args.sizes:
        points = make_clustered_points(n_samples, args.clusters, args.seed)
        square = square_distances(points)
        distance_matrix = dense_distance_matrix(square)

        if not args.skip_parity:
            check_dbscan_parity(square, distance_matrix, eps=0.2, min_samples=5)
            check_hdbscan_parity(
                square,
                distance_matrix,
                min_cluster_size=5,
                min_samples=5,
            )

        results.append(
            benchmark_dbscan(
                square,
                distance_matrix,
                args.repeats,
                args.warmups,
                args.num_threads,
            )
        )
        results.append(
            benchmark_hdbscan(
                square,
                distance_matrix,
                args.repeats,
                args.warmups,
                args.num_threads,
            )
        )

    if args.as_json:
        payload = [asdict(result) | {"speedup": result.speedup} for result in results]
        print(json.dumps(payload, indent=2))
        return

    print("| algorithm | n | native_s | sklearn_s | speedup |")
    print("| --- | ---: | ---: | ---: | ---: |")
    for result in results:
        print(
            f"| {result.algorithm} | {result.n_samples} | "
            f"{result.native_seconds:.6f} | {result.sklearn_seconds:.6f} | "
            f"{result.speedup:.2f}x |"
        )


if __name__ == "__main__":
    main()
