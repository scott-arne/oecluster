# Benchmarks

OECluster ships two opt-in benchmark scripts under `benchmarks/`. They are
developer tools for tracking performance against reference implementations, not
part of the normal test suite. Each prints a Markdown timing table by default
and accepts `--json` for machine-readable output.

## Clustering Algorithms

```bash
python benchmarks/cluster_algorithms.py
```

`cluster_algorithms.py` times the native DBSCAN and HDBSCAN implementations
against scikit-learn on precomputed distance matrices. It builds synthetic
two-dimensional clusters once per size, converts them to a dense OECluster
distance matrix, and times each algorithm repeatedly. Unless `--skip-parity` is
given, it first checks that the native labels match scikit-learn's.

Flags:

- `--sizes N [N ...]` — sample counts to benchmark (default `250 500 1000`).
- `--clusters N` — number of synthetic clusters (default `10`).
- `--repeats N` / `--warmups N` — timing repeats and warmup runs.
- `--seed N` — RNG seed for the synthetic data.
- `--num-threads N` — worker threads (`0` = auto).
- `--skip-parity` — skip the scikit-learn label-parity checks.
- `--json` — emit JSON instead of the Markdown table.

This benchmark imports scikit-learn for the reference timings.

## BitBirch

```bash
python benchmarks/bitbirch.py
```

`bitbirch.py` times the native BitBirch workflows against a local reference
BitBirch implementation. It generates random binary fingerprint batches and
benchmarks the selected workflows across sizes.

Flags:

- `--workflows {cluster,recluster,reassign,prune,prune_reassign} [...]` —
  workflows to benchmark (default `cluster`).
- `--mode {strict_parity,fast}` — execution mode. Fast runs the parallel
  partition-merge path for `cluster` and `recluster` (deterministic and
  quality-equivalent to strict); refine always runs `strict_parity`.
- `--sizes N [N ...]` — sample counts (default `250 500 1000`).
- `--bits N` / `--density F` — fingerprint width and on-bit density.
- `--regime {random,duplicate_blocks}` and `--prototype-count N` — control the
  synthetic fingerprint structure.
- `--threshold`, `--second-threshold`, `--second-tolerance`,
  `--branching-factor`, `--reassign-top-clusters` — algorithm parameters.
- `--num-threads N` — worker threads (`0` = auto).
- `--native-only` — skip the reference implementation and time only the native
  path.
- `--repeats N` / `--warmups N` / `--seed N` — timing controls.
- `--compare-modes` — report strict-vs-fast native speedup, quality delta, and
  fast-mode determinism across thread counts.
- `--json` — emit JSON instead of the Markdown table.

The reference comparison loads a local BitBirch checkout; use `--native-only`
when that reference is not available so the script times only the native
implementation.
