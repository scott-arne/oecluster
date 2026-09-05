# C++ Core

This page summarizes the C++ API for users who want to embed OECluster directly
in a C++ application or extend the library. Include the umbrella header in
user-facing code:

```cpp
#include <oecluster/oecluster.h>
```

All public types live in the `OECluster` namespace. The generated, per-symbol
reference is in [](api/cpp).

## Distance Computation

A pairwise distance computation has three pieces: a `PairwiseComparison` that
scores item pairs, a `StorageBackend` that holds the results, and the `pdist()`
engine that drives the work across threads.

```cpp
#include <oecluster/oecluster.h>

std::vector<OEChem::OEMolBase*> mols = /* load molecules */;

OECluster::FingerprintOptions opts;
opts.fp_type = "morgan";
opts.metric = "tanimoto";
OECluster::FingerprintComparison comparison(mols, opts);

OECluster::DenseStorage storage(mols.size());
OECluster::PDistOptions pdist_opts;
pdist_opts.num_threads = 8;
OECluster::pdist(comparison, storage, pdist_opts);

for (size_t i = 0; i < mols.size(); ++i) {
    for (size_t j = i + 1; j < mols.size(); ++j) {
        std::cout << storage.Get(i, j) << "\n";
    }
}
```

`PDistOptions` controls threading (`num_threads`, `0` for auto-detect),
`chunk_size`, a sparse `cutoff`, and an optional progress callback. `CDist.h`
provides the matching cross-distance engine for the rectangular NxM case.

## Pairwise Comparisons

`PairwiseComparison` is the abstract interface implemented by every comparison
method. Implementations own their data and scorer state and must be cloneable so
each worker thread evaluates pairs without sharing scorer internals. A
comparison can optionally provide a batched `TryPDist()`/`TryCDist()` kernel;
returning `false` asks the engine to fall back to per-pair `Compare()` calls.

The concrete comparisons are:

- `FingerprintComparison` — OEFP fingerprint similarity/distance, configured by
  `FingerprintOptions` (`fp_type`, `storage`, `metric`, `numbits`, the
  per-family shape parameters, the per-metric parameters, `use_chirality`,
  `similarity`).
- `ROCSComparison` — ROCS shape/color overlay, configured by `ROCSOptions`
  (`score_type`, `color_ff_type`, `similarity`).
- `SuperposeComparison` — protein superposition and binding-site comparison,
  configured by `SuperposeOptions` (`method`, `score_type`, oeselect
  `predicate`/`ref_predicate`/`fit_predicate`, `similarity`).
- `DescriptorComparison` — distance in descriptor space, standardized only
  under `standardized_euclidean`, `seuclidean` and `mahalanobis`, the three
  metrics that fit a scale from the input; configured by `DescriptorOptions`
  (`sources`, `columns`, `groups`, `metric`, `variances`,
  `inverse_covariance`, `missing`, `p`).
- `RMSDComparison` — coordinate RMSD between poses of one molecule, configured
  by `RMSDOptions` (`overlay`, `automorph`, `heavy_only`).

## Storage Backends

`StorageBackend` is the abstract interface for distance storage. All backends
use scipy-compatible condensed indexing.

- `DenseStorage` — full condensed matrix held in contiguous memory. Best when
  the matrix fits in RAM and most distances are accessed.
- `MMapStorage` — full condensed matrix backed by a memory-mapped file, for
  out-of-core datasets larger than RAM. The file persists and can be reloaded.
- `SparseStorage` — stores only entries below a cutoff, for large sparse
  distance graphs. It cannot serve complete pairwise distances, so workflows
  that need every distance must use dense or memory-mapped storage.

## Clustering

The clustering headers under `oecluster/clustering` provide the algorithms and
their option/result types:

- `Butina.h` — threshold neighbor-count clustering (`ButinaOptions`,
  `ButinaResult`).
- `DBSCAN.h` — density-based clustering (`DBSCANOptions`, `DBSCANResult`).
- `HDBSCAN.h` — hierarchical density clustering (`HDBSCANOptions`,
  `HDBSCANResult`).
- `Agglomerative.h` — bottom-up linkage clustering (`AgglomerativeOptions`,
  `AgglomerativeResult`).
- `BitBirch.h` — Birch-style clustering of binary fingerprint batches
  (`BitBirchOptions`, `BitBirchReclusteringOptions`,
  `BitBirchRefinementOptions`, `BitBirchResult`).

`ClusterTypes.h` defines the shared cluster representation, and `ClusterReport.h`
defines the method-agnostic quality scorecard exposed in Python as
`cluster_report()`/`compare_reports()`.

## Representatives

`Representative.h` provides representative selection over a cluster and distance
matrix: the medoid, minimax (radius) center, highest-neighborhood
(Butina-style) representative, and the weighted medoid. It also computes the
per-representative quality metrics and supports ranked and k-representative
selection.

## Errors

`Error.h` defines the OECluster exception hierarchy. Invalid arguments and
misuse (unknown method names, a sparse matrix where complete distances are
required, and similar) raise descriptive exceptions; the Python layer surfaces
these as `ValueError`/`RuntimeError`.

## Python Use

Python users normally work through the top-level `oecluster` package, which
wraps these C++ classes with OpenEye-molecule input handling, NumPy-friendly
distance matrices, and result objects. The same comparison and option types are
re-exported there for advanced use. See the [Python API](python-api.md) guide.
