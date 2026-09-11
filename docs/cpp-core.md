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

### Internal validity indices

`ClusterReport` carries seven internal cluster-validity indices, all `double`
and all NaN when undefined. Five are always computed:
`calinski_harabasz_medoid` (higher is better), `davies_bouldin_medoid` (lower is
better), `dunn_mean_separation_mean_diameter` and
`dunn_medoid_separation_medoid_spread` (higher is better), and `point_biserial`
(higher is better, positive meaning that between-cluster pairs are the more
distant). `dunn_medoid_separation_medoid_spread` ignores
`representative_method` and always uses the true medoid, so a report requested
with `representative_method="minimax"` still reports medoid-based values here.

> These are **medoid-substituted** indices. The published Calinski-Harabasz and
> Davies-Bouldin definitions use centroids, which do not exist for a distance
> matrix; each cluster's medoid stands in for its centroid, and the global
> medoid stands in for the grand mean. The values are therefore not comparable
> with published figures or with scikit-learn's. Both ignore
> `representative_method` and always use the true medoid, so a report requested
> with `representative_method="minimax"` still reports medoid-based values here.

The remaining two, `c_index` (lower is better) and `baker_hubert_gamma` (higher
is better), are computed only under `ClusterReportOptions::compute_pair_rank_indices`.

### Optional stages and their cost

`ClusterReportOptions` carries two flags, both `false` by default:

- `compute_pair_rank_indices` fills `c_index` and `baker_hubert_gamma`. The
  stage materialises every pairwise distance among the `Nc` clustered points as
  two sortable arrays, `Nc * (Nc - 1) / 2` doubles in total -- roughly 400 MB at
  `Nc = 10,000` and 10 GB at `Nc = 50,000` -- which is why it is opt-in. One
  flag covers both indices because both are read off the same sorted arrays.
- `compute_per_cluster_records` fills `ClusterReport::records`, a
  `std::vector<ClusterRecord>` holding one row per cluster in member-list order:
  label, size, representative, intra-distance mean and median, radius, diameter,
  mean representative distance, nearest cluster and its distance, silhouette,
  and boundary-violation count. It buffers the largest cluster's `n(n-1)/2`
  pairwise distances to take their median, and `detail::median_distance` copies
  that buffer to sort it, so roughly 400 MB for the buffer and 400 MB again for
  the copy, transiently, at `n = 10,000`. A record's `boundary_violations`
  counts pairs involving that cluster, so the sum over records is twice the
  report's own count, which counts each pair once.

`ClusterReport::noise_coverage_at` is the coverage curve restricted to the noise
points. Its length matches `coverage_at`: the threshold count when there is at
least one cluster, and zero when there is none. Every entry is NaN when the
clustering has no noise, rather than 0.0, which would read as "no noise point is
covered" instead of "there is nothing to cover".

`ClusterReport::requested` is a `ClusterReportRequested` recording which of the
two optional computations the caller asked for. It records the request and not
the outcome, so a NaN can be read unambiguously: false means nobody asked, true
with NaN means asked and undefined. The Python comparison table carries the same
distinction into its cells, rendering an unasked cell as `None` -- printed as
`--` -- and reserving `nan` for asked and undefined. A coverage cell is also
`None` when the report carries the threshold but answered nothing at it, as when
the clustering has no clusters and both coverage curves come back empty.

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
