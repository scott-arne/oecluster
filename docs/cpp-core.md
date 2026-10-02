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
`chunk_size`, and an optional progress callback; its `cutoff` field is not
read, because a `SparseStorage` applies its own cutoff. `CDist.h`
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
- `MCSComparison` — maximum common substructure, scored as Tanimoto over
  matched bonds with each pair searched in both directions and the larger match
  taken; configured by `MCSOptions` (`search_mode`, `match_level`,
  `max_matches`, `similarity`).

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

Both accessors refuse a pair they cannot address rather than answering about
it. `Get` and `Set` throw `std::out_of_range` for an index at or beyond
`NumSamples()`, and `Set` additionally throws `std::invalid_argument` for an
in-range diagonal, which owns no stored slot: `Get` answers `i == j` from a
shortcut rather than from memory. The range check runs first, so an
out-of-range diagonal is reported as the range error it also is, and it runs
ahead of `SparseStorage`'s cutoff shortcut, so whether a bad call is diagnosed
does not depend on the value it carried.

`pdist()` checks the storage size against the comparison size separately.
`Set`'s own range check cannot see a storage larger than the comparison,
because every index that loop produces is in range for the oversized storage
and the distances would simply land in the wrong slots.

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
- `KMedoids.h` — k-medoids (PAM) clustering over precomputed distances
  (`KMedoidsOptions`, `KMedoidsResult`).
- `BitBirch.h` — Birch-style clustering of binary fingerprint batches
  (`BitBirchOptions`, `BitBirchReclusteringOptions`,
  `BitBirchRefinementOptions`, `BitBirchResult`).
- `MurckoScaffold.h` — Bemis-Murcko scaffold assignment and scaffold-identity
  clustering of molecules (`ScaffoldType`, `MurckoOptions`, `MurckoResult`).
- `DiversitySelection.h` — farthest-first subset selection and the #Circles
  coverage measure (`MaxMinSeed`, `MaxMinStop`, `CirclesMethod`,
  `MaxMinOptions`, `MaxMinSelection`, `CirclesOptions`, `CirclesResult`).
- `SetDiversity.h` — the Vendi score and log-determinant diversity of a set
  (`DiversityKernel`, `VendiOptions`, `VendiResult`, `LogDetOptions`,
  `LogDetResult`).
- `SphereExclusion.h` — sphere-exclusion clustering under a leader, Butina
  or DISE seed order (`SphereOrder`, `SphereAssignment`,
  `SphereExclusionOptions`, `SphereExclusionResult`).

`ClusterTypes.h` defines the shared cluster representation, `ClusterReport.h`
defines the method-agnostic quality scorecard exposed in Python as
`cluster_report()`/`compare_reports()`, `PartitionAgreement.h` defines the
labels-only agreement metrics exposed as
`partition_agreement()`/`scaffold_agreement()`, and `SARCoherence.h` defines
the structure-activity coherence metrics exposed as
`sar_coherence()`/`activity_landscape()`/`modelability()`. `SARCoherence.h`
is the one of the four that reads activity data as well as structure:
`sar_coherence` decomposes an activity vector across a labeling, while
`activity_landscape` and `modelability` sweep a distance matrix, or a
comparison, directly and need no clustering at all. `cluster_report`,
`activity_landscape` and `modelability` each have a `PairwiseComparison&`
overload beside the `StorageBackend` one; see
[Comparison overloads](#comparison-overloads-for-the-reports).

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

> Two of these five are **medoid-substituted** indices. The published
> Calinski-Harabasz and Davies-Bouldin definitions use centroids, which do not
> exist for a distance matrix, so each cluster's medoid stands in for its
> centroid. Calinski-Harabasz also has a grand-mean term, and the global medoid
> stands in for that; Davies-Bouldin has no such term. The values are therefore
> not comparable with published figures or with scikit-learn's. Both ignore
> `representative_method` and always use the true medoid, so a report requested
> with `representative_method="minimax"` still reports medoid-based values here.

The remaining two, `c_index` (lower is better) and `baker_hubert_gamma` (higher
is better), are computed only under `ClusterReportOptions::compute_pair_rank_indices`.

### Optional stages and their cost

`ClusterReportOptions` carries two flags, both `false` by default:

- `compute_pair_rank_indices` fills `c_index` and `baker_hubert_gamma`. Both are
  read off two sorted arrays that between them hold every pairwise distance
  among the `Nc` clustered points, `Nc * (Nc - 1) / 2` doubles -- roughly 400 MB
  at `Nc = 10,000` and 10 GB at `Nc = 50,000` -- and both are built only under
  this flag. Without it no overload holds a
  pair-sized array: `median_intra_distance` is taken by an exact selector
  whose state is fixed size. One flag covers both indices because both are
  read off the same sorted arrays.
- `compute_per_cluster_records` fills `ClusterReport::records`, a
  `std::vector<ClusterRecord>` holding one row per cluster in member-list order:
  label, size, representative, intra-distance mean and median, radius, diameter,
  mean representative distance, nearest cluster and its distance, silhouette,
  and boundary-violation count. A cluster with at most 2^20 pairs has its
  median taken from a buffer of its distances; a larger cluster's median is
  selected by re-walking its pairs, so the stage buffers at most 2^20 doubles
  (8 MB). Under `compute_pair_rank_indices`, which already holds every pair,
  each cluster is buffered whole. A record's `boundary_violations`
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

### SAR coherence

The three functions are four entry points. `sar_coherence` is overloaded on
how the labels arrive, taking `const ClusteringResult&` or `const
std::vector<ClusterLabel>&`, then `const std::vector<double>& activity`.
`activity_landscape` takes `const StorageBackend& storage` and the same
activity vector; `modelability` takes that storage and `const
std::vector<std::string>& activity_classes`. Each also takes a defaulted
options struct by const reference and returns by value. All four raise
`std::invalid_argument` on an empty or mismatched input vector, the two taking
storage also on anything short of complete pairwise distances: `SparseStorage`
derives from `StorageBackend`, so the signature admits what the runtime
refuses.

`SARCoherenceOptions::noise_handling` defaults to `NoiseHandling::Excluded`,
where `PartitionAgreementOptions`, whose header declares the enum, defaults to
`Singletons`. `ActivityLandscapeOptions` carries `distance_threshold` 0.30,
`activity_threshold` 1.0, `rmodi_delta` 0.625 (the RMODI band half-width, in
activity standard deviations) and `num_threads` 0; `ModelabilityOptions`
carries `num_threads` 0 alone, and zero selects the hardware concurrency.

Each result carries `num_samples`, the input length, and `num_scored`, how
much entered the metric. `SARCoherence` adds `num_clusters`, `eta_squared`,
`omega_squared` and `std::vector<ClusterActivity> clusters`: `label`,
`num_scored`, `mean_activity`, `stddev_activity`. `ActivityLandscape` adds
`num_pairs_scored`, `num_cliffs`, `cliff_density`, `num_zero_distance_pairs`,
`max_sali`, `mean_sali`, `rmodi`, `activity_stddev`. `Modelability` adds
`num_classes`, `modi` and `std::vector<ClassConcordance> classes`: `label`,
`num_members`, `fraction_same_class`.

Every `double` above is NaN where undefined, never substituted: `eta_squared`
when `num_scored < 2` or the activity has no variance, `omega_squared` also
when every scored sample is its own cluster, `cliff_density` when
`num_pairs_scored` is 0, `max_sali` and `mean_sali` when no scored pair has a
nonzero distance, `rmodi`, `activity_stddev` and a row's `stddev_activity` when
the `num_scored` each reads is below 2, `modi` when `num_classes < 2`, and a
row's `fraction_same_class` when its class is the only one scored.

### Comparison overloads for the reports

`cluster_report(const ClusteringResult&, PairwiseComparison&, const
ClusterReportOptions&)`, `activity_landscape(PairwiseComparison&, const
std::vector<double>&, const ActivityLandscapeOptions&)` and
`modelability(PairwiseComparison&, const std::vector<std::string>&, const
ModelabilityOptions&)` score a comparison without materializing a matrix,
with one `Clone()` per concurrently running chunk. Each options struct
gains `chunk_size`, default 4096 pairwise distances per unit. The
comparison overloads refuse a zero `chunk_size`; the storage overloads
ignore the field, and Python validates it on every path. A unit is whole
rows; `activity_landscape` never makes one smaller than 64 rows, the floor
its storage path already applies. Memory is O(N + K) plus O(N + C) per
worker, for K clusters, where C is one comparison clone's size: O(N) for an
`MCSComparison`, which deep-copies its molecules, and O(1) for clones that
share their data. Workers are `num_threads`, or the hardware concurrency
when it is 0, capped at N. No term grows with the pair count. A `cluster_report` fill block
holds at most 2^20 distances, or one row if a row is longer.

The comparison contract is the one `maxmin_select`, `knn_graph` and
`leiden` already document: the engines call `Compare(min(i, j), max(i,
j))`, never a self-pair, and `Compare` must return the same value, or
throw the same exception, for the same pair on every call and on every
clone. `cluster_report`'s exact median rereads pairs above its 2^20-value
budget, so a comparison that breaks this voids the exactness. It is a
documented precondition, not a checked one. Given it, each result is
bit-identical to the storage overload over a matrix holding the values
`Compare` returns, except that a zero median is always `+0.0` on both
overloads.

All three comparison overloads refuse a zero `chunk_size` with
`std::invalid_argument`, then
check the comparison's `Facts()` before any `Compare` call.
`activity_landscape` and `modelability` call `validate_comparison_facts`,
refusing a similarity, a nonzero self-distance, `NaNPresent` and
`SubsetScored` with `ComparisonError`. `cluster_report` refuses the first
three only, as its storage overload has no metric gate, and then refuses
`compute_pair_rank_indices` with `std::invalid_argument`. Length
mismatches name "the comparison". Undeclared non-finite or negative
distances are caught per pair as they are scored, with the storage
overloads' messages; only scored pairs are checked. `cluster_report`
reports the earliest bad pair in row order for every thread count;
`activity_landscape` and `modelability` keep their storage overloads'
selection.

### Diversity selection

`maxmin_select` and `circles` are each overloaded on where distances come
from: `const StorageBackend& storage`, which must hold complete dense or
memory-mapped distances, or `PairwiseComparison& comparison`, evaluated
lazily with one `Clone()` per concurrently running chunk. The comparison
overloads compare only the pairs the run needs, always as
`Compare(min(i, j), max(i, j))` and never a self-pair, and refuse a comparison
whose `Facts()` report a similarity, a nonzero self-distance, `NaNPresent` or
`SubsetScored` with `ComparisonError`. Every other refusal is
`std::invalid_argument`.

`MaxMinOptions` carries `count` 0 (no count limit), `threshold` NaN (unset),
`seed_mode` `MaxMinSeed::Index`, `seed` 0, an empty `initial`, `num_threads` 0
and `chunk_size` 256; at least one of `count` and `threshold` must be set.
`MaxMinSeed::Medoid` is accepted by the storage overload only, and a non-empty
`initial` only with the default seed. `MaxMinSelection` returns `indices`,
`pick_distances` (NaN for the seed and every `initial` entry) and a
`MaxMinStop`. `circles` takes the threshold as its own argument, finite and
non-negative, with `CirclesOptions` carrying `method` `CirclesMethod::MaxMin`,
`num_threads` 0 and `chunk_size` 256; `CirclesResult` returns `count`,
`members`, `threshold` and `method`. An explicit `num_threads` is capped at the
item count.

Native code checks finiteness only on the distances it reads, plus a full
scan for the medoid seed, so a non-finite matrix entry the run never reaches
is not reported. The Python matrix path refuses one up front.

Both entry points share `src/clustering/MaxMinKernel.h`, a farthest-first
kernel templated on a row provider. k-medoids' FarthestFirst initialization
runs on the same kernel through a provider that, as before, does not validate
finiteness, and its global-medoid helper moved there too; its output is
unchanged.

### Set diversity scores

`vendi_score` and `logdet_diversity` are overloaded on the same two distance
sources as `maxmin_select`, with the same `ComparisonError` fact refusals.
Every pair is compared once as `Compare(i, j)` with `i < j`. Each distance
becomes a kernel entry through `DiversityKernel::Complement` (`1 - d`, which
refuses d outside [0, 1]) or `DiversityKernel::Laplacian`
(`exp(-d / bandwidth)`, which refuses a negative d). Either kernel refuses a
non-finite d, and every refusal is `std::invalid_argument`.

`VendiOptions` carries `order` 1, `kernel` `Complement`, `bandwidth` NaN
(unset, required finite and positive for `Laplacian`), `max_exact` 2048,
`num_threads` 0 and `chunk_size` 256. `LogDetOptions` swaps `order` for
`ridge` 0.0, which must be finite and non-negative. Options are validated
first, then storage or facts, then the item count, then the `max_exact`
ceiling, on the exact paths only. `VendiResult` returns `score`, `order`,
`size`, `kernel`, `min_eigenvalue` and `negative_mass` (both NaN at order 2).
`LogDetResult` returns `score` (`-infinity` when `K + ridge I` is not
numerically positive definite), `ridge`, `size`, `kernel`, `min_eigenvalue`
and `nonpositive_count`.

Order-2 Vendi sums each row's squared kernel entries in one slot and then
sums the slots in row order. Chunks own whole rows, so the score is
bit-identical for every `num_threads` and `chunk_size`, and between the two
sources.

The exact scores decompose the dense kernel with a private solver,
`src/clustering/SymmetricEigen.h`: Householder tridiagonalization, then
implicit QL, eigenvalues only, single-threaded and deterministic. It is in
the tree rather than taken from LAPACK or Eigen, because the build's
FetchContent cannot reach either from behind the firewall. The solver throws
`std::runtime_error` if an eigenvalue fails to converge. The shared
validators and the clone-leasing comparison runner that E1 and E2 both use
live in `src/clustering/DiversityValidation.h` and
`src/clustering/ChunkedComparisons.h`.

### Sphere exclusion

`sphere_exclusion` is overloaded on a `StorageBackend` and on a
`PairwiseComparison`:
- `SphereOrder::Input` takes centers in index order.
- `SphereOrder::Permutation` takes them in
  `SphereExclusionOptions::permutation` order, which must be a complete
  permutation.
- `SphereOrder::Neighbors` takes them in Butina's descending neighbor-count
  order, with optional `reordering`, and is available on the storage overload
  only.

`SphereAssignment::Nearest` reassigns non-center items to their nearest center
after the centers are fixed. `SphereExclusionResult::Centers()` returns one
center per cluster.

`butina_cluster` runs on the same engine (`src/clustering/SphereExclusionEngine.h`)
and its outputs are unchanged. Under the neighbor order with
`SphereAssignment::First`, `sphere_exclusion` returns exactly what
`butina_cluster` returns.

The comparison overload evaluates lazily through the clone-leasing
`ChunkedComparisons` runner. Workers write only their own distance slots, and
claims are committed on the calling thread. As a result, results do not depend
on `num_threads` or `chunk_size`. Its fact refusals throw `ComparisonError`.

Other refusals throw `std::invalid_argument`. A non-finite distance throws
`std::runtime_error`. The storage overload scans every distance first under
the neighbor order, so a NaN cannot pass as "not a neighbor".

### k-nearest-neighbor graph and Jarvis-Patrick

`include/oecluster/clustering/KNNGraph.h` declares `KNNGraphOptions{k,
num_threads, chunk_size}`, the `KNNGraph` type, and `knn_graph` overloads for
`const StorageBackend&` and `PairwiseComparison&`.
`include/oecluster/clustering/JarvisPatrick.h` declares
`JarvisPatrickOptions{k, kmin, num_threads, chunk_size}`,
`JarvisPatrickResult` (`K()`, `KMin()`, `Method() == "jarvis_patrick"`), and
`jarvis_patrick` overloads for a `KNNGraph` plus `kmin`, for storage, and for
a comparison.

```cpp
OECluster::KNNGraphOptions options;
options.k = 6;
const OECluster::KNNGraph graph = OECluster::knn_graph(storage, options);
const auto result = OECluster::jarvis_patrick(graph, 3);
```

- `Indices()` and `Distances()` are row-major with `NumItems() * K()`
  entries. Row i never names i and is ordered by ascending (distance,
  index). The values are raw distances, not affinities.
- Every row is an independent bounded selection. Matrix paths distribute
  rows with `ThreadPool::ParallelFor`, and the comparison path uses one
  clone per running unit. A unit covers `max(1, chunk_size / (n - 1))` rows,
  and the result is identical for every `num_threads` and `chunk_size`.
- The comparison overload calls `Compare(min(i, j), max(i, j))` N(N-1)
  times. `Compare` must be repeatable across calls and clones.
- Sparse storage must hold every pair at or within its cutoff; the builder
  cannot verify that, and a missing nearer pair silently changes a row. An
  item with fewer than `k` distinct stored neighbors is refused with
  `std::invalid_argument`. Duplicate entries count once, with the value
  `Get()` reports.
- Validation order on every storage and comparison entry point: `chunk_size`,
  then zero items (empty result whatever `k` and `kmin` are), then
  `1 <= k <= n - 1`, then `kmin < k`, then the input checks. The
  `jarvis_patrick(const KNNGraph&, kmin)` overload checks zero items then
  `kmin < k`. A NaN or infinite distance read raises `std::runtime_error`
  (sparse: stored entries only). Comparison-facts refusals raise
  `ComparisonError`.
- The `KNNGraph` constructor validates its arrays: size, `k` range, index
  range, no self, no repeats, finite distances and row order. It throws
  `std::invalid_argument` otherwise, including when `num_items * k`
  overflows.
- Jarvis-Patrick links i and j when each is in the other's row and the rows
  share at least `kmin` items. Clusters are ordered by smallest member. A
  self-inclusive formulation's `k` is this `k + 1`, and its `kmin` is this
  `kmin + 2`.

### Leiden community detection

`include/oecluster/clustering/Leiden.h` declares `LeidenObjective`
(`Modularity`, `CPM`), `LeidenOptions{k, objective, resolution, prune,
theta, n_iterations, seed, num_threads, chunk_size}`, `LeidenResult`
(`Quality()`, `Iterations()`, `Objective()`, `Resolution()`, `K()`,
`Method() == "leiden"`), and `leiden` overloads for a `KNNGraph`, for
storage, and for a comparison.

```cpp
OECluster::LeidenOptions options;
options.k = 15;
options.objective = OECluster::LeidenObjective::CPM;
options.resolution = 0.05;
const auto result = OECluster::leiden(storage, options);
```

- The graph is reweighted internally by shared-nearest-neighbor Jaccard
  (`src/clustering/SNNWeights.h`): with N+(i) row i plus i, the edge {i, j}
  weighs s / (2(k + 1) - s) for s = |N+(i) ∩ N+(j)|. Only kNN arcs become
  edges, each pair once whichever rows name it, and a weight below `prune`
  is dropped. The weights are computed in parallel with
  `ThreadPool::ParallelFor` and do not depend on `num_threads`.
- Seurat's `k.param` is this `k + 1`. Seurat also links pairs that share
  neighbors without either naming the other; this graph does not.
- Modularity is Reichardt-Bornholdt, sum_c [e_c / m - resolution *
  (K_c / 2m)^2]; CPM is sum_c [e_c - resolution * N_c (N_c - 1) / 2]. 1.0
  is the modularity default; CPM over Jaccard weights typically needs a
  resolution well below the median weight. `Quality()` is the objective of
  the returned partition on the weighted graph.
- The optimizer (`src/clustering/LeidenEngine.h`) is serial: local moving
  with a queue, refinement that merges singletons only into well-connected
  subsets, and aggregation, repeated until a level does not shrink. Every
  returned cluster is connected. `n_iterations == -1` repeats passes until
  one leaves the canonical labels unchanged and counts that pass;
  otherwise exactly `n_iterations` passes run.
- Randomness comes from `std::mt19937_64` seeded with `seed`, with bounded
  integers, uniform doubles and shuffles written out rather than taken from
  the standard distributions, so the stream is the same on every platform.
  The same input, options and build give identical results; across
  platforms, a one-ulp `std::exp` difference or a different floating-point
  contraction can change a draw or a tie.
- Validation order on the storage and comparison overloads: `chunk_size`,
  then the options (`objective` a known enumerator, `resolution` finite and
  non-negative, `prune` finite and in [0, 1), `theta` finite and positive,
  `n_iterations >= -1`), then zero items (empty result), then at most
  `INT_MAX` items, then `1 <= k <= n - 1`, then everything `knn_graph`
  checks. The graph overload checks the options, zero items, then the item
  count. These refusals are `std::invalid_argument`, raised before any
  comparison runs. The raw overloads also raise what `knn_graph` raises:
  `ComparisonError` when the comparison's facts rule out ranking, and
  `std::runtime_error` when a distance is NaN or infinite.
- Memory at 1e6 items with k = 50: the kNN graph is 800 MB, the weighted
  CSR at most 1.2 GB, and building it needs another 0.6 GB that is freed
  before optimization. The raw overloads free the kNN graph before
  optimizing, which then holds the CSR, O(N) work arrays and the next
  level's CSR while aggregating.

## Partition Agreement

`PartitionAgreement.h` scores two labelings of the same samples against each
other. It reads labels and nothing else -- no distance matrix, unlike
`ClusterReport.h` -- so comparing two clustering methods does not require
owning the matrix either of them was built from. There are four entry points,
two per function, differing only in whether the labels arrive inside a
`ClusteringResult`:

```cpp
PartitionAgreement partition_agreement(
    const ClusteringResult& a, const ClusteringResult& b,
    const PartitionAgreementOptions& options = PartitionAgreementOptions());

PartitionAgreement partition_agreement(
    const std::vector<ClusterLabel>& a, const std::vector<ClusterLabel>& b,
    const PartitionAgreementOptions& options = PartitionAgreementOptions());

PartitionAgreement scaffold_agreement(
    const ClusteringResult& result,
    const std::vector<std::string>& scaffold_labels,
    const PartitionAgreementOptions& options = PartitionAgreementOptions());

PartitionAgreement scaffold_agreement(
    const std::vector<ClusterLabel>& labels,
    const std::vector<std::string>& scaffold_labels,
    const PartitionAgreementOptions& options = PartitionAgreementOptions());
```

The `ClusteringResult` overloads read `Labels()` only. Unlike `cluster_report`,
they never cross-check `Labels()` against `Members()`, so a result whose two
views disagree is not rejected here. All four raise `std::invalid_argument`
when the two sides differ in length or either is empty.

`PartitionAgreement` carries seven `double` metrics. Two are pair-counting --
`adjusted_rand_index` (Hubert-Arabie; 1.0 is exact agreement, 0.0 is the value
expected by chance, negative is worse than chance) and `fowlkes_mallows` (the
geometric mean of pair precision and pair recall). Four are entropy-based:
`normalized_mutual_information`, `homogeneity`, `completeness`, and
`v_measure`. The seventh, `adjusted_mutual_information`, is opt-in.

`v_measure` and `normalized_mutual_information` are assigned from one computed
value, `2*MI/(H(a)+H(b))`, and are bitwise equal on every input. That is the
definition rather than the harmonic mean of the two components, because the
harmonic form is 0/0 both when MI is zero with two positive entropies and when
one entropy is zero, and the composite is well defined in each case.

### The positional rule

Side A is always the first argument. `homogeneity` is `MI / H(a)` and
`completeness` is `MI / H(b)`, so swapping the arguments exchanges that pair --
and with it the counts `num_clusters_a` and `num_clusters_b`. Every other
metric is symmetric, `adjusted_mutual_information` only to within rounding: a
swap transposes the contingency table and exchanges the two marginal values
within each expected-MI term, so those terms evaluate in a different
floating-point order and the last bits can move.

`scaffold_agreement` puts the clustering on side A and the scaffold annotation
on side B, so `completeness` carries the scaffold-purity reading -- whether
each cluster's members share a single scaffold, the same question
`RepresentativeMetrics::scaffold_purity` asks per cluster -- and `homogeneity`
carries its transpose, whether each scaffold landed in a single cluster.

### Noise handling

A sample is noise on a side when its label there is negative -- not only
`NOISE_LABEL`, matching how `ClusterReport` and `labels_to_clusters` already
read labels. In `scaffold_agreement` an empty scaffold string is noise on the
string side. `PartitionAgreementOptions::noise_handling` chooses one of three
readings:

- `NoiseHandling::Singletons` (the default) makes each noise sample its own
  cluster, so `num_clusters_a` and `num_clusters_b` include one entry per noise
  sample.
- `NoiseHandling::Grouped` collapses each side's noise into one cluster, which
  is how scikit-learn reads a -1 label.
- `NoiseHandling::Excluded` drops a sample noisy on either side from both
  partitions, so `num_samples` can be smaller than the input length. Every
  other mode leaves it equal to the input length.

### Adjusted mutual information

`PartitionAgreementOptions::compute_adjusted_mutual_information` is `false` by
default. The other six metrics come essentially free once the contingency table
is built, and AMI is the only one that costs more than that pass: its
expected-MI correction first builds an O(N) table of log factorials -- `8*(n+1)`
bytes, about 800 KB at `n = 100,000` -- and then sums over pairs of distinct
marginal values, that is over distinct cluster sizes rather than over clusters,
each term walking the hypergeometric support. That is why it is opt-in.
`PartitionAgreement::requested` records the request and not the outcome, the
same convention `ClusterReportRequested` uses: `false` means nobody asked, and
`true` with NaN means asked and undefined.

AMI's accuracy is limited where its denominator -- the mean entropy minus the
expected mutual information -- approaches zero, which happens when both
partitions are close to all-singleton. Numerator and denominator are then each
a difference of nearly equal sums over N terms, and the quotient loses
significance: two 1.5-million-sample partitions differing by one merged pair
have a true AMI of zero and report about 0.035. Partitions whose denominator is
order one are unaffected. The limit is inherent to computing the correction in
double precision -- scikit-learn shares it -- rather than a property of this
implementation.

### Divergences from scikit-learn

An undefined metric is NaN here rather than a convention. scikit-learn
substitutes a value in five cases, and this table is the complete list:

| Case | Fixture | scikit-learn 1.9.1 | OECluster |
| --- | --- | --- | --- |
| Fewer than two surviving samples | `a = b = {0}` | every metric `1.0`, except `fowlkes_mallows = 0.0` | every metric NaN |
| Side A is one cluster, partitions differ | `{0,0,0,0}` vs `{0,0,1,1}` | `homogeneity = 1.0` | NaN |
| Side B is one cluster, partitions differ | `{0,0,1,1}` vs `{0,0,0,0}` | `completeness = 1.0` | NaN |
| One side all singletons, partitions differ | `{0,1,2,3}` vs `{0,0,1,1}` | `fowlkes_mallows = 0.0` | NaN |
| Both sides all singletons, N >= 2 | `{0,1,2,3}` vs `{3,2,1,0}` | `fowlkes_mallows = 0.0` | `1.0` |

The last row is not a NaN case: two all-singleton partitions are the same
partition, so the six always-computed metrics are 1.0, `fowlkes_mallows` among
them, and `adjusted_mutual_information` is 1.0 too whenever it was requested.

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
