# Python API

The top-level `oecluster` package is the recommended way to use OECluster from
notebooks, scripts, and data-processing pipelines. It accepts OpenEye molecule
objects, computes pairwise and cross distance matrices, clusters them with
several algorithms, selects representatives, and reports cluster quality.

Advanced users can still reach the lower-level SWIG bindings at
`oecluster.oecluster` (and the compiled extension `oecluster._oecluster`) when
they need direct control over the C++ classes.

The generated reference for every public symbol is in
[](api/python). This page is the narrative guide.

## Distance Matrices

`pdist()` computes the symmetric within-set pairwise distance matrix. The first
argument is the item set (OpenEye molecules), the second names the comparison
method, and method-specific keywords are forwarded to that comparison.

```python
import oecluster

dm = oecluster.pdist(
    mols,
    "fingerprint",
    fp_type="morgan",
    metric="tanimoto",
    num_threads=8,
)

print(dm.condensed)     # scipy-compatible condensed vector
print(dm.squareform())  # full NxN matrix
dm.to_file("distances.npz")
```

`pdist()` returns a `SymmetricDistanceMatrix`. The supported comparison names
are `"fingerprint"`, `"rocs"`, `"superpose"`, and `"sitehopper"` (the last is a
superpose mode). See [Comparison Methods](#comparison-methods) for the keyword
arguments each one accepts.

Storage is chosen by keyword:

```python
# Sparse: keep only short distances, for threshold graph algorithms.
dm = oecluster.pdist(mols, "fingerprint", metric="tanimoto", cutoff=0.35)

# Memory-mapped: dense distances that fit on disk but not in RAM.
dm = oecluster.pdist(mols, "fingerprint", output="distances.mmap")
```

`cdist()` computes the rectangular NxM cross distance matrix between two
distinct item sets and returns a `CrossDistanceMatrix`:

```python
cross_dm = oecluster.cdist(queries, targets, "fingerprint", metric="tanimoto")
print(cross_dm.matrix[0, :])  # distances from first query to all targets
```

Use `cdist()` when you need distances from one set (for example virtual
screening hits) to another (for example an in-house collection), rather than
the symmetric within-set distances `pdist()` computes.

Both functions accept `similarity=True` to return similarities instead of
distances where the comparison supports it, plus shared parallelism controls
(`num_threads`, `chunk_size`) and a `progress` callback.

Reload a saved matrix with `load_distance_matrix()`:

```python
dm = oecluster.load_distance_matrix("distances.npz")
```

## Clustering

All clustering functions take a `DistanceMatrix` (except BitBirch, which takes
an OEFP fingerprint batch) and return a result that subclasses
`ClusteringResult`.

```python
result = oecluster.butina(dm, threshold=0.35, reordering=False)

# Length-n, scikit-learn-style assignment array.
labels = result.labels

# Tuple of clusters, each a tuple of member indices.
for cluster_id, cluster in enumerate(result.clusters):
    print(cluster_id, len(cluster))

len(result)        # number of clusters
result[0]          # first cluster's member indices
result.method      # "butina"
```

The available algorithms and their key parameters:

| Function | Input | Key parameters |
|----------|-------|----------------|
| `butina(dm, threshold, ...)` | `DistanceMatrix` | `threshold`, `reordering` |
| `dbscan(dm, eps, ...)` | `DistanceMatrix` | `eps`, `min_samples` |
| `hdbscan(dm, ...)` | `DistanceMatrix` | `min_cluster_size`, `min_samples`, `cluster_selection_method` |
| `agglomerative(dm, ...)` | `DistanceMatrix` | `n_clusters`, `distance_threshold`, `linkage` |
| `bitbirch(fingerprints, ...)` | `oefp.OEFPBatch` | `threshold`, `branching_factor`, `merge_criterion` |

Algorithm-specific outputs live on the specific result subclass — for example
`DBSCANResult.core_sample_indices` or `BitBirchResult.centroids` — so a result
never carries fields that do not apply to its algorithm.

### BitBirch Variants

`bitbirch()` incrementally inserts binary fingerprints into a Birch-style tree
and merges leaf subclusters. Two specialized strategies build on it:

`bitbirch_recluster()` runs a two-stage pass. The first stage fits the
fingerprints at `initial_threshold`, then re-clusters the leaf summaries at
`second_threshold` with an optional `second_tolerance` penalty:

```python
result = oecluster.bitbirch_recluster(
    fingerprints,
    initial_threshold=0.65,
    second_threshold=0.7,
    branching_factor=50,
    mode="strict_parity",
)
```

`bitbirch_refine()` fits a tree and then applies refinement passes. Enable
`redistribute_largest_cluster` to redistribute the largest cluster, or set
`reassign_top_clusters` (a count of two or more) to reassign members of the
top-K largest clusters against cluster centroids:

```python
result = oecluster.bitbirch_refine(
    fingerprints,
    threshold=0.65,
    branching_factor=50,
    redistribute_largest_cluster=True,
    reassign_top_clusters=3,
)
```

All three return a `BitBirchResult` with `labels`, `clusters`, `centroids`, and
`cluster_sizes`. The `mode` parameter accepts `"strict_parity"` (exact
reference-implementation parity) or `"fast"` (partition-merge parallelism with
deterministic output). Fast mode applies to `bitbirch()` and
`bitbirch_recluster()`; `bitbirch_refine()` always runs in strict parity
because its prune/reassign passes are order-sensitive, so its `mode` argument is
accepted for API symmetry but does not change behavior.

## Representatives

Representatives are cluster members chosen to summarize a cluster. They are
distinct from BitBirch centroid fingerprints, which are synthetic summaries and
may not correspond to a real molecule.

`representative()` returns the single best member index, `rank_representatives()`
returns a ranked list of `ClusterRepresentative` objects, and
`select_representatives()` returns more than one.

```python
medoid = oecluster.representative(cluster, dm, method="medoid")
minimax = oecluster.representative(cluster, dm, method="minimax")

ranked = oecluster.rank_representatives(cluster, dm, method="medoid")
best = ranked[0]
print(best.member, best.metrics.cluster_radius)

selected = oecluster.select_representatives(
    cluster, dm, k=3, method="medoid", selection="diversity",
)
```

The supported `method` values are `"medoid"`, `"minimax"`,
`"highest_neighborhood"`, and `"weighted_medoid"`. The
`"highest_neighborhood"` method requires a `threshold`. The weighted medoid
biases centrality by project metadata:

```python
ranked = oecluster.rank_representatives(
    cluster,
    dm,
    method="weighted_medoid",
    alpha=1.0,
    beta=0.5,
    gamma=0.4,
    liability_penalties=liabilities,
    priority_scores=priorities,
    scaffold_labels=scaffolds,
)
```

`select_representatives()` takes `selection="score"` (top-K by representative
score) or `selection="diversity"` (spread across the cluster). Each
`ClusterRepresentative` carries a `RepresentativeMetrics` object exposing
`mean_distance_to_cluster`, `max_distance_to_cluster`,
`neighbor_fraction_at_threshold`, `cluster_radius`, `cluster_diameter`,
`silhouette_like_score`, `scaffold_purity`, `representative_rank`, and related
quality fields.

## Cluster Quality Reports

`cluster_report()` computes a method-agnostic scorecard for any clustering
result, and `compare_reports()` aligns two or more scorecards side by side:

```python
butina_result = oecluster.butina(dm, threshold=0.35)
dbscan_result = oecluster.dbscan(dm, eps=0.35, min_samples=5)

butina_report = oecluster.cluster_report(butina_result, dm)
dbscan_report = oecluster.cluster_report(dbscan_result, dm)

print(butina_report)
print(oecluster.compare_reports(butina_report, dbscan_report))
```

Distances are Tanimoto/Jaccard (`distance = 1 - similarity`), so smaller means
more similar. Choose a threshold preset to match the use case, or override the
individual thresholds:

```python
report = oecluster.cluster_report(butina_result, dm, preset="tight")

report = oecluster.cluster_report(
    butina_result,
    dm,
    coverage_thresholds=[0.2, 0.3],
    boundary_threshold=0.25,
)
```

The presets are `"tight"`, `"default"`, and `"diversity"`. Metrics that are
undefined for a given clustering (for example separation metrics when there is
only one cluster) are reported as `NaN`. `num_noise` (HDBSCAN-style unclustered
points, label `-1`) and `num_singletons` (size-1 clusters) are always reported
separately; `treat_noise_as_singletons=True` (the default) folds noise into the
singleton interpretation. The report requires complete pairwise distances
(dense or memory-mapped storage); a sparse (`cutoff`) matrix raises.

## Comparison Methods

The comparison keyword arguments are forwarded by `pdist()` and `cdist()` to the
selected method.

### Fingerprint

| Parameter | Values | Default |
|-----------|--------|---------|
| `fp_type` | `morgan`, `atom_pair` | `morgan` |
| `metric` | `tanimoto`, `dice`, `manhattan` | `tanimoto` |
| `numbits` | Fingerprint size | `2048` |
| `min_distance` | Minimum Atom Pair graph distance | `0` |
| `max_distance` | Morgan radius or maximum Atom Pair graph distance | `2` |

Distance mode maps Tanimoto to Jaccard distance, uses OEFP's Dice distance for
Dice, and returns raw Manhattan distance. Similarity mode is supported for
Tanimoto.

### ROCS

| Parameter | Values | Default |
|-----------|--------|---------|
| `score_type` | `combo_norm`, `combo`, `shape`, `color` | `combo_norm` |
| `color_ff_type` | `1` (ImplicitMillsDean), `2` (ExplicitMillsDean) | `1` |

Distance ranges: combo_norm [0, 1], combo [0, 2], shape [0, 1], color [0, 1].

### Superpose

| Parameter | Values | Default |
|-----------|--------|---------|
| `method` | `global_carbon_alpha`, `global`, `ddm`, `weighted`, `sse`, `sitehopper` | `global_carbon_alpha` |
| `score_type` | `auto`, `rmsd`, `tanimoto`, `patch_score` | `auto` |
| `predicate` | oeselect atom selection expression | empty |
| `ref_predicate` | Reference atom selection override | empty |
| `fit_predicate` | Fit atom selection override | empty |

Accepts molecules directly or design units; design units are converted to their
protein components automatically. Use `method="sitehopper"` for binding-site
patch-score comparison.

## Storage Backends

`pdist()` selects the backend from its keywords: dense by default, sparse when
`cutoff` is set, and memory-mapped when `output` is a path. The backends are
also available directly as `DenseStorage`, `MMapStorage`, and `SparseStorage`
for advanced use through the lower-level API.

| Backend | Use case | Memory |
|---------|----------|--------|
| `DenseStorage` | Default full pairwise matrix | O(N^2) in memory |
| `MMapStorage` | Full pairwise matrix backed by a file | O(N^2) on disk |
| `SparseStorage` | Cutoff-filtered results | Stored entries only |

All backends use scipy-compatible condensed distance-matrix indexing.

## Advanced C++ Binding Access

Most users do not need this section. The generated SWIG wrapper is available as
`oecluster.oecluster` and the compiled extension as `oecluster._oecluster` for
users who need direct access to the C++ options and classes. The
`FingerprintComparison`, `ROCSComparison`, and `SuperposeComparison` Python
wrappers, and the option structs (`PDistOptions`, `ButinaOptions`, and the
rest), are re-exported on the top-level package.

## Exceptions

Invalid arguments raise standard Python exceptions: unknown comparison or
representative method names, a `"highest_neighborhood"` request without a
`threshold`, a sparse matrix passed where complete distances are required, and
similar misuse raise `ValueError` or `RuntimeError` with a descriptive message.
