# Python API

The top-level `oecluster` package is the recommended way to use OECluster from
notebooks, scripts, and data-processing pipelines. It accepts OpenEye molecule
objects, computes pairwise and cross distance matrices, clusters them with
several algorithms, selects representatives, and reports cluster quality.

Advanced users can still reach the lower-level SWIG bindings at
`oecluster.oecluster` (and the compiled extension `oecluster._oecluster`) when
they need direct control over the C++ classes.

A generated reference for the classes and functions defined in the `oecluster`
package is in [](api/python). The option and storage structs that come straight
from the SWIG bindings -- `FingerprintOptions`, `ButinaOptions` and
`DenseStorage` among them -- are outside it. This page is the narrative guide.

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
are `"descriptor"`, `"fingerprint"`, `"rmsd"`, `"rocs"`, `"sitehopper"`, and
`"superpose"` (`"sitehopper"` is a superpose mode). See
[Comparison Methods](#comparison-methods) for the keyword arguments each one
accepts, and [Metric Requirements](#metric-requirements) for which of them
produce a matrix the clustering algorithms will accept.

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

### Internal validity indices

Seven internal cluster-validity indices sit alongside the original scorecard
fields. Five of them are always computed:

| Metric | Direction | What it reports |
|--------|-----------|-----------------|
| `calinski_harabasz_medoid` | higher is better | Between-cluster scatter divided by within-cluster scatter, each scaled by its degrees of freedom. `NaN` when there are fewer than two clusters, when every cluster is a singleton, or when the within-cluster scatter is zero. |
| `davies_bouldin_medoid` | lower is better | Mean over clusters of the worst ratio of two clusters' spreads to the distance between their centers. `NaN` when there are fewer than two clusters, and infinite when two centers coincide, which is reported rather than divided by. |
| `dunn_mean_separation_mean_diameter` | higher is better | Smallest mean between-cluster distance over the largest mean within-cluster distance. Averaging both terms makes it far less outlier-sensitive than `dunn_index`, which takes the extreme of each. |
| `dunn_medoid_separation_medoid_spread` | higher is better | Smallest medoid-to-medoid distance over the largest medoid spread, the spread being twice a cluster's mean medoid-to-member distance, averaged over all `n_k` members and so counting the medoid's own zero. This field ignores `representative_method` and always uses the true medoid, so a report requested with `representative_method="minimax"` still reports medoid-based values here. |
| `point_biserial` | higher is better | Correlation between the pairwise distances and the within-versus-between split. Positive means between-cluster pairs are the more distant ones. The sign convention is stated because published sources differ on it. |

> Two of these five are **medoid-substituted** indices. The published
> Calinski-Harabasz and Davies-Bouldin definitions use centroids, which do not
> exist for a distance matrix; each cluster's medoid stands in for its centroid,
> and the global medoid stands in for the grand mean. The values are therefore
> not comparable with published figures or with scikit-learn's. Both ignore
> `representative_method` and always use the true medoid, so a report requested
> with `representative_method="minimax"` still reports medoid-based values here.

The remaining two of the seven are read off the ranked pairwise distances and
are computed only when `compute_pair_rank_indices=True`:

| Metric | Direction | What it reports |
|--------|-----------|-----------------|
| `c_index` | lower is better | Where the within-cluster distance sum falls between the smallest and largest sums the same number of pairs could have taken. Zero is perfect. `NaN` when there are no within-pairs, no between-pairs, or no spread between those two bounds. |
| `baker_hubert_gamma` | higher is better | Rank correlation over couples of one within-pair and one between-pair: the concordant count minus the discordant, over their total. Ranges from -1 to 1. `NaN` when no couple is either. |

### Optional report stages

Two options, both off by default, turn on the report's optional stages:

```python
report = oecluster.cluster_report(
    butina_result,
    dm,
    compute_pair_rank_indices=True,
    compute_per_cluster_records=True,
)
```

`compute_pair_rank_indices` is off by default because it is the one optional
stage whose pair-scaled cost grows with all `Nc` clustered points rather than
with the largest single cluster: it materialises every pairwise distance among
them as sortable arrays, `Nc * (Nc - 1) / 2` doubles in total, which is roughly
400 MB at `Nc = 10,000` and 10 GB at `Nc = 50,000`. A failed allocation raises
`MemoryError`. One option covers both indices rather than two because they come
off the same sorted arrays; once those are paid for, the second index is nearly
free.

`compute_per_cluster_records` populates `report.records`, one `ClusterRecord`
per cluster in member-list order. Each record carries the cluster's `label`,
`size`, `representative`, `mean_intra_distance`, `median_intra_distance`,
`radius`, `diameter`, `mean_representative_distance`, `nearest_cluster`,
`nearest_cluster_distance`, `silhouette` and `boundary_violations`. A record's
`boundary_violations` counts pairs *involving* that cluster, so summing the
column gives twice the scorecard's `boundary_violations`, which counts each
pair once. `ClusterRecord` is a `typing.NamedTuple`, so
`pandas.DataFrame(report.records)` works without a conversion step; pandas is
not a dependency. This stage is pair-scaled too, but bounded by the largest
cluster rather than by the whole clustering: it buffers that cluster's
`n * (n - 1) / 2` distances to take their median, and the median is taken over a
copy of the buffer, so roughly 400 MB for the buffer and 400 MB again for the
copy, transiently, at `n = 10,000`.

`report.noise_coverage_at` restricts the coverage curve to the noise points,
parallel to `coverage_thresholds` in the same way `coverage_at` is. Its length
always matches `coverage_at`: the threshold count when the clustering has at
least one cluster, and empty when it has none. Every entry is `NaN` when the
clustering has no noise, rather than 0.0, which would read as "no noise point is
covered" instead of "there is nothing to cover".

### Asked-for versus undefined

`report.requested` is a `ClusterReportRequested` naming the two optional
computations the caller asked for, `pair_rank_indices` and
`per_cluster_records`. It records the request and not the outcome, which is what
makes a `NaN` readable: `False` means nobody asked, and `True` with `NaN` means
asked and undefined.

`compare_reports(...).to_table()` applies the same distinction to its cells. A
cell is `None` when that report never asked the question -- an opt-in metric it
did not request, a threshold it did not use, or a threshold it does carry but
answered nothing at, as when the clustering has no clusters and both coverage
curves come back empty -- while `nan` keeps its single meaning of asked and
undefined. `__repr__` renders `None` as `--`. A caller reading cells as floats
must test for `None` before doing arithmetic on them.

## Metric Requirements

`butina()`, `dbscan()`, `hdbscan()`, `agglomerative()`, and `cluster_report()`
assume their input is a metric: an item's distance to itself is zero, and the
triangle inequality holds. Those five are the entry points that check.
`representative()`, `rank_representatives()`, and `select_representatives()`
also take a distance matrix, and consult none of these facts. Not every
comparison produces a metric, so each matrix records what its comparison
actually guarantees:

```python
dm = oecluster.pdist(mols, "fingerprint", metric="dice")
dm.is_distance           # True
dm.metric_capabilities   # {'zero_self': True, 'triangle': False}
dm.data_integrity        # 'complete'
dm.metric_probe          # 'not_run'
```

`is_distance` is kept separate from the capabilities because it answers a
different question: not whether the numbers form a metric, but which direction
they run. A similarity is not a distance however its diagonal behaves.

Three refusals cannot be overridden, because the algorithms would return
plausible-looking wrong answers rather than fail: a matrix of similarities, a
matrix whose measure does not score an item as identical to itself, and a
matrix stamped `data_integrity == 'nan_present'` or holding a non-finite
distance. The remedies are, respectively, to recompute with
`similarity=False`, or to choose a metric that has a distance form at all --
`tversky` does not, and refuses `similarity=False`; to choose a comparison or
configuration whose measured diagonal vanishes; and to recompute with
`missing='complete_case'` if the values came from descriptors, or otherwise to
remove the offending items. The third check does not trust the stamp alone --
it also scans the stored distances, so a matrix edited through `.condensed`
after it was stamped is still caught.

`zero_self` is measured rather than assumed. A ROCS comparison scores every
molecule against itself at construction and stamps the capability from what it
finds, so the same `score_type` can pass over one molecule set and refuse over
another. As a distance, two independent things make a diagonal fail to vanish.
A shape self-overlay does not always saturate: on small compact molecules the
best self-overlay comes back a little short, leaving a `score_type="shape"`
self-distance around `1e-2` rather than `0`. Methane, water, Cl2 and Br2 have
all been measured there, while ethane, benzene and ethanol seat on themselves
exactly -- and whether a given small molecule lands on zero depends on the
conformer it was embedded with, which is why the stamp is measured per set and
not tabulated per `score_type`. Separately, a colour term contributes nothing
for a molecule that carries no colour features, so `color` self-scores `1.0`
as a distance there; methane and Cl2 do that, water does not. `combo_norm`
(the default) and `combo` average the two terms and inherit both causes:
methane's `combo_norm` self-distance was measured at `5.1e-01`, and water's at
`5.6e-03` even though water's colour term vanishes. No `score_type` is
guaranteed to vanish for every molecule, so read
`dm.metric_capabilities['zero_self']` for the set in hand rather than choosing
a `score_type` in the hope of one. None of that carries over to
`similarity=True`, where a self-score sits at or just under 1.0 instead of
vanishing -- the same shortfall as above, read from the other end: methane's
shape self-similarity measures `0.9896` and Cl2's `0.9863` where ethane's is
exactly `1.0`. Either way the value is nowhere near zero, so on a non-empty set
`shape` stamps `zero_self` false, and so do `combo_norm`, `combo` and `color`
on coloured molecules. The one exception among non-empty sets is a colour
similarity over molecules with no colour features, which self-scores 0.0 and is
stamped `zero_self` true. It is refused anyway: a similarity fails
`is_distance` before `zero_self` is consulted, which is why the two facts are
kept apart. An empty set is stamped true as well, and vacuously so: `pdist`
accepts it, and with no diagonal to measure the stamp comes back true for every
`score_type` in either mode.

The remaining checks are soundness warnings that `allow_nonmetric=True`
overrides: a measure known to violate the triangle inequality, distances
scored on a per-pair feature subset (`missing='ignore'`), and violations found
by the ingress probe.

```python
oecluster.butina(dm, 0.4, allow_nonmetric=True)
```

A matrix built with `SymmetricDistanceMatrix.from_condensed()` is checked at
ingress: values must be finite and non-negative, a square input must
additionally be symmetric to within `rtol=1e-5, atol=1e-8` and have a diagonal
zero to within `3e-2` of the largest distance, and a sample of triples is
tested against the triangle inequality. The probe can only disprove -- finding
no violation leaves `triangle` unknown, and the gate stays out of the way.

```python
dm = oecluster.SymmetricDistanceMatrix.from_condensed(scipy_condensed)
dm.metric_probe      # 'no_violations_found'
dm.is_distance       # 'unknown'
```

The three capability facts -- `is_distance`, `zero_self` and `triangle` --
stay `'unknown'` on that path. There is no comparison to interrogate, and a
caller's assurance is not evidence, so the probe result is the only capability
claim made. What the checks do record is `data_integrity`, stamped
`'complete'`, and the probe's own `metric_probe`, `probe_violations` and
`probe_sampled`. Pass `probe_triples=0` to skip the probe alone, or
`check=False` to skip the value checks and the probe together, which also
leaves `data_integrity` at `'unknown'` since nothing measured it. The refusals
that do not read the numbers still run either way, among them the shape and
length arithmetic, the label count, the masked-array refusal and the
complex-dtype refusal.

Distance matrices written before 5.0.0 record no facts at all and load with
`is_distance`, `zero_self`, `triangle` and `data_integrity` every one
`'unknown'`. An unknown fact never refuses on its own, so those files cluster
as they always did -- except for the non-finite scan, which reads the numbers
rather than the stamp and so applies to them too.

## Comparison Methods

The comparison keyword arguments are forwarded by `pdist()` and `cdist()` to the
selected method.

### Fingerprint

| Parameter | Values | Default |
|-----------|--------|---------|
| `fp_type` | `morgan`, `atom_pair`, `topological_atom_pair`, `topological_torsions` | `morgan` |
| `storage` | `binary`, `count`, `sparse`, `sparse_count` | `binary` |
| `metric` | `jaccard`, `tanimoto`, `dice`, `sokal_sneath`, `matching`, `rogers_tanimoto`, `russell_rao`, `kulsinski`, `sokal_michener`, `euclidean`, `manhattan`, `chebyshev`, `hamming`, `canberra`, `bray_curtis`, `minkowski`, `tversky` | `tanimoto` |
| `numbits` | Fingerprint size in bits | `2048` |
| `radius` | Morgan radius | `2` |
| `min_distance` | Minimum atom-pair graph distance | `1` |
| `max_distance` | Maximum atom-pair graph distance | `30` |
| `torsion_atom_count` | Torsion path length | `4` |
| `use_chirality` | Distinguish stereocenters | `False` |
| `p` | Minkowski order | `2.0` |
| `tversky_alpha` | Tversky reference weight | `0.5` |
| `tversky_beta` | Tversky fit weight | `0.5` |

Every metric in that list works under the default `similarity=False` except
`tversky`, which has no distance form and requires `similarity=True`.

Family, storage and metric names are matched case-insensitively.
`topological_atom_pair` is an alias of `atom_pair`: OEFP has one atom-pair
generator and it is the 2D graph-distance model. Every family and storage
combination works except `topological_torsions` with `sparse_count`, for which
OEFP provides no batch kernel.

Binary and sparse storage record which features are present; count and
sparse-count storage record how many times each occurs. The seven numeric
metrics -- `euclidean`, `manhattan`, `chebyshev`, `hamming`, `canberra`,
`bray_curtis` and `minkowski` -- read those counts. The other ten are bit-set
coefficients that would binarize their input first, so a counted storage
refuses them and names the numeric alternatives. Counts matter when molecules
differ mainly in how often a feature repeats: binary Morgan fingerprints of
hexadecane and triacontane are identical, putting their Tanimoto distance at
`0.0`, while count fingerprints under `bray_curtis` put them `0.313` apart.

Similarity mode (`similarity=True`) is available for `tanimoto` and `tversky`,
the two metrics with a similarity form. That is the only place `jaccard` and
`tanimoto` differ: their distances are equal to the bit, but only `tanimoto`
can be asked for a similarity. `tversky` goes the other way -- it has no
distance form, so it requires `similarity=True` and refuses the default
`similarity=False` rather than quietly returning a similarity anyway. A
similarity is not a distance, and the clustering entry points refuse one
outright; see [Metric Requirements](#metric-requirements).

Naming an option that the rest of the configuration would ignore raises
`TypeError` naming the ignored option and a way out -- a replacement option
where one exists, and otherwise dropping the option or changing the setting
that made it inapplicable -- rather than accepting a value that has no effect:

```python
oecluster.pdist(mols, "fingerprint", max_distance=2)
# TypeError: max_distance does not apply to fp_type='morgan'; it belongs to
# 'atom_pair'. Use radius instead, or select one of those families.

oecluster.pdist(mols, "fingerprint", p=3.0)
# TypeError: p does not apply to metric='tanimoto'; it belongs to
# 'minkowski'. Drop p, or select metric='minkowski'.
```

The rule covers `radius`, `min_distance`, `max_distance` and
`torsion_atom_count` against the selected `fp_type`; `numbits` against a
sparse `storage`; and `p`, `tversky_alpha` and `tversky_beta` against the
selected `metric`. `use_chirality` reaches every family and is never rejected.
The C++ `FingerprintOptions` struct applies none of these rules: it receives
values rather than argument names and cannot tell a default from a choice.
C++ still speaks first either way -- a family, storage or metric it refuses
outright raises before any of these rules run, so the message names the
problem that has to be fixed first.

**Changed in 5.0.0.** `max_distance` no longer sets the Morgan radius. It is
the atom-pair window only, the Morgan radius is `radius`, and passing
`max_distance` with `fp_type="morgan"` now raises rather than quietly doing
something else. The atom-pair window itself now defaults to OEFP's `1`-`30`
instead of the previous `0`-`2`. That old window discarded every pair more
than two bonds apart, at a cost that grows with the molecule: it leaves
ethanol's three on-bits untouched, but takes aspirin from 68 on-bits to 27,
caffeine from 79 to 30, and hexadecane from 70 to 12. Code that passed
`max_distance=2` for Morgan must pass `radius=2`; code that wanted the old
atom-pair window must ask for it with
`fp_type="atom_pair", min_distance=0, max_distance=2`.

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

### Descriptor

| Parameter | Values | Default |
|-----------|--------|---------|
| `sources` | `openeye`, `mordred`, `rdkit` | `["openeye"]` |
| `columns` | Explicit column names | every numeric column of the selected sources |
| `groups` | Descriptor group names | none |
| `metric` | `euclidean`, `manhattan`, `chebyshev`, `hamming`, `canberra`, `bray_curtis`, `minkowski`, `standardized_euclidean`, `seuclidean`, `mahalanobis` | `standardized_euclidean` |
| `variances` | Explicit per-column variances | fitted from the input |
| `inverse_covariance` | Explicit inverse covariance for `mahalanobis` | fitted from the input |
| `missing` | `complete_case`, `propagate`, `ignore` | `complete_case` |
| `p` | Minkowski order | `2.0` |

`seuclidean` is an alias of `standardized_euclidean`. Non-numeric columns are
dropped from any selection, on every metric. A zero-variance column is dropped
only where a variance is being fitted -- `standardized_euclidean`, `seuclidean`
and `mahalanobis` -- and only when that fit actually runs: over `["CCO", "CCC"]`
those three metrics lose `HeavyAtomCount`, `FractionCsp3`, `AromaticRingCount`
and `RotatableBondCount` that way, while `euclidean`, `hamming` and the other
raw metrics keep all eleven openeye columns. Supplying `variances=` or
`inverse_covariance=` skips the fit, so the selection is used exactly as given
and nothing is dropped for variance even on a fitted metric.

`descriptor_statistics()` has no metric to key on and fits regardless, so it
reports those same four as `zero-variance` over that pair and leaves them out
of its `columns`. Its `dropped` list and a raw-metric comparison's
`DroppedColumns()` therefore disagree over the same molecules, and the
difference reaches the distances rather than staying in the metadata:
`pdist(mols, "descriptor", metric="hamming")` scores `["CCO", "CCC"]` at
`0.636364` over all eleven columns, and the same call restricted to the
`columns` that `descriptor_statistics()` hands back scores `1.0` over the seven
the fit left.

The drops are readable as `(name, reason)` pairs from `descriptor_statistics()`
under its `dropped` key, and from a prebuilt `DescriptorComparison` through
`DroppedColumns()` and `DroppedReasons()`; `pdist()` does not carry them in
`params`.

Descriptor columns span wildly different scales, so the default metric
standardizes each by the variance fitted over the molecules passed in. That
fit is data-dependent: the same pair of molecules gets a different distance in
a different set. Pass `columns=` and `variances=` from
`descriptor_statistics()` to fix the scaling across runs. Both have to travel
together, because the values are matched to columns by position, and so does
the same `sources=` the statistics were fitted over.

```python
stats = oecluster.descriptor_statistics(mols)
dm = oecluster.pdist(mols, "descriptor",
                     columns=stats["columns"], variances=stats["variance"])
```

`missing` controls what happens to a molecule whose descriptor value is
absent. `complete_case` (the default) drops those molecules before computing
anything and records them in `params["excluded_items"]` as
`[index, "missing-descriptor"]` pairs, indexed against the list as passed. On
the `cdist()` path each side is filtered against its own indices and reported
separately, under `params["excluded_items_a"]` and
`params["excluded_items_b"]`; a side that lost nothing gets no key at all.
`propagate` keeps every molecule and lets the absence flow into the distances
as NaN; it stamps `data_integrity` as `'nan_present'` on the strength of the
policy, whether or not a NaN actually reached the matrix, and the clustering
entry points refuse that with no override available. `ignore` scores each pair
over the features both molecules have, which produces distances that are not
mutually comparable; it stamps `'subset_scored'`, which the clustering entry
points refuse unless `allow_nonmetric=True`. An observed NaN outranks that
stamp: a value that is present but not finite is not skipped, and takes
`data_integrity` to `'nan_present'` under `ignore` too. `ignore` is not
available with `standardized_euclidean` or `mahalanobis` at all, since their
fitted transforms mix columns and a per-pair subset of them is incoherent.

### RMSD

| Parameter | Values | Default |
|-----------|--------|---------|
| `overlay` | Superpose before measuring | `False` |
| `automorph` | Minimize over graph automorphisms | `True` |
| `heavy_only` | Ignore hydrogens | `True` |
| `expand_conformers` | Expand multi-conformer inputs into one item per pose | `True` |

Compares poses of one molecule. Every input must share a topology and carry
coordinates, and a mismatch raises rather than returning the `-1.0` that
`OERMSD` reports for it. Use `rocs` to compare different molecules by shape.
With `expand_conformers` left on, a multi-conformer `OEMol` becomes one item
per conformer, labeled `"<title>:conf<n>"` with `mol_<index>` standing in for
an empty title; the caller's molecules are never modified.
`expand_conformers` is a `pdist()`/`cdist()` keyword only -- the
`RMSDComparison` factory takes the molecules exactly as passed.

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
users who need direct access to the C++ options and classes. The comparison
wrappers -- `DescriptorComparison`, `FingerprintComparison`, `RMSDComparison`,
`ROCSComparison` and `SuperposeComparison` -- are on the top-level package, as
are most of the option structs (`PDistOptions`, `ButinaOptions`,
`FingerprintOptions` and the rest). Four are not: `ClusterReportOptions`,
`DescriptorOptions`, `DescriptorStatisticsOptions` and `RMSDOptions`. Reach
those through `oecluster.oecluster`, or let the wrapper build them from
keywords. To check the split against the version you have installed:

```python
import oecluster

raw = oecluster.oecluster
print([n for n in dir(raw)
       if n.endswith("Options") and not hasattr(oecluster, n)])
```

## Exceptions

Invalid arguments raise standard Python exceptions with a descriptive message.
Among the causes: an unknown comparison or representative method name and a
`"highest_neighborhood"` request without a `threshold` raise `ValueError`; a
sparse matrix passed where complete distances are required raises
`RuntimeError`, which is also how a refusal from the C++ layer usually
surfaces.

`TypeError` covers argument misuse, and is the usual outcome of the
explicitness rules described under [Fingerprint](#fingerprint): an argument
named where the rest of the configuration would ignore it, an argument of the
wrong type, and an unrecognized keyword all raise it.
