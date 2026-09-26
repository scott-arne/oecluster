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
are `"descriptor"`, `"fingerprint"`, `"mcs"`, `"rmsd"`, `"rocs"`,
`"sitehopper"`, and `"superpose"` (`"sitehopper"` is a superpose mode). See
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

Most clustering functions take a `DistanceMatrix`. Two do not: BitBirch takes an
OEFP fingerprint batch, and Murcko takes molecules directly, because it
partitions on chemical structure rather than on distance. All of them return a
result that subclasses `ClusteringResult`.

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
| `k_medoids(dm, ...)` | `DistanceMatrix` | `n_clusters`, `init`, `initial_medoids` |
| `bitbirch(fingerprints, ...)` | `oefp.OEFPBatch` | `threshold`, `branching_factor`, `merge_criterion` |
| `murcko(mols, ...)` | `list[OEMolBase]` | `scaffold` |

Algorithm-specific outputs live on the specific result subclass — for example
`DBSCANResult.core_sample_indices` or `BitBirchResult.centroids` — so a result
never carries fields that do not apply to its algorithm.

### k-medoids

`k_medoids()` places exactly `n_clusters` centers, each of which is a real
member of the input rather than a synthetic average, and chooses them by
minimizing the sum of every item's distance to its assigned center. It is the
algorithm to reach for when the cluster count is a requirement and the centers
have to be orderable compounds.

```python
result = oecluster.k_medoids(dm, n_clusters=10)

result.medoids        # (17, 42, ...) one item index per cluster, ascending
result.cost           # sum of each item's distance to its medoid
result.n_iterations   # swap iterations performed
result.converged      # True when no single swap lowers the cost
```

`init` selects the seeding strategy: `"build"` (the default; greedy PAM BUILD),
`"farthest_first"` (deterministic MaxMin from the global medoid, for
spread-out seeds), or `"explicit"` with `initial_medoids`. Passing a non-empty
`initial_medoids` with any other `init` raises, rather than silently deciding
which one you meant.

Output is byte-identical across runs, `num_threads` values and `chunk_size`
values; those two options change runtime and memory use but not the result.

`converged` is a real guarantee, not a loop-exit flag: when it is `True`, no
single medoid swap lowers the cost the result reports. That claim is verified
by recomputing each candidate total rather than by trusting the optimizer's
incremental arithmetic (or holds vacuously when `n_clusters` equals the item
count, which returns the identity partition). Reaching `max_iterations` is not
an error -- the partition and its medoids are valid -- but `converged` is then
`False` and no optimality claim is made.

That guarantee is not free. The verification pass recomputes a full objective
for each of the `n_clusters * (n - n_clusters)` candidate swaps, so it costs on
the order of `n^2 * n_clusters` distance lookups and on large inputs can take
longer than the swap phase it verifies. The cheaper incremental check was
rejected deliberately: recomputing from scratch is what makes the comparison
independent of the optimizer's own accumulated arithmetic.

One consequence worth knowing: at a converged solution every medoid
minimizes the distance sum within its own cluster, so
`representative(result.clusters[i], dm, method="medoid")` returns
`result.medoids[i]` whenever that minimum is unique. Where two members of a
cluster tie on within-cluster distance sum, both are minimizers and the two
functions may name different ones; the medoid the result reports is still a
minimizer.

### Murcko scaffolds

`murcko()` is the one clustering entry point that partitions on chemical
structure rather than on distance, so it takes molecules directly. Two molecules
share a cluster exactly when their Bemis-Murcko scaffolds canonicalize to the
same SMILES.

```python
result = oecluster.murcko(mols)

result.scaffolds          # ('c1ccccc1', 'c1ccccc1', '', ...) one per molecule
result.cluster_scaffolds  # ('c1ccccc1', 'c1ccncc1') sorted, indexed by label
result.labels             # -1 for a molecule with no ring system
```

`scaffold="generic"` reduces each framework to its topology -- every heavy atom
carbon, every bond single -- so scaffolds differing only in element or bond
order collapse together.

`murcko_scaffolds()` returns the same per-molecule strings without clustering,
which is what closes the loop with the functions that consume a scaffold
annotation:

```python
scaffolds = oecluster.murcko_scaffolds(mols)
dm = oecluster.pdist(mols, "fingerprint", metric="tanimoto")
clustering = oecluster.butina(dm, threshold=0.35)
agreement = oecluster.scaffold_agreement(clustering.labels, scaffolds)
```

`scaffolds` is per-item and in input order, so it is also what the
representative functions take as `scaffold_labels`:

```python
index = oecluster.representative(clustering.clusters[0], dm,
                                 method="weighted_medoid",
                                 scaffold_labels=scaffolds)
```

Molecules are taken as given: there is no salt stripping and no largest-component
selection. This matters less than it sounds, because an acyclic counter-ion has
no ring system and so contributes no scaffold region -- a hydrochloride, a
sodium salt and a mesylate all give the same scaffold as the free base. Only a
counter-ion that carries its own ring, such as a tosylate, adds a `.`-joined
component that will not match the free base. Stereochemistry is dropped and
explicit hydrogens are suppressed, so scaffold identity does not depend on how a
molecule was read, and hydrogen counts and formal charges are recomputed on the
atoms a sidechain cut touched, so diphenyl sulfone, diphenyl sulfoxide and
diphenyl sulfide share one framework scaffold. Only the cut atoms: a charged
ring atom whose bonds all survive keeps its charge, so an N-methylpyridinium
reduces to pyridine while a pyridinium stays a pyridinium.

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

One shipped combination is refused. `bitbirch_refine()` can return an emptied
leaf subcluster as an empty member list -- reference-parity behaviour for its
prune pass -- and `cluster_report()` refuses any result carrying an empty
cluster, with a `RuntimeError` reading `Cluster must contain at least one
member`. Nothing else the library produces trips that check.

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
> exist for a distance matrix, so each cluster's medoid stands in for its
> centroid. Calinski-Harabasz also has a grand-mean term, and the global medoid
> stands in for that; Davies-Bouldin has no such term. The values are therefore
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
with the largest single cluster. The two sorted arrays it reads hold every
pairwise distance among those points, `Nc * (Nc - 1) / 2` doubles in total,
which is roughly 400 MB at `Nc = 10,000` and 10 GB at `Nc = 50,000` -- but the
flag pays for only the between-cluster array. The within-cluster one is built on
every call, with or without the flag, because `median_intra_distance` is taken
over it, and taking that median copies the array transiently. So the flag adds
nothing to a single-cluster result and nearly the whole figure to one with small
clusters, which is the case it is off by default for. A failed allocation raises
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

## Partition Agreement

`partition_agreement()` scores two labelings of the same samples against each
other, and `scaffold_agreement()` scores a clustering against a per-sample
scaffold annotation. Neither takes a distance matrix, so comparing two methods
costs nothing beyond the labels they already produced:

```python
butina_result = oecluster.butina(dm, threshold=0.35)
dbscan_result = oecluster.dbscan(dm, eps=0.35, min_samples=5)

agreement = oecluster.partition_agreement(
    butina_result,
    dbscan_result,
    noise="excluded",
    adjusted_mutual_information=True,
)
print(agreement)
print(agreement.adjusted_rand_index, agreement.v_measure)
```

`scaffold_agreement()` takes one scaffold string per sample, in the same order
as the molecules the clustering was built from -- `scaffolds` below is the same
per-molecule list `rank_representatives()` accepts, so it has one entry for
every molecule in `mols`. A length mismatch raises `ValueError`:

```python
print(oecluster.scaffold_agreement(butina_result, scaffolds).completeness)
```

Either function accepts a clustering result, a list or tuple of ints, or a
numpy integer array on each label side. Labels are held natively as 32-bit
signed ints, so one outside that range raises `ValueError` naming the argument
rather than being truncated. Side A is the first argument:
`homogeneity` is `MI / H(a)` and `completeness` is `MI / H(b)`, and swapping
the arguments exchanges that pair, along with `num_clusters_a` and
`num_clusters_b`; every other metric is symmetric, with
`adjusted_mutual_information` symmetric only to within rounding, since a swap
transposes the contingency table and trades the two marginal values inside each
expected-MI term, changing the order those terms evaluate in. For
`scaffold_agreement()` the clustering is side A, so `completeness` is the
scaffold-purity reading -- whether each cluster's members share a single
scaffold -- and `homogeneity` is its transpose, whether each scaffold landed in
a single cluster.

`noise=` takes `"singletons"` (the default; each negatively-labelled sample
becomes its own cluster), `"grouped"` (each side's noise forms one cluster,
which is how scikit-learn reads a -1 label), or `"excluded"` (a sample noisy on
either side is dropped from both, which can make `num_samples` smaller than the
input length). An empty scaffold string is missing data, not a category, and
follows `noise=` exactly as a negative label does.

`adjusted_mutual_information=True` adds the seventh metric. It is opt-in
because the other six come essentially free once the contingency table is
built, while the expected-MI correction pays for a pass of its own: an O(N)
table of log factorials, then a sum over pairs of distinct marginal values --
distinct cluster sizes, not clusters -- each term walking the hypergeometric
support. `agreement.requested` records the request and not the outcome, the
same convention `cluster_report()` uses: `False` means nobody asked, `True`
with `nan` means asked and undefined. An opt-in metric nobody asked for reads
`None` from `to_table()` and `--` from `repr()`, while `nan` keeps its single
meaning of asked and undefined.

AMI's accuracy is limited where its denominator -- the mean entropy
minus the expected mutual information -- approaches zero, which happens when
both partitions are close to all-singleton. Numerator and denominator are then
each a difference of nearly equal sums over N terms, and the quotient loses
significance: two 1.5-million-sample partitions differing by one merged pair
have a true value of zero and report about 0.035. Partitions whose denominator
is order one are unaffected. The limit is inherent to computing the correction
in double precision -- scikit-learn shares it -- rather than a property of this
implementation.

An undefined metric is `nan` rather than a substituted value, which diverges
from scikit-learn on five degenerate inputs. This table is the complete list:

| Case | Fixture | scikit-learn 1.9.1 | OECluster |
| --- | --- | --- | --- |
| Fewer than two surviving samples | `a = b = [0]` | every metric `1.0`, except `fowlkes_mallows = 0.0` | every metric `nan` |
| Side A is one cluster, partitions differ | `[0,0,0,0]` vs `[0,0,1,1]` | `homogeneity = 1.0` | `nan` |
| Side B is one cluster, partitions differ | `[0,0,1,1]` vs `[0,0,0,0]` | `completeness = 1.0` | `nan` |
| One side all singletons, partitions differ | `[0,1,2,3]` vs `[0,0,1,1]` | `fowlkes_mallows = 0.0` | `nan` |
| Both sides all singletons, N >= 2 | `[0,1,2,3]` vs `[3,2,1,0]` | `fowlkes_mallows = 0.0` | `1.0` |

The last row is not a `nan` case: two all-singleton partitions are the same
partition, so the six always-computed metrics are 1.0 -- and
`adjusted_mutual_information` with them if it was asked for.

## SAR Coherence

Three functions ask whether a structure-activity relationship is there at all,
and whether a clustering captured it. `sar_coherence()` decomposes an activity
vector across a labeling you already have. `activity_landscape()` and
`modelability()` take a distance matrix instead of a labeling, so they answer
the question without clustering first:

```python
import math

activity = [...]   # one measurement per molecule, in the matrix's sample order

coherence = oecluster.sar_coherence(butina_result, activity)
print(coherence.eta_squared, coherence.omega_squared)

landscape = oecluster.activity_landscape(
    dm, activity, distance_threshold=0.30, activity_threshold=1.0)
print(landscape.num_cliffs, landscape.cliff_density, landscape.max_sali)

# The nan arm is not optional: nan >= 6.0 is False, so a comprehension
# without it labels every missing measurement "inactive".
classes = ["" if math.isnan(a) else "active" if a >= 6.0 else "inactive"
           for a in activity]
print(oecluster.modelability(dm, classes).modi)
```

All three treat a `nan` activity -- or, for `modelability()`, an empty class
string -- as missing rather than as a value. The sample is dropped, and every
result reports `num_samples`, the length of the input, beside `num_scored`,
how many of those samples entered the metric. For `activity_landscape()` and
`modelability()` that is exactly the count that carried a usable value.
`sar_coherence()` drops whatever `noise=` excludes as well, so under its
`"excluded"` default a gap between the two figures can be noise rather than
missing data: `sar_coherence([-1, 0, 0, 1, 1], [1.0, 2.0, 3.0, 4.0, 5.0])`
reports `num_samples` 5 and `num_scored` 4 with every activity usable. An
undefined metric is `nan`, never a substituted number, on the same convention
`cluster_report()` and `partition_agreement()` follow.

### `sar_coherence()`

`eta_squared` is the fraction of activity variance that falls between clusters
rather than within them, and `omega_squared` is the same quantity with the
variance a random labeling would explain removed. η² rises with the number of
clusters whether or not the clusters mean anything, so ω² is the one to compare
across labelings of different granularity. That comparison only means something
across labelings that score the same samples, and under the `"excluded"` default
two labelings of the same molecules need not: over activities
`[0., 0., 10., 10., 0., 10.]`, the labelings `[0, 0, 1, 1, 0, -1]` and
`[0, 0, 1, 1, -1, 0]` both report `num_scored` 5 and two clusters yet read ω²
1.0 and 0.21875, having scored different fives. Equal `num_scored` does not
establish that two labelings scored the same samples. ω² can go below zero,
which says the labeling explains less than chance would; η² cannot. Read
"chance" there as approximately rather than exactly zero. Two conditions
together lift the chance level: the clusters are mostly singletons, so
`num_scored - num_clusters` is small, and the activity variance is concentrated
in a few samples. Over activities `[1.0] + [0.0] * 99`, a labeling of one
11-member cluster plus 89 singletons has an exact chance expectation of 0.075;
five actives among 95 inactives give 0.014 at the same shape. Neither condition
acts alone -- under balanced clusterings of two, four, ten, twenty and fifty
clusters those same two activities have exact expectations of zero and of
0.0000017 through 0.00016. Distinct values are not the protective property: the
100 distinct activities `[1.0] + [i * 1e-6 for i in range(1, 100)]` reach 0.075
again at the singleton-heavy shape, because what matters is where the variance
sits, not whether the values are tied or discrete. A tight `butina()` threshold
produces just such a singleton-heavy labeling, so on screening data read a
small positive ω² as unresolved rather than as weak signal. Both are
`nan` when fewer than two samples are scored, and both are `nan` when the
scored activities have zero variance -- there is nothing to apportion. ω² has
one further `nan` case of its own: when every scored sample is its own cluster
there are no within-cluster degrees of freedom for the chance correction to
divide by. η² does not share that case. On an all-singleton labeling it reads
1.0 whenever the scored activities vary at all, and `nan` when they do not, the
zero-variance rule above outranking everything else.

`coherence.clusters` is a tuple of `ClusterActivity` rows -- `label`,
`num_scored`, `mean_activity`, `stddev_activity` -- and it feeds
`pandas.DataFrame(coherence.clusters)` directly. There is a row for every
cluster that kept at least one scored sample, and no row for the others: a
cluster whose every member was dropped as noise or as a missing measurement is
absent from the table rather than present with `nan` statistics, and
`num_clusters` counts the rows that are there. A row that is present is not
therefore a row with complete statistics -- `stddev_activity` is `nan` for
every cluster down to one scored member, a spread over one sample not being
defined -- so `dropna()` over this table deletes every cluster with one scored
member, which need not be a cluster with one member, and not the absent ones.
Nor is `label` a key: under `noise="singletons"` each promoted noise point
becomes its own row carrying the negative label it arrived with, and nothing
stops two of those from being equal, so
`sar_coherence([-1, -1, 0, 0], [1.0, 2.0, 3.0, 4.0], noise="singletons")`
returns three rows labelled `[-1, -1, 0]`, where the same call over
`[-1, -2, 0, 0]` returns `[-1, -2, 0]` and repeats nothing. A dict
comprehension keyed by `label` silently keeps the last of any repeat and
reports two clusters where there are three, and `set_index("label").loc[-1]`
hands back a frame rather than the row the caller was reaching for. Position
distinguishes those rows; the label does not. The rows follow the order in
which their labels first appear among the scored samples, not sorted label
order, so a clustering result and a reordering of the same labels can list the
same clusters in different positions. "Among the scored samples" is the
operative part: with labels `[1, 0, 1]` and activities `[nan, 2.0, 3.0]` the
rows come out `[0, 1]`, because the first `1` carries no measurement and so
does not get to claim the first position for its cluster.

The first argument is a clustering result, a list or tuple of ints, or a numpy
integer array, on the same terms as `partition_agreement()`. `noise=` takes the
same three spellings -- `"singletons"`, `"grouped"`, `"excluded"` -- but here
the default is `"excluded"`: noise is not a structural hypothesis, and counting
each noise point as its own cluster inflates η² for a reason that has nothing
to do with the labeling under test.

### `activity_landscape()`

`activity_landscape()` sweeps every pair of scored samples once.
A pair is an activity cliff when the two are close and their activities are
far apart: `distance <= distance_threshold` and `abs(delta) >=
activity_threshold`, both inclusive, both in the caller's own units.
`cliff_density` is `num_cliffs / num_pairs_scored`, where `num_pairs_scored` is
every pair of scored samples, so the reading is per-pair and comparable across
sets of different size.

`max_sali` and `mean_sali` summarise the SALI ratio `abs(delta) / distance`
over the pairs where that ratio is defined. At zero distance it is not: the
ratio is infinite where the two activities differ and `0/0` where they agree.
Zero-distance pairs are therefore excluded from both figures either way and
counted in `num_zero_distance_pairs` instead -- a non-zero count there is worth
reading before the SALI figures, because a zero-distance pair whose activities
differ sharply is the sharpest cliff there is and neither figure sees it. Such
a pair is still counted as a cliff, on the same
`abs(delta) >= activity_threshold` test as any other -- a zero distance clears
every permitted `distance_threshold` -- so `num_cliffs` and `cliff_density`
include it.

`rmodi` scores the landscape without a distance threshold at all. For each
scored sample it finds the nearest neighbour whose activity is inside the
sample's band and the nearest one outside it. A neighbour is inside the band
when its activity differs from the sample's by at most `rmodi_delta *
activity_stddev`, the comparison being inclusive; `rmodi_delta` is thus a
half-width in standard deviations, reaching that far either side of the
sample's own activity, and the band spans twice it. `activity_stddev` is the
population standard
deviation of the scored activities and is reported on the result, so the band
is auditable. `rmodi` is the fraction of samples whose in-band neighbour is
strictly nearer than the out-of-band one. The default `rmodi_delta=0.625` is
the published one.

Both extremes are reachable from the band width alone, so read `rmodi` against
`activity_stddev` rather than as a verdict on the landscape. A sample with no
in-band neighbour has its in-band minimum left at infinity and counts as
discordant, so if every pair's activity difference exceeds the band then every
sample is discordant and `rmodi` is 0.0 however faithfully distance tracks
activity. Five points whose pairwise distances equal their activity
differences exactly do that when they are also equally spaced: activities
`[1, 2, 3, 4, 5]` score 0.0 at the default `rmodi_delta` and 1.0 at
`rmodi_delta=1.0`. Read that as one five-point fixture rather than as a rule
about even spacing. The band is scaled to the whole spread, so lengthening the
run widens it while leaving the gap alone, and `[1, 2, 3, 4, 5, 6]` already
reads 1.0 at the default; the same construction over `[1, 2, 4, 8, 16]` reads
0.6 and 0.8 at those two widths. At the other end, flat
activity collapses the band to zero width but puts every pair inside it, no
sample has an out-of-band neighbour, and `rmodi` is 1.0 from two scored samples
up -- correctly, since nothing about the activity contradicts the distances --
and `nan` below that. Sweep `rmodi_delta` before concluding anything from a
value at either end.

### `modelability()`

`modelability()` asks the same question of a class annotation: `modi` is the
mean, over the classes, of the fraction of each class's members whose nearest
scored neighbour carries the same class. Averaging per class rather than per
sample is what keeps a large class from burying a small one. `report.classes`
is a tuple of `ClassConcordance` rows -- `label`, `num_members`,
`fraction_same_class` -- and also feeds `pandas.DataFrame` directly. As with
`coherence.clusters`, there is a row for every class that kept at least one
scored sample and no row for the others, ordered by first appearance among the
scored samples rather than sorted by label. The distinction that bites
`coherence.clusters` does not arise here: the only unscored sample is one whose
class string is empty, and that is no class at all, so this is plain input
order over the classes that exist. `modi` is `nan` when fewer than two classes
are scored: with one class no neighbour could carry a different one, so the
comparison has no content. The concordance sweep is skipped rather than run to
a foregone answer -- the distances are still validated -- so the single row's
`fraction_same_class` is `nan` as well, not 1.0.

Nearest-neighbour ties resolve to the lowest scored index, so `modi` and the
`fraction_same_class` rows behind it depend on the order the samples arrive in
wherever distances tie: over a matrix whose off-diagonal distances are all
equal, the classes `["A", "A", "B", "B"]` read 0.5 where `["A", "B", "A", "B"]`
read 0.25 on the same four samples, and across the six orderings of that data
`modi` spans 0.25 to 0.5 where each class's fraction spans the whole 0.0 to
1.0. Two duplicate molecules are enough to produce such a tie, being
equidistant from every third sample, so reordering rows can move these figures
on a real matrix as well.

`modi` is not chance-corrected, and its chance level is not zero: random labels
over N scored samples in K classes average about `(N - K) / (K * (N - 1))`,
which approaches `1/K` as N grows and does not move with the class balance.
Read a `modi` against that figure rather than against 0.0.

The annotation is a sequence of strings, one per sample, and an empty string is
missing data rather than a category -- the same reading `scaffold_agreement()`
gives it. A bare `str` is rejected rather than iterated, so `"AAB"` raises
`TypeError` instead of quietly becoming three annotations.

### What the matrix functions require

`activity_landscape()` and `modelability()` check their distance matrix, but
against a weaker standard than the clustering entry points use. They rank
distances against one another and compare them to a threshold; neither
operation needs the triangle inequality, so a matrix a clustering algorithm
would refuse is usually fine here and there is no `allow_nonmetric` parameter
to pass. What they do refuse is sparse storage, a similarity, a measure whose
self-distance does not vanish, a non-finite distance, and -- alone among the
entry points in refusing it outright, with no override -- a matrix stamped
`data_integrity == 'subset_scored'`.

Sparse storage is the refusal most callers meet first, because `pdist(...,
cutoff=...)` produces it as a matter of course. Both functions read every
pairwise distance among the scored samples, which a matrix that kept only the
distances below a cutoff cannot supply, so both raise `ValueError` from the
Python layer before the gate runs. The rest are the gate's:

```python
dm = oecluster.pdist(mols, "descriptor", missing="ignore")
oecluster.activity_landscape(dm, activity)   # ValueError
```

A matrix built with `missing='ignore'` scores each pair on whatever features
that pair happens to share, so two of its distances answer different questions
and the smaller one is not necessarily the nearer pair. Ranking is the whole
operation here, so this one is not overridable; recompute with
`missing='complete_case'`.

### Which exception you get

Four exception types are in play, not two. The order they are described in
below is the order of this prose and not a precedence rule, and the three
functions do not share one order between them either. `activity_landscape` and
`modelability` check their first argument's type, then whether its storage is
sparse, then everything else, so `activity_landscape(sparse_dm, "abc")` reports
the sparse `ValueError` and never looks at the activity. The `subset_scored`
refusal is not in that early position: the activity and the options are read
first, so against a `missing='ignore'` matrix a bad one of those is what you
see, and the type does not tell them apart either -- that refusal, a length
mismatch and a rejected option all raise `ValueError`. `sar_coherence` checks
the activity's type, whether its elements convert to a double, and whether it
is empty, all before it judges `result`, so `sar_coherence(3.5, "bad")` names
the activity rather than the result and `sar_coherence(3.5, [])` raises
`ValueError` over a `result` `TypeError` already pending. The length check runs
the other way about, after `result` rather than before it. Catch on the
exception type rather than on the check you expect to run first.

`TypeError` is for an argument of the wrong kind rather than the wrong value:
a first argument that is not a `SymmetricDistanceMatrix`, a `result` that is
neither a clustering result nor an iterable of ints, an activity element that
is not a number, and a class element that is not a string. A bare `str` and a
mapping are both refused rather than iterated in any of the three iterable
positions -- `result`, the activity, the class annotation. What is not checked
in any of them is ordering: a `set` of labels is accepted, and its members are
paired to the activities in the set's own iteration order, which need not be
the order they were written in. `{5, 3, 1, 0}` gives label 0 the first activity
where `[5, 3, 1, 0]` gives it the last, and neither is an error, so pass a
sequence when the pairing matters. The same silence in the activity position
costs more than a pairing: `sar_coherence([0, 0, 1, 1], [10., 20., 30., 40.])`
reports an `eta_squared` of 0.8, and the same call with those four values as a
`set`, which iterates `40.0, 10.0, 20.0, 30.0`, reports 0.0. All four numeric
options behave alike here, each being coerced before it is range-checked:
`None` or a list raises `TypeError` from that coercion for any of
`distance_threshold`, `activity_threshold`, `rmodi_delta` and `num_threads`,
while a string the coercion cannot parse raises `ValueError` from the same
call.

Each of the three raises `ValueError` for whatever its own signature lets
Python see for itself: a length mismatch, an empty labeling, an empty activity,
an out-of-range label, an unknown `noise=` spelling, a non-finite or negative
threshold, a `num_threads` of -1 or below (`-0.5` truncates to 0 and is
accepted, selecting the hardware concurrency rather than raising), and every
refusal above from the matrix check. The two overloads of `sar_coherence()`
agree on all of the conditions they share, which is worth saying because they
reach the native layer by different routes: an empty activity and a label too
wide for the native 32-bit type are refused in Python on both routes, not left
to whichever one happens to notice. The emptiness check runs ahead of the
length check in all three functions, and that ordering is the whole of what
makes the promise true: a zero-sample matrix and an empty clustering are both
constructible, and against an empty annotation a length check compares 0 with 0
and agrees. Leave it out of any one of the three and that function alone
reports the empty case as `RuntimeError`. Conditions only the C++ sweep can
reach arrive as `RuntimeError` instead -- an infinite activity value, a
negative stored distance, a mean or an accumulator that overflows to
infinity -- because the bindings map every native exception to `RuntimeError`,
the same way `cluster_report()` already does. The negative distance is the one
worth knowing about: the matrix check tests that entries are finite, which a
negative number is, so nothing catches it until the sweep reads it.

The fourth type is `OverflowError`, and it derives from `ArithmeticError`, so
an `except ValueError` does not catch it. Only `num_threads` reaches it, in
both functions that take one, by two routes: `float("inf")` fails in `int()`
with "cannot convert float infinity to integer", and a finite value past the
native `size_t` -- `1e300`, or `10**30`, though not `2**64 - 1` -- fails where
the bindings assign it to the options struct. The three `float()` options
cannot reach it at all, each being finiteness-guarded before use, so
`distance_threshold=float("inf")` is
`ValueError: distance_threshold must be finite`; and a `nan` is a `ValueError`
for all four.

## Metric Requirements

`butina()`, `dbscan()`, `hdbscan()`, `agglomerative()`, and `cluster_report()`
assume their input is a metric: an item's distance to itself is zero, and the
triangle inequality holds. Those five are the entry points that check.
`k_medoids()`, `activity_landscape()` and `modelability()` check a weaker
standard described below -- they rank, threshold and add distances but never
assume the triangle inequality. `k_medoids()` is the one *clustering* algorithm
in that weaker group, which is not an oversight: PAM's objective is a sum of
distances and its swap step compares two such sums, so nothing in it appeals to
the triangle inequality and it takes no `allow_nonmetric` parameter because
there is no assumption for a flag to override.

That produces one asymmetry worth expecting. A Dice matrix that `k_medoids()`
clusters without complaint will be refused by `cluster_report()` unless you pass
`allow_nonmetric=True`, because the internal validity indices do lean on metric
behavior where PAM does not. Both calls are behaving correctly; it is the
sequence that surprises.

`representative()`, `rank_representatives()`, and
`select_representatives()` also take a distance matrix, and consult none of
these facts. Not every comparison produces a metric, so each matrix records
what its comparison actually guarantees:

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

`activity_landscape()` and `modelability()` take the three refusals above as
they stand, waive the two triangle-inequality warnings without being asked --
which is why they have no `allow_nonmetric` parameter -- and promote
`missing='ignore'` from a warning to a refusal they will not waive. See
[SAR Coherence](#sar-coherence) for why ranking incomparable distances has no
correct reading.

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

### MCS

| Parameter | Values | Default |
|-----------|--------|---------|
| `search_mode` | `"approximate"`, `"exhaustive"` | `"approximate"` |
| `match_level` | `"default"`, `"exact"`, `"loose"` | `"default"` |
| `max_matches` | Matches one directed search may enumerate | `1024` |

Scores the maximum common substructure as Tanimoto over matched bonds:
`c / (|A| + |B| - c)` with `c` the matched-bond count and `|A|`, `|B|` the two
molecules' heavy-atom bond counts. The distance is one minus that.
`similarity=True` is supported, unlike `rmsd`.

Coordinates are never read, so molecules parsed from SMILES need no embedding
step and a multi-conformer `OEMol` is scored once rather than once per pose.
Hydrogens are always suppressed, isotopic ones included, so a deuterated
analogue scores as identical to its parent. Molecules with no bonds after
suppression -- methane, water, argon -- are refused at construction, because
bond Tanimoto has a zero denominator for them rather than an extreme value.

`match_level` chooses how strictly atoms and bonds must correspond: `"loose"` is
atomic number only with bonds unconstrained, so benzene matches cyclohexane
completely; `"default"` is OEChem's own pair of expressions, under which those
two share nothing; `"exact"` adds hydrogen count, charge, degree and bond order.

`search_mode` deliberately inverts the toolkit's own default. Exhaustive search
genuinely finds larger matches on rigid polycyclic and sugar-like input -- on
eleven of 120 measured pairs it did, by one to three bonds -- but it costs one
to three orders of magnitude, and it is not uniformly better: on a 53-bond
against 54-bond macrolide pair it took 16.2 s and matched 50 bonds where
approximate took 8.5 ms and matched 51. Exhaustive mode also prints
`Warning: MCS search truncated` to OpenEye's process-global error stream, once
per truncated pair, which at 500,000 pairs is unusable output; the library does
not redirect that stream, because doing so would silence warnings from the
caller's own unrelated OpenEye code.

There is **no metric guarantee**. No triangle-inequality violation appeared in
74,400 ordered triples across three molecule sets, one of them built to stress
transitivity, and the thinnest observed margin was 1.0000 against 1.1444. But
the inclusion-exclusion bound that a genuine set intersection satisfies was
violated 66 times over 59,280 triples, so the Jaccard metric proof is
unavailable rather than merely unattempted. The matrix therefore reports
`triangle` as `"unknown"`, which every clustering entry point accepts without
`allow_nonmetric=True`. To check the property empirically on your own data, use
`from_array(..., probe_triples=N)`, which samples triples and refuses on a
violation.

**Cost.** An order-of-magnitude planning estimate, not a measurement. No
`pdist` run at this size has been executed.

| quantity | 1,000 molecules, approximate, drug-like |
|----------|------------------------------------------|
| pairs / searches | 499,500 / 999,000 |
| search time | ~400 s |
| search construction, 0.287 ms each | ~287 s |
| total, single-threaded | ~11.5 min |
| total, 8 threads | ~90 s |
| peak memory added by per-clone molecule copies | ~79 MB |

Two caveats travel with those numbers. The per-search figures come from small
drug-like inputs and scale steeply with molecule size -- the erythromycin
against azithromycin pair above is 8.5 ms in approximate mode, about twelve
times the 0.691 ms mean the microbenchmark measured. And the 8-thread row
assumes linear scaling, which was never measured. Use `fingerprint` for large
sets and `mcs` for focused series.

Memory grows as `num_threads x n x 10.4 KB`, because each worker holds private
copies of the molecules rather than sharing one set. That is the one cost here
that grows with thread count instead of shrinking.

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

Reached directly, `Get` and `Set` refuse a pair they cannot address: an index
at or beyond `NumSamples()` on either, and the diagonal on `Set`, which owns no
stored slot because `Get` answers `i == j` from a shortcut rather than from
memory. Both surface as `RuntimeError`.

## Advanced C++ Binding Access

Most users do not need this section. The generated SWIG wrapper is available as
`oecluster.oecluster` and the compiled extension as `oecluster._oecluster` for
users who need direct access to the C++ options and classes. The comparison
wrappers -- `DescriptorComparison`, `FingerprintComparison`, `MCSComparison`,
`RMSDComparison`, `ROCSComparison` and `SuperposeComparison` -- are on the
top-level package, and
so are many of the option structs, among them `PDistOptions`, `ButinaOptions`
and `FingerprintOptions`. Others are not, `ClusterReportOptions` and
`RMSDOptions` among them. Reach those through `oecluster.oecluster`, or let the
wrapper build them from keywords. Which structs fall on which side moves as new
metrics are added, so read the split off the installed package rather than off
a list here:

```python
import oecluster

raw = oecluster.oecluster
print([n for n in dir(raw)
       if n.endswith("Options") and not hasattr(oecluster, n)])
```

## Exceptions

Invalid arguments raise standard Python exceptions with a descriptive message.
Among the causes: an unknown comparison or representative method name and a
`"highest_neighborhood"` request without a `threshold` raise `ValueError`. A
sparse matrix passed where complete distances are required splits by where the
refusal lives. `hdbscan`, `agglomerative`, `cluster_report`,
`activity_landscape` and `modelability` check the storage
in the Python layer and raise `ValueError`; `rank_representatives` and
`select_representatives` have no such pre-check and refuse from the C++ layer
with `RuntimeError`, which is also how a refusal from that layer usually
surfaces.

`cluster_report` raises `MemoryError` where the C++ layer reports a memory
condition rather than a bad argument: a failed allocation in the pair-rank
stage, and a pair or couple count that would exceed the integer holding it.

`TypeError` covers argument misuse, and is the usual outcome of the
explicitness rules described under [Fingerprint](#fingerprint): an argument
named where the rest of the configuration would ignore it, an argument of the
wrong type, and an unrecognized keyword all raise it.
