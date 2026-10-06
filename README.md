# OECluster

Molecular clustering, distance-matrix computation, and representative selection
for cheminformatics workflows built on [OpenEye Toolkits](https://www.eyesopen.com/).

OECluster gives Python, C++, and command-line users a fast path from molecules to
clusters, ranked representatives, and reusable pairwise distance matrices. It is
designed for medicinal chemists and cheminformaticians who need practical
cluster summaries as well as lower-level control over distance computation.

---

## Table Of Contents

- [What You Can Do](#what-you-can-do)
- [Requirements](#requirements)
- [Installation](#installation)
- [Quickstart: Cluster Molecules From SMILES](#quickstart-cluster-molecules-from-smiles)
- [Core Python Workflow](#core-python-workflow)
- [Choosing A Clustering Algorithm](#choosing-a-clustering-algorithm)
  - [From The Command Line](#from-the-command-line)
- [Choosing Representatives](#choosing-representatives)
- [Assessing And Comparing Clustering Quality](#assessing-and-comparing-clustering-quality)
- [Scaling Guidance](#scaling-guidance)
- [Comparison Methods](#comparison-methods)
- [Storage Backends](#storage-backends)
- [Command-Line Tool](#command-line-tool)
- [C++ API](#c-api)
- [Troubleshooting](#troubleshooting)
- [Examples](#examples)
- [License](#license)

---

## What You Can Do

- **Cluster molecular collections** with Butina, DBSCAN, HDBSCAN,
  agglomerative clustering, k-medoids, BitBirch, or Murcko scaffold
  clustering.
- **Cluster without a matrix**: Butina, DBSCAN and Butina-order sphere
  exclusion run straight from a comparison and keep only the threshold
  neighbor graph, behind an exact memory guard, so a large fingerprint set
  never needs its N x N distance matrix.
- **Choose representatives** with true medoids, minimax/radius centers,
  highest-neighborhood Butina-style representatives, weighted medoids, ranked
  representative lists, and k-representative selection.
- **Choose a clustering parameter from data** with `select_parameter`, which
  sweeps one parameter over a grid, scores every partition, and picks the
  best validity index under noise and cluster-count bounds.
- **Measure cluster stability** with `cluster_stability`, which reruns a
  spec on resampled subsets and reports per cluster how often it dissolved
  or was recovered (Hennig's clusterboot statistics), plus one adjusted Rand
  index per resample.
- **Combine an ensemble of partitions** with `consensus`, which builds the
  co-association matrix behind a sweep or a stability run and reads one
  agreed partition off it, with Monti's per-cluster evidence scores.
- **Drive all of that from a shell** with the `oecluster` command, which
  clusters, sweeps a parameter, scores stability and reaches a consensus over
  a precomputed distance matrix, and writes the result as CSV or JSON.
- **Compute molecular distances** for fingerprints, ROCS shape/color overlay,
  protein superposition, and binding-site comparison.
- **Scale distance storage** with dense in-memory arrays, memory-mapped files,
  or sparse cutoff-filtered storage.
- **Use the same core from Python, C++, or CLI workflows**.

---

## Requirements

- **Python** 3.11+ with NumPy 1.20 or later.
- **OpenEye Toolkits** 2026.1 or later.
- **A valid OpenEye license** at build time and runtime.
- **OEFP** 0.3.0 exactly, for fingerprint generation and comparison. oecluster
  compiles OEFP's core into its own extension and exchanges raw
  fingerprint-batch pointers with the separately compiled `oefp` wheel, so the
  compiled-against and installed versions must be identical. The extension
  compares the two the first time such a pointer crosses -- passing an
  `oefp.OEFPBatch` to `bitbirch`, for instance -- and raises `ImportError` on
  any difference.
- **pyarrow** 25.x. The extension loads Arrow and Parquet out of the installed
  `pyarrow` package by versioned library name rather than bundling its own
  copy, so the pyarrow major is an ABI dependency, the same one `oefp` carries.
- **C++17**, **CMake** 3.21+, and **SWIG** 4.0+ when building from source.

---

## Installation

Install a wheel when one is available:

```bash
pip install oecluster
```

Verify that Python can import OECluster and the OpenEye runtime can parse a
molecule:

```bash
python - <<'PY'
from openeye import oechem
import oecluster

mol = oechem.OEGraphMol()
assert oechem.OESmilesToMol(mol, "CCO")
print(f"oecluster {oecluster.__version__} is ready")
PY
```

Build a Python wheel from source when a wheel is not available:

```bash
pip install scikit-build-core
python scripts/build_python.py --openeye-root /path/to/openeye/toolkits
```

Build the C++ library, Python extension, tests, and CLI directly with CMake:

```bash
cmake -B build \
  -DCMAKE_BUILD_TYPE=Release \
  -DOPENEYE_ROOT=/path/to/openeye/toolkits
cmake --build build
cd build && ctest --output-on-failure
cmake --install . --prefix /usr/local
```

Useful build options:

| Option                    | Default | Description                      |
|---------------------------|---------|----------------------------------|
| `OECLUSTER_BUILD_TESTS`   | ON      | Build C++ tests                  |
| `OECLUSTER_BUILD_PYTHON`  | ON      | Build Python SWIG bindings       |
| `OECLUSTER_BUILD_TOOLS`   | ON      | Build the `oepdist` CLI          |
| `OECLUSTER_UNIVERSAL2`    | OFF     | macOS universal2 binary build    |

---

## Quickstart: Cluster Molecules From SMILES

This example needs no input files. It creates molecules from SMILES, computes
Morgan/Tanimoto distances, runs Butina clustering, and prints one medoid per
cluster.

```python
from openeye import oechem
import oecluster

records = [
    ("benzene", "c1ccccc1"),
    ("toluene", "Cc1ccccc1"),
    ("phenol", "Oc1ccccc1"),
    ("ethanol", "CCO"),
    ("propanol", "CCCO"),
    ("acetic_acid", "CC(=O)O"),
]

mols = []
for title, smiles in records:
    mol = oechem.OEGraphMol()
    if not oechem.OESmilesToMol(mol, smiles):
        raise ValueError(f"Could not parse {title}: {smiles}")
    mol.SetTitle(title)
    mols.append(mol)

dm = oecluster.pdist(
    mols,
    "fingerprint",
    fp_type="morgan",
    metric="tanimoto",
)
result = oecluster.butina(dm, threshold=0.55)

for cluster_id, cluster in enumerate(result.clusters):
    ranked = oecluster.rank_representatives(cluster, dm, method="medoid")
    medoid = ranked[0]
    print(
        cluster_id,
        mols[medoid.member].GetTitle(),
        len(cluster),
        medoid.metrics.cluster_radius,
    )
```

Run the same workflow from the tested example file:

```bash
python examples/quickstart_smiles.py
```

Expected output resembles:

```text
clusters: 5
representatives:
  cluster 0: size=2 medoid=propanol radius=...
```

---

## Core Python Workflow

### 1. Load Molecules

```python
from openeye import oechem

mols = []
ifs = oechem.oemolistream("molecules.sdf")
mol = oechem.OEGraphMol()
while oechem.OEReadMolecule(ifs, mol):
    mols.append(oechem.OEGraphMol(mol))
```

### 2. Compute A Distance Matrix

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

Use `cutoff` when you only need short distances and want sparse storage:

```python
dm = oecluster.pdist(mols, "fingerprint", metric="tanimoto", cutoff=0.35)
```

Use `output` when dense distances are large enough to keep on disk:

```python
dm = oecluster.pdist(mols, "fingerprint", output="distances.mmap")
```

**Cross-distance** (`cdist`) computes the rectangular NxM distance matrix between two distinct item sets:

```python
queries = [...]  # N query molecules
targets = [...]  # M target molecules
cross_dm = oecluster.cdist(queries, targets, "fingerprint", metric="tanimoto")
# cross_dm is a CrossDistanceMatrix with shape (N, M)
print(cross_dm.matrix[0, :])  # distances from first query to all targets
```

Use `cdist` when you need distances from one set (e.g., virtual screening hits) to another (e.g., in-house compounds), rather than the symmetric within-set distances that `pdist` computes.

### 3. Cluster And Select Representatives

```python
result = oecluster.butina(dm, threshold=0.35, reordering=False)

# result.labels is a length-n, scikit-learn-style assignment array, ready to
# drop into a DataFrame column (e.g. df["cluster"] = result.labels).
for cluster in result.clusters:
    medoid = oecluster.representative(cluster, dm, method="medoid")
    minimax = oecluster.representative(cluster, dm, method="minimax")
    ranked = oecluster.rank_representatives(cluster, dm, method="medoid")
    selected = oecluster.select_representatives(
        cluster,
        dm,
        k=2,
        method="medoid",
        selection="diversity",
    )
    print(medoid, minimax, ranked[0].metrics.cluster_radius, selected)
```

### 4. Export Cluster Assignments

```python
import csv

with open("clusters.csv", "w", newline="") as handle:
    writer = csv.writer(handle)
    writer.writerow(["cluster_id", "member_index", "title", "representative_rank"])
    for cluster_id, cluster in enumerate(result.clusters):
        ranks = {
            item.member: item.metrics.representative_rank
            for item in oecluster.rank_representatives(cluster, dm, method="medoid")
        }
        for member in cluster:
            writer.writerow([
                cluster_id,
                member,
                mols[member].GetTitle(),
                ranks[member],
            ])
```

---

## Choosing A Clustering Algorithm

| Algorithm       | Input                      | Behavior                                                                                                                                                                                                                                            | Best For                                                                                                                                                              | Key Parameters                                                |
| --------------- | -------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------- |
| `butina`        | `oecluster.DistanceMatrix`, comparison, items | Finds the molecule with the largest number of neighbors within the distance threshold, forms a cluster around that molecule, removes those molecules, and repeats. With `reordering`, neighbor counts are updated after each cluster is formed.     | Datasets with compact threshold-neighborhoods, interpretable local similarity groups, and a moderate number of singleton or near-singleton compounds.                 | `threshold`, `reordering`                                     |
| `dbscan`        | `oecluster.DistanceMatrix`, comparison, items | Identifies core molecules that have at least `min_samples` neighbors within `eps`, then expands clusters through density-connected neighborhoods. Molecules that are not density-connected to any cluster are labeled as noise.                     | Datasets with dense clusters separated by sparse regions, where a single global distance threshold is meaningful across the dataset.                                  | `eps`, `min_samples`                                          |
| `hdbscan`       | `oecluster.DistanceMatrix` | Builds a hierarchy of density-connected clusters using mutual-reachability distance, then extracts the most stable clusters from that hierarchy. Molecules in low-density or ambiguous regions can be assigned as noise or given weaker membership. | Datasets with variable-density chemical series, uneven cluster sizes, and a meaningful population of sparse or ambiguous molecules.                                   | `min_cluster_size`, `min_samples`, `cluster_selection_method` |
| `agglomerative` | `oecluster.DistanceMatrix` | Starts with each molecule as its own cluster and repeatedly merges the closest clusters according to the selected linkage rule. Clustering stops when either `n_clusters` is reached or the `distance_threshold` is exceeded.                       | Small-to-medium datasets where all molecules should be assigned to compact clusters and the desired granularity is controlled by cluster count or distance threshold. | `n_clusters`, `distance_threshold`, `linkage`                 |
| `k_medoids`     | `oecluster.DistanceMatrix` | Places exactly `n_clusters` medoids -- real members of the input, never synthetic averages -- and minimizes the sum of every item's distance to its assigned medoid, using PAM BUILD seeding and the FastPAM1 swap. Does not assume a metric.       | Datasets where the cluster count is a requirement and every center must be an orderable compound.                                                                     | `n_clusters`, `init`, `initial_medoids`                       |
| `sphere_exclusion` | `oecluster.DistanceMatrix`, comparison, items | Leader clustering or Directed Sphere Exclusion (DISE) with a fixed threshold. Each center claims every unclaimed item within the threshold, so centers are farther apart than the threshold. Runs lazily over a comparison for every order.         | Datasets where a threshold-based seed order is desired (leader, DISE, or Butina's neighbor order), and centers must be farther apart than the threshold.            | `threshold`, `order`, `assignment`                            |
| `jarvis_patrick` | `oecluster.DistanceMatrix` (dense or sparse), comparison, items, `KNNGraph` | Classic Jarvis-Patrick: two items link when each is among the other's k nearest neighbors and they share at least kmin of them; clusters are the connected components. Builds the graph with `knn_graph`, which is also public. | Density-adaptive clustering without a distance threshold, where cluster membership should follow shared neighborhoods. | `k`, `kmin` |
| `leiden` | `oecluster.DistanceMatrix` (dense or sparse), comparison, items, `KNNGraph` | Leiden community detection on a shared-nearest-neighbor graph: kNN edges weighted by the Jaccard overlap of the two neighborhoods, then partitioned by modularity or the Constant Potts Model. Every cluster is connected, and a seed makes runs repeatable. | Graph-based clustering at a chosen resolution, the approach single-cell tools such as Seurat and Scanpy use. | `k`, `objective`, `resolution`, `seed` |
| `bitbirch`      | `oefp.OEFPBatch`           | Incrementally inserts binary fingerprints into a Birch-style tree that summarizes nearby fingerprints in feature space. Leaf subclusters can then be merged according to the selected merge criterion to produce final clusters.                    | Large binary fingerprint datasets with many locally similar molecules where scalable feature-space clustering is preferred over a full pairwise distance matrix.      | `threshold`, `branching_factor`, `merge_criterion`            |
| `murcko`        | `list[OEMolBase]`          | Assigns each molecule its Bemis-Murcko scaffold -- ring systems plus the linkers that connect them -- and clusters by scaffold identity. Molecules with no ring system are labeled as noise. Takes molecules rather than a distance matrix, because the partition is on structure rather than on distance.                        | Chemical-series partitioning, scaffold-diversity accounting, and producing the `scaffold_labels` that `scaffold_agreement` and the weighted-medoid representative consume.                                                             | `scaffold`                                                    |

Start with **Butina** or **agglomerative** for familiar fingerprint threshold clustering. Use
**k-medoids** when the cluster count is fixed in advance and each center has to be a real molecule. Use
**DBSCAN/HDBSCAN** when noise and density matter. Use **BitBirch** when the workflow already has
OEFP dense binary fingerprints and you want feature-space clustering without
materializing an initial pairwise distance matrix (good for VERY large datasets).
Use **Murcko** when the grouping you want is chemical series rather than
fingerprint neighborhood -- it is the only algorithm here that partitions on
structure, and its clusters come with the scaffold string that names each one.

All clustering functions return a result that subclasses `ClusteringResult`,
which exposes read-only `labels` (a length-n, scikit-learn-style assignment
array) and `clusters` (a tuple of member-index tuples). Results support
`len(result)` (number of clusters), iteration over clusters, and `result[i]`
indexing. Algorithm-specific outputs live on the specific subclass (for
example `DBSCANResult.core_sample_indices` or `BitBirchResult.centroids`), so a
result never carries fields that do not apply to its algorithm.

### BitBirch Variants

Beyond the standard `bitbirch` function, two specialized clustering strategies are available for binary fingerprint batches:

**`bitbirch_recluster`** applies a two-stage reclustering pass. The first pass fits the fingerprints at `initial_threshold`, then the second pass re-clusters the leaf summaries at `second_threshold` with an optional `second_tolerance` penalty. This approach is useful when you want to first group locally similar molecules and then merge those groups at a coarser level.

```python
result = oecluster.bitbirch_recluster(
    fingerprints,
    initial_threshold=0.65,
    second_threshold=0.7,
    branching_factor=50,
    mode="strict_parity",
)
```

**`bitbirch_refine`** fits a BitBirch tree and then applies refinement passes to improve cluster quality. You can enable `redistribute_largest_cluster` to redistribute molecules from the largest cluster, or set `reassign_top_clusters` to a count (≥2) to reassign molecules from the top-K largest clusters by comparing them to cluster centroids. Refinement is helpful when the initial fit produces one or a few oversized clusters.

```python
result = oecluster.bitbirch_refine(
    fingerprints,
    threshold=0.65,
    branching_factor=50,
    redistribute_largest_cluster=True,
    reassign_top_clusters=3,
)
```

Both functions return a `BitBirchResult` with `labels`, `clusters`, `centroids`, and `cluster_sizes`. The `mode` parameter accepts `"strict_parity"` (exact reference-implementation parity) or `"fast"` (partition-merge parallelism with deterministic output). Fast mode applies to `bitbirch` and `bitbirch_recluster`; `bitbirch_refine` always runs in strict parity because its prune/reassign passes are order-sensitive (the `mode` argument is accepted for API symmetry but does not change refine's behavior).

### Choosing A Parameter

`select_parameter` sweeps one parameter over a grid you supply, scores
every partition with `cluster_report` (or `isim_report` for fingerprints),
and picks the best value under a validity index, with bounds that keep
degenerate partitions out of the running:

```python
spec = oecluster.ClusteringSpec("butina", num_threads=4)
selection = oecluster.select_parameter(
    spec, dm, "threshold", [0.2, 0.3, 0.4, 0.5, 0.6],
    criterion="silhouette", min_clusters=2, max_clusters=dm.num_samples // 2)
print(selection)                  # the scored table, winner marked
if selection.winner is not None:  # None when no threshold met the bounds
    result = selection.winner.result
```

`ClusteringSpec` names any clustering function plus its fixed options and
runs it on any input, so the same object describes the algorithm to later
workflow steps. See
[Parameter Selection](docs/python-api.md#parameter-selection) for the
criterion table, the bounds and the fingerprint path.

### Measuring Cluster Stability

`cluster_stability` reruns a spec on repeated half-size subsamples of the
items and matches every reference cluster to its best counterpart in each
resample by Jaccard overlap:

```python
spec = oecluster.ClusteringSpec("butina", threshold=0.3)
stability = oecluster.cluster_stability(spec, dm, resamples=100, seed=0)
print(stability)          # mean Jaccard, dissolved and recovered fractions
stability.mean_agreement  # adjusted Rand index, reference versus resamples
```

`take(items, indices)` is the row-subset primitive beneath it, for a
`SymmetricDistanceMatrix` or an `oefp.OEFPBatch`. See
[Cluster Stability](docs/python-api.md#cluster-stability) for the scoring
rules, the noise handling and the memory figures.

### Reaching A Consensus

`consensus` turns an ensemble of partitions into one. It counts how often
each pair of items was clustered together, divides by how often the pair was
seen together, and extracts a partition from the resulting distances:

```python
stability = oecluster.cluster_stability(spec, dm, resamples=100, seed=0)
agreed = oecluster.consensus(stability)
print(agreed)                 # clusters with their consensus scores
agreed.matrix                 # the co-association SymmetricDistanceMatrix
```

The default merges pairs supported by at least half the members that saw
them; `method=` runs any matrix-consuming `ClusteringSpec` on the matrix
instead. See
[Consensus Clustering](docs/python-api.md#consensus-clustering) for the
ensemble kinds, the scoring rules and the memory figures.

### From The Command Line

The `oecluster` command runs the same four features over a distance matrix
computed earlier, as a `.npz` from Python or as `oepdist`'s `.npy` or `.bin`
beside its JSON sidecar. The algorithm is named by its roster name, with its
options as repeated `--set key=value` — or, for a cross-algorithm `consensus`
ensemble, inside one `--member` per member:

```bash
oecluster algorithms butina
oecluster cluster distances.npz --algorithm butina --set threshold=1.5
```

```bash
oecluster select-parameter distances.npz --algorithm butina \
  --parameter threshold --values 0.5,1.0,1.5,2.0 --criterion silhouette
```

```bash
oecluster consensus distances.npz \
  --member 'butina;threshold=1.5' \
  --member 'dbscan;eps=1.5;min_samples=3' --output agreed.csv
```

Options are validated against the roster before the file is read, so a
misspelling is answered with a suggestion rather than with a `TypeError` from
inside the library. `--output` writes CSV or JSON, and `--help` documents
every command and option. See
[Command Line](docs/python-api.md#command-line) for the accepted inputs, the
output schema and the exit codes.

---

## Choosing Representatives

Representatives are cluster members chosen to summarize a cluster. They are
different from BitBirch centroid fingerprints, which are synthetic fingerprint
summaries and may not correspond to a real molecule.

| Method | What It Optimizes | Use When |
|--------|-------------------|----------|
| `medoid` | Lowest mean distance to other cluster members | You want the mathematically central real molecule. |
| `minimax` | Lowest maximum distance to any cluster member | You want the smallest cluster radius / best worst-case coverage. |
| `highest_neighborhood` | Largest fraction of neighbors within `threshold` | You want the Butina-style high-neighborhood representative. |
| `weighted_medoid` | `alpha * mean_distance + beta * liability_penalty - gamma * priority_score` | You want central molecules biased by project metadata. |

Representative quality metrics:

| Metric | Meaning |
|--------|---------|
| `mean_distance_to_cluster` | Average centrality. |
| `max_distance_to_cluster` | Worst-case coverage from this representative. |
| `median_distance_to_cluster` | Robust centrality. |
| `neighbor_fraction_at_threshold` | Fraction of cluster members within the supplied distance threshold. |
| `nearest_external_distance` | Distance to the closest molecule outside the cluster. |
| `cluster_radius` | Same as max distance from representative to cluster members. |
| `cluster_diameter` | Maximum pairwise distance within the cluster. |
| `silhouette_like_score` | Separation proxy comparing in-cluster fit to nearest external molecule. |
| `scaffold_purity` | Fraction of cluster members with the representative's scaffold label. |
| `representative_rank` | Rank under the selected scoring function. |

Weighted representative example:

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

Select more than one representative:

```python
top_by_score = oecluster.select_representatives(
    cluster,
    dm,
    k=3,
    method="medoid",
    selection="score",
)
diverse = oecluster.select_representatives(
    cluster,
    dm,
    k=3,
    method="medoid",
    selection="diversity",
)
```

---

## Assessing And Comparing Clustering Quality

`cluster_report` computes a method-agnostic scorecard for any clustering result,
and `compare_reports` aligns two scorecards side by side (e.g. Butina vs DBSCAN).

```python
butina_result = oecluster.butina(dm, threshold=0.35)
dbscan_result = oecluster.dbscan(dm, eps=0.35, min_samples=5)

butina_report = oecluster.cluster_report(butina_result, dm)
dbscan_report = oecluster.cluster_report(dbscan_result, dm)

print(butina_report)               # ClusterReport(method='butina', num_clusters=..., ...)
print(oecluster.compare_reports(butina_report, dbscan_report))
```

`cluster_report` can also score a comparison directly, without building a
distance matrix: pass the molecules with `comparison=` (or a prebuilt
comparison). The report is exact and memory stays O(N).

```python
report = oecluster.cluster_report(butina_result, mols, comparison="fingerprint")
```

Every clustering result and its report expose a read-only `.method` name
(`"butina"`, `"dbscan"`, `"hdbscan"`, `"agglomerative"`, `"k_medoids"`,
`"bitbirch"`, `"murcko"`, `"sphere_exclusion"`, `"jarvis_patrick"`, or
`"leiden"`).
`compare_reports` accepts two or more reports and labels each column by method:

```python
agglomerative_report = oecluster.cluster_report(
    oecluster.agglomerative(dm, n_clusters=10), dm)

comparison = oecluster.compare_reports(
    butina_report, dbscan_report, agglomerative_report)
print(comparison)
#   metric              butina   dbscan   agglomerative
#   num_clusters             5        7              10
#   ...
```

Distances are Tanimoto/Jaccard (`distance = 1 - similarity`). Choose a threshold
preset to match the use case:

| Preset | `coverage_thresholds` | `boundary_threshold` | Tanimoto similarity |
|--------|-----------------------|----------------------|---------------------|
| `"tight"` | `[0.20, 0.30, 0.40]` | `0.25` | cover ≥0.80/0.70/0.60; boundary ≥0.75 |
| `"default"` | `[0.25, 0.35, 0.45]` | `0.30` | cover ≥0.75/0.65/0.55; boundary ≥0.70 |
| `"diversity"` | `[0.40, 0.50, 0.60]` | `0.40` | cover ≥0.60/0.50/0.40; boundary ≥0.60 |

```python
report = oecluster.cluster_report(butina_result, dm, preset="tight")
# or override individual thresholds:
report = oecluster.cluster_report(
    butina_result, dm, coverage_thresholds=[0.2, 0.3], boundary_threshold=0.25)
# or opt in to the two stages that are off by default:
report = oecluster.cluster_report(
    butina_result, dm, compute_pair_rank_indices=True,
    compute_per_cluster_records=True)
```

All distances below are Tanimoto/Jaccard distances, so **smaller means more
similar** (distance `0.2` ≈ Tanimoto similarity `0.8`). Metrics that are
undefined for a given clustering (for example separation metrics when there is
only one cluster) are reported as `NaN`.

### Basic profile

How many clusters there are and how the molecules are distributed across them —
pure bookkeeping on the cluster labels, no distances required.

| Metric | What it tells you |
|--------|-------------------|
| `num_samples` | Total number of molecules clustered. |
| `num_clusters` | Number of clusters found (noise points are not counted as a cluster). |
| `num_noise` | Number of molecules left unassigned. HDBSCAN labels these `-1`; methods like Butina assign everything, so this is `0` for them. |
| `num_singletons` | Number of clusters containing exactly one molecule. A high count signals over-fragmentation — many compounds that didn't group with anything. |
| `noise_fraction` | `num_noise / num_samples` — the share of the library left unclustered. |
| `singleton_fraction` | Share of clusters that are singletons. With `treat_noise_as_singletons=True` (the default) noise points are counted as their own singletons here; set it `False` to base this only on real size-1 clusters. |
| `largest_cluster_fraction` | Fraction of all molecules sitting in the single biggest cluster. A value near 1.0 means one giant cluster swallowed most of the library (over-merging). |
| `cluster_size_median` | Median cluster size — the "typical" number of molecules per cluster, robust to a few very large clusters. |
| `cluster_size_p90` | 90th-percentile cluster size — most clusters are at or below this; useful for spotting a heavy tail of large clusters. |
| `size_gini` | Gini coefficient of the cluster sizes (0 = all clusters equal in size, approaching 1 = highly uneven). A quick read on whether sizes are balanced or skewed. |
| `size_entropy` | Shannon entropy (in bits) of the cluster-size distribution. Higher means sizes are more evenly spread; lower means a few clusters dominate. |

### Compactness and separation

Whether clusters are chemically tight inside and well separated from each other —
the core "are these good clusters?" view. Needs the full distance matrix.

| Metric | What it tells you |
|--------|-------------------|
| `mean_intra_distance` | Average distance between all pairs of molecules within the same cluster. Lower means members are, on average, more similar to one another. |
| `median_intra_distance` | Median of those within-cluster pairwise distances — the same idea as above but less sensitive to a few outlier pairs. |
| `median_radius` | For each cluster, the distance from its medoid (most central member) to its farthest member; reported as the median across clusters. A small radius means a typical cluster is tightly packed around its center. |
| `p95_diameter` | Cluster diameter is the largest distance between any two members of a cluster; this reports the 95th percentile across clusters. It surfaces the worst-case internal spread — how heterogeneous the loosest clusters get. |
| `silhouette` | Mean silhouette score over all clustered molecules (range −1 to 1). For each molecule it compares how close it sits to its own cluster versus the nearest other cluster; values near 1 mean tight, well-separated clusters, near 0 mean overlapping clusters, and negative means molecules may be in the wrong cluster. A singleton has no within-cluster distance to compare against, so it scores 0 rather than the 1.0 the formula would otherwise give it — Rousseeuw's convention, which scikit-learn also follows. Fragmenting a clustering therefore cannot inflate this number. |
| `dunn_index` | Smallest between-cluster distance divided by the largest cluster diameter. Higher is better: it rewards clusters that are far apart relative to how wide they are. Sensitive to outliers, so read it alongside the other metrics. |
| `boundary_violations` | Count of cross-cluster molecule pairs at or below `boundary_threshold` — i.e. pairs that look like near-neighbors yet were split into different clusters. The comparison is inclusive, so a pair sitting exactly on the threshold counts. A high count suggests the method is cutting through groups of similar compounds. |

### Representatives and coverage

If you pick one representative molecule (medoid) per cluster, how well do those
representatives stand in for the whole library? Relevant for compound selection,
purchasing, and diversity triage.

| Metric | What it tells you |
|--------|-------------------|
| `median_medoid_member_distance` | For each cluster, the average distance from its medoid to the other members; reported as the median across clusters. A small value means the chosen representative is genuinely typical of its cluster. |
| `representative_redundancy` | Median nearest-neighbor distance among the medoids themselves. A small value warns that different clusters' representatives are near-duplicates of each other (redundant chemotypes); a larger value means the representatives are diverse. |
| `coverage_thresholds` | The distance cutoffs at which coverage is evaluated (set by the preset or your override). |
| `coverage_at` | For each threshold, the fraction of **all** molecules that fall within that distance of some medoid (noise molecules are included in the denominator). It answers "if I only kept the representatives, what fraction of the library would still have a close analog?" |

`num_noise` (HDBSCAN-style unclustered points, label `-1`) and `num_singletons`
(size-1 clusters) are always reported separately. By default
`treat_noise_as_singletons=True` folds noise into the `singleton_fraction`
interpretation; set it `False` to keep them distinct. The report requires
complete pairwise distances (dense or memory-mapped storage); a sparse
(`cutoff`) matrix raises. `bitbirch_refine` can return an emptied leaf
subcluster as an empty member list, and `cluster_report` refuses a result
carrying one; nothing else the library produces trips that check.

### Internal validity indices and optional stages

Seven internal cluster-validity indices round out the scorecard. Five are always
computed: `calinski_harabasz_medoid` (higher is better),
`davies_bouldin_medoid` (lower is better),
`dunn_mean_separation_mean_diameter` and `dunn_medoid_separation_medoid_spread`
(higher is better), and `point_biserial` (higher is better). The first two are
medoid-substituted — a distance matrix has no centroids, so each cluster's
medoid stands in for one — which makes them incomparable with published or
scikit-learn figures; those two and `dunn_medoid_separation_medoid_spread`
always use the true medoid and ignore `representative_method`.
`compute_pair_rank_indices=True` adds `c_index` (lower is better) and
`baker_hubert_gamma` (higher is better). Both are read off two sorted arrays
holding every pairwise distance among the `Nc` clustered points,
`Nc * (Nc - 1) / 2` doubles between them, or roughly 400 MB at 10,000 clustered
points and 10 GB at 50,000. Both are built only under the flag; without it the
report holds no pair-sized array, which is why it is off by default.

`compute_per_cluster_records=True` fills `report.records` with one
`ClusterRecord` per cluster — size, representative, spread, nearest cluster,
silhouette and boundary violations — as a `NamedTuple` that feeds
`pandas.DataFrame(report.records)` directly. `report.noise_coverage_at` is the
coverage curve restricted to the noise points, and `report.requested` names
which optional stages the caller asked for, so a `NaN` reads unambiguously:
nobody asked, or asked and undefined. `compare_reports(...).to_table()` draws
the same line, publishing an unasked cell as `None`, printed `--`. See
[docs/python-api.md](docs/python-api.md#cluster-quality-reports) for the field
list and the memory notes.

### Agreement between two labelings

`partition_agreement` scores two labelings of the same samples against each
other, and `scaffold_agreement` scores a clustering against a per-sample
scaffold annotation. Neither needs a distance matrix.

```python
agreement = oecluster.partition_agreement(
    butina_result, dbscan_result, noise="excluded")
print(agreement.adjusted_rand_index, agreement.v_measure)
```

| Metric | Range | A low value means |
|--------|-------|-------------------|
| `adjusted_rand_index` | −0.5 to 1.0, 0.0 by chance | The two labelings put pairs together and apart no better than chance |
| `fowlkes_mallows` | 0.0 to 1.0, or `nan` | Few of the pairs grouped by one side are grouped by the other. `nan` when either side is all singletons and the partitions differ |
| `normalized_mutual_information` | 0.0 to 1.0 | Knowing one labeling tells you little about the other |
| `homogeneity` | 0.0 to 1.0, or `nan` | Side B's clusters each span many of side A's. `MI / H(a)`, so it is side A's information that side B has to explain. `nan` when side A is a single cluster and the partitions differ |
| `completeness` | 0.0 to 1.0, or `nan` | Side A's clusters each span many of side B's. `MI / H(b)`, the transpose. `nan` when side B is a single cluster and the partitions differ |
| `v_measure` | 0.0 to 1.0 | Same number as `normalized_mutual_information`, reported under both names |
| `adjusted_mutual_information` | Below 0.0 to 1.0, 0.0 by chance, or `nan` | As NMI, but corrected for the agreement many small clusters produce by chance. Opt in with `adjusted_mutual_information=True`; `nan` until you do |

Every metric is `nan` when fewer than two samples survive noise handling. See
[docs/python-api.md](docs/python-api.md#partition-agreement) for the `noise=`
readings and the divergences from scikit-learn.

### Structure-activity coherence

`sar_coherence` decomposes an activity vector across a labeling.
`activity_landscape` and `modelability` score the structure-activity
relationship straight from a distance matrix or a comparison
(`activity_landscape(mols, activity, comparison="fingerprint")`), with no
clustering in between.

```python
coherence = oecluster.sar_coherence(butina_result, activity)
landscape = oecluster.activity_landscape(dm, activity)
print(coherence.omega_squared, landscape.cliff_density, landscape.rmodi)
print(oecluster.modelability(dm, classes).modi)
```

| Metric | Range | A low value means |
|--------|-------|-------------------|
| `eta_squared` | 0.0 to 1.0, or `nan` | Activity varies as much inside the clusters as between them. Rises with the cluster count on its own, so compare it only across labelings of the same granularity |
| `omega_squared` | Below 0.0 to 1.0, about 0.0 by chance, or `nan` | The labeling explains no more activity variance than a random one of the same shape. Negative means less than chance. Read that chance level as approximately rather than exactly zero: it lifts when the clusters are mostly singletons and the activity variance is concentrated in a few samples -- conditions a tight `butina()` threshold on screening data meets together. One active among 99 inactives under one 11-member cluster plus 89 singletons has an exact chance expectation of 0.075, and five actives among 95 have 0.014, where balanced clusterings of two to fifty clusters stay below 0.0002 on those same activities. Distinct values are no protection: 100 distinct activities spanning 1e-6 to 1.0 reach 0.075 again |
| `cliff_density` | 0.0 to 1.0, or `nan` | Few of the scored pairs are both near and sharply different in activity. The denominator is every scored pair, not the near ones, so a low value is no evidence of a smooth neighbourhood: 100 samples whose only near pair is a cliff read 0.000202 |
| `max_sali`, `mean_sali` | 0.0 upward, or `nan` | For `max_sali`, no pair is both structurally close and far apart in activity. For `mean_sali`, only that the typical scored pair is not: it averages over every pair with a defined ratio, so a lone cliff contributes just its own ratio divided by that count and a low mean can sit beside a high `max_sali`. Zero-distance pairs are excluded from both and counted in `num_zero_distance_pairs` |
| `rmodi` | 0.0 to 1.0, or `nan` | A sample's nearest neighbour is usually outside its activity band -- which reaches `rmodi_delta * activity_stddev` either side. Either end of the range is reachable from the band width alone, so read the figure against the reported `activity_stddev` and sweep `rmodi_delta` before concluding anything from it |
| `modi` | 0.0 to 1.0, or `nan` | Some class's members usually have a nearest scored neighbour of another class. It is the mean over classes of each class's same-class fraction, not a per-sample rate, so one small class holds it down however well the rest separate: 100 evenly spaced class-A points with one class-B point half a step past the end read 0.495 -- about what random labels at those class counts would score -- while 98% of samples do have a same-class nearest neighbour. `report.classes` carries the per-class fractions |

A `nan` activity is missing data rather than a value, as is an empty class
string, and every result reports `num_scored` beside `num_samples` so the
difference is visible. Only for `activity_landscape` and `modelability` is that
difference entirely missing data: `sar_coherence` counts what `noise=` drops in
it as well, so under its `"excluded"` default
`sar_coherence([-1, 0, 0, 1, 1], [1.0, 2.0, 3.0, 4.0, 5.0])` reports
`num_scored` 4 of 5 with every activity usable. See
[docs/python-api.md](docs/python-api.md#sar-coherence) for the per-cluster and
per-class tables, the `noise=` readings, and the one distance-matrix refusal
these two entry points will not waive.

`isim_report()` is an approximate scorecard for a fingerprint clustering,
built directly from the fingerprint batch rather than an N x N distance
matrix. It is built on `isim()`, the set-similarity primitive underneath it:
the union-weighted Tanimoto similarity of a set of binary fingerprints. The
default core is linear in the number of fingerprints; an optional
`compute_centroid_indices=True` stage adds silhouette, nearest cluster,
medoid Davies-Bouldin/Dunn and coverage at O(N K + K^2) cost.

```python
result = oecluster.bitbirch(fps, threshold=0.65)
report = oecluster.isim_report(result, fps)
print(report.isim_intra_distance)
```

See [docs/python-api.md](docs/python-api.md#approximate-reports-from-fingerprints)
for which fields are exact, which are iSIM ratios, and when to reach for
`cluster_report()` instead.

---

## Scaling Guidance

Dense pairwise storage uses scipy-compatible condensed indexing and stores
`N * (N - 1) / 2` doubles. Memory is approximately:

| Molecules | Dense Distance Memory |
|-----------|-----------------------|
| 10,000    | 0.4 GB                |
| 50,000    | 10 GB                 |
| 100,000   | 40 GB                 |

Practical guidance:

- Use `DenseStorage` for small to medium datasets and full representative
  metrics.
- Use `MMapStorage` through `pdist(..., output="distances.mmap")` when dense
  distances fit on disk but not comfortably in RAM.
- Use `SparseStorage` through `pdist(..., cutoff=...)` for threshold graph
  algorithms such as Butina/DBSCAN when you do not need complete distances.
- Use complete dense or memory-mapped distances for representative ranking,
  HDBSCAN, and agglomerative workflows that need all pairwise distances.
- Use BitBirch for dense binary OEFP batches when pairwise precomputation is
  not the right first step.

`threshold` means an algorithm-level neighbor or merge distance. `cutoff`
means a storage-level distance filter used while computing pairwise distances.
If `cutoff` is lower than an algorithm threshold, the algorithm cannot see all
required neighbors.

---

## Comparison Methods

### Fingerprint

| Parameter | Values | Default |
|-----------|--------|---------|
| `fp_type` | `morgan`, `atom_pair`, `topological_atom_pair`, `topological_torsions` | `morgan` |
| `storage` | `binary`, `count`, `sparse`, `sparse_count` | `binary` |
| `metric` | 17 OEFP scalar metrics, from `jaccard` to `tversky` | `tanimoto` |
| `numbits` | Fingerprint size in bits | `2048` |
| `radius` | Morgan radius | `2` |
| `min_distance` | Minimum atom-pair graph distance | `1` |
| `max_distance` | Maximum atom-pair graph distance | `30` |
| `torsion_atom_count` | Torsion path length | `4` |
| `use_chirality` | Distinguish stereocenters | `False` |
| `p` | Minkowski order | `2.0` |
| `tversky_alpha`, `tversky_beta` | Tversky weights | `0.5` |

Similarity mode is supported for `tanimoto` and required by `tversky`, which
has no distance form and refuses the default `similarity=False`. A counted
storage needs a metric defined on counts, such as `bray_curtis` or `manhattan`,
and naming an option the rest of the configuration would ignore raises rather
than accepting a value that has no effect. **Changed in 5.0.0:** `max_distance`
no longer sets the Morgan radius; use `radius` for that. See
[docs/python-api.md](docs/python-api.md#fingerprint) for the full metric list
and the migration notes.

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

Accepts molecules directly or design units. Design units are converted to their
protein components automatically.

### SiteHopper

Use `method="sitehopper"` in the superpose comparison. It accepts design units,
generates patch surfaces when needed, and compares binding-site patch scores.

### Descriptor

| Parameter | Values | Default |
|-----------|--------|---------|
| `sources` | `openeye`, `mordred`, `rdkit` | `["openeye"]` |
| `columns` | Explicit column names | every numeric column of the selected sources |
| `metric` | 10 metrics, from `euclidean` to `mahalanobis` | `standardized_euclidean` |
| `missing` | `complete_case`, `propagate`, `ignore` | `complete_case` |

Computes distances over molecular descriptor columns. The default metric
standardizes each column by the variance fitted over the molecules passed in,
so the same pair gets a different distance in a different set; pass `sources=`,
`columns=` and `variances=` from `descriptor_statistics()` together to fix the
scaling across runs, `sources=` included because column names are resolved
against it and an `rdkit` or `mordred` name is not in the default `openeye`
schema.
See [docs/python-api.md](docs/python-api.md#descriptor) for the full parameter
table, the column-selection rules, and what each `missing` policy stamps on the
matrix.

### RMSD

| Parameter | Values | Default |
|-----------|--------|---------|
| `overlay` | Superpose before measuring | `False` |
| `automorph` | Minimize over graph automorphisms | `True` |
| `heavy_only` | Ignore hydrogens | `True` |
| `expand_conformers` | Expand multi-conformer inputs into one item per pose | `True` |

Compares poses of one molecule: every input must share a topology and carry
coordinates. Use `rocs` to compare different molecules by shape. See
[docs/python-api.md](docs/python-api.md#rmsd) for the conformer labeling and
the `expand_conformers` caveat.

### MCS

| Parameter | Values | Default |
|-----------|--------|---------|
| `search_mode` | `approximate`, `exhaustive` | `approximate` |
| `match_level` | `default`, `exact`, `loose` | `default` |
| `max_matches` | Matches one directed search may enumerate | `1024` |

Maximum common substructure scored as Tanimoto over matched bonds. Topological,
so coordinates are never read and hydrogens are suppressed wherever they can be
folded into a heavy atom; molecules with no bonds after suppression are refused.
`exhaustive` is one to three orders of magnitude slower than `approximate` and
is not reliably better. There is no metric guarantee, so the matrix reports
`triangle` as `unknown`. See [docs/python-api.md](docs/python-api.md#mcs) for
the cost estimate and the measurements behind both caveats.

---

## Storage Backends

| Backend | Use Case | Memory | Thread-Safe |
|---------|----------|--------|-------------|
| `DenseStorage` | Default full pairwise matrix | O(N^2) in memory | Yes |
| `MMapStorage` | Full pairwise matrix backed by a file | O(N^2) on disk | Yes |
| `SparseStorage` | Cutoff-filtered results | Stored entries only | Yes |

All backends use scipy-compatible condensed distance matrix indexing:

```text
n * i - i * (i + 1) / 2 + j - i - 1
```

Reached directly, `Get` and `Set` refuse a pair they cannot address: an index
at or beyond `NumSamples()` on either, and the diagonal on `Set`, which owns no
stored slot. Python sees `RuntimeError`.

---

## Command-Line Tool

`oepdist` computes pairwise and cross-distance matrices from molecular files.

```bash
oepdist fp molecules.sdf -o distances.npy \
  --fp-type morgan --metric tanimoto --threads 8
```

```bash
oepdist rocs conformers.oeb -o shape_dist.npy --score combo_norm
```

```bash
oepdist superpose structures.oeb -o rmsd.npy --method global_carbon_alpha
```

Cross-distance uses two input files:

```bash
oepdist fp queries.sdf targets.sdf -o cross_dist.npy
```

Output format is determined by extension:

| Extension | Format |
|-----------|--------|
| `.npy` | NumPy array plus JSON sidecar |
| `.csv` | Labeled comma-separated values |
| `.bin` | Raw double array plus JSON sidecar |

`oecluster` clusters what `oepdist` wrote. It reads the `.npy` and `.bin`
output through the sidecar and refuses the `.csv`, which records neither
provenance nor full precision; see
[From The Command Line](#from-the-command-line).

---

## C++ API

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

---

## Troubleshooting

| Symptom | Likely Cause | Fix |
|---------|--------------|-----|
| `ImportError` for `_oecluster` | Python extension was not built or cannot find runtime libraries | Rebuild with `scripts/build_python.py` or ensure the wheel matches your Python and platform. |
| OpenEye import or license failure | OpenEye Toolkits or license is missing at runtime | Install OpenEye Toolkits and configure your OpenEye license before importing or running examples. |
| `Unknown comparison` | The comparison string is misspelled | Use one of `descriptor`, `fingerprint`, `mcs`, `rmsd`, `rocs`, `sitehopper`, or `superpose`; the message itself lists the valid names. |
| `Unknown representative method` | Representative method name is misspelled | Use `medoid`, `minimax`, `highest_neighborhood`, or `weighted_medoid`. |
| `highest_neighborhood representative requires a threshold` | The method needs a neighbor cutoff | Pass `threshold=<distance>`. |
| Sparse storage cutoff error | `cutoff` is lower than the clustering threshold | Recompute distances with a cutoff at least as large as the clustering threshold, or use dense/mmap storage. |
| `SparseStorage cannot provide complete distances` | The workflow needs every pairwise distance | Use dense or memory-mapped storage for representative ranking, HDBSCAN, or agglomerative clustering. |
| ROCS or superpose warnings about coordinates | Input molecules lack required 3D coordinates | Generate conformers or use a 2D/fingerprint workflow instead. |
| Process runs out of memory | Dense pairwise matrix is too large | Use `output=...` for mmap storage, `cutoff=...` for sparse threshold workflows, or BitBirch for OEFP batches. |

---

## Examples

Runnable examples live in `examples/`:

| Example | What It Shows |
|---------|---------------|
| `quickstart_smiles.py` | End-to-end clustering from inline SMILES. |
| `rank_representatives.py` | Weighted ranking plus score/diversity k-representative selection. |
| `getting-started.ipynb` | Notebook walkthrough over a real dataset: distances, four of the clustering algorithms, and a side-by-side quality comparison. |

Run the scripts with the local package on your `PYTHONPATH`:

```bash
python examples/quickstart_smiles.py
python examples/rank_representatives.py
```

The notebook additionally imports `oepandas` and `cnotebook`, which are not
OECluster dependencies; it reads its input from `examples/assets/Bender.csv`
and must be run from the `examples/` directory.

## License

MIT License. See [LICENSE](LICENSE) for details.
