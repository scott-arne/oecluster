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
- [Choosing Representatives](#choosing-representatives)
- [Scaling Guidance](#scaling-guidance)
- [Comparison Methods](#comparison-methods)
- [Storage Backends](#storage-backends)
- [Command-Line Tool](#command-line-tool)
- [C++ API](#c-api)
- [Troubleshooting](#troubleshooting)
- [Examples](#examples)
- [Project Structure](#project-structure)
- [License](#license)

---

## What You Can Do

- **Cluster molecular collections** with Butina, DBSCAN, HDBSCAN,
  agglomerative clustering, or BitBirch.
- **Choose representatives** with true medoids, minimax/radius centers,
  highest-neighborhood Butina-style representatives, weighted medoids, ranked
  representative lists, and k-representative selection.
- **Compute molecular distances** for fingerprints, ROCS shape/color overlay,
  protein superposition, and binding-site comparison.
- **Scale distance storage** with dense in-memory arrays, memory-mapped files,
  or sparse cutoff-filtered storage.
- **Use the same core from Python, C++, or CLI workflows**.

---

## Requirements

- **Python** 3.10+ with NumPy.
- **OpenEye Toolkits** 2025.2 or later.
- **A valid OpenEye license** at build time and runtime.
- **OEFP** 0.2.4 or later for fingerprint generation and comparison.
- **C++17**, **CMake** 3.16+, and **SWIG** 4.0+ when building from source.

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

| Algorithm | Input | Best For | Key Parameters | Notes |
|-----------|-------|----------|----------------|-------|
| `butina` | `DistanceMatrix` | Classic fingerprint clustering and threshold-neighborhood summaries | `threshold`, `reordering` | Returns a `ButinaResult` (subclass of `ClusteringResult`); the first member of each cluster is the highest-neighborhood representative chosen by the Butina pass. |
| `dbscan` | `DistanceMatrix` | Density clustering with explicit noise points | `eps`, `min_samples` | Returns a `DBSCANResult` (subclass of `ClusteringResult`) with `core_sample_indices`. |
| `hdbscan` | `DistanceMatrix` | Variable-density clusters and soft membership | `min_cluster_size`, `min_samples`, `cluster_selection_method` | Requires complete precomputed distances. Returns an `HDBSCANResult` (subclass of `ClusteringResult`) with `probabilities`. |
| `agglomerative` | `DistanceMatrix` | Hierarchical clustering and dendrogram-style analysis | `n_clusters`, `distance_threshold`, `linkage` | Supports single, complete, average, and weighted linkage. Returns an `AgglomerativeResult` (subclass of `ClusteringResult`) with `children`, `distances`, and `cluster_sizes`. |
| `bitbirch` | `oefp.OEFPBatch` | Large dense binary fingerprint sets where feature-space clustering is preferred | `threshold`, `branching_factor`, `merge_criterion` | `mode="fast"` is currently a reserved strict-parity alias. Returns a `BitBirchResult` (subclass of `ClusteringResult`) with `centroids` and `cluster_sizes`. |

Start with **Butina** for familiar fingerprint threshold clustering. Use
**DBSCAN/HDBSCAN** when noise and density matter. Use **agglomerative** when
you need linkage structure. Use **BitBirch** when the workflow already has
OEFP dense binary fingerprints and you want feature-space clustering without
materializing an initial pairwise distance matrix.

All clustering functions return a result that subclasses `ClusteringResult`,
which exposes read-only `labels` (a length-n, scikit-learn-style assignment
array) and `clusters` (a tuple of member-index tuples). Results support
`len(result)` (number of clusters), iteration over clusters, and `result[i]`
indexing. Algorithm-specific outputs live on the specific subclass (for
example `DBSCANResult.core_sample_indices` or `BitBirchResult.centroids`), so a
result never carries fields that do not apply to its algorithm.

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

print(butina_report)               # ClusterReport(num_clusters=..., ...)
print(oecluster.compare_reports(butina_report, dbscan_report))
```

Every clustering result and its report expose a read-only `.method` name
(`"butina"`, `"dbscan"`, `"hdbscan"`, `"agglomerative"`, or `"bitbirch"`).
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
```

Reported metrics:

| Family | Metrics |
|--------|---------|
| Basic profile | `num_samples`, `num_clusters`, `num_noise`, `num_singletons`, `noise_fraction`, `singleton_fraction`, `largest_cluster_fraction`, `cluster_size_median`, `cluster_size_p90`, `size_gini`, `size_entropy` |
| Compactness / separation | `mean_intra_distance`, `median_intra_distance`, `median_radius`, `p95_diameter`, `silhouette` (true mean per-point), `dunn_index`, `boundary_violations` |
| Representative / coverage | `median_medoid_member_distance`, `representative_redundancy`, `coverage_thresholds`, `coverage_at` |

`num_noise` (HDBSCAN-style unclustered points, label `-1`) and `num_singletons`
(size-1 clusters) are always reported separately. By default
`treat_noise_as_singletons=True` folds noise into the `singleton_fraction`
interpretation; set it `False` to keep them distinct. The report requires
complete pairwise distances (dense or memory-mapped storage); a sparse
(`cutoff`) matrix raises.

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
| `score_type` | `ComboNorm`, `Combo`, `Shape`, `Color` | `ComboNorm` |
| `color_ff_type` | `1` (ImplicitMillsDean), `2` (ExplicitMillsDean) | `1` |

Distance ranges: ComboNorm [0, 1], Combo [0, 2], Shape [0, 1], Color [0, 1].

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
| `Unknown comparison` | The comparison string is misspelled | Use `fingerprint`, `rocs`, `superpose`, or `sitehopper` where supported. |
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

Run them with the local package on your `PYTHONPATH`:

```bash
PYTHONPATH=python python examples/quickstart_smiles.py
PYTHONPATH=python python examples/rank_representatives.py
```

---

## Project Structure

```text
oecluster/
  include/oecluster/         C++ public headers
    comparisons/             Comparison-specific headers
    clustering/              Clustering and representative APIs
  src/                        C++ implementation
  tools/                      CLI tool
  swig/                       SWIG interface for Python bindings
  python/oecluster/           Python package
  examples/                   Runnable examples
  tests/
    cpp/                      C++ tests
    python/                   Python tests
  scripts/                    Build and packaging utilities
  CMakeLists.txt              Root build configuration
  pyproject.toml              Python package configuration
```

---

## License

MIT License. See [LICENSE](LICENSE) for details.
