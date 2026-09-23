# Quickstart

This guide shows the usual OECluster workflow: load molecules, compute a
pairwise distance matrix, cluster the molecules, choose representatives, and
review cluster quality.

## Cluster Molecules From SMILES

This example needs no input files. It builds molecules from SMILES, computes
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

dm = oecluster.pdist(mols, "fingerprint", fp_type="morgan", metric="tanimoto")
result = oecluster.butina(dm, threshold=0.55)

for cluster_id, cluster in enumerate(result.clusters):
    ranked = oecluster.rank_representatives(cluster, dm, method="medoid")
    medoid = ranked[0]
    print(cluster_id, mols[medoid.member].GetTitle(), len(cluster),
          medoid.metrics.cluster_radius)
```

Run the same workflow from the tested example file:

```bash
python examples/quickstart_smiles.py
```

## Loading Molecules From Files

OECluster works with any OpenEye molecule objects. Read them with the OpenEye
toolkit and pass the list to `pdist()`:

```python
from openeye import oechem

mols = []
ifs = oechem.oemolistream("molecules.sdf")
mol = oechem.OEGraphMol()
while oechem.OEReadMolecule(ifs, mol):
    mols.append(oechem.OEGraphMol(mol))
```

## Computing A Distance Matrix

`pdist()` computes the symmetric within-set pairwise distances. Choose the
storage backend by keyword: dense by default, sparse with `cutoff`, or
memory-mapped with `output`.

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

# Sparse storage for threshold graph algorithms:
sparse_dm = oecluster.pdist(mols, "fingerprint", metric="tanimoto", cutoff=0.35)

# Memory-mapped storage for large dense matrices:
mmap_dm = oecluster.pdist(mols, "fingerprint", output="distances.mmap")
```

Use `cdist()` for the rectangular cross distance matrix between two distinct
item sets, for example screening hits against an in-house collection:

```python
cross_dm = oecluster.cdist(queries, targets, "fingerprint", metric="tanimoto")
print(cross_dm.matrix[0, :])  # distances from first query to all targets
```

## Clustering And Representatives

Most clustering functions take a `DistanceMatrix`. Two do not: BitBirch takes an
OEFP fingerprint batch, and Murcko takes molecules directly. All of them return
a result that subclasses `ClusteringResult`, which exposes a scikit-learn-style
`labels` array and a tuple of `clusters`.

```python
result = oecluster.butina(dm, threshold=0.35, reordering=False)

# Drop the assignment array straight into a dataframe column.
labels = result.labels

for cluster in result.clusters:
    medoid = oecluster.representative(cluster, dm, method="medoid")
    minimax = oecluster.representative(cluster, dm, method="minimax")
    ranked = oecluster.rank_representatives(cluster, dm, method="medoid")
    selected = oecluster.select_representatives(
        cluster, dm, k=2, method="medoid", selection="diversity",
    )
    print(medoid, minimax, ranked[0].metrics.cluster_radius, selected)
```

The other algorithms follow the same shape:

```python
db = oecluster.dbscan(dm, eps=0.35, min_samples=5)
hdb = oecluster.hdbscan(dm, min_cluster_size=5)
agg = oecluster.agglomerative(dm, n_clusters=10)
km = oecluster.k_medoids(dm, n_clusters=10)
```

BitBirch clusters OEFP binary fingerprint batches directly, without
materializing a pairwise distance matrix first:

```python
result = oecluster.bitbirch(fingerprints, threshold=0.65, branching_factor=50)
print(result.centroids, result.cluster_sizes)
```

Murcko clusters molecules by scaffold identity, with no distance matrix at all:

```python
scaffolded = oecluster.murcko(mols)
print(scaffolded.cluster_scaffolds[0], scaffolded.labels[:5])
```

See [Choosing representatives](python-api.md#representatives) for the weighted
medoid and k-representative selection options shown in
`examples/rank_representatives.py`.

## Reviewing Cluster Quality

`cluster_report()` computes a scorecard for any clustering result, and
`compare_reports()` aligns several scorecards side by side:

```python
butina_report = oecluster.cluster_report(result, dm)
dbscan_report = oecluster.cluster_report(db, dm)

print(butina_report)
print(oecluster.compare_reports(butina_report, dbscan_report))
```

Distances are Tanimoto/Jaccard, so smaller means more similar. Choose a preset
(`"tight"`, `"default"`, or `"diversity"`) or override individual thresholds.
The report needs complete pairwise distances, so dense or memory-mapped storage
is required; a sparse (`cutoff`) matrix raises.

Alongside the original scorecard the report carries seven internal validity
indices. Five are always computed: `calinski_harabasz_medoid` (higher is
better), `davies_bouldin_medoid` (lower is better),
`dunn_mean_separation_mean_diameter` and `dunn_medoid_separation_medoid_spread`
(higher is better), and `point_biserial` (higher is better).
`dunn_medoid_separation_medoid_spread` ignores `representative_method` and
always uses the true medoid, so a report requested with
`representative_method="minimax"` still reports medoid-based values here.

> Two of these five are **medoid-substituted** indices. The published
> Calinski-Harabasz and Davies-Bouldin definitions use centroids, which do not
> exist for a distance matrix, so each cluster's medoid stands in for its
> centroid. Calinski-Harabasz also has a grand-mean term, and the global medoid
> stands in for that; Davies-Bouldin has no such term. The values are therefore
> not comparable with published figures or with scikit-learn's. Both ignore
> `representative_method` and always use the true medoid, so a report requested
> with `representative_method="minimax"` still reports medoid-based values here.

Two options, both off by default, turn on the optional stages:

```python
report = oecluster.cluster_report(
    result,
    dm,
    compute_pair_rank_indices=True,    # adds c_index and baker_hubert_gamma
    compute_per_cluster_records=True,  # populates report.records
)
```

`compute_pair_rank_indices` adds `c_index` (lower is better) and
`baker_hubert_gamma` (higher is better). Both are read off two sorted arrays
holding every pairwise distance among the `Nc` clustered points,
`Nc * (Nc - 1) / 2` doubles between them -- roughly 400 MB at `Nc = 10,000` and
10 GB at `Nc = 50,000`. Only the between-cluster array is the flag's own cost,
since the within-cluster one is built on every call, so the flag adds nothing to
a single-cluster result and nearly the whole figure to one with small clusters.
It is off by default for the second case. `compute_per_cluster_records` fills `report.records` with one
`ClusterRecord` per cluster, a `NamedTuple` that feeds
`pandas.DataFrame(report.records)` directly. `report.noise_coverage_at` is the
coverage curve restricted to the noise points, parallel to `coverage_at`.

`report.requested` names which of the two optional computations the caller asked
for, so a `NaN` can be read unambiguously: false means nobody asked, true with
`NaN` means asked and undefined. `compare_reports(...).to_table()` follows the
same rule -- a cell is `None` when that report never asked, including a coverage
threshold the report carries but answered nothing at, and `nan` only when the
answer is undefined. `__repr__` renders `None` as `--`.

## Exporting Cluster Assignments

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
            writer.writerow([cluster_id, member, mols[member].GetTitle(), ranks[member]])
```

## Command-Line Use

The `oepdist` command computes pairwise and cross distance matrices from
molecular files without writing a Python script:

```bash
oepdist fp molecules.sdf -o distances.npy --fp-type morgan --metric tanimoto
```

See [CLI](cli.md) for the full command surface.
