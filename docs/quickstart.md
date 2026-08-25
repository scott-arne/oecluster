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

All clustering functions take a `DistanceMatrix` (except BitBirch, which takes
an OEFP fingerprint batch) and return a result that subclasses
`ClusteringResult`. Results expose a scikit-learn-style `labels` array and a
tuple of `clusters`.

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
```

BitBirch clusters OEFP binary fingerprint batches directly, without
materializing a pairwise distance matrix first:

```python
result = oecluster.bitbirch(fingerprints, threshold=0.65, branching_factor=50)
print(result.centroids, result.cluster_sizes)
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
