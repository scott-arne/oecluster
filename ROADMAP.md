# OECluster Clustering Expansion Roadmap

The clustering expansion is delivered as a sequence of sub-projects, each with
its own design spec, implementation plan, and release. This file records the
decomposition and the order.

It exists because the original decomposition was never written to a tracked
file. It lived in conversation and in a single order string copy-pasted into
four design specs under `docs/hyperpowers/specs/`, which is gitignored. By the
time sub-project A closed, three of the eight letters could no longer be
resolved to any content. Sections marked *re-derived* below were reconstructed
on 2026-09-18 from the surviving cross-references and the state of the library;
they are not recoveries of the original intent.

## Order

    F -> A -> D1 -> C -> E -> D2a -> D2b -> D2c -> B -> D3/D4

| Piece | Content | Status |
| --- | --- | --- |
| F | Comparison layer expansion | shipped 5.0.0 |
| A1 | Internal indices, per-cluster record table, `noise_coverage_at` | shipped 5.1.0, 5.2.0 |
| A2 | Partition agreement: ARI/AMI/NMI/V-measure/Fowlkes-Mallows, scaffold-ARI | shipped 5.3.0 |
| A3 | SAR coherence: eta-squared, omega-squared, SALI, cliff density, MODI/RMODI | shipped 5.4.0 |
| D1 | k-medoids/PAM: exactly `k` clusters with real-member centers | shipped 5.5.0 |
| C1 | Murcko scaffold assignment and scaffold-identity clustering | shipped 5.6.0 |
| C2 | MCS-based comparison: Tanimoto over matched bonds | shipped 5.7.0 |
| E1 | Diversity selection: `maxmin_select` and #Circles | shipped 5.8.0 |
| E2 | Set diversity scores: Vendi score and log-determinant diversity | shipped 5.9.0 |
| D2a | Sphere exclusion: leader, Butina and DISE as one engine | shipped 5.10.0 |
| D2b | k-nearest-neighbor graph and Jarvis-Patrick clustering | shipped 5.11.0 |
| D2c | Leiden community detection over the D2b graph | shipped 5.12.0 |
| B1 | Exact O(N)-memory paths for the A metrics over a comparison | shipped 5.13.0 |
| B2 | Approximate fingerprint-native (iSIM) A metrics | shipped 5.14.0 |
| D3a | Parameter selection: `ClusteringSpec`, `select_parameter` | shipped 5.15.0 |
| D3b | Stability resampling: subsample-Jaccard over a `ClusteringSpec` | shipped 5.16.0 |
| D3c | Consensus clustering over an ensemble of partitions | shipped 5.17.0 |
| D3d | Clustering CLI over oepdist output | shipped 5.18.0 |
| D4a | Streaming threshold-graph clustering: Butina, DBSCAN, neighbor-order sphere exclusion from a comparison | shipped 5.19.0 |
| D4b | Streaming MST clustering: HDBSCAN, single-linkage agglomerative | shipped 5.20.0 |
| D4c | Algorithms that cannot stream exactly: k-medoids, complete/average/weighted linkage | planned |

Sub-project A was originally scoped as one piece covering seven metric families
across five input shapes. It was decomposed into A1, A2 and A3 during
execution, following the rule that the pieces follow the input shapes. That
decomposition is the reason the remaining sub-projects are kept to a size that
fits one spec.

## D1 - k-medoids

*Re-derived.* Adds k-medoids/PAM: exactly `k` clusters whose centers are real
members chosen by minimizing total within-cluster distance.

Nothing in the library produces a fixed cluster count except
`agglomerative(n_clusters=...)`, and nothing chooses centers by optimizing a
criterion. Medoids are real members rather than synthetic summaries, so they
compose with the existing `Representative` machinery instead of duplicating it.

One of the supported initialization strategies is deterministic farthest-first
selection, better known as MaxMin. That kernel is the standard primitive for
diverse-subset selection, so E reuses it for a public library-scale
`maxmin_select()`. Having a consumer inside D1 is what keeps it from being
speculative.

Two algorithms were considered and excluded:

- **OPTICS** belongs to the density family alongside DBSCAN and HDBSCAN rather
  than the exemplar family, and HDBSCAN already answers most of what it adds.
- **Sphere exclusion** is already in the library under another name.
  `butina_cluster()` seeds by descending neighbor count, claims every unseen
  neighbor within the threshold, and repeats - the leader/sphere-exclusion
  loop, with `ButinaOptions::reordering` as the Taylor-Butina refinement. The
  family is assigned to D2 as "leader/DISE" (DISE being Directed Sphere
  Exclusion), where it can generalize the shipped algorithm rather than ship a
  near-duplicate beside it.

## C - Chemistry-native clustering

*Re-derived.* Partitions on chemical structure rather than on distance: Murcko
scaffold assignment and MCS-based comparison.

Before C, every algorithm in the library was distance-driven. `scaffold_labels`
is consumed as an input in three places - scaffold-ARI in A2, the
weighted-medoid representative, and scaffold purity in the report - but nothing
in the library produced them. C closes that asymmetry and gives E a chemically
meaningful axis to diversify along.

C1 shipped in 5.6.0 as `murcko_scaffolds` and `murcko`, closing the producer
side of that asymmetry. C2 shipped in 5.7.0 as the `mcs` comparison. It is
scoped to a comparison method rather than a clustering entry point: a dedicated
MCS clustering algorithm and common-core reporting were both offered and
declined in favour of feeding the existing algorithms, either of which
remains available as an additive follow-on.

## E - Diversity metrics and library-scale selection

*Pinned by the A1 and A3 specs.* Adds #Circles, Vendi score and log-determinant
diversity, plus a selection module that operates over a whole library rather
than within a single cluster. The existing `select_representatives()` selects
within one cluster; E selects across a collection.

These metrics score a set rather than a partition, which is why they were moved
out of sub-project A.

E was split in two. E1 shipped in 5.8.0 as `maxmin_select`, deterministic
farthest-first selection, and `circles`, the #Circles coverage measure; both
run on one kernel and read a precomputed matrix or evaluate a comparison
lazily, so a library too large for an O(N^2) matrix can still be subset. E2
shipped in 5.9.0 as `vendi_score` and `logdet_diversity`. Both read a
set's similarity kernel. The exact scores decompose it with an in-tree
symmetric eigenvalue solver (Householder tridiagonalization and implicit QL),
chosen over Eigen, LAPACK through numpy, and a Jacobi solver. It adds no
dependency behind the build firewall, keeps the exact scores in the C++ core,
and is several times faster than Jacobi. Order-2 Vendi needs no spectrum and
runs lazily at any size.

## D2 - Graph and leader algorithms

*Pinned by the A1 and A2 specs.* Adds Leiden, Jarvis-Patrick and leader/DISE to
the algorithm roster, in three slices:

- **D2a** (shipped 5.10.0): `sphere_exclusion()`, one engine with input
  (leader), neighbor-count (Butina) and caller-permutation (DISE) seed orders
  and first-claim or nearest assignment. `butina_cluster()` became an adapter
  over it with unchanged outputs.
- **D2b** (shipped 5.11.0): `knn_graph()` and `KNNGraph`, a public k-nearest-neighbor graph over dense, memory-mapped, sparse and lazy input, and `jarvis_patrick()` with the classic mutual shared-neighbor rule, over a graph or raw input.
- **D2c** (shipped 5.12.0): `leiden()`, Leiden community detection with modularity and CPM objectives over a shared-nearest-neighbor Jaccard weighting of the D2b graph, seeded and with connected clusters guaranteed.

Cluster stability by resampling shipped in 5.16.0 as D3b (subsampling without
replacement; see the D3 section below), once the D2 roster was final.

## B - Fingerprint-native metrics

*Pinned by the A1 and A3 specs.* Adds counterparts of the metrics defined in A
that avoid the materialized O(N^2) distance matrix the current forms require.
A owns the definitions; B owns the scale path.

- **B1 (shipped 5.13.0).** `cluster_report`, `activity_landscape` and
  `modelability` score a comparison directly. Results are exact -- identical
  to the matrix forms over a matrix filled through `Compare` -- with O(N)
  memory and O(N^2) time.
- **B2 (shipped 5.14.0).** `isim()` and `isim_report()` score binary
  fingerprint batches directly in O(N) time, trading exactness for scale:
  `isim()` is the union-weighted Tanimoto similarity of a fingerprint set,
  and `isim_report()` is a `cluster_report`-shaped scorecard built from it,
  with an opt-in O(N K) centroid stage for silhouette, nearest cluster,
  medoid Davies-Bouldin/Dunn and coverage. `activity_landscape` and
  `modelability` stay exact-only; B2 does not add approximate forms of
  either.

## D3 - Workflow layer

*Re-derived.* Four gaps that share a consumer rather than an implementation,
shipped as four slices:

- **D3a (shipped 5.15.0).** Parameter selection. `ClusteringSpec` describes
  an algorithm plus its fixed options and runs it on any input;
  `select_parameter` sweeps one keyword over an explicit grid, scores every
  partition with `cluster_report` (a matrix or a prebuilt comparison) or
  `isim_report` (fingerprints), and picks the winner under a named validity
  index with optional noise and cluster-count bounds. Pure Python over the
  existing entry points: the first slice with no native core.
- **D3b (shipped 5.16.0).** Stability resampling. `cluster_stability` reruns
  a `ClusteringSpec` on subsamples drawn without replacement and reports
  Hennig's per-cluster Jaccard statistics (mean, dissolved, recovered) plus
  one adjusted Rand index per resample; `take` is the row-subset primitive
  over `SymmetricDistanceMatrix` and `oefp.OEFPBatch`, backed by the native
  `take_pairs` and `take_fingerprints`, which D3c and D3d reuse.
- **D3c (shipped 5.17.0).** Consensus clustering. `consensus` accumulates a
  co-association matrix over an ensemble of full or partial partitions, which
  D3a sweeps and D3b resamples both produce, extracts a partition by
  majority-evidence components or by any matrix-consuming `ClusteringSpec`,
  and scores every member against the result with A2's agreement machinery.
- **D3d (shipped 5.18.0).** The `oecluster` command line. `cluster`,
  `select-parameter`, `stability` and `consensus` expose the roster and the
  three slices above as subcommands, with `algorithms` deriving the schema
  of each entry from its signature so a new roster entry appears without an
  edit. It reads a `.npz` written by the Python API, and `oepdist`'s `.npy`
  and `.bin` through their JSON sidecar; `.csv` is refused because oepdist
  writes it unquoted at eight significant digits with no provenance, and a
  file the sidecar proves holds similarities is refused outright, while an
  unproven orientation warns and proceeds. `oepdist` computed distances and
  nothing drove clustering from the command line before this.

Sequenced after D2 because consensus and stability are both defined over the
algorithm roster, and that roster was not final until D2 landed.

## D4 - Out-of-core and streaming clustering

*Re-derived.* When this section was re-derived, every algorithm except BitBirch
and Murcko required a materialized distance matrix. D2 changed that before D4
began: sphere exclusion's input and permutation orders, `jarvis_patrick` and
`leiden` already ran from a comparison. The out-of-core half was met as well,
since `pdist(..., output=path)` writes a memory-mapped matrix that every
algorithm taking a matrix accepts. What remained was never materializing the
matrix at all. The remaining algorithms divide by the structure they need, so
D4 is three slices:

- **D4a (shipped 5.19.0).** Threshold-graph streaming. `butina`, `dbscan` and
  `sphere_exclusion(order="neighbors")` need no pairwise structure beyond the
  threshold neighbor graph, which a comparison now builds in two passes over
  every pair, one to count each item's neighbors and one to record them. The
  graph's exact size is known before it is allocated, so a too-dense run is
  refused rather than attempted. The clustering engines are unchanged, and
  the result equals the matrix path's over a matrix filled through `Compare`.
- **D4b (shipped 5.20.0).** MST streaming. `hdbscan` and single-linkage
  `agglomerative` build a minimum spanning tree with Prim's algorithm on a
  persistent thread team, from a matrix or straight from a comparison.
  HDBSCAN's core distances come from one pass over every pair on a schedule
  whose concurrent tiles share no item. Matrix single linkage moved onto the
  same tree, so its merges at a tied height now follow the tree's order.
- **D4c (planned).** The algorithms that cannot stream exactly: `k_medoids`
  and complete, average and weighted linkage.

## Deferred defect backlog

`docs/hyperpowers/carry-forward-comparison-layer-expansion.md` holds findings
accepted as real during sub-project F and deliberately deferred. It stood at 73
items when F closed and has since grown to 80. The standing disposition is that
each item is fixed by whichever sub-project next touches its area, rather than
as a tranche of its own.

That file is gitignored, for the same reason the specs are. Items that outlive
their sub-project should be promoted to issues rather than left there.
