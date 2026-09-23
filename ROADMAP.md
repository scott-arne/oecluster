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

    F -> A -> D1 -> C -> E -> D2 -> B -> D3/D4

| Piece | Content | Status |
| --- | --- | --- |
| F | Comparison layer expansion | shipped 5.0.0 |
| A1 | Internal indices, per-cluster record table, `noise_coverage_at` | shipped 5.1.0, 5.2.0 |
| A2 | Partition agreement: ARI/AMI/NMI/V-measure/Fowlkes-Mallows, scaffold-ARI | shipped 5.3.0 |
| A3 | SAR coherence: eta-squared, omega-squared, SALI, cliff density, MODI/RMODI | shipped 5.4.0 |
| D1 | k-medoids/PAM: exactly `k` clusters with real-member centers | shipped 5.5.0 |
| C1 | Murcko scaffold assignment and scaffold-identity clustering | shipped 5.6.0 |
| C2 | MCS-based clustering | next |
| E | Diversity metrics and library-scale selection | planned |
| D2 | Graph and leader algorithms: Leiden, Jarvis-Patrick, leader/DISE | planned |
| B | Fingerprint-native O(N) counterparts of the A metrics | planned |
| D3 | Workflow layer: clustering CLI, parameter selection, consensus, stability | planned |
| D4 | Out-of-core and streaming clustering | planned |

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
scaffold assignment and MCS-based clustering.

Every algorithm in the library today is distance-driven. `scaffold_labels` is
consumed as an input in three places - scaffold-ARI in A2, the weighted-medoid
representative, and scaffold purity in the report - but the library has never
produced them. C closes that asymmetry and gives E a chemically meaningful axis
to diversify along.

C1 shipped in 5.6.0 as `murcko_scaffolds` and `murcko`, closing the producer
side of that asymmetry. C2, MCS-based clustering, is the remaining half.

## E - Diversity metrics and library-scale selection

*Pinned by the A1 and A3 specs.* Adds #Circles, Vendi score and log-determinant
diversity, plus a selection module that operates over a whole library rather
than within a single cluster. The existing `select_representatives()` selects
within one cluster; E selects across a collection.

These metrics score a set rather than a partition, which is why they were moved
out of sub-project A.

## D2 - Graph and leader algorithms

*Pinned by the A1 and A2 specs.* Adds Leiden, Jarvis-Patrick and leader/DISE to
the algorithm roster.

Bootstrap-Jaccard cluster stability is deferred until after D2: it must re-run
the clustering algorithm, which inverts the current layering, and it is not
worth defining against a roster that is still growing.

## B - Fingerprint-native metrics

*Pinned by the A1 and A3 specs.* Adds O(N) fingerprint-direct counterparts of
the metrics defined in A, avoiding the materialized O(N^2) distance matrix that
the current forms require. A owns the definitions; B owns the scale path.

## D3 - Workflow layer

*Re-derived.* Four gaps that share a consumer rather than an implementation:

- A clustering CLI. `oepdist` computes distances; nothing drives clustering from
  the command line.
- Automatic parameter selection. A1's quality metrics exist, but nothing uses
  them to choose a Butina threshold or a DBSCAN `eps`.
- Consensus clustering over A2's agreement machinery.
- Stability resampling, including the bootstrap-Jaccard work deferred from A.

Sequenced after D2 because consensus and stability are both defined over the
algorithm roster, and that roster is not final until D2 lands.

## D4 - Out-of-core and streaming clustering

*Re-derived.* Every algorithm except BitBirch requires a materialized distance
matrix. D4 lifts that constraint at the algorithm level, building on B's
fingerprint-native path and the completed roster.

## Deferred defect backlog

`docs/hyperpowers/carry-forward-comparison-layer-expansion.md` holds findings
accepted as real during sub-project F and deliberately deferred. It stood at 73
items when F closed and has since grown to 80. The standing disposition is that
each item is fixed by whichever sub-project next touches its area, rather than
as a tranche of its own.

That file is gitignored, for the same reason the specs are. Items that outlive
their sub-project should be promoted to issues rather than left there.
