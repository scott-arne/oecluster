# Changelog

This file starts at 5.0.0; earlier releases are not recorded here.

## [5.7.0] - 2026-09-26

### Added

- `mcs`, a maximum-common-substructure comparison for `pdist` and `cdist`, with
  the C++ class `MCSComparison` configured through `MCSOptions`. The score is
  Tanimoto over matched bonds, `c / (|A| + |B| - c)`, searched with complete
  cycles. It is topological: coordinates are never read, hydrogens are
  suppressed wherever they can be folded into a heavy atom (isotopic ones
  included; bridging and charged hydrogens survive), and a multi-conformer
  molecule is scored once rather than once per pose. Molecules with no bonds
  after suppression are refused at construction, because bond Tanimoto has a
  zero denominator there.
- `search_mode` selects `approximate` (the default) or `exhaustive`. This
  inverts the toolkit's own default deliberately: exhaustive search costs one to
  three orders of magnitude and is not reliably better, having returned a
  smaller match than approximate on a large macrolide pair while taking 16.2 s
  against 8.5 ms.
- `match_level` selects `default`, `exact` or `loose` matching strictness, and
  `max_matches` bounds how many matches one directed search enumerates.
- Each pair is searched in both directions and the larger match is taken.
  Approximate search is asymmetric -- the two directions can return different
  match sizes -- and `pdist` fills only one triangle, so a one-direction score
  would depend on input order.
- The matrix reports `triangle` as `unknown`: no violation appeared in 74,400
  ordered triples, but the inclusion-exclusion bound a genuine set intersection
  satisfies was violated 66 times over 59,280 triples, so the metric proof is
  unavailable. The five entry points that gate on metric facts -- `butina`,
  `dbscan`, `hdbscan`, `agglomerative` and `cluster_report` -- accept an MCS
  matrix without `allow_nonmetric=True`, as do `k_medoids`,
  `activity_landscape` and `modelability`, which do not gate at all.
- `MCSComparison::Clone()` deep-copying its molecule snapshots is now asserted
  rather than assumed. Two clones are taken, and their snapshot addresses are
  required to be disjoint from each other and from the parent's. This needed a
  test-only accessor, because the property has no consequence any score can
  show: an aliasing clone returns identical numbers, keeps the same molecules
  alive and reports the same `Size()`. Disjointness *between the clones* is the
  operative half -- `pdist` and `cdist` give each worker its own clone and
  never dereference the parent inside the parallel region, so a `Clone()` that
  deep-copied once and then handed every later caller the same snapshots would
  put one molecule set under every thread. Nothing in the suite had
  distinguished any of this: with `Clone()` neutralised to alias the parent's
  snapshots, all 685 other C++ tests still passed, so the documented
  thread-safety guarantee had rested on reading the code.
- A `tsan` CMake preset builds the C++ tests under ThreadSanitizer into their
  own `build-tsan/` tree, with the procedure documented in
  `docs/developer.md`. It covers the concurrency this project owns -- the
  thread pool, the storage backends, the progress callback and the
  clone-distribution loops -- and reported no races. The OpenEye libraries are
  prebuilt and uninstrumented, so it can say nothing about the toolkit's own
  internals, and it does not verify the `Clone()` deep copy above; the test
  does that.

### Fixed

- `pdist` now raises on `cutoff > 0` together with `similarity=True` instead of
  silently discarding the highest scores. Sparse storage zeroes values *above*
  the cutoff, which drops the far pairs of a distance matrix but the near pairs
  of a similarity matrix. The result was not an obviously broken matrix of
  zeros -- everything below the cutoff came back untouched, so the call looked
  plausible with exactly the most-similar pairs replaced by `0.0`. **Calls that
  would have used sparse storage -- a named comparison with no `output=` -- now
  raise `ValueError`**, matching the guard `cdist` has always had. This affected
  every similarity-capable comparison, `fingerprint` included, not only the
  `mcs` comparison new in this release.

### Changed

- `MCSComparison` and `ROCSComparison` now accept `OEGraphMol` as well as
  `OEMol`, snapshotting either into their own storage. Every example in the
  documentation builds `OEGraphMol`, and these were the two comparisons in
  scope for that change. `rmsd` still requires `OEMol` and continues to reject
  an `OEGraphMol` list at construction. A multi-conformer `OEMol` still binds
  to the `OEMol` overload and keeps its whole ensemble, so no Python call that
  worked before returns anything different now.

  C++ callers have to spell the conversion. The new overload takes
  `std::vector<OEChem::OEMolBase*>`, and an `OEGraphMol*` does not convert to
  an `OEMolBase*` implicitly; a pointer `static_cast` is rejected too. Write
  `&static_cast<OEChem::OEMolBase&>(graph_mol)`. Python callers are unaffected,
  the bindings converting for them.
- C++ callers keep source compatibility, but it was not free. A second vector
  overload made `MCSComparison({})` and `MCSComparison({nullptr})` ambiguous --
  either vector type can be brace-initialized from those, so neither candidate
  wins -- and the same held for `ROCSComparison`. Both classes therefore carry
  a `std::initializer_list<std::shared_ptr<OEMol>>` constructor that restores
  those spellings. It delegates to the `OEMol` overload and is hidden from the
  bindings, so a braced list of molecules -- which binds it now rather than the
  vector overload -- keeps every conformer and scores exactly what it scored
  before.
- Passing `rocs` a list whose **first** element is an `OEGraphMol` now selects
  the permissive overload for the entire list, so an `OEMol` later in that list
  is reduced to its active conformer with no error raised. The reverse ordering
  still raises. See the `rocs` section of `docs/python-api.md`.

## [5.6.0] - 2026-09-23

### Added

- `murcko_scaffolds`, in C++ and Python, assigning each molecule its
  Bemis-Murcko scaffold as a canonical SMILES. `scaffold="framework"` keeps ring
  systems plus their linkers; `scaffold="generic"` reduces that framework to its
  topology, every heavy atom carbon and every bond single. A molecule with no
  ring system yields the empty string, which is the "missing scaffold"
  convention `scaffold_agreement` already consumes -- the library can now
  produce the `scaffold_labels` it has consumed since 5.3.0 in three places.
- `murcko`, clustering molecules by scaffold identity; the C++ entry point is
  `murcko_cluster`, configured through `MurckoOptions`. It is the first
  partition in the library computed from chemical structure rather than from a
  distance matrix, so it takes molecules directly. Labels are the rank of each
  scaffold in the sorted distinct set, making the labeling independent of input
  order; acyclic molecules are noise. `MurckoResult` adds `scaffolds` and
  `cluster_scaffolds`.

### Notes

- Molecules are taken as given: no salt stripping and no largest-component
  selection. An acyclic counter-ion carries no ring system and so contributes
  no scaffold region, which makes the usual salt forms come out identical to
  the free base; only a ring-bearing counter-ion adds a `.`-joined component.
  Stereochemistry is dropped and explicit hydrogens are suppressed, so scaffold
  identity does not depend on how a molecule was read. Hydrogen counts and
  formal charges are recomputed on the atoms a sidechain cut touched, so
  diphenyl sulfone, diphenyl sulfoxide and diphenyl sulfide share one framework
  scaffold; a charged ring atom whose bonds all survive the cut keeps its
  charge.
- Extraction parallelizes across molecules only when the OpenEye memory-pool
  mode reports a thread-safe setting; otherwise it runs serially and returns the
  same answer. `num_threads` is clamped to the molecule count and to a multiple
  of the hardware concurrency.
- `oecluster` now links `OpenEye::OEMedChem`.

## [5.5.0] - 2026-09-20

### Added

- `k_medoids`, in C++ and Python, placing exactly `n_clusters` medoids -- real
  members of the input rather than synthetic averages -- by minimizing the sum
  of every item's distance to its assigned medoid. Initialization is greedy PAM
  BUILD (the default), deterministic farthest-first from the global medoid, or
  caller-supplied indices; the swap phase is FastPAM1, an exact algebraic
  reformulation of textbook PAM that gets all `k` deltas for a candidate in one
  scan rather than one scan per pair. `KMedoidsResult` carries `medoids`,
  `cost`, `n_iterations` and `converged`.
- `k_medoids()` is the first clustering entry point that does not assume a
  metric. Its objective is a sum of distances and its swap step compares two
  such sums, so no step appeals to the triangle inequality; it uses the same
  weaker gate as `activity_landscape()` and `modelability()` and takes no
  `allow_nonmetric` parameter, because there is no assumption for a flag to
  override. `cluster_report()` still requires the stronger gate, so a Dice
  matrix that clusters here may need `allow_nonmetric=True` there.
- `converged` on a k-medoids result is a verified claim rather than a loop-exit
  flag: it is only ever `True` when no single medoid swap would lower the cost
  -- verified for non-trivial cases (`n_clusters < n`) by recomputing every
  candidate total from scratch, or holding vacuously when the identity
  partition is returned (`n_clusters` equals the item count). Reaching
  `max_iterations` reports `converged=False` and does not raise. Output is
  byte-identical across runs, thread counts and chunk sizes.

## [5.4.0] - 2026-09-17

### Added

- `sar_coherence`, in C++ and Python, decomposing an activity vector across a
  labeling: `eta_squared` is the fraction of activity variance that falls
  between clusters rather than within them, and `omega_squared` is the same
  quantity with the variance a random labeling of the same shape would explain
  removed, so it can be compared across labelings of different granularity.
  The input is a clustering result or a bare labeling, and a `nan` activity is
  missing data rather than a value. `SARCoherence.clusters` carries a
  per-cluster mean and standard deviation. `noise=` takes the same three
  spellings as the 5.3.0 agreement metrics, but the default here is
  `"excluded"` and not the `Singletons` reading that entry called "matching the
  rest of the library": noise is not a structural hypothesis, and promoting
  each noise point to a cluster of its own inflates eta_squared for a reason
  that has nothing to do with the labeling under test.
- `activity_landscape`, in C++ and Python, scoring the structure-activity
  landscape from a distance matrix with no clustering in between: the count and
  per-pair density of activity cliffs, the maximum and mean SALI over the pairs
  where the ratio is defined, and RMODI. Pairs at zero distance are excluded
  from both SALI figures -- the ratio is infinite where the two activities
  differ and `0/0` where they agree -- and reported as
  `num_zero_distance_pairs`, while still counting as cliffs when their activity
  difference qualifies.
- `modelability`, in C++ and Python, reporting MODI over a per-sample class
  annotation: the mean over classes of the fraction of each class's members
  whose nearest scored neighbour shares their class, with a per-class
  breakdown. An empty class string is missing data rather than a category,
  matching `scaffold_agreement`.
- `OECluster::detail` activity helpers in `src/clustering/ActivityMetrics.h`,
  and the public surface in `include/oecluster/clustering/SARCoherence.h`,
  reachable through the `oecluster.h` umbrella header.
- A distance-matrix gate for the two matrix entry points. They rank and
  threshold distances rather than assuming a metric, so they waive the
  triangle-inequality warnings and take no `allow_nonmetric` parameter, and
  they refuse a matrix stamped `data_integrity == 'subset_scored'` outright: a
  nearest neighbour chosen among distances scored on different feature subsets
  is not a nearest neighbour.

## [5.3.0] - 2026-09-14

### Added

- `partition_agreement` and `scaffold_agreement`, in C++ and Python, scoring
  the agreement between two labelings of the same samples: adjusted Rand index,
  Fowlkes-Mallows, normalized mutual information, homogeneity, completeness,
  V-measure, and adjusted mutual information behind an opt-in flag. The input
  is labels and nothing else -- no distance matrix -- so two clustering methods
  can be compared without owning one. `scaffold_agreement` scores a clustering
  against a per-sample scaffold annotation, where an empty string is missing
  data rather than a category.
- `NoiseHandling` selects how negatively-labelled samples enter the contingency
  table: `Singletons` (the default, matching the rest of the library),
  `Grouped` (matching scikit-learn's reading of a -1 label), or `Excluded`.
  Python callers pass `noise="singletons"`, `"grouped"` or `"excluded"`.

### Fixed

- `include/oecluster/oecluster.h` now includes `ClusterReport.h`. The umbrella
  header the documentation names as the entry point had never carried it, so
  the whole cluster-quality surface was unreachable through it.
- `scikit-learn` is now declared in the `dev` extra. Three test modules had
  imported it and CI had installed it since before 5.0.0 -- this release brings
  a fourth -- but a fresh `uv pip install -e '.[dev]'` produced a checkout
  whose parity suites failed to collect.

## [5.2.0] - 2026-09-12

### Changed

- A singleton cluster's silhouette is now 0 rather than 1.0, following
  Rousseeuw's definition and matching scikit-learn's `silhouette_score`. **This
  changes the number `cluster_report` returns for any clustering that has a
  size-1 cluster**, in both `ClusterReport::silhouette` and the affected
  `ClusterRecord::silhouette`; nothing raises and no other field moves. A
  singleton's within-cluster mean arrives at the formula as 0.0 for want of an
  own-cluster pair, not because its neighbours are coincident, so the general
  `(b - a) / max(a, b)` read it as `b / b` and handed a cluster of one the
  perfect score. Averaged into the scalar, that pulled the mean **up** on
  exactly the fragmented clusterings the scorecard exists to discriminate: on a
  four-point fixture split `{0,1}`, `{2}`, both real members scored 0.75 and the
  report read 0.833, above two of the three terms it averaged. It now reads 0.5.
  The old value was never a defensible reading of the index -- no published
  definition awards a lone point the maximum -- so this is corrected rather than
  made configurable. Reports produced before this release are not comparable
  with ones produced after it whenever `num_singletons` is non-zero; the two are
  identical when it is zero. The correction applies to the silhouette alone:
  `silhouette_like_score` on `RepresentativeMetrics` is a different quantity
  over representatives and is unchanged.

## [5.1.0] - 2026-09-10

### Added

- Seven internal cluster-validity indices on `ClusterReport`. Always computed:
  `calinski_harabasz_medoid`, `davies_bouldin_medoid`,
  `dunn_mean_separation_mean_diameter`,
  `dunn_medoid_separation_medoid_spread` and `point_biserial`. Behind the new
  `compute_pair_rank_indices` option: `c_index` and `baker_hubert_gamma`. The
  two are behind one flag rather than two because both come off the same sorted
  pair arrays -- once those are paid for, the second index is nearly free, and a
  separate flag would advertise a saving that does not exist. Those arrays hold
  `Nc * (Nc - 1) / 2` doubles between them, roughly 400 MB at 10,000 clustered
  points and 10 GB at 50,000, but only the between-cluster one is the flag's own
  cost: the within-cluster array is built on every call, with or without the
  flag, because `median_intra_distance` is taken over it. The flag's share of
  the figure runs from nothing on a single-cluster result to nearly all of it
  when the clusters are small, and it is off by default for the second case.
- The Calinski-Harabasz and Davies-Bouldin indices are **medoid-substituted**:
  the published definitions use centroids, which do not exist for a distance
  matrix, so each cluster's medoid stands in for its centroid.
  Calinski-Harabasz also has a grand-mean term, and the global medoid stands in
  for that; Davies-Bouldin has no such term. The values are not comparable with
  published or scikit-learn figures. These two fields and
  `dunn_medoid_separation_medoid_spread` always use the true medoid and ignore
  `representative_method`, so a minimax-configured report still reports
  medoid-based values for them.
- A per-cluster record table, `ClusterReport.records`, behind the new
  `compute_per_cluster_records` option. Each `ClusterRecord` carries the
  cluster's label, size, representative, intra-distance mean and median, radius,
  diameter, mean representative distance, nearest cluster and distance,
  silhouette, and boundary-violation count. `boundary_violations` on a record
  counts pairs *involving* that cluster, so the sum over records is twice the
  scorecard's, which counts each pair once.
- `ClusterReport.requested`, a `ClusterReportRequested` recording which optional
  computations the caller asked for. It records the request, not the outcome, so
  a NaN can be read unambiguously: false means nobody asked, true with NaN means
  asked and undefined.
- `ClusterReport.noise_coverage_at`, the coverage curve restricted to noise
  points. Its length always matches `coverage_at`: the threshold count when the
  clustering has at least one cluster, and empty when it has none. Every entry
  is NaN when the clustering has no noise -- not 0.0, which would read as "no
  noise point is covered" rather than "there is nothing to cover".
- Python exports `ClusterRecord` and `ClusterReportRequested`, both
  `typing.NamedTuple` subclasses, so `pandas.DataFrame(report.records)` works
  without a conversion step. pandas is not a dependency.

### Changed

- `cluster_report` now refuses a `ClusteringResult` whose labels and members
  describe different partitions. Newly `std::invalid_argument`: a sample in more
  than one cluster; a member whose label entry does not match the cluster
  holding it; a clustered sample omitted from every member list; a non-noise
  label naming no cluster; a noise-labelled sample sitting in a cluster. Newly
  `std::out_of_range`: a member index inside the storage range but past the end
  of a shorter label vector, which previously read out of bounds. The
  empty-cluster, duplicate-within-a-cluster and beyond-`NumSamples()` refusals
  are unchanged and keep their present types and messages. The check now also
  runs when the member list is empty, so a result with a non-noise label and no
  clusters refuses instead of returning an all-NaN report. The rules added here
  are satisfied by every algorithm the library ships, so the new refusals can
  only fire on a hand-built result. That is not true of the unchanged
  empty-cluster refusal: `bitbirch_refine` can return an emptied leaf subcluster
  as an empty member list, and `cluster_report` rejects such a result.
  `ClusterReportRefusesEmptiedRefinementSubclusters` in
  `tests/cpp/test_bitbirch_clustering.cpp` pins that gap.
- `cluster_report` now refuses a non-finite distance that reaches a reported
  value with `std::invalid_argument`. Previously a NaN reached `std::sort` through the
  intra-distance median, which is undefined behaviour, so this removes a hazard
  rather than a defined result. A caller who loaded a matrix with holes will see
  it.
- `cluster_report` now raises `ValueError` for a NaN coverage threshold
  (`coverage thresholds must not be NaN`) or a NaN `boundary_threshold`
  (`boundary_threshold must not be NaN`). Both were previously accepted, and
  each failed in its own way. A NaN coverage threshold produced a computed value
  that `compare_reports(...).to_table()` could not match back to its threshold,
  so the cell rendered as though the question had never been asked. A NaN
  `boundary_threshold` failed worse: no distance compares against it, so every
  pair cleared the boundary and the report returned `boundary_violations = 0`.
  That is a plain scalar row in the table, published as an ordinary `0` -- a
  clean bill of health for a question that was never answerable, with nothing
  marking it as suspect, where the coverage case at least rendered `--`.
  Positive infinity remains accepted for both and is unaffected. Negative
  infinity is refused, as it always was, by the non-negative check rather than
  by anything added here.
- The C++ `cluster_report` now throws `std::invalid_argument` for a NaN
  `boundary_threshold` (`cluster_report: boundary_threshold must not be NaN`) or
  a NaN entry in `coverage_thresholds` (`cluster_report: coverage threshold <i>
  must not be NaN`, naming the offending index). **This refuses C++ calls that
  5.0.0 accepted.** A NaN threshold never produced a NaN and never raised: every
  comparison against it is false, so the call returned `boundary_violations = 0`
  and `coverage_at` entries of `0.0` -- a plausible, in-range answer to a
  question that was never answerable. Positive infinity is still accepted and
  still means what it always meant: every cross pair is a boundary violation,
  and every sample is covered. The checks sit inside the members-non-empty
  block, behind only the cluster count and the unsupported-method refusal and
  ahead of every distance read, so a partition with no clusters is accepted and
  from the first cluster onwards every call is refused. At `K == 0` no
  comparison consumes either option -- `coverage_at` stays empty and
  `boundary_violations` is zero for want of a pair to count, not for want of a
  comparison that held -- so there is no wrong number for the NaN to hide
  behind, and refusing it would be over-refusal. `coverage_thresholds` is still
  echoed back verbatim in the report, so a NaN passed there is returned in that
  field; it is the caller's own value coming back, not a computed one. That is a
  placement rule and not a test of whether the particular call would have
  consumed the value -- a single-cluster result never reads `boundary_threshold`
  at all, and a NaN one is refused there anyway. Which of the two refusals a
  call that is NaN in both places receives is deliberately not a guarantee, and
  no test pins it. **Python callers see no change.** The wrapper already refuses
  both with its own `ValueError` before the call reaches C++, so this closes a
  C++/Python asymmetry rather than altering Python behaviour.
- `cluster_report` now raises `TypeError` for a `treat_noise_as_singletons` that
  is not a `bool` or a `numpy.bool_`. **This refuses calls that 5.0.0 accepted.**
  The keyword predates this release and was coerced with a bare `bool(...)`, so
  `treat_noise_as_singletons="no"` read to a caller as off while switching the
  folding on, and `0`, `1` and `None` were silently reinterpreted the same way.
  All four now raise. `numpy.bool_` is still accepted, on the same terms as the
  new `compute_pair_rank_indices` and `compute_per_cluster_records` keywords and
  as `allow_nonmetric`; the message carries the offending value as well as its
  type, because `numpy.bool_.__name__` is itself `bool`. The check is placed in
  signature order, behind the `representative_method` refusal and ahead of the
  `num_threads` one, so a call that is wrong in two ways names the keyword the
  caller wrote first.
- `compare_reports(...).to_table()` renders a cell as `None` rather than `nan`
  when the report never asked the question, in any of three ways: it did not
  request an opt-in metric (`c_index`, `baker_hubert_gamma`); it does not carry
  that coverage threshold at all; or it carries the threshold but answered
  nothing at it, as when the clustering has no clusters and both coverage curves
  come back empty. `nan` keeps its single meaning of asked-and-undefined; the
  two states were previously indistinguishable, which is why the change was
  made. A caller who consumed the table as floats and used `math.isnan` to
  detect absence must now also test for `None`, and arithmetic on a cell without
  that test raises `TypeError` where it used to propagate a NaN. `__repr__`
  renders `None` as `--`. The table also gains rows for the new metrics, for
  `noise_coverage_at`, and a `requested_pair_rank_indices` row stating in one
  line what the `None` in those two metric rows means. There is deliberately no
  matching `requested_per_cluster_records` row: that flag governs `records`,
  which the table does not carry.
- The default path's cost changes. It grows no pairwise-sized allocation and no
  additional distance sweep, and pass fusion removes two of the three
  cross-cluster walks, one of the two medoid selections, and a factor of the
  cluster count from the coverage loop, so the expected direction is faster.
  Unchanged cost is not claimed.
- Python `MemoryError` now surfaces from `cluster_report` for `std::bad_alloc`
  and `std::length_error`. Both previously reached Python as `RuntimeError`,
  indistinguishable from a validation failure.
- `StorageBackend::Set` now refuses a pair it cannot store, on all three
  backends: `std::out_of_range` for an index at or beyond `NumSamples()`, and
  `std::invalid_argument` for an in-range diagonal. Both reach Python as
  `RuntimeError`. **This refuses calls that 5.0.0 accepted.** `Get` has been
  range-checked since 5.0.0 and the write side never was, and it fails worse
  than the read side: `Get` answered about a pair that does not exist, `Set`
  destroyed one that does. Out of range the condensed index lands on an
  unrelated pair -- `CondensedIndex(5, 2, 5) == 9 == CondensedIndex(5, 3, 4)` --
  or, further out, leaves the allocation entirely: `Set(0, 5000, v)` on a
  five-sample `DenseStorage` wrote roughly 40 KiB past the buffer without
  faulting, and the `MMapStorage` equivalent runs past the mapping. On the
  diagonal both indices are in range yet the pair owns no slot at all, since
  `Get` answers `i == j` from a shortcut rather than from memory, so
  `CondensedIndex(5, 2, 2) == 6` handed the write the `(1, 4)` pair. The
  `i != j` precondition was a bare `assert`, compiled out of release builds. The
  range check runs ahead of the diagonal check, so an out-of-range diagonal is
  reported as the range error it also is, and ahead of `SparseStorage`'s cutoff
  shortcut, so whether a bad call is diagnosed does not depend on the value it
  carried. `pdist` keeps its separate storage-size check: `Set`'s own range
  check cannot see a storage larger than the comparison, because every index
  that loop produces is in range for the oversized storage.

## [5.0.0] - 2026-09-06

### Changed

- The OEFP requirement moved from `oefp>=0.2.4` to `oefp==0.3.0`, and the
  vendored source tag from `v0.2.4` to `v0.3.0`. The pin is exact rather than a
  range: oecluster compiles OEFP's core into its own extension and exchanges
  raw fingerprint-batch pointers with the separately compiled wheel, so the
  compiled-against and installed versions must be identical. The extension
  compares all three version components the first time such a pointer crosses
  -- passing an `oefp.OEFPBatch` to `bitbirch`, for instance -- and raises
  `ImportError` on any difference. That check does not run on `import
  oecluster`, and the fingerprint comparisons never reach it because they build
  their own fingerprints, so a range pin left a mismatched wheel working right
  up to the first batch call.
- `max_distance` no longer sets the Morgan radius. It previously meant both the
  Morgan radius and the atom-pair maximum graph distance. A separate `radius=`
  now carries the Morgan meaning, `max_distance=` applies to the atom-pair
  family only, and naming `max_distance` together with `fp_type="morgan"`
  raises `TypeError`. This is not a rename: `max_distance` still exists, with a
  narrowed meaning.
- The atom-pair window defaults changed from 0-2 to 1-30, matching OEFP's own
  defaults (`min_distance` 0 to 1, `max_distance` 2 to 30). A fingerprint built
  with the defaults is therefore not comparable to one built under 4.x.
- `ROCSComparison` refuses input it used to accept. A molecule whose dimension
  attribute, recomputed from its coordinates, is below three now raises at
  construction: `ComparisonError` in C++, which the bindings render as
  `RuntimeError`. 4.2.3 had no such check and ran on 2D input. Separately, the
  constructor now measures the diagonal, and linear species such as N#N, O=C=O
  and C#N score nonzero against themselves under the default `combo_norm`
  (0.500, 0.507 and 0.0067), so they stamp `zero_self = No` and the metric gate
  refuses to cluster them with no override available.
- ROCS scores changed, silently. 4.2.3 named a color force field on the overlay
  options but never assigned color atoms to the molecules, so
  `GetColorTanimoto()` answered 0.0 for every pair -- a molecule against itself
  included, which put a floor of 0.5 under every `combo_norm` self-distance. The
  constructor now runs `OEOverlapPrep` over its copies. On Omega-embedded phenol
  against catechol the `color` distance moves from 1.0 to 0.400135 and
  `combo_norm` from 0.519059 to 0.220440. `shape` is exempt: keeping hydrogens
  reproduces the shape numbers this comparison returned before color prep
  existed. Phenol's own `combo_norm` self-distance falls from 0.500 to 0.000. A
  4.2.3 script produces different numbers here with no error and no warning.
- `ROCSComparison` no longer mutates the caller's molecules. It deep-copies each
  input, so the dimension refresh and the newly added color-atom preparation,
  both of which write into the molecule, act on the copies rather than on
  objects the caller still owns.
- An index past the end is refused rather than answered. Neither
  `PairwiseComparison::Compare` nor `StorageBackend::Get` range-checked its
  arguments in 4.2.3. `Compare` read off the end of its container, and what it
  did next varied with how far past the end the index landed rather than with
  which comparison class was in use: a two-molecule fingerprint comparison
  quietly returned 0.0 for `Compare(1000000, 1000001)`, and that same class took
  the process down with SIGSEGV once the index was far enough out to leave the
  mapping. `Get` returned 0.0 for any `i == j` because the diagonal shortcut ran
  before anything looked at the bounds -- `Get(1000, 1000)` called a nonexistent
  sample identical to itself on a three-sample matrix. All five comparison
  classes now throw `ComparisonError` and all three storage backends
  `std::out_of_range`, both of which the bindings render as `RuntimeError`.
- `oepdist fp` refuses an option the selected family or storage does not read,
  where it used to accept and discard it. `--max-distance 4 --fp-type morgan`
  produced a matrix and a JSON sidecar byte-identical to the same run without
  the flag; 4.2.3 had no `--radius` and documented `--max-distance` as "Morgan
  radius or maximum Atom Pair graph distance", so a 4.x script carried straight
  over gets a different fingerprint than it asks for. `--radius`,
  `--min-distance`, `--max-distance`, `--torsion-atom-count`, and `--numbits`
  under `--storage sparse` or `--storage sparse_count` now exit 1 with a message
  naming the replacement. These are the rules the Python surface has applied
  since they were introduced; the CLI never ran them.
- `butina`, `dbscan`, `hdbscan`, `agglomerative` and `cluster_report` refuse
  matrices they used to accept. Each now calls the metric gate described under
  Added, so a similarity matrix -- which 4.2.3 clustered without complaint --
  and any matrix whose recorded facts fail the gate raise `ValueError`.

### Added

- Descriptor comparisons (`DescriptorComparison`, `descriptor_statistics`) and
  coordinate RMSD comparisons (`RMSDComparison`).
- The descriptor options that take a sequence -- `sources`, `columns`, `groups`,
  `variances` and `inverse_covariance` on the comparison paths, and `sources`,
  `columns` and `groups` on `descriptor_statistics` -- refuse an empty sequence
  with `ValueError`. The C++ layer cannot tell an empty sequence from an omitted
  option, so accepting one would silently resolve to the default.
- A metric-capability gate. Comparisons record what their distances guarantee,
  and the clustering entry points raise `ValueError` on a matrix whose recorded
  facts violate their assumptions. `allow_nonmetric=True` overrides the second
  tier of checks -- among them the triangle inequality and subset-scored
  distances -- but not the first, which covers distance orientation, the zero
  self-distance and non-finite entries.
- `SymmetricDistanceMatrix.from_condensed`.
- The `metric_probe`, `probe_violations` and `probe_sampled` facts.
- A `storage=` axis on `FingerprintComparison`, accepting `binary`, `count`,
  `sparse` and `sparse_count`, and a `topological_torsions` fingerprint family
  alongside `morgan` and `atom_pair`. `topological_atom_pair` is accepted as a
  further spelling of `atom_pair`.
- A shared metric table behind both comparison surfaces. The fingerprint surface
  resolves 17 scalar metric names where 4.2.3 recognized four -- `tanimoto`,
  `jaccard`, `dice` and `manhattan` -- and refused `euclidean` outright; the
  descriptor surface resolves 10, seven of them shared with fingerprints.
  `tversky` is the one similarity-only entry: it has no distance form, so it
  requires `similarity=True` and refuses `similarity=False` rather than
  answering a distance request with a similarity. That is the mirror of the
  refusal `jaccard` gives `similarity=True`.
- `SymmetricDistanceMatrix.from_file` reads sparse-storage matrices, and both
  matrix classes read the recorded facts payload. A file with no `storage_kind`
  key is read as dense, so matrices written before sparse serialization existed
  still load unchanged.
