# Changelog

This file starts at 5.0.0; earlier releases are not recorded here.

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
