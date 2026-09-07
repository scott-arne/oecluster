# Changelog

This file starts at 5.0.0; earlier releases are not recorded here.

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
  and C#N score nonzero against themselves under the
  default `combo_norm` (0.500, 0.507 and 0.0067), so they stamp
  `zero_self = No` and the metric gate refuses to cluster them with no override
  available.
- ROCS scores changed, silently. 4.2.3 named a color force field on the overlay
  options but never assigned color atoms to the molecules, so
  `GetColorTanimoto()` answered 0.0 for every pair -- a molecule against itself
  included, which put a floor of 0.5 under every `combo_norm` self-distance. The
  constructor now runs `OEOverlapPrep` over its copies. On Omega-embedded phenol
  against catechol the `color` distance moves from 1.0 to 0.400135 and
  `combo_norm` from 0.519059 to 0.220440; `shape` is not exempt, moving from
  0.038117 to 0.040745 on that same pair. Phenol's own `combo_norm`
  self-distance falls from 0.500 to 0.000. A 4.2.3 script produces different
  numbers here with no error and no warning.
- `ROCSComparison` no longer mutates the caller's molecules. It deep-copies each
  input, so the dimension refresh and the newly added color-atom preparation,
  both of which write into the molecule, act on the copies rather than on
  objects the caller still owns.
- An index past the end is refused rather than answered. Neither
  `PairwiseComparison::Compare` nor `StorageBackend::Get` range-checked its
  arguments in 4.2.3: `Compare` read off the end of its container, quietly
  returning 0.0 from the descriptor comparison and crashing the process in
  others, and `Get` returned 0.0 for any `i == j` because the diagonal shortcut
  ran before anything looked at the bounds -- `Get(1000, 1000)` called a
  nonexistent sample identical to itself on a three-sample matrix. All five
  comparison classes now throw `ComparisonError` and all three storage backends
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
