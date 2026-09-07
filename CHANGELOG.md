# Changelog

This file starts at 5.0.0; earlier releases are not recorded here.

## [5.0.0] - 2026-09-06

### Changed

- The OEFP requirement moved from `oefp>=0.2.4` to `oefp>=0.3.0,<0.4`, and the
  vendored source tag from `v0.2.4` to `v0.3.0`. An 0.2.x wheel no longer
  satisfies the install, and 0.4 is excluded ahead of its release: oecluster
  compiles OEFP's core into its own extension and exchanges batch pointers with
  the installed wheel, so the two must share a minor series.
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
  attribute, recomputed from its coordinates, is below three now raises
  `ComparisonError` at construction; 4.2.3 had no such check and ran on 2D
  input. Separately, the constructor now measures the diagonal, and linear
  species such as N#N, O=C=O and C#N score nonzero against themselves under the
  default `combo_norm` (0.500, 0.507 and 0.0067), so they stamp
  `zero_self = No` and the metric gate refuses to cluster them with no override
  available.
- `ROCSComparison` no longer mutates the caller's molecules. It deep-copies each
  input, so the color-atom preparation and the dimension refresh, which write
  into the molecule, act on the copies rather than on objects the caller still
  owns.
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
- `SymmetricDistanceMatrix.from_file` reads sparse-storage matrices, and both
  matrix classes read the recorded facts payload. A file with no `storage_kind`
  key is read as dense, so matrices written before sparse serialization existed
  still load unchanged.
