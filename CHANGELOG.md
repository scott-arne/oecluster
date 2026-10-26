# Changelog

This file starts at 5.0.0; earlier releases are not recorded here.

## [5.28.0] - 2026-10-26

### Changed

- **`butina`, `dbscan`, `hdbscan`, `agglomerative`, `bitbirch`,
  `bitbirch_recluster`, `bitbirch_refine`, `cluster_report`, `isim_report`,
  `activity_landscape` and `modelability` raise `TypeError` for a
  `num_threads` that is a float, a string or `None`.** They now coerce it with
  `operator.index`, as `k_medoids`, `pdist` and `cdist` already did, where
  they used `int()`. The first four coerce `chunk_size` the same way, as
  `k_medoids` does. A float used to be truncated, so `num_threads=2.5` was
  read as 2 and `num_threads=-0.5` as 0, which chooses the thread count
  automatically; an infinity raised `OverflowError` and a NaN `ValueError`. A
  string holding an integer, such as `"8"`, used to be parsed, and any other
  string raised `ValueError`. `None` already raised `TypeError`. Python bools
  and numpy integers are still accepted. On NumPy 2.3 and later a
  `numpy.bool_` now raises `TypeError` too, as it already did for the
  `num_threads` of `pdist`, `cdist` and `k_medoids`; earlier NumPy reads it
  as 0 or 1 with a `DeprecationWarning`.
- When the metric gate scans a distance matrix for NaN and infinite
  distances, which it does unless the matrix is stamped `nan_present` or an
  earlier check has refused it, the scan no longer allocates memory in
  proportion to the matrix. On dense and memory-mapped storage it tests
  1,048,576 values at a time, where it used to allocate a one-byte flag for
  every distance at once. On sparse storage it walks the stored entries in
  C++ and stops at the first non-finite one, where it used to build a Python
  list holding an `(i, j, value)` tuple for every stored entry. `pdist` uses
  the same scan when it checks a `nan_present` stamp against a dense or
  memory-mapped result, and so does the load-time correction under Fixed: in
  chunks for a dense or cross matrix, entry by entry for a sparse one.

### Fixed

- `SymmetricDistanceMatrix.from_file` and `CrossDistanceMatrix.from_file`,
  which `load_distance_matrix` calls, used to report a file's
  `data_integrity` stamp as saved, so a matrix holding a NaN could load as
  `complete`. A `complete` or `unknown` stamp now reads back as `nan_present`
  when the loaded matrix holds a NaN or an infinity. A `subset_scored` stamp
  is kept even then, and a scan that finds nothing changes no stamp, so
  `nan_present` over finite values stays `nan_present`. A cross matrix stored
  with a dtype other than integer or floating point is not scanned and keeps
  its stamp. The metric gate measured the values, not only the stamp, before
  this release too, so a call that runs it already refused such a matrix.
- `consensus_components` and `consensus_strength`, called from C++ or
  through the raw bindings in `oecluster.oecluster`, refuse a matrix of fewer
  than two items with "requires at least 2 items" whatever its storage. On a
  dense, memory-mapped or sparse matrix of 0 or 1 items they used to raise
  the refusal written for sparse storage, "requires contiguous storage; a
  SparseStorage matrix holds only the pairs below its cutoff", and
  `consensus_strength` had no item-count check at all. `consensus()` refuses
  fewer than two items itself and is unchanged.
- The `oecluster` CLI reads a `.npy` or `.bin` input that holds no values as
  a matrix of 0 items when its sidecar records `n_rows` as the integer 0. It
  used to refuse that input with "0 values describe 1 items, but the sidecar
  says 0". With a sidecar recording 1, such an input still reads as 1 item.
- `SparseStorage` refuses a sample count for which `n * (n - 1)` does not fit
  in `size_t`, which on a 64-bit platform is any count above 4294967296, with
  `std::overflow_error`; through the bindings that arrives as `RuntimeError`.
  Its pair count used to wrap silently, and its condensed index could:
  `SparseStorage(SIZE_MAX, 0.5).NumPairs()` returned 1. This supersedes
  5.25.0's note that `SparseStorage` was not covered.
- Documentation only: the `MMapStorage` constructor's comment now says how
  an existing file is reused, extended and zero-filled, or truncated.

## [5.27.0] - 2026-10-25

### Added

- The public header `oecluster/PoolMode.h`, which decides how many workers
  may run OpenEye code under the memory-pool mode of a call, and a Concurrency
  section in the Python API docs. The thread counts,
  timings and per-worker memory figures in the docs now say they assume a
  thread-safe pool mode.

### Changed

- **Calling any wrapped accessor, method or destructor with `None` as `self`
  raises `TypeError` ("invalid null reference").** SWIG accepts `None` as a
  null pointer, and getters and methods then read through it. Setters and
  `delete_*` used to return silently and now raise too, which reverses what
  5.26.1 documented for the `progress` setter: `PDistOptions.progress.fset(None,
  callback)` now raises.
- A `ThreadPool` with one thread runs its work on the calling thread and
  starts no worker. C++ callers see this.
- `PairwiseComparison` gains the virtual `RunsOpenEyeCode()`, which changes
  its C++ ABI. A `TryPDist` or `TryCDist` override must now honour
  `options.num_threads` and call `progress` only on the calling thread.

### Fixed

- **Under an OpenEye memory-pool mode without the Mutexed or ThreadLocal bit,
  parallel `pdist`, `cdist` and the comparison clustering paths (for every
  comparison but fingerprint and descriptor, which run no OpenEye code while
  scoring), descriptor calculation and Murcko extraction ran OpenEye code on
  several threads and could crash.** On every build they now run that work on
  the calling thread.
  Under such a mode a call keeps the GIL while its native work runs, except
  while a progress callback runs Python code, so two oecluster calls never run
  native work at once. That does not stop native code, OpenEye calls included,
  that another thread is already running. The Windows wheel and statically
  linked builds previously ignored a mode set through `openeye.oechem`,
  Murcko extraction included.
- `pdist` and `cdist` build at most one comparison clone per worker that
  runs, not one per requested thread, and an absurd `num_threads` no longer
  throws `length_error`.
- `pdist` and `cdist` refuse a negative `num_threads` with
  `ValueError: num_threads must be non-negative`, and a non-integer one (a
  float, a numeric string, `None`) with `TypeError`. The value goes through
  `operator.index`, as in `k_medoids`. Both checks run before any native call.
  The binding used to refuse such a value only after `pdist` had created its
  storage, and created or resized any `output` file, or after `cdist` had
  built its comparison.
- `pdist` converts its items before it reads any label, so an item that only
  claims to be a molecule is refused by oecluster. `pdist` and `cdist` with
  RMSD conformer expansion verify each item with the conversion's own check
  before the expansion reads it, and refuse an item that check catches with
  oecluster's own error rather than letting it reach OEChem's binding.
- `SuperposeComparison` accepts an empty list.
- `ThreadPool::ParallelFor` computed its chunk count and chunk end in
  arithmetic that overflowed for a `chunk_size` near `SIZE_MAX`, so a `pdist`
  or `cdist` with such a chunk size returned success having computed nothing.
  It now covers the range.
- With fewer than two items there are no pairs. `MMapStorage` used to empty
  an existing file and then fail to map it, on POSIX and Win32 alike (on POSIX
  with `Failed to mmap file ...: Invalid argument`), so
  `pdist(..., output=path)` failed. It now sizes
  the file to 0 bytes, maps nothing, and returns the matrix.
  `consensus(..., output=path)` still refuses fewer than two items before it
  opens the path.
- `SymmetricDistanceMatrix.squareform()` of a matrix with no items returns an
  empty `(0, 0)` array. It returned a 1x1 array, because SciPy reads an empty
  condensed array as one sample. A `pdist` over no items reached this too.

## [5.26.1] - 2026-10-21

### Added

- **`docs/api/swig-classes.rst`, a reference page for the 18 SWIG option and
  storage classes that the `oecluster` package exports**, such as
  `PDistOptions` and `DenseStorage`. The documentation build mocks the compiled
  extension, so autodoc cannot describe these classes. `docs/swig_reference.py`
  writes the page from the built extension, and
  `python -m invoke docs-swig-reference` runs it. The page describes whichever
  `oecluster` the interpreter imports, so rebuild `oecluster_python` before
  regenerating it. `tests/python/test_swig_reference_page.py` fails when the
  committed page differs from what the current build generates. Classes
  reached only through `oecluster.oecluster`, such as `ClusterReportOptions`
  and `RMSDOptions`, are not on the page.
- The wheel workflow has a `lint` job that runs ruff 0.17.0 over `python`,
  `tests` and `benchmarks`. The `publish` job, which runs only for a `v*` tag,
  now requires it, so a release cannot publish while ruff fails.

### Fixed

- **Reading `PDistOptions.progress` or `CDistOptions.progress` now returns
  the callable that was set, or `None` when none is.** It used to return an
  opaque `std::function` object, set or not, that Python could not call: a
  heap copy of the C++ callback that SWIG had no destructor for. Every read
  leaked one, and SWIG printed
  `swig/python detected a memory leak of type 'std::function< void (size_t,size_t) > *', no destructor found.`
  as each returned object was freed. A copy of a set callback also held the
  callable, so a single read kept the callable alive after the options were
  cleared and deleted. Setting the member is unchanged: it takes a callable
  or `None`, refuses anything else with
  `TypeError: Expected callable or None for progress callback`, and a call on
  a `None` self, such as `PDistOptions.progress.fset(None, callback)`, still
  does nothing. `pdist()` and `cdist()` only set the member, so they behave as
  before.
- An allocation failure while the setter copies a progress callback into the
  options now raises `MemoryError: out of memory setting a progress callback`,
  where it would have left the wrapper as an uncaught C++ exception and
  terminated the interpreter. The `MemoryError` that 5.26.0 added for a
  progress callback caught only the conversion, and SWIG applies the
  extension's `%exception` to no member accessor, so nothing caught the
  setter's assignment or the getter's copy. The getter now copies nothing.
  Found by reading the wrapper, not reproduced.
- The `leiden()` docstring's overlap formula renders as written. docutils
  parsed `|N+(i) & N+(j)|` as a reference to an undefined substitution, so
  every documentation build since 5.12.0, which added `leiden()`, reported
  `ERROR: Undefined substitution referenced: "N+(i) & N+(j)".` and marked the
  formula as an error on the Python API page. It is now an inline literal.
- The C++ API page now includes its class hierarchy. Exhale writes the
  hierarchy as an include whose longest line is 17,432 characters, over
  docutils' default `line_length_limit` of 10,000, so docutils rejected the
  include and the page was built without it. A new `docs/docutils.conf`
  raises the limit to 100,000.
- **Re-running CMake on an existing build directory now picks up a SWIG
  upgrade.** FindSWIG trusts its cache: with an executable already cached, it
  asks that executable for its library directory only when `SWIG_DIR` is
  unset or false in CMake's sense, and for its version only when
  `SWIG_VERSION` is, so a build directory went on generating the wrapper from
  the SWIG library it was first configured with, through any number of
  configures. A Homebrew upgrade that removes the old version's directory
  broke the build that way. `cmake/DropStaleSwigCache.cmake` now drops
  `SWIG_DIR` and `SWIG_VERSION` from the cache on every configure, so
  `find_package(SWIG)` derives both again. It drops a cached
  `SWIG_EXECUTABLE` only when it is an absolute path that no longer exists,
  so that SWIG is searched for again. A live executable, including one passed
  with `-DSWIG_EXECUTABLE`, is kept, and so is a bare program name such as
  `swig`. Re-passing `-DSWIG_EXECUTABLE` over a cache that holds an older
  SWIG now works without deleting `CMakeCache.txt`, and the configure error
  for an unsuitable SWIG no longer says to delete it. **`SWIG_DIR` can no
  longer be set in the cache: a value passed with `-DSWIG_DIR` is dropped on
  the very configure that passes it, and FindSWIG's own search sets it
  again.** To point that search at another SWIG library, pass
  `-DCMAKE_INCLUDE_PATH=<dir>`, naming the directory that holds `swig.swg`;
  it stays in the cache. A normal variable named `SWIG_DIR`, set by a
  toolchain file or a parent project, is not touched and still wins. The
  guard runs only when CMake configures, and a plain `cmake --build` after a
  SWIG upgrade does not re-run CMake, so run `cmake <build-dir>` first.

## [5.26.0] - 2026-10-20

### Changed

- **ROCS now reads every conformer of an `OEMol` wherever it sits in a
  molecule list, so a ROCS matrix built from a list that put an `OEGraphMol`
  ahead of an `OEMol` can change under this release. Recompute it rather than
  compare against a stored one.** The bindings exposed two `ROCSComparison`
  constructors: one copied each molecule with every conformer and accepted
  only `OEMol`, the other copied each molecule as its active conformer alone,
  and the list's first element chose between them. A first element that was
  a molecule other than an `OEMol` sent the whole list to the second, so every
  `OEMol` after it was scored on its active conformer, with no error and no
  warning. The second is no longer wrapped, and the first now takes every
  list through one conversion that copies each element by its own type: an
  `OEMol` with every conformer, anything else as the one conformer it has.
  `ROCSComparison`, `pdist()`, `cdist()` and every function that takes
  `comparison="rocs"` build that same native class, so they change together.
  The C++ constructors are unchanged.
- **Lists that were all `OEMol`, or all `OEGraphMol`, are unaffected.** Each
  element reaches the copy it reached before -- `OEMol(const OEMCMolBase&)`
  for an `OEMol`, keeping every conformer, and `OEMol(const OEMolBase&)` for
  anything else -- so only a list that put a molecule other than an `OEMol`
  ahead of an `OEMol` can score differently. A list with an `OEMol` first and
  an `OEGraphMol` later used to raise instead; see below.
- The size of the change, on one fixture: the one behind
  `test_rocs_reads_every_conformer_whatever_the_list_order` in
  `tests/python/test_rocs_conformer_order.py`. It is one Omega conformer of
  `CCCCCC(=O)NCCc1ccccc1` as an `OEMol` (`query`) and as an `OEGraphMol` copy
  (`query_graph`), against thirty Omega conformers of the same SMILES as an
  `OEMol` (`flex`) whose active conformer is set to the one that overlays
  `query_graph` worst. `pdist([query_graph, flex], "rocs", num_threads=1)`
  reads 0.0 as a `combo_norm` distance, and so does the reverse order. The
  active-conformer reading, which is what the previous release returned for
  `[query_graph, flex]`, is 0.604896; the test computes it under this release
  as `active_only`, from an `OEGraphMol` copy of `flex`. Measured on
  2026-10-20 with OpenEye Toolkits 2026.1.0. Omega's conformers can differ
  between toolkit versions, so those numbers belong to this fixture on that
  toolkit, not to the molecule.
- **`cdist()` inherited the order dependence through set A.** It builds one
  comparison over `items_a + items_b`, so set A's first element decided how
  both sets were read: an `OEGraphMol` there reduced every `OEMol` in either
  set to its active conformer, and two sets that were each of one type, but
  not the same one, still interacted. Each element is now read by its own
  type, whichever set holds it. On the fixture above,
  `cdist([query_graph], [flex], "rocs", num_threads=1)` reads 0.0, where the
  active-conformer reading is 0.604896.
- **Mixed lists that used to raise `TypeError` now construct, for ROCS and
  MCS alike.** A first element that was an `OEMol` sent the list to the
  `OEMol`-only constructor, which refused a later `OEGraphMol` with
  `TypeError: List item is not an OEMol object`; through `cdist()` the
  refusal reached across the sets, so `cdist([oemol], [graph_mol], ...)`
  raised it too. `MCSComparison` had the same pair of constructors and takes
  the same conversion now. MCS scores are topological and never read
  coordinates, so for MCS a mixed list changes only in being accepted: benzene
  as an `OEGraphMol` against toluene as an `OEMol` reads 1/7 through `pdist()`
  in either order and through `cdist()` in either direction, as
  `test_mcs_accepts_a_mixed_list_in_either_order` pins. `RMSDComparison`
  still accepts only `OEMol`.
- **Reversing a ROCS `cdist()` asks a different question, and the answer
  depends on direction, not type.** `cdist([flex], [query_graph], "rocs")`,
  refused before, reads 0.604908 on the fixture above at `num_threads=1`, as
  does `cdist([flex], [query], "rocs")`. Entry `[i, j]` of a `cdist()` block is
  what `items_a[i]` finds when it searches `items_b[j]`'s conformers, so a
  multi-conformer `OEMol` in set A is overlaid as its active conformer, here
  the worst one, whatever set B holds. `pdist()` scores both orientations and
  keeps the better, and reads 0.0 for the same pair.
- **The list conversions now report a refused element by its position, its
  type and the type expected, and the message text changed.** Each molecule
  and design-unit list conversion used to refuse an element with
  "List item is not an OEMolBase object",
  "List item is not an OEMol object" or
  "List item is not an OEDesignUnit object", naming no element. Where the bad
  element came first on one of the six comparison constructors, the overload
  dispatcher refused the list before any conversion ran, with
  "Wrong number or type of arguments for overloaded function
  'new_FingerprintComparison'" or its counterpart, naming a private symbol.
  The messages now read
  `molecule list item 2 is a int, not an OEMolBase` for
  `FingerprintComparison`, `DescriptorComparison`, `ROCSComparison`,
  `MCSComparison`, `SuperposeComparison` given molecules,
  `descriptor_statistics()`, `murcko_scaffolds()` and `murcko()`;
  `molecule list item 1 is a OEGraphMol, not an OEMol` for `RMSDComparison`;
  and `design unit list item 1 is a str, not an OEDesignUnit` for
  `SuperposeComparison` given design units. An instance that passes the type
  check but wraps nothing, as one made by `__new__` without `__init__` does,
  gets `molecule list item 1 is a OEGraphMol that wraps no molecule`, or
  `... that wraps no design unit`; `pdist()` reads molecule titles for its
  labels before building the comparison, and on such a molecule that read can
  raise OEChem's own `TypeError` first. The index is the element's position
  in the list the caller passed; `cdist()` names the set and counts within it,
  as in `items_b item 1 is a int, not an OEMolBase`, rather than within the
  `items_a + items_b` it builds over. The type is read from the type's
  `__name__`, or given as `<unnamed type>` when that read raises an
  `Exception` or returns something other than a string; any other exception
  from the read, such as `KeyboardInterrupt` or `SystemExit`, propagates. A
  caller matching on "List item is not" or on the dispatcher's message for
  these inputs has to change.
- `pdist()`, `cdist()` and the other functions that take a comparison name,
  such as `vendi_score()` and `maxmin_select()`, check every item with
  `comparison="rmsd"` before expanding conformers, which they do by default,
  so a refused item gets the message above, indexed before any expansion. An
  item could fail inside the expansion before: an object with a `NumConfs`
  method reporting more than one conformer raised `AttributeError` for its
  missing `GetTitle`, and an `OEMol` made by `__new__` raised OEChem's own
  `TypeError` from `NumConfs`.

### Fixed

- **`murcko_scaffolds()` and `murcko()` could crash the interpreter when
  another thread removed molecules from the list during the call.** Both
  release the GIL while they extract scaffolds, and the list conversion
  borrowed each molecule from the caller's list without taking a reference,
  so emptying the list mid-call freed any molecule nothing else referenced
  while native code was still reading it. Reproduced as a segmentation fault
  inside scaffold extraction at `num_threads=1`:
  `test_clearing_the_list_mid_call_does_not_crash` in
  `tests/python/test_binding_lifetime.py` builds 20,000 `OEGraphMol`s and has
  a second thread call `mols.clear()` after 0.2 s, and its child process died
  on signal 11 under 5.25.1. The conversion now holds a reference to every
  element until the call returns. No result changes. The other molecule and
  design-unit list conversions hold references the same way, but the two
  Murcko functions are the only entry points taking a molecule or design-unit
  list that release the GIL.
- **The list conversions now propagate an exception raised by an element's
  type check, for a check that raises each time it is asked.** The native
  molecule and design-unit list conversions and the `OEFPBatch` conversion
  used to replace such an exception with a `TypeError` about the element's
  type, and now let it reach the caller, as does the item check `pdist()` and
  `cdist()` now run for RMSD. An `AttributeError` from `__class__` does not
  count: Python's `isinstance` treats it as no match, so such an element is
  refused with the positional `TypeError` above. An exception other than
  `AttributeError` raised while reading an object's wrapped pointer now
  propagates too, where the conversion used to discard it. Some Python entry
  points check item types themselves first, for example `pdist()`'s label
  extraction and the public `isim()` and `bitbirch()`, so what a caller sees
  from those depends on that check.
- An allocation failure inside the conversion of a molecule list, a
  design-unit list or a progress callback now raises `MemoryError`, where it
  used to leave the wrapper as an uncaught C++ exception and terminate the
  interpreter. Found by reading the conversions, not reproduced.
- **An object that passes oecluster's Python type check for an OpenEye or
  oefp type, but whose SWIG type says it wraps a C++ object that is not of
  that type, is now refused, except that `ROCSComparison` and
  `MCSComparison` read a molecule as the kind of molecule its SWIG type
  records.** Such an object can be one whose `__class__` misreports its
  type, or one that holds another object's `this`. The molecule list,
  design-unit list and `OEFPBatch` conversions decided an object's C++ type
  with `isinstance` and then used the pointer it wraps as that type,
  unchecked. Each now finds the SWIG type stored with that pointer in the
  SWIG runtime that created it, OpenEye's or the one oecluster shares with
  oefp, and converts the pointer through that runtime's own cast list. The
  `isinstance` checks still run before this and keep their messages. A
  molecule or design unit whose SWIG type does not convert to the type
  required is refused with the positional `TypeError` described under
  Changed, and a batch with a `TypeError` reading
  `Expected an oefp.OEFPBatch.`, the batch conversion's existing message.
  A molecule or design-unit list item whose `this` is not a SwigPyObject at
  all is refused with that positional `TypeError` too, where 5.25.1 read the
  object's memory as a pointer. ROCS and MCS take any molecule, so a
  molecule that claims one kind of molecule but whose SWIG type records
  another is read as the one its SWIG type records: an `OEGraphMol` whose
  `__class__` reports `OEMol` is copied as its one conformer, where
  `RMSDComparison` refuses it. An object that
  passes the type check but whose SWIG type is in no runtime oecluster can
  read raises `RuntimeError`, as
  `molecule list item is a <type> whose SWIG type is in no SWIG runtime oecluster can read (loaded: <runtimes>)`,
  with `design unit` in place of `molecule` for a design-unit list and
  `argument` in place of `molecule list item` for a batch. `<runtimes>` lists
  the `swig_runtime_data` modules imported, such as
  `swig_runtime_data4, swig_runtime_data5`, and the message names no
  position. These conversions also keep the `this` they verified alive for
  as long as they use its pointer, so when `this` is a new object on each
  read, the C++ object it owns is no longer freed while oecluster still uses
  it. These hazards existed in 5.25.1 and earlier. Arguments of oecluster's
  own types still convert through SWIG's standard conversion, which drops
  its reference to `this` before the pointer is used. Some Python entry
  points check or read items themselves first, for example `pdist()`'s
  label extraction and the public `isim()` and `bitbirch()`, so what a
  caller sees from those depends on that step.

### Notes

- The `%exception` registration for `OECluster::cdist` is removed, and that
  changes nothing: SWIG generated no wrapper for that function, and
  `cdist()` calls `cdist_into_address`, whose registration is unchanged. A
  new `tests/python/test_gil_release.py` checks the build rather than the
  interface file. It reads the generated wrapper and requires each of the 36
  registered functions to reach its C++ call with the GIL released, skipping
  with the reason when it cannot tell which build tree the imported module
  came from; and it times `butina_cluster`, `k_medoids_cluster` and
  `murcko_scaffolds` while a second thread tries to run.

## [5.25.1] - 2026-10-17

### Fixed

- **`bitbirch_refine(reassign_top_clusters=k)` could return the same sample in
  two clusters.** The pass pooled the members of the `k` largest leaf
  subclusters, assigned each member to its nearest of the `k` centroids, then
  wrote the centroids that won members into the selected subclusters in order,
  taking one slot per winning centroid. A centroid that won no members took no
  slot, so the selected subclusters past the last slot taken were never
  written and kept whatever members they held before the pass -- any of which
  this same pass had already written into an earlier subcluster. Each centroid
  now writes back to its own subcluster, and a centroid that won no members
  empties its subcluster instead, which `Result()` drops as of 5.25.0.
- **What duplicated a sample was more starved centroids than already-empty
  selected subclusters.** A starved centroid is a selected centroid that is
  the nearest centroid of none of the pooled members. With `s` of them the old
  writer filled only the first `k - s` selected slots, so the last `s` went
  unwritten whichever centroids had starved -- not the starved centroids' own
  subclusters. The selected subclusters are ordered by member count, so any
  already empty sort last: when at least `s` of them were empty the unwritten
  slots held nothing and no sample was duplicated, and otherwise at least one
  unwritten slot still held members this pass had already written into
  another cluster. `redistribute_largest_cluster` is how a selected subcluster
  comes to be empty already. On the four-row fixture behind
  `test_bitbirch_refine_reassign_leaves_an_already_empty_selected_slot_alone`
  one centroid is starved and one selected subcluster is already empty, and
  the result is identical under both writers.
- **Starving a centroid needs neither duplicate fingerprints nor
  near-identical centroids.** Two routes are measured. A tie resolves to the
  lower centroid index, so the later of two equal centroids is starved, and
  duplicate fingerprints are the ordinary way to make two centroids equal;
  that is the route the sweep below takes. A centroid no member prefers is
  starved just as surely, with every fingerprint distinct: over 24 distinct
  32-bit fingerprints (NumPy `default_rng(2)`, density 0.18) at threshold 0.1,
  branching factor 2 and `reassign_top_clusters=5`, the ten-member cluster's
  majority centroid is
  all-zero, so it scores 0 against every fingerprint and wins none of the 23
  pooled members, while the greatest Tanimoto between any two of the five
  selected centroids is 0.0909. **Distinct fingerprints are not on their own
  a reason to think a result was unaffected.**
- **Where the duplication occurred the result loses a cluster -- one per
  unwritten slot that still held members.** Over the 450-point sweep now committed as
  `test_bitbirch_refine_reassign_partitions_duplicate_bearing_inputs` -- 25
  seeds of four 24-bit prototypes repeated three times each, by thresholds
  0.3/0.5/0.7, by branching factors 2/3, by `reassign_top_clusters` 3/4/5 --
  8 points duplicated a sample before this release and none do now. Each of
  the 8 lost exactly one cluster, in every case a trailing one holding a
  member another cluster also held; results from the other 442 points are
  byte-identical to 5.25.0's. On the twelve-fingerprint fixture behind
  `test_bitbirch_refine_reassign_drops_a_starved_cluster_instead_of_duplicating`
  five clusters become four, the duplicate `(9,)` is gone, and `num_clusters`,
  `len(result)`, `cluster_sizes`, `centroids` and the duplicated sample's
  `labels` entry change with it; on the distinct-fingerprint fixture behind
  `test_bitbirch_refine_reassign_starves_a_centroid_with_distinct_fingerprints`
  six become five. `bitbirch()` and `bitbirch_recluster()` never run this
  pass, and neither does `bitbirch_refine()` with `reassign_top_clusters=0`.
- A `select_parameter()` grid point that reached the defect raised
  `RuntimeError: isim_report: sample 9 appears in clusters 1 and 4` and ended
  the sweep; sweeping `reassign_top_clusters` over `[3, 4, 5]` on that fixture
  now scores all three points and picks 5.
- **Reference parity is deliberately given up here, as it was for the emptied
  subcluster in 5.25.0.** The Python reference has the same defect: run on
  2026-10-17 from the sibling checkout the parity tests import, it put sample
  9 in two clusters on the twelve-fingerprint fixture and samples 4 and 12 in
  two on the distinct-fingerprint one, which is what the native pass returned
  before this release. No parity test changed. Every one that reassigns uses
  `top=2`, and none of their fixtures reaches the defect: each still agrees with
  the defective reference under the corrected writer, which it could not if a
  sample were duplicated there. `top=2` itself rules nothing out -- two tied
  centroids starve the second one there as well.

## [5.25.0] - 2026-10-16

### Changed

- **ROCS distances change, and a ROCS matrix computed before 5.25.0 will not
  reproduce under this release -- whichever function produced it.** Two things
  moved. `pdist()` now scores an unordered pair in both orientations and keeps
  the better overlay, where it scored one orientation and stored that as the
  pair's distance; and every search -- `pdist()`'s two, and the one behind
  `cdist()` and a direct `Compare` -- now runs on an `OEOverlay` built for that
  search and discarded after it, where one overlay was built per clone and
  reused. `pdist()` values move for any pair whose two orientations disagree,
  and both measured fixtures show a gap: flexible molecules at twenty
  conformers most, but single-conformer input too -- on the fixture behind
  `test_pdist_rocs_end_to_end` the gap reaches 0.00028 on `combo_norm` and
  0.00056 on `shape`. `cdist()` values move wherever the reused overlay's
  state was reaching the score; the case measured is a cell that read a
  saturated 1.0 where a first-use overlay reads 0.546160012. Recompute rather
  than compare against a stored matrix.
- The reuse was not a rounding-scale effect, and four separate behaviours went
  with it. The orientation bias itself: measured 2026-10-15 over three
  flexible molecules at twenty conformers each, the two orientations of a pair
  differ by up to 0.0355 on `combo_norm`, 0.0710 on `combo`, 0.0660 on `shape`
  and 0.0010 on `color`, each on its own scale and on one fixture rather than
  as a bound.
  A **0.45** order dependence on accepted input: a phenol with one atom
  displaced to 200 angstroms -- a fiftieth of the 10,000 angstrom
  coordinate-extent guard, and accepted -- scored against benzene and octane,
  where three of the six input permutations moved the `pdist()` matrix by
  0.4538. That is the phenol-benzene entry moving between a saturated 1.0 and
  0.546160012. Scoring both orientations did not remove that one; only the
  fresh overlay did. A `cdist()` cell that depended on the thread count and on
  where an unrelated molecule sat in the list: the same phenol-against-benzene
  cell read 1.0 at one thread and 0.546160012 at two, and 0.546160012 again
  with benzene listed first. The phenol-against-octane entry of the `pdist()`
  matrix above now reads **0.772408**, the orientation with octane as the
  overlay reference, where `cdist()` of the same two molecules reports
  0.8297612816095352 in the phenol-as-reference direction it is asked for.
  That 0.772408 is a measurement of one fixture under this release, not a
  property of the input: the fixture is
  `test_rocs_pdist_survives_an_input_that_saturates_a_reused_overlay` in
  `tests/python/test_pdist.py` -- phenol, benzene and octane from SMILES at
  one Omega conformer each, phenol's atom 0 (a ring carbon) moved to x = 200
  angstroms -- read at `combo_norm` as a distance from `pdist(num_threads=1)`.
  Displacing any one of the other twelve atoms instead, and nothing else,
  moves that entry anywhere in 0.747319 to 0.777447 on the same run.
- What holds after the change, on the fixtures that measure it: the `pdist()`
  matrix over the three-molecule flexible fixture agrees with the one computed
  from the reversed input to a residual of 0.0 on all four score types; all
  six permutations of the displaced-phenol fixture agree at 0.0, against
  0.4538 before; and the `cdist()` cell above reads 0.546160012 under all
  three arrangements that used to disagree. The tests assert `1e-6` rather
  than equality, because OEShape promises no bit-reproducibility: these are
  measurements of this toolkit on these fixtures, not a property the library
  can guarantee. Clustering stability is a separate question and this does not
  settle it -- a reproducible matrix is necessary for reordering not to change
  labels, not sufficient, since ties are the algorithm's to break.
- The cost. `pdist()` performs two overlay searches per off-diagonal pair
  instead of one -- n*(n-1) searches over n molecules rather than n*(n-1)/2,
  with the n self-overlays the constructor takes unchanged -- and building an
  overlay per search costs on top of that doubling. `benchmarks/rocs_overlay_reuse.py`
  is committed so that second part can be re-measured: it runs both strategies
  in one process, warmed up, nine repetitions, the two arms in randomised
  order within each repetition. On 2026-10-16, over the eight flexible
  molecules it names, the median of fresh over shared was 1.20 at one
  conformer, 1.09 to 1.11 at twenty and 1.06 at fifty. That is one panel on
  one machine; re-run it rather than carrying the numbers.
- **ROCS results from the comparison forms move too, and land on what
  `pdist()` gives.** Every entry point that consumes a pair as a symmetric
  distance reads the same seam `pdist()` fills through, where each called
  `Compare(min(i, j), max(i, j))` and so disagreed with a `pdist()` matrix of
  the same molecules: `cluster_report()`, `activity_landscape()`,
  `modelability()`, `maxmin_select()`, `circles()`, `knn_graph()` -- and with
  it `jarvis_patrick()` and `leiden()` -- `sphere_exclusion()` away from the
  neighbor order, `vendi_score()` and `logdet_diversity()`. Measured on the
  eight SMILES of `benchmarks/rocs_overlay_reuse.py` at one Omega conformer
  each, `combo_norm` as a distance, `num_threads=1`, Butina at threshold 0.5:
  `cluster_report(result, mols, comparison="rocs").silhouette` was
  **0.15030862789190502** against the `pdist()` matrix form's
  **0.15076074788278154**, and both now read 0.15076074788278154. Only ROCS
  moves; every other comparison leaves the seam at its default, which is the
  `Compare(min, max)` these paths already read.
- `cdist()` and a direct `Compare` keep their direction, deliberately. Entry
  `[i, j]` of a `cdist()` block is what `items_a[i]` finds when it searches
  `items_b[j]`'s conformers, which is the question a cross-distance block is
  asked, and keeping the better of the two orientations there would answer a
  different one. The seam stays a protected member of `PairwiseComparison`
  that defaults to `Compare`, and it does not appear in the bindings: callers
  reach it through one internal accessor, `detail::compare_unordered` in the
  private `src/CompareUnordered.h`, which is the class's only friend. The
  paths that still call `Compare(min, max)` themselves are the ones that
  refuse a ROCS comparison outright -- Butina, DBSCAN, sphere exclusion's
  neighbor order, agglomerative, HDBSCAN and k-medoids -- so no comparison in
  this library can tell the two sources apart there. Lifting one of those
  refusals means routing its path too; `src/clustering/DiversityValidation.h`
  says so beside the refusal.
- The refusals that keep ROCS off the clustering and streaming paths are
  **unchanged in behaviour** -- every one of them still refuses a ROCS
  comparison -- but the reason they give is not. They asserted that a ROCS
  score depends on what its overlay scored before, which three new C++ tests
  now contradict on the fixtures they cover; they now say ROCS is not
  certified for the repeated scoring those paths do. The whole-graph half of
  the repeatability instrument still cannot run for ROCS, because these
  refusals stop it before it starts, so "uncertified" is the accurate word
  rather than "unrepeatable". **The message text changed**: "a ROCS score
  depends on what its overlay scored before" became "ROCS is not certified for
  the repeated scoring these paths do", and the k-medoids and single-linkage
  clauses were reworded with it. A caller matching on "cannot cluster a ROCS
  comparison without a matrix" is unaffected.
- `DenseStorage(n)` and both platform-conditional `MMapStorage(path, n)`
  constructors refuse a sample count whose buffer length cannot be
  represented, where they computed it unchecked and carried on. Both
  multiplications are tested -- `n * (n - 1) / 2`, then
  `pairs * sizeof(double)` -- and `std::overflow_error` names `n` and which
  limit it crossed; through the bindings that arrives as `RuntimeError`
  carrying the same message. The smallest refused counts are 6074001001 for
  the pair count and 2147483649 for the byte count. `MMapStorage` also refuses
  a byte count past this platform's file size limit -- `off_t`'s maximum on
  POSIX, `INT64_MAX` on Windows, first refused at n = 1518500251 -- and every
  refusal happens before the file is opened, so a refused size leaves nothing
  on disk. **`SparseStorage` is not covered and is not fixed**: it allocates
  nothing up front, so its wrapped pair count mis-sizes no buffer, but it
  still reports one -- `SparseStorage(SIZE_MAX, 0.5).NumPairs()` returns 1,
  measured under this release. What the fix replaces:
  `DenseStorage(SIZE_MAX)` used to construct successfully, reporting
  `NumSamples() == 18446744073709551615` with `NumPairs() == 1`, after which
  `Set`, `Get` and every reader that walks `Data()` indexed past a
  one-element allocation.
- `bitbirch_refine()` drops a leaf subcluster its prune pass emptied instead
  of reporting it as an empty cluster. Where that happens, `num_clusters`,
  `len(result)`, `cluster_sizes` and `centroids` are each one shorter -- five
  clusters become four on the four-fingerprint fixture in
  `tests/python/test_bitbirch.py` -- while `labels` is unchanged: member lists
  are ordered by member count, largest first, so an emptied subcluster is
  trailing and no label is renumbered. `cluster_report()` and `isim_report()`
  stop raising `RuntimeError: Cluster must contain at least one member` on
  such a result, and a `select_parameter()` grid point that reached it now
  scores instead of propagating that error. Reference parity is deliberately
  given up here: the Python reference reports the empty cluster, together with
  the zero-sample centroid it computes for it. `bitbirch()` and
  `bitbirch_recluster()` are unaffected -- neither can empty a subcluster -- and
  the empty-cluster refusal itself is unchanged, still firing on a result a
  caller or a foreign clusterer builds.
  `RefinePruneMatchesZeroSampleReferenceCentroid` and
  `ClusterReportRefusesEmptiedRefinementSubclusters` were flipped and renamed
  to `RefinePruneDropsTheSubclusterItEmptied` and
  `BothScorecardsAcceptARefinementThatEmptiedASubcluster`.
- The `Generic` scaffold level stops running the cut-atom normalization whose
  result it discards. `OEUncolorMol` under the `BemisMurcko` strategy
  recomputes formal charge and implicit hydrogen count over every atom, so the
  charge and hydrogen-count transfer running ahead of it was overwritten at
  that level; the transfer and the re-perception that followed it are now
  guarded to `Framework`. **No molecule's scaffold changes and no failure mode
  is fixed.** Over the 937 molecules of `examples/assets/Bender.csv`, the
  `Framework` and `Generic` strings are byte-identical before and after, and
  the three `std::nullopt` exits inside the guarded region were each counted
  at zero at `Generic` on the unguarded build, so no molecule on that panel
  ever reached one. This is work removed rather than output changed. No timing
  was measured, and that panel is the whole of the evidence.

## [5.24.0] - 2026-10-14

### Changed

- A dropped descriptor column now reports one of three reasons where
  `"zero-variance"` was the only one any selection could produce.
  `"zero-variance"` keeps the narrowest of the three: a variance that was
  measured and came out zero or negative, so the column is constant over
  these molecules. `"insufficient-data"` is new, for a column with a value
  in fewer than two of the rows the fit ran over -- no variance was ever
  measured, so nothing is yet established about the column itself.
  `"non-finite-variance"` is new, for a measured variance that is not
  finite. **A caller matching `"zero-variance"` therefore stops matching
  two populations it matched in 5.23.1**: the too-few-rows columns, and the
  non-finite ones, which is where a column with a present NaN lands when it
  has a value in two or more rows -- OpenEye's `FractionCsp3` over a set
  containing water, and all eight of RDKit's `BCUT2D_*` columns over a
  selection containing a bare metal ion, the common real case. The reasons
  are read from `descriptor_statistics()['dropped']` and from
  `DescriptorComparison.DroppedReasons()`; both surfaces now call one shared
  rule, where each carried its own copy of it.
- `"non-finite-variance"` names the variance and claims nothing about the
  values. A present NaN is one way to reach it, but an accumulation over
  perfectly finite values spanning a wide magnitude range can overflow on
  its own, and the variance alone does not separate the two. The guide gives
  a narrowing `descriptor_excluded_indices` call for the dropped column, and
  what its answer does and does not settle. The three tests are ordered --
  row count, then finiteness, then sign -- so a `-inf` variance over two or
  more rows reports `"non-finite-variance"` rather than being read as a
  negative one.
- The refusal both surfaces raise when the selection empties names the
  causes present instead of asserting one. `descriptor_statistics()` and the
  `DescriptorComparison` constructor each said "Every selected descriptor
  column has zero variance over these molecules"; they now say "Every
  selected descriptor column was dropped over these molecules (<causes>)",
  with the comparison's trailing "; there is no distance to compute"
  unchanged. The causes are deduplicated and ordered `insufficient-data`,
  `non-finite-variance`, `zero-variance`, so one unusable selection gives
  one message however the drops were collected. This message is the whole
  account the caller gets: the constructor throws before `DroppedReasons()`
  can be read, and `descriptor_statistics()` returns no drop report when it
  refuses.

### Fixed

- The descriptor section of `docs/python-api.md` presented the
  `"non-numeric"` drop as current behaviour -- "Non-numeric columns are
  dropped from any selection, on every metric". oecluster counts Bool, Int
  and Float columns as numeric, and every descriptor source it can build --
  `openeye`, `mordred` and `rdkit` -- emits only those kinds, so no
  selection can produce that reason. It is now documented as reserved,
  against OEFP adding a source that emits another kind. The same passage
  gained what each live reason asks of the caller, and the header comments
  on `DescriptorStatisticsResult::dropped_reasons` and
  `DescriptorComparison::DroppedReasons`, which described a single
  `"zero-variance"` covering "zero, negative, or not finite" alike, were
  corrected with it.

## [5.23.1] - 2026-10-14

### Fixed

- A test fixture added in 5.23.0 asserted that an ill-conditioned inverse's
  asymmetry exceeded `1e-8`. That quantity is round-off, so it depends on the
  platform's LAPACK -- 4.8e-8 on macOS, 8.1e-9 and 9.2e-9 on the two Linux
  builds -- and the wheel build for 5.23.0 failed on Linux and Windows while
  passing locally and on macOS. The fixture is now bracketed rather than
  pinned: comfortably above pure round-off, and below the band it is proving
  silence over. No library code changed; 5.23.0 was never published.

## [5.23.0] - 2026-10-13

### Changed

Six argument surfaces replace silent acceptance with an error or a warning.
Each took an input, discarded or reinterpreted it, and returned a result
anyway, so a caller who relied on that silence sees new behaviour: a
`TypeError`, a `ValueError` or a `UserWarning` where 5.22.0 returned normally.

- `pdist()` refuses the arguments its prebuilt-comparison branch would
  discard. A prebuilt comparison carries its own items and its own options, so
  any comparison keyword passed beside one now raises `TypeError` naming the
  keywords, and an item count that disagrees with `comparison.Size()` raises
  `ValueError` naming both counts. The keywords are judged first, because no
  edit to `items` makes a dropped option take effect.
- On that same branch, `items` must now be **sized**. Anything without
  `__len__` -- a generator, an iterator, a lazily built view -- raises
  `TypeError` rather than being accepted, because there is no count to check
  the comparison against. This narrows more than the mismatch rule does: a
  call that passed a generator beside a prebuilt comparison ran in 5.22.0 and
  raises in 5.23.0. Pass a list. `similarity=` is unaffected -- it is a named
  parameter rather than a comparison option, and a prebuilt comparison still
  ignores it and reads its own orientation.
- `p` is refused with `TypeError` on any descriptor metric but `minkowski`,
  on every Python entry point that builds a descriptor comparison. That is
  not only `pdist()`, `cdist()` and the `DescriptorComparison` constructor:
  the clustering, selection and diversity entry points route through the same
  descriptor validation, so the rule was measured reaching `butina()`,
  `k_medoids()`, `maxmin_select()`, `circles()`, `sphere_exclusion()`,
  `knn_graph()`, `jarvis_patrick()`, `leiden()`, `vendi_score()`,
  `logdet_diversity()`, `modelability()` and `activity_landscape()`, and the
  clustering entry points `dbscan()`, `hdbscan()`, `agglomerative()` and
  `cluster_report()` as well. The list is given rather than a count, because
  a count drifts the moment another entry point takes descriptor options.
  It was accepted and inert before --
  measured, the condensed values were bit-identical with and without it. The
  rule is case-insensitive and covers the `seuclidean` alias, and C++ still
  speaks first: a metric that is unknown or unavailable on this surface raises
  `RuntimeError` before the `p` rule runs. No precedence against the
  schema-dependent checks is promised, because supplying `variances` or
  `inverse_covariance` resolves the schema earlier and moves them. The rule
  does **not** reach the raw `oecluster.oecluster` bindings, the documented
  escape hatch to the C++ options: `DescriptorOptions::p` is a plain `double`
  with a default, so C++ cannot tell "not supplied" from "supplied as the
  default", and `p` stays inert on that route.
- An asymmetric `inverse_covariance` now produces a `UserWarning` instead of
  nothing. The `mahalanobis` metric has always scored `(A + A.T) / 2`, so
  `[[2, 1], [0, 2]]`, its transpose and the symmetrization returned
  bit-identical distances; the defect was the silence, not the asymmetry. The
  Python layer now substitutes `(A + A.T) / 2` itself and warns, naming one
  offending entry pair with both of its values. It is **not** refused, and the
  C++ validator does not refuse for asymmetry at all: no tolerance separates a
  caller's mistake from the round-off of a real inverse, since a perturbation
  below any band still moves a distance arbitrarily far in relative terms.
  Both the named paths and the constructor warn once per call, and only once
  the comparison has been built -- after every descriptor-option refusal,
  after the input normalization that can empty the item list, and after the
  C++ constructor's own verdicts. It is **not** ordered against validation
  that happens later in the calling entry point: under `-W error` a warning
  about `inverse_covariance` will stop a call that would otherwise have
  failed on, say, `num_threads`, because that filter makes the first warning
  fatal by design. What is guaranteed is that the advisory never displaces a
  refusal about the options it is itself describing. The band is `1e-6` of an entry pair's own
  Cauchy-Schwarz scale, plus a triangular-fill rule for one side exactly zero
  and the other above `1e-6/n` of that scale. Measured over 30 draws per
  decade and 3 to 20 variables, a legitimate inverse stays silent to about
  condition number `1e10` and can warn from `1e11` up, which is acceptable
  because the message reports a substitution rather than an error.
- `SymmetricDistanceMatrix.from_condensed()` requires `probe_triples` to be an
  integer, raising `TypeError` naming the parameter. The check is
  unconditional rather than inside `check=True`, because the downstream guard
  only reached the value at `n >= 3`, so a nonsense `probe_triples` passed
  with a two-item matrix was accepted in silence.
- `murcko()` and `murcko_scaffolds()` decide emptiness with
  `list.__len__` rather than truthiness, matching the `PyList_Size` the SWIG
  typemap reads. An empty `list` subclass whose `__len__` reports otherwise is
  now refused as empty with the documented `ValueError`; before, it reached
  the native call and surfaced a `RuntimeError`.
- `k_medoids()` refuses an empty `initial_medoids` under every `init` with
  `ValueError`, rather than reading it as the option never having been passed.
  Under `init="explicit"` that message comes ahead of the count check. A
  generator `initial_medoids` is still accepted: the emptiness is tested after
  materializing, not with `len()`.

### Fixed

- `missing="propagate"` no longer reports `nan_present` for a matrix that
  contains no non-finite entry. The comparison stamps `data_integrity` from
  its missing-value policy before any pair is scored, which is the only honest
  answer at that point but not after the sweep, so `pdist()` and `cdist()` now
  let the finished matrix have the last word. Only a stamp reading
  `nan_present` is ever rewritten; `subset_scored` answers a different
  question and is untouched. The conservative stamp is kept where the output
  is not a witness to what was computed: sparse `pdist()` storage discards
  values above the cutoff, `+inf` among them, and a `cdist()` run with a
  positive finite cutoff rewrites them to `0.0`. An infinite, negative or
  `NaN` `cdist()` cutoff filters nothing, so those runs are judged on the
  data. A matrix cleared this way is no longer refused on that ground by the
  clustering entry points, which have no override for `nan_present`.
- `cdist()` converts `cutoff` to a C `double` once, at the top, and every
  later reader works from that double rather than from the caller's object. A
  float subclass with an unstable `__gt__` could previously answer one way to
  the `similarity` guard and another to the stamping test while native
  `cdist` filtered with a value neither saw. A non-numeric `cutoff` was
  already a `TypeError` before the matrix was computed in 5.22.0, and still
  is; what moved is how much runs first. The conversion now sits at the top of
  the function, so the refusal precedes building the comparison, materializing
  the inputs and allocating the output, where it used to follow all three. A
  `cutoff` that is not even orderable against a float, such as a string, died
  at the `similarity` guard in 5.22.0 and is unaffected.
- `from_condensed()` names `params=` when it cannot build a mapping from it.
  **Nothing about what `params` accepts or raises has changed.** `dict()`'s
  own exception propagates unaltered -- same class, same args, same object,
  same traceback -- with a note appended through `BaseException.add_note`, and
  only when the exception is a built-in type, so no user code can run on the
  way. The class a caller catches, the `args` they read and the object's
  identity are unchanged, and `params=[]` and `params=""` remain an empty
  mapping. What does change is the rendered traceback and `__notes__`: a
  note-aware matcher -- `pytest.raises(match=...)` from pytest 9 onward
  searches notes -- can match text that only appears in the note, where in
  5.22.0 it could not. Only the diagnostic for
  a case such as `params=0`, which reported `'int' object is not iterable` and
  named neither the parameter nor the entry point, gained context.
- The fingerprint keyword refusal no longer says "select one of those
  families" when it names exactly one. Every current entry names a single
  family, so the plural was always wrong; the sentence is derived from the
  count and stays correct if an entry grows a second family.
- Five docstrings described the `allow_nonmetric` refusal as firing when the
  argument is not a `bool`, but the metric gate accepts `numpy.bool_` as well,
  so each documented the refusal of an input the code takes. They now match
  the wording `cluster_report` already carried.
- The descriptor section of `docs/python-api.md` left two rules unstated:
  `variances` and `inverse_covariance` each belong to one metric and are
  refused on any other, and an empty sequence is refused rather than read as
  an omitted option. It now also carries the `p` rule, the
  `inverse_covariance` symmetrization warning, and -- in the `pdist()`
  section -- the prebuilt-comparison contract.

## [5.22.0] - 2026-10-12

### Added

- `k_medoids()` clusters straight from a comparison, holding no N x N
  matrix. It takes `items` first, in the three forms `butina()` takes: a
  `SymmetricDistanceMatrix`, a prebuilt comparison such as
  `FingerprintComparison(mols)`, or a sequence of items with `comparison=`.
  `distance_matrix=` remains a keyword alias that accepts only a matrix, so
  existing calls behave as before. A comparison form returns exactly what the
  matrix form returns over a matrix filled through the same `Compare`, for
  every `num_threads` and `chunk_size` and for a repeatable comparison; the
  test suite verifies the two forms bitwise on labels, medoids, cost and
  iteration count. ROCS is refused on this path, as on the other comparison
  paths, because its overlay keeps state between calls. Against `pdist()`,
  whose batched kernels agree with `Compare` only to about 1e-12, a cost that
  close to a tie can be decided differently.
- Native overload
  `k_medoids_cluster(PairwiseComparison&, const KMedoidsOptions&)`, and
  `KMedoidsOptions::max_comparisons`, which the storage overload ignores.
- `max_comparisons` bounds the work a streamed run may do, counting one pass
  over every pair as `N * (N - 1)`. The default is 300,000,000,000, which is
  about 2.97 hours of single-threaded comparison *work* at the measured rate
  below -- not of elapsed time, which is that work divided by the parallelism
  the run achieves. The rate is a property of the comparison rather than of
  k-medoids: a SiteHopper `Compare` costs around 0.9 s, where the same count
  is geological, so a caller whose comparison is far dearer or far cheaper
  should set a value. The passes no run can avoid -- the initializer, the
  medoid column fill and one verification pass -- are priced before any pair
  is scored, and a call that cannot afford them is refused with a
  `ValueError` that prices the floor for the `n_clusters` asked for and,
  under `init="build"`, the cheaper `init="farthest_first"` floor beside it.
  Within the budget the swap loop stops between passes once it is spent,
  which reads as `converged=False` with `n_iterations` below
  `max_iterations`. The bound is a brake rather than a hard ceiling: it is
  accurate to within an additive `O(n_clusters * N)`, the per-pass row
  bookkeeping each initializer does being left unpriced, so a run can finish
  slightly over it. `max_comparisons` is refused with a matrix, which leaves
  no pair to score.
- `benchmarks/kmedoids_streaming.py`, which runs the two paths against each
  other, one subprocess per cell, and records the core count, the resolved
  worker count and the 1-minute load average beside every row.

### Changed

- What streaming k-medoids costs is passes, not exactness. A converged run
  makes `init passes + (iterations + 1) + n_clusters` passes over every pair,
  each of `N * (N - 1)` comparisons, plus `N * n_clusters` for the one-off
  medoid column fill, where the matrix form reads each of those pairs from
  storage instead. `init="farthest_first"` is the lever: it spends one
  initializing pass against BUILD's `n_clusters`, making the floors
  `n_clusters + 1` passes against `2 * n_clusters`, so switching saves under
  half however large `n_clusters` grows -- the `n_clusters` verification
  passes that make `converged` a guarantee are what neither initializer can
  shed. Streaming is therefore slower in time and far smaller in memory than
  the matrix path; it is not a speedup.
- Measured with `benchmarks/kmedoids_streaming.py` on an Apple M3 Max
  (14 CPUs, 36 GiB of memory) over seeded library molecules at
  `n_clusters=20`, one subprocess per cell, each figure a single measurement.
  Process peak RSS, which does not depend on machine load: at 10,000 items
  the streamed path held 198 MB against the matrix path's 651 MB, a 3.28x
  reduction; at 5,000 items 155 MB against 267 MB, 1.72x; at 2,000 items
  125 MB against 141 MB, 1.13x. Both paths returned the same cost, the same
  iteration count and the same medoids in all six pairs. That ratio is a
  process ratio and understates the difference: both paths carry about
  110 MB of interpreter, oecluster and OpenEye plus the cost of building the
  molecules (about 88 MB at 10,000 items) before clustering starts, and what
  the clustering call itself adds is 0.0 to 4.6 MB streamed against 2.4 to
  52.7 MB from a matrix. In place of the matrix the streamed path keeps a
  medoid column cache of `N * n_clusters * 8` bytes, 1.6 MB at 10,000 items
  and `n_clusters=20` against the condensed matrix's 400 MB. Both paths now
  build that cache, so the matrix form is untouched except for an
  `N * n_clusters` table it holds beside its matrix -- 16 MB at 100,000 items
  with `n_clusters=20`, 0.04% of the condensed matrix there -- which is what
  keeps the assignment build one code path.
- A single-threaded in-kernel `Compare` on a 2048-bit fingerprint measured
  0.0356 microseconds, flat to within 2% over a 4x range of item count. A
  comparison count multiplied by that rate is **comparison work, not elapsed
  time**: at 45,000 items with `n_clusters=20`, the `init="build"` floor of
  8.1e10 comparisons is about 48 minutes of comparison work, which the
  default chunking is predicted to run in about 4.4 minutes of wall time,
  and the `init="farthest_first"` floor of 4.2e10 is about 25 minutes of
  work, predicted at about 2.3 minutes. Those predictions are a model, not a
  measurement. **No elapsed time for a clustering run is published here as a
  validated measurement**: every such timing was taken at a 1-minute load
  average of 12.1 to 30.2 against the benchmark's validity limit of 7. The
  per-comparison rate above is a measured rate rather than a duration -- it
  is what the default budget rests on -- and it was taken in the same window.
  Within that window a streamed run took 11.3x to 24.2x the matrix path's
  clustering time, narrowing as the item count grows, which is indicative
  rather than a validated timing.
- `chunk_size` is the tuning lever for a streamed run. The default
  `chunk_size=4096` divides a pass into `ceil(N / 4096)` work units, and the
  binding quantity is the largest unit rather than their number, because a
  short tail unit finishes early while one worker carries a full 4096: the
  parallelism a pass achieves is about `N / max(unit size)`, 1.00x at 2,000
  items, 2.44x at 10,000, 10.99x at 45,000, saturating at the core count of a
  14-core machine at 57,344 items. Setting `chunk_size = ceil(N / num_threads)`
  uses the machine; at 2,000 items, where the default runs a pass serially, a
  14-unit chunk measured 6.4x the default's single-threaded per-comparison
  rate, while at 10,000 items the default already delivers 2.44x and the same
  tuning buys roughly 2.6x more. This chunking is unchanged and shared with
  the matrix path, which has always behaved this way; it costs more here
  because a streamed pass is bound by `Compare` calls where a matrix pass
  reads memory.
- A streamed run checks each distance as it scores it and raises
  `RuntimeError` naming the two items when it is not finite, where the matrix
  form relies on the metric gate refusing a non-finite matrix with a
  `ValueError` before the native call. This is consistent with what 5.20.0
  established for the streamed MST paths: a path with no stored matrix for the
  gate to scan checks what it scores. 5.21.0 extended the same named-pair
  error to complete, average and weighted linkage at the C++ layer, where
  those paths do hold a matrix and the Python gate still refuses a non-finite
  one first, so a Python caller never reaches it there.

## [5.21.0] - 2026-10-10

### Changed

- Complete, average and weighted linkage are built from one condensed
  matrix over the live clusters, in place of the `2N-1` node distance table
  and the all-pairs heap. Peak memory falls from about 14x the input matrix
  to about 2x for the algorithm itself, and from about 15x to about 2.4x
  end to end through Python, with the results unchanged for a finite,
  non-negative input. Measured in C++ alone on a Red Hat Enterprise Linux
  8.10 server (128 cores, 1.5 TB of memory) at `num_threads=16`, which is
  not what a default call resolves to, with both paths in one binary: at
  20,000 items process peak RSS falls from 22,402 MB (14.00x the matrix)
  to 3,205 MB (2.00x), and at 10,000 items from 5,603 MB (14.01x) to
  804 MB (2.01x). What the algorithm itself allocates is 1.00x the input
  matrix against the heap's 13.00x. End to end through Python on an Apple
  M3 Max (14 CPUs, 36 GiB of memory), where the interpreter, numpy and a
  different allocator add to the total, 10,000 items fall from
  6,137-6,142 MB (15.35-15.36x) to 947-950 MB (2.37x). Peak RSS does not
  depend on machine load. Wall time at 20,000 items on that server falls
  from 316 / 271 / 248 s to 5.57 / 5.34 / 5.12 s for complete / average /
  weighted; both sides ran at 16 threads and the old path's merge loop is
  serial where the new one's is not, so part of that ratio is parallelism
  rather than algorithm. One shape regresses, and the cost there is cubic:
  a synthetic all-equal matrix measured 1.98x to 5.31x slower than 5.20.0
  from 500 to 4,000 items, while a duplicate-heavy library stayed 2.6x to
  11x faster. Those degenerate timings are indicative only, several having
  been taken at a 1-minute load above the benchmark's validity limit of 7.
  What governs the cost is how often a merge invalidates the cached nearest
  neighbours: each row remembers its nearest cluster and must rescan when
  that cluster is merged away, so when the cluster being built stays the
  nearest neighbour of most rows, nearly every row rescans at every step.
  No simple property of the distances decides this. Equal distances can
  produce it, as the all-equal fixture does, but concentration alone does
  not: a few tight, well-separated groups share one inter-group distance
  across most of the matrix and stay fast, because each row's nearest stays
  inside its own group. A geometry with entirely distinct distances in
  which one growing cluster remains everyone's nearest produces it without
  any repeated value at all. Deduplicating identical items does not rule
  the shape out either -- the all-equal fixture holds no identical items,
  only mutually equidistant ones. Measure before committing to a large run.
- A non-finite distance in the input matrix raises `std::runtime_error`
  naming the two items from complete, average and weighted linkage, where
  5.20.0 produced an unspecified result. Single linkage has behaved this
  way since 5.20.0. This is a C++ API change only: Python callers are
  unaffected, because the metric gate refuses a non-finite matrix with a
  `ValueError` before the native call, as it already did.

### Fixed

- A `chunk_size` large enough to overflow the chunk-count arithmetic no
  longer skips the initial copy. In 5.20.0 `chunk_size = SIZE_MAX` returned
  a tree built from an all-infinity table.

## [5.20.0] - 2026-10-08

### Added

- `hdbscan()` and `agglomerative(..., linkage="single")` cluster straight
  from a comparison, holding no N x N matrix. Both take `items` first, in the
  three forms `butina()` takes; `distance_matrix=` remains a keyword alias
  that accepts only a matrix, so existing calls behave as before. Single
  linkage compares every pair once; HDBSCAN compares every pair once for its
  core distances and at most once more for the tree. The result equals the
  matrix form's over a matrix filled through the same `Compare`, for a
  repeatable comparison; ROCS is refused on these paths.
- Native overloads `hdbscan_cluster(PairwiseComparison&, const HDBSCANOptions&)`
  and `agglomerative_cluster(PairwiseComparison&, const AgglomerativeOptions&)`.
- `benchmarks/streaming_mst.py`. Measured on an Apple M3 Max (14 CPUs, 36 GiB of memory), each path in its own process,
  with the 1-minute load average recorded before each row (6.4 to 7.4 for
  5.20.0 unless stated, 41.9 to 44.1 for the 5.19.0 baseline). Matrix single
  linkage's process peak RSS fell from 5.19.0 to 5.20.0: 4.88 to 0.65 GB at
  10,000 molecules with one thread, 5.85 to 0.65 GB at 10,000 with the
  default threads, 8.91 to 2.09 GB at 20,000 with one thread, and 12.11 to
  2.10 GB at 20,000 with the default threads. Single linkage from a
  comparison at 20,000 molecules held 0.29 to 0.30 GB where the matrix path
  held 2.10 GB. The clustering step of matrix single linkage measured 17.7x
  to 37.7x faster than 5.19.0 on that contended machine; the load difference
  between the two runs inflates this, so treat it as an upper bound (a
  core-share correction, an estimate and not a measurement, suggests 6x to
  11x). At 100,000 molecules, HDBSCAN and single linkage from a comparison
  peaked at 1.04 to 1.07 GB of process RSS, where the condensed matrix alone
  would need 40.0 GB and is never allocated; those four runs were taken at
  loads of 7.26 to 66.90, so they carry no timing claim.

### Changed

- Single linkage, from a matrix too, is built from the minimum spanning tree.
  Merge heights and `distance_threshold` cuts are unchanged, but merges at a
  tied height come in the tree's order: children and cluster-size order
  within a tie level can differ from 5.19.0, and an `n_clusters` cut that
  falls inside a tie level can return a different, equally valid partition.
  Matrix single linkage no longer allocates the heap algorithm's cluster
  distances and heap.
- An allocation failure in `hdbscan()` or `agglomerative()` raises
  `MemoryError` on every path, where it raised `RuntimeError`.
- `hdbscan(chunk_size=)` defaults to 64 and sets the matrix core-distance
  work unit, where it was unused.
- A NaN or infinite distance read from a comparison raises `RuntimeError`
  naming the two items in `hdbscan()` and single-linkage `agglomerative()`,
  and a negative distance in `hdbscan()`, from a matrix or a comparison,
  raises `RuntimeError` naming the two items; a distance whose quotient by
  `alpha` is not finite raises `RuntimeError` naming `alpha`; and
  `hdbscan(alpha=nan)` raises `ValueError`. 5.19.0 returned undefined results
  for the negative and alpha cases. A non-finite distance in a matrix passed
  to `hdbscan()` or `agglomerative()` is unchanged: the metric gate refuses it
  with a `ValueError` before the native call. The native matrix overloads
  `hdbscan_cluster` and single-linkage `agglomerative_cluster` have no such
  gate and now raise `std::runtime_error` naming the pair.
- A zero merge distance, or a zero core distance inside HDBSCAN, is +0.0
  even when the input held -0.0.
- `oecluster algorithms` lists `items` as the input of `hdbscan` and
  `agglomerative`. Their `--set` options are unchanged.

## [5.19.0] - 2026-10-05

### Added

- `butina()`, `dbscan()` and `sphere_exclusion(..., order="neighbors")`
  cluster straight from a comparison, holding no N x N matrix. `butina()` and
  `dbscan()` take `items` first, in the three forms `cluster_report()` takes:
  a `SymmetricDistanceMatrix`, a prebuilt comparison, or a sequence of items
  with `comparison=`. `distance_matrix=` remains a keyword alias that accepts
  only a matrix, so existing calls need no change. The comparison forms build
  the threshold graph by comparing every pair twice, once to count each
  item's neighbors and once to record them, and hold no pairwise structure
  but that graph: 16 bytes per within-threshold pair plus 16 per item on a
  64-bit platform, where a condensed matrix costs 8 bytes for every pair.
  `sphere_exclusion()`'s nearest assignment then compares each non-center
  item with every center, as it always has. The result equals the matrix
  form's over a matrix filled through the same `Compare`.
- An exact memory guard for those graphs. The first pass gives the graph's
  size before it is allocated, and a graph above the limit raises
  `MemoryError` naming the item count, the edge count and the bytes needed.
  The default limit is the larger of the condensed matrix and 1 GiB;
  `max_graph_bytes=` sets an explicit budget, which large runs should do.
- ROCS is refused on those paths, named or prebuilt: a ROCS score depends on
  what its overlay scored before, so the graph's two passes can disagree.
  Cluster a ROCS matrix from `pdist()` instead. Every other built-in
  comparison family passed the new repeatability tests, which run one or two
  configurations each.
- Native overloads `butina_cluster(PairwiseComparison&, const ButinaOptions&)`
  and `dbscan_cluster(PairwiseComparison&, const DBSCANOptions&)`, and the
  `Neighbors` order on `sphere_exclusion`'s comparison overload.
  `ButinaOptions`, `DBSCANOptions` and `SphereExclusionOptions` gain
  `max_graph_bytes`.
- `benchmarks/streaming_threshold.py`, which measures both paths in separate
  processes. On an Apple M3 Max (14 CPUs, 36 GiB of memory), Butina over
  20,000 molecules at threshold 0.4519 took 1.54 s with a process peak RSS of
  293 MB from a comparison, and 0.68 s with a process peak RSS of 2.10 GB
  through `pdist()` and the 1.6 GB matrix it builds: about 2.3 times as long,
  for about a seventh of the peak memory. At 100,000 molecules Butina from a
  comparison took 51.8 s with a process peak RSS of 1.08 GB, and did not
  build the 40 GB matrix the matrix path would.

### Changed

- An allocation failure in `butina()`, `dbscan()` or `sphere_exclusion()` now
  raises `MemoryError` on every path, the matrix path included, where it
  raised `RuntimeError`.
- `oecluster algorithms` lists `items` as the input of `butina` and `dbscan`.
  Their `--set` options are unchanged.

## [5.18.1] - 2026-10-05

### Fixed

- Wheels build again. 5.18.0's wheel tests failed on every platform, so it
  never reached PyPI: the command-line tests read rendered terminal text, and
  that text changed shape under CI. rich-click forces a colour terminal when
  `GITHUB_ACTIONS` is set, and the escape codes split an expected message
  around a styled option name; rich also wraps at 80 columns, so beside
  Windows' long temporary paths an expected message split across two lines.
  The tests now pin plain, unwrapped rendering. The command line itself is
  unchanged.

## [5.18.0] - 2026-10-04

### Added

- An `oecluster` command line over a precomputed distance matrix, with
  `cluster`, `select-parameter`, `stability` and `consensus` commands plus
  `algorithms` for discovery. `cluster`, `select-parameter`, `stability` and
  `consensus` in bootstrap mode take `--algorithm NAME` with repeatable
  `--set key=value`; `consensus` in cross-algorithm mode takes one
  `--member 'NAME;key=value'` per ensemble member and refuses `--algorithm`
  and `--set`; `algorithms` takes neither, only an optional algorithm name.
  Options are validated against the roster before any file is read, so an
  unknown one is reported with a suggestion rather than as a `TypeError` from
  deep in the library. The workflow options are checked before the read too
  — `--criterion`, `--min-clusters`, `--max-clusters`,
  `--max-noise-fraction`, `--resamples`, `--fraction`, `--seed` and
  `--threshold` — by calling the library's own validators, so a bad bound
  never costs a file read. `--threads` and `--allow-nonmetric` are
  first-class flags and route to every place the underlying call needs them.
- The command line reads `.npz` matrices and `oepdist`'s `.npy` and `.bin`
  output through its JSON sidecar. A file holding similarities rather than
  distances is refused, one whose orientation neither the sidecar nor the
  archive records warns on stderr and proceeds, matching the library's own
  rule, and one the sidecar proves holds distances is reported as such. A
  `.npy` on the sidecar path must be one-dimensional: `oepdist` writes a
  condensed matrix as 1-D and a cross-distance matrix as 2-D, and a 2-D
  file beside a mismatched sidecar would otherwise cluster as nonsense.
  `.csv` is not accepted: `oepdist` writes titles unquoted and values at
  eight significant digits, and records no provenance.
- `--output` writes results as CSV or JSON. The JSON carries
  `schema_version` 1; NaN encodes as `null` and the infinities as `"inf"`
  and `"-inf"`, since collapsing them together would erase a real
  distinction. Every destination is checked before the run: one that would
  clobber the input, its sidecar or another destination is refused, decided
  by filesystem identity rather than by comparing names.
- `--quiet` suppresses the terminal summary, `--traceback` (or
  `OECLUSTER_CLI_TRACEBACK`) shows the original exception instead of a
  one-line message, and the exit code is 0 for success, 1 for a failed run
  and 2 for a bad invocation. The unproven-orientation warning goes to
  stderr unconditionally and `--quiet` does not silence it.

## [5.17.0] - 2026-10-03

### Added

- `consensus()` combines an ensemble of partitions into one. It accumulates
  how often each pair of items shared a cluster, divides by how often the
  pair was observed together, and extracts a partition from the resulting
  co-association distances: by default the connected components of the pairs
  supported by at least `threshold` (0.5) of the members that saw them, or
  any matrix-consuming `ClusteringSpec` passed as `method=`. The ensemble is
  a `ClusterStability`, a `ParameterSelection`, or a sequence of
  `ClusteringResult` and `(positions, labels)` members, so partial partitions
  are normalized per pair. The `ConsensusResult` is itself a
  `ClusteringResult` and carries the co-association matrix, one
  `ConsensusRecord` per cluster with Monti's cluster consensus, the per-item
  consensus, and one adjusted Rand index per member.
- The matrix is a `DenseStorage` by default and an `MMapStorage` when
  `output=` names a file, the `pdist` convention, so a consensus over tens of
  thousands of items can live on disk.

## [5.16.0] - 2026-10-03

### Added

- `cluster_stability()` reruns a `ClusteringSpec` on repeated subsamples of a
  `SymmetricDistanceMatrix` or an `oefp.OEFPBatch` (drawn without
  replacement, half the items by default) and matches every reference
  cluster to its best counterpart in each resample by Jaccard overlap. The
  read-only `ClusterStability` carries one `ClusterStabilityRecord` per
  cluster (mean best Jaccard, dissolved and recovered fractions at Hennig's
  0.5 and 0.75 thresholds, resamples evaluated), the `(K, R)` Jaccard
  matrix, one adjusted Rand index per resample through
  `partition_agreement()`, and optionally every resample's positions and
  labels. The Jaccard matching never treats noise as a cluster; `noise`
  governs the agreement call only. A callable that returns a partition of the
  wrong size, hands back the reference result object again, relabels the
  reference while resampling, or returns labels that are not one-dimensional
  is refused with `ValueError`.
- `take()` builds the sub-matrix or sub-batch over selected item positions,
  in the given order, through the new native `take_pairs` (dense,
  memory-mapped and sparse storage, parallel over rows for the dense kinds)
  and `take_fingerprints`. The comparison facts carry over; a permutation
  of every item inherits every fact; a proper subset of a dense or
  memory-mapped source re-runs the triangle probe when the source's had
  run (a sparse subset inherits the probe fields), and a source stamped
  NaN-present is re-measured on the subset.
- Descriptor distance matrices record the `missing` policy in `params`, so
  `take()` can tell an ignore-scored matrix from a propagated one after
  excluding the NaN items.

## [5.15.2] - 2026-10-03

### Fixed

- Wheels build and import again. oefp 0.3.0 (since 5.0.0) compiles Arrow into
  the extension, which no tag since v4.2.2 had exercised: the RHEL 8 build now
  uses `gcc-toolset-12` for pyarrow's C++20 headers (the wheel stays
  `manylinux_2_28`), and the Linux and macOS wheels load Arrow and Parquet out
  of the installed `pyarrow` package instead of trying to vendor it, with the
  needed sonames validated against that pyarrow at build time.
- `pyarrow>=25,<26` is declared as a runtime dependency: the extension links
  libarrow/libparquet by versioned soname, the same ABI coupling `oefp`
  carries.
- The `oepdist` command puts pyarrow's library directories on the loader path
  before starting the binary.

## [5.15.1] - 2026-10-02

### Fixed

- `select_parameter()` refuses a stage flag in `report_options`
  (`compute_pair_rank_indices`, `compute_centroid_indices`) that is not a
  `bool` with `TypeError` before the first clustering, instead of letting
  the scorer refuse it after a run.
- `ClusteringSpec` equality short-circuits on identity like Python's
  container equality, so a spec whose option is NaN (a float or an array)
  is equal to itself.
- The `ClusteringSpec` docstring scopes "not validated at construction" to
  the options, as the documentation already did.

## [5.15.0] - 2026-10-02

### Added

- `ClusteringSpec` describes a clustering entry point (a roster name, its
  function, or any callable returning a `ClusteringResult`) plus the
  keyword options it is run with, and runs it on any input with
  `spec.run(items, **overrides)`.
- `select_parameter()` sweeps one keyword of a spec over an explicit grid,
  scores every partition with `cluster_report()` (a `SymmetricDistanceMatrix`
  or a prebuilt comparison) or `isim_report()` (an `oefp.OEFPBatch`), and
  returns a read-only `ParameterSelection` of `SweepRow` records with the
  winner under a named validity index. Optional `max_noise_fraction`,
  `min_clusters` and `max_clusters` bounds exclude rows before ranking; NaN
  scores never win; ties keep grid order. The opt-in scorer stage a
  criterion needs is switched on automatically. Pure Python; no native
  change.

## [5.14.0] - 2026-10-02

### Added

- `isim()` returns the iSIM Tanimoto similarity of an `oefp.OEFPBatch`
  (summed pairwise intersections over summed unions) in linear time, in C++
  and Python.
- `isim_report()` scores a clustering from binary fingerprints without a
  distance matrix. The core is linear in the number of fingerprints: the
  exact size profile, iSIM intra- and inter-cluster distances, per-cluster
  iSIM separation, and radius, mean medoid distance and Calinski-Harabasz
  exact relative to iSIM medoids. `compute_centroid_indices=True` adds an
  iSIM silhouette, iSIM nearest cluster, medoid Davies-Bouldin and Dunn, and
  coverage at O(N K) cost. Integer accumulation makes results independent of
  the thread count. Tanimoto only; `metric=` is reserved.
- `isim()` and `isim_report()` raise `MemoryError` in Python when the native
  call fails with `std::bad_alloc` or `std::length_error`, matching
  `cluster_report()`. Other native refusals raise `RuntimeError`.

## [5.13.0] - 2026-10-01

### Added

- `cluster_report`, `activity_landscape` and `modelability` score a
  comparison directly, holding no distance matrix, in C++ and Python.
  - The first argument (after `result` for `cluster_report`) is now `items`:
    a `SymmetricDistanceMatrix`, a prebuilt comparison such as
    `FingerprintComparison`, or a sequence of items with `comparison=` and
    its options as keyword arguments. `distance_matrix=` is a keyword alias
    accepting only a matrix; passing both is a `TypeError`.
  - Results are exact: identical to the matrix form over a matrix filled
    through `Compare`. Individual fingerprint distances agree with `pdist()`
    to about 1e-12, so a thresholded count or an argmin-derived field
    (medoids, representatives, nearest cluster, CH/DB, RMODI, MODI) can
    differ when inputs sit within that margin of a tie or threshold. Memory
    is O(N) per worker; time stays O(N^2).
  - New `chunk_size` option (default 4096 distances per work unit) on all
    three, validated on every Python path (the C++ storage overloads ignore
    it) and used only on the comparison paths.
  - Comparisons are gated on their declared facts as the matrix gate gates
    a matrix; `cluster_report` honours `allow_nonmetric`. Pair-rank indices
    are matrix-only, and items that normalizing drops or expands are refused.
  - `modelability` makes `N * (N - 1)` comparisons, twice
    `activity_landscape`'s, because its nearest-neighbour search covers each
    sample's whole row; with a single scored class it only validates each
    pair once.
  - The C++ overloads take `PairwiseComparison&` and follow the
    `maxmin_select` / `knn_graph` comparison contract; SWIG dispatches on the
    argument type, so Python reaches them through the existing functions.

### Changed

- A zero median in `cluster_report` is always `+0.0`, on the matrix path as
  well. An input holding `-0.0` could previously report `-0.0`.
- The comparison paths check only the pairs they score. A non-finite value
  on a pair no metric reads is refused by the matrix gate but never seen on
  a comparison path; no result depends on it.

## [5.12.0] - 2026-09-30

### Added

- `leiden`, Leiden community detection (Traag, Waltman and van Eck, 2019)
  on a shared-nearest-neighbor graph, in C++ and Python. It takes a
  `KNNGraph` or any input `knn_graph` accepts.
  - Each kNN edge {i, j} is weighted by the Jaccard overlap of the two
    neighborhoods (each item's `k` neighbors plus itself), and edges below
    `prune` (default 1/15) are dropped. Seurat's `k.param` is this `k + 1`.
  - `objective` is `"modularity"` (with `resolution`, default 1.0) or
    `"cpm"` (Constant Potts Model). `n_iterations=-1` iterates until a pass
    changes nothing.
  - Every returned cluster is connected. The same `seed`, inputs and build
    give the same result, whatever `num_threads` is.
  - `LeidenResult` reports `quality`, `iterations`, `objective`,
    `resolution` and `k`; `.method` is `"leiden"`.

## [5.11.1] - 2026-09-30

### Fixed

- Python functions that return a result object (every `*Result`, `KNNGraph`,
  `ClusterReport`, `PartitionAgreement` and the rest) no longer abort the
  interpreter when memory runs out while the result is handed to Python. That
  step used to happen outside the binding's exception handler; it now raises
  `MemoryError`. The result is also moved into its Python wrapper instead of
  copied, so a large result no longer needs twice its memory for a moment.

## [5.11.0] - 2026-09-30

### Added

- `knn_graph`, a public k-nearest-neighbor graph. Each row lists an item's
  `k` nearest other items and their distances, ordered by (distance,
  position), with ties going to the lower position. The result is a
  `KNNGraph` with `indices`, `distances` and `positions` arrays.
  - It takes a `SymmetricDistanceMatrix` (dense, memory-mapped or sparse), a
    prebuilt comparison, or a list of items with a comparison name. The lazy
    paths keep only the O(N·k) graph, plus one comparison clone per running
    worker during the build.
  - Sparse matrices must hold every pair within their cutoff, as `pdist()`
    writes them; the builder cannot verify that, and a missing nearer pair
    silently changes a row. An item with fewer than `k` stored neighbors is
    refused.
- `jarvis_patrick`, classic Jarvis-Patrick clustering. Two items link when
  each is among the other's `k` nearest neighbors and they share at least
  `kmin` of them. Clusters are the connected components. It takes a
  `KNNGraph` or any input `knn_graph` accepts, and returns
  `JarvisPatrickResult` with `k` and `kmin`.
- The C++ entry points of the same names, configured through
  `KNNGraphOptions` and `JarvisPatrickOptions`.

## [5.10.0] - 2026-09-29

### Added

- `sphere_exclusion`, sphere-exclusion clustering with a pluggable seed order.
  Each center claims every unclaimed item at or within the threshold, and the
  result lists one center per cluster.
  - `order="input"` is leader clustering. `order="neighbors"` is
    Taylor-Butina, and with the default first-claim assignment equals
    `butina()` with the same `threshold` and `reordering`. A sequence of
    positions is Directed Sphere Exclusion (DISE), such as an `argsort` by
    distance to a reference compound.
  - `assignment="nearest"` keeps the centers and moves each other item to its
    nearest center, with ties going to the earlier center.
  - It takes a precomputed `SymmetricDistanceMatrix`, a prebuilt comparison, or
    a list of items with a comparison name, as `maxmin_select` does. The input
    and permutation orders run lazily in O(N·k) comparisons. The neighbor
    order needs a matrix.
  - The C++ entry point of the same name is configured through
    `SphereExclusionOptions` and returns `SphereExclusionResult`.

### Fixed

- `butina()`, `butina_cluster()`, `dbscan()`, and `dbscan_cluster()` no longer
  skip every pair when `chunk_size` is near `SIZE_MAX`, which made Butina
  return every item as a singleton and DBSCAN label every item noise. The
  threshold-graph build wrapped its chunk count to zero; it now clamps the
  chunk size to the pair count.

### Changed

- `butina_cluster()` runs on the shared sphere-exclusion engine. Its outputs
  are unchanged, and its candidate loop no longer erases from the front of a
  vector, which removes a quadratic cost on large inputs.

## [5.9.0] - 2026-09-29

### Added

- `vendi_score`, the Vendi score of Friedman and Dieng (TMLR 2023): the
  effective number of distinct items in a set, from the similarity kernel its
  distances define. The C++ entry point of the same name is configured
  through `VendiOptions`.
  - `order=1` is exact, from the kernel's eigenvalues. It drops non-positive
    eigenvalues without renormalizing, as the reference implementation does,
    and reports `min_eigenvalue` and `negative_mass` so a kernel that is not
    positive semidefinite can be judged. It is capped at `max_exact` items
    (default 2048).
  - `order=2` is `n^2 / ||K||_F^2`. It needs no spectrum, runs at any size,
    and on a comparison holds only O(N) memory.
- `logdet_diversity`, `log det(K + ridge I)` from the same spectrum. It is a
  positive-definite log-determinant: a singular or indefinite ridged kernel
  scores `-inf`, and `nonpositive_count` and `min_eigenvalue` tell the two
  apart. The C++ entry point of the same name is configured through
  `LogDetOptions`.
- Both build the kernel as `1 - d` (`kernel="complement"`, distances in
  [0, 1]) or `exp(-d / bandwidth)` (`kernel="laplacian"`). They take a
  precomputed `SymmetricDistanceMatrix`, a prebuilt comparison, or a list of
  items with a comparison name, as `maxmin_select` does. Results are
  bit-identical for every `num_threads` and `chunk_size`.
- The exact scores use an in-tree symmetric eigenvalue solver, so the project
  gains no linear-algebra dependency.

## [5.8.0] - 2026-09-28

### Added

- `maxmin_select`, farthest-first (MaxMin) selection of a diverse subset, with
  the C++ entry point of the same name configured through `MaxMinOptions`.
  Each pick is the unselected item farthest from the selection so far, ties to
  the smaller index, so results never depend on `num_threads`. It stops at a
  `count`, at a distance `threshold`, or at whichever comes first, starts from
  an item, the item farthest from the first item that survived normalization,
  or (matrix input only) the global medoid, and can extend an existing
  `initial` selection. The result reports the picks, each pick's distance to
  the earlier selection, and why it stopped.
- `circles`, the #Circles coverage measure of Xie et al. (ICLR 2023): the size
  of a packing whose members are pairwise strictly farther apart than a
  threshold. `method="maxmin"` (the default) packs farthest-first;
  `method="sequential"` is the paper's reference greedy pass over input order,
  without its shuffle. Either count is a lower bound on the packing number.
- Both take a precomputed `SymmetricDistanceMatrix`, a prebuilt comparison, or
  a list of items with a comparison name, as `pdist` does. The last two are
  evaluated lazily, comparing only the pairs the run needs, so a library too
  large for an O(N^2) matrix can still be subset. Returned positions always
  refer to the caller's items, including when normalization dropped some,
  which are reported in `excluded`. A comparison whose distances cannot be
  ranked -- a similarity, `missing='propagate'` or `missing='ignore'` -- is
  refused before any pair is scored.

### Changed

- k-medoids' FarthestFirst initialization now runs on the farthest-first
  kernel the new entry points share. Its output is identical, including on
  matrices holding NaN, which it still does not validate; a regression test
  recorded against the previous kernel pins that.

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
  `max_matches` bounds how many matches one directed search enumerates. A bool
  is refused: `bool` is an `int` subclass, so `True` would otherwise be read as
  a budget of 1, the most destructive value the option has, and would change
  scores silently rather than raise.
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
  rather than assumed. Four clones are taken; their snapshot addresses are
  required to be pairwise disjoint and disjoint from the parent's, and each
  pointer is dereferenced to confirm it reaches a molecule with the bond count
  it should. This needed a test-only accessor, because the property has no
  consequence any score can show: an aliasing clone returns identical numbers,
  keeps the same molecules alive and reports the same `Size()`. Disjointness
  *between the clones* is the operative half -- `pdist` and `cdist` give each
  worker its own clone and never dereference the parent inside the parallel
  region, so a `Clone()` that deep-copied once and then handed every later
  caller the same snapshots would put one molecule set under every thread.
  Nothing in the suite had distinguished any of this: with `Clone()`
  neutralised to alias the parent's snapshots, all 685 other C++ tests still
  passed, so the documented thread-safety guarantee had rested on reading the
  code.
- A `tsan` CMake preset builds the C++ tests under ThreadSanitizer into their
  own `build-tsan/` tree, with the procedure documented in
  `docs/developer.md`. It exercises the concurrency this project owns -- the
  thread pool, the storage backends and the progress callback -- and reported
  no races across the 101 distinct tests the two documented filters select.
  The OpenEye libraries are prebuilt and uninstrumented, so it can say nothing
  about the toolkit's own internals, and it does not verify the `Clone()` deep
  copy above; the test does that.

### Fixed

- `pdist` now raises on `cutoff > 0` together with `similarity=True` instead of
  silently discarding the highest scores. Sparse storage zeroes values *above*
  the cutoff, which drops the far pairs of a distance matrix but the near pairs
  of a similarity matrix. The result was not an obviously broken matrix of
  zeros -- everything below the cutoff came back untouched, so the call looked
  plausible with exactly the most-similar pairs replaced by `0.0`. **Calls with
  a named comparison that would have used sparse storage -- no `output=` -- now
  raise `ValueError`**, matching the guard `cdist` has always had. This affected
  every similarity-capable comparison, `fingerprint` included, not only the
  `mcs` comparison new in this release.
- `pdist` now also raises on `cutoff > 0` with no `output=` for a prebuilt
  comparison object that reports similarity values, such as
  `MCSComparison(mols, similarity=True)`. The `similarity` argument says
  nothing about a prebuilt object, so the guard above never saw these calls
  and the same corruption went through: the benzene-toluene MCS score of
  0.857 came back as `0.0`. The orientation is read from the object's own
  facts before any pair is scored. An object that reports no orientation is
  still accepted, and the Python `cdist` needs no change because it refuses
  prebuilt objects outright.
- The native `pdist` and `cdist` drivers now refuse the same combination with
  `ComparisonError`, which the Python bindings raise as `RuntimeError`. The
  guards above live in the Python wrappers, so C++ callers and callers of the
  raw SWIG bindings still reached the corruption. Native `pdist` refuses a
  comparison whose facts report similarities when the storage is a
  `SparseStorage`, whatever that storage's cutoff; native `cdist` refuses one
  when `CDistOptions::cutoff > 0`. Both refuse before any pair is scored, and
  a comparison that reports no orientation is still accepted.
  `PDistOptions::cutoff` is not read by the native driver at all -- the
  storage applies its own cutoff -- and is now documented that way.

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
  still raises. `cdist` builds one comparison over `items_a + items_b`, so the
  first element of `items_a` decides for both sets. `mcs` follows the same
  rule for which lists are accepted, but never reads coordinates, so its
  scores are unchanged. See the `rocs` and `mcs` sections of
  `docs/python-api.md`.

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
  `BothScorecardsAcceptARefinementThatEmptiedASubcluster` in
  `tests/cpp/test_bitbirch_clustering.cpp` pinned that gap under its former
  name. See 5.25.0, which closed it.
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
