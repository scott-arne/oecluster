# MCS Comparison Measurements

The `mcs` comparison's documentation cites specific measured numbers: a
16.2-second exhaustive search, 66 inclusion-exclusion violations, eleven
discriminating pairs out of 118, 10.4 KB per molecule per clone. This page is
where those numbers come from.

It is a record of measurements taken during development, not a benchmark suite.
Nothing here re-runs automatically, and the figures are not tracked over time.
For runnable performance tools see [Benchmarks](benchmarks.md).

## How to read this page

Each entry gives the claim as it appears in the shipped documentation, where it
is cited, what was actually measured, and a **provenance** verdict:

- **Pinned** — the exact inputs live in a tracked test file, so anyone can
  re-derive the figure today.
- **Described** — the procedure is recorded but the input set was not
  preserved. A re-run would measure something similar, not the same thing.
- **Lost** — the inputs were named only by common name. The figure stands as a
  historical observation and cannot be checked against anything.

Most of the headline figures are **Described** or **Lost**. That is the honest
state of this evidence, and stating it is the reason this page exists.

## Environment

All figures come from probes against OpenEye Toolkits 2026.1.0 on a single
Apple Silicon machine, single-threaded unless the entry says otherwise.
Molecules were parsed from SMILES and hydrogen-suppressed with
`OESuppressHydrogens(mol, false, false, false)`. Except where an entry says
otherwise, searches used the `OEMCSMaxBondsCompleteCycles(1.0)` ranking functor
with `max_matches = 1024`.

Timings are wall-clock on one machine and were never repeated across hardware.
Treat every ratio as sound and every absolute millisecond figure as indicative.

## Pinned inputs

These SMILES are tracked, and every **Pinned** entry below uses them.
`tests/cpp/test_mcs_comparison.cpp` holds every molecule in this table;
`tests/python/test_mcs.py` holds all but benzene-d1, methane and testosterone.
Where both files define a molecule the strings are byte-identical.

| name | SMILES |
| --- | --- |
| benzene | `c1ccccc1` |
| toluene | `Cc1ccccc1` |
| cyclohexane | `C1CCCCC1` |
| benzene-d1 | `[2H]c1ccccc1` |
| methane | `C` |
| morphine | `CN1CC[C@]23c4c5ccc(O)c4O[C@H]2[C@@H](O)C=C[C@H]3[C@H]1C5` |
| penicillin G | `CC1(C)S[C@@H]2[C@H](NC(=O)Cc3ccccc3)C(=O)N2[C@H]1C(=O)O` |

Testosterone, sucrose and the macrolide fragment are also pinned; their SMILES
are long and are best read from the test files directly.

## Search mode: approximate against exhaustive

**Claim.** "On a 53-bond against 54-bond macrolide pair, it took 16.2 s and
matched 50 bonds where approximate took 8.5 ms and matched 51."

**Cited in.** `MCSComparison.h` (the `MCSSearchMode` enum documentation),
`CHANGELOG.md` under 5.7.0, `docs/python-api.md`, and the
`MCSComparison.__new__` docstring in `python/oecluster/__init__.py`.

**Measured.** Erythromycin against azithromycin, 53 and 54 bonds. Exhaustive
took 16,206 ms, printed `Warning: MCS search truncated`, and returned 50
matched bonds. Approximate took 8.5 ms and returned 51. Exhaustive was about
1,900 times slower and returned the smaller match.

**Provenance: Lost.** The two molecules were recorded by name only. No SMILES
for either appears anywhere in this repository.

The supporting sweep -- 15 drug-like pairs of 28 to 44 bonds, exhaustive 1 to
5 ms against approximate 0.1 to 1.4 ms, with identical bond counts on all 15 --
is **Lost** for the same reason. So are sucrose against raffinose (39 ms
against 0.7 ms) and a C30 against C28 alkane (4.6 ms against 1.0 ms, 27 bonds
both ways).

## Search mode: how often exhaustive actually wins

**Claim.** "In a 120-pair scan it won on eleven of the 118 pairs that
completed, by one to three bonds."

**Cited in.** `docs/python-api.md`.

**Measured.** Sixteen molecules ranging from benzene to a 56-bond vancomycin
fragment, all 120 pairs, symmetrized bond counts in both modes. Two pairs
exceeded a 4-second ceiling and were abandoned, leaving 118 completed. Eleven
of those discriminated, and exhaustive returned the larger match on every one.
Its cost on those eleven ran from 6.4 ms to 1,839.1 ms; approximate was 10 to
100 times faster on every pair.

The eleven, exhaustive against approximate bond counts:

| pair | exhaustive | approximate | exhaustive cost |
| --- | --- | --- | --- |
| sucrose / macrolide fragment | 17 | 15 | 6.4 ms |
| raffinose / macrolide fragment | 17 | 15 | 10.0 ms |
| sucrose / paclitaxel fragment | 16 | 13 | 11.3 ms |
| morphine / testosterone | 9 | 8 | 19.1 ms |
| raffinose / paclitaxel fragment | 17 | 16 | 19.5 ms |
| morphine / digoxin fragment | 11 | 10 | 32.6 ms |
| testosterone / macrolide fragment | 19 | 18 | 95.9 ms |
| macrolide / digoxin fragment | 21 | 19 | 370.2 ms |
| macrolide / paclitaxel fragment | 22 | 21 | 411.7 ms |
| testosterone / paclitaxel fragment | 19 | 18 | 857.9 ms |
| cholesterol / testosterone | 20 | 19 | 1,839.1 ms |

**Provenance: Lost** for the scan as a whole: the sixteen-molecule set was not
recorded, so neither the 120-pair nor the 118-pair denominator can be
reconstructed. Five of the sixteen do survive as pinned fixtures -- benzene,
morphine, sucrose, testosterone and the macrolide fragment -- while raffinose,
cholesterol and the paclitaxel, digoxin and vancomycin fragments do not.

Three rows are better off than the verdict suggests, in decreasing order.
Sucrose against the macrolide fragment, the first row, is asserted in both
suites and re-derived on every run: `tests/cpp/test_mcs_comparison.cpp` and
`tests/python/test_mcs.py` each pin 0.605263 approximate against 0.527778
exhaustive, which is 15 matched bonds against 17, and the C++ comment quotes
the 6.4 ms. Morphine against testosterone has its directed approximate counts
of 8 and 7 asserted in the C++ suite, whose maximum is the 8 in this table,
but nothing asserts the exhaustive 9. Testosterone against the macrolide
fragment is built from pinned fixtures as well, so its inputs survive, but no
test pairs those two and neither of its bond counts is asserted anywhere.

One caveat travels with the scan and was recorded at the time: one of the
sixteen SMILES, a vancomycin fragment, parsed with an `Unclosed ring` warning
and so was not the structure it was meant to be. It still yielded a well-formed
56-bond molecule and the scan treated it as another input. It appears in none
of the eleven discriminating pairs and in no test fixture, so the table above
and the pinned pairs are unaffected. It is counted in the 120 pairs, in the
118 that completed, and in the 84 that change their score at `max_matches = 1`
under "Ranking functor and `max_matches`" below.

## Metric properties

**Claim.** "No violation appeared in 74,400 ordered triples, but the
inclusion-exclusion bound that would prove the Jaccard metric property was
violated 66 times over 59,280 triples."

**Cited in.** `MCSComparison.h` (class documentation), `CHANGELOG.md` under
5.7.0, `docs/python-api.md`. It is the evidence behind `triangle = Unknown` on
the gate facts.

**Measured.** Three molecule sets, every ordered triple in each:

| set | molecules | ordered triples | triangle violations | inclusion-exclusion violations |
| --- | --- | --- | --- | --- |
| drug-like | 12 | 1,320 | 0 | not measured |
| fragment-like, transitivity-stressing | 25 | 13,800 | 0 | not measured |
| mixed | 40 | 59,280 | 0 | 66 |

The 74,400 in the shipped claim is the sum of the three triple counts. The 66
violations come from the mixed set alone; the inclusion-exclusion bound was not
evaluated on the other two. Self-distances were exactly zero for all 12
molecules of the first set, and the thinnest triangle margin observed anywhere
was 1.0000 against 1.1444.

**Provenance: Described.** The set sizes and their selection intent are
recorded; the members are not. A re-run on fresh sets of the same sizes would
be a new experiment, not a check of this one. The shipped reason for reporting
`triangle` as `unknown` is the one stated in `MCSComparison.h`,
`CHANGELOG.md` and `docs/python-api.md`: the inclusion-exclusion bound was
violated, so the proof is unavailable rather than merely unattempted. The thin
1.0000 against 1.1444 margin is a further reason not to lean on the empirical
result -- a different 40 molecules could plausibly find a triangle violation
where these did not -- but it is a supporting observation, not the stated
ground.

## Symmetry and reproducibility

**Claim.** Approximate search is asymmetric often enough that `Compare` must
search both directions, which is why the class documentation says so.

**Measured.** Over 45 pairs times 3 repeats, approximate search was
deterministic on all 45 -- 0 nondeterministic -- but asymmetric on 8. Over a
separate 190-pair drug-like scan, 21 were asymmetric.

**Provenance: Described** for both sweeps; the pair sets were not preserved.

**Pinned** for the three pairs that matter. A 20-molecule drug-like scan over
all 190 pairs found 14 asymmetric pairs and searched every 3-subset of them for
a cycle in the winner relation. It found exactly one, and those three pairs are
pinned as C++ fixtures precisely because an acyclic set of asymmetric pairs is
consistent with some total order and so cannot refute a ranking-based
implementation:

| pair | bonds | A as pattern | B as pattern | max | other |
| --- | --- | --- | --- | --- | --- |
| testosterone / morphine | 24 / 25 | 8 | 7 | 0.804878 | 0.833333 |
| morphine / penicillin G | 25 / 25 | 11 | 10 | 0.717949 | 0.750000 |
| penicillin G / testosterone | 25 / 24 | 5 | 4 | 0.886364 | 0.911111 |

The winner relation is cyclic: testosterone beats morphine beats penicillin G
beats testosterone. These figures are asserted in
`tests/cpp/test_mcs_comparison.cpp`, so the suite re-derives them on every run.

## Ranking functor and `max_matches`

**Claim.** "Quality saturated by 256 on all 45 pairs measured", and "small
values are destructive rather than merely faster: at 1, 22 of those 45 pairs
scored wrong."

**Cited in.** The `MCSOptions::max_matches` documentation in
`MCSComparison.h`, and the `MCSComparison.__new__` docstring.

**Measured.** Symmetrized bond counts over 45 pairs, first varying the ranking
functor and then the match budget. The source records "45 pairs" for this
sweep and for the symmetry sweep above without stating that they are the same
45, so do not read the two as one set.

| functor | total bonds | pairs differing from CC(1.0) |
| --- | --- | --- |
| `OEMCSMaxBonds` | 293 | 1 |
| `CompleteCycles(1.0)` | 292 | — |
| `CompleteCycles(1.5)` | 287 | 2 |
| `CompleteCycles(2.0)` | 281 | 4 |

| `max_matches` | total bonds | pairs differing from saturated |
| --- | --- | --- |
| 1 | 147 | 22 of 45 |
| 16 | 290 | 2 |
| 256 | 292 | 0 |
| 1024 | 292 | 0 |
| 4096 | 292 | 0 |

The default of 1024 is the toolkit's own, kept for its fourfold headroom over
the observed saturation point rather than because 1024 was measured to be
necessary.

A second scan agrees on the direction: of the 120 pairs in the search-mode
sweep above, 84 change their score at `max_matches = 1` against the saturated
value.

**Provenance: Described** for the sweeps. **Pinned** for the destructive
effect: morphine against penicillin G scores 0.717949 at the default and
0.958333 at `max_matches = 1`, both from pinned fixtures, and
`tests/python/test_mcs.py` asserts both.

## Construction cost and search reuse

**Claim.** The cost table in `docs/python-api.md` uses 0.287 ms per search
construction.

**Measured.** Twenty targets against one pattern:

| | time |
| --- | --- |
| construction only | 5.74 ms, 0.287 ms each |
| fresh search per target | 13.82 ms |
| one reused search | 8.90 ms |

Construction is 41.5% of the fresh-per-target cost on that fast workload and
about 3% of the 8.5 ms erythromycin pair. Reuse produced identical bond counts,
and results were independent of target order. Breaking out of the match
iterator early does not avoid a slow search: 15,935 ms against 16,679 ms on the
erythromycin pair.

The shipped 0.691 ms mean is the 13.82 ms fresh-per-target row divided by the
twenty targets. That mean includes the 0.287 ms construction, so the
search-only share is about 0.404 ms. The cost table for 1,000 molecules runs
999,000 directed searches, two per pair: its `~400 s` search-time row is
0.404 ms times 999,000 searches, its `~287 s` construction row is 0.287 ms
times 999,000 searches, and the two sum to the `~11.5 min` single-threaded
total.

**Provenance: Lost.** Neither the pattern nor the twenty targets were recorded.

## Per-clone molecule copies

**Claim.** "Memory grows as `num_threads x n x 10.4 KB`", and the cost table's
"peak memory added by per-clone molecule copies, ~79 MB".

**Cited in.** `MCSComparison.h` class documentation and `docs/python-api.md`.

**Measured.** Copy-constructing 1,000 drug-like molecules moved peak RSS by
10,128 KiB. Copying those 1,000 molecules took 5.51 ms, so building 8 clones
costs about 44 ms -- against the roughly 90 s of search time the cost table
projects for that same input.

A note on the units, because the division does not work otherwise: 10,128 KiB
across 1,000 molecules is 10,371 bytes each, which is 10.4 kB in decimal
kilobytes and 10.13 KiB in binary ones. The shipped "10.4 KB per molecule" is
the decimal reading. The "~79 MB" for 8 clones is the binary one --
10,128 KiB x 8 = 79.1 MiB. The two shipped figures are consistent with each
other only under that mixed reading, which is worth knowing before anyone tries
to reconcile them.

Correctness was checked alongside cost: running a directed search against a
private copy rather than the original returned the same bond count on all three
pairs tried -- morphine/penicillin G 11, sucrose/macrolide 15,
cholesterol/digoxin 24.

**Provenance: Described.** The 1,000-molecule set was not preserved. The per-
molecule figure is an average over drug-like input and should not be applied to
peptides or other large molecules.

## Match-level observability

**Claim.** The three match-level presets are distinguishable, and ring
membership is the surprising one.

**Measured.** Over 153 pairs of small molecules, 132 distinguish at least two
levels: loose differs from default on 85, default differs from exact on 81.

**Provenance: Described** for the sweep. **Pinned** for four worked examples,
two given as distances and two as similarities:

Distances: benzene against cyclohexane gives 0.000 under loose and 1.000 under
default; benzene against toluene gives 0.143 under default and 0.556 under
exact. Both pairs are asserted in the test suites.

Similarities: ethane against ethene (`CC` against `C=C`) gives 1.0 under loose
and 0.0 under default, asserted in `tests/python/test_mcs.py` as the
measurement behind "`default` already constrains bond order." Cyclohexane
against hexane (`C1CCCCC1` against `CCCCCC`) gives 0.833333 under default and
0.0 under exact, asserted in the same file as the measurement behind "ring
membership is the addition that surprises." Both pairs are pinned fixtures.

**Lost** for the one isotope claim: 13C-butane against butane scoring 1.0 at
every level, cited in `docs/python-api.md` and `MCSComparison.h`. Neither
molecule appears in any test and neither SMILES is written down anywhere, so
this is the weakest-evidenced of the match-level claims the shipped docs make.

## Hydrogen suppression and isotopes

**Claim.** Isotopic hydrogens are suppressed along with the rest, so a
deuterated analogue scores as identical to its parent.

**Measured.** `OESuppressHydrogens` under SDK defaults, against the same call
with `retainIsotope=false`, as (atoms, bonds):

| SMILES | parsed | default suppression | `retainIsotope=false` |
| --- | --- | --- | --- |
| `c1ccccc1` | (6, 6) | (6, 6) | (6, 6) |
| `[2H]c1ccccc1` | (7, 7) | (7, 7) | (6, 6) |
| `[2H]C([2H])([2H])c1ccccc1` | (10, 10) | (10, 10) | (7, 7) |
| `Cc1ccccc1` | (7, 7) | (7, 7) | (7, 7) |

This is why the constructor passes `retainIsotope=false` explicitly rather than
taking the SDK default. Under the default, benzene-d1 carries 7 bonds against
benzene's 6 and the pair scores `6/(6+7-6) = 0.857` instead of 1.0 -- and
benzene-d1 becomes indistinguishable from toluene by bond count alone.

**Provenance: Pinned** for three of the four rows. `c1ccccc1`,
`[2H]c1ccccc1` and `Cc1ccccc1` are all fixtures in
`tests/cpp/test_mcs_comparison.cpp`, which asserts benzene against benzene-d1
at a distance of exactly 0. The trideuterated row is a fixture nowhere, but
its SMILES is written out above, so it is the one entry on this page that is
re-derivable from the page itself rather than from the suite.

## Zero-bond molecules

**Claim.** Methane, water and argon are refused at construction because bond
Tanimoto has a zero denominator for them.

**Measured.** Each has 1 atom and 0 bonds after hydrogen suppression. Methane
against water yields no MCS matches, so the denominator is `0 + 0 - 0`.

**Provenance: Pinned** for the refusal, which both suites assert -- the C++
suite through methane, and the Python suite through the hydride cases
`[H][Li]` and `[H][H]`, each expecting a `RuntimeError` naming the missing
bonds. **Lost** for water, argon, and the methane-against-water calculation --
all three were probed once and recorded by name only.

## Toolkit behaviour this design depends on

Not performance figures, but direct probes against 2026.1.0 that the
implementation relies on. These are probes of the SDK itself rather than
measurements of the library, so none of them has its inputs in a tracked test
file and **Pinned** is not available for any of them. Each bullet gives its own
verdict.

- `OEMCSSearch` copies its pattern. After constructing from a 7-atom molecule
  and calling `Clear()` on the source, `GetPattern().NumAtoms()` still reads 7
  and `Match` still returns 7 bonds. **Provenance: Lost** -- the molecule is
  not named anywhere, not even by common name.

- `umatch` does not change the score. Over all 120 pairs of the sixteen
  molecules in the search-mode scan, symmetrized bond counts with
  `umatch=false` and `umatch=true` were identical on all 120. **Provenance:
  Lost** -- it ran over the unrecorded sixteen-molecule set from the
  search-mode sweep above.

- `Match` left the target's observable state unchanged across five repeated
  matches against morphine, comparing canonical SMILES and every atom's and
  bond's ring and aromatic flags. This is consistent with thread-safety but is
  not a proof of race-freedom, and the design does not treat it as one.
  **Provenance: Described** -- morphine is a pinned fixture in both suites, so
  the input survives, but no test performs the comparison this bullet describes.

- `SetMCSFunc` and `SetMaxMatches(1024)` both return `true`.
  **`SetMaxMatches(0)` also returns `true`**, and `GetMaxMatches()` then reads
  0 -- the SDK accepts a zero budget silently, which is why the constructor
  refuses it. **Provenance: Described** -- `tests/cpp/test_mcs_comparison.cpp`
  calls the first two at lines 140-141 without asserting their return values,
  and `SetMaxMatches(0)` and `GetMaxMatches()` appear in no test at all.

- There is no `OEFindRingAtomAndBond` in this SDK's OEChem, in the C++ headers
  or the Python module; `OEAssignAromaticFlags` and `OEPerceiveSymmetry` exist.
  This is why the snapshot is warmed with a read traversal rather than an eager
  perception call. This is a symbol-existence fact about the installed SDK,
  re-checkable at any time by grepping the headers and the Python module, so
  the three-verdict vocabulary does not apply.

## What this record cannot support

Stated plainly, so nobody has to infer it:

- The erythromycin/azithromycin pair behind the most widely cited figure in the
  MCS documentation cannot be reconstructed. Neither molecule's SMILES was
  recorded.
- The 66 inclusion-exclusion violations, the saturation point of 256, the
  0.287 ms construction cost and the 10.4 kB per molecule all rest on input
  sets that were not preserved.
- Every timing is from one machine and one toolkit build. Ratios should
  transfer; absolute milliseconds should not be quoted as portable.
- The 8-thread row of the cost table assumes linear scaling, which was never
  measured at all.

The figures are reported here as they were observed. Anyone who needs one of
them to be load-bearing for a new decision should re-measure it on inputs they
record.
