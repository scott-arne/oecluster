# SuperposeMetric Redesign: OESuperpose Integration

## Summary

Replace `SuperposeMetric` and `SiteHopperMetric` with a single `SuperposeMetric`
backed by the `oespruce::OESuperpose` API. Support all `OESuperposeMethod` values
(except SiteSequence and Undefined), expose score type selection, a global
similarity flag, and oeselect-based atom predicates. Add the `similarity` flag to
FingerprintMetric and ROCSMetric as well.

## Motivation

The current implementation uses low-level OEBio APIs (`OEGetAlignment` + `OERMSD`)
and duplicates logic between `SuperposeMetric` and `SiteHopperMetric`. The
`OESuperpose` API from oespruce provides a high-level interface that supports
multiple superposition methods (Global, DDM, Weighted, SSE, SiteHopper) with
built-in atom predicate support.

## Decisions

- **Single class**: `SiteHopperMetric` is deleted. `"sitehopper"` becomes an alias
  that routes to `SuperposeMetric` with `method=SiteHopper`.
- **Skip SiteSequence**: The `SiteSequence` method requires per-structure site
  residues, which is inherently asymmetric and incompatible with pairwise distance
  matrices.
- **Score semantics per method**: Methods return either RMSD or a similarity score.
  No method returns both non-sentinel RMSD and Tanimoto simultaneously.
  - RMSD group (Global, GlobalCarbonAlpha, DDM, Weighted): RMSD
  - Tanimoto group (SSE): Tanimoto [0,1]
  - PatchScore group (SiteHopper): Patch Score [0,4]
- **Distance/similarity handled per-metric**: Each metric has a `similarity` option
  in its options struct. No base class changes.
- **RMSD has no similarity form**: `similarity=true` with RMSD methods returns raw
  RMSD unchanged (no-op). No error is raised.
- **oeselect for predicates**: Predicate strings compiled once at construction,
  passed through to `OESuperpose::SetupRef` and `OESuperpose::Superpose`.

## C++ API

### Enums

```cpp
enum class SuperposeMethod {
    GlobalCarbonAlpha,  // Default
    Global,
    DDM,
    Weighted,
    SSE,
    SiteHopper
};

enum class SuperposeScoreType {
    Auto,        // Selects natural score for the method
    RMSD,
    Tanimoto,
    PatchScore
};
```

### SuperposeOptions

```cpp
struct SuperposeOptions {
    SuperposeMethod method = SuperposeMethod::GlobalCarbonAlpha;
    SuperposeScoreType score_type = SuperposeScoreType::Auto;
    bool similarity = false;
    std::string predicate;      // oeselect expression for both ref and fit
    std::string ref_predicate;  // Overrides predicate for ref
    std::string fit_predicate;  // Overrides predicate for fit
};
```

### SuperposeMetric

```cpp
class SuperposeMetric : public DistanceMetric {
public:
    using Options = SuperposeOptions;

    explicit SuperposeMetric(
        const std::vector<std::shared_ptr<OEBio::OEDesignUnit>>& dus,
        const Options& opts = Options());

    explicit SuperposeMetric(
        const std::vector<OEChem::OEMolBase*>& mols,
        const Options& opts = Options());

    double Distance(size_t i, size_t j) override;
    std::unique_ptr<DistanceMetric> Clone() const override;
    size_t Size() const override;
    std::string Name() const override;
};
```

`Name()` returns `"superpose:<method>"`, e.g. `"superpose:global_carbon_alpha"`,
`"superpose:sitehopper"`. This appears in output metadata and `DistanceMatrix`.

### Error Handling

**Construction time:**
- If a predicate string fails to compile via oeselect, throw `MetricError` with
  the invalid expression.
- If structures are empty or have no atoms, throw `MetricError`.

**Distance() call time:**
- If `SetupRef()` returns false, throw `MetricError` (e.g. no atoms match the
  predicate for structure `i`).
- If `Superpose()` returns false, throw `MetricError` (e.g. no common atoms
  between ref and fit).
- The `OESuperpose` API returns sentinel values (-1.0 for RMSD/Tanimoto) when it
  cannot compute a result. If the extracted score is a sentinel, throw
  `MetricError` rather than returning a meaningless value.

### Distance() Behavior

1. Get thread-local `OESuperpose` instance (created during `Clone()`)
2. Call `SetupRef(structure[i], ref_pred)`
3. Call `Superpose(results, structure[j], fit_pred)`
4. Extract score per resolved `score_type`
5. Apply distance or similarity conversion

### Score Resolution Table

| Method            | Auto Score | Raw Range | Distance         | Similarity        |
|-------------------|-----------|-----------|------------------|-------------------|
| GlobalCarbonAlpha | RMSD      | [0, +inf) | raw RMSD         | raw RMSD (no-op)  |
| Global            | RMSD      | [0, +inf) | raw RMSD         | raw RMSD (no-op)  |
| DDM               | RMSD      | [0, +inf) | raw RMSD         | raw RMSD (no-op)  |
| Weighted          | RMSD      | [0, +inf) | raw RMSD         | raw RMSD (no-op)  |
| SSE               | Tanimoto  | [0, 1]   | 1.0 - tanimoto   | raw tanimoto      |
| SiteHopper        | PatchScore| [0, 4]   | 4.0 - patch_score| raw patch_score   |

Requesting an incompatible score type (e.g., RMSD from SSE) throws `MetricError`.

### Predicate Resolution

- If `ref_predicate` is non-empty, compile it via oeselect for ref; else if
  `predicate` is non-empty, use that; else `OEIsTrueAtom()`.
- Same logic for `fit_predicate`.
- Predicate strings compiled once during construction, stored in `SharedData`,
  shared across clones.

## Cross-Metric Similarity Flag

Add `bool similarity = false` to `FingerprintOptions` and `ROCSOptions`.

### FingerprintMetric

| Similarity Function | Distance             | Similarity         |
|---------------------|----------------------|--------------------|
| Tanimoto            | 1.0 - sim            | raw sim [0,1]      |
| Dice                | 1.0 - sim            | raw sim [0,1]      |
| Cosine              | 1.0 - sim            | raw sim [0,1]      |
| Manhattan           | raw value            | error              |
| Euclidean           | raw value            | error              |

### ROCSMetric

| Score Type | Distance               | Similarity           |
|------------|------------------------|----------------------|
| ComboNorm  | 1.0 - combo/2.0       | combo/2.0 [0,1]      |
| Combo      | 2.0 - combo            | raw combo [0,2]      |
| Shape      | 1.0 - tanimoto         | raw tanimoto [0,1]   |
| Color      | 1.0 - tanimoto         | raw tanimoto [0,1]   |

## CLI Changes (`oepdist`)

### `superpose` subcommand (replaces both `superpose` and `sitehopper`)

```
oepdist superpose structures.oedu -o dist.npy \
  --method global_carbon_alpha \
  --score-type auto \
  --sim \
  --predicate "protein and backbone" \
  --ref-predicate "chain A" \
  --fit-predicate "chain A"
```

Options:
- `--method`: `global_carbon_alpha` (default), `global`, `ddm`, `weighted`,
  `sse`, `sitehopper`
- `--score-type`: `auto` (default), `rmsd`, `tanimoto`, `patch_score`
- `--sim`: return similarity instead of distance
- `--predicate`: oeselect expression for both ref and fit
- `--ref-predicate`: override predicate for ref structures
- `--fit-predicate`: override predicate for fit structures

### `sitehopper` subcommand

Removed. Users should use `superpose --method sitehopper`.

### `fp` and `rocs` subcommands

Add `--sim` flag to both.

## Python API

### High-level `pdist()`

```python
dm = oecluster.pdist(structures, "superpose", method="ddm")
dm = oecluster.pdist(structures, "superpose", method="sitehopper", similarity=True)

# "sitehopper" as top-level alias
dm = oecluster.pdist(structures, "sitehopper")  # routes to method="sitehopper"

# Predicates
dm = oecluster.pdist(structures, "superpose", predicate="backbone")

# Similarity mode for other metrics
dm = oecluster.pdist(mols, "fingerprint", similarity=True)
dm = oecluster.pdist(mols, "rocs", similarity=True)
```

## SWIG Bindings

- Remove all `%include`/`%ignore` for `SiteHopperMetric` and `SiteHopperOptions`.
- Add `%include "oecluster/metrics/SuperposeMetric.h"` (already present, but the
  header content changes).
- SWIG parses `enum class` via `%include` of the header. The `SuperposeMethod` and
  `SuperposeScoreType` enums will be exposed automatically.
- The `SuperposeOptions` struct with `std::string` members for predicates is
  SWIG-friendly (already have `std_string.i`).
- Ignore the private clone constructor as with other metrics.

## Python `pdist()` Signature Changes

The `pdist()` function currently accepts `metric` as a string and constructs the
metric internally. To support metric-specific kwargs:

```python
def pdist(items, metric, *, similarity=False, **kwargs):
```

The `similarity` kwarg is common to all metrics. Remaining `**kwargs` are routed
based on metric type:

- `"fingerprint"`: `fp_type`, `numbits`, `similarity_func`, etc. mapped to
  `FingerprintOptions`
- `"rocs"`: `score_type`, `color_ff_type`, etc. mapped to `ROCSOptions`
- `"superpose"`: `method`, `score_type`, `predicate`, `ref_predicate`,
  `fit_predicate` mapped to `SuperposeOptions`
- `"sitehopper"`: alias — routes to `"superpose"` with `method="sitehopper"`,
  remaining kwargs passed through

Unknown kwargs raise `TypeError`. The `similarity` flag is set on the
metric-specific options struct before construction.

## Dependencies

### New

- **oespruce**: `OESuperpose`, `OESuperposeOptions`, `OESuperposeResults`,
  `OESuperposeMethod`. Added to CMakeLists.txt via `_find_oe_dep` pattern.
- **oeselect**: Compiles predicate strings to `OEUnaryPredicate<OEAtomBase>`.
  Added via `FetchContent` from `https://github.com/scott-arne/oeselect.git`.

### Existing (unchanged)

- oechem, oebio, oeshape, oegraphsim

## Deletions

- `include/oecluster/metrics/SiteHopperMetric.h`
- `src/metrics/SiteHopperMetric.cpp`
- `tests/cpp/test_bio_metrics.cpp` (SiteHopper portions; SuperposeMetric tests
  rewritten)
- All SWIG `%include`/`%ignore` for SiteHopperMetric and SiteHopperOptions
- Python `SiteHopperMetric` wrapper class
- CLI `sitehopper` subcommand in `oepdist.cpp`

## Testing

### C++ Tests

- **Score validation**: For each `SuperposeMethod`, run against 5FQD/8G66 test
  assets, verify the expected score type is non-sentinel and unexpected ones are
  sentinel.
- **Distance/similarity modes**: Verify distance conversion matches formula per
  method; `similarity=true` returns raw scores.
- **Score type override**: Verify `Auto` selects correctly per method; explicit
  override works; incompatible combinations throw `MetricError`.
- **Predicates**: Test with oeselect predicate (e.g., `"backbone"`), verify
  different RMSD than default.
- **Clone**: Verify `Clone()` produces independent instances for parallel use.
- **Integration with pdist**: Run `pdist` with 2+ structures.

### Python Tests

- `pdist(..., "superpose", method=X)` for each method.
- `"sitehopper"` alias routing.
- `similarity=True` for superpose, fingerprint, and rocs.
- Predicate passthrough.

### CLI Tests

- `oepdist superpose` with `--method`, `--sim`, `--predicate`.
- `oepdist fp --sim` and `oepdist rocs --sim`.
