# CLI

The `oepdist` command computes pairwise and cross-distance matrices from
molecular files without writing a Python script. It exposes one subcommand per
comparison method: `fp` (fingerprint), `rocs` (shape overlay), and `superpose`
(protein superposition, including the SiteHopper patch-score mode).

```text
oepdist — pairwise and cross-distance matrices
```

Every subcommand requires exactly one of itself, so run `oepdist fp ...`,
`oepdist rocs ...`, or `oepdist superpose ...`.

## Common Options

These options apply to all subcommands:

| Option | Meaning |
|--------|---------|
| `input` (positional) | Input structure file (required). |
| `input2` (positional) | Optional second input file; switches to cross-distance (NxM) mode. |
| `-o, --output` | Output file path (required). The extension selects the format. |
| `-t, --threads` | Worker threads (`0` = auto-detect). |
| `-c, --cutoff` | Distance cutoff for sparse output. |
| `--chunk-size` | Pairs processed per work unit. |
| `--no-progress` | Disable the progress display. |
| `-v, --verbose` | Verbose logging. |

When a single input file is given, `oepdist` computes the symmetric pairwise
matrix. When a second input file is given, it computes the rectangular
cross-distance matrix between the two sets.

## Fingerprint Distance

```bash
oepdist fp molecules.sdf -o distances.npy \
  --fp-type morgan --metric tanimoto --threads 8
```

| Option | Values | Default |
|--------|--------|---------|
| `--fp-type` | `morgan`, `atom_pair` | `morgan` |
| `--metric` | `tanimoto`, `dice`, `manhattan` | `tanimoto` |
| `--numbits` | Fingerprint size | `2048` |
| `--min-distance` | Minimum Atom Pair graph distance | `0` |
| `--max-distance` | Morgan radius or maximum Atom Pair graph distance | `2` |
| `--sim` | Return similarity instead of distance | off |

## ROCS Distance

```bash
oepdist rocs conformers.oeb -o shape_dist.npy --score combo_norm
```

| Option | Values | Default |
|--------|--------|---------|
| `--score` | `combo_norm`, `combo`, `shape`, `color` | `combo_norm` |
| `--color-ff` | Color force field name | `implicit-mills-dean` |
| `--sim` | Return similarity instead of distance | off |

ROCS reads multi-conformer molecules, so supply conformer files (for example
`.oeb`).

## Superpose Distance

```bash
oepdist superpose structures.oeb -o rmsd.npy --method global_carbon_alpha
```

| Option | Values | Default |
|--------|--------|---------|
| `--method` | `global_carbon_alpha`, `global`, `ddm`, `weighted`, `sse`, `sitehopper` | `global_carbon_alpha` |
| `--score-type` | `auto`, `rmsd`, `tanimoto`, `patch_score` | `auto` |
| `--predicate` | oeselect expression applied to both reference and fit | empty |
| `--ref-predicate` | Override predicate for reference structures | empty |
| `--fit-predicate` | Override predicate for fit structures | empty |
| `--sim` | Return similarity instead of distance | off |

Superpose accepts molecules or design units; design units are converted to
their protein components automatically. Use `--method sitehopper` for
binding-site patch-score comparison.

## Cross-Distance Mode

Pass a second input file to any subcommand to compute the NxM cross-distance
matrix between the two sets:

```bash
oepdist fp queries.sdf targets.sdf -o cross_dist.npy
```

## Output Formats

The output format is determined by the file extension:

| Extension | Format |
|-----------|--------|
| `.npy` | NumPy array plus a JSON sidecar with mode, method, parameters, and labels |
| `.csv` | Labeled comma-separated values |
| `.bin` | Raw double array plus a JSON sidecar |

Any other extension is rejected with an "Unknown output format" error.

## Running The CLI From Python

The console script is installed as `oepdist`. The Python wrapper locates the
bundled binary, configures the OpenEye runtime library path, and execs the
binary, so `oepdist ...` behaves identically whether OECluster was installed
from a wheel or built in place.
