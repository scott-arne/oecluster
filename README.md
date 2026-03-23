# OECluster

High-performance pairwise distance matrix computation for molecular datasets,
built on [OpenEye Toolkits](https://www.eyesopen.com/).

OECluster provides parallel distance computation across four molecular
comparison methods -- fingerprint similarity, 3D shape overlay (ROCS),
protein superposition RMSD, and binding site comparison -- with a C++ core
library, Python bindings, and a command-line interface.

---

## Table of Contents

- [Features](#features)
- [Requirements](#requirements)
- [Installation](#installation)
  - [Python Package](#python-package)
  - [Building from Source](#building-from-source)
- [Usage](#usage)
  - [Python API](#python-api)
  - [Command-Line Tool (oepdist)](#command-line-tool-oepdist)
  - [C++ API](#c-api)
- [Distance Metrics](#distance-metrics)
- [Storage Backends](#storage-backends)
- [Project Structure](#project-structure)
- [License](#license)

---

## Features

- **Four distance metrics** for molecular comparison:
  - **Fingerprint** -- Tanimoto, Dice, Cosine, Manhattan, or Euclidean distance
    between molecular fingerprints (circular, tree, path, MACCS, lingo).
  - **ROCS** -- 3D shape and color overlay using OEShape with configurable
    scoring (ComboNorm, Combo, Shape, Color).
  - **Superpose** -- Protein superposition RMSD via sequence alignment and
    optional structural overlay.
  - **SiteHopper** -- Binding site RMSD comparison from design units.

- **Parallel computation** with automatic hardware concurrency detection,
  dynamic chunk-based scheduling, and per-thread metric cloning for
  lock-free execution.

- **Flexible storage**: dense in-memory, memory-mapped files for out-of-core
  datasets, and sparse storage with distance cutoff filtering.

- **Two computation modes**: pairwise self-distance (NxN, `pdist`) and
  cross-distance between two sets (NxM, `cdist`).

- **Python integration** with zero-copy NumPy arrays, scipy sparse matrix
  support, and serialization to `.npz` files.

- **CLI tool** (`oepdist`) for batch distance computation with NPY, CSV,
  and binary output formats.

---

## Requirements

- **C++17** compiler (GCC 9+, Clang 10+, MSVC 2019+)
- **CMake** 3.16 or later
- **OpenEye C++ Toolkits** 2025.2 or later (OEChem, OEGraphSim, OEShape, OEBio)
- **SWIG** 4.0+ (for Python bindings)
- **Python** 3.10+ with NumPy (for Python bindings)

A valid OpenEye license is required at both build time and runtime.

---

## Installation

### Python Package

Install from a pre-built wheel (if available):

```bash
pip install oecluster
```

Or build a wheel from source:

```bash
pip install scikit-build-core
python scripts/build_python.py --openeye-root /path/to/openeye/toolkits
```

### Building from Source

**1. Configure**

```bash
cmake -B build \
  -DCMAKE_BUILD_TYPE=Release \
  -DOPENEYE_ROOT=/path/to/openeye/toolkits
```

Or use the provided CMake preset (edit `CMakePresets.json` to set your
`OPENEYE_ROOT` path first):

```bash
cmake --preset default
```

**2. Build**

```bash
cmake --build build
```

This produces:

| Target             | Description                      |
|--------------------|----------------------------------|
| `liboecluster.a`   | Core C++ static library          |
| `oepdist`          | Command-line distance tool       |
| `_oecluster.so`    | Python extension module (SWIG)   |
| `oecluster_tests`  | C++ test suite                   |

**3. Run tests**

```bash
cd build && ctest --output-on-failure
```

**4. Install**

```bash
cmake --install build --prefix /usr/local
```

**Build options:**

| Option                    | Default | Description                       |
|---------------------------|---------|-----------------------------------|
| `OECLUSTER_BUILD_TESTS`  | ON      | Build C++ test suite              |
| `OECLUSTER_BUILD_PYTHON` | ON      | Build Python SWIG bindings        |
| `OECLUSTER_BUILD_TOOLS`  | ON      | Build the `oepdist` CLI tool      |
| `OECLUSTER_UNIVERSAL2`   | OFF     | macOS universal (x86_64 + arm64)  |

---

## Usage

### Python API

```python
from openeye import oechem
import oecluster

# Load molecules
mols = []
ifs = oechem.oemolistream("molecules.sdf")
mol = oechem.OEGraphMol()
while oechem.OEReadMolecule(ifs, mol):
    mols.append(oechem.OEGraphMol(mol))

# Compute pairwise fingerprint distances
dm = oecluster.pdist(mols, "fingerprint")

# Access results
print(dm.condensed)         # NumPy array (condensed form)
print(dm.squareform())      # Full square distance matrix
print(dm.labels)            # Molecule titles

# Save to disk
dm.to_file("distances.npz")
```

**With options:**

```python
dm = oecluster.pdist(
    mols,
    "fingerprint",
    fp_type="circular",
    similarity="tanimoto",
    num_threads=8,
    cutoff=0.7,
    progress=lambda done, total: print(f"{done}/{total}")
)
```

**Cross-distance (NxM) via C++ bindings:**

```python
from oecluster.oecluster import (
    FingerprintMetric, FingerprintOptions,
    CDistOptions, cdist
)
import numpy as np

opts = FingerprintOptions()
opts.fp_type = "circular"

# Concatenate both sets: [set_a... , set_b...]
all_mols = list_a + list_b
metric = FingerprintMetric([m for m in all_mols], opts)

output = np.zeros(len(list_a) * len(list_b))
cdist(metric, len(list_a), output)
result = output.reshape(len(list_a), len(list_b))
```

### Command-Line Tool (oepdist)

```
oepdist -- pairwise and cross-distance matrices

Usage: oepdist [OPTIONS] SUBCOMMAND

Subcommands:
  fp            Fingerprint distance
  rocs          ROCS shape overlay distance
  superpose     Protein superposition RMSD
  sitehopper    Binding site RMSD
```

**Fingerprint distance (NxN):**

```bash
oepdist fp molecules.sdf -o distances.npy \
  --fp-type circular --similarity tanimoto --threads 8
```

**ROCS shape overlay (NxN):**

```bash
oepdist rocs conformers.oeb -o shape_dist.npy \
  --score-type combonorm
```

**Protein superposition (NxN):**

```bash
oepdist superpose structures.oeb -o rmsd.npy \
  --alignment-method pam250 --only-calpha
```

**Cross-distance (NxM) -- supply two input files:**

```bash
oepdist fp queries.sdf targets.sdf -o cross_dist.npy
```

**Output formats** (determined by file extension):

| Extension | Format                                 |
|-----------|----------------------------------------|
| `.npy`    | NumPy binary array + JSON sidecar      |
| `.csv`    | Labeled comma-separated values         |
| `.bin`    | Raw double array + JSON sidecar        |

### C++ API

```cpp
#include <oecluster/oecluster.h>

// Load molecules into a vector<OEChem::OEMolBase*>
std::vector<OEChem::OEMolBase*> mols = /* ... */;

// Create metric
OECluster::FingerprintOptions opts;
opts.fp_type = "circular";
opts.similarity = "tanimoto";
OECluster::FingerprintMetric metric(mols, opts);

// Allocate storage and compute
OECluster::DenseStorage storage(mols.size());
OECluster::PDistOptions pdist_opts;
pdist_opts.num_threads = 8;
OECluster::pdist(metric, storage, pdist_opts);

// Access results
for (size_t i = 0; i < mols.size(); ++i)
    for (size_t j = i + 1; j < mols.size(); ++j)
        std::cout << storage.Get(i, j) << "\n";
```

---

## Distance Metrics

### Fingerprint

| Parameter      | Values                                              | Default    |
|----------------|-----------------------------------------------------|------------|
| `fp_type`      | circular, tree, path, maccs, lingo                  | circular   |
| `similarity`   | tanimoto, dice, cosine, manhattan, euclidean         | tanimoto   |
| `numbits`      | Fingerprint size                                    | 2048       |
| `min_distance`  | Minimum radius/path length                         | 0          |
| `max_distance`  | Maximum radius/path length                         | 0 (auto)   |
| `atom_type_mask` | Pipe-delimited atom features                      | (per type) |
| `bond_type_mask` | Pipe-delimited bond features                      | (per type) |

Distance: `1.0 - similarity` for Tanimoto/Dice/Cosine; raw value for Manhattan/Euclidean.

### ROCS

| Parameter      | Values                                      | Default         |
|----------------|---------------------------------------------|-----------------|
| `score_type`   | ComboNorm, Combo, Shape, Color              | ComboNorm       |
| `color_ff_type`| 1 (ImplicitMillsDean), 2 (ExplicitMillsDean)| 1               |

Distance ranges: ComboNorm [0,1], Combo [0,2], Shape [0,1], Color [0,1].

### Superpose

| Parameter          | Values                                   | Default  |
|--------------------|------------------------------------------|----------|
| `alignment_method` | PAM250 (2), BLOSUM62 (3), Gonnet, Identity| PAM250  |
| `gap_penalty`      | Sequence alignment gap penalty           | -10      |
| `extend_penalty`   | Gap extension penalty                    | -2       |
| `only_calpha`      | Use only C-alpha atoms                   | true     |
| `overlay`          | Superpose structures before RMSD         | true     |

Accepts molecules directly or design units (protein component extracted automatically).

### SiteHopper

Same alignment parameters as Superpose. Accepts design units only; always
performs overlay. Extracts and compares binding site protein components.

---

## Storage Backends

| Backend        | Use Case                      | Memory           | Thread-Safe |
|----------------|-------------------------------|------------------|-------------|
| `DenseStorage` | Default, small-medium datasets| O(N^2) in-memory | Yes         |
| `MMapStorage`  | Large datasets, persistence   | File-backed      | Yes         |
| `SparseStorage`| Cutoff-filtered results       | Entries only     | Yes         |

All backends use scipy-compatible condensed distance matrix indexing
(upper triangle, `n*i - i*(i+1)/2 + j - i - 1`).

---

## Project Structure

```
oecluster/
  include/oecluster/         C++ public headers
    metrics/                  Metric-specific headers
  src/                        C++ implementation
    metrics/                  Metric implementations
  tools/                      CLI tool (oepdist)
  swig/                       SWIG interface for Python bindings
  python/oecluster/           Python package
  tests/
    cpp/                      C++ tests (GoogleTest)
    python/                   Python tests (pytest)
  scripts/                    Build and packaging utilities
  CMakeLists.txt              Root build configuration
  pyproject.toml              Python package configuration
```

---

## License

MIT License. See [LICENSE](LICENSE) for details.
