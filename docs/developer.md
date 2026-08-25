# Developer Guide

## Documentation Build

The documentation build uses Sphinx, MyST, autodoc, Doxygen, Breathe, and
Exhale. Generated files under `docs/_build/`, `docs/_doxygen/`, and
`docs/cpp-api/` are regenerated on every build and are ignored by git.

The build needs the documentation dependency set. Install it into the active
interpreter with:

```bash
python -m invoke docs-deps
```

or directly:

```bash
uv pip install -r docs/requirements.txt
```

Build the HTML documentation:

```bash
python -m invoke docs
```

Build strictly, with warnings treated as errors (used for CI and release
checks):

```bash
python -m invoke docs-check
```

Build and serve the HTML tree on a local port:

```bash
python -m invoke serve-docs
```

Serve with live rebuilds when `sphinx-autobuild` is installed:

```bash
python -m invoke serve-docs --watch
```

The same targets are available through the Makefile in `docs/` (`make html`,
`make check`, `make clean`).

### How The C++ Reference Is Generated

Exhale runs Doxygen over `include/oecluster` during the Sphinx build and emits
the reStructuredText that becomes the C++ API pages. The Python autodoc pass
imports the `oecluster` package and introspects its public API. On a
documentation-only host (such as Read the Docs) the compiled SWIG extension and
the OpenEye toolkits are not installed, so `conf.py` mocks the SWIG layers
(`oecluster._oecluster` and `oecluster.oecluster`); NumPy must still be
importable because the package imports it at module load.

## Building The Extension

The C++/SWIG extension is built and tested with CMake. The committed presets
hard-code an interpreter path that may not exist, so override it on the command
line and use the micromamba `main` interpreter:

```bash
cmake --build build --target oecluster_python
cd build && ctest --output-on-failure
```

Run the Python suite with the same interpreter:

```bash
/Users/johnss51/Applications/micromamba/envs/main/bin/python -m pytest tests/python
```

The package version is read live from `python/oecluster/__init__.py` in an
editable install, so a documentation-only change to docstrings does not require
rebuilding the extension before the autodoc pass picks it up.

## Static Analysis

`ruff` and `mypy` cover the Python package and benchmarks. The repository's
`pyproject.toml` configures both but defines no default target, so pass the
edited files explicitly:

```bash
ruff check <edited files>
mypy <edited files>
```

The generated SWIG wrapper `python/oecluster/oecluster.py` and the build-info
module are excluded from both tools.

## Project Layout

| Path | Contents |
|------|----------|
| `include/oecluster/` | Public C++ headers (the Doxygen/Exhale input). |
| `src/` | C++ implementations and internal headers. |
| `tools/` | The `oepdist` CLI and its molecule readers/output writers. |
| `python/oecluster/` | The Python package and SWIG-generated wrapper. |
| `swig/` | The SWIG interface definitions. |
| `tests/` | C++ and Python test suites. |
| `benchmarks/` | Opt-in performance benchmark scripts. |
| `examples/` | Runnable example scripts. |
| `docs/` | This documentation. |
