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

### Thread Sanitizer

The `tsan` preset builds the C++ test binary with ThreadSanitizer into its own
`build-tsan/` tree, so the `build/` and `build-debug/` caches are untouched. It
turns off the Python bindings and the CLI tools, because the target is the gtest
binary and neither is needed to reach the concurrency code.

```bash
cmake --preset tsan
cmake --build build-tsan --target oecluster_tests -j8
```

TSan is slow, so run the concurrency-relevant suites rather than all 686 tests:

```bash
./build-tsan/tests/oecluster_tests \
    --gtest_filter='*ThreadPool*:*Storage*:*PDist*:*CDist*:*Parallel*:*Progress*'
./build-tsan/tests/oecluster_tests --gtest_filter='*MCS*'
```

gtest filters are case-sensitive: `*PDist*` matches `PDistTest`, `*Pdist*`
matches nothing.

No suppression file is committed, because these runs produce no warnings. If
one is ever needed, note that the toolkit is linked statically here
(`OPENEYE_USE_SHARED=OFF`), so a `called_from_lib:` entry has no shared object
to name and the suppression would have to be written against a symbol instead.

**What this covers, and what it cannot.** TSan instruments the concurrency this
project owns -- `ThreadPool`, `ParallelFor`, the storage backends, the progress
callback, and the serial clone-distribution loops in `pdist` and `cdist` -- all
of which are compiled from this tree. The OpenEye libraries are prebuilt and
uninstrumented, so TSan cannot see inside `OEMol`, `OEMCSSearch`, or any other
toolkit type: it can neither report a race there nor rule one out. A clean run
therefore says the machinery around the toolkit is race-free on the paths the
tests exercise, and says nothing about the toolkit itself.

A clean sanitizer run proves nothing until the instrumentation has been shown to
speak. Before trusting one, introduce a deliberate unsynchronised write inside a
`ParallelFor` body, confirm TSan reports it, and revert.

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
