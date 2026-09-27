# OECluster Documentation

OECluster is an OpenEye-based toolkit for molecular clustering, pairwise
distance computation, and representative selection. It is designed for
cheminformatics and medicinal-chemistry workflows where you want to go from a
set of molecules to clusters, ranked representatives, and reusable distance
matrices.

Most users should start with the quickstart and Python API guide. The
command-line tool is useful for quick file-based distance jobs, and the C++
reference is available for users embedding OECluster directly in C++
applications or extending the library.

```{toctree}
:maxdepth: 2
:caption: User Guides

quickstart
python-api
cli
benchmarks
```

```{toctree}
:maxdepth: 2
:caption: API And Developer References

cpp-core
api/python
api/cpp
developer
mcs-measurements
```
