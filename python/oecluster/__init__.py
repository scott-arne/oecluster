"""
oecluster: High-performance pairwise distance computation for molecular datasets.

This package provides efficient computation of pairwise distance matrices for
molecular and protein structure datasets using OpenEye toolkits. It supports
multiple comparison methods including fingerprint similarity, ROCS shape
overlay, protein superposition, and binding site comparison.
"""

import abc
import json
import hashlib
import importlib.machinery
import importlib.util
import os
import re
import shutil
import sys
import warnings
from importlib import metadata
from pathlib import Path
from typing import Any

import numpy as np

__version__ = "4.2.3"
__version_info__ = (4, 2, 3)


_OPENEYE_COMPAT_PRELOAD_PATHS: list[str] = []
_OPENEYE_COMPAT_EXTENSION_DIR: Path | None = None
__all__ = [
    "__version__",
    "__version_info__",
    "DenseStorage",
    "MMapStorage",
    "SparseStorage",
    "PDistOptions",
    "DistanceMatrix",
    "SymmetricDistanceMatrix",
    "CrossDistanceMatrix",
    "load_distance_matrix",
    "ClusteringResult",
    "ClusterReport",
    "ClusterReportComparison",
    "cluster_report",
    "compare_reports",
    "ButinaResult",
    "DBSCANResult",
    "HDBSCANResult",
    "AgglomerativeResult",
    "BitBirchResult",
    "RepresentativeMetrics",
    "ClusterRepresentative",
    "pdist",
    "cdist",
    "butina",
    "representative",
    "rank_representatives",
    "select_representatives",
    "dbscan",
    "hdbscan",
    "agglomerative",
    "bitbirch",
    "bitbirch_recluster",
    "bitbirch_refine",
    "ButinaOptions",
    "RepresentativeOptions",
    "RepresentativeWeights",
    "DBSCANOptions",
    "HDBSCANOptions",
    "AgglomerativeOptions",
    "BitBirchOptions",
    "BitBirchReclusteringOptions",
    "BitBirchRefinementOptions",
    "FingerprintComparison",
    "ROCSComparison",
    "SuperposeComparison",
    "DescriptorComparison",
    "descriptor_statistics",
]



def _user_cache_root():
    """Return the per-user cache root for OpenEye compatibility aliases."""
    cache_home = os.environ.get("XDG_CACHE_HOME")
    if cache_home:
        return Path(cache_home) / "oecluster"
    return Path.home() / ".cache" / "oecluster"


def _runtime_openeye_version():
    """Return the installed OpenEye toolkit distribution version if available."""
    try:
        return metadata.version("openeye-toolkits")
    except metadata.PackageNotFoundError:
        return "unknown"


def _cache_key(oe_lib_dir, expected_libs, build_version, runtime_version):
    """Build a stable cache key for one OpenEye runtime library set."""
    key_data = "\n".join(
        [
            os.path.realpath(oe_lib_dir),
            build_version or "unknown",
            runtime_version or "unknown",
            *sorted(expected_libs),
        ]
    )
    return hashlib.sha256(key_data.encode("utf-8")).hexdigest()[:16]


def _runtime_shared_library_names(lib_names):
    """Return filenames that can participate in runtime dynamic loading."""
    return [
        lib_name
        for lib_name in lib_names
        if ".so" in lib_name
        or lib_name.endswith(".dylib")
        or lib_name.endswith(".dll")
    ]


def _is_openeye_runtime_library_name(lib_name):
    """Return whether a dependency belongs to the OpenEye runtime set."""
    return lib_name.startswith("liboe") or lib_name.startswith("libzstd.")


def _find_openeye_runtime_lib_dir(expected_libs=()):
    """Find the OpenEye runtime library directory without importing oechem."""
    search_locations = []
    openeye_module = sys.modules.get("openeye")
    openeye_path = getattr(openeye_module, "__path__", None)
    if openeye_path is not None:
        search_locations.extend(openeye_path)

    if not search_locations:
        try:
            openeye_spec = importlib.util.find_spec("openeye")
        except (ImportError, ValueError):
            openeye_spec = None
        if (
            openeye_spec is not None
            and openeye_spec.submodule_search_locations is not None
        ):
            search_locations.extend(openeye_spec.submodule_search_locations)

    expected_libs = set(_runtime_shared_library_names(expected_libs or ()))
    fallback_dir = None
    for package_root in search_locations:
        libs_root = Path(package_root) / "libs"
        if not libs_root.is_dir():
            continue

        # Importing openeye.libs eagerly imports oechem in some environments.
        # The runtime libraries are shipped below openeye/libs, so filesystem
        # discovery preserves the fresh-import condition.
        for root, _, files in os.walk(libs_root):
            file_set = set(files)
            if expected_libs and expected_libs.intersection(file_set):
                return root
            if fallback_dir is None and any(
                ".dylib" in lib_name or ".so" in lib_name or ".dll" in lib_name
                for lib_name in files
            ):
                fallback_dir = root

    return fallback_dir


def _library_family(lib_name):
    """Return the stable library family name for a versioned shared library."""
    match = re.match(r"(lib\w+?)(-[\d.]+)?(\.[\d.]*\w+)$", lib_name)
    if match is None:
        return None
    return match.group(1)


def _candidate_runtime_libraries(oe_lib_dir, expected_name):
    """Find runtime libraries with the same family as an expected filename."""
    family = _library_family(expected_name)
    if family is None:
        return []
    candidates = []
    for file_name in os.listdir(oe_lib_dir):
        candidate_path = os.path.join(oe_lib_dir, file_name)
        if not os.path.isfile(candidate_path):
            continue
        if file_name.startswith(f"{family}-") or file_name.startswith(f"{family}."):
            candidates.append(candidate_path)
    return sorted(candidates)


def _compatible_library_path(oe_lib_dir, expected_name):
    """Return a runtime library path and whether it needs an expected-name alias."""
    exact_path = os.path.join(oe_lib_dir, expected_name)
    if os.path.isfile(exact_path):
        return exact_path, False

    candidates = _candidate_runtime_libraries(oe_lib_dir, expected_name)
    if len(candidates) != 1:
        candidate_names = ", ".join(os.path.basename(path) for path in candidates)
        raise ImportError(
            f"Could not find a compatible OpenEye runtime library for "
            f"{expected_name!r} in {oe_lib_dir!r}. "
            f"Candidates: {candidate_names or 'none'}."
        )
    return candidates[0], True


def _extension_runtime_library_names(pkg_dir):
    """Return OpenEye runtime library names recorded by the extension."""
    extension_path = _find_extension_module_path(pkg_dir)
    if extension_path is None:
        return []

    if sys.platform == "darwin":
        return _mach_o_runtime_library_names(extension_path)
    if sys.platform.startswith("linux"):
        return _elf_runtime_library_names(extension_path)
    return []


def _mach_o_runtime_library_names(extension_path):
    """Return OpenEye dylib dependencies recorded in a Mach-O extension."""
    import subprocess

    try:
        result = subprocess.run(
            ["otool", "-L", str(extension_path)],
            check=True,
            capture_output=True,
            text=True,
        )
    except (FileNotFoundError, OSError, subprocess.CalledProcessError):
        return []

    dependencies = []
    for line in result.stdout.splitlines()[1:]:
        dependency = line.strip().split(" ", 1)[0]
        lib_name = os.path.basename(dependency)
        if _is_openeye_runtime_library_name(lib_name):
            dependencies.append(lib_name)
    return dependencies


def _elf_runtime_library_names(extension_path):
    """Return OpenEye shared-library dependencies recorded in an ELF extension."""
    import subprocess

    try:
        result = subprocess.run(
            ["readelf", "-d", str(extension_path)],
            check=True,
            capture_output=True,
            text=True,
        )
    except (FileNotFoundError, OSError, subprocess.CalledProcessError):
        return []

    dependencies = []
    for match in re.finditer(r"Shared library: \[(?P<name>[^\]]+)\]", result.stdout):
        lib_name = match.group("name")
        if _is_openeye_runtime_library_name(lib_name):
            dependencies.append(lib_name)
    return dependencies


def _ensure_cache_alias(cache_dir, expected_name, target_path):
    """Create or refresh an expected-name symlink in the user cache."""
    alias_path = cache_dir / expected_name
    if alias_path.is_symlink():
        if alias_path.resolve() == Path(target_path).resolve():
            return alias_path
        alias_path.unlink()
    elif alias_path.exists():
        raise ImportError(
            f"Cannot create OpenEye compatibility alias {alias_path}: "
            "a non-symlink file already exists at that path."
        )

    try:
        alias_path.symlink_to(target_path)
    except OSError as exc:
        raise ImportError(
            f"Could not create OpenEye compatibility alias "
            f"{alias_path} -> {target_path}: {exc}"
        ) from exc
    return alias_path


def _ensure_library_compat():
    """Prepare compatibility aliases when OpenEye library filenames drift.

    When oecluster is built with shared OpenEye libraries, the compiled extension
    records the exact versioned library filenames (e.g., liboechem-4.3.0.1.dylib).
    If the user upgrades openeye-toolkits, these filenames change and the dynamic
    linker fails to load the extension.

    This function creates expected-name aliases in a user-writable cache instead
    of mutating the installed package directory. When aliases are needed, the
    extension is later loaded from the same cache directory so its $ORIGIN lookup
    can find those aliases.
    """
    global _OPENEYE_COMPAT_EXTENSION_DIR, _OPENEYE_COMPAT_PRELOAD_PATHS

    _OPENEYE_COMPAT_PRELOAD_PATHS = []
    _OPENEYE_COMPAT_EXTENSION_DIR = None

    try:
        from . import _build_info
    except ImportError:
        return False

    if getattr(_build_info, 'OPENEYE_LIBRARY_TYPE', 'STATIC') != 'SHARED':
        return False

    expected_libs = set(_runtime_shared_library_names(
        getattr(_build_info, 'OPENEYE_EXPECTED_LIBS', [])
    ))
    expected_libs.update(_extension_runtime_library_names(os.path.dirname(__file__)))
    expected_libs = sorted(expected_libs)
    if not expected_libs:
        return False

    oe_lib_dir = _find_openeye_runtime_lib_dir(expected_libs)
    if oe_lib_dir is None:
        return False

    if not os.path.isdir(oe_lib_dir):
        return False

    build_version = getattr(_build_info, 'OPENEYE_BUILD_VERSION', None)
    runtime_version = _runtime_openeye_version()
    cache_dir = (
        _user_cache_root()
        / "openeye-libs"
        / _cache_key(oe_lib_dir, expected_libs, build_version, runtime_version)
    )

    preload_paths = []
    needs_cached_origin = False
    for expected_name in expected_libs:
        actual_path, needs_alias = _compatible_library_path(oe_lib_dir, expected_name)
        if needs_alias:
            try:
                cache_dir.mkdir(parents=True, exist_ok=True)
            except OSError as exc:
                raise ImportError(
                    f"Could not create OpenEye compatibility cache directory "
                    f"{cache_dir}: {exc}"
                ) from exc
            alias_path = _ensure_cache_alias(cache_dir, expected_name, actual_path)
            preload_paths.append(str(alias_path))
            needs_cached_origin = True
        else:
            preload_paths.append(actual_path)

    _OPENEYE_COMPAT_PRELOAD_PATHS = preload_paths
    if needs_cached_origin:
        _OPENEYE_COMPAT_EXTENSION_DIR = cache_dir

    return needs_cached_origin


def _extension_suffixes():
    """Return extension-module suffixes for the active Python interpreter."""
    return tuple(importlib.machinery.EXTENSION_SUFFIXES)


def _find_extension_module_path(pkg_dir):
    """Find the installed _oecluster extension file."""
    for suffix in _extension_suffixes():
        candidate = Path(pkg_dir) / f"_oecluster{suffix}"
        if candidate.is_file():
            return candidate
    for candidate in Path(pkg_dir).glob("_oecluster*"):
        if candidate.is_file() and str(candidate).endswith(_extension_suffixes()):
            return candidate
    return None


def _copy_if_stale(source_path, target_path):
    """Copy a file into the cache when size or mtime changed."""
    if (
        target_path.exists()
        and target_path.stat().st_size == source_path.stat().st_size
        and target_path.stat().st_mtime_ns == source_path.stat().st_mtime_ns
    ):
        return
    shutil.copy2(source_path, target_path)


def _copy_package_shared_sidecars(pkg_dir, cache_dir, extension_path):
    """Copy package-local shared library sidecars needed by cached extension."""
    for candidate in Path(pkg_dir).iterdir():
        name = candidate.name
        if not candidate.is_file() or candidate == extension_path:
            continue
        if (
            ".so" not in name
            and not name.endswith(".dylib")
            and not name.endswith(".dll")
            and not name.endswith(".pyd")
        ):
            continue
        _copy_if_stale(candidate, cache_dir / name)


def _load_cached_extension_if_needed():
    """Load _oecluster from the cache when OpenEye aliases live there."""
    cache_dir = _OPENEYE_COMPAT_EXTENSION_DIR
    if cache_dir is None:
        return

    module_name = f"{__name__}._oecluster"
    if module_name in sys.modules:
        return

    pkg_dir = os.path.dirname(__file__)
    extension_path = _find_extension_module_path(pkg_dir)
    if extension_path is None:
        return

    cached_extension_path = cache_dir / extension_path.name
    try:
        cache_dir.mkdir(parents=True, exist_ok=True)
        _copy_if_stale(extension_path, cached_extension_path)
        _copy_package_shared_sidecars(pkg_dir, cache_dir, extension_path)
    except OSError as exc:
        raise ImportError(
            f"Could not prepare cached oecluster extension in {cache_dir}: {exc}"
        ) from exc

    spec = importlib.util.spec_from_file_location(module_name, cached_extension_path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Could not create import spec for {cached_extension_path}")

    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    try:
        spec.loader.exec_module(module)
    except Exception:
        sys.modules.pop(module_name, None)
        raise


def _preload_shared_libs():
    """Preload OpenEye shared libraries so the C extension can find them.

    On Linux, the extension's RUNPATH (set at build time) normally handles
    dependency resolution, but preloading ensures libraries are available
    even if RUNPATH is stripped (e.g. by certain packaging tools).
    On macOS, @rpath references may not resolve without preloading.

    Only the libraries recorded in ``OPENEYE_EXPECTED_LIBS`` are loaded,
    and they are loaded with ``RTLD_GLOBAL`` so that cross-module C++
    symbol references resolve correctly. Loading the entire OpenEye
    library directory (which can contain 70+ unrelated shared objects)
    would pollute the global symbol namespace and cause segfaults in
    unrelated C extensions such as ``_sqlite3``.
    """
    import ctypes
    import sys
    if sys.platform not in ('linux', 'darwin'):
        return

    try:
        from . import _build_info
    except ImportError:
        return

    if getattr(_build_info, 'OPENEYE_LIBRARY_TYPE', 'STATIC') != 'SHARED':
        return

    expected_libs = _runtime_shared_library_names(
        getattr(_build_info, 'OPENEYE_EXPECTED_LIBS', [])
    )
    if not expected_libs:
        return

    oe_lib_dir = _find_openeye_runtime_lib_dir(expected_libs)
    if oe_lib_dir is None:
        return

    if not os.path.isdir(oe_lib_dir):
        return

    paths = _OPENEYE_COMPAT_PRELOAD_PATHS
    if not paths:
        paths = [
            os.path.join(oe_lib_dir, lib_name)
            for lib_name in expected_libs
            if os.path.exists(os.path.join(oe_lib_dir, lib_name))
        ]

    for path in paths:
        if os.path.exists(path) or os.path.islink(path):
            try:
                ctypes.CDLL(path, mode=ctypes.RTLD_GLOBAL)
            except OSError:
                pass

def _preload_bundled_libs():
    """Preload libraries bundled by auditwheel from the .libs directory.

    auditwheel repair bundles non-manylinux dependencies (e.g. libraries
    from FetchContent or system packages) into a ``<package>.libs/``
    directory next to the package. The bundled copies have hashed filenames
    and must be loaded before the C extension to satisfy its DT_NEEDED
    entries.

    Libraries may have inter-dependencies, so we do multiple passes
    until no new libraries can be loaded. Libraries are loaded without
    ``RTLD_GLOBAL`` to avoid polluting the global symbol namespace.
    """
    import sys
    if sys.platform != 'linux':
        return

    import ctypes
    pkg_name = __name__
    pkg_dir = os.path.dirname(os.path.abspath(__file__))
    site_dir = os.path.dirname(pkg_dir)
    for libs_name in (f'{pkg_name}.libs', f'.{pkg_name}.libs'):
        libs_dir = os.path.join(site_dir, libs_name)
        if not os.path.isdir(libs_dir):
            continue
        remaining = [
            os.path.join(libs_dir, f)
            for f in sorted(os.listdir(libs_dir))
            if '.so' in f
        ]
        while remaining:
            failed = []
            for lib_path in remaining:
                try:
                    ctypes.CDLL(lib_path)
                except OSError:
                    failed.append(lib_path)
            if len(failed) == len(remaining):
                break
            remaining = failed


def _check_openeye_version():
    """Check that the OpenEye version matches what was used at build time."""
    try:
        from . import _build_info
    except ImportError:
        return

    if getattr(_build_info, 'OPENEYE_LIBRARY_TYPE', 'STATIC') != 'SHARED':
        return

    build_version = getattr(_build_info, 'OPENEYE_BUILD_VERSION', None)
    if not build_version:
        return

    try:
        runtime_version = metadata.version("openeye-toolkits")
    except metadata.PackageNotFoundError:
        warnings.warn(
            "openeye-toolkits package not found. "
            "This package requires openeye-toolkits to be installed.",
            ImportWarning
        )
        return

    build_parts = build_version.split('.')[:2]
    runtime_parts = runtime_version.split('.')[:2]
    if build_parts != runtime_parts:
        warnings.warn(
            f"OpenEye version mismatch: oecluster was built with OpenEye Toolkits {build_version} "
            f"but runtime has OpenEye Toolkits {runtime_version}. "
            f"This may cause compatibility issues.",
            RuntimeWarning
        )


# Initialize compatibility checks
_ensure_library_compat()
_preload_shared_libs()
_preload_bundled_libs()
_load_cached_extension_if_needed()
_check_openeye_version()

# Import C++ bindings from SWIG module
try:
    from .oecluster import (
        DenseStorage,
        MMapStorage,
        SparseStorage,
        PDistOptions,
        CDistOptions,
        ButinaOptions,
        RepresentativeOptions,
        RepresentativeWeights,
        DBSCANOptions,
        HDBSCANOptions,
        AgglomerativeOptions,
        BitBirchOptions,
        BitBirchReclusteringOptions,
        BitBirchRefinementOptions,
        pdist as _cpp_pdist,
        cdist_into_address as _cpp_cdist_into_address,
        cluster_report as _cluster_report,
        butina_cluster as _butina_cluster,
        cluster_representative as _cluster_representative,
        rank_representatives as _rank_representatives,
        select_representatives as _select_representatives,
        dbscan_cluster as _dbscan_cluster,
        hdbscan_cluster as _hdbscan_cluster,
        agglomerative_cluster as _agglomerative_cluster,
        bitbirch_cluster as _bitbirch_cluster,
        bitbirch_recluster as _bitbirch_recluster,
        bitbirch_refine as _bitbirch_refine,
    )
except ImportError as e:
    raise ImportError(
        f"Failed to import _oecluster C++ extension: {e}. "
        "The package may not be built correctly."
    ) from e

from . import oecluster as _oecluster

# Docstrings for SWIG-imported storage backend classes
DenseStorage.__doc__ = """
Default in-memory storage backend for dense distance matrices.

Stores the full condensed distance matrix in contiguous memory. Best for
datasets where most distances will be accessed and the matrix fits in RAM.
"""

MMapStorage.__doc__ = """
Memory-mapped file storage backend for large distance matrices.

Backs the condensed distance matrix with a memory-mapped file, allowing
out-of-core computation for datasets larger than available RAM. The file
persists after computation and can be reloaded.
"""

SparseStorage.__doc__ = """
Sparse storage backend that stores only distances below a cutoff.

Stores only non-zero entries (distances below a threshold) in a sparse
representation. Efficient for large, sparse distance graphs where most
pairwise distances exceed the cutoff.
"""

PDistOptions.__doc__ = """
Options for parallel pairwise-distance computation.

:ivar num_threads: Number of threads (0 = auto-detect).
:ivar chunk_size: Pairs processed per work unit.
:ivar cutoff: Distance cutoff for sparse storage (0 = store all).
:ivar progress: Optional callback(completed, total) for progress reporting.
"""

from .oecluster import DescriptorComparison as _DescriptorComparison
from .oecluster import FingerprintComparison as _FingerprintComparison
from .oecluster import FingerprintOptions
from .oecluster import ROCSComparison as _ROCSComparison
from .oecluster import ROCSOptions
from .oecluster import SuperposeComparison as _SuperposeComparison
from .oecluster import SuperposeOptions

from . import _comparisons
from . import _gate


class _StorageView:
    """Expose a storage backend's buffer to numpy while owning a reference.

    A zero-copy view is only safe if the buffer outlives every array over it.
    Building the array from a raw pointer does not establish that: the pointer
    carries no ownership, so the array stays valid only for as long as
    something else happens to hold the storage. When the array outlives its
    matrix -- ``arr = pdist(...).condensed`` discards the matrix immediately --
    the buffer is freed and the array reads whatever the allocator puts there
    next. That fails silently, returning plausible distances rather than
    crashing.

    Handing numpy an ``__array_interface__`` on an object that holds the
    storage makes the array's ``base`` this view, so the storage is reachable
    for exactly as long as the array is.

    :ivar _storage: The storage backend whose buffer is exposed. Held solely
        to keep it alive; never read.
    """

    def __init__(self, storage, ptr, num_pairs):
        """
        Wrap a storage buffer for numpy without copying it.

        :param storage: Storage backend owning the buffer.
        :param ptr: Address of the first ``double`` in the buffer.
        :param num_pairs: Number of ``double`` values in the buffer.
        """
        self._storage = storage
        self.__array_interface__ = {
            'data': (int(ptr), False),
            'shape': (num_pairs,),
            'typestr': np.dtype(np.float64).str,
            'version': 3,
        }


class DistanceMatrix(abc.ABC):
    """
    Abstract base for distance matrices computed from a set of items.

    Subclasses model specific shapes: :class:`SymmetricDistanceMatrix` for the
    symmetric within-set result of :func:`pdist`, and
    :class:`CrossDistanceMatrix` for the rectangular cross-set result of
    :func:`cdist`. The base holds only the shared comparison metadata.
    """

    def __init__(self, comparison_name, params=None, facts=None):
        """
        Initialize shared distance-matrix metadata.

        :param comparison_name: Name of the comparison method used.
        :param params: Optional dictionary of comparison parameters.
        :param facts: Optional capability facts recorded by the comparison.
            Anything omitted stays at its "unknown" default, which the
            clustering gate never refuses.
        """
        self._comparison_name = comparison_name
        self._params = params if params is not None else {}
        self._facts = _gate.default_facts()
        if facts:
            self._facts.update(facts)

    @property
    def comparison_name(self):
        """Get the name of the comparison method."""
        return self._comparison_name

    @property
    def params(self):
        """Get the comparison parameters dictionary."""
        return self._params

    @property
    def facts(self):
        """
        Get a copy of the capability facts recorded for this matrix.

        :returns: Dict with ``is_distance``, ``zero_self``, ``triangle``,
                  ``data_integrity``, ``metric_probe``, ``probe_violations``,
                  ``probe_sampled``.
        """
        return dict(self._facts)

    @property
    def is_distance(self):
        """
        Get whether a larger entry means "further apart".

        This is the matrix's orientation, and it is independent of the two
        metric axioms: a similarity is not a distance however its diagonal
        behaves.

        :returns: True, False, or the string ``"unknown"``.
        """
        return self._facts['is_distance']

    @property
    def metric_capabilities(self):
        """
        Get the two metric axioms the clustering gate checks.

        Orientation is not one of them -- read :attr:`is_distance` for that.

        :returns: Dict with ``zero_self`` and ``triangle``, each True, False,
                  or the string ``"unknown"``.
        """
        return {'zero_self': self._facts['zero_self'],
                'triangle': self._facts['triangle']}

    @property
    def data_integrity(self):
        """Get "complete", "nan_present", "subset_scored", or "unknown"."""
        return self._facts['data_integrity']

    @property
    def metric_probe(self):
        """Get "not_run", "no_violations_found", or "violations_found"."""
        return self._facts['metric_probe']

    @property
    def probe_violations(self):
        """Get the number of triple violations the probe found."""
        return self._facts['probe_violations']

    @property
    def probe_sampled(self):
        """Get the number of triples the probe sampled."""
        return self._facts['probe_sampled']

    @abc.abstractmethod
    def to_file(self, path):
        """Save the distance matrix to a compressed ``.npz`` file."""
        raise NotImplementedError


class SymmetricDistanceMatrix(DistanceMatrix):
    """
    A symmetric (within-set) distance matrix produced by :func:`pdist`.

    Provides zero-copy access to condensed distance matrix data via numpy arrays,
    conversion to scipy sparse matrices, and serialization support.
    """

    def __init__(self, storage, comparison_name, labels=None, params=None,
                 facts=None):
        """
        Construct a SymmetricDistanceMatrix wrapper.

        :param storage: C++ StorageBackend instance.
        :param comparison_name: Name of the comparison method used.
        :param labels: Optional list of labels for items.
        :param params: Optional dictionary of comparison parameters.
        :param facts: Optional capability facts recorded by the comparison.
        """
        super().__init__(comparison_name, params, facts)
        self._storage = storage
        self._labels = labels if labels is not None else []
        self._condensed_cache = None

    @property
    def storage(self):
        """Get the underlying storage backend."""
        return self._storage

    @property
    def labels(self):
        """Get the item labels."""
        return self._labels

    @property
    def num_samples(self):
        """Get the number of items in the dataset."""
        return self._storage.NumSamples()

    @property
    def num_pairs(self):
        """Get the number of pairwise distances."""
        return self._storage.NumPairs()

    @property
    def shape(self):
        """Get the shape of the full distance matrix (N, N)."""
        n = self.num_samples
        return (n, n)

    @property
    def condensed(self):
        """
        Get the condensed distance array as a numpy array (zero-copy when possible).

        :returns: 1D numpy array of pairwise distances.
        """
        if self._condensed_cache is not None:
            return self._condensed_cache

        # Check if this is sparse storage
        if isinstance(self._storage, SparseStorage):
            # Sparse storage requires dense conversion
            n = self.num_samples
            num_pairs = self.num_pairs
            condensed = np.zeros(num_pairs, dtype=np.float64)

            # Get sparse entries
            entries = self._storage._entries()
            for i, j, val in entries:
                # Convert (i,j) to condensed index
                if i < j:
                    idx = n * i + j - ((i + 2) * (i + 1)) // 2
                else:
                    idx = n * j + i - ((j + 2) * (j + 1)) // 2
                condensed[idx] = val

            self._condensed_cache = condensed
            return condensed

        # Dense or MMap storage: zero-copy access. The view owns a reference to
        # the storage, so the buffer cannot be freed while this array is alive
        # -- caching on the matrix is not enough, because the array is routinely
        # kept after the matrix is dropped.
        arr = np.asarray(
            _StorageView(self._storage, self._storage._data_ptr(),
                         self.num_pairs))

        self._condensed_cache = arr
        return arr

    @property
    def sparse(self):
        """
        Get the distance matrix as a scipy sparse COO matrix.

        Only available if scipy is installed. For sparse storage, only non-zero
        entries are included. For dense storage, all entries are converted.

        :returns: scipy.sparse.coo_matrix in square form.
        :raises ImportError: If scipy is not installed.
        """
        try:
            from scipy.sparse import coo_matrix
        except ImportError as e:
            raise ImportError(
                "scipy is required for sparse matrix support. "
                "Install it with: pip install scipy"
            ) from e

        n = self.num_samples

        if isinstance(self._storage, SparseStorage):
            # Get sparse entries directly
            entries = self._storage._entries()
            if not entries:
                return coo_matrix((n, n), dtype=np.float64)

            rows, cols, data = zip(*entries)
            return coo_matrix(
                (
                    np.asarray(data, dtype=np.float64),
                    (
                        np.asarray(rows, dtype=np.intp),
                        np.asarray(cols, dtype=np.intp),
                    ),
                ),
                shape=(n, n),
                dtype=np.float64,
            )

        # Dense/MMap: convert condensed to COO
        condensed = self.condensed
        rows = []
        cols = []
        data = []

        idx = 0
        for i in range(n):
            for j in range(i + 1, n):
                val = condensed[idx]
                if val != 0.0:  # Only store non-zero
                    rows.append(i)
                    cols.append(j)
                    data.append(val)
                idx += 1

        return coo_matrix(
            (
                np.asarray(data, dtype=np.float64),
                (
                    np.asarray(rows, dtype=np.intp),
                    np.asarray(cols, dtype=np.intp),
                ),
            ),
            shape=(n, n),
            dtype=np.float64,
        )

    def squareform(self):
        """
        Convert the condensed distance matrix to square form.

        :returns: 2D numpy array of shape (n, n) with zeros on diagonal.
        :raises ImportError: If scipy is not installed.
        """
        try:
            from scipy.spatial.distance import squareform
        except ImportError as e:
            raise ImportError(
                "scipy is required for squareform conversion. "
                "Install it with: pip install scipy"
            ) from e

        return squareform(self.condensed)

    def to_file(self, path):
        """
        Save the distance matrix to a compressed .npz file.

        Sparse storage is written as its entry list rather than as a condensed
        array. ``.condensed`` densifies a sparse matrix by filling every pair
        the cutoff omitted with ``0.0``, and those are the *farthest* pairs, so
        saving it would hand every algorithm "identical" for the pairs the
        cutoff called most distant.

        :param path: Output file path.
        :raises ValueError: If sparse storage holds an entry that could not be
            loaded back from the file.
        """
        arrays = {
            'comparison_name': np.array(self._comparison_name),
            'params_json': np.array(json.dumps(self._params)),
            'facts_json': np.array(json.dumps(self._facts)),
            'labels': np.array(self._labels),
            'num_samples': np.array(self.num_samples),
        }

        if isinstance(self._storage, SparseStorage):
            # Written verbatim, duplicates included: ThresholdGraph tests every
            # tuple Entries() returns, so a "tidied" list would cluster
            # differently from the matrix that was saved.
            entries = self._storage._entries()
            cutoff = self._storage.Cutoff()
            rows = np.array([e[0] for e in entries], dtype=np.int64)
            cols = np.array([e[1] for e in entries], dtype=np.int64)
            values = np.array([e[2] for e in entries], dtype=np.float64)
            # Checked before anything is written, not only on load.
            # ``SparseStorage.Set`` enforces ``i != j`` with a bare ``assert``,
            # compiled out of release builds, so a poked-at matrix reached here
            # and wrote a file ``from_file`` then refused -- the caller lost the
            # data and only found out on the next load.
            self._validate_sparse_entries(
                rows, cols, values, self.num_samples, cutoff,
                "Refusing to write a malformed sparse matrix")
            arrays['storage_kind'] = np.array("sparse")
            arrays['sparse_cutoff'] = np.array(cutoff, dtype=np.float64)
            arrays['sparse_i'] = rows
            arrays['sparse_j'] = cols
            arrays['sparse_v'] = values
        else:
            arrays['storage_kind'] = np.array("dense")
            arrays['condensed'] = self.condensed

        np.savez_compressed(path, **arrays)

    @classmethod
    def from_file(cls, path):
        """
        Load a symmetric distance matrix from a .npz file.

        :param path: Input file path.
        :returns: SymmetricDistanceMatrix instance.
        :raises ValueError: If the file is a cross-distance matrix or malformed.
        """
        data = np.load(path, allow_pickle=False)
        if 'matrix_kind' in data:
            kind = str(data['matrix_kind'])
            if kind == "cross":
                raise ValueError(
                    "File is a cross-distance matrix; use "
                    "CrossDistanceMatrix.from_file or load_distance_matrix")
            raise ValueError(
                f"Unknown matrix_kind {kind!r}; a symmetric matrix file must not "
                f"carry a matrix_kind key")
        # Absence of ``storage_kind`` means dense, so every file written before
        # sparse serialization existed still loads exactly as it did.
        storage_kind = str(data['storage_kind']) if 'storage_kind' in data \
            else "dense"
        if storage_kind not in ("dense", "sparse"):
            raise ValueError(
                f"Unknown storage_kind {storage_kind!r}; expected 'dense' or "
                f"'sparse'")

        if storage_kind == "sparse":
            required_keys = ('comparison_name', 'sparse_cutoff', 'sparse_i',
                             'sparse_j', 'sparse_v')
        else:
            required_keys = ('condensed', 'comparison_name')
        for required in required_keys:
            if required not in data:
                raise ValueError(
                    f"Malformed symmetric matrix: missing required key {required!r}")

        comparison_name = str(data['comparison_name'])
        labels = list(data.get('labels', np.array([])))

        try:
            params = json.loads(str(data['params_json']))
        except (KeyError, json.JSONDecodeError):
            params = {}

        try:
            facts = json.loads(str(data['facts_json']))
        except (KeyError, json.JSONDecodeError):
            facts = None

        # Accept the legacy ``num_items`` key so distance matrices saved before
        # the rename to ``num_samples`` still load.
        if 'num_samples' in data:
            num_samples = int(data['num_samples'])
        else:
            num_samples = int(data['num_items'])

        if storage_kind == "sparse":
            storage = cls._sparse_storage_from_file(data, num_samples)
        else:
            condensed = data['condensed']
            expected = num_samples * (num_samples - 1) // 2
            if condensed.shape[0] != expected:
                raise ValueError(
                    f"Malformed symmetric matrix: condensed length "
                    f"{condensed.shape[0]} != expected {expected} for "
                    f"{num_samples} samples")
            storage = DenseStorage(num_samples)

            idx = 0
            for i in range(num_samples):
                for j in range(i + 1, num_samples):
                    storage.Set(i, j, float(condensed[idx]))
                    idx += 1

        return cls(storage, comparison_name, labels, params, facts)

    @staticmethod
    def _validate_sparse_entries(rows, cols, values, num_samples, cutoff,
                                 subject):
        """
        Refuse a sparse entry list that would not replay into the same matrix.

        Shared by :meth:`to_file` and :meth:`from_file` so the write and read
        sides cannot drift: anything ``to_file`` accepts, ``from_file``
        accepts. The rules are what a replay needs. ``SparseStorage.Set``
        silently drops a value above the cutoff, its ``i != j`` guard is a bare
        ``assert`` compiled out of release builds, and ``ThresholdGraph``
        indexes an ``n``-element neighbour vector with the entry indices.

        Non-finite values are deliberately allowed: a sparse matrix may
        legitimately carry them, and the metric gate is what refuses them at
        clustering time.

        :param rows: First index of each entry.
        :param cols: Second index of each entry.
        :param values: Distance of each entry.
        :param num_samples: Number of samples the indices must fall below.
        :param cutoff: Sparse storage cutoff the values must not exceed.
        :param subject: Prefix naming what is being refused.
        :raises ValueError: If any entry would not survive a replay.
        """
        for name, array in (('sparse_i', rows), ('sparse_j', cols),
                            ('sparse_v', values)):
            if array.ndim != 1:
                raise ValueError(
                    f"{subject}: {name} must be a 1-D array, not "
                    f"{array.ndim}-D")
        # Ahead of the comparisons below, which only mean what they say for an
        # integer index. A float64 index array satisfies every one of them and
        # is then truncated by the ``int(i)`` in the replay loop, so a file
        # holding (0.7, 1.9) loads as the pair (0, 1) with nothing raised.
        for name, array in (('sparse_i', rows), ('sparse_j', cols)):
            if not np.issubdtype(array.dtype, np.integer):
                raise ValueError(
                    f"{subject}: {name} must hold integer indices, not "
                    f"dtype {array.dtype}")
        if not (rows.shape[0] == cols.shape[0] == values.shape[0]):
            raise ValueError(
                f"{subject}: sparse entry arrays have lengths "
                f"{rows.shape[0]}, {cols.shape[0]}, {values.shape[0]}")

        # Elementwise, not ``values.max() > cutoff``: max() propagates NaN and
        # ``nan > cutoff`` is False, so one NaN entry switched the whole check
        # off and every genuinely over-cutoff value was then silently dropped
        # by Set -- the exact reload-a-different-matrix failure this prevents.
        above = np.flatnonzero(values > cutoff)
        if above.size:
            first = int(above[0])
            raise ValueError(
                f"{subject}: sparse entry {first} has value "
                f"{float(values[first])}, above the cutoff {cutoff}")

        bad = np.flatnonzero(
            (rows < 0) | (cols < 0) | (rows >= num_samples)
            | (cols >= num_samples) | (rows == cols))
        if bad.size:
            first = int(bad[0])
            raise ValueError(
                f"{subject}: sparse entry {first} is "
                f"({int(rows[first])}, {int(cols[first])}), not a pair of "
                f"distinct indices below {num_samples}")

    @staticmethod
    def _sparse_storage_from_file(data, num_samples):
        """
        Rebuild sparse storage from a saved entry list.

        :param data: Open ``.npz`` archive.
        :param num_samples: Sample count recorded in the file.
        :returns: A finalized :class:`SparseStorage`.
        :raises ValueError: If the saved entries cannot be replayed faithfully.
        """
        # Checked before the conversion, not after: on the numpy in use,
        # ``float()`` raises TypeError on a ``(1,)``-shaped array, and TypeError
        # escapes the ValueError ``from_file`` documents for a malformed file.
        # Hence ``ndim``, not ``size`` -- a one-element array is a single value
        # but still not the 0-d scalar ``to_file`` writes.
        stored_cutoff = data['sparse_cutoff']
        if stored_cutoff.ndim != 0:
            raise ValueError(
                f"Malformed symmetric matrix: sparse_cutoff must be a scalar, "
                f"not a {stored_cutoff.ndim}-D array")
        cutoff = float(stored_cutoff)
        rows = data['sparse_i']
        cols = data['sparse_j']
        values = data['sparse_v']
        # Kept on the load side as well as in to_file: from_file reads
        # untrusted input. It is the same rule, so the two cannot disagree
        # about which files are legal.
        SymmetricDistanceMatrix._validate_sparse_entries(
            rows, cols, values, num_samples, cutoff,
            "Malformed symmetric matrix")

        storage = SparseStorage(num_samples, cutoff)
        for i, j, value in zip(rows, cols, values):
            storage.Set(int(i), int(j), float(value))
        storage.Finalize()
        return storage

    def __array__(self):
        """Support numpy array interface."""
        return self.condensed

    def __len__(self):
        """Return the number of pairwise distances."""
        return self.num_pairs

    def __repr__(self):
        return (f"SymmetricDistanceMatrix(comparison={self._comparison_name!r}, "
                f"num_samples={self.num_samples}, num_pairs={self.num_pairs})")


class CrossDistanceMatrix(DistanceMatrix):
    """
    A rectangular cross-distance matrix produced by :func:`cdist`.

    Holds an ``(n_a, n_b)`` numpy array of distances between two item sets,
    where entry ``[i, j]`` compares ``items_a[i]`` (reference) against
    ``items_b[j]`` (fit). Unlike :class:`SymmetricDistanceMatrix`, it has no
    condensed or squareform representation — cross-distances are not symmetric.
    """

    def __init__(self, matrix, comparison_name, labels_a=None, labels_b=None,
                 params=None, facts=None):
        """
        Construct a CrossDistanceMatrix wrapper.

        :param matrix: 2D numpy array of shape (n_a, n_b), dtype float64.
        :param comparison_name: Name of the comparison method used.
        :param labels_a: Optional labels for set A (the rows / reference side).
        :param labels_b: Optional labels for set B (the columns / fit side).
        :param params: Optional dictionary of comparison parameters.
        :param facts: Optional capability facts recorded by the comparison.
        """
        super().__init__(comparison_name, params, facts)
        self._matrix = matrix
        self._labels_a = labels_a if labels_a is not None else []
        self._labels_b = labels_b if labels_b is not None else []

    @property
    def matrix(self):
        """Get the rectangular cross-distance array of shape (n_a, n_b)."""
        return self._matrix

    @property
    def shape(self):
        """Get the shape of the cross-distance matrix (n_a, n_b)."""
        return tuple(self._matrix.shape)

    @property
    def labels_a(self):
        """Get the labels for set A (rows / reference side)."""
        return self._labels_a

    @property
    def labels_b(self):
        """Get the labels for set B (columns / fit side)."""
        return self._labels_b

    def to_file(self, path):
        """
        Save the cross-distance matrix to a compressed .npz file.

        :param path: Output file path.
        """
        n_a, n_b = self._matrix.shape
        np.savez_compressed(
            path,
            matrix_kind=np.array("cross"),
            matrix=self._matrix,
            n_a=np.array(n_a),
            n_b=np.array(n_b),
            labels_a=np.array(self._labels_a),
            labels_b=np.array(self._labels_b),
            comparison_name=np.array(self._comparison_name),
            params_json=np.array(json.dumps(self._params)),
            facts_json=np.array(json.dumps(self._facts)),
        )

    @classmethod
    def from_file(cls, path):
        """
        Load a cross-distance matrix from a .npz file.

        :param path: Input file path.
        :returns: CrossDistanceMatrix instance.
        :raises ValueError: If the file is not a cross-distance matrix or is malformed.
        """
        data = np.load(path, allow_pickle=False)
        if 'matrix_kind' not in data or str(data['matrix_kind']) != "cross":
            raise ValueError(
                "File is a symmetric distance matrix; use "
                "SymmetricDistanceMatrix.from_file or load_distance_matrix")
        for required in ('matrix', 'n_a', 'n_b', 'comparison_name'):
            if required not in data:
                raise ValueError(
                    f"Malformed cross matrix: missing required key {required!r}")
        matrix = data['matrix']
        n_a = int(data['n_a'])
        n_b = int(data['n_b'])
        if matrix.shape != (n_a, n_b):
            raise ValueError(
                f"Malformed cross matrix: matrix shape {tuple(matrix.shape)} "
                f"!= ({n_a}, {n_b})")
        labels_a = list(data.get('labels_a', np.array([])))
        labels_b = list(data.get('labels_b', np.array([])))
        if labels_a and len(labels_a) != n_a:
            raise ValueError(
                f"Malformed cross matrix: labels_a length {len(labels_a)} != {n_a}")
        if labels_b and len(labels_b) != n_b:
            raise ValueError(
                f"Malformed cross matrix: labels_b length {len(labels_b)} != {n_b}")
        try:
            params = json.loads(str(data['params_json']))
        except (KeyError, json.JSONDecodeError):
            params = {}
        try:
            facts = json.loads(str(data['facts_json']))
        except (KeyError, json.JSONDecodeError):
            facts = None
        return cls(matrix, str(data['comparison_name']), labels_a, labels_b,
                   params, facts)

    def __array__(self):
        """Support numpy array interface."""
        return self._matrix

    def __len__(self):
        """Return the number of rows (size of set A)."""
        return self._matrix.shape[0]

    def __repr__(self):
        n_a, n_b = self._matrix.shape
        return (f"CrossDistanceMatrix(comparison={self._comparison_name!r}, "
                f"shape=({n_a}, {n_b}))")


def load_distance_matrix(path):
    """
    Load a distance matrix from a .npz file, dispatching on its stored kind.

    Symmetric files (written by :class:`SymmetricDistanceMatrix`) have no
    ``matrix_kind`` key; cross files carry ``matrix_kind="cross"``.

    :param path: Input file path.
    :returns: A :class:`SymmetricDistanceMatrix` or :class:`CrossDistanceMatrix`.
    """
    data = np.load(path, allow_pickle=False)
    if 'matrix_kind' in data:
        kind = str(data['matrix_kind'])
        if kind == "cross":
            return CrossDistanceMatrix.from_file(path)
        raise ValueError(f"Unknown matrix_kind {kind!r} in {path}")
    # No matrix_kind key => legacy/symmetric format.
    return SymmetricDistanceMatrix.from_file(path)


class ClusteringResult:
    """Base clustering result: per-item labels and grouped clusters.

    Results are read-only. ``labels`` is a length-n NumPy ``intp`` array where
    ``-1`` marks noise; ``clusters`` is a tuple of member-index tuples.
    """

    def __init__(self, labels, clusters, *, native_owner=None):
        """
        Construct a clustering result.

        :param labels: Per-item integer labels. Noise is labeled -1.
        :param clusters: Iterable of cluster member iterables.
        :param native_owner: Native C++ result object that owns borrowed arrays
            (e.g. BitBirch centroids). Kept alive by this Python reference to
            prevent premature deallocation of zero-copy numpy views.
        """
        self._labels = np.asarray(list(labels), dtype=np.intp)
        self._clusters = tuple(
            tuple(int(member) for member in cluster) for cluster in clusters
        )
        # Held only to keep native-owned storage (e.g. BitBirch centroids) alive.
        self._native_owner = native_owner

    @property
    def labels(self):
        """Per-item integer labels as a NumPy ``intp`` array; -1 is noise."""
        return self._labels

    @property
    def clusters(self):
        """Tuple of clusters, each a tuple of member indices."""
        return self._clusters

    @property
    def method(self):
        """Clustering method that produced this result; '' if unknown."""
        return ""

    @property
    def num_clusters(self):
        """Number of clusters."""
        return len(self._clusters)

    @property
    def num_samples(self):
        """Number of items (length of ``labels``)."""
        return len(self._labels)

    def __len__(self):
        """Number of clusters."""
        return self.num_clusters

    def __iter__(self):
        """Iterate over clusters (each a tuple of member indices)."""
        return iter(self._clusters)

    def __getitem__(self, index):
        """Return the cluster at ``index``."""
        return self._clusters[index]

    def __repr__(self):
        return (f"{type(self).__name__}(num_clusters={self.num_clusters}, "
                f"num_samples={self.num_samples})")


class ButinaResult(ClusteringResult):
    """Butina clustering result.

    Butina clusters in descending order by representative neighborhood size.
    The first member of each cluster is the highest-neighborhood representative,
    selected as the molecule with the largest number of unassigned neighbors at
    the time it was chosen.
    """

    @property
    def method(self):
        return "butina"


class DBSCANResult(ClusteringResult):
    """DBSCAN clustering result with core sample indices.

    Core samples have at least ``min_samples`` neighbors within ``eps`` and
    form the skeletons of clusters. Border points (non-core members) are
    assigned to clusters via core samples; noise points have no core neighbor
    within eps.
    """

    def __init__(self, labels, clusters, *, core_sample_indices=(),
                 native_owner=None):
        super().__init__(labels, clusters, native_owner=native_owner)
        self._core_sample_indices = tuple(int(i) for i in core_sample_indices)

    @property
    def core_sample_indices(self):
        """Tuple of core sample indices (items with >= min_samples neighbors)."""
        return self._core_sample_indices

    @property
    def method(self):
        return "dbscan"


class HDBSCANResult(ClusteringResult):
    """HDBSCAN clustering result with membership probabilities.

    Probabilities reflect the stability of each item's cluster assignment,
    ranging from 0 (weakly assigned / noise) to 1 (strongly assigned core
    member). Noise points typically have probabilities near 0.
    """

    def __init__(self, labels, clusters, *, probabilities=None,
                 native_owner=None):
        super().__init__(labels, clusters, native_owner=native_owner)
        self._probabilities = (
            None if probabilities is None
            else np.asarray(list(probabilities), dtype=np.float64)
        )

    @property
    def probabilities(self):
        """Per-item cluster membership probabilities (0=noise, 1=core)."""
        return self._probabilities

    @property
    def method(self):
        return "hdbscan"


class AgglomerativeResult(ClusteringResult):
    """Agglomerative clustering result with merge-tree metadata.

    The merge tree encodes the hierarchical dendrogram: ``children[i]`` holds
    the ``(left, right)`` child-node indices merged at step ``i``, with
    ``distances[i]`` as the linkage distance and ``cluster_sizes[i]`` as the
    resulting cluster size. Node indices < n are leaf samples; indices >= n
    are internal merge nodes.
    """

    def __init__(self, labels, clusters, *, children=(), distances=None,
                 cluster_sizes=(), native_owner=None):
        super().__init__(labels, clusters, native_owner=native_owner)
        self._children = tuple(
            (int(left), int(right)) for left, right in children
        )
        self._distances = (
            None if distances is None
            else np.asarray(list(distances), dtype=np.float64)
        )
        self._cluster_sizes = tuple(int(size) for size in cluster_sizes)

    @property
    def children(self):
        """Merge tree: tuple of ``(left, right)`` child-node indices per merge."""
        return self._children

    @property
    def distances(self):
        """Linkage distance per merge, as a NumPy ``float64`` array."""
        return self._distances

    @property
    def cluster_sizes(self):
        """Cluster size after each merge, as a tuple of ints."""
        return self._cluster_sizes

    @property
    def method(self):
        return "agglomerative"


class BitBirchResult(ClusteringResult):
    """BitBirch clustering result with centroid fingerprints.

    Centroids are arithmetic-mean binary fingerprints computed by averaging
    the binary vectors of each cluster's members. They can be used for
    representative selection, reclustering, or refinement without reloading
    the input.
    """

    def __init__(self, labels, clusters, *, centroids=None, cluster_sizes=(),
                 native_owner=None):
        super().__init__(labels, clusters, native_owner=native_owner)
        self._centroids = centroids
        self._cluster_sizes = tuple(int(size) for size in cluster_sizes)

    @property
    def centroids(self):
        """Cluster centroid fingerprints as an ``oefp.OEFPBatch``."""
        return self._centroids

    @property
    def cluster_sizes(self):
        """Tuple of cluster sizes (member counts), aligned with centroids."""
        return self._cluster_sizes

    @property
    def method(self):
        return "bitbirch"




def pdist(items,
          comparison,
          *,
          similarity=False,
          num_threads=0,
          chunk_size=256,
          cutoff=0.0,
          output=None,
          progress=None,
          **kwargs) -> "SymmetricDistanceMatrix":
    """
    Compute pairwise distances for a collection of items using a comparison.

    :param items: List of molecules, design units, or other items.
    :param comparison: Comparison method: "fingerprint", "rocs", "superpose",
                       "sitehopper", "descriptor", or a C++ comparison object.
    :param similarity: Return similarities instead of distances.
    :param num_threads: Number of threads (0 = auto).
    :param chunk_size: Pairs per work unit.
    :param cutoff: Distance cutoff for sparse storage (0 = store all).
    :param output: Optional file path for memory-mapped storage.
    :param progress: Optional callback(completed, total).
    :param kwargs: Comparison-specific options.
    :returns: SymmetricDistanceMatrix with computed distances/similarities.
    :raises TypeError: If unknown kwargs are passed.
    """
    if isinstance(comparison, str):
        items, excluded = _comparisons.normalize_items(
            comparison, items, kwargs)
        # Refuse only when normalization is what emptied the list; an input
        # that arrived empty is left for the comparison to answer.
        if excluded and len(items) == 0:
            raise ValueError(
                "pdist requires a non-empty input set, but normalizing the "
                f"inputs for the {comparison!r} comparison left 0 items")
        labels = _comparisons.extract_labels(items)
        comparison_obj, comparison_name, params = _comparisons.build_comparison(
            items, comparison, similarity, kwargs, symmetric=True)
        if excluded:
            params['excluded_items'] = excluded
    else:
        comparison_obj = comparison
        comparison_name = comparison_obj.ComparisonName()
        labels = []
        params = {}

    n = comparison_obj.Size()

    storage: Any
    if output is not None:
        storage = MMapStorage(output, n)
    elif cutoff > 0.0:
        storage = SparseStorage(n, cutoff)
    else:
        storage = DenseStorage(n)

    options = PDistOptions()
    options.num_threads = num_threads
    options.chunk_size = chunk_size
    options.cutoff = cutoff
    if progress is not None:
        options.progress = progress

    _cpp_pdist(comparison_obj, storage, options)
    # After the computation, never before it: a comparison can only discover a
    # non-finite distance while it scores pairs, so an earlier read would stamp
    # ``complete`` on a matrix that turned out to contain NaN.
    facts = _gate.facts_from_comparison(comparison_obj)
    return SymmetricDistanceMatrix(storage, comparison_name, labels, params,
                                   facts)


def cdist(items_a, items_b, comparison, *,
          similarity=False,
          num_threads=0,
          chunk_size=256,
          cutoff=0.0,
          progress=None,
          **kwargs) -> "CrossDistanceMatrix":
    """
    Compute cross-distances between two item sets using a comparison.

    Entry ``[i, j]`` of the returned matrix compares ``items_a[i]`` (reference)
    against ``items_b[j]`` (fit). For asymmetric comparisons (e.g. Superpose with
    distinct reference/fit predicates), ``cdist(A, B)`` is not guaranteed to equal
    ``cdist(B, A).T``.

    :param items_a: Reference items (rows of the result).
    :param items_b: Fit items (columns of the result).
    :param comparison: Comparison method name: "fingerprint", "rocs", "superpose",
                       "sitehopper", or "descriptor". Prebuilt comparison objects
                       are not supported.
    :param similarity: Return similarities instead of distances.
    :param num_threads: Number of threads (0 = auto).
    :param chunk_size: Pairs per work unit.
    :param cutoff: Distance cutoff (0 = no cutoff); values above are zeroed.
    :param progress: Optional callback(completed, total).
    :param kwargs: Comparison-specific options.
    :returns: CrossDistanceMatrix of shape (len(items_a), len(items_b)).
    :raises TypeError: If a prebuilt comparison object is passed, or unknown kwargs.
    :raises ValueError: If either input is empty, or cutoff > 0 with similarity=True.
    """
    if not isinstance(comparison, str):
        raise TypeError(
            "cdist requires a string comparison name (e.g. 'fingerprint'); "
            "prebuilt comparison objects are not supported because cdist must "
            "construct the combined A+B comparison itself")

    if cutoff > 0.0 and similarity:
        raise ValueError(
            "cutoff > 0 is not supported with similarity=True: the cutoff zeroes "
            "values above the threshold, which would discard high similarities")

    # Materialize to lists (SWIG comparison constructors require Python lists) and
    # validate non-empty before allocating or calling into C++.
    a = list(items_a)
    b = list(items_b)
    if not a or not b:
        raise ValueError("cdist requires non-empty input sets (set A or B is empty)")

    a, excluded_a = _comparisons.normalize_items(comparison, a, kwargs)
    b, excluded_b = _comparisons.normalize_items(comparison, b, kwargs)
    n_a = len(a)
    n_b = len(b)
    if n_a == 0 or n_b == 0:
        raise ValueError(
            "cdist requires non-empty input sets, but normalizing the inputs "
            f"for the {comparison!r} comparison left {n_a} item(s) in set A "
            f"and {n_b} in set B")

    comparison_obj, comparison_name, params = _comparisons.build_comparison(
        a + b, comparison, similarity, kwargs, symmetric=False)
    if excluded_a:
        params['excluded_items_a'] = excluded_a
    if excluded_b:
        params['excluded_items_b'] = excluded_b

    # Guard the raw-pointer write with an explicit runtime check (not assert,
    # which -O strips).
    if comparison_obj.Size() != n_a + n_b:
        raise RuntimeError(
            f"Combined comparison size {comparison_obj.Size()} != "
            f"n_a + n_b ({n_a + n_b})")

    output = np.empty((n_a, n_b), dtype=np.float64)

    options = CDistOptions()
    options.num_threads = num_threads
    options.chunk_size = chunk_size
    options.cutoff = cutoff
    if progress is not None:
        options.progress = progress

    _cpp_cdist_into_address(comparison_obj, n_a, output.ctypes.data, options)
    facts = _gate.facts_from_comparison(comparison_obj)

    return CrossDistanceMatrix(
        output, comparison_name,
        labels_a=_comparisons.extract_labels(a),
        labels_b=_comparisons.extract_labels(b), params=params, facts=facts)


def butina(distance_matrix, threshold, *, reordering=False,
           num_threads=0, chunk_size=4096, allow_nonmetric=False):
    """
    Cluster a precomputed distance matrix using the Butina algorithm.

    :param distance_matrix: SymmetricDistanceMatrix returned by :func:`pdist`.
    :param threshold: Maximum distance for two items to be neighbors.
    :param reordering: Recompute candidate neighbor counts after each cluster.
    :param num_threads: Thread count for threshold graph construction.
    :param chunk_size: Condensed-distance pairs per work unit.
    :param allow_nonmetric: Cluster anyway when the distances are known not to
        satisfy the triangle inequality. Does not override the refusals for
        similarity-valued or non-finite matrices.
    :returns: ButinaResult with per-item labels and grouped clusters. The
        first member of each cluster is the highest-neighborhood representative,
        and each member's label equals its cluster index.
    :raises TypeError: If distance_matrix is not a SymmetricDistanceMatrix, or
        allow_nonmetric is not a bool.
    :raises ValueError: If threshold or a size argument is negative, or the
        matrix is not a metric.
    """
    if threshold < 0.0:
        raise ValueError("Butina threshold must be non-negative")
    if not isinstance(distance_matrix, SymmetricDistanceMatrix):
        raise TypeError("butina() expects a SymmetricDistanceMatrix")

    # Coerce caller arguments before the gate so that an invalid type is reported
    # ahead of an advisory refusal that names a remedy which cannot rescue it.
    num_threads_int = int(num_threads)
    chunk_size_int = int(chunk_size)

    # Both reach size_t option fields, where a negative value raises
    # OverflowError below the gate. Zero stays legal: it means "choose for me".
    if num_threads_int < 0:
        raise ValueError("num_threads must be non-negative")
    if chunk_size_int < 0:
        raise ValueError("chunk_size must be non-negative")

    _gate.require_metric(distance_matrix, "butina",
                         allow_nonmetric=allow_nonmetric)

    options = ButinaOptions()
    options.distance_threshold = float(threshold)
    options.reordering = bool(reordering)
    options.num_threads = num_threads_int
    options.chunk_size = chunk_size_int

    result = _butina_cluster(distance_matrix.storage, options)
    return ButinaResult(result.Labels(), result.Members())


class RepresentativeMetrics:
    """Quality metrics for one cluster representative.

    :ivar mean_distance_to_cluster: Mean distance to cluster members.
    :ivar max_distance_to_cluster: Maximum distance to any cluster member.
    :ivar median_distance_to_cluster: Median distance to cluster members.
    :ivar neighbor_fraction_at_threshold: Fraction of cluster members within
        the threshold distance.
    :ivar nearest_external_distance: Distance to the nearest non-cluster item.
    :ivar cluster_radius: Maximum distance from the representative to any
        cluster member (same as max_distance_to_cluster).
    :ivar cluster_diameter: Maximum pairwise distance within the cluster.
    :ivar silhouette_like_score: Silhouette-like separation score
        (higher = better separation).
    :ivar scaffold_purity: Fraction of cluster members sharing the
        representative's scaffold (when scaffold labels are provided).
    :ivar representative_rank: Zero-based rank of this representative by score
        (0 = best).
    """

    __slots__ = (
        "mean_distance_to_cluster",
        "max_distance_to_cluster",
        "median_distance_to_cluster",
        "neighbor_fraction_at_threshold",
        "nearest_external_distance",
        "cluster_radius",
        "cluster_diameter",
        "silhouette_like_score",
        "scaffold_purity",
        "representative_rank",
    )

    def __init__(self, native_metrics):
        """
        Construct metrics from the native representative result.

        :param native_metrics: Native metrics object returned by the extension.
        """
        self.mean_distance_to_cluster = float(
            native_metrics.mean_distance_to_cluster)
        self.max_distance_to_cluster = float(native_metrics.max_distance_to_cluster)
        self.median_distance_to_cluster = float(
            native_metrics.median_distance_to_cluster)
        self.neighbor_fraction_at_threshold = float(
            native_metrics.neighbor_fraction_at_threshold)
        self.nearest_external_distance = float(
            native_metrics.nearest_external_distance)
        self.cluster_radius = float(native_metrics.cluster_radius)
        self.cluster_diameter = float(native_metrics.cluster_diameter)
        self.silhouette_like_score = float(native_metrics.silhouette_like_score)
        self.scaffold_purity = float(native_metrics.scaffold_purity)
        self.representative_rank = int(native_metrics.representative_rank)


class ClusterRepresentative:
    """A scored cluster representative and its quality metrics.

    :ivar member: Item index of the representative (int).
    :ivar score: Representative score (float). Lower is better for
        distance-based methods; interpretation depends on the scoring method.
    :ivar metrics: Quality metrics as a :class:`RepresentativeMetrics` instance.
    """

    __slots__ = ("member", "score", "metrics")

    def __init__(self, native_representative):
        """
        Construct a representative from the native result.

        :param native_representative: Native representative object.
        """
        self.member = int(native_representative.member)
        self.score = float(native_representative.score)
        self.metrics = RepresentativeMetrics(native_representative.metrics)


def _cpp_cluster(cluster, function_name):
    cpp_cluster = _oecluster.SizeTVector()
    for member in cluster:
        cpp_cluster.push_back(int(member))
    if len(cpp_cluster) == 0:
        raise ValueError(f"{function_name}() requires a non-empty cluster")
    return cpp_cluster


def _representative_method(method):
    method_map = {
        "medoid": _oecluster.RepresentativeMethod_Medoid,
        "minimax": _oecluster.RepresentativeMethod_Minimax,
        "highest_neighborhood": _oecluster.RepresentativeMethod_HighestNeighborhood,
        "weighted_medoid": _oecluster.RepresentativeMethod_WeightedMedoid,
    }
    method_key = str(method).lower()
    if method_key not in method_map:
        raise ValueError(f"Unknown representative method: {method!r}")
    return method_key, method_map[method_key]


def _representative_selection(selection):
    selection_map = {
        "score": _oecluster.RepresentativeSelection_Score,
        "diversity": _oecluster.RepresentativeSelection_Diversity,
    }
    selection_key = str(selection).lower()
    if selection_key not in selection_map:
        raise ValueError(f"Unknown representative selection: {selection!r}")
    return selection_map[selection_key]


def _optional_float_vector(values, name, distance_matrix):
    vector = _oecluster.DoubleVector()
    if values is None:
        return vector
    converted = [float(value) for value in values]
    if converted and len(converted) != distance_matrix.num_samples:
        raise ValueError(
            f"{name} must be empty or the same length as the distance matrix")
    for value in converted:
        vector.push_back(value)
    return vector


def _optional_string_vector(values, name, distance_matrix):
    vector = _oecluster.StringVector()
    if values is None:
        return vector
    converted = [str(value) for value in values]
    if converted and len(converted) != distance_matrix.num_samples:
        raise ValueError(
            f"{name} must be empty or the same length as the distance matrix")
    for value in converted:
        vector.push_back(value)
    return vector


def _representative_options(
    distance_matrix,
    *,
    method,
    threshold,
    selection="score",
    alpha=1.0,
    beta=1.0,
    gamma=1.0,
    liability_penalties=None,
    priority_scores=None,
    scaffold_labels=None,
):
    if not isinstance(distance_matrix, SymmetricDistanceMatrix):
        raise TypeError("representative functions expect a SymmetricDistanceMatrix")

    method_key, native_method = _representative_method(method)
    if threshold is None:
        threshold_value = -1.0
    else:
        threshold_value = float(threshold)
        if threshold_value < 0.0:
            raise ValueError("Representative threshold must be non-negative")
    if method_key == "highest_neighborhood" and threshold is None:
        raise ValueError(
            "highest_neighborhood representative requires a threshold")

    weights = RepresentativeWeights()
    weights.alpha = float(alpha)
    weights.beta = float(beta)
    weights.gamma = float(gamma)

    options = RepresentativeOptions()
    options.method = native_method
    options.selection = _representative_selection(selection)
    options.neighbor_threshold = threshold_value
    options.weights = weights
    options.liability_penalties = _optional_float_vector(
        liability_penalties,
        "liability_penalties",
        distance_matrix,
    )
    options.priority_scores = _optional_float_vector(
        priority_scores,
        "priority_scores",
        distance_matrix,
    )
    options.scaffold_labels = _optional_string_vector(
        scaffold_labels,
        "scaffold_labels",
        distance_matrix,
    )
    return options


def _wrap_representatives(native_representatives):
    return tuple(
        ClusterRepresentative(representative)
        for representative in native_representatives
    )


def representative(cluster, distance_matrix, *, method="medoid", threshold=None,
                   alpha=1.0, beta=1.0, gamma=1.0,
                   liability_penalties=None, priority_scores=None,
                   scaffold_labels=None):
    """
    Select the best representative member from a cluster.

    :param cluster: Iterable of item indices, such as one cluster returned by
        :func:`butina`.
    :param distance_matrix: SymmetricDistanceMatrix used to compute representative scores.
    :param method: Scoring method: "medoid", "minimax",
        "highest_neighborhood", or "weighted_medoid".
    :param threshold: Distance threshold for highest-neighborhood scoring and
        neighbor-fraction metrics.
    :param alpha: Weight applied to mean distance for weighted medoids.
    :param beta: Weight applied to liability penalties for weighted medoids.
    :param gamma: Weight applied to priority scores for weighted medoids.
    :param liability_penalties: Optional per-item penalty vector.
    :param priority_scores: Optional per-item priority vector.
    :param scaffold_labels: Optional per-item scaffold label vector.
    :returns: Selected member index.
    :raises TypeError: If distance_matrix is not a SymmetricDistanceMatrix.
    :raises ValueError: If cluster, method, threshold, or metadata is invalid.
    """
    cpp_cluster = _cpp_cluster(cluster, "representative")
    options = _representative_options(
        distance_matrix,
        method=method,
        threshold=threshold,
        alpha=alpha,
        beta=beta,
        gamma=gamma,
        liability_penalties=liability_penalties,
        priority_scores=priority_scores,
        scaffold_labels=scaffold_labels,
    )
    return int(_cluster_representative(cpp_cluster, distance_matrix.storage, options))


def rank_representatives(cluster, distance_matrix, *, method="medoid",
                         threshold=None, alpha=1.0, beta=1.0, gamma=1.0,
                         liability_penalties=None, priority_scores=None,
                         scaffold_labels=None):
    """
    Rank all cluster members as representatives.

    :param cluster: Iterable of item indices.
    :param distance_matrix: SymmetricDistanceMatrix used to compute representative scores.
    :param method: Scoring method: "medoid", "minimax",
        "highest_neighborhood", or "weighted_medoid".
    :param threshold: Distance threshold for highest-neighborhood scoring and
        neighbor-fraction metrics.
    :param alpha: Weight applied to mean distance for weighted medoids.
    :param beta: Weight applied to liability penalties for weighted medoids.
    :param gamma: Weight applied to priority scores for weighted medoids.
    :param liability_penalties: Optional per-item penalty vector.
    :param priority_scores: Optional per-item priority vector.
    :param scaffold_labels: Optional per-item scaffold label vector.
    :returns: Tuple of ClusterRepresentative objects sorted by score.
    :raises TypeError: If distance_matrix is not a SymmetricDistanceMatrix.
    :raises ValueError: If cluster, method, threshold, or metadata is invalid.
    """
    cpp_cluster = _cpp_cluster(cluster, "rank_representatives")
    options = _representative_options(
        distance_matrix,
        method=method,
        threshold=threshold,
        alpha=alpha,
        beta=beta,
        gamma=gamma,
        liability_penalties=liability_penalties,
        priority_scores=priority_scores,
        scaffold_labels=scaffold_labels,
    )
    return _wrap_representatives(
        _rank_representatives(cpp_cluster, distance_matrix.storage, options))


def select_representatives(cluster, distance_matrix, *, k, method="medoid",
                           selection="score", threshold=None, alpha=1.0,
                           beta=1.0, gamma=1.0, liability_penalties=None,
                           priority_scores=None, scaffold_labels=None):
    """
    Select up to k representatives from a cluster.

    :param cluster: Iterable of item indices.
    :param distance_matrix: SymmetricDistanceMatrix used to compute representative scores.
    :param k: Maximum number of representatives to return.
    :param method: Scoring method: "medoid", "minimax",
        "highest_neighborhood", or "weighted_medoid".
    :param selection: "score" for top-k ranking or "diversity" for greedy
        coverage after the top-scoring representative.
    :param threshold: Distance threshold for highest-neighborhood scoring and
        neighbor-fraction metrics.
    :param alpha: Weight applied to mean distance for weighted medoids.
    :param beta: Weight applied to liability penalties for weighted medoids.
    :param gamma: Weight applied to priority scores for weighted medoids.
    :param liability_penalties: Optional per-item penalty vector.
    :param priority_scores: Optional per-item priority vector.
    :param scaffold_labels: Optional per-item scaffold label vector.
    :returns: Tuple of selected ClusterRepresentative objects.
    :raises TypeError: If distance_matrix is not a SymmetricDistanceMatrix.
    :raises ValueError: If cluster, k, method, selection, threshold, or metadata
        is invalid.
    """
    if k < 0:
        raise ValueError("k must be non-negative")

    cpp_cluster = _cpp_cluster(cluster, "select_representatives")
    options = _representative_options(
        distance_matrix,
        method=method,
        threshold=threshold,
        selection=selection,
        alpha=alpha,
        beta=beta,
        gamma=gamma,
        liability_penalties=liability_penalties,
        priority_scores=priority_scores,
        scaffold_labels=scaffold_labels,
    )
    return _wrap_representatives(
        _select_representatives(
            cpp_cluster,
            distance_matrix.storage,
            int(k),
            options,
        ))


def dbscan(distance_matrix, eps, *, min_samples=5, num_threads=0,
           chunk_size=4096, allow_nonmetric=False):
    """
    Cluster a precomputed distance matrix using DBSCAN.

    :param distance_matrix: SymmetricDistanceMatrix returned by :func:`pdist`.
    :param eps: Maximum distance for two items to be neighbors.
    :param min_samples: Minimum self-inclusive neighbor count for a core sample.
    :param num_threads: Thread count for threshold graph construction.
    :param chunk_size: Condensed-distance pairs per work unit.
    :param allow_nonmetric: Cluster anyway when the distances are known not to
        satisfy the triangle inequality. Does not override the refusals for
        similarity-valued or non-finite matrices.
    :returns: DBSCANResult with labels, clusters, and core sample indices.
    :raises TypeError: If distance_matrix is not a SymmetricDistanceMatrix, or
        allow_nonmetric is not a bool.
    :raises ValueError: If eps, min_samples, or a size argument is invalid, or
        the matrix is not a metric.
    """
    if eps < 0.0:
        raise ValueError("DBSCAN eps must be non-negative")
    if min_samples < 1:
        raise ValueError("DBSCAN min_samples must be at least one")
    if not isinstance(distance_matrix, SymmetricDistanceMatrix):
        raise TypeError("dbscan() expects a SymmetricDistanceMatrix")

    # Coerce caller arguments before the gate so that an invalid type is reported
    # ahead of an advisory refusal that names a remedy which cannot rescue it.
    num_threads_int = int(num_threads)
    chunk_size_int = int(chunk_size)

    # Both reach size_t option fields, where a negative value raises
    # OverflowError below the gate. Zero stays legal: it means "choose for me".
    if num_threads_int < 0:
        raise ValueError("num_threads must be non-negative")
    if chunk_size_int < 0:
        raise ValueError("chunk_size must be non-negative")

    _gate.require_metric(distance_matrix, "dbscan",
                         allow_nonmetric=allow_nonmetric)

    options = DBSCANOptions()
    options.eps = float(eps)
    options.min_samples = int(min_samples)
    options.num_threads = num_threads_int
    options.chunk_size = chunk_size_int

    result = _dbscan_cluster(distance_matrix.storage, options)
    return DBSCANResult(
        result.Labels(),
        result.Members(),
        core_sample_indices=result.CoreSampleIndices(),
    )


def hdbscan(distance_matrix, *, min_cluster_size=5, min_samples=None,
            cluster_selection_epsilon=0.0, max_cluster_size=None, alpha=1.0,
            cluster_selection_method="eom", allow_single_cluster=False,
            num_threads=0, chunk_size=4096, allow_nonmetric=False):
    """
    Cluster a precomputed distance matrix using HDBSCAN.

    :param distance_matrix: SymmetricDistanceMatrix returned by :func:`pdist`.
    :param min_cluster_size: Minimum size for selected clusters.
    :param min_samples: Self-inclusive core-distance neighbor count. Defaults
        to min_cluster_size when omitted.
    :param cluster_selection_epsilon: Epsilon threshold for merging selected
        clusters.
    :param max_cluster_size: Optional maximum selected cluster size.
    :param alpha: Mutual-reachability distance scaling.
    :param cluster_selection_method: Cluster selection method, "eom" or "leaf".
    :param allow_single_cluster: Whether the root cluster may be selected.
    :param num_threads: Thread count for core-distance computation.
    :param chunk_size: Reserved for parity with other clustering wrappers.
    :param allow_nonmetric: Cluster anyway when the distances are known not to
        satisfy the triangle inequality. Does not override the refusals for
        similarity-valued or non-finite matrices.
    :returns: HDBSCANResult with labels, clusters, and probabilities.
    :raises TypeError: If distance_matrix is not a SymmetricDistanceMatrix, or
        allow_nonmetric is not a bool.
    :raises ValueError: If options are invalid, the matrix uses sparse storage,
        or the matrix is not a metric.
    """
    if min_cluster_size < 2:
        raise ValueError("HDBSCAN min_cluster_size must be at least two")
    if min_samples is not None and min_samples < 1:
        raise ValueError("HDBSCAN min_samples must be at least one")
    if cluster_selection_epsilon < 0.0:
        raise ValueError("HDBSCAN cluster_selection_epsilon must be non-negative")
    if max_cluster_size is not None and max_cluster_size < 1:
        raise ValueError("HDBSCAN max_cluster_size must be at least one")
    if alpha <= 0.0:
        raise ValueError("HDBSCAN alpha must be positive")
    if not isinstance(distance_matrix, SymmetricDistanceMatrix):
        raise TypeError("hdbscan() expects a SymmetricDistanceMatrix")

    method_map = {
        "eom": _oecluster.HDBSCANClusterSelectionMethod_EOM,
        "leaf": _oecluster.HDBSCANClusterSelectionMethod_Leaf,
    }
    method_key = str(cluster_selection_method).lower()
    if method_key not in method_map:
        raise ValueError(
            f"Unknown HDBSCAN cluster_selection_method: {cluster_selection_method!r}"
        )

    # Coerce caller arguments before the gate so that an invalid type is reported
    # ahead of an advisory refusal that names a remedy which cannot rescue it.
    num_threads_int = int(num_threads)
    chunk_size_int = int(chunk_size)

    # Both reach size_t option fields, where a negative value raises
    # OverflowError below the gate. Zero stays legal: it means "choose for me".
    if num_threads_int < 0:
        raise ValueError("num_threads must be non-negative")
    if chunk_size_int < 0:
        raise ValueError("chunk_size must be non-negative")

    # HDBSCAN bounds the min_samples it actually uses, which is min_cluster_size
    # when the caller leaves min_samples unset, so the mirror substitutes the
    # same way the native code does.
    effective_min_samples = (int(min_cluster_size) if min_samples is None
                             else int(min_samples))
    if effective_min_samples > distance_matrix.num_samples:
        raise ValueError(
            "HDBSCAN min_samples must be at most "
            f"the item count ({distance_matrix.num_samples})")
    # ValueError, not TypeError: the argument's type is right, its storage is not.
    if isinstance(distance_matrix.storage, SparseStorage):
        raise ValueError(  # noqa: TRY004
            "HDBSCAN requires complete pairwise "
            "distances; SparseStorage is not supported")

    # Local argument validation first: allow_nonmetric cannot rescue an unknown
    # cluster_selection_method, an out-of-range min_samples, or sparse storage,
    # so the gate must not pre-empt those messages. The bound and the storage
    # check run in the order hdbscan_cluster() applies them, so the fix-first
    # reason is the same whichever layer reports it.
    _gate.require_metric(distance_matrix, "hdbscan",
                         allow_nonmetric=allow_nonmetric)

    options = HDBSCANOptions()
    options.min_cluster_size = int(min_cluster_size)
    options.min_samples = 0 if min_samples is None else int(min_samples)
    options.cluster_selection_epsilon = float(cluster_selection_epsilon)
    options.max_cluster_size = 0 if max_cluster_size is None else int(max_cluster_size)
    options.alpha = float(alpha)
    options.cluster_selection_method = method_map[method_key]
    options.allow_single_cluster = bool(allow_single_cluster)
    options.num_threads = num_threads_int
    options.chunk_size = chunk_size_int

    result = _hdbscan_cluster(distance_matrix.storage, options)
    return HDBSCANResult(
        result.Labels(),
        result.Members(),
        probabilities=result.Probabilities(),
    )


def agglomerative(distance_matrix, *, n_clusters=2, distance_threshold=None,
                  linkage="average", compute_full_tree=True,
                  num_threads=0, chunk_size=4096, allow_nonmetric=False):
    """
    Cluster a precomputed distance matrix using hierarchical agglomerative clustering.

    :param distance_matrix: SymmetricDistanceMatrix returned by :func:`pdist`.
    :param n_clusters: Number of flat clusters when distance_threshold is omitted.
    :param distance_threshold: Optional merge-distance cutoff for flat clusters.
    :param linkage: Linkage method: "single", "complete", "average", or "weighted".
    :param compute_full_tree: Whether to request full-tree computation.
    :param num_threads: Thread count for initial distance materialization.
    :param chunk_size: Rows per work unit during distance materialization.
    :param allow_nonmetric: Cluster anyway when the distances are known not to
        satisfy the triangle inequality. Does not override the refusals for
        similarity-valued or non-finite matrices.
    :returns: AgglomerativeResult with labels, clusters, children, distances, and cluster sizes.
    :raises TypeError: If distance_matrix is not a SymmetricDistanceMatrix, or
        allow_nonmetric is not a bool.
    :raises ValueError: If options are invalid, the matrix uses sparse storage,
        or the matrix is not a metric.
    """
    if not isinstance(distance_matrix, SymmetricDistanceMatrix):
        raise TypeError("agglomerative() expects a SymmetricDistanceMatrix")
    if distance_threshold is None and n_clusters < 1:
        raise ValueError("Agglomerative n_clusters must be at least one")
    if distance_threshold is not None and distance_threshold < 0.0:
        raise ValueError("Agglomerative distance_threshold must be non-negative")
    if distance_threshold is not None and np.isnan(distance_threshold):
        raise ValueError("Agglomerative distance_threshold must not be NaN")

    linkage_map = {
        "single": _oecluster.AgglomerativeLinkageMethod_Single,
        "complete": _oecluster.AgglomerativeLinkageMethod_Complete,
        "average": _oecluster.AgglomerativeLinkageMethod_Average,
        "weighted": _oecluster.AgglomerativeLinkageMethod_Weighted,
    }
    linkage_key = str(linkage).lower()
    if linkage_key not in linkage_map:
        raise ValueError(f"Unknown agglomerative linkage: {linkage!r}")

    # Coerce caller arguments before the gate so that an invalid type is reported
    # ahead of an advisory refusal that names a remedy which cannot rescue it.
    n_clusters_int = int(n_clusters)
    num_threads_int = int(num_threads)
    chunk_size_int = int(chunk_size)

    # All three reach size_t option fields, where a negative value raises
    # OverflowError below the gate. n_clusters needs its own check because the
    # "at least one" guard above only runs when distance_threshold is omitted.
    # Zero stays legal for num_threads: it means "choose for me".
    if n_clusters_int < 0:
        raise ValueError("Agglomerative n_clusters must be non-negative")
    if num_threads_int < 0:
        raise ValueError("num_threads must be non-negative")
    if chunk_size_int < 0:
        raise ValueError("chunk_size must be non-negative")

    # ValueError, not TypeError: the argument's type is right, its storage is not.
    if isinstance(distance_matrix.storage, SparseStorage):
        raise ValueError(  # noqa: TRY004
            "Agglomerative clustering requires complete pairwise "
            "distances; SparseStorage is not supported")

    # Native rule: chunk_size == 0 is rejected. Zero is legal for butina/dbscan.
    if chunk_size_int == 0:
        raise ValueError("Agglomerative chunk_size must be at least one")

    # The native bound on n_clusters applies only when no distance_threshold is
    # given, because a threshold cut ignores n_clusters entirely. Mirroring it
    # unconditionally would refuse calls that work today.
    if (distance_threshold is None
            and n_clusters_int > distance_matrix.num_samples):
        raise ValueError(
            "Agglomerative n_clusters must be at most "
            f"the item count ({distance_matrix.num_samples})")

    # Local argument validation first: allow_nonmetric cannot rescue a bad
    # n_clusters, distance_threshold, linkage, or sparse storage, so the gate
    # must not pre-empt those messages. The storage, chunk_size, and n_clusters
    # checks run in the order agglomerative_cluster() applies them, so the
    # fix-first reason is the same whichever layer reports it.
    _gate.require_metric(distance_matrix, "agglomerative",
                         allow_nonmetric=allow_nonmetric)

    options = AgglomerativeOptions()
    options.n_clusters = n_clusters_int
    options.distance_threshold = (
        -1.0 if distance_threshold is None else float(distance_threshold)
    )
    options.linkage = linkage_map[linkage_key]
    options.compute_full_tree = bool(compute_full_tree)
    options.num_threads = num_threads_int
    options.chunk_size = chunk_size_int

    result = _agglomerative_cluster(distance_matrix.storage, options)
    children = zip(result.ChildrenLeft(), result.ChildrenRight())
    return AgglomerativeResult(
        result.Labels(),
        result.Members(),
        children=children,
        distances=result.Distances(),
        cluster_sizes=result.ClusterSizes(),
    )


def _bitbirch_merge_criterion(value):
    criterion_map = {
        "radius": _oecluster.BitBirchMergeCriterion_Radius,
        "diameter": _oecluster.BitBirchMergeCriterion_Diameter,
        "tolerance": _oecluster.BitBirchMergeCriterion_Tolerance,
        "tolerance_tough": _oecluster.BitBirchMergeCriterion_ToleranceTough,
    }
    key = str(value).lower()
    if key not in criterion_map:
        raise ValueError(f"Unknown BitBirch merge_criterion: {value!r}")
    return criterion_map[key]


def _bitbirch_mode(value):
    mode_map = {
        "strict_parity": _oecluster.BitBirchMode_StrictParity,
        "fast": _oecluster.BitBirchMode_Fast,
    }
    key = str(value).lower()
    if key not in mode_map:
        raise ValueError(f"Unknown BitBirch mode: {value!r}")
    return mode_map[key]


def _require_oefp_batch(fingerprints, function_name):
    try:
        import oefp as _oefp_api
    except ImportError as exc:
        raise TypeError(
            f"{function_name}() expects an oefp.OEFPBatch and oefp is not importable"
        ) from exc
    if not isinstance(fingerprints, _oefp_api.OEFPBatch):
        raise TypeError(f"{function_name}() expects an oefp.OEFPBatch")


def _bitbirch_centroids(native_centroids):
    import oefp as _oefp_api

    return _oefp_api.OEFPBatch._from_native(native_centroids)


def bitbirch(fingerprints, *, threshold=0.65, branching_factor=50,
             merge_criterion="diameter", tolerance=0.05, singly=True,
             mode="strict_parity", num_threads=0):
    """
    Cluster dense binary OEFP fingerprints using BitBirch.

    :param fingerprints: `oefp.OEFPBatch` containing dense binary fingerprints.
    :param threshold: Similarity threshold used by the merge criterion.
    :param branching_factor: Maximum number of subclusters per tree node.
    :param merge_criterion: "radius", "diameter", "tolerance", or
        "tolerance_tough".
    :param tolerance: Tolerance penalty for tolerance-based criteria.
    :param singly: Whether to skip parent-pointer maintenance for the faster
        single-pass reference behavior.
    :param mode: "strict_parity" (default, exact reference parity) or "fast".
        Fast partitions the input, fits chunk trees in parallel, and merges them;
        its output is deterministic and independent of num_threads, may differ
        from strict cluster shapes, and is byte-identical to strict for small
        inputs.
    :param num_threads: Thread count for parallel-safe phases.
    :returns: BitBirchResult with labels, clusters, centroids, and cluster sizes.
    :raises TypeError: If fingerprints is not an `oefp.OEFPBatch`.
    :raises ValueError: If an option is invalid.
    """
    _require_oefp_batch(fingerprints, "bitbirch")
    if threshold < 0.0:
        raise ValueError("BitBirch threshold must be non-negative")
    if branching_factor < 1:
        raise ValueError("BitBirch branching_factor must be at least one")
    if tolerance < 0.0:
        raise ValueError("BitBirch tolerance must be non-negative")

    options = BitBirchOptions()
    options.threshold = float(threshold)
    options.branching_factor = int(branching_factor)
    options.merge_criterion = _bitbirch_merge_criterion(merge_criterion)
    options.tolerance = float(tolerance)
    options.singly = bool(singly)
    options.mode = _bitbirch_mode(mode)
    options.num_threads = int(num_threads)

    result = _bitbirch_cluster(fingerprints, options)
    return BitBirchResult(
        result.Labels(),
        result.Members(),
        cluster_sizes=result.ClusterSizes(),
        centroids=_bitbirch_centroids(result.Centroids()),
        native_owner=result,
    )


def bitbirch_recluster(fingerprints, *, initial_threshold=0.65,
                       second_threshold=0.7, second_tolerance=0.0,
                       branching_factor=50, mode="strict_parity",
                       num_threads=0):
    """
    Cluster dense binary OEFP fingerprints using two-stage BitBirch reclustering.

    :param fingerprints: `oefp.OEFPBatch` containing dense binary fingerprints.
    :param initial_threshold: Diameter threshold for the first pass.
    :param second_threshold: Tolerance threshold for the second pass.
    :param second_tolerance: Tolerance penalty for the second pass.
    :param branching_factor: Maximum number of subclusters per tree node.
    :param mode: "strict_parity" (default, exact reference parity) or "fast".
        Fast partitions the input, fits chunk trees in parallel, and merges them;
        its output is deterministic and independent of num_threads, may differ
        from strict cluster shapes, and is byte-identical to strict for small
        inputs.
    :param num_threads: Thread count for parallel-safe phases.
    :returns: BitBirchResult with labels, clusters, centroids, and cluster sizes.
    """
    _require_oefp_batch(fingerprints, "bitbirch_recluster")
    if initial_threshold < 0.0 or second_threshold < 0.0:
        raise ValueError("BitBirch reclustering thresholds must be non-negative")
    if second_tolerance < 0.0:
        raise ValueError("BitBirch reclustering tolerance must be non-negative")
    if branching_factor < 1:
        raise ValueError("BitBirch branching_factor must be at least one")

    options = BitBirchReclusteringOptions()
    options.initial_threshold = float(initial_threshold)
    options.second_threshold = float(second_threshold)
    options.second_tolerance = float(second_tolerance)
    options.branching_factor = int(branching_factor)
    options.mode = _bitbirch_mode(mode)
    options.num_threads = int(num_threads)

    result = _bitbirch_recluster(fingerprints, options)
    return BitBirchResult(
        result.Labels(),
        result.Members(),
        cluster_sizes=result.ClusterSizes(),
        centroids=_bitbirch_centroids(result.Centroids()),
        native_owner=result,
    )


def bitbirch_refine(fingerprints, *, threshold=0.65, branching_factor=50,
                    merge_criterion="diameter", tolerance=0.05, singly=False,
                    redistribute_largest_cluster=False, reassign_top_clusters=0,
                    mode="strict_parity", num_threads=0):
    """
    Fit BitBirch and apply requested refinement passes.

    :param fingerprints: `oefp.OEFPBatch` containing dense binary fingerprints.
    :param threshold: Similarity threshold used by the merge criterion.
    :param branching_factor: Maximum number of subclusters per tree node.
    :param merge_criterion: "radius", "diameter", "tolerance", or
        "tolerance_tough".
    :param tolerance: Tolerance penalty for tolerance-based criteria.
    :param singly: Whether to skip parent-pointer maintenance during the fit.
    :param redistribute_largest_cluster: Whether to redistribute molecules from
        the largest cluster after fitting.
    :param reassign_top_clusters: Number of largest clusters to reassign by
        centroid similarity. Use zero to disable reassignment.
    :param mode: "strict_parity" or "fast"; refinement always runs the
        strict-parity path. The fast partition-merge path is not applied to
        refine because its prune/reassign passes are order-sensitive.
    :param num_threads: Thread count for parallel-safe refinement scoring.
    :returns: BitBirchResult with labels, clusters, centroids, and cluster sizes.
    """
    _require_oefp_batch(fingerprints, "bitbirch_refine")
    if threshold < 0.0:
        raise ValueError("BitBirch threshold must be non-negative")
    if branching_factor < 1:
        raise ValueError("BitBirch branching_factor must be at least one")
    if tolerance < 0.0:
        raise ValueError("BitBirch tolerance must be non-negative")
    if reassign_top_clusters == 1:
        raise ValueError("BitBirch reassign_top_clusters must be zero or at least two")
    if redistribute_largest_cluster and singly:
        raise ValueError("BitBirch redistribute_largest_cluster requires singly=False")

    fit_options = BitBirchOptions()
    fit_options.threshold = float(threshold)
    fit_options.branching_factor = int(branching_factor)
    fit_options.merge_criterion = _bitbirch_merge_criterion(merge_criterion)
    fit_options.tolerance = float(tolerance)
    fit_options.singly = bool(singly)
    fit_options.mode = _bitbirch_mode(mode)
    fit_options.num_threads = int(num_threads)

    options = BitBirchRefinementOptions()
    options.fit_options = fit_options
    options.redistribute_largest_cluster = bool(redistribute_largest_cluster)
    options.reassign_top_clusters = int(reassign_top_clusters)
    options.num_threads = int(num_threads)

    result = _bitbirch_refine(fingerprints, options)
    return BitBirchResult(
        result.Labels(),
        result.Members(),
        cluster_sizes=result.ClusterSizes(),
        centroids=_bitbirch_centroids(result.Centroids()),
        native_owner=result,
    )


def _native_clustering_result(result):
    """Rebuild a native ClusteringResult from a Python result's labels/clusters."""
    label_vec = _oecluster.IntVector()
    for label in result.labels:
        label_vec.push_back(int(label))
    member_vec = _oecluster.ClusterVector()
    for cluster in result.clusters:
        inner = _oecluster.SizeTVector()
        for member in cluster:
            inner.push_back(int(member))
        member_vec.push_back(inner)
    return _oecluster.ClusteringResult(label_vec, member_vec)


class ClusterReport:
    """Read-only clustering-quality scorecard.

    Provides 20 scalar quality metrics plus threshold-coverage curves for
    evaluating clustering results. Key metrics include silhouette (separation),
    Dunn index (compactness vs. isolation), size Gini coefficient (imbalance),
    median radius/diameter, and coverage fractions at user-specified thresholds.

    Users should construct reports via the :func:`cluster_report` function,
    not by calling ``__init__`` directly.

    All metrics are exposed as read-only properties. Undefined metrics are NaN.
    Vector metrics (``coverage_thresholds``, ``coverage_at``) are tuples.
    """

    _SCALAR_FIELDS = (
        "num_samples", "num_clusters", "num_noise", "num_singletons",
        "noise_fraction", "singleton_fraction", "largest_cluster_fraction",
        "cluster_size_median", "cluster_size_p90", "size_gini", "size_entropy",
        "mean_intra_distance", "median_intra_distance", "median_radius",
        "p95_diameter", "silhouette", "dunn_index", "boundary_violations",
        "median_medoid_member_distance", "representative_redundancy",
    )

    def __init__(self, native_report, method=""):
        """Capture every field from the native report into Python values.

        Internal constructor; users should call :func:`cluster_report` instead.

        :param native_report: Native C++ ClusterReport object.
        :param method: Clustering method name.
        """
        for name in self._SCALAR_FIELDS:
            object.__setattr__(self, f"_{name}", getattr(native_report, name))
        object.__setattr__(
            self, "_coverage_thresholds",
            tuple(float(v) for v in native_report.coverage_thresholds))
        object.__setattr__(
            self, "_coverage_at",
            tuple(float(v) for v in native_report.coverage_at))
        object.__setattr__(self, "_method", str(method))

    def __setattr__(self, name, value):
        """ClusterReport is read-only; reject external attribute assignment."""
        raise AttributeError(
            f"ClusterReport is read-only; cannot set {name!r}")

    def __getattr__(self, name):
        if name in ClusterReport._SCALAR_FIELDS:
            return object.__getattribute__(self, f"_{name}")
        raise AttributeError(name)

    @property
    def method(self):
        """Clustering method that produced the result this report describes."""
        return self._method

    @property
    def coverage_thresholds(self):
        """Distance thresholds for coverage (tuple of floats)."""
        return self._coverage_thresholds

    @property
    def coverage_at(self):
        """Coverage fraction at each threshold (tuple, aligned with coverage_thresholds)."""
        return self._coverage_at

    def __repr__(self):
        return (f"ClusterReport(method={self.method!r}, "
                f"num_clusters={self.num_clusters}, "
                f"num_samples={self.num_samples}, "
                f"silhouette={self.silhouette:.4f})")


class ClusterReportComparison:
    """Two or more ClusterReports aligned for side-by-side reading."""

    def __init__(self, reports):
        self._reports = tuple(reports)

    @property
    def reports(self):
        """The compared reports, in input order (tuple)."""
        return self._reports

    def _column_labels(self):
        """Per-report column labels: the method name, or report{i} if empty."""
        labels = []
        for index, report in enumerate(self._reports):
            labels.append(report.method if report.method else f"report{index}")
        return labels

    def to_table(self):
        """Return rows ``(metric_name, *values)`` -- one value per report.

        Scalar-metric rows come first, then coverage rows aligned by threshold
        value across all reports; a report lacking a given threshold shows NaN.
        """
        rows = []
        for name in ClusterReport._SCALAR_FIELDS:
            rows.append((name, *(getattr(r, name) for r in self._reports)))

        def _coverage_lookup(report, threshold):
            for i, t in enumerate(report.coverage_thresholds):
                if t == threshold and i < len(report.coverage_at):
                    return report.coverage_at[i]
            return float("nan")

        thresholds = sorted(
            set().union(*(r.coverage_thresholds for r in self._reports)))
        for threshold in thresholds:
            rows.append((
                f"coverage_at[{threshold}]",
                *(_coverage_lookup(r, threshold) for r in self._reports),
            ))
        return rows

    def __repr__(self):
        def _fmt(value):
            if isinstance(value, float):
                return f"{value:.4g}"
            return str(value)

        labels = self._column_labels()
        rows = self.to_table()
        metric_width = max([len("metric")] + [len(r[0]) for r in rows])
        col_width = max([12] + [len(label) for label in labels])
        header = f"{'metric':<{metric_width}}"
        for label in labels:
            header += f"  {label:>{col_width}}"
        lines = [header]
        for row in rows:
            line = f"{row[0]:<{metric_width}}"
            for value in row[1:]:
                line += f"  {_fmt(value):>{col_width}}"
            lines.append(line)
        return "\n".join(lines)


def _cluster_threshold(preset):
    preset_map = {
        "default": _oecluster.ClusterThreshold_Default,
        "tight": _oecluster.ClusterThreshold_Tight,
        "diversity": _oecluster.ClusterThreshold_Diversity,
    }
    key = str(preset).lower()
    if key not in preset_map:
        raise ValueError(f"Unknown cluster report preset: {preset!r}")
    return preset_map[key]


def cluster_report(result, distance_matrix, *, preset="default",
                   coverage_thresholds=None, boundary_threshold=None,
                   representative_method="medoid",
                   treat_noise_as_singletons=True, num_threads=0,
                   allow_nonmetric=False):
    """
    Compute a method-agnostic clustering-quality report.

    :param result: A clustering result (e.g. from :func:`butina`/:func:`dbscan`).
    :param distance_matrix: Complete SymmetricDistanceMatrix for the same items.
    :param preset: Threshold preset: "default", "tight", or "diversity".
    :param coverage_thresholds: Optional override for coverage distances.
    :param boundary_threshold: Optional override for the boundary-violation distance.
    :param representative_method: Centrality method for medoid selection:
        "medoid" (default), "minimax", or "weighted_medoid". The
        "highest_neighborhood" method is not supported.
    :param treat_noise_as_singletons: Fold noise into singleton accounting.
    :param num_threads: Reserved for parallel-safe computation.
    :param allow_nonmetric: Score anyway when the distances are known not to
        satisfy the triangle inequality. Does not override the refusals for
        similarity-valued or non-finite matrices.
    :returns: A ClusterReport.
    :raises TypeError: If result/distance_matrix have the wrong type, or
        allow_nonmetric is not a bool.
    :raises ValueError: If a preset/method/threshold is invalid, the result and
        the matrix cover different numbers of samples, the matrix uses sparse
        storage, or the matrix is not a metric.
    :raises RuntimeError: If the distance matrix cannot provide complete distances.
    """
    if not isinstance(result, ClusteringResult):
        raise TypeError("cluster_report() expects a ClusteringResult")
    if not isinstance(distance_matrix, SymmetricDistanceMatrix):
        raise TypeError("cluster_report() expects a SymmetricDistanceMatrix")

    # Every caller argument is resolved into a local before the gate runs, so
    # that an invalid preset, threshold, or representative method is reported
    # ahead of an advisory refusal that names a remedy which cannot rescue it.
    # The options object is still built only after the gate, as it was.
    threshold = _cluster_threshold(preset)

    coverage_vector = None
    if coverage_thresholds is not None:
        coverage_vector = _oecluster.DoubleVector()
        for value in coverage_thresholds:
            v = float(value)
            if v < 0.0:
                raise ValueError("coverage thresholds must be non-negative")
            coverage_vector.push_back(v)

    bt = None
    if boundary_threshold is not None:
        bt = float(boundary_threshold)
        if bt < 0.0:
            raise ValueError("boundary_threshold must be non-negative")

    method_key, native_method = _representative_method(representative_method)
    if method_key == "highest_neighborhood":
        raise ValueError(
            "cluster_report does not support the 'highest_neighborhood' "
            "representative method; use 'medoid', 'minimax', or 'weighted_medoid'")

    num_threads_int = int(num_threads)

    # num_threads reaches a size_t option field, where a negative value raises
    # OverflowError below the gate. Zero stays legal: it means "choose for me".
    if num_threads_int < 0:
        raise ValueError("num_threads must be non-negative")

    # Nothing else ties the result to the matrix: the native reporter reads
    # cluster members as storage indices, so a result scored against a larger
    # unrelated matrix stays in range and returns a confident, wrong scorecard.
    # This outranks the storage and metric refusals below -- which matrix it is,
    # sparse or dense, metric or not, cannot be the caller's first problem when
    # it is the wrong matrix. ValueError, not TypeError: both argument types are
    # right, their pairing is not.
    if result.num_samples != distance_matrix.num_samples:
        raise ValueError(
            f"cluster_report requires a result and a distance matrix over the "
            f"same items, but the result covers {result.num_samples} samples "
            f"and the matrix {distance_matrix.num_samples}")

    # ValueError, not TypeError: the argument's type is right, its storage is not.
    if isinstance(distance_matrix.storage, SparseStorage):
        raise ValueError(  # noqa: TRY004
            "cluster_report requires complete pairwise "
            "distances; SparseStorage is not supported")

    # Local argument validation first: allow_nonmetric cannot rescue a bad
    # preset, representative method, or sparse storage, so the gate must not
    # pre-empt those messages.
    _gate.require_metric(distance_matrix, "cluster_report",
                         allow_nonmetric=allow_nonmetric)

    options = _oecluster.ClusterReportOptions(threshold)
    if coverage_vector is not None:
        options.coverage_thresholds = coverage_vector
    if bt is not None:
        options.boundary_threshold = bt
    options.representative_method = native_method
    options.treat_noise_as_singletons = bool(treat_noise_as_singletons)
    options.num_threads = num_threads_int

    native = _cluster_report(
        _native_clustering_result(result), distance_matrix.storage, options)
    return ClusterReport(native, method=result.method)


def compare_reports(*reports):
    """
    Compare two or more ClusterReports side by side.

    :param reports: Two or more ClusterReport objects.
    :returns: A ClusterReportComparison.
    :raises TypeError: If any argument is not a ClusterReport.
    :raises ValueError: If fewer than two reports are given.
    """
    if len(reports) < 2:
        raise ValueError("compare_reports() requires at least two reports")
    for report in reports:
        if not isinstance(report, ClusterReport):
            raise TypeError("compare_reports() expects ClusterReport objects")
    return ClusterReportComparison(reports)


def descriptor_statistics(mols, *, sources=None, columns=None, groups=None,
                          inverse_covariance=False):
    """
    Compute per-column descriptor statistics over a molecule set.

    The statistics are the ones the descriptor comparison fits internally, so
    computing them here and passing ``variances=`` back to :func:`pdist` gives
    a reusable, explicitly scoped standardization.

    :param mols: List of OEMolBase molecules.
    :param sources: Descriptor source names: "openeye" (default), "mordred",
        or "rdkit".
    :param columns: Optional explicit column names to restrict to.
    :param groups: Optional descriptor group names to restrict to.
    :param inverse_covariance: Also compute the pseudo-inverse covariance over
        the surviving columns, for use with ``metric="mahalanobis"``.
    :returns: Dict with keys ``columns``, ``mean``, ``variance``, ``minimum``,
        ``maximum``, ``present_count``, ``dropped`` (a list of
        ``(name, reason)`` pairs), ``num_rows``, ``inverse_covariance``
        (a ``(k, k)`` array or None), ``inverse_covariance_rank``, and
        ``inverse_covariance_rows``. The last is the row count the covariance
        was actually fitted over, which can be smaller than ``num_rows``:
        covariance uses listwise deletion, so a molecule missing any selected
        descriptor still reaches the per-column statistics but not the
        covariance.
    :raises RuntimeError: If the descriptor layer refuses the request. Among
        the reasons: an unknown source, column, or group name; fewer than two
        molecules, which is too few to fit a variance; and a selection whose
        columns are all constant.
    """
    options = _oecluster.DescriptorStatisticsOptions()
    if sources is not None:
        options.sources = _comparisons._string_vector(sources)
    if columns is not None:
        options.columns = _comparisons._string_vector(columns)
    if groups is not None:
        options.groups = _comparisons._string_vector(groups)
    options.inverse_covariance = bool(inverse_covariance)

    native = _oecluster.descriptor_statistics(mols, options)

    names = list(native.columns)
    result = {
        'columns': names,
        'mean': list(native.mean),
        'variance': list(native.variance),
        'minimum': list(native.minimum),
        'maximum': list(native.maximum),
        'present_count': [int(count) for count in native.present_count],
        'dropped': list(zip(list(native.dropped_columns),
                            list(native.dropped_reasons))),
        'num_rows': int(native.num_rows),
        'inverse_covariance': None,
        'inverse_covariance_rank': int(native.inverse_covariance_rank),
        'inverse_covariance_rows': int(native.inverse_covariance_rows),
    }

    flat = list(native.inverse_covariance)
    if flat:
        k = len(names)
        result['inverse_covariance'] = np.array(flat).reshape(k, k)
    return result


# Python wrapper classes for comparison construction
class FingerprintComparison:
    """Fingerprint-based comparison using OEFP scalar metrics."""

    def __new__(cls, mols, *, fp_type=None, metric=None, numbits=None,
                min_distance=None, max_distance=None, similarity=False):
        """
        Construct a FingerprintComparison.

        :param mols: List of OEMolBase molecules.
        :param fp_type: Fingerprint type.
        :param metric: OEFP scalar metric name.
        :param numbits: Fingerprint size in bits.
        :param min_distance: Minimum Atom Pair graph distance.
        :param max_distance: Morgan radius or maximum Atom Pair graph distance.
        :param similarity: Return similarity instead of distance.
        :returns: C++ FingerprintComparison object.
        """
        opts = FingerprintOptions()
        opts.similarity = similarity
        if fp_type is not None:
            opts.fp_type = fp_type
        if metric is not None:
            opts.metric = metric
        if numbits is not None:
            opts.numbits = numbits
        if min_distance is not None:
            opts.min_distance = min_distance
        if max_distance is not None:
            opts.max_distance = max_distance
        return _FingerprintComparison(mols, opts)


class ROCSComparison:
    """ROCS-style shape overlay comparison."""

    def __new__(cls, mols, *, similarity=False):
        """
        Construct a ROCSComparison.

        :param mols: List of OEMol molecules with 3D coordinates.
        :param similarity: Return similarity instead of distance.
        :returns: C++ ROCSComparison object.
        """
        opts = ROCSOptions()
        opts.similarity = similarity
        return _ROCSComparison(mols, opts)


class SuperposeComparison:
    """Protein superposition comparison using oespruce OESuperpose."""

    def __new__(cls, items, *, method="global_carbon_alpha", similarity=False,
                predicate=None, ref_predicate=None, fit_predicate=None):
        """
        Construct a SuperposeComparison.

        :param items: List of OEDesignUnit or OEMolBase objects.
        :param method: Superposition method name.
        :param similarity: Return similarity instead of distance.
        :param predicate: oeselect expression for both ref and fit.
        :param ref_predicate: Override predicate for ref.
        :param fit_predicate: Override predicate for fit.
        :returns: C++ SuperposeComparison object.
        """
        opts = SuperposeOptions()
        method_map = {
            'global_carbon_alpha': _oecluster.SuperposeMethod_GlobalCarbonAlpha,
            'global': _oecluster.SuperposeMethod_Global,
            'ddm': _oecluster.SuperposeMethod_DDM,
            'weighted': _oecluster.SuperposeMethod_Weighted,
            'sse': _oecluster.SuperposeMethod_SSE,
            'sitehopper': _oecluster.SuperposeMethod_SiteHopper,
        }
        if method not in method_map:
            raise ValueError(f"Unknown superpose method: {method}")
        opts.method = method_map[method]
        opts.similarity = similarity
        if predicate is not None:
            opts.predicate = predicate
        if ref_predicate is not None:
            opts.ref_predicate = ref_predicate
        if fit_predicate is not None:
            opts.fit_predicate = fit_predicate
        return _SuperposeComparison(items, opts)


class DescriptorComparison:
    """Distance in standardized molecular-descriptor space.

    The molecules are taken as given. Unlike ``pdist(mols, "descriptor")``,
    which resolves the complete-case mask before it fixes the item list, this
    constructor applies no filtering: a prebuilt comparison object fixes its
    own size.
    """

    def __new__(cls, mols, *, sources=None, columns=None, groups=None,
                metric=None, variances=None, inverse_covariance=None,
                missing=None, p=None):
        """
        Construct a DescriptorComparison.

        :param mols: List of OEMolBase molecules.
        :param sources: Descriptor source names; defaults to ["openeye"].
        :param columns: Optional explicit column names.
        :param groups: Optional descriptor group names.
        :param metric: Descriptor metric name; defaults to
            "standardized_euclidean".
        :param variances: Explicit per-column variances, bypassing the pooled
            fit. The length must equal the number of selected columns, and a
            ``columns`` list passed alongside must already be in ascending
            schema order, which is the order ``descriptor_statistics`` returns.
        :param inverse_covariance: Explicit inverse covariance for
            "mahalanobis", bypassing the pooled fit.
        :param missing: Missing-value policy: "complete_case" (default),
            "propagate", or "ignore".
        :param p: Minkowski order.
        :returns: C++ DescriptorComparison object.
        :raises RuntimeError: If the C++ layer rejects an option value, or if
            the policy is "complete_case" and a molecule has an absent or
            non-finite value for a selected descriptor.
        """
        kwargs = {
            'sources': sources,
            'columns': columns,
            'groups': groups,
            'metric': metric,
            'variances': variances,
            'inverse_covariance': inverse_covariance,
            'missing': missing,
            'p': p,
        }
        return _DescriptorComparison(mols,
                                     _comparisons.descriptor_options(kwargs))
