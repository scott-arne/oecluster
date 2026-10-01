"""
oecluster: High-performance pairwise distance computation for molecular datasets.

This package provides efficient computation of pairwise distance matrices for
molecular and protein structure datasets using OpenEye toolkits. It supports
multiple comparison methods including fingerprint similarity, ROCS shape
overlay, protein superposition, and binding site comparison.
"""

import abc
import collections.abc
import ctypes
import hashlib
import importlib.machinery
import importlib.util
import json
import math
import operator
import os
import re
import shutil
import sys
import warnings
from importlib import metadata
from pathlib import Path
from typing import Any, ClassVar, NamedTuple

import numpy as np

__version__ = "5.11.1"
__version_info__ = (5, 11, 1)


_OPENEYE_COMPAT_PRELOAD_PATHS: list[str] = []
_OPENEYE_COMPAT_EXTENSION_DIR: Path | None = None
# Grouped by role -- version, storage, matrices, results, then the callables --
# so the export list reads as a tour of the API. Alphabetical order would
# interleave those groups, so RUF022 is suppressed rather than applied.
__all__ = [  # noqa: RUF022
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
    "ClusterRecord",
    "ClusterReport",
    "ClusterReportComparison",
    "ClusterReportRequested",
    "PartitionAgreement",
    "PartitionAgreementRequested",
    "SARCoherence",
    "ClusterActivity",
    "ActivityLandscape",
    "Modelability",
    "ClassConcordance",
    "MaxMinSelection",
    "CirclesResult",
    "VendiResult",
    "LogDetResult",
    "cluster_report",
    "compare_reports",
    "partition_agreement",
    "scaffold_agreement",
    "sar_coherence",
    "activity_landscape",
    "modelability",
    "maxmin_select",
    "circles",
    "vendi_score",
    "logdet_diversity",
    "sphere_exclusion",
    "knn_graph",
    "jarvis_patrick",
    "leiden",
    "ButinaResult",
    "DBSCANResult",
    "HDBSCANResult",
    "AgglomerativeResult",
    "BitBirchResult",
    "KMedoidsResult",
    "SphereExclusionResult",
    "KNNGraph",
    "JarvisPatrickResult",
    "LeidenResult",
    "MurckoResult",
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
    "k_medoids",
    "murcko",
    "murcko_scaffolds",
    "ButinaOptions",
    "RepresentativeOptions",
    "RepresentativeWeights",
    "DBSCANOptions",
    "HDBSCANOptions",
    "AgglomerativeOptions",
    "BitBirchOptions",
    "BitBirchReclusteringOptions",
    "BitBirchRefinementOptions",
    "KMedoidsOptions",
    "FingerprintComparison",
    "ROCSComparison",
    "SuperposeComparison",
    "DescriptorComparison",
    "MCSComparison",
    "RMSDComparison",
    "descriptor_statistics",
]

# Maximum value for size_t fields forwarded to C++
_SIZE_T_MAX = (1 << (8 * ctypes.sizeof(ctypes.c_size_t))) - 1


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
        if ".so" in lib_name or lib_name.endswith((".dylib", ".dll"))
    ]


def _is_openeye_runtime_library_name(lib_name):
    """Return whether a dependency belongs to the OpenEye runtime set."""
    return lib_name.startswith(("liboe", "libzstd."))


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
        if file_name.startswith((f"{family}-", f"{family}.")):
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

# Import C++ bindings from SWIG module. Ordered storage backends, then options,
# then functions, mirroring the C++ headers; isort would both alphabetize that
# away and split the aliased names into a dozen separate statements.
try:
    from .oecluster import (  # noqa: I001
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
        KMedoidsOptions,
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
        k_medoids_cluster as _k_medoids_cluster,
        murcko_scaffolds as _murcko_scaffolds,
        murcko_cluster as _murcko_cluster,
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
:ivar cutoff: Not read by the native pdist; a SparseStorage applies its own
    cutoff.
:ivar progress: Optional callback(completed, total) for progress reporting.
"""

# Import order in this module is load-bearing -- everything here runs only after
# the preload sequence above -- so the block is ordered deliberately, with the
# sibling modules that pull in the extension themselves left until last. isort
# would hoist them to the front and detach the comment below from its import.
from .oecluster import DescriptorComparison as _DescriptorComparison  # noqa: I001
# Re-exported only. The wrapper class below builds its options through
# _comparisons, but oecluster.FingerprintOptions is a documented name and has
# to keep resolving here; the redundant alias says so rather than leaving the
# import looking dead.
from .oecluster import FingerprintOptions as FingerprintOptions
from .oecluster import MCSComparison as _MCSComparison
from .oecluster import RMSDComparison as _RMSDComparison
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


def _flag(value, argument_name):
    """Coerce a caller-supplied boolean option, refusing strings.

    Every non-empty str is truthy, so ``bool("false")`` and ``bool("no")`` are
    both ``True`` and a flag the caller meant to clear would silently be set.
    Only str is rejected: ints, numpy bools and None all have an unambiguous
    truth value, and callers do pass ``1`` and ``0``.

    :func:`cluster_report` deliberately does not route through this. Its flags
    gate expensive stages, so it admits bool and numpy.bool_ and nothing else,
    refusing the ``1`` this function allows.

    :param value: The option as the caller passed it.
    :param argument_name: Parameter name, for the error message.
    :returns: The value as a native bool.
    :raises TypeError: If ``value`` is a str.
    """
    if isinstance(value, str):
        raise TypeError(
            f"{argument_name} must be a bool, not a str: every non-empty "
            f"string is true, so {value!r} would enable it")
    return bool(value)


def _fill_dense_storage(storage, condensed):
    """
    Copy a condensed distance array into a dense storage buffer in one pass.

    Writing through ``storage.Set`` costs one SWIG call per pair, which is
    quadratic in the item count; the storage buffer is contiguous float64, so
    a single copy replaces the loop.

    :param storage: DenseStorage whose ``NumPairs()`` matches ``condensed``.
    :param condensed: 1-D array of condensed distances.
    :raises ValueError: If the input is complex, or if the shape does not
        match the storage capacity.
    """
    num_pairs = storage.NumPairs()
    # Ahead of the conversion, which is what discards the imaginary part.
    # ``from_condensed`` refuses complex before it ever reaches here, so this
    # guard exists for ``from_file``: the .npz payload is caller-supplied and
    # crosses into the buffer with no dtype check of its own.
    if np.iscomplexobj(condensed):
        raise ValueError(
            f"expected real distances, got a complex input (dtype "
            f"{np.asarray(condensed).dtype}); a distance is a real number")
    values = np.ascontiguousarray(condensed, dtype=np.float64)
    # ``ndim`` and ``size``, not ``shape[0]``: an array whose leading
    # dimension happens to match reaches ``np.copyto`` and is refused there
    # with a broadcast message that names neither this helper nor the storage.
    if values.ndim != 1 or values.size != num_pairs:
        raise ValueError(
            f"condensed shape {values.shape} does not match the storage "
            f"capacity {num_pairs}")
    if num_pairs == 0:
        return

    buffer = np.ctypeslib.as_array(
        ctypes.cast(storage._data_ptr(), ctypes.POINTER(ctypes.c_double)),
        shape=(num_pairs,))
    np.copyto(buffer, values)


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
        """Get the number of distinct violating triples the probe found."""
        return self._facts['probe_violations']

    @property
    def probe_sampled(self):
        """Get the number of distinct triples the probe tested."""
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
            # Checked before anything is written, not only on load. This is a
            # backstop now: ``SparseStorage.Set`` refuses an out-of-range or
            # diagonal pair outright, where it once enforced ``i != j`` with a
            # bare ``assert`` compiled out of release builds -- so a poked-at
            # matrix reached here and wrote a file ``from_file`` then refused,
            # and the caller lost the data and only found out on the next load.
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
            _fill_dense_storage(storage, condensed)

        return cls(storage, comparison_name, labels, params, facts)

    @classmethod
    def from_condensed(cls, values, *, labels=None,
                       comparison_name="precomputed", params=None, check=True,
                       probe_triples=100000, seed=0):
        """
        Build a distance matrix from distances computed elsewhere.

        Accepts either a 1-D condensed array (upper triangle, row-major, as
        ``scipy.spatial.distance.pdist`` returns) or a 2-D square matrix,
        whose strict upper triangle is taken.

        Ingress measures what the numbers alone can show: values must be
        finite and non-negative, a square input must additionally be symmetric
        with a diagonal near zero, and a sample of triples is tested against
        the triangle inequality. A violation found here is recorded and later
        refused by the clustering entry points unless ``allow_nonmetric=True``.
        A computed matrix instead reads a ``triangle`` capability off its
        comparison object, which can be a positive claim; there is no object
        to ask here, which is why the probe exists.

        All three capability facts -- ``is_distance``, ``zero_self`` and
        ``triangle`` -- are stamped ``"unknown"``. There is no metric object
        to interrogate, and the caller asserting a value is not evidence, so
        the probe result is the only claim this path makes.

        :param values: 1-D condensed array or 2-D square distance matrix. A
            square input must be symmetric to within ``rtol=1e-5,
            atol=1e-8`` -- loose enough for a float32-computed matrix, which
            is asymmetric in its last bit -- and its strict upper triangle is
            the half kept, so a disagreement below that tolerance is
            discarded unread. Its diagonal must likewise be zero only to
            within ``3e-2`` of the largest distance, which the same float32
            producers miss by roundoff; the diagonal is discarded too, so
            this catches a matrix that is not a distance matrix rather than
            an imprecise one. Neither array is accepted masked: the values
            behind an active mask would be read as distances. The
            two empty shapes differ: a length-0 condensed array is one item,
            which is what the single-item ``pdist`` round trip needs, while a
            0x0 square is no items.
        :param labels: Optional item labels; must match the item count.
        :param comparison_name: Name recorded on the matrix.
        :param params: Optional provenance recorded on the matrix and echoed
            by ``repr``; copied, so later mutation of the caller's dict does
            not change the matrix.
        :param check: When False, skip the value checks and the probe
            together, and leave ``data_integrity`` at ``"unknown"`` since
            nothing measured it. The refusals that do not read the numbers
            still run -- among them the shape and length arithmetic, the
            label count, the masked-array refusal and the complex-dtype
            refusal: an input that maps to no item count, or whose numbers
            cannot be stored as given, cannot be stored at all.
        :param probe_triples: Triples to sample for the triangle probe; a
            non-positive count skips it.
        :param seed: Seed for the probe sampler.
        :returns: SymmetricDistanceMatrix.
        :raises TypeError: If check is not a bool or ``numpy.bool_``.
        :raises ValueError: If the shape, dtype, length, or values are
            invalid.
        """
        # Rejected before the input is even read: a malformed switch is a call
        # the caller must fix whatever the numbers look like. Truthiness would
        # be the wrong rule, since check=None reads to a caller as
        # "unspecified" while turning off every check that reads the numbers.
        if not isinstance(check, (bool, np.bool_)):
            raise TypeError(
                f"check must be True or False, not {type(check).__name__} "
                f"({check!r}). A falsy value would silently disable the value "
                f"checks and the probe.")

        # Ahead of ``np.asarray``, because that call is what loses the mask:
        # it hands back the data buffer, so the entries the caller marked
        # invalid are read as distances and stamped ``complete``. Refused
        # rather than filled, since there is no distance to substitute. A
        # masked array with nothing actually masked is accepted -- the caller
        # loses nothing, and refusing it would be refusing a container.
        if np.ma.isMaskedArray(values):
            masked = int(np.ma.getmaskarray(values).sum())
            if masked:
                raise ValueError(
                    f"expected a plain array, got a masked array with "
                    f"{masked} masked entries; the values behind the mask "
                    f"would be read as distances")

        # Ahead of the float64 conversion, which keeps the real part behind a
        # ComplexWarning the caller may have filtered. Not gated on check=:
        # no assertion by the caller makes half of a complex number a
        # distance. The message states the contract rather than a loss, since
        # an all-zero imaginary part loses nothing and is still refused.
        array = np.asarray(values)
        if np.iscomplexobj(array):
            raise ValueError(
                f"expected real distances, got a complex input (dtype "
                f"{array.dtype}); a distance is a real number")
        array = array.astype(np.float64, copy=False)

        square = None
        if array.ndim == 2:
            if array.shape[0] != array.shape[1]:
                raise ValueError(
                    f"a 2-D input must be square, got shape "
                    f"{tuple(array.shape)}")
            n = array.shape[0]
            square = array
            condensed = np.ascontiguousarray(array[np.triu_indices(n, k=1)])
        elif array.ndim == 1:
            length = array.shape[0]
            n = round((1.0 + np.sqrt(1.0 + 8.0 * length)) / 2.0)
            if n * (n - 1) // 2 != length:
                raise ValueError(
                    f"{length} is not a valid condensed length: no item count "
                    f"n satisfies n * (n - 1) / 2 == {length}")
            condensed = np.ascontiguousarray(array)
        else:
            raise ValueError(
                f"expected a 1-D or 2-D array, got {array.ndim} dimensions")

        # Ahead of the value checks, because no number in the input can make a
        # miscounted labels= valid. Leaving it below them made the answer to
        # one bad argument depend on the data and on check=.
        if labels is not None:
            labels = list(labels)
            if len(labels) != n:
                raise ValueError(
                    f"labels length {len(labels)} != item count {n}")

        if check:
            # The whole input, so a 2-D diagonal is measured too: a NaN there
            # is a missing self-distance, and calling it a non-zero one would
            # send the caller to recompute rather than to impute or drop.
            # Every check below can then assume finite numbers.
            if not np.all(np.isfinite(array)):
                raise ValueError(
                    "distance matrix contains non-finite values; remove or "
                    "impute them before clustering")
            if square is not None:
                # A stated tolerance that happens to equal np.allclose's
                # defaults, not an inherited one. A tighter rtol=1e-9 was
                # tried and refused ordinary float32 output: the Gram-trick
                # euclidean that torch.cdist, faiss and cuML use adds the two
                # squared norms in sequence, so d(i, j) and d(j, i) differ in
                # the association order alone -- measured up to 1.99e-7
                # relative, under two float32 ulp, across 150 configurations.
                # 1e-5 is roughly 84 of those ulp, which leaves room for a
                # worse-conditioned pipeline than any measured here.
                #
                # The cost is real and accepted: a float64 producer whose
                # halves genuinely disagree at 1e-6 passes, and the lower one
                # is then discarded unread. The input's dtype cannot separate
                # the two cases -- the float32 matrix above arrives as
                # float64, having been widened by a caller for unrelated
                # reasons -- so a dtype-dependent tolerance would only make
                # the same data refuse or pass by accident.
                if not np.allclose(square, square.T, rtol=1e-5, atol=1e-8):
                    raise ValueError("a 2-D distance matrix must be symmetric")
                # The diagonal is discarded with the rest of the lower
                # triangle, so this check is here to catch a matrix that is
                # not a distance matrix at all -- not to police roundoff.
                # Exact zero did the latter, and refused the same float32
                # Gram-trick output the tolerance above was widened to accept:
                # self-distances of 2.8e-3 against distances of 10.8 on the
                # set the paired test builds.
                #
                # 3e-2 relative to the largest distance is the geometric
                # midpoint of the two regimes, roughly 4.6x clear of each.
                # Measured roundoff tops out at 6.6e-3 (offset 10, d = 128;
                # it grows with the offset, because the trick cancels the
                # squared norms against themselves). The cheapest real
                # confusion measured is a caller who forgot to zero a diagonal
                # of 1.0 against distances of order 7, at 1.4e-1; a similarity
                # matrix sits at 1.2 and an RBF kernel at 5.5.
                #
                # Widening costs no safety here, and that is worth stating:
                # zeroing the diagonal already admits the same numbers at any
                # offset, so this check was never guarding precision. It only
                # ever refused callers who had not pre-zeroed.
                zero_diagonal_rtol = 3e-2
                # ``if n`` because a 0x0 square has no diagonal to reduce over
                # and ``max`` has no identity for an empty array.
                diagonal = (float(np.abs(np.diagonal(square)).max())
                            if n else 0.0)
                # ``abs`` on the scale too, because the negative-value check
                # runs below this one: a matrix of negative distances has to
                # be reported as negative, not as a diagonal that failed a
                # negative tolerance. With no distance to measure against
                # (n < 2, or every distance zero) the only defensible
                # tolerance is zero, which is what a zero scale gives.
                scale = float(np.abs(condensed).max()) if condensed.size else 0.0
                if diagonal > zero_diagonal_rtol * scale:
                    raise ValueError(
                        f"a 2-D distance matrix must have a zero diagonal, "
                        f"within {zero_diagonal_rtol:g} of the largest "
                        f"distance ({scale:g}); the largest self-distance is "
                        f"{diagonal:g}")
            if np.any(condensed < 0.0):
                raise ValueError("distance matrix contains negative values")

        # is_distance, zero_self and triangle all stay "unknown" from
        # default_facts(). A zero diagonal, where there was one to check, is
        # the caller's diagonal rather than a property of whatever produced
        # the numbers, so it proves nothing durable -- and calling the
        # argument a "distance matrix" is not evidence that it is one.
        facts = _gate.default_facts()
        if check:
            # Only now, having measured it. Under check=False nothing looked
            # at the numbers, and ``to_file`` would write the claim into the
            # .npz for ``from_file`` to read back as evidence.
            facts['data_integrity'] = "complete"
            facts.update(_gate.probe_triangle(condensed, n,
                                              samples=probe_triples,
                                              seed=seed))

        storage = DenseStorage(n)
        _fill_dense_storage(storage, condensed)
        # ``is not None`` rather than truthiness, matching
        # ``DistanceMatrix.__init__``: an explicitly passed value must not be
        # silently equivalent to passing nothing. ``dict`` because this path
        # promises the caller a copy.
        copied = dict(params) if params is not None else {}
        return cls(storage, comparison_name, labels, copied, facts)

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


class KMedoidsResult(ClusteringResult):
    """k-medoids clustering result with medoids and the objective value.

    Medoids are real members of the input, one per cluster, sorted ascending by
    item index so that label ``i`` always belongs to ``medoids[i]``. ``cost`` is
    the sum of every item's distance to its assigned medoid, recomputed from the
    returned assignment rather than accumulated from the optimizer's deltas.
    """

    def __init__(self, labels, clusters, *, medoids=(), cost=0.0,
                 n_iterations=0, converged=False, native_owner=None):
        super().__init__(labels, clusters, native_owner=native_owner)
        self._medoids = tuple(int(medoid) for medoid in medoids)
        self._cost = float(cost)
        self._n_iterations = int(n_iterations)
        self._converged = bool(converged)

    @property
    def medoids(self):
        """Medoid item index per cluster, ascending and aligned with labels."""
        return self._medoids

    @property
    def cost(self):
        """Sum of each item's distance to its assigned medoid."""
        return self._cost

    @property
    def n_iterations(self):
        """Swap iterations performed."""
        return self._n_iterations

    @property
    def converged(self):
        """True when no single medoid swap lowers the reported cost.

        False means ``max_iterations`` was reached; the partition is valid but
        no local-optimality claim is made about it.
        """
        return self._converged

    @property
    def method(self):
        return "k_medoids"


class SphereExclusionResult(ClusteringResult):
    """Sphere-exclusion clustering result with one center per cluster.

    Clusters are in center order and list their center first, then the other
    members in ascending position. ``centers[i]`` is ``clusters[i][0]``.
    Positions refer to the caller's items; an item that normalization dropped
    has the label -1, appears in no cluster, and is listed in ``excluded``.
    """

    def __init__(self, labels, clusters, *, centers=(), excluded=(),
                 native_owner=None):
        super().__init__(labels, clusters, native_owner=native_owner)
        self._centers = tuple(int(center) for center in centers)
        self._excluded = [list(entry) for entry in excluded]

    @property
    def centers(self):
        """Center position per cluster, in cluster order."""
        return self._centers

    @property
    def excluded(self):
        """``[position, reason]`` for each item normalization dropped."""
        return [list(entry) for entry in self._excluded]

    @property
    def method(self):
        return "sphere_exclusion"


class KNNGraph:
    """Each item's ``k`` nearest other items, with their raw distances.

    Built only by :func:`knn_graph`. Rows are the items kept after
    normalization, in native order; ``positions[r]`` is row ``r``'s caller
    position. ``indices`` holds caller positions, and each row is ordered by
    ascending (distance, position), so equal distances go to the lower
    position. The values are distances, not affinities. Excluded items have
    no row and appear in no neighbor list.
    """

    def __init__(self):
        raise TypeError("KNNGraph is built by knn_graph()")

    @classmethod
    def _from_native(cls, native, source):
        """
        Wrap a native graph built over a dispatched input.

        :param native: The native KNNGraph.
        :param source: The :class:`_DiversitySource` it was built from.
        :returns: A :class:`KNNGraph`.
        """
        graph = cls.__new__(cls)
        graph._native = native
        graph._k = int(native.K())
        rows = int(native.NumItems())
        count = rows * graph._k
        positions = (np.arange(rows, dtype=np.int64) if source.positions is None
                     else np.asarray(source.positions, dtype=np.int64))
        native_indices = np.fromiter(native.Indices(), dtype=np.int64,
                                     count=count).reshape(rows, graph._k)
        graph._indices = positions[native_indices]
        graph._distances = np.fromiter(native.Distances(), dtype=np.float64,
                                       count=count).reshape(rows, graph._k)
        graph._row_positions = positions
        graph._positions = source.positions
        graph._num_positions = source.num_positions
        graph._excluded = [list(entry) for entry in source.excluded]
        return graph

    @property
    def k(self):
        """Neighbors per row."""
        return self._k

    @property
    def indices(self):
        """``(rows, k)`` ``int64`` array of neighbor caller positions."""
        return self._indices.copy()

    @property
    def distances(self):
        """``(rows, k)`` ``float64`` array of neighbor distances."""
        return self._distances.copy()

    @property
    def positions(self):
        """``(rows,)`` ``int64`` array: the caller position of each row."""
        return self._row_positions.copy()

    @property
    def excluded(self):
        """``[position, reason]`` for each item normalization dropped."""
        return [list(entry) for entry in self._excluded]

    def __len__(self):
        """Number of rows."""
        return len(self._row_positions)

    def __repr__(self):
        return f"KNNGraph(num_rows={len(self)}, k={self._k})"


class JarvisPatrickResult(ClusteringResult):
    """Jarvis-Patrick clustering result with the ``k`` and ``kmin`` used.

    Clusters are ordered by their smallest member and list members in
    ascending position; an item that links to nothing is a singleton.
    Positions refer to the caller's items; an item that normalization dropped
    has the label -1, appears in no cluster, and is listed in ``excluded``.
    """

    def __init__(self, labels, clusters, *, k=0, kmin=0, excluded=(),
                 native_owner=None):
        super().__init__(labels, clusters, native_owner=native_owner)
        self._k = int(k)
        self._kmin = int(kmin)
        self._excluded = [list(entry) for entry in excluded]

    @property
    def k(self):
        """Neighbors per item in the graph."""
        return self._k

    @property
    def kmin(self):
        """Shared neighbors required to link a mutual pair."""
        return self._kmin

    @property
    def excluded(self):
        """``[position, reason]`` for each item normalization dropped."""
        return [list(entry) for entry in self._excluded]

    @property
    def method(self):
        return "jarvis_patrick"


class LeidenResult(ClusteringResult):
    """Leiden community detection result with the optimization settings.

    Clusters are ordered by their smallest member and list members in
    ascending position. Positions refer to the caller's items; an item that
    normalization dropped has the label -1, appears in no cluster, and is
    listed in ``excluded``.
    """

    def __init__(self, labels, clusters, *, quality=0.0, iterations=0,
                 objective="modularity", resolution=0.0, k=0, excluded=(),
                 native_owner=None):
        super().__init__(labels, clusters, native_owner=native_owner)
        self._quality = float(quality)
        self._iterations = int(iterations)
        self._objective = str(objective)
        self._resolution = float(resolution)
        self._k = int(k)
        self._excluded = [list(entry) for entry in excluded]

    @property
    def quality(self):
        """Objective value of the partition on the SNN graph."""
        return self._quality

    @property
    def iterations(self):
        """Full passes run, including a final pass that changed nothing."""
        return self._iterations

    @property
    def objective(self):
        """``"modularity"`` or ``"cpm"``."""
        return self._objective

    @property
    def resolution(self):
        """Resolution the objective used."""
        return self._resolution

    @property
    def k(self):
        """Neighbors per item in the graph."""
        return self._k

    @property
    def excluded(self):
        """``[position, reason]`` for each item normalization dropped."""
        return [list(entry) for entry in self._excluded]

    @property
    def method(self):
        return "leiden"


class MurckoResult(ClusteringResult):
    """Murcko scaffold clustering result with the scaffold strings.

    Clusters are scaffold identity classes: two molecules share a cluster
    exactly when their scaffolds canonicalize to the same non-empty SMILES.
    Labels come from sorting the distinct non-empty scaffold strings, so
    ``cluster_scaffolds`` is ascending and ``cluster_scaffolds[i]`` names
    ``clusters[i]``. A molecule with no ring system has an empty scaffold
    string, which names no cluster: it carries ``-1`` however many other
    acyclic molecules the input holds.
    """

    def __init__(self, labels, clusters, *, scaffolds=(), cluster_scaffolds=(),
                 native_owner=None):
        super().__init__(labels, clusters, native_owner=native_owner)
        self._scaffolds = tuple(str(scaffold) for scaffold in scaffolds)
        self._cluster_scaffolds = tuple(
            str(scaffold) for scaffold in cluster_scaffolds)

    @property
    def scaffolds(self):
        """Per-item scaffold SMILES in input order; '' for an acyclic molecule."""
        return self._scaffolds

    @property
    def cluster_scaffolds(self):
        """Scaffold naming each cluster, ascending; entry i belongs to label i."""
        return self._cluster_scaffolds

    @property
    def method(self):
        return "murcko"


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

    Warning: the first element of ``items`` decides how the whole list is
    read. For ``rocs`` an ``OEGraphMol`` there silently reduces every later
    ``OEMol`` to its active conformer and the scores move with it; the reverse
    ordering raises ``TypeError`` instead. Convert up front with
    ``oechem.OEMol(...)`` when a list would otherwise be mixed. ``mcs`` follows
    the same rule for which lists are accepted, but its scores are unchanged.
    See the ROCS and MCS sections of ``docs/python-api.md``.

    :param items: List of molecules, design units, or other items.
    :param comparison: Comparison method: "fingerprint", "rocs", "superpose",
                       "sitehopper", "descriptor", "rmsd", "mcs", or a C++
                       comparison object.
    :param similarity: Return similarities instead of distances.
    :param num_threads: Number of threads (0 = auto).
    :param chunk_size: Pairs per work unit.
    :param cutoff: Distance cutoff for sparse storage (0 = store all).
    :param output: Optional file path for memory-mapped storage.
    :param progress: Optional callback(completed, total).
    :param kwargs: Comparison-specific options.
    :returns: SymmetricDistanceMatrix with computed distances/similarities.
    :raises TypeError: If unknown kwargs are passed.
    :raises ValueError: If cutoff > 0 with no ``output`` and either
        similarity=True on a named comparison or a prebuilt comparison that
        reports similarities, or if normalizing the inputs leaves no items.
    """
    if isinstance(comparison, str):
        # Before normalization: filtering can empty the list, and the refusal
        # below would then answer for an argument no input could rescue.
        _comparisons.validate_request(comparison, similarity, kwargs)

        # Decided once and reused at the storage branch below: testing
        # ``cutoff > 0.0`` in both places would let a value whose comparison
        # is not stable answer differently there, selecting sparse storage
        # for a call this guard had already cleared. ``bool`` is what makes
        # the decision final. ``and`` yields its operand, so without it the
        # binding holds whatever ``__gt__`` returned, and the guard and the
        # storage branch each convert that to a truth value again -- the same
        # divergence, one level down.
        #
        # Below validate_request because the cutoff message names dropping
        # the cutoff as the remedy, which cannot fix a misspelled argument.
        # That only reorders the two for a comparison that validates its
        # keywords in validate_request; most register no validator, and for
        # those this refusal comes first regardless. cdist orders them the
        # same way in both cases.
        #
        # An mmap output never consults the cutoff, so this refusal skips it.
        # The prebuilt branch ignores the ``similarity`` argument, so this
        # refusal cannot guard it; that branch asks the object for its
        # orientation instead.
        sparse = output is None and bool(cutoff > 0.0)
        if sparse and similarity:
            raise ValueError(
                "cutoff > 0 is not supported with similarity=True: the cutoff "
                "zeroes values above the threshold, which would discard high "
                "similarities")

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
        # Deliberately a second copy rather than a hoist above the branch:
        # the string branch evaluates the cutoff only after validate_request,
        # and a hoist would reverse that order.
        sparse = output is None and bool(cutoff > 0.0)
        # Read ahead of the storage allocation and the computation, unlike
        # the facts read below: orientation is known before any pair is
        # scored, while data integrity is not. Only a reported ``False`` is
        # refused; an object without facts reports "unknown" and keeps
        # working.
        if (sparse and _gate.facts_from_comparison(
                comparison_obj)['is_distance'] is False):
            raise ValueError(
                "cutoff > 0 is not supported for a prebuilt comparison that "
                "reports similarities: the cutoff zeroes values above the "
                "threshold, which would discard high similarities. Drop the "
                "cutoff or build the comparison with similarity=False")

    n = comparison_obj.Size()

    storage: Any
    if output is not None:
        storage = MMapStorage(output, n)
    elif sparse:
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

    Warning: one comparison is built over the concatenation ``items_a +
    items_b``, so the first element of ``items_a`` decides how *both* sets are
    read, set A's own later elements included. For ``rocs`` an ``OEGraphMol``
    there silently reduces every later ``OEMol`` to its active conformer and
    the scores move with it. ``mcs`` is reduced the same way but scores
    identically, never reading coordinates, so for it only which calls are
    *accepted* differs. The reverse ordering raises ``TypeError`` instead. See
    the ROCS and MCS sections of ``docs/python-api.md`` for the measured
    effect.

    :param items_a: Reference items (rows of the result).
    :param items_b: Fit items (columns of the result).
    :param comparison: Comparison method name: "fingerprint", "rocs", "superpose",
                       "sitehopper", "descriptor", "rmsd", or "mcs". Prebuilt
                       comparison objects are not supported.
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

    # Above the two guards below, not merely above normalization as in
    # ``pdist``: both of them can describe the wrong problem. The cutoff
    # message names dropping the cutoff as the remedy, which cannot work when
    # ``similarity=True`` is itself the invalid argument, and the arrived-empty
    # message answers a shape question the caller did not ask. Hoisting also
    # makes ``pdist`` and ``cdist`` give one answer to one bad argument.
    _comparisons.validate_request(comparison, similarity, kwargs)

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
    options.reordering = _flag(reordering, "reordering")
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

    # Field order mirrors the C++ RepresentativeMetrics struct, the :ivar: list
    # above, and the assignments below; alphabetizing would desynchronize all
    # three from the native layout they document.
    __slots__ = (  # noqa: RUF023
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

    # Ordered as the native result presents them, matching the :ivar: list above.
    __slots__ = ("member", "score", "metrics")  # noqa: RUF023

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
    options.allow_single_cluster = _flag(
        allow_single_cluster, "allow_single_cluster")
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
    options.compute_full_tree = _flag(compute_full_tree, "compute_full_tree")
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


def k_medoids(distance_matrix, *, n_clusters=2, init="build",
              initial_medoids=None, max_iterations=100,
              num_threads=0, chunk_size=4096):
    """
    Cluster a precomputed distance matrix with k-medoids (PAM).

    Places exactly ``n_clusters`` medoids, each a real member of the input,
    minimizing the sum of every item's distance to its assigned medoid. Every
    item receives a label; ``-1`` is never emitted.

    Unlike :func:`butina`, :func:`dbscan`, :func:`hdbscan` and
    :func:`agglomerative`, this function does **not** require a metric and takes
    no ``allow_nonmetric`` parameter. PAM's objective is a sum of distances and
    its swap step compares two such sums, so no step appeals to the triangle
    inequality and there is nothing for a flag to override. A Dice or Tanimoto
    matrix clusters here with no override at all.

    That creates one asymmetry worth knowing about in advance:
    :func:`cluster_report` *does* assume a metric, so a matrix this function
    accepted may be refused by the report unless you pass
    ``allow_nonmetric=True`` there.

    :param distance_matrix: SymmetricDistanceMatrix returned by :func:`pdist`.
    :param n_clusters: Number of medoids to place; must be in [1, item count].
    :param init: Initialization strategy: "build" (greedy PAM BUILD),
        "farthest_first" (deterministic MaxMin), or "explicit"
        (use ``initial_medoids``). Case-insensitive.
    :param initial_medoids: Starting medoid indices; a non-empty sequence is
        required when ``init="explicit"`` and refused otherwise.
    :param max_iterations: Swap iterations before giving up. Reaching the cap
        is not an error: the result is valid and ``converged`` is False.
    :param num_threads: Worker threads; 0 auto-detects hardware concurrency.
    :param chunk_size: Items per work unit in the parallelized swap scan.
    :returns: KMedoidsResult with labels, clusters, medoids, cost, iteration
        count, and convergence.
    :raises TypeError: If distance_matrix is not a SymmetricDistanceMatrix, or
        an integer argument is not an integer.
    :raises ValueError: If options are invalid (including negative or oversized
        integers for max_iterations/num_threads/chunk_size), the matrix uses
        sparse storage, or the matrix is not comparable.
    :raises IndexError: If an explicit medoid index is outside the matrix.

    Example::

        dm = oecluster.pdist(mols, "fingerprint", metric="dice")
        result = oecluster.k_medoids(dm, n_clusters=10)
        for label, medoid in enumerate(result.medoids):
            print(label, mols[medoid].GetTitle())
    """
    if not isinstance(distance_matrix, SymmetricDistanceMatrix):
        raise TypeError("k_medoids() expects a SymmetricDistanceMatrix")

    # An unrecognized init is deliberately NOT rejected here. Its row is the
    # last one in the native table, so rejecting it first would make
    # k_medoids(dm, init="bogus", chunk_size=0) report the init while the
    # native call reports the chunk size -- the divergence the complete mirror
    # exists to prevent. An unrecognized name is carried as "not Explicit",
    # which is exactly what an unrecognized enumerator is to the native switch,
    # and reported below in its own place.
    init_map = {
        "build": _oecluster.KMedoidsInit_Build,
        "farthest_first": _oecluster.KMedoidsInit_FarthestFirst,
        "explicit": _oecluster.KMedoidsInit_Explicit,
    }
    init_key = str(init).lower()

    # Coerce caller arguments before the gate so that an invalid type is reported
    # ahead of an advisory refusal that names a remedy which cannot rescue it.
    # operator.index() accepts int, bool and numpy integers and rejects
    # float/str/None; int() would silently truncate 2.5 to two clusters.
    n_clusters_int = operator.index(n_clusters)
    max_iterations_int = operator.index(max_iterations)
    num_threads_int = operator.index(num_threads)
    chunk_size_int = operator.index(chunk_size)

    if initial_medoids is None:
        seeds = []
    else:
        try:
            seeds = [operator.index(index) for index in initial_medoids]
        except TypeError as error:
            raise TypeError(
                "k_medoids() initial_medoids must be a sequence of ints"
            ) from error

    # All four reach size_t option fields, where a negative value or an oversized
    # positive raises OverflowError below the gate. Zero stays legal for
    # num_threads: it means "choose for me". The upper-bound checks are
    # representation constraints about what a size_t can hold, so they live here
    # with the negativity checks rather than in the mirror, and they fire before
    # the gate so an oversized value cannot hide behind an advisory refusal.
    if n_clusters_int < 0:
        raise ValueError("K-medoids n_clusters must be non-negative")
    if max_iterations_int < 0:
        raise ValueError("K-medoids max_iterations must be non-negative")
    if max_iterations_int > _SIZE_T_MAX:
        raise ValueError("K-medoids max_iterations exceeds size_t maximum")
    if num_threads_int < 0:
        raise ValueError("num_threads must be non-negative")
    if num_threads_int > _SIZE_T_MAX:
        raise ValueError("num_threads exceeds size_t maximum")
    if chunk_size_int < 0:
        raise ValueError("chunk_size must be non-negative")
    if chunk_size_int > _SIZE_T_MAX:
        raise ValueError("chunk_size exceeds size_t maximum")
    if any(index < 0 for index in seeds):
        raise ValueError("K-medoids initial_medoids must be non-negative")

    # Every native check below is mirrored here, in the order
    # k_medoids_cluster() applies them, and before the gate, except the
    # null-data guard (no Python-reachable way to construct storage with null
    # data today, so that would surface as RuntimeError from the GIL wrapper
    # rather than ValueError). That is stronger than agglomerative()'s partial
    # mirror and deliberately so: the GIL wrapper collapses every native
    # exception to RuntimeError, so a check left to the native layer loses both
    # its type and its place in the order, and a caller could not write one
    # except clause for "you passed me something invalid".

    # ValueError, not TypeError: the argument's type is right, its storage is not.
    if isinstance(distance_matrix.storage, SparseStorage):
        raise ValueError(  # noqa: TRY004
            "K-medoids clustering requires complete pairwise "
            "distances; SparseStorage is not supported")
    if chunk_size_int == 0:
        raise ValueError("K-medoids chunk_size must be at least one")
    if max_iterations_int == 0:
        raise ValueError("K-medoids max_iterations must be at least one")
    if n_clusters_int == 0:
        raise ValueError("K-medoids n_clusters must be at least one")
    if n_clusters_int > distance_matrix.num_samples:
        raise ValueError(
            "K-medoids n_clusters must be at most "
            f"the item count ({distance_matrix.num_samples})")
    if init_key != "explicit":
        if seeds:
            raise ValueError(
                "K-medoids initial_medoids requires an explicit initialization")
    else:
        if len(seeds) != n_clusters_int:
            raise ValueError(
                "K-medoids initial_medoids must hold exactly n_clusters indices")
        for index in seeds:
            if index >= distance_matrix.num_samples:
                raise IndexError(
                    "K-medoids initial_medoids index is outside "
                    "the storage range")
        if len(set(seeds)) != len(seeds):
            raise ValueError("K-medoids initial_medoids must be unique")
    # Last, mirroring the last row of the native table.
    if init_key not in init_map:
        raise ValueError(
            f"Unknown k-medoids init: {init!r}; expected 'build', "
            "'farthest_first', or 'explicit'")

    # require_comparable, not require_metric: PAM never appeals to the triangle
    # inequality, so refusing a matrix for violating it would be over-refusal.
    # The subset_scored refusal it keeps matters more here than for the ranking
    # metrics that motivated it, because PAM adds the incomparable numbers up.
    _gate.require_comparable(distance_matrix, "k_medoids")

    options = KMedoidsOptions()
    options.n_clusters = n_clusters_int
    options.init = init_map[init_key]
    if seeds:
        native_seeds = _oecluster.SizeTVector()
        for index in seeds:
            native_seeds.push_back(index)
        options.initial_medoids = native_seeds
    options.max_iterations = max_iterations_int
    options.num_threads = num_threads_int
    options.chunk_size = chunk_size_int

    result = _k_medoids_cluster(distance_matrix.storage, options)
    return KMedoidsResult(
        result.Labels(),
        result.Members(),
        medoids=result.Medoids(),
        cost=result.Cost(),
        n_iterations=result.NumIterations(),
        converged=result.Converged(),
    )


# Module level rather than function local, because both entry points share it.
_SCAFFOLD_TYPES = {
    "framework": _oecluster.ScaffoldType_Framework,
    "generic": _oecluster.ScaffoldType_Generic,
}


def _murcko_options(mols, scaffold, num_threads, caller):
    """Validate the Murcko keywords and build the native options struct.

    The SWIG layer collapses every native exception to ``RuntimeError``, so the
    native checks are mirrored here -- before the call -- to give callers the
    exception type that actually describes what happened. The scaffold level is
    judged before the input, matching the native ordering; ``num_threads`` has
    no native check to mirror -- the native layer clamps rather than refuses --
    and is only bounded here to what a ``size_t`` can hold.

    :param mols: The caller's molecule list, judged only for type and emptiness.
    :param scaffold: Scaffold level name; case-insensitive.
    :param num_threads: Worker thread request.
    :param caller: Public function name, for the messages.
    :returns: A populated native ``MurckoOptions``.
    :raises TypeError: If ``scaffold`` is not a string, ``mols`` is not a list,
        or ``num_threads`` is not index-coercible.
    :raises ValueError: If ``scaffold`` is not a known level, ``mols`` is empty,
        or ``num_threads`` is negative or larger than a ``size_t``.
    """
    if not isinstance(scaffold, str):
        raise TypeError(f"{caller}() scaffold must be a string")
    key = scaffold.lower()
    if key not in _SCAFFOLD_TYPES:
        raise ValueError(
            f"Unknown Murcko scaffold type: {scaffold!r}; "
            "expected 'framework' or 'generic'")
    # The typemap is a strict PyList_Check, so anything else is rejected there
    # anyway -- but both entry points are overloaded, and SWIG's overload
    # dispatcher discards the typemap's message in favour of "Wrong number or
    # type of arguments for overloaded function 'murcko_cluster'", naming a
    # symbol murcko() never told the caller about. Mirroring the check exactly
    # keeps a tuple, which has a length and so used to slip through, from
    # reaching that dispatcher.
    if not isinstance(mols, list):
        raise TypeError(f"{caller}() requires a list of molecules")
    if not mols:
        raise ValueError(f"{caller}() requires at least one molecule")

    num_threads_int = operator.index(num_threads)
    if num_threads_int < 0:
        raise ValueError("num_threads must be non-negative")
    if num_threads_int > _SIZE_T_MAX:
        raise ValueError("num_threads exceeds size_t maximum")

    options = _oecluster.MurckoOptions()
    options.scaffold = _SCAFFOLD_TYPES[key]
    options.num_threads = num_threads_int
    return options


def murcko_scaffolds(mols, *, scaffold="framework", num_threads=0):
    """
    Compute the Bemis-Murcko scaffold of each molecule.

    Returns one canonical SMILES per molecule, in input order, suitable for
    passing straight to :func:`scaffold_agreement` or to the
    ``scaffold_labels=`` argument of the representative functions.

    A molecule with no ring system yields ``''``, which is the "missing
    scaffold" convention those consumers already use. A molecule that *has*
    rings but from which no framework can be extracted is an error rather than
    an empty string.

    Three normalizations are applied, all of them so that scaffold identity is
    string identity. Stereochemistry is dropped, because deleting sidechains can
    orphan a stereocenter and leave a configuration that means nothing;
    enantiomers and diastereomers therefore share a scaffold. Explicit hydrogens
    are suppressed, so an SD-file molecule and the same molecule read from
    SMILES agree. Hydrogen counts and formal charges are recomputed on the atoms
    the cut touched, because removing a sidechain takes its bond order with it;
    diphenyl sulfone, diphenyl sulfoxide and diphenyl sulfide therefore share
    one framework scaffold, as their carbon analogues already did, while a
    charged ring atom whose bonds all survive keeps its charge. Nothing else is
    standardized: there is no salt stripping and no largest-component
    selection. An acyclic counter-ion contributes no framework, so most salt
    forms already match the free base; a counter-ion that *also* carries a ring
    -- a tosylate or besylate -- yields one ``.``-joined scaffold that will
    not. Strip salts first if you want the parent scaffold of those.

    Input molecules are never modified.

    :param mols: List of OEMolBase molecules.
    :param scaffold: ``"framework"`` for the classic Bemis-Murcko scaffold --
        ring systems plus their linkers -- or ``"generic"`` to reduce that
        framework to its topology, every heavy atom carbon and every bond
        single. Case-insensitive.
    :param num_threads: Worker threads; 0 auto-detects hardware concurrency. An
        over-large value is clamped, not refused. Extraction runs serially
        whatever this value is when the OpenEye memory-pool mode is not
        thread-safe.
    :returns: Scaffold SMILES, one per molecule, in input order.
    :raises TypeError: If ``scaffold`` is not a string, ``num_threads`` is not
        an integer, or ``mols`` is not a list of OEMolBase molecules.
    :raises ValueError: If ``scaffold`` is not ``"framework"`` or ``"generic"``,
        ``mols`` is empty, or ``num_threads`` is negative or exceeds a
        ``size_t``.
    :raises RuntimeError: If a ring-containing molecule yields no framework.
        The message names the first such molecule's index.
    """
    options = _murcko_options(mols, scaffold, num_threads, "murcko_scaffolds")
    return list(_murcko_scaffolds(mols, options))


def murcko(mols, *, scaffold="framework", num_threads=0):
    """
    Cluster molecules by Bemis-Murcko scaffold identity.

    Two molecules share a cluster exactly when their scaffolds canonicalize to
    the same non-empty SMILES. Unlike :func:`butina` and the other
    distance-driven algorithms, this one takes molecules rather than a distance
    matrix -- there is no distance in it -- following the same rule as
    :func:`bitbirch`, whose first argument is a fingerprint batch.

    Labels come from sorting the distinct non-empty scaffold strings, so two
    runs over permuted inputs give the same label to the same scaffold. Acyclic
    molecules are noise: they carry ``-1``, appear in ``scaffolds`` as ``''``,
    and belong to no cluster.

    Every normalization and caveat on :func:`murcko_scaffolds` applies here
    unchanged, including the salt-stripping one.

    :param mols: List of OEMolBase molecules.
    :param scaffold: ``"framework"`` (default) or ``"generic"``;
        case-insensitive. See :func:`murcko_scaffolds`.
    :param num_threads: Worker threads; 0 auto-detects hardware concurrency.
    :returns: MurckoResult with labels, clusters, per-item ``scaffolds`` and the
        sorted ``cluster_scaffolds``.
    :raises TypeError: If ``scaffold`` is not a string, ``num_threads`` is not
        an integer, or ``mols`` is not a list of OEMolBase molecules.
    :raises ValueError: If ``scaffold`` is not ``"framework"`` or ``"generic"``,
        ``mols`` is empty, or ``num_threads`` is negative or exceeds a
        ``size_t``.
    :raises RuntimeError: If a ring-containing molecule yields no framework.
    """
    options = _murcko_options(mols, scaffold, num_threads, "murcko")
    result = _murcko_cluster(mols, options)
    return MurckoResult(
        result.Labels(),
        result.Members(),
        scaffolds=result.Scaffolds(),
        cluster_scaffolds=result.ClusterScaffolds(),
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
    options.singly = _flag(singly, "singly")
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
    fit_options.singly = _flag(singly, "singly")
    fit_options.mode = _bitbirch_mode(mode)
    fit_options.num_threads = int(num_threads)

    options = BitBirchRefinementOptions()
    options.fit_options = fit_options
    options.redistribute_largest_cluster = _flag(
        redistribute_largest_cluster, "redistribute_largest_cluster")
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


class ClusterRecord(NamedTuple):
    """One row of a :class:`ClusterReport`'s per-cluster table.

    ``label`` and ``nearest_cluster`` are ordinals into the clustering's member
    lists, not values read out of its label vector. The two are the same integer
    under ``cluster_report``'s partition precondition; naming which one is
    authoritative removes the ambiguity rather than relying on the coincidence.

    A ``NamedTuple`` rather than a plain class so the table is immutable without
    a hand-written ``__setattr__`` and feeds ``pandas.DataFrame(report.records)``
    directly, which is what a per-cluster table is for. pandas is not a
    dependency.

    Every field defaults, mirroring the C++ struct value for value, so that a
    default-constructed record reads the same from either side. The four
    undefined-valued fields default to NaN rather than 0.0: a zero
    ``nearest_cluster_distance`` would say another cluster sits at zero
    distance, and 0.0 is the real singleton value only for ``radius`` and
    ``diameter``.
    """

    label: int = 0
    size: int = 0
    representative: int = 0
    mean_intra_distance: float = float("nan")
    median_intra_distance: float = float("nan")
    radius: float = 0.0
    diameter: float = 0.0
    mean_representative_distance: float = 0.0
    nearest_cluster: int = -1
    nearest_cluster_distance: float = float("nan")
    silhouette: float = float("nan")
    boundary_violations: int = 0


class ClusterReportRequested(NamedTuple):
    """Which optional computations a :class:`ClusterReport`'s caller asked for.

    Records the request, not the outcome. A requested metric whose value is
    undefined still reports ``True`` here, so a NaN can be read unambiguously:
    ``False`` means nobody asked, ``True`` with NaN means asked and undefined.
    """

    pair_rank_indices: bool
    per_cluster_records: bool


class ClusterReport:
    """Read-only clustering-quality scorecard.

    Provides 27 scalar quality metrics plus threshold-coverage curves for
    evaluating clustering results. Key metrics include silhouette (separation),
    Dunn index (compactness vs. isolation), size Gini coefficient (imbalance),
    median radius/diameter, and coverage fractions at user-specified thresholds.

    Users should construct reports via the :func:`cluster_report` function,
    not by calling ``__init__`` directly.

    All metrics are exposed as read-only properties. Undefined metrics are NaN.
    ``c_index`` and ``baker_hubert_gamma`` are computed only when
    ``cluster_report`` is called with ``compute_pair_rank_indices``, so their
    NaN carries a second meaning; read :attr:`requested` to tell "nobody asked"
    apart from "asked and undefined".
    Vector metrics (``coverage_thresholds``, ``coverage_at``,
    ``noise_coverage_at``, ``records``) are tuples.
    """

    _SCALAR_FIELDS = (
        "num_samples", "num_clusters", "num_noise", "num_singletons",
        "noise_fraction", "singleton_fraction", "largest_cluster_fraction",
        "cluster_size_median", "cluster_size_p90", "size_gini", "size_entropy",
        "mean_intra_distance", "median_intra_distance", "median_radius",
        "p95_diameter", "silhouette", "dunn_index", "boundary_violations",
        "median_medoid_member_distance", "representative_redundancy",
        "calinski_harabasz_medoid",
        "davies_bouldin_medoid",
        "dunn_mean_separation_mean_diameter",
        "dunn_medoid_separation_medoid_spread",
        "point_biserial",
        "c_index",
        "baker_hubert_gamma",
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
        object.__setattr__(
            self, "_noise_coverage_at",
            tuple(float(v) for v in native_report.noise_coverage_at))
        object.__setattr__(
            self, "_records",
            tuple(
                ClusterRecord(
                    label=record.label,
                    size=record.size,
                    representative=record.representative,
                    mean_intra_distance=record.mean_intra_distance,
                    median_intra_distance=record.median_intra_distance,
                    radius=record.radius,
                    diameter=record.diameter,
                    mean_representative_distance=(
                        record.mean_representative_distance),
                    nearest_cluster=record.nearest_cluster,
                    nearest_cluster_distance=record.nearest_cluster_distance,
                    silhouette=record.silhouette,
                    boundary_violations=record.boundary_violations,
                )
                for record in native_report.records
            ))
        object.__setattr__(
            self, "_requested",
            ClusterReportRequested(
                pair_rank_indices=bool(
                    native_report.requested.pair_rank_indices),
                per_cluster_records=bool(
                    native_report.requested.per_cluster_records),
            ))
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

    @property
    def noise_coverage_at(self) -> tuple[float, ...]:
        """Coverage over noise points only, parallel to ``coverage_thresholds``.

        Same length as :attr:`coverage_at`: the threshold count when the
        clustering has at least one cluster, and empty when it has none. Every
        entry is NaN when the clustering has no noise -- not 0.0, which would
        read as "no noise point is covered".

        :returns: One fraction per coverage threshold.
        """
        return self._noise_coverage_at

    @property
    def records(self) -> tuple["ClusterRecord", ...]:
        """The per-cluster table, in member-list order.

        Empty unless ``compute_per_cluster_records`` was set. Read
        :attr:`requested` to tell "nobody asked" apart from "no clusters".

        :returns: One :class:`ClusterRecord` per cluster.
        """
        return self._records

    @property
    def requested(self) -> "ClusterReportRequested":
        """What the caller asked for, independent of what was computable.

        :returns: The two opt-in flags as passed.
        """
        return self._requested

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

    # Scalar metric -> the ClusterReportRequested flag that gates it. A metric
    # named here reads as None when its flag is False, which is what keeps
    # "nobody asked" apart from the NaN of "asked and undefined".
    _OPT_IN_SCALARS: ClassVar[dict[str, str]] = {
        "c_index": "pair_rank_indices",
        "baker_hubert_gamma": "pair_rank_indices",
    }

    def to_table(self):
        """Return rows ``(metric_name, *values)`` -- one value per report.

        Scalar-metric rows come first, then ``requested_pair_rank_indices``,
        then coverage rows aligned by threshold value across all reports, then
        the matching noise-coverage rows.

        A cell is ``None`` when that report never asked the question: an opt-in
        metric it did not request, a threshold it did not use, or a threshold it
        does carry but answered nothing at, as when the clustering has no
        clusters and both coverage curves come back empty. ``nan`` keeps
        its single meaning of asked-and-undefined. The two were previously
        indistinguishable, which is the collision this separates.
        """
        def _scalar_cell(report, name):
            flag = self._OPT_IN_SCALARS.get(name)
            if flag is not None and not getattr(report.requested, flag):
                return None
            return getattr(report, name)

        rows = []
        for name in ClusterReport._SCALAR_FIELDS:
            rows.append((name, *(_scalar_cell(r, name) for r in self._reports)))

        # Placed directly beneath the two rows it summarises, c_index and
        # baker_hubert_gamma being the last two scalar fields. Those two now
        # carry the unasked state themselves, as None, so this row states in
        # one line what their None means rather than being the only way to see
        # it. There is deliberately no matching row for per_cluster_records:
        # that flag governs ``records``, which this table does not carry, so
        # the row would be inert.
        rows.append((
            "requested_pair_rank_indices",
            *(r.requested.pair_rank_indices for r in self._reports),
        ))

        def _coverage_lookup(report, threshold, curve_name):
            """This report's value at ``threshold`` on ``curve_name``, or None.

            The ``i < len(curve)`` bound is what turns an empty curve -- a
            report with no clusters, which carries its thresholds but has no
            answers -- into "not asked" rather than a stale value.
            """
            curve = getattr(report, curve_name)
            for i, t in enumerate(report.coverage_thresholds):
                if t == threshold and i < len(curve):
                    return curve[i]
            return None

        thresholds = sorted(
            set().union(*(r.coverage_thresholds for r in self._reports)))
        for curve_name in ("coverage_at", "noise_coverage_at"):
            for threshold in thresholds:
                rows.append((
                    f"{curve_name}[{threshold}]",
                    *(_coverage_lookup(r, threshold, curve_name)
                      for r in self._reports),
                ))
        return rows

    def __repr__(self):
        return _format_metric_table(self.to_table(), self._column_labels())


def _format_metric_table(rows, labels):
    """Render metric rows as the fixed-width table both scorecard reprs print.

    Shared by :class:`ClusterReportComparison` and :class:`PartitionAgreement`,
    which print the same table at different widths rather than being related
    types. A row is its metric name followed by one value per column, and None
    reads as "--" so a metric nobody asked for is distinguishable from one that
    came back undefined.

    :param rows: Sequence of ``(name, *values)``, one value per label.
    :param labels: Column headers, left to right.
    :returns: The rendered table as a single newline-joined str.
    """
    def _fmt(value):
        if value is None:
            return "--"
        if isinstance(value, float):
            return f"{value:.4g}"
        return str(value)

    metric_width = max([len("metric")] + [len(row[0]) for row in rows])
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


class PartitionAgreementRequested(NamedTuple):
    """Which optional computations a :class:`PartitionAgreement`'s caller asked for.

    Records the request, not the outcome. A requested metric whose value is
    undefined still reports ``True`` here, so a NaN can be read unambiguously:
    ``False`` means nobody asked, ``True`` with NaN means asked and undefined.
    """

    adjusted_mutual_information: bool


class PartitionAgreement:
    """Agreement between two labelings of the same samples.

    Construct these with :func:`partition_agreement` or
    :func:`scaffold_agreement` rather than by calling ``__init__`` directly.
    The attributes are read-only by convention; assigning to one changes this
    object and nothing else.

    All seven metrics are NaN when fewer than two samples survive noise
    handling; the counts keep their real values. The six always-computed
    metrics are 1.0 when the two partitions are identical after noise
    handling -- ``adjusted_mutual_information`` joins them only when it was
    requested, and otherwise stays NaN with ``requested`` false. The per-field
    notes below cover the remaining cases.

    :ivar num_samples: Samples entering the contingency table. Equals the input
        length except under ``noise="excluded"``, which can shrink it.
    :ivar num_clusters_a: Distinct clusters on side A after noise handling.
    :ivar num_clusters_b: Distinct clusters on side B after noise handling.
    :ivar adjusted_rand_index: Hubert-Arabie adjusted Rand index. 1.0 is exact
        agreement, 0.0 is the value expected by chance, negative is worse than
        chance.
    :ivar fowlkes_mallows: Geometric mean of pair precision and pair recall.
        NaN when either side is all singletons and the partitions differ;
        scikit-learn reports 0.0 there.
    :ivar normalized_mutual_information: Mutual information over the arithmetic
        mean of the two entropies. Bitwise equal to ``v_measure``.
    :ivar homogeneity: ``MI / H(a)``. Asymmetric. NaN when side A is a single
        cluster and the partitions differ; scikit-learn reports 1.0 there.
    :ivar completeness: ``MI / H(b)``. Asymmetric. NaN when side B is a single
        cluster and the partitions differ; scikit-learn reports 1.0 there.
    :ivar v_measure: The beta = 1 V-measure, assigned from the same value as
        ``normalized_mutual_information``. Stays defined where the harmonic
        mean of homogeneity and completeness does not.
    :ivar adjusted_mutual_information: Mutual information corrected for chance.
        NaN unless ``adjusted_mutual_information=True`` was passed; read
        :attr:`requested` to tell "nobody asked" from "asked and undefined".
        Its accuracy is limited where the denominator -- mean entropy minus
        expected mutual information -- approaches zero, which happens when both
        partitions are close to all-singleton. Two 1.5-million-sample
        partitions differing by one merged pair have a true value of zero and
        report about 0.035. Partitions whose denominator is order one are
        unaffected. The limit is inherent to computing the correction in double
        precision -- scikit-learn shares it -- rather than a property of this
        implementation.
    :ivar requested: A :class:`PartitionAgreementRequested`.
    """

    # Field order mirrors the C++ PartitionAgreement struct and the :ivar: list
    # above; alphabetizing would desynchronize both from the native layout they
    # document.
    __slots__ = (  # noqa: RUF023
        "num_samples",
        "num_clusters_a",
        "num_clusters_b",
        "adjusted_rand_index",
        "fowlkes_mallows",
        "normalized_mutual_information",
        "homogeneity",
        "completeness",
        "v_measure",
        "adjusted_mutual_information",
        "requested",
    )

    # The rows to_table() emits, in order. Excludes `requested`, which gates
    # the rows rather than being one.
    _TABLE_FIELDS: ClassVar[tuple[str, ...]] = (
        "num_samples",
        "num_clusters_a",
        "num_clusters_b",
        "adjusted_rand_index",
        "fowlkes_mallows",
        "normalized_mutual_information",
        "homogeneity",
        "completeness",
        "v_measure",
        "adjusted_mutual_information",
    )

    # Scalar metric -> the PartitionAgreementRequested flag that gates it. A
    # metric named here reads as None when its flag is False, which is what
    # keeps "nobody asked" apart from the NaN of "asked and undefined".
    _OPT_IN_SCALARS: ClassVar[dict[str, str]] = {
        "adjusted_mutual_information": "adjusted_mutual_information",
    }

    def __init__(self, native_agreement):
        """
        Construct an agreement scorecard from the native result.

        :param native_agreement: Native agreement object returned by the
            extension.
        """
        self.num_samples = int(native_agreement.num_samples)
        self.num_clusters_a = int(native_agreement.num_clusters_a)
        self.num_clusters_b = int(native_agreement.num_clusters_b)
        self.adjusted_rand_index = float(native_agreement.adjusted_rand_index)
        self.fowlkes_mallows = float(native_agreement.fowlkes_mallows)
        self.normalized_mutual_information = float(
            native_agreement.normalized_mutual_information)
        self.homogeneity = float(native_agreement.homogeneity)
        self.completeness = float(native_agreement.completeness)
        self.v_measure = float(native_agreement.v_measure)
        self.adjusted_mutual_information = float(
            native_agreement.adjusted_mutual_information)
        self.requested = PartitionAgreementRequested(
            adjusted_mutual_information=bool(
                native_agreement.requested.adjusted_mutual_information),
        )

    def to_table(self):
        """Return ``(metric_name, value)`` rows, one per reported quantity.

        A value is ``None`` when nobody asked the question -- an opt-in metric
        whose :attr:`requested` flag is ``False`` -- so ``nan`` keeps its single
        meaning of asked-and-undefined.

        :returns: A list of ``(str, value)`` pairs.
        """
        rows = []
        for name in self._TABLE_FIELDS:
            flag = self._OPT_IN_SCALARS.get(name)
            if flag is not None and not getattr(self.requested, flag):
                rows.append((name, None))
            else:
                rows.append((name, getattr(self, name)))
        return rows

    def __repr__(self):
        return _format_metric_table(self.to_table(), ("value",))


def _noise_handling(noise):
    noise_map = {
        "singletons": _oecluster.NoiseHandling_Singletons,
        "grouped": _oecluster.NoiseHandling_Grouped,
        "excluded": _oecluster.NoiseHandling_Excluded,
    }
    key = str(noise).lower()
    if key not in noise_map:
        raise ValueError(
            f"Unknown noise handling: {noise!r}; expected 'singletons', "
            f"'grouped' or 'excluded'")
    return noise_map[key]


def _agreement_options(noise, adjusted_mutual_information):
    options = _oecluster.PartitionAgreementOptions()
    options.noise_handling = _noise_handling(noise)
    options.compute_adjusted_mutual_information = _flag(
        adjusted_mutual_information, "adjusted_mutual_information")
    return options


def _agreement_labels(value, argument_name):
    """Coerce a clustering result or a sequence of ints to a native IntVector."""
    labels = getattr(value, "labels", value)
    # A Mapping iterates its keys, so {0: "a", 1: "b"} would score the keys and
    # report a plausible number for a labeling the caller never passed. The
    # values are the likelier intent, but guessing between the two is worse
    # than refusing.
    if isinstance(labels, collections.abc.Mapping):
        raise TypeError(
            f"{argument_name} must be a clustering result or a sequence of "
            f"ints, not a mapping")
    vector = _oecluster.IntVector()
    try:
        for label in labels:
            # operator.index() accepts int, bool, numpy integers, and rejects
            # float/Decimal/str/None. int() would silently truncate floats and
            # coerce strings, publishing a wrong metric instead of raising.
            vector.push_back(operator.index(label))
    except OverflowError as error:
        # operator.index() admits any Python int, so a label too wide for the
        # native std::vector<int> only fails inside push_back. Folding it into
        # the clause below would report a type problem, which is false.
        raise ValueError(
            f"{argument_name} must contain labels that fit a 32-bit signed "
            f"int") from error
    except (TypeError, ValueError) as error:
        raise TypeError(
            f"{argument_name} must be a clustering result or a sequence of "
            f"ints") from error
    return vector


def _string_vector(value, argument_name, noun):
    """Coerce a sequence of strings to a native StringVector.

    :param value: The caller's sequence.
    :param argument_name: Name of the argument, for the messages.
    :param noun: Plural noun for what the strings are, e.g. ``"scaffold
        strings"``. Two entry points annotate samples with strings and each
        wants its own word in the refusal; the rules are identical.
    :returns: A native StringVector.
    :raises TypeError: If the value is a bare str, a mapping, not iterable, or
        yields a non-str.
    """
    # A bare str is iterable, so without this guard a single annotation string
    # would silently become one annotation per character.
    if isinstance(value, str):
        raise TypeError(
            f"{argument_name} must be a sequence of {noun}, not a single str")
    # As in _agreement_labels: a Mapping's iteration yields its keys, which
    # would score something the caller did not pass.
    if isinstance(value, collections.abc.Mapping):
        raise TypeError(
            f"{argument_name} must be a sequence of {noun}, not a mapping")
    vector = _oecluster.StringVector()
    try:
        iterator = iter(value)
    except TypeError as error:
        raise TypeError(
            f"{argument_name} must be a sequence of {noun}") from error
    for label in iterator:
        # Require actual strings. str() would coerce None, numbers, etc.,
        # publishing a wrong metric instead of raising.
        if not isinstance(label, str):
            raise TypeError(f"{argument_name} must be a sequence of {noun}")
        try:
            # A str SWIG cannot encode to UTF-8 -- a lone surrogate, which
            # surrogateescape decoding of a mis-encoded file produces -- is
            # rejected here rather than by the isinstance check above, and
            # would otherwise escape naming the internal container type. The
            # check stays outside this guard so its raise is not self-caught.
            vector.push_back(label)
        except TypeError as error:
            raise TypeError(
                f"{argument_name} must be a sequence of {noun}") from error
    return vector


def _agreement_scaffolds(value, argument_name):
    """Coerce a sequence of scaffold strings to a native StringVector."""
    return _string_vector(value, argument_name, "scaffold strings")


def cluster_report(result, distance_matrix, *, preset="default",
                   coverage_thresholds=None, boundary_threshold=None,
                   representative_method="medoid",
                   treat_noise_as_singletons=True, num_threads=0,
                   compute_pair_rank_indices=False,
                   compute_per_cluster_records=False,
                   allow_nonmetric=False):
    """
    Compute a method-agnostic clustering-quality report.

    :param result: A clustering result (e.g. from :func:`butina`/:func:`dbscan`).
    :param distance_matrix: Complete SymmetricDistanceMatrix for the same items.
    :param preset: Threshold preset: "default", "tight", or "diversity".
    :param coverage_thresholds: Optional override for coverage distances; NaN
        is refused, because no reported value can be matched back to it.
    :param boundary_threshold: Optional override for the boundary-violation
        distance; NaN is refused, because no distance compares against it.
    :param representative_method: Centrality method for medoid selection:
        "medoid" (default), "minimax", or "weighted_medoid". The
        "highest_neighborhood" method is not supported.
    :param treat_noise_as_singletons: Fold noise into singleton accounting.
    :param num_threads: Reserved for parallel-safe computation.
    :param compute_pair_rank_indices: Compute ``c_index`` and
        ``baker_hubert_gamma``. Both are read off two sorted arrays holding
        every pairwise distance among clustered points,
        ``Nc * (Nc - 1) / 2`` doubles between them -- roughly 400 MB at
        ``Nc = 10,000`` and 10 GB at ``Nc = 50,000``. Only the between-cluster
        array is this flag's own cost; the within-cluster one is built on every
        call, because ``median_intra_distance`` is taken over it. So the flag
        adds nothing to a single-cluster result and nearly the whole figure to
        one with small clusters, and it is off by default for the second case.
        Raises ``MemoryError`` if the allocation fails.
    :param compute_per_cluster_records: Populate :attr:`ClusterReport.records`.
        Off by default: the stage buffers the largest cluster's pairwise
        distances to take their median, ``n * (n - 1) / 2`` doubles, and the
        median is taken over a copy of that buffer -- roughly 400 MB for the
        buffer and 400 MB again for the copy, transiently, at ``n = 10,000``.
    :param allow_nonmetric: Score anyway when the distances are known not to
        satisfy the triangle inequality. Does not override the refusals for
        similarity-valued or non-finite matrices.
    :returns: A ClusterReport.
    :raises TypeError: If result/distance_matrix have the wrong type, or
        treat_noise_as_singletons, compute_pair_rank_indices,
        compute_per_cluster_records or allow_nonmetric is not a bool or
        ``numpy.bool_``.
    :raises ValueError: If a preset/method/threshold is invalid, the result and
        the matrix cover different numbers of samples, the matrix uses sparse
        storage, or the matrix is not a metric -- and additionally: any
        non-finite entry anywhere in the matrix is refused before the report is
        computed, whether or not that pair reaches a reported value.
    :raises RuntimeError: If the distance matrix cannot provide complete
        distances, or if the result's labels and members describe different
        partitions.
    :raises MemoryError: If the pair-rank stage cannot allocate its distance
        arrays, or would exceed its couple counter.
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
            # NaN slips past the negative check, and every later comparison
            # against it is false: the native computes a coverage value, the
            # comparison table's threshold match never finds it, and the cell
            # publishes as None, which in that table means nobody asked.
            # Refusing here is the only place the caller can still be told.
            if math.isnan(v):
                raise ValueError("coverage thresholds must not be NaN")
            if v < 0.0:
                raise ValueError("coverage thresholds must be non-negative")
            coverage_vector.push_back(v)

    bt = None
    if boundary_threshold is not None:
        bt = float(boundary_threshold)
        if math.isnan(bt):
            raise ValueError("boundary_threshold must not be NaN")
        if bt < 0.0:
            raise ValueError("boundary_threshold must be non-negative")

    method_key, native_method = _representative_method(representative_method)
    if method_key == "highest_neighborhood":
        raise ValueError(
            "cluster_report does not support the 'highest_neighborhood' "
            "representative method; use 'medoid', 'minimax', or 'weighted_medoid'")

    # Gated on the same terms as the two bool keywords below -- see that comment
    # for why numpy.bool_ is admitted and why the message carries the value --
    # but checked here rather than beside them, because this block reports its
    # arguments in signature order. A call that is wrong in two ways names the
    # keyword the caller wrote first, and treat_noise_as_singletons is declared
    # ahead of num_threads.
    if not isinstance(treat_noise_as_singletons, (bool, np.bool_)):
        raise TypeError(
            "treat_noise_as_singletons must be True or False, "
            f"not {type(treat_noise_as_singletons).__name__} "
            f"({treat_noise_as_singletons!r})."
        )

    num_threads_int = int(num_threads)

    # num_threads reaches a size_t option field, where a negative value raises
    # OverflowError below the gate. Zero stays legal: it means "choose for me".
    if num_threads_int < 0:
        raise ValueError("num_threads must be non-negative")

    # numpy.bool_ is admitted on the same terms as allow_nonmetric, which the
    # gate below has always taken: one call must not apply two admissibility
    # rules to its bool keywords, and np.bool_ is what arr.any() and every
    # comparison of numpy scalars returns. Truthiness stays refused, because
    # compute_pair_rank_indices="no" reads to a caller as off while switching
    # an expensive stage on. The message carries the value as well as the type
    # for the same reason the gate's does: numpy.bool_.__name__ is itself
    # "bool", so a bare type name would read as "must be a bool, got bool".
    if not isinstance(compute_pair_rank_indices, (bool, np.bool_)):
        raise TypeError(
            "compute_pair_rank_indices must be True or False, "
            f"not {type(compute_pair_rank_indices).__name__} "
            f"({compute_pair_rank_indices!r})."
        )
    if not isinstance(compute_per_cluster_records, (bool, np.bool_)):
        raise TypeError(
            "compute_per_cluster_records must be True or False, "
            f"not {type(compute_per_cluster_records).__name__} "
            f"({compute_per_cluster_records!r})."
        )

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
    # Coerced like treat_noise_as_singletons above, and for a harder reason:
    # the SWIG bool setter takes only a Python bool, so an admitted numpy.bool_
    # would otherwise fail here with a message naming a generated setter.
    options.compute_pair_rank_indices = bool(compute_pair_rank_indices)
    options.compute_per_cluster_records = bool(compute_per_cluster_records)

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


def partition_agreement(a, b, *, noise="singletons",
                        adjusted_mutual_information=False):
    """
    Score the agreement between two labelings of the same samples.

    Takes labels and nothing else -- no distance matrix, unlike
    :func:`cluster_report`. Side A is always the first argument:
    ``homogeneity`` is ``MI / H(a)`` and ``completeness`` is ``MI / H(b)``, and
    those two exchange values when the arguments are swapped, as do
    ``num_clusters_a`` and ``num_clusters_b``; every other metric is
    symmetric, though ``adjusted_mutual_information`` only to within rounding,
    because a swap exchanges the two marginal values inside each expected-MI
    term and so changes the floating-point order they evaluate in.

    Undefined metrics are NaN rather than a convention, which diverges from
    scikit-learn on five degenerate inputs; ``docs/python-api.md`` tabulates
    them.

    :param a: A clustering result, or a sequence of ints fitting the native
        32-bit signed label type (list, tuple, or a numpy integer array). The
        reference labeling.
    :param b: The candidate labeling, in the same forms.
    :param noise: How negatively-labelled samples enter the contingency table:
        ``"singletons"`` (the default, each noise sample is its own cluster),
        ``"grouped"`` (each side's noise forms one cluster, which is how
        scikit-learn reads a -1 label), or ``"excluded"`` (a sample noisy on
        either side is dropped from both).
    :param adjusted_mutual_information: Also compute AMI. Off by default: the
        other six metrics are essentially free once the contingency table is
        built, while the expected-MI correction adds an O(N) log-factorial
        table and a sum over pairs of distinct cluster sizes.
    :returns: A :class:`PartitionAgreement`.
    :raises TypeError: If either argument is neither a clustering result nor a
        sequence of ints.
    :raises ValueError: If the labelings differ in length, either is empty, a
        label does not fit the native 32-bit signed label type, or ``noise`` is
        not one of the three accepted strings.
    """
    labels_a = _agreement_labels(a, "a")
    labels_b = _agreement_labels(b, "b")
    if len(labels_a) == 0 or len(labels_b) == 0:
        raise ValueError("partition_agreement() requires non-empty labelings")
    if len(labels_a) != len(labels_b):
        raise ValueError(
            f"b has {len(labels_b)} entries but a has {len(labels_a)} samples")
    options = _agreement_options(noise, adjusted_mutual_information)
    return PartitionAgreement(
        _oecluster.partition_agreement(labels_a, labels_b, options))


def scaffold_agreement(result, scaffold_labels, *, noise="singletons",
                       adjusted_mutual_information=False):
    """
    Score a clustering against a per-sample scaffold annotation.

    The clustering is side A and the annotation is side B, following the same
    positional rule as :func:`partition_agreement`. So ``completeness`` carries
    the scaffold-purity reading -- whether each cluster's members share a
    single scaffold -- and ``homogeneity`` carries its transpose, whether each
    scaffold landed in a single cluster.

    An empty scaffold string is missing data, not a category, and follows
    ``noise`` exactly as a negative label does on the integer side.

    :param result: A clustering result, or a sequence of ints fitting the
        native 32-bit signed label type.
    :param scaffold_labels: One scaffold string per sample.
    :param noise: As for :func:`partition_agreement`.
    :param adjusted_mutual_information: As for :func:`partition_agreement`.
    :returns: A :class:`PartitionAgreement`.
    :raises TypeError: If ``result`` is not a clustering result or sequence of
        ints, or ``scaffold_labels`` is not a sequence of strings.
    :raises ValueError: If the two differ in length, either is empty, a label
        does not fit the native 32-bit signed label type, or ``noise`` is not
        one of the three accepted strings.
    """
    labels = _agreement_labels(result, "result")
    scaffolds = _agreement_scaffolds(scaffold_labels, "scaffold_labels")
    if len(labels) == 0 or len(scaffolds) == 0:
        raise ValueError(
            "scaffold_agreement() requires a non-empty clustering and a "
            "non-empty scaffold annotation")
    if len(labels) != len(scaffolds):
        raise ValueError(
            f"scaffold_labels has {len(scaffolds)} entries but the clustering "
            f"has {len(labels)} samples")
    options = _agreement_options(noise, adjusted_mutual_information)
    return PartitionAgreement(
        _oecluster.scaffold_agreement(labels, scaffolds, options))


class _Scorecard:
    """Fixed-width table rendering shared by the three SAR-coherence results.

    :class:`PartitionAgreement` is deliberately not folded onto this. Its
    ``to_table`` consults ``_OPT_IN_SCALARS`` to keep "nobody asked" apart from
    "asked and undefined", and none of these three has an opt-in scalar, so
    sharing a base would put a branch in it for a case one subclass has.
    """

    __slots__ = ()

    # The rows to_table() emits, in order. Each subclass overrides it.
    _TABLE_FIELDS: ClassVar[tuple[str, ...]] = ()

    def to_table(self):
        """Return ``(metric_name, value)`` rows, one per reported scalar.

        The per-row tables -- :attr:`SARCoherence.clusters` and
        :attr:`Modelability.classes` -- are not rows here. They are read as
        sequences; a table nested in a table cell does not render.

        :returns: A list of ``(str, value)`` pairs.
        """
        return [(name, getattr(self, name)) for name in self._TABLE_FIELDS]

    def __repr__(self):
        return _format_metric_table(self.to_table(), ("value",))


class ClusterActivity(NamedTuple):
    """One row of a :class:`SARCoherence`'s per-cluster table.

    A ``NamedTuple`` for the same reasons as :class:`ClusterRecord`: immutable
    without a hand-written ``__setattr__``, and it feeds
    ``pandas.DataFrame(coherence.clusters)`` directly. pandas is not a
    dependency.

    Every field defaults, mirroring the C++ struct value for value. The two
    statistics default to NaN rather than 0.0, because a zero mean activity is
    a real measurement rather than the absence of one.

    :ivar label: The cluster's label. Under ``noise="grouped"`` the merged noise
        row reports -1; under ``noise="singletons"`` each noise sample keeps its
        own label, so several rows may share a label and are distinguished by
        position. Keying the table by label -- ``{row.label: row for row in
        coherence.clusters}`` -- silently discards all but the last of them.
    :ivar num_scored: Samples in this cluster with a finite activity.
    :ivar mean_activity: Their mean activity.
    :ivar stddev_activity: Their population standard deviation, NaN for every
        row whose num_scored is below 2, whatever kind of row it is: a spread
        over one sample is not defined. That covers the singleton noise rows
        and both routes an ordinary cluster takes there -- holding a single
        member, or having missing activity thin it down to one.
    """

    label: int = 0
    num_scored: int = 0
    mean_activity: float = float("nan")
    stddev_activity: float = float("nan")


class SARCoherence(_Scorecard):
    """How much of an activity's variance a clustering explains.

    Construct these with :func:`sar_coherence` rather than by calling
    ``__init__`` directly. The attributes are read-only by convention;
    assigning to one changes this object and nothing else.

    :ivar num_samples: Length of the input activity vector.
    :ivar num_scored: Samples with a finite activity that survived noise
        handling. Both effect sizes are NaN when this is below two.
    :ivar num_clusters: Clusters retaining at least one scored sample.
    :ivar eta_squared: ``SS_between / SS_total``: the share of activity
        variance the clustering accounts for. It rises with the cluster count
        even when activity is independent of the labels, by roughly
        ``(K - 1) / (n - 1)``, so it does not compare across clusterings that
        differ in cluster count. NaN when the activity has no variance.
    :ivar omega_squared: The chance-corrected counterpart, which does compare
        across cluster counts. Negative when the clustering explains less than
        chance would; that is a reading, not an error. NaN wherever
        :attr:`eta_squared` is, and also when every scored sample is its own
        cluster, which leaves no within-cluster variance to correct against.
    :ivar clusters: A tuple of :class:`ClusterActivity`, one per cluster with a
        scored member, ordered by where the cluster's label first appears among
        the scored samples rather than in the raw input. Under
        ``noise="singletons"`` several rows carry the label -1; see
        :class:`ClusterActivity`.
    """

    # Field order mirrors the C++ SARCoherence struct and the :ivar: list
    # above; alphabetizing would desynchronize both from the native layout.
    __slots__ = (  # noqa: RUF023
        "num_samples",
        "num_scored",
        "num_clusters",
        "eta_squared",
        "omega_squared",
        "clusters",
    )

    _TABLE_FIELDS: ClassVar[tuple[str, ...]] = (
        "num_samples",
        "num_scored",
        "num_clusters",
        "eta_squared",
        "omega_squared",
    )

    def __init__(self, native_coherence):
        """
        Construct a coherence scorecard from the native result.

        :param native_coherence: Native object returned by the extension.
        """
        self.num_samples = int(native_coherence.num_samples)
        self.num_scored = int(native_coherence.num_scored)
        self.num_clusters = int(native_coherence.num_clusters)
        self.eta_squared = float(native_coherence.eta_squared)
        self.omega_squared = float(native_coherence.omega_squared)
        # Every row is copied out into a plain Python record while the native
        # result is still alive. The member vector is owned by that result, so
        # a scorecard holding the vector itself would read as empty the moment
        # the caller let the native object go -- and for the common
        # ``sar_coherence(...).clusters`` spelling, that is before the caller
        # ever sees it.
        self.clusters = tuple(
            ClusterActivity(
                label=int(row.label),
                num_scored=int(row.num_scored),
                mean_activity=float(row.mean_activity),
                stddev_activity=float(row.stddev_activity),
            )
            for row in native_coherence.clusters
        )


class ActivityLandscape(_Scorecard):
    """Activity-cliff structure over a precomputed distance matrix.

    Construct these with :func:`activity_landscape`. The attributes are
    read-only by convention.

    :ivar num_samples: Length of the input activity vector.
    :ivar num_scored: Samples with a finite activity.
    :ivar num_pairs_scored: ``num_scored * (num_scored - 1) / 2``; zero below
        two scored samples.
    :ivar num_cliffs: Pairs at or below ``distance_threshold`` whose activity
        differs by at least ``activity_threshold``.
    :ivar cliff_density: ``num_cliffs / num_pairs_scored``. NaN when no pair
        was scored.
    :ivar num_zero_distance_pairs: Pairs at exactly zero distance. Their SALI
        is undefined, so they are counted here and left out of
        :attr:`max_sali` and :attr:`mean_sali` -- but they are still eligible to
        count as cliffs, and do whenever their activity difference reaches
        ``activity_threshold``, because two identical structures with different
        activities are the sharpest cliff there is. The comparison is ``>=``,
        so a zero-distance pair whose activities agree is counted here and is
        also a cliff exactly when ``activity_threshold`` is zero.
    :ivar max_sali: Largest ``|activity difference| / distance`` over the pairs
        at non-zero distance. NaN when there are none.
    :ivar mean_sali: Their mean. NaN on the same condition.
    :ivar rmodi: The regression modelability index: the fraction of scored
        samples whose nearest neighbour inside the activity band is strictly
        closer than their nearest neighbour outside it. The band reaches
        ``rmodi_delta`` standard deviations either side of a molecule's own
        activity, so it spans twice that. NaN below two scored samples.
    :ivar activity_stddev: Population standard deviation of the scored
        activity, the quantity the band is measured in. NaN below two scored
        samples.
    """

    __slots__ = (  # noqa: RUF023
        "num_samples",
        "num_scored",
        "num_pairs_scored",
        "num_cliffs",
        "cliff_density",
        "num_zero_distance_pairs",
        "max_sali",
        "mean_sali",
        "rmodi",
        "activity_stddev",
    )

    # Exact rather than a second copy: every one of this scorecard's fields is
    # a scalar metric row, unlike the other two, whose row tables are slots and
    # not rows.
    _TABLE_FIELDS: ClassVar[tuple[str, ...]] = __slots__

    def __init__(self, native_landscape):
        """
        Construct a landscape scorecard from the native result.

        :param native_landscape: Native object returned by the extension.
        """
        self.num_samples = int(native_landscape.num_samples)
        self.num_scored = int(native_landscape.num_scored)
        self.num_pairs_scored = int(native_landscape.num_pairs_scored)
        self.num_cliffs = int(native_landscape.num_cliffs)
        self.cliff_density = float(native_landscape.cliff_density)
        self.num_zero_distance_pairs = int(
            native_landscape.num_zero_distance_pairs)
        self.max_sali = float(native_landscape.max_sali)
        self.mean_sali = float(native_landscape.mean_sali)
        self.rmodi = float(native_landscape.rmodi)
        self.activity_stddev = float(native_landscape.activity_stddev)


class ClassConcordance(NamedTuple):
    """One row of a :class:`Modelability`'s per-class table.

    A ``NamedTuple`` on the same terms as :class:`ClusterActivity`.

    :ivar label: The class string, as it appears in the input.
    :ivar num_members: Scored samples carrying it -- those with a non-empty
        class string.
    :ivar fraction_same_class: The fraction of those whose nearest scored
        neighbour shares the class. NaN when this is the only scored class,
        where no molecule has a neighbour that could differ.
    """

    label: str = ""
    num_members: int = 0
    fraction_same_class: float = float("nan")


class Modelability(_Scorecard):
    """Whether a descriptor separates the activity classes at all.

    Construct these with :func:`modelability`. The attributes are read-only by
    convention.

    :ivar num_samples: Length of the input class vector.
    :ivar num_scored: Samples with a non-empty class string.
    :ivar num_classes: Distinct non-empty class strings.
    :ivar modi: The unweighted mean of :attr:`ClassConcordance.
        fraction_same_class` over the classes. A low value says no classifier
        is likely to learn this dataset from this descriptor. NaN below two
        classes, where no molecule has a neighbour that could differ.
    :ivar classes: A tuple of :class:`ClassConcordance`, ordered by where each
        class string first appears among the scored samples. Only the empty
        string is unscored, so over the classes that exist this is input order.
    """

    __slots__ = (  # noqa: RUF023
        "num_samples",
        "num_scored",
        "num_classes",
        "modi",
        "classes",
    )

    _TABLE_FIELDS: ClassVar[tuple[str, ...]] = (
        "num_samples",
        "num_scored",
        "num_classes",
        "modi",
    )

    def __init__(self, native_modelability):
        """
        Construct a modelability scorecard from the native result.

        :param native_modelability: Native object returned by the extension.
        """
        self.num_samples = int(native_modelability.num_samples)
        self.num_scored = int(native_modelability.num_scored)
        self.num_classes = int(native_modelability.num_classes)
        self.modi = float(native_modelability.modi)
        # Copied out while the native result is alive, for the reason given in
        # SARCoherence.__init__.
        self.classes = tuple(
            ClassConcordance(
                label=str(row.label),
                num_members=int(row.num_members),
                fraction_same_class=float(row.fraction_same_class),
            )
            for row in native_modelability.classes
        )


def _activity_classes(value, argument_name):
    """Coerce a sequence of activity-class strings to a native StringVector."""
    return _string_vector(value, argument_name, "class strings")


def _activity_values(value, argument_name):
    """Coerce a sequence of activity measurements to a native DoubleVector.

    :param value: The caller's sequence; NaN marks a missing measurement.
    :param argument_name: Name of the argument, for the messages.
    :returns: A native DoubleVector.
    :raises TypeError: If the value is a bare str, a ``bytes``, a ``bytearray``
        or a ``memoryview``, a mapping, not iterable, or yields one of those
        four or an element ``float()`` refuses. Those three binary containers
        are named one by one rather than tested for as buffers, so a
        byte-format ``array.array`` or a numpy ``uint8`` column is iterated
        like any other sequence and read as the numbers it holds:
        ``array("B", b"12")`` scores ``[49.0, 50.0]``. Exporting a buffer is
        neither what admits a value nor what refuses one -- an ``mmap`` is
        refused because it yields ``bytes``, and a memoryview is refused
        whatever its format, so a float buffer is refused with the rest; pass
        ``np.asarray(view)`` to score one.
    :raises ValueError: If it yields a number too large to convert to a double,
        such as ``10 ** 1000``. A magnitude fault rather than a type fault, so
        it is not folded into the TypeError above. An infinity converts
        perfectly well and is refused in C++ instead.
    """
    if isinstance(value, str):
        raise TypeError(
            f"{argument_name} must be a sequence of floats, not a single str")
    # Refused whatever the buffer's format. bytes and bytearray iterate as
    # ints, so b"12" would score as the activities [49.0, 50.0]; a numeric
    # memoryview would convert correctly, but it is refused with them rather
    # than branching on .format -- pass np.asarray(view) to score a float
    # buffer. Three concrete types are named here rather than the buffer
    # protocol tested for, and that closes the reinterpretation for these three
    # only: a buffer test cannot tell a uint8 measurement column from misread
    # text, so array("B", b"12") and np.frombuffer(b"12", "u1") are accepted and
    # read as the numbers they hold. The three named types earn the refusal by
    # being the text-and-binary containers, which a caller holding one has
    # almost certainly not meant as a column of measurements.
    if isinstance(value, (bytes, bytearray, memoryview)):
        raise TypeError(
            f"{argument_name} must be a sequence of floats, not a bytes-like "
            "object")
    # As in _agreement_labels: a Mapping iterates its keys, so
    # {0: 5.4, 1: 6.1} would score the indices and report a plausible number.
    if isinstance(value, collections.abc.Mapping):
        raise TypeError(
            f"{argument_name} must be a sequence of floats, not a mapping")
    vector = _oecluster.DoubleVector()
    try:
        iterator = iter(value)
    except TypeError as error:
        raise TypeError(
            f"{argument_name} must be a sequence of floats") from error
    for item in iterator:
        # float() rather than operator.index(), which _agreement_labels uses:
        # an activity is a measurement, so 7 and 7.0 are the same input and a
        # numpy float has to pass. str, bytes, bytearray and memoryview are
        # refused explicitly because float() converts str, bytes, bytearray and
        # a byte-format memoryview -- a column read from a CSV without
        # conversion, or one whose entries survived a single layer of
        # deserialization, would otherwise score as numbers. These four are
        # named concretely rather than tested for as buffers, for the reason
        # the column guard gives, and they are the whole of the type test:
        # every other element is handed to float(), which reads a numpy scalar
        # numerically, parses a nested byte-format array.array as text -- such
        # an element of array("B", b"1") scores 1.0, the digit its byte spells,
        # rather than the 49 the array holds -- and refuses the rest. Among the
        # four, a numeric memoryview is the one spelling float() does not
        # convert; it is refused with them for consistency with that guard.
        if isinstance(item, (str, bytes, bytearray, memoryview)):
            raise TypeError(f"{argument_name} must be a sequence of floats")
        try:
            vector.push_back(float(item))
        except OverflowError as error:
            # float() admits any Python int until the cast itself, so a value
            # beyond double range only fails inside it, never at the type check
            # above. Folding it into the clause below would report a type
            # problem, which is false: the magnitude is the only thing wrong.
            raise ValueError(
                f"{argument_name} must contain values that fit a "
                f"double") from error
        except (TypeError, ValueError) as error:
            raise TypeError(
                f"{argument_name} must be a sequence of floats") from error
    return vector


def sar_coherence(result, activity, *, noise="excluded"):
    """
    Score how much of an activity's variance a clustering explains.

    Takes labels and measurements and nothing else -- no distance matrix,
    unlike :func:`activity_landscape` -- so it answers the question for a
    clustering computed any way at all, including one read in from elsewhere.

    Read ``omega_squared`` when comparing clusterings that differ in cluster
    count, and ``eta_squared`` only within a fixed count: the raw ratio rises
    with the number of clusters even when activity is independent of the
    labels.

    :param result: A clustering result, or a sequence of ints fitting the
        native 32-bit signed label type. Negative labels are noise.
    :param activity: One measurement per sample. NaN marks a missing
        measurement, and those samples are dropped rather than refused.
    :param noise: How negatively-labelled samples are treated: ``"excluded"``
        (the default, dropped), ``"grouped"`` (all noise forms one cluster) or
        ``"singletons"`` (each noise sample is its own cluster). The default
        differs from :func:`partition_agreement`'s, because a noise point
        promoted to a singleton has a group mean equal to its own value and so
        reads as perfectly explained variance.
    :returns: A :class:`SARCoherence`.
    :raises TypeError: If ``result`` is neither a clustering result nor a
        sequence of ints, or ``activity`` is not a sequence of floats.
    :raises ValueError: If ``activity`` is empty, differs in length from the
        labeling, or holds a number too large to convert to a double; if a
        label does not fit the native 32-bit signed label type,
        or a clustering result carries a member index that does not fit a
        native ``size_t``, a negative one being the reachable case; or if
        ``noise`` is not one of the three accepted strings.
    :raises RuntimeError: If an activity value is infinite, in which case the
        message names the offending index; or if the magnitudes are large enough
        that a mean or a sum of squares overflows, in which case it names the
        quantity that overflowed rather than an index, no single sample being
        responsible. These are refusals raised in C++, and SWIG maps every
        native exception to ``RuntimeError``.

    Example::

        result = oecluster.butina(dm, threshold=0.35)
        coherence = oecluster.sar_coherence(result, activity)
        print(coherence.eta_squared, coherence.omega_squared)
        print(coherence.clusters[0].mean_activity)
    """
    values = _activity_values(activity, "activity")
    # Ahead of the overload split, because emptiness is the one refusal both
    # branches would otherwise disagree about. A length check alone lets
    # ClusteringResult([], []) with an empty activity through -- 0 == 0 -- and
    # the native refusal that catches it downstream arrives as RuntimeError,
    # contradicting the ValueError this function documents.
    if len(values) == 0:
        raise ValueError("sar_coherence() requires a non-empty activity")

    options = _oecluster.SARCoherenceOptions()
    options.noise_handling = _noise_handling(noise)

    # A ClusteringResult takes the native overload that reads the result;
    # anything else goes through the shared label coercion to the overload that
    # reads a label vector. The split exists so a caller holding a result does
    # not have to take it apart, and one holding bare labels does not have to
    # build a result around them.
    if isinstance(result, ClusteringResult):
        if len(values) != result.num_samples:
            raise ValueError(
                f"activity has {len(values)} entries but the clustering has "
                f"{result.num_samples} samples")
        # ClusteringResult stores its labels and its member indices as intp
        # arrays without range validation, so a value too wide for either native
        # vector -- an IntVector of labels, a SizeTVector of members -- reaches
        # _native_clustering_result and surfaces as OverflowError. The bare-label
        # branch already reports the label case as a ValueError, and the two
        # overloads must not disagree about which exception a caller catches.
        # The message names both fields because OverflowError does not say which
        # vector rejected the value, and a negative member index is the easier
        # of the two to hit.
        try:
            native_result = _native_clustering_result(result)
        except OverflowError as error:
            raise ValueError(
                "result must contain labels that fit a 32-bit signed int and "
                "member indices that fit a native size_t") from error
        native = _oecluster.sar_coherence(native_result, values, options)
    else:
        labels = _agreement_labels(result, "result")
        # No separate emptiness refusal here. The activity is already known to
        # be non-empty, so an empty labeling is always a length mismatch, and
        # the comparison below reports it with the same message the
        # ClusteringResult branch gives for the same input shape.
        if len(values) != len(labels):
            raise ValueError(
                f"activity has {len(values)} entries but the clustering has "
                f"{len(labels)} samples")
        native = _oecluster.sar_coherence(labels, values, options)
    # Bound to a local, never scored inline: the scorecard reads the result's
    # member vector, which the result owns and frees with itself.
    return SARCoherence(native)


def activity_landscape(distance_matrix, activity, *, distance_threshold=0.30,
                       activity_threshold=1.0, rmodi_delta=0.625,
                       num_threads=0):
    """
    Measure the activity cliffs in a precomputed distance matrix.

    Answers a different question from :func:`sar_coherence`: not whether a
    clustering groups molecules that behave alike, but whether the descriptor
    itself puts similar activities near one another. A high cliff density says
    small structural changes swing the activity, which is what makes a series
    hard to model and interesting to a chemist.

    :param distance_matrix: Complete SymmetricDistanceMatrix. SparseStorage is
        refused: a nearest neighbour read off a partial matrix is not one.
    :param activity: One measurement per sample; NaN marks a missing one.
    :param distance_threshold: Pairs at or below this are structurally near.
        Shares :func:`cluster_report`'s boundary default of 0.30, because it
        encodes the same judgement about fingerprint distance.
    :param activity_threshold: Activity differences at or above this are
        sharp. One log unit by default.
    :param rmodi_delta: Half-width of the RMODI activity band, in standard
        deviations. 0.625 is the published value.
    :param num_threads: 0 selects the hardware concurrency. The result does not
        depend on this value, bit for bit. Truncated toward zero, so 1.9 selects
        one thread and -0.5 truncates to 0 and therefore selects the hardware
        concurrency, like every other value that truncates there.
        ``OverflowError`` has two sources, with different causes and different
        moments: an infinite value has no integer to truncate to and fails
        inside ``int()``, before the non-negative check runs at all, while a
        finite but oversized value coerces cleanly there and fails later, in the
        binding layer's ``size_t`` assignment.
    :returns: An :class:`ActivityLandscape`.
    :raises TypeError: If ``distance_matrix`` is not a SymmetricDistanceMatrix,
        ``activity`` is not a sequence of floats, or ``num_threads`` is a value
        ``int()`` cannot accept at all, such as None or a complex. That
        coercion runs before the non-negative check below.
    :raises ValueError: If the matrix uses sparse storage, the activity is
        empty, the activity holds a number too large to convert to a double,
        the activity and the matrix cover different numbers of samples,
        any of the three thresholds is non-finite or negative,
        ``num_threads`` is a str ``int()`` cannot parse or a NaN, which the
        same coercion refuses,
        ``num_threads`` truncates toward zero to a negative integer,
        or the gate refuses the matrix --
        similarity-valued, a non-zero self-distance, a non-finite entry, or
        scored on a per-pair feature subset.
    :raises RuntimeError: If an activity value is infinite, or a stored distance
        is negative; or if the activity magnitudes are large enough that
        ``activity_stddev`` or a SALI accumulator overflows to infinity, in
        which case the message names the quantity that overflowed rather than an
        index. These are refusals raised in C++, and SWIG maps every native
        exception to ``RuntimeError``. A negative distance reaches C++ because
        the gate measures finiteness, not sign.
    :raises OverflowError: If ``num_threads`` is an infinity or coerces to an
        integer too large for a ``size_t``.

    Example::

        landscape = oecluster.activity_landscape(
            dm, activity, distance_threshold=0.30, activity_threshold=1.0)
        print(landscape.num_cliffs, landscape.cliff_density)
        print(landscape.max_sali, landscape.rmodi)
    """
    if not isinstance(distance_matrix, SymmetricDistanceMatrix):
        raise TypeError(
            "activity_landscape() expects a SymmetricDistanceMatrix")

    # ValueError, not TypeError: the argument's type is right, its storage is
    # not. Ahead of the gate, whose remedies cannot rescue a sparse matrix.
    if isinstance(distance_matrix.storage, SparseStorage):
        raise ValueError(  # noqa: TRY004
            "activity_landscape requires complete pairwise distances; "
            "SparseStorage is not supported")

    values = _activity_values(activity, "activity")
    # Ahead of the length check for the same reason as in sar_coherence(): a
    # zero-sample matrix is constructible, so 0 == 0 agrees and the emptiness
    # would only be caught in C++, arriving as RuntimeError.
    if len(values) == 0:
        raise ValueError("activity_landscape() requires a non-empty activity")

    if len(values) != distance_matrix.num_samples:
        raise ValueError(
            f"activity has {len(values)} entries but the matrix covers "
            f"{distance_matrix.num_samples} samples")

    # Mirrors validate_landscape_options() in src/clustering/SARCoherence.cpp:
    # the same three thresholds in the same order, each tested for finiteness
    # before sign. Refused here rather than left to that function because SWIG
    # maps every native exception to RuntimeError, and a threshold Python can
    # inspect for itself belongs in the ValueError this signature documents.
    # isfinite rather than isnan: an infinite distance_threshold would call
    # every pair structurally near instead of failing.
    checked = []
    for name, value in (("distance_threshold", distance_threshold),
                        ("activity_threshold", activity_threshold),
                        ("rmodi_delta", rmodi_delta)):
        try:
            coerced = float(value)
        except OverflowError as error:
            # An int beyond double range is the same condition isfinite() is
            # there to refuse; it just fails one step earlier, in the cast.
            raise ValueError(f"{name} must be finite") from error
        if not math.isfinite(coerced):
            raise ValueError(f"{name} must be finite")
        if coerced < 0.0:
            raise ValueError(f"{name} must be non-negative")
        checked.append(coerced)
    distance_value, activity_value, rmodi_value = checked

    num_threads_int = int(num_threads)
    # num_threads reaches a size_t option field, where a negative value raises
    # OverflowError below the gate. Zero stays legal: it means "choose for me".
    if num_threads_int < 0:
        raise ValueError("num_threads must be non-negative")

    _gate.require_comparable(distance_matrix, "activity_landscape")

    options = _oecluster.ActivityLandscapeOptions()
    options.distance_threshold = distance_value
    options.activity_threshold = activity_value
    options.rmodi_delta = rmodi_value
    options.num_threads = num_threads_int
    # Bound rather than scored inline, on the same terms as sar_coherence().
    # This result carries no member vector today; keeping the three entry
    # points identical means adding one later cannot quietly reintroduce a read
    # off a freed parent.
    native = _oecluster.activity_landscape(
        distance_matrix.storage, values, options)
    return ActivityLandscape(native)


def modelability(distance_matrix, activity_classes, *, num_threads=0):
    """
    Score how well a descriptor separates activity classes.

    The MODI index of Golbraikh et al. (2014), generalized to K classes: the
    mean over classes of the fraction of members whose nearest neighbour shares
    the class. Run it before fitting a classifier, to find out whether the
    descriptor carries the signal at all.

    :param distance_matrix: Complete SymmetricDistanceMatrix. SparseStorage is
        refused, on the same grounds as :func:`activity_landscape`.
    :param activity_classes: One class string per sample. An empty string is a
        missing annotation rather than a category, and its sample is dropped.
    :param num_threads: 0 selects the hardware concurrency. Nearest-neighbour
        ties resolve to the lowest scored index, so the result does not depend
        on this value. Truncated toward zero on :func:`activity_landscape`'s
        terms: 1.9 selects one thread, -0.5 truncates to 0 and therefore
        selects the hardware concurrency, and ``OverflowError`` arrives from
        the same two places -- an infinite value out of ``int()``, a finite but
        oversized one out of the binding layer's ``size_t`` assignment.
    :returns: A :class:`Modelability`.
    :raises TypeError: If ``distance_matrix`` is not a SymmetricDistanceMatrix,
        ``activity_classes`` is not a sequence of strings, or ``num_threads``
        is a value ``int()`` cannot accept at all, on
        :func:`activity_landscape`'s terms.
    :raises ValueError: On the same conditions as :func:`activity_landscape`,
        reading ``activity_classes`` for ``activity``, less the three threshold
        refusals and the oversized-value refusal: this function has no
        thresholds, and a class string has no magnitude to overflow. Every
        ``num_threads`` condition carries over unchanged, including the ones
        the coercion raises: a str ``int()`` cannot parse and a NaN as the
        ValueError this clause describes, an infinite value as the
        ``OverflowError`` the parameter above sets out.
    :raises RuntimeError: If a stored distance is negative, on the same terms
        as :func:`activity_landscape`.
    :raises OverflowError: If ``num_threads`` is an infinity or coerces to an
        integer too large for a ``size_t``.

    Example::

        # nan >= 6.0 is False, so without the nan arm every missing
        # measurement would be annotated "inactive" rather than dropped.
        classes = ["" if math.isnan(a) else "active" if a >= 6.0 else "inactive"
                   for a in activity]
        report = oecluster.modelability(dm, classes)
        print(report.modi, report.num_classes)
        print(report.classes[0].label, report.classes[0].fraction_same_class)
    """
    if not isinstance(distance_matrix, SymmetricDistanceMatrix):
        raise TypeError("modelability() expects a SymmetricDistanceMatrix")

    if isinstance(distance_matrix.storage, SparseStorage):
        raise ValueError(  # noqa: TRY004
            "modelability requires complete pairwise distances; "
            "SparseStorage is not supported")

    classes = _activity_classes(activity_classes, "activity_classes")
    # See activity_landscape(): the length check alone lets a zero-sample
    # matrix with an empty annotation through.
    if len(classes) == 0:
        raise ValueError(
            "modelability() requires a non-empty activity_classes")

    if len(classes) != distance_matrix.num_samples:
        raise ValueError(
            f"activity_classes has {len(classes)} entries but the matrix "
            f"covers {distance_matrix.num_samples} samples")

    num_threads_int = int(num_threads)
    if num_threads_int < 0:
        raise ValueError("num_threads must be non-negative")

    _gate.require_comparable(distance_matrix, "modelability")

    options = _oecluster.ModelabilityOptions()
    options.num_threads = num_threads_int
    # Bound rather than scored inline: the scorecard copies this result's
    # per-class rows, which the result owns and frees with itself.
    native = _oecluster.modelability(
        distance_matrix.storage, classes, options)
    return Modelability(native)


# Stands for the default seed=0 while letting maxmin_select() tell an explicit
# seed apart from none: an explicit seed together with initial is a conflict to
# report, not one to resolve silently in either argument's favour.
_SEED_UNSET = object()

_MAXMIN_STOP_NAMES = {
    _oecluster.MaxMinStop_Count: "count",
    _oecluster.MaxMinStop_Threshold: "threshold",
    _oecluster.MaxMinStop_Exhausted: "exhausted",
}

_CIRCLES_METHODS = {
    "maxmin": _oecluster.CirclesMethod_MaxMin,
    "sequential": _oecluster.CirclesMethod_Sequential,
}

_DIVERSITY_KERNELS = {
    "complement": _oecluster.DiversityKernel_Complement,
    "laplacian": _oecluster.DiversityKernel_Laplacian,
}


class MaxMinSelection:
    """A farthest-first selection, returned by :func:`maxmin_select`.

    :ivar indices: Selected items in selection order, initial entries first, as
        positions in the caller's input.
    :ivar pick_distances: Each pick's distance to the earlier selection when it
        was picked; NaN for the seed and for every initial entry. Never
        increases after them.
    :ivar stop: Why the selection stopped: ``"count"``, ``"threshold"`` or
        ``"exhausted"``.
    :ivar excluded: ``[original_index, reason]`` pairs for the items
        normalization dropped; empty outside the named-comparison path.
    """

    # Field order mirrors the :ivar: list above.
    __slots__ = ("indices", "pick_distances", "stop", "excluded")  # noqa: RUF023

    def __init__(self, indices, pick_distances, stop, excluded):
        """
        Construct a selection from values already copied out of native memory.

        :param indices: Selected caller positions.
        :param pick_distances: One distance per index.
        :param stop: Stop reason name.
        :param excluded: ``[original_index, reason]`` pairs.
        """
        self.indices = indices
        self.pick_distances = pick_distances
        self.stop = stop
        self.excluded = excluded

    def __repr__(self):
        return (f"MaxMinSelection(indices={self.indices!r}, "
                f"stop={self.stop!r}, excluded={len(self.excluded)})")


class CirclesResult:
    """A #Circles packing, returned by :func:`circles`.

    :ivar count: Number of members; the #Circles value. A lower bound on the
        packing number under either method.
    :ivar members: Members as positions in the caller's input, in pick order
        (``"maxmin"``) or input order (``"sequential"``). Pairwise strictly
        farther apart than the threshold.
    :ivar threshold: The distance threshold the packing was built at.
    :ivar method: ``"maxmin"`` or ``"sequential"``.
    :ivar excluded: ``[original_index, reason]`` pairs for the items
        normalization dropped; empty outside the named-comparison path.
    """

    # Field order mirrors the :ivar: list above.
    __slots__ = ("count", "members", "threshold", "method", "excluded")  # noqa: RUF023

    def __init__(self, count, members, threshold, method, excluded):
        """
        Construct a packing from values already copied out of native memory.

        :param count: Number of members.
        :param members: Member caller positions.
        :param threshold: The distance threshold.
        :param method: Method name.
        :param excluded: ``[original_index, reason]`` pairs.
        """
        self.count = count
        self.members = members
        self.threshold = threshold
        self.method = method
        self.excluded = excluded

    def __repr__(self):
        return (f"CirclesResult(count={self.count}, "
                f"threshold={self.threshold!r}, method={self.method!r}, "
                f"excluded={len(self.excluded)})")


class VendiResult:
    """A Vendi score, returned by :func:`vendi_score`.

    :ivar score: The Vendi score; between 1 and ``size`` for a PSD kernel.
    :ivar order: 1 or 2.
    :ivar size: Number of items scored, after normalization.
    :ivar kernel: ``"complement"`` or ``"laplacian"``.
    :ivar min_eigenvalue: The smallest kernel eigenvalue, on the scale of K;
        None for order 2.
    :ivar negative_mass: Sum of ``|lambda| / size`` over the eigenvalues order 1
        dropped as negative; None for order 2.
    :ivar excluded: ``[original_index, reason]`` pairs for the items
        normalization dropped; empty outside the named-comparison path.
    """

    # Field order mirrors the :ivar: list above.
    __slots__ = ("score", "order", "size", "kernel", "min_eigenvalue",  # noqa: RUF023
                 "negative_mass", "excluded")

    def __init__(self, score, order, size, kernel, min_eigenvalue,
                 negative_mass, excluded):
        """
        Construct a result from values already copied out of native memory.

        :param score: The score.
        :param order: The order.
        :param size: Items scored.
        :param kernel: Kernel name.
        :param min_eigenvalue: Smallest eigenvalue, or None.
        :param negative_mass: Dropped negative mass, or None.
        :param excluded: ``[original_index, reason]`` pairs.
        """
        self.score = score
        self.order = order
        self.size = size
        self.kernel = kernel
        self.min_eigenvalue = min_eigenvalue
        self.negative_mass = negative_mass
        self.excluded = excluded

    def __repr__(self):
        return (f"VendiResult(score={self.score!r}, order={self.order}, "
                f"size={self.size}, kernel={self.kernel!r}, "
                f"excluded={len(self.excluded)})")


class LogDetResult:
    """A log-determinant diversity, returned by :func:`logdet_diversity`.

    :ivar score: ``log det(K + ridge I)``, or ``-inf`` when the ridged kernel
        is not numerically positive definite.
    :ivar ridge: The ridge added to every eigenvalue.
    :ivar size: Number of items scored, after normalization.
    :ivar kernel: ``"complement"`` or ``"laplacian"``.
    :ivar min_eigenvalue: The smallest kernel eigenvalue, before the ridge.
    :ivar nonpositive_count: Number of ridged eigenvalues at or below the
        tolerance; nonzero exactly when ``score`` is ``-inf``.
    :ivar excluded: ``[original_index, reason]`` pairs for the items
        normalization dropped; empty outside the named-comparison path.
    """

    # Field order mirrors the :ivar: list above.
    __slots__ = ("score", "ridge", "size", "kernel", "min_eigenvalue",  # noqa: RUF023
                 "nonpositive_count", "excluded")

    def __init__(self, score, ridge, size, kernel, min_eigenvalue,
                 nonpositive_count, excluded):
        """
        Construct a result from values already copied out of native memory.

        :param score: The score.
        :param ridge: The ridge.
        :param size: Items scored.
        :param kernel: Kernel name.
        :param min_eigenvalue: Smallest eigenvalue before the ridge.
        :param nonpositive_count: Ridged eigenvalues at or below tolerance.
        :param excluded: ``[original_index, reason]`` pairs.
        """
        self.score = score
        self.ridge = ridge
        self.size = size
        self.kernel = kernel
        self.min_eigenvalue = min_eigenvalue
        self.nonpositive_count = nonpositive_count
        self.excluded = excluded

    def __repr__(self):
        return (f"LogDetResult(score={self.score!r}, ridge={self.ridge!r}, "
                f"size={self.size}, "
                f"nonpositive_count={self.nonpositive_count}, "
                f"excluded={len(self.excluded)})")


class _DiversitySource:
    """What one diversity call runs on, after its input has been dispatched.

    ``positions`` maps a native index to a caller position and is None when
    the two coincide; ``num_positions`` is the caller's item count.
    """

    # A plain class rather than a NamedTuple: in this module, a NamedTuple
    # with methods drives mypy into an internal error.
    def __init__(self, target, size, num_positions, positions, excluded,
                 matrix):
        self.target = target
        self.size = size
        self.num_positions = num_positions
        self.positions = positions
        self.excluded = excluded
        self.matrix = matrix

    def caller_position(self, index):
        """Map a native index back to the caller's position."""
        return index if self.positions is None else self.positions[index]

    def native_index(self, position, what):
        """
        Map a caller position to a native index.

        :param position: Non-negative caller position.
        :param what: Argument name for the messages.
        :returns: The native index.
        :raises ValueError: If the position is out of range, or names an item
            normalization dropped.
        """
        if position >= self.num_positions:
            raise ValueError(
                f"{what} {position} is outside the item range "
                f"(0 to {self.num_positions - 1})")
        if self.positions is None:
            return position
        for index, kept in enumerate(self.positions):
            if kept == position:
                return index
        reason = next(r for p, r in self.excluded if p == position)
        raise ValueError(
            f"{what} {position} names an item that normalization dropped "
            f"({reason})")


def _refuse_comparison_facts(comparison_obj, caller):
    """
    Refuse a comparison whose declared facts rule out ranking its distances.

    Mirrors validate_comparison_facts in src/clustering/DiversityValidation.h,
    ahead of it, because SWIG turns the native ComparisonError into
    RuntimeError. Every fact read here is declared before any pair is scored,
    so a count-limited call cannot pass merely by never reaching a bad pair.
    "unknown" is accepted throughout.

    :param comparison_obj: Native comparison about to be run.
    :param caller: Entry point name for the messages.
    :raises ValueError: If a fact refuses.
    """
    facts = _gate.facts_from_comparison(comparison_obj)
    if facts['is_distance'] is False:
        raise ValueError(
            f"{caller} requires distances, but the comparison reports "
            "similarities; build it with similarity=False")
    if facts['zero_self'] is False:
        raise ValueError(
            f"{caller} requires a zero self-distance, but the comparison "
            "reports that d(x, x) is not zero")
    if facts['data_integrity'] == "nan_present":
        raise ValueError(
            f"{caller} cannot rank distances the comparison declares may be "
            "non-finite (missing='propagate'); use "
            "missing='complete_case'")
    if facts['data_integrity'] == "subset_scored":
        raise ValueError(
            f"{caller} cannot rank distances scored on per-pair feature "
            "subsets (missing='ignore'); they are not mutually comparable. "
            "Use missing='complete_case'")


def _diversity_source(items, comparison, kwargs, caller, *,
                      allow_sparse=False, defer_comparable_check=False):
    """
    Dispatch a diversity entry point's input onto one of its three paths.

    :param items: A SymmetricDistanceMatrix, a native PairwiseComparison, or a
        sequence of items for a named comparison.
    :param comparison: Comparison name; the named path only.
    :param kwargs: Comparison options; the named path only. Consumed.
    :param caller: Entry point name for the messages.
    :param allow_sparse: Accept SparseStorage; for callers whose native entry
        points read sparse entries directly.
    :param defer_comparable_check: Skip the up-front comparison fact check;
        for callers that need to inspect source.size before gating.
    :returns: A :class:`_DiversitySource`.
    :raises TypeError: If the input and the comparison arguments do not fit
        one path.
    :raises ValueError: If the storage is sparse and allow_sparse is False,
        normalization expanded or emptied the item list, or the comparison's
        facts refuse it (unless defer_comparable_check defers the check).
    """
    if isinstance(items, CrossDistanceMatrix):
        raise TypeError(
            f"{caller}() requires a SymmetricDistanceMatrix; a "
            "CrossDistanceMatrix is rectangular and has no storage backend")

    if isinstance(items, (SymmetricDistanceMatrix,
                          _oecluster.PairwiseComparison)):
        if comparison is not None or kwargs:
            raise TypeError(
                f"{caller}() takes no comparison or comparison options with a "
                "distance matrix or a prebuilt comparison, which already fix "
                "the distances")
        if isinstance(items, SymmetricDistanceMatrix):
            # ValueError, not TypeError: the argument's type is right, its
            # storage is not.
            if isinstance(items.storage, SparseStorage) and not allow_sparse:
                raise ValueError(
                    f"{caller} requires complete pairwise distances; "
                    "SparseStorage is not supported")
            size = items.num_samples
            return _DiversitySource(items.storage, size, size, None, [], items)
        if not defer_comparable_check:
            _refuse_comparison_facts(items, caller)
        size = items.Size()
        return _DiversitySource(items, size, size, None, [], None)

    if not isinstance(comparison, str):
        raise TypeError(
            f"{caller}() requires comparison= to name a comparison, such as "
            "'fingerprint', when items is a sequence of items")

    items = list(items)
    _comparisons.validate_request(comparison, False, kwargs)
    kept, excluded = _comparisons.normalize_items(comparison, items, kwargs)
    # Normalization may drop items but must not add them: conformer expansion
    # would map one caller position onto several selectable items, and no
    # index this function returns could then name what was selected.
    if len(kept) + len(excluded) != len(items):
        raise ValueError(
            f"normalizing the inputs for the {comparison!r} comparison turned "
            f"{len(items)} items into {len(kept)}; {caller} reports caller "
            "positions and cannot map expanded items back to them. Expand "
            "conformers up front, or pass expand_conformers=False")
    if not kept:
        raise ValueError(f"{caller}() requires at least one item")
    dropped = {index for index, _ in excluded}
    positions = [p for p in range(len(items)) if p not in dropped]
    comparison_obj, _, _ = _comparisons.build_comparison(
        kept, comparison, False, kwargs, symmetric=True)
    _refuse_comparison_facts(comparison_obj, caller)
    return _DiversitySource(comparison_obj, comparison_obj.Size(), len(items),
                            positions, [list(entry) for entry in excluded],
                            None)


def _diversity_int(value, name, minimum):
    """
    Coerce an integer option, refusing a non-integer with ValueError.

    :param value: Caller value.
    :param name: Option name for the messages.
    :param minimum: Smallest accepted value.
    :returns: The coerced int.
    :raises ValueError: If the value is not an integer, is below the minimum,
        or does not fit a size_t.
    """
    # operator.index() accepts int and numpy integers and refuses
    # float/str/None; int() would silently truncate 2.5. The design's error
    # table makes a non-integer a ValueError, so the TypeError is translated.
    # bool is an int subclass that index() also accepts, and count=True is a
    # caller mistake rather than a count of one, so it is refused first.
    if isinstance(value, bool):
        raise ValueError(f"{name} must be an integer, got {value!r}")  # noqa: TRY004
    try:
        coerced = operator.index(value)
    except TypeError:
        raise ValueError(f"{name} must be an integer, got {value!r}") from None
    if coerced < minimum:
        raise ValueError(f"{name} must be at least {minimum}")
    if coerced > _SIZE_T_MAX:
        raise ValueError(f"{name} exceeds size_t maximum")
    return coerced


def _diversity_threshold(value):
    """
    Coerce a distance threshold, refusing NaN, infinity and negatives.

    :param value: Caller value.
    :returns: The coerced float.
    :raises ValueError: If the value is NaN, infinite or negative.
    """
    try:
        coerced = float(value)
    except OverflowError as error:
        # An int beyond double range is the infinity refused below, failing
        # one step earlier.
        raise ValueError("threshold must be finite") from error
    if math.isnan(coerced):
        raise ValueError("threshold must be a number, not NaN")
    if math.isinf(coerced):
        raise ValueError("threshold must be finite")
    if coerced < 0.0:
        raise ValueError("threshold must be non-negative")
    return coerced


def _maxmin_seed(seed):
    """
    Resolve maxmin_select()'s seed argument.

    :param seed: ``_SEED_UNSET``, a non-negative int, ``"medoid"`` or
        ``"farthest"``.
    :returns: ``(seed_mode, position)``; position is None unless the mode is
        Index with an explicit seed.
    :raises ValueError: If the seed is none of those.
    """
    if seed is _SEED_UNSET:
        return _oecluster.MaxMinSeed_Index, None
    if isinstance(seed, str):
        modes = {"medoid": _oecluster.MaxMinSeed_Medoid,
                 "farthest": _oecluster.MaxMinSeed_Farthest}
        mode = modes.get(seed.lower())
        if mode is None:
            raise ValueError(
                f"Unknown seed: {seed!r}; expected a non-negative int, "
                "'medoid' or 'farthest'")
        return mode, None
    return _oecluster.MaxMinSeed_Index, _diversity_int(seed, "seed", 0)


def maxmin_select(items, *, count=None, threshold=None, seed=_SEED_UNSET,
                  initial=None, comparison=None, similarity=False,
                  num_threads=0, chunk_size=256,
                  **kwargs) -> "MaxMinSelection":
    """
    Select a diverse subset by farthest-first (MaxMin) picking.

    Each pick is the unselected item farthest from everything already
    selected, ties to the smaller index, so the result depends only on the
    distances and the arguments. With a threshold, a candidate at or within
    that distance of the selection is not added, and every pick after the
    seed and ``initial`` is strictly farther than the threshold from all
    earlier picks.

    ``items`` selects the path: a :class:`SymmetricDistanceMatrix` reads a
    precomputed matrix; a native comparison object (for example the one
    :class:`FingerprintComparison` returns) is evaluated lazily; anything else
    is a sequence of items and ``comparison`` names how to compare them, as in
    :func:`pdist`, again evaluated lazily. The lazy paths compare only the
    pairs the selection needs, O(N·k) for k picks, and never build the O(N^2)
    matrix.

    Every position argument and result -- ``seed``, ``initial``, ``indices``,
    ``excluded`` -- refers to the caller's ``items``, even when normalization
    dropped some of them. Left unset, ``seed`` starts from the first item that
    survived normalization.

    :param items: Matrix, prebuilt comparison, or sequence of items.
    :param count: Total selection size, ``initial`` included; a positive int.
    :param threshold: Stop before a candidate at or within this distance.
        At least one of ``count`` and ``threshold`` is required.
    :param seed: Starting position (default 0), ``"medoid"`` for the item
        with the smallest distance sum (matrix path only), or
        ``"farthest"`` for the item farthest from the first item that
        survived normalization.
    :param initial: An existing selection to extend, reported first. Cannot
        be combined with an explicit ``seed``.
    :param comparison: Comparison name, required with a sequence of items.
    :param similarity: Refused when True; selection runs on distances.
    :param num_threads: Worker threads for the lazy paths; 0 selects the
        hardware concurrency. Each worker holds one clone of the comparison,
        which for comparisons that copy their items (MCS) costs O(N) apiece.
    :param chunk_size: Items per work unit, at least one.
    :param kwargs: Comparison options for a named comparison.
    :returns: A :class:`MaxMinSelection`.
    :raises TypeError: If the arguments fit none of the three paths, a
        comparison option is unknown, an ``initial`` entry is not an int, or
        ``initial`` comes with an explicit ``seed``.
    :raises ValueError: On an invalid ``count``, ``threshold``, ``seed``,
        ``initial``, ``num_threads`` or ``chunk_size``; ``similarity=True``;
        sparse storage; an empty input; normalization that expanded the
        input or dropped a named position; ``seed="medoid"`` without a
        matrix; or a matrix or comparison whose distances cannot be ranked.
    :raises RuntimeError: If a distance read during the selection is NaN or
        infinite, or a medoid seed's distance sum overflows.

    Example::

        picked = oecluster.maxmin_select(mols, comparison="fingerprint",
                                         count=50)
        subset = [mols[i] for i in picked.indices]
    """
    if similarity:
        raise ValueError(
            "maxmin_select() selects on distances; similarity=True is not "
            "supported")
    count_value = None if count is None else _diversity_int(count, "count", 1)
    threshold_value = (None if threshold is None
                       else _diversity_threshold(threshold))
    if count_value is None and threshold_value is None:
        raise ValueError(
            "maxmin_select() requires count, threshold, or both")
    num_threads_value = _diversity_int(num_threads, "num_threads", 0)
    chunk_size_value = _diversity_int(chunk_size, "chunk_size", 1)
    seed_mode, seed_position = _maxmin_seed(seed)
    if initial is None:
        initial_positions = []
    else:
        # bool passes operator.index(); it is refused as _diversity_int
        # refuses it, so initial=[True] is not read as position 1.
        # Materialize the iterable once so one-shot iterators work.
        try:
            initial_list = list(initial)
            if any(isinstance(p, bool) for p in initial_list):
                raise TypeError("bool entry")
            initial_positions = [operator.index(p) for p in initial_list]
        except TypeError as error:
            raise TypeError(
                "maxmin_select() initial must be a sequence of ints"
            ) from error
    if initial_positions and seed is not _SEED_UNSET:
        raise TypeError(
            "maxmin_select() takes initial or seed, not both: an initial "
            "selection replaces the seed")
    if any(p < 0 for p in initial_positions):
        raise ValueError("initial entries must be non-negative")
    if len(set(initial_positions)) != len(initial_positions):
        raise ValueError("initial entries must be unique")
    if (count_value is not None and len(initial_positions) > count_value):
        raise ValueError("initial holds more entries than count")

    source = _diversity_source(items, comparison, kwargs, "maxmin_select")
    if source.size == 0:
        raise ValueError("maxmin_select() requires at least one item")
    # Tested on the caller's string rather than seed_mode: mypy reports a
    # comparison against a wrapper enum attribute as a cyclic definition.
    wants_medoid = isinstance(seed, str) and seed.lower() == "medoid"
    if wants_medoid and source.matrix is None:
        raise ValueError(
            "seed='medoid' requires a distance matrix: on a comparison it "
            "would cost every pairwise comparison, which the lazy path exists "
            "to avoid")
    if count_value is not None and count_value > source.size:
        raise ValueError(
            f"count must be at most the item count ({source.size})")
    seed_index = (0 if seed_position is None
                  else source.native_index(seed_position, "seed"))
    initial_indices = [source.native_index(p, "initial entry")
                       for p in initial_positions]
    if source.matrix is not None:
        _gate.require_comparable(source.matrix, "maxmin_select")

    options = _oecluster.MaxMinOptions()
    options.count = 0 if count_value is None else count_value
    options.threshold = (math.nan if threshold_value is None
                         else threshold_value)
    options.seed_mode = seed_mode
    options.seed = seed_index
    if initial_indices:
        native_initial = _oecluster.SizeTVector()
        for index in initial_indices:
            native_initial.push_back(index)
        options.initial = native_initial
    options.num_threads = num_threads_value
    options.chunk_size = chunk_size_value

    native = _oecluster.maxmin_select(source.target, options)
    # Copied out while the native result is alive: its vectors are owned by
    # it and read as empty once it is gone.
    return MaxMinSelection(
        indices=[source.caller_position(i) for i in native.indices],
        pick_distances=[float(d) for d in native.pick_distances],
        stop=_MAXMIN_STOP_NAMES[native.stop],
        excluded=source.excluded,
    )


def circles(items, *, threshold, method="maxmin", comparison=None,
            similarity=False, num_threads=0, chunk_size=256,
            **kwargs) -> "CirclesResult":
    """
    Count the #Circles coverage of a set: a packing at a distance threshold.

    #Circles (Xie et al., ICLR 2023) is the size of a set of items that are
    pairwise strictly farther apart than ``threshold``. The paper's headline
    threshold is a Tanimoto distance of 0.75. Any valid packing is a lower
    bound on the true packing number, so ``count`` is a lower bound under
    either method, and the two methods can disagree.

    ``method="maxmin"`` packs farthest-first from the first item that
    survived normalization: it is :func:`maxmin_select` with this threshold,
    no count, and the default seed.
    ``method="sequential"`` is the paper's reference greedy pass over input
    order, accepting an item when it is farther than ``threshold`` from every
    member so far. The paper's implementation also shuffles and repeats that
    pass in chunks; this one does not.

    ``items`` selects the path exactly as for :func:`maxmin_select`, and
    ``members`` and ``excluded`` refer to the caller's positions.

    :param items: Matrix, prebuilt comparison, or sequence of items.
    :param threshold: Distance threshold; finite and non-negative.
    :param method: ``"maxmin"`` (the default) or ``"sequential"``.
    :param comparison: Comparison name, required with a sequence of items.
    :param similarity: Refused when True; the threshold is a distance.
    :param num_threads: Worker threads for the lazy paths; 0 selects the
        hardware concurrency.
    :param chunk_size: Items per work unit, at least one.
    :param kwargs: Comparison options for a named comparison.
    :returns: A :class:`CirclesResult`.
    :raises TypeError: If the arguments fit none of the three paths, a
        comparison option is unknown, or ``seed`` or ``initial`` is passed.
    :raises ValueError: On an invalid ``threshold``, ``method``,
        ``num_threads`` or ``chunk_size``; ``similarity=True``; sparse
        storage; an empty input; normalization that expanded the input; or a
        matrix or comparison whose distances cannot be ranked.
    :raises RuntimeError: If a distance read during the packing is NaN or
        infinite.

    Example::

        packing = oecluster.circles(mols, comparison="fingerprint",
                                    threshold=0.75)
        print(packing.count)
    """
    if "seed" in kwargs or "initial" in kwargs:
        raise TypeError(
            "circles() takes no seed or initial: the packing always starts "
            "from the first item that survived normalization")
    if similarity:
        raise ValueError(
            "circles() packs on distances; similarity=True is not supported")
    threshold_value = _diversity_threshold(threshold)
    method_key = method.lower() if isinstance(method, str) else None
    if method_key not in _CIRCLES_METHODS:
        raise ValueError(
            f"Unknown circles method: {method!r}; expected 'maxmin' or "
            "'sequential'")
    num_threads_value = _diversity_int(num_threads, "num_threads", 0)
    chunk_size_value = _diversity_int(chunk_size, "chunk_size", 1)

    source = _diversity_source(items, comparison, kwargs, "circles")
    if source.size == 0:
        raise ValueError("circles() requires at least one item")
    if source.matrix is not None:
        _gate.require_comparable(source.matrix, "circles")

    options = _oecluster.CirclesOptions()
    options.method = _CIRCLES_METHODS[method_key]
    options.num_threads = num_threads_value
    options.chunk_size = chunk_size_value

    native = _oecluster.circles(source.target, threshold_value, options)
    return CirclesResult(
        count=int(native.count),
        members=[source.caller_position(i) for i in native.members],
        threshold=float(native.threshold),
        method=method_key,
        excluded=source.excluded,
    )


def _diversity_number(value, name):
    """
    Coerce a real-valued option, refusing a non-number with TypeError.

    :param value: Caller value.
    :param name: Option name for the messages.
    :returns: The coerced float; an int beyond double range becomes infinity,
        which every caller refuses as non-finite.
    :raises TypeError: If the value is a bool, a string, or not a real number.
    """
    if isinstance(value, bool) or not isinstance(
            value, (int, float, np.integer, np.floating)):
        raise TypeError(f"{name} must be a number, got {value!r}")
    try:
        return float(value)
    except OverflowError:
        return math.inf


def _diversity_kernel(kernel, bandwidth):
    """
    Resolve the kernel name and its bandwidth.

    :param kernel: ``"complement"`` or ``"laplacian"``.
    :param bandwidth: None for complement; a positive finite number for
        laplacian.
    :returns: ``(kernel_key, bandwidth_value)``; the bandwidth is NaN (the
        native "unset") for complement.
    :raises ValueError: If the kernel is unknown or the bandwidth does not fit
        it.
    :raises TypeError: If the bandwidth is not a number.
    """
    key = kernel.lower() if isinstance(kernel, str) else None
    if key not in _DIVERSITY_KERNELS:
        raise ValueError(
            f"Unknown kernel: {kernel!r}; expected 'complement' or "
            "'laplacian'")
    if key == "complement":
        if bandwidth is not None:
            raise ValueError("bandwidth applies only to kernel='laplacian'")
        return key, math.nan
    if bandwidth is None:
        raise ValueError("kernel='laplacian' requires a bandwidth")
    value = _diversity_number(bandwidth, "bandwidth")
    if not math.isfinite(value) or value <= 0.0:
        raise ValueError("bandwidth must be positive and finite")
    return key, value


def _condensed_pair(n, index):
    """Map a condensed index back to its (row, column) pair, row < column."""
    row = 0
    while index >= n - 1 - row:
        index -= n - 1 - row
        row += 1
    return row, row + 1 + index


def _require_kernel_distances(matrix, kernel_key, caller):
    """
    Refuse a matrix holding a distance the kernel cannot use.

    Native code would refuse it too, but as a RuntimeError, after reading up
    to the offending pair; checking the whole matrix here makes it a
    ValueError, as the gate's refusal of a non-finite entry is.

    :param matrix: A SymmetricDistanceMatrix that already passed the gate.
    :param kernel_key: ``"complement"`` or ``"laplacian"``.
    :param caller: Entry point name for the messages.
    :raises ValueError: On a complement distance outside [0, 1], or a negative
        Laplacian distance, naming the first such pair.
    """
    condensed = np.asarray(matrix.condensed)
    # A boolean mask and argmax, not flatnonzero: an input that is wrong
    # throughout (distances above 1 under the complement kernel) would
    # otherwise allocate an int64 index for every entry.
    mask = condensed < 0.0
    if kernel_key == "complement":
        mask |= condensed > 1.0
    if not mask.any():
        return
    index = int(np.argmax(mask))
    row, column = _condensed_pair(matrix.num_samples, index)
    value = float(condensed[index])
    if kernel_key == "complement":
        raise ValueError(
            f"{caller}() kernel='complement' requires distances in [0, 1], "
            f"but d({row}, {column}) = {value!r}; use kernel='laplacian' for "
            "other distances")
    raise ValueError(
        f"{caller}() kernel='laplacian' requires non-negative distances, "
        f"but d({row}, {column}) = {value!r}")


def _require_exact_size(size, max_exact, caller):
    """
    Refuse an input above the exact-score ceiling, ahead of native code.

    :param size: Items that survived normalization.
    :param max_exact: The ceiling.
    :param caller: Entry point name for the messages.
    :raises ValueError: If size exceeds max_exact.
    """
    if size > max_exact:
        raise ValueError(
            f"{caller}() computes an exact spectrum of at most "
            f"max_exact={max_exact} items, but the input has {size}; raise "
            "max_exact (memory grows as 8n^2 bytes and time as n^3) or use "
            "vendi_score(order=2), which needs no spectrum")


def vendi_score(items, *, order=1, kernel="complement", bandwidth=None,
                max_exact=2048, comparison=None, similarity=False,
                num_threads=0, chunk_size=256, **kwargs) -> "VendiResult":
    """
    Score a set's diversity as its Vendi score.

    The Vendi score (Friedman and Dieng, TMLR 2023) is the effective number of
    distinct items in a set: 1 when every item is identical, and the item
    count when every pair is maximally dissimilar. It reads the similarity
    kernel K built from the distances by ``kernel``: ``"complement"`` is
    ``1 - d`` and needs distances in [0, 1]; ``"laplacian"`` is
    ``exp(-d / bandwidth)``.

    ``order=1`` is ``exp(-sum p log p)`` over ``p = lambda / n`` for the
    eigenvalues of K above ``n * eps * max|lambda|``. The rest are dropped
    without renormalizing, as the reference implementation does, and
    ``min_eigenvalue`` and ``negative_mass`` report what a non-PSD kernel
    lost. It decomposes an n x n kernel, so it refuses more than
    ``max_exact`` items. ``order=2`` is ``n^2 / ||K||_F^2``: it needs no
    spectrum, has no ceiling, and on a comparison runs in O(N) memory. For a
    non-PSD kernel it includes the negative eigenvalues' squares, and so
    differs from the reference.

    ``items`` selects the path exactly as for :func:`maxmin_select`.

    :param items: Matrix, prebuilt comparison, or sequence of items.
    :param order: 1 or 2.
    :param kernel: ``"complement"`` (the default) or ``"laplacian"``.
    :param bandwidth: Laplacian bandwidth, positive and finite; None for
        complement.
    :param max_exact: Largest input order 1 decomposes. Memory grows as
        ``8 n^2`` bytes and time as ``n^3``.
    :param comparison: Comparison name, required with a sequence of items.
    :param similarity: Refused when True; the kernel is built from distances.
    :param num_threads: Worker threads for the lazy paths; 0 selects the
        hardware concurrency. Each worker holds one clone of the comparison,
        which for comparisons that copy their items (MCS) costs O(N) apiece.
    :param chunk_size: Rows per work unit, at least one.
    :param kwargs: Comparison options for a named comparison.
    :returns: A :class:`VendiResult`.
    :raises TypeError: If the arguments fit none of the three paths, a
        comparison option is unknown, or ``bandwidth`` is not a number.
    :raises ValueError: On an invalid ``order``, ``kernel``, ``bandwidth``,
        ``max_exact``, ``num_threads`` or ``chunk_size``; ``similarity=True``;
        sparse storage; an empty input; normalization that expanded the
        input; more than ``max_exact`` items at order 1; a refused matrix or
        comparison; or a matrix distance the kernel cannot use.
    :raises RuntimeError: If a distance read from a comparison is NaN,
        infinite, or one the kernel cannot use, or the eigenvalue solver does
        not converge.

    Example::

        result = oecluster.vendi_score(mols, comparison="fingerprint")
        print(result.score)
    """
    if similarity:
        raise ValueError(
            "vendi_score() scores distances; similarity=True is not "
            "supported")
    if (isinstance(order, bool)
            or not isinstance(order, (int, np.integer))
            or order not in (1, 2)):
        raise ValueError(f"vendi_score() order must be 1 or 2, got {order!r}")
    order_value = int(order)
    kernel_key, bandwidth_value = _diversity_kernel(kernel, bandwidth)
    max_exact_value = _diversity_int(max_exact, "max_exact", 1)
    num_threads_value = _diversity_int(num_threads, "num_threads", 0)
    chunk_size_value = _diversity_int(chunk_size, "chunk_size", 1)

    source = _diversity_source(items, comparison, kwargs, "vendi_score")
    if source.size == 0:
        raise ValueError("vendi_score() requires at least one item")
    if order_value == 1:
        _require_exact_size(source.size, max_exact_value, "vendi_score")
    if source.matrix is not None:
        _gate.require_comparable(source.matrix, "vendi_score")
        _require_kernel_distances(source.matrix, kernel_key, "vendi_score")

    options = _oecluster.VendiOptions()
    options.order = order_value
    options.kernel = _DIVERSITY_KERNELS[kernel_key]
    options.bandwidth = bandwidth_value
    options.max_exact = max_exact_value
    options.num_threads = num_threads_value
    options.chunk_size = chunk_size_value

    native = _oecluster.vendi_score(source.target, options)
    exact = order_value == 1
    return VendiResult(
        score=float(native.score),
        order=order_value,
        size=int(native.size),
        kernel=kernel_key,
        min_eigenvalue=float(native.min_eigenvalue) if exact else None,
        negative_mass=float(native.negative_mass) if exact else None,
        excluded=source.excluded,
    )


def logdet_diversity(items, *, ridge=0.0, kernel="complement",
                     bandwidth=None, max_exact=2048, comparison=None,
                     similarity=False, num_threads=0, chunk_size=256,
                     **kwargs) -> "LogDetResult":
    """
    Score a set's diversity as the log-determinant of its kernel.

    The score is ``log det(K + ridge I)`` over the similarity kernel K that
    ``kernel`` builds from the distances, as for :func:`vendi_score`. It is a
    positive-definite log-determinant: when any eigenvalue of ``K + ridge I``
    is at or below ``n * eps * max|mu|``, the score is ``-inf``. That covers a
    singular kernel (duplicate items) and an indefinite one, even when an
    even number of negative eigenvalues leaves the determinant positive.
    ``nonpositive_count`` and ``min_eigenvalue`` tell the two apart. A ridge
    makes a singular PSD kernel finite only when it clears that tolerance; for
    n identical items the smallest useful ridge is about ``n * eps * n``.

    It decomposes an n x n kernel, so it refuses more than ``max_exact``
    items.

    :param items: Matrix, prebuilt comparison, or sequence of items.
    :param ridge: Added to every eigenvalue; finite and non-negative.
    :param kernel: ``"complement"`` (the default) or ``"laplacian"``.
    :param bandwidth: Laplacian bandwidth, positive and finite; None for
        complement.
    :param max_exact: Largest input decomposed. Memory grows as ``8 n^2``
        bytes and time as ``n^3``.
    :param comparison: Comparison name, required with a sequence of items.
    :param similarity: Refused when True; the kernel is built from distances.
    :param num_threads: Worker threads for the lazy paths; 0 selects the
        hardware concurrency.
    :param chunk_size: Rows per work unit, at least one.
    :param kwargs: Comparison options for a named comparison.
    :returns: A :class:`LogDetResult`.
    :raises TypeError: If the arguments fit none of the three paths, a
        comparison option is unknown, or ``ridge`` or ``bandwidth`` is not a
        number.
    :raises ValueError: As :func:`vendi_score`, plus a negative or non-finite
        ``ridge``.
    :raises RuntimeError: As :func:`vendi_score`.

    Example::

        result = oecluster.logdet_diversity(mols, comparison="fingerprint",
                                            ridge=1e-6)
        print(result.score, result.nonpositive_count)
    """
    if similarity:
        raise ValueError(
            "logdet_diversity() scores distances; similarity=True is not "
            "supported")
    ridge_value = _diversity_number(ridge, "ridge")
    if not math.isfinite(ridge_value) or ridge_value < 0.0:
        raise ValueError("ridge must be finite and non-negative")
    kernel_key, bandwidth_value = _diversity_kernel(kernel, bandwidth)
    max_exact_value = _diversity_int(max_exact, "max_exact", 1)
    num_threads_value = _diversity_int(num_threads, "num_threads", 0)
    chunk_size_value = _diversity_int(chunk_size, "chunk_size", 1)

    source = _diversity_source(items, comparison, kwargs, "logdet_diversity")
    if source.size == 0:
        raise ValueError("logdet_diversity() requires at least one item")
    _require_exact_size(source.size, max_exact_value, "logdet_diversity")
    if source.matrix is not None:
        _gate.require_comparable(source.matrix, "logdet_diversity")
        _require_kernel_distances(source.matrix, kernel_key,
                                  "logdet_diversity")

    options = _oecluster.LogDetOptions()
    options.ridge = ridge_value
    options.kernel = _DIVERSITY_KERNELS[kernel_key]
    options.bandwidth = bandwidth_value
    options.max_exact = max_exact_value
    options.num_threads = num_threads_value
    options.chunk_size = chunk_size_value

    native = _oecluster.logdet_diversity(source.target, options)
    return LogDetResult(
        score=float(native.score),
        ridge=ridge_value,
        size=int(native.size),
        kernel=kernel_key,
        min_eigenvalue=float(native.min_eigenvalue),
        nonpositive_count=int(native.nonpositive_count),
        excluded=source.excluded,
    )


_SPHERE_ORDERS = {
    "input": _oecluster.SphereOrder_Input,
    "neighbors": _oecluster.SphereOrder_Neighbors,
}

_SPHERE_ASSIGNMENTS = {
    "first": _oecluster.SphereAssignment_First,
    "nearest": _oecluster.SphereAssignment_Nearest,
}

_SPHERE_ORDER_TYPE = (
    "sphere_exclusion() order must be 'input', 'neighbors' or a sequence of "
    "ints")


def _sphere_order(order):
    """
    Resolve sphere_exclusion()'s order argument.

    :param order: ``"input"``, ``"neighbors"``, or an iterable of caller
        positions.
    :returns: ``(native_order, positions)``; positions is None unless the
        order is a sequence.
    :raises ValueError: If a string names no order.
    :raises TypeError: If the order is neither a string nor an iterable of
        ints, or an entry is a bool.
    """
    if isinstance(order, str):
        native = _SPHERE_ORDERS.get(order.lower())
        if native is None:
            raise ValueError(
                f"Unknown sphere_exclusion order: {order!r}; expected "
                "'input', 'neighbors' or a sequence of positions")
        return native, None
    try:
        entries = list(order)
    except TypeError:
        raise TypeError(_SPHERE_ORDER_TYPE) from None
    positions = []
    for entry in entries:
        # bool is an int subclass; order=[0, True] is a mistake, not item 1.
        if isinstance(entry, bool):
            raise TypeError(_SPHERE_ORDER_TYPE)
        try:
            positions.append(operator.index(entry))
        except TypeError:
            raise TypeError(_SPHERE_ORDER_TYPE) from None
    return _oecluster.SphereOrder_Permutation, positions


def _sphere_native_order(positions, source):
    """
    Map a permutation of caller positions to native indices.

    :param positions: Caller positions from :func:`_sphere_order`.
    :param source: The call's :class:`_DiversitySource`.
    :returns: Native indices in the caller's order, dropped positions skipped.
    :raises ValueError: If ``positions`` is not a permutation of every caller
        position, dropped ones included.
    """
    count = source.num_positions
    if len(positions) != count:
        raise ValueError(
            f"sphere_exclusion order has {len(positions)} entries for "
            f"{count} items")
    seen = set()
    for position in positions:
        if position < 0 or position >= count:
            raise ValueError(
                f"sphere_exclusion order entry {position} is outside the "
                f"item range (0 to {count - 1})")
        if position in seen:
            raise ValueError(
                f"sphere_exclusion order repeats position {position}")
        seen.add(position)
    if source.positions is None:
        return positions
    native = {kept: index for index, kept in enumerate(source.positions)}
    return [native[position] for position in positions if position in native]


def sphere_exclusion(items, threshold, *, order="input", reordering=False,
                     assignment="first", comparison=None, similarity=False,
                     num_threads=0, chunk_size=4096,
                     **kwargs) -> "SphereExclusionResult":
    """
    Cluster by sphere exclusion: leader, Butina or DISE, by seed order.

    Centers are taken in ``order``. Each center claims every unclaimed item
    at or within ``threshold`` of it. ``order="input"`` is leader clustering.
    ``order="neighbors"`` is Butina: descending neighbor count, with equal
    counts going to the larger index, and optional ``reordering``. With
    ``assignment="first"`` it equals :func:`butina` on every matrix both
    accept; nearest assignment keeps Butina's centers but may move members.
    A sequence of caller positions is Directed Sphere Exclusion (DISE): for
    example, ``argsort`` of the distances to a reference item. The sequence
    must name every caller position, including positions that normalization
    dropped.

    ``assignment="nearest"`` keeps the centers and moves each other item to
    its nearest center, with ties going to the earlier center.

    ``items`` selects the path exactly as for :func:`maxmin_select`. Labels,
    clusters, centers and ``excluded`` refer to the caller's positions.

    :param items: Matrix, prebuilt comparison, or sequence of items.
    :param threshold: Distance threshold; finite and non-negative.
    :param order: ``"input"`` (the default), ``"neighbors"``, or a sequence
        of caller positions.
    :param reordering: Butina's reordering; ``order="neighbors"`` only.
    :param assignment: ``"first"`` (the default) or ``"nearest"``.
    :param comparison: Comparison name, required with a sequence of items.
    :param similarity: Refused when True; the threshold is a distance.
    :param num_threads: Worker threads for the lazy paths and the neighbor
        graph; 0 selects the hardware concurrency.
    :param chunk_size: Items or pairs per work unit, at least one.
    :param kwargs: Comparison options for a named comparison.
    :returns: A :class:`SphereExclusionResult`.
    :raises TypeError: If the arguments fit none of the three paths, a
        comparison option is unknown, ``order`` is not a string or a sequence
        of ints, or ``reordering`` is a string.
    :raises ValueError: On an invalid ``threshold``, ``order``,
        ``assignment``, ``num_threads`` or ``chunk_size``; ``reordering``
        without the neighbor order; ``similarity=True``; sparse storage; an
        empty sequence; the neighbor order on a lazy path; or a matrix or
        comparison whose distances cannot be ranked.
    :raises RuntimeError: If a comparison returns a NaN or infinite distance.

    Example::

        reference = distances_to_reference  # one distance per item
        result = oecluster.sphere_exclusion(
            mols, 0.6, comparison="fingerprint",
            order=numpy.argsort(reference, kind="stable"))
        print(result.centers)
    """
    if similarity:
        raise ValueError(
            "sphere_exclusion() clusters on distances; similarity=True is not "
            "supported")
    threshold_value = _diversity_threshold(threshold)
    order_native, order_positions = _sphere_order(order)
    reordering_value = _flag(reordering, "reordering")
    if reordering_value and order_native != _oecluster.SphereOrder_Neighbors:
        raise ValueError(
            "sphere_exclusion() reordering requires order='neighbors'")
    assignment_key = (assignment.lower() if isinstance(assignment, str)
                      else None)
    if assignment_key not in _SPHERE_ASSIGNMENTS:
        raise ValueError(
            f"Unknown sphere_exclusion assignment: {assignment!r}; expected "
            "'first' or 'nearest'")
    num_threads_value = _diversity_int(num_threads, "num_threads", 0)
    chunk_size_value = _diversity_int(chunk_size, "chunk_size", 1)

    source = _diversity_source(items, comparison, kwargs, "sphere_exclusion")
    if (order_native == _oecluster.SphereOrder_Neighbors
            and source.matrix is None):
        raise ValueError(
            "sphere_exclusion() with order='neighbors' needs every pairwise "
            "distance; pass a precomputed distance matrix from pdist()")
    native_order = (None if order_positions is None
                    else _sphere_native_order(order_positions, source))
    if source.matrix is not None:
        _gate.require_comparable(source.matrix, "sphere_exclusion")

    options = _oecluster.SphereExclusionOptions()
    options.distance_threshold = threshold_value
    options.order = order_native
    if native_order is not None:
        permutation = _oecluster.SizeTVector()
        for index in native_order:
            permutation.push_back(index)
        options.permutation = permutation
    options.reordering = reordering_value
    options.assignment = _SPHERE_ASSIGNMENTS[assignment_key]
    options.num_threads = num_threads_value
    options.chunk_size = chunk_size_value

    native = _oecluster.sphere_exclusion(source.target, options)
    labels = [-1] * source.num_positions
    for index, label in enumerate(native.Labels()):
        labels[source.caller_position(index)] = int(label)
    return SphereExclusionResult(
        labels,
        [[source.caller_position(i) for i in cluster]
         for cluster in native.Members()],
        centers=[source.caller_position(i) for i in native.Centers()],
        excluded=source.excluded,
    )


_INTP_MAX = int(np.iinfo(np.intp).max)


def _knn_check_k(k, size, caller):
    """
    Refuse a neighbor count that does not fit the item count.

    Mirrors validate_knn_k in src/clustering/KNNGraphBuild.h ahead of it,
    because SWIG turns the native invalid_argument into RuntimeError.

    :param k: Coerced neighbor count.
    :param size: Native item count.
    :param caller: Entry point name for the messages.
    :raises ValueError: If there is one item, or k is outside 1..size - 1.
    """
    if size == 0:
        return
    if size == 1:
        raise ValueError(
            f"{caller}() needs at least two items: a single item has no "
            "neighbors")
    if not 1 <= k <= size - 1:
        raise ValueError(
            f"{caller}() k must be between 1 and {size - 1} for {size} items, "
            f"got {k}")


def knn_graph(items, k, *, comparison=None, similarity=False, num_threads=0,
              chunk_size=4096, **kwargs) -> "KNNGraph":
    """
    Find each item's ``k`` nearest other items.

    Each row lists ``k`` items other than its own, ordered by ascending
    (distance, position); equal distances go to the lower position. The
    result does not depend on ``num_threads`` or ``chunk_size``.

    ``items`` selects the path as for :func:`sphere_exclusion`, and sparse
    storage is also accepted. A sparse matrix must hold every pair at or
    within its cutoff, as :func:`pdist` with a cutoff produces, and every item
    needs at least ``k`` stored neighbors; raise the cutoff or lower ``k``
    otherwise. The lazy paths evaluate every ordered pair, N(N-1)
    comparisons, and keep only the graph.

    :param items: Matrix, prebuilt comparison, or sequence of items.
    :param k: Neighbors per item, 1 to n - 1.
    :param comparison: Comparison name, required with a sequence of items.
    :param similarity: Refused when True; the graph ranks distances.
    :param num_threads: Worker threads; 0 selects the hardware concurrency.
    :param chunk_size: Pairwise distances per work unit, at least one.
    :param kwargs: Comparison options for a named comparison.
    :returns: A :class:`KNNGraph`.
    :raises TypeError: If the arguments fit none of the three paths, or a
        comparison option is unknown.
    :raises ValueError: On an invalid ``k``, ``num_threads`` or
        ``chunk_size`` (including one item, ``k`` outside 1 to n - 1, or
        ``k`` above the largest NumPy array dimension);
        ``similarity=True``; an empty sequence; or a matrix or comparison
        whose distances cannot be ranked.
    :raises RuntimeError: If a sparse item has fewer than ``k`` stored
        neighbors, or a comparison returns a NaN or infinite distance.

    Example::

        graph = oecluster.knn_graph(mols, 5, comparison="fingerprint")
        print(graph.indices[0], graph.distances[0])
    """
    if similarity:
        raise ValueError(
            "knn_graph() ranks distances; similarity=True is not supported")
    k_value = _diversity_int(k, "k", 0)
    # The arrays are shaped (rows, k) even with zero rows, and NumPy refuses
    # a dimension above intp's maximum, so a zero-item call cannot accept
    # every size_t k the native layer does.
    if k_value > _INTP_MAX:
        raise ValueError(
            f"knn_graph() k must be at most {_INTP_MAX}, the largest NumPy "
            f"array dimension, got {k_value}")
    num_threads_value = _diversity_int(num_threads, "num_threads", 0)
    chunk_size_value = _diversity_int(chunk_size, "chunk_size", 1)

    source = _diversity_source(items, comparison, kwargs, "knn_graph",
                               allow_sparse=True, defer_comparable_check=True)
    _knn_check_k(k_value, source.size, "knn_graph")
    if source.size > 0:
        if source.matrix is not None:
            _gate.require_comparable(source.matrix, "knn_graph")
        else:
            _refuse_comparison_facts(source.target, "knn_graph")

    options = _oecluster.KNNGraphOptions()
    options.k = k_value
    options.num_threads = num_threads_value
    options.chunk_size = chunk_size_value
    native = _oecluster.knn_graph(source.target, options)
    return KNNGraph._from_native(native, source)


def _jarvis_patrick_check_kmin(kmin, k, size):
    """
    Refuse a shared-neighbor count no mutual pair can reach.

    :param kmin: Coerced shared-neighbor count.
    :param k: Neighbors per item.
    :param size: Item count; nothing is refused for zero items.
    :raises ValueError: If kmin >= k with items.
    """
    if size and kmin >= k:
        raise ValueError(
            f"jarvis_patrick() kmin must be less than k = {k}, got {kmin}; a "
            "mutual pair shares at most k - 1 neighbors")


def _jarvis_patrick_result(native, positions, num_positions, excluded):
    """
    Map a native result's labels and clusters back to caller positions.

    :param native: The native JarvisPatrickResult.
    :param positions: Caller position per native index, or None for identity.
    :param num_positions: The caller's item count.
    :param excluded: ``[position, reason]`` entries normalization dropped.
    :returns: A :class:`JarvisPatrickResult`.
    """
    def caller(index):
        return index if positions is None else positions[index]

    labels = [-1] * num_positions
    for index, label in enumerate(native.Labels()):
        labels[caller(index)] = int(label)
    return JarvisPatrickResult(
        labels,
        [[caller(i) for i in cluster] for cluster in native.Members()],
        k=native.K(), kmin=native.KMin(), excluded=excluded)


def jarvis_patrick(items, *, kmin, k=None, comparison=None, similarity=False,
                   num_threads=0, chunk_size=4096,
                   **kwargs) -> "JarvisPatrickResult":
    """
    Cluster by the classic Jarvis-Patrick shared-nearest-neighbor rule.

    Items i and j are linked when each is among the other's ``k`` nearest
    items and the two neighbor lists share at least ``kmin`` items; clusters
    are the connected components of the links, and an unlinked item is a
    singleton. Neighbor lists exclude the item itself, so a formulation that
    counts the item among its own neighbors uses ``k + 1`` for this ``k`` and
    ``kmin + 2`` for this ``kmin``.

    ``items`` is a :class:`KNNGraph` from :func:`knn_graph`, or any input
    :func:`knn_graph` accepts, in which case ``k`` is required and the graph
    is built first. ``k`` and ``kmin`` are keyword-only so they cannot be
    swapped by position. With a graph, ``k`` must be omitted or equal
    ``graph.k``, and ``num_threads`` and ``chunk_size`` are validated but
    unused.

    :param items: A KNNGraph, matrix, prebuilt comparison, or sequence of
        items.
    :param kmin: Shared neighbors required to link a mutual pair; below k.
    :param k: Neighbors per item, 1 to n - 1; required unless items is a
        KNNGraph.
    :param comparison: Comparison name, required with a sequence of items.
    :param similarity: Refused when True; the graph ranks distances.
    :param num_threads: Worker threads for the graph; 0 selects the hardware
        concurrency.
    :param chunk_size: Pairwise distances per work unit, at least one.
    :param kwargs: Comparison options for a named comparison.
    :returns: A :class:`JarvisPatrickResult`.
    :raises TypeError: If the arguments fit no path, comparison arguments
        accompany a KNNGraph, or ``k`` is missing with raw input.
    :raises ValueError: On an invalid ``k``, ``kmin``, ``num_threads`` or
        ``chunk_size`` (including ``kmin >= k``); a graph whose ``k``
        differs from ``k``; ``similarity=True``; an empty sequence; or a
        matrix or comparison whose distances cannot be ranked.
    :raises RuntimeError: If a sparse item has fewer than ``k`` stored
        neighbors, or a comparison returns a NaN or infinite distance.

    Example::

        result = oecluster.jarvis_patrick(mols, k=6, kmin=3,
                                          comparison="fingerprint")
        print(result.clusters)
    """
    if similarity:
        raise ValueError(
            "jarvis_patrick() clusters on distances; similarity=True is not "
            "supported")
    kmin_value = _diversity_int(kmin, "kmin", 0)
    num_threads_value = _diversity_int(num_threads, "num_threads", 0)
    chunk_size_value = _diversity_int(chunk_size, "chunk_size", 1)

    if isinstance(items, KNNGraph):
        if comparison is not None or kwargs:
            raise TypeError(
                "jarvis_patrick() takes no comparison or comparison options "
                "with a KNNGraph, which already fixes the neighbors")
        if k is not None:
            k_value = _diversity_int(k, "k", 0)
            if k_value != items.k:
                raise ValueError(
                    f"jarvis_patrick() k={k_value} does not match the graph's "
                    f"k={items.k}")
        _jarvis_patrick_check_kmin(kmin_value, items.k, len(items))
        native = _oecluster.jarvis_patrick(items._native, kmin_value)  # type: ignore[attr-defined]
        return _jarvis_patrick_result(native, items._positions,  # type: ignore[attr-defined]
                                      items._num_positions, items._excluded)  # type: ignore[attr-defined]

    if k is None:
        raise TypeError(
            "jarvis_patrick() requires k= unless items is a KNNGraph")
    k_value = _diversity_int(k, "k", 0)
    source = _diversity_source(items, comparison, kwargs, "jarvis_patrick",
                               allow_sparse=True, defer_comparable_check=True)
    _knn_check_k(k_value, source.size, "jarvis_patrick")
    _jarvis_patrick_check_kmin(kmin_value, k_value, source.size)
    if source.size > 0:
        if source.matrix is not None:
            _gate.require_comparable(source.matrix, "jarvis_patrick")
        else:
            _refuse_comparison_facts(source.target, "jarvis_patrick")

    options = _oecluster.JarvisPatrickOptions()
    options.k = k_value
    options.kmin = kmin_value
    options.num_threads = num_threads_value
    options.chunk_size = chunk_size_value
    native = _oecluster.jarvis_patrick(source.target, options)
    return _jarvis_patrick_result(native, source.positions,
                                  source.num_positions, source.excluded)


# leiden labels clusters with int and indexes nodes with uint32_t; this is
# the native INT_MAX limit in src/clustering/SNNWeights.h.
_LEIDEN_MAX_ITEMS = 2147483647


def _leiden_check_size(size):
    """
    Refuse more items than leiden can label.

    Mirrors validate_leiden_item_count in src/clustering/SNNWeights.h ahead
    of it, because SWIG turns the native invalid_argument into RuntimeError.

    :param size: Native item count.
    :raises ValueError: If size exceeds INT_MAX.
    """
    if size > _LEIDEN_MAX_ITEMS:
        raise ValueError(
            f"leiden() supports at most {_LEIDEN_MAX_ITEMS} items, got {size}")


def _leiden_int(value, name, minimum, maximum):
    """
    Coerce a bounded integer option of leiden().

    :param value: Caller value.
    :param name: Option name for the messages.
    :param minimum: Smallest accepted value.
    :param maximum: Largest accepted value.
    :returns: The coerced int.
    :raises ValueError: If the value is not an integer or is out of range.
    """
    # Bounded here rather than by _diversity_int because n_iterations is
    # signed and seed spans the full uint64_t range; SWIG must never see a
    # value it cannot convert.
    if isinstance(value, bool):
        raise ValueError(f"leiden() {name} must be an integer, got {value!r}")  # noqa: TRY004
    try:
        coerced = operator.index(value)
    except TypeError:
        raise ValueError(
            f"leiden() {name} must be an integer, got {value!r}") from None
    if not minimum <= coerced <= maximum:
        raise ValueError(
            f"leiden() {name} must be between {minimum} and {maximum}, "
            f"got {coerced}")
    return coerced


def _leiden_real(value, name, requirement, accepts):
    """
    Coerce a finite real option of leiden().

    :param value: Caller value.
    :param name: Option name for the messages.
    :param requirement: Phrase completing "must be" in the message.
    :param accepts: Predicate a finite value must also satisfy.
    :returns: The coerced float.
    :raises ValueError: If the value is NaN, infinite or refused by accepts.
    """
    try:
        coerced = float(value)
    except OverflowError:
        # An int beyond double range is the infinity refused below, failing
        # one step earlier.
        coerced = math.inf
    if not (math.isfinite(coerced) and accepts(coerced)):
        raise ValueError(
            f"leiden() {name} must be {requirement}, got {value!r}")
    return coerced


def _leiden_objective(objective):
    """
    Map an objective name to the native enum, matched exactly.

    :param objective: ``"modularity"`` or ``"cpm"``.
    :returns: The native LeidenObjective constant.
    :raises ValueError: For any other value.
    """
    objectives = {"modularity": _oecluster.LeidenObjective_Modularity,
                  "cpm": _oecluster.LeidenObjective_CPM}
    if not isinstance(objective, str) or objective not in objectives:
        raise ValueError(
            "leiden() objective must be 'modularity' or 'cpm', got "
            f"{objective!r}")
    return objectives[objective]


def _leiden_result(native, positions, num_positions, excluded):
    """
    Map a native result's labels and clusters back to caller positions.

    :param native: The native LeidenResult.
    :param positions: Caller position per native index, or None for identity.
    :param num_positions: The caller's item count.
    :param excluded: ``[position, reason]`` entries normalization dropped.
    :returns: A :class:`LeidenResult`.
    """
    def caller(index):
        return index if positions is None else positions[index]

    labels = [-1] * num_positions
    for index, label in enumerate(native.Labels()):
        labels[caller(index)] = int(label)
    objective = ("cpm" if native.Objective() == _oecluster.LeidenObjective_CPM
                 else "modularity")
    return LeidenResult(
        labels,
        [[caller(i) for i in cluster] for cluster in native.Members()],
        quality=native.Quality(), iterations=native.Iterations(),
        objective=objective, resolution=native.Resolution(), k=native.K(),
        excluded=excluded)


def leiden(items, *, k=None, objective="modularity", resolution=1.0,
           prune=1 / 15, theta=0.01, n_iterations=-1, seed=0,
           comparison=None, similarity=False, num_threads=0,
           chunk_size=4096, **kwargs) -> "LeidenResult":
    """
    Partition items by Leiden community detection on a shared-neighbor graph.

    The k-nearest-neighbor graph is reweighted by shared-nearest-neighbor
    Jaccard: with N+(i) the ``k`` neighbors of i plus i itself, the edge
    {i, j} weighs s / (2(k + 1) - s) for s = |N+(i) & N+(j)|. Only pairs in
    which one item names the other get an edge, and edges below ``prune``
    are dropped. Leiden (Traag, Waltman and van Eck, 2019) then optimizes
    ``objective`` on that graph; every returned cluster is connected.

    ``items`` is a :class:`KNNGraph` from :func:`knn_graph`, or any input
    :func:`knn_graph` accepts, in which case ``k`` is required and the graph
    is built first. With a graph, ``k`` must be omitted or equal
    ``graph.k``. The same inputs, options and build give the same result,
    whatever ``num_threads`` is.

    :param items: A KNNGraph, matrix, prebuilt comparison, or sequence of
        items.
    :param k: Neighbors per item, 1 to n - 1; required unless items is a
        KNNGraph. Seurat's ``k.param`` is this ``k + 1``.
    :param objective: ``"modularity"`` (Reichardt-Bornholdt, with
        ``resolution``) or ``"cpm"`` (Constant Potts Model).
    :param resolution: Finite and non-negative. Higher values give smaller
        clusters. 1.0 is the modularity default; CPM with Jaccard weights
        typically needs a value well below the median edge weight.
    :param prune: Jaccard weights below this are dropped; in [0, 1).
    :param theta: Randomness of the refinement step; finite and positive.
    :param n_iterations: -1 repeats passes until one changes nothing;
        otherwise exactly that many passes, up to 2**63 - 1.
    :param seed: Seed for the random stream, 0 to 2**64 - 1.
    :param comparison: Comparison name, required with a sequence of items.
    :param similarity: Refused when True; the graph ranks distances.
    :param num_threads: Worker threads for the graph and the weights; 0
        selects the hardware concurrency.
    :param chunk_size: Pairwise distances per work unit, at least one.
    :param kwargs: Comparison options for a named comparison.
    :returns: A :class:`LeidenResult`.
    :raises TypeError: If the arguments fit no path, comparison arguments
        accompany a KNNGraph, or ``k`` is missing with raw input.
    :raises ValueError: On an invalid option; a graph whose ``k`` differs
        from ``k``; ``similarity=True``; more than 2147483647 items; an
        empty sequence; or a matrix or comparison whose distances cannot be
        ranked.
    :raises RuntimeError: If a sparse item has fewer than ``k`` stored
        neighbors, or a comparison returns a NaN or infinite distance.

    Example::

        result = oecluster.leiden(mols, k=15, comparison="fingerprint")
        print(result.clusters, result.quality)
    """
    if similarity:
        raise ValueError(
            "leiden() clusters on distances; similarity=True is not "
            "supported")
    options = _oecluster.LeidenOptions()
    options.objective = _leiden_objective(objective)
    options.resolution = _leiden_real(resolution, "resolution",
                                      "finite and non-negative",
                                      lambda value: value >= 0.0)
    options.prune = _leiden_real(prune, "prune", "finite and in [0, 1)",
                                 lambda value: 0.0 <= value < 1.0)
    options.theta = _leiden_real(theta, "theta", "finite and positive",
                                 lambda value: value > 0.0)
    options.n_iterations = _leiden_int(n_iterations, "n_iterations", -1,
                                       2**63 - 1)
    options.seed = _leiden_int(seed, "seed", 0, 2**64 - 1)
    options.num_threads = _diversity_int(num_threads, "num_threads", 0)
    options.chunk_size = _diversity_int(chunk_size, "chunk_size", 1)

    if isinstance(items, KNNGraph):
        if comparison is not None or kwargs:
            raise TypeError(
                "leiden() takes no comparison or comparison options with a "
                "KNNGraph, which already fixes the neighbors")
        if k is not None:
            k_value = _diversity_int(k, "k", 0)
            if k_value != items.k:
                raise ValueError(
                    f"leiden() k={k_value} does not match the graph's "
                    f"k={items.k}")
        _leiden_check_size(len(items))
        native = _oecluster.leiden(items._native, options)  # type: ignore[attr-defined]
        return _leiden_result(native, items._positions,  # type: ignore[attr-defined]
                              items._num_positions, items._excluded)  # type: ignore[attr-defined]

    if k is None:
        raise TypeError("leiden() requires k= unless items is a KNNGraph")
    k_value = _diversity_int(k, "k", 0)
    source = _diversity_source(items, comparison, kwargs, "leiden",
                               allow_sparse=True, defer_comparable_check=True)
    _leiden_check_size(source.size)
    _knn_check_k(k_value, source.size, "leiden")
    if source.size > 0:
        if source.matrix is not None:
            _gate.require_comparable(source.matrix, "leiden")
        else:
            _refuse_comparison_facts(source.target, "leiden")

    options.k = k_value
    native = _oecluster.leiden(source.target, options)
    return _leiden_result(native, source.positions, source.num_positions,
                          source.excluded)


def descriptor_statistics(mols, *, sources=None, columns=None, groups=None,
                          inverse_covariance=False):
    """
    Compute per-column descriptor statistics over a molecule set.

    The statistics are the ones the descriptor comparison fits internally, so
    computing them here and passing them back to :func:`pdist` gives a
    reusable, explicitly scoped standardization. Two arguments have to travel
    with ``variances=`` or ``inverse_covariance=``: ``columns=stats['columns']``,
    because the values are matched to columns by position and any dropped
    column leaves them shorter than the selection :func:`pdist` would otherwise
    make, and the same ``sources=`` these statistics were fitted over, because
    a column name from one source is not in another source's schema.

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
        ``inverse_covariance_rows``. The last is the number of rows behind the
        fitted matrix: ``0`` when no matrix was requested, and not always
        ``num_rows``.
    :raises TypeError: Among the reasons: ``sources``, ``columns`` or
        ``groups`` given a value that cannot be iterated. All three are
        converted before any of them is judged empty, so a conversion failure
        precedes the empty-sequence refusal whichever option carries which
        fault.
    :raises ValueError: If ``sources``, ``columns`` or ``groups`` is passed as
        an empty sequence, which C++ cannot tell apart from an omitted option.
        Decided before the request is computed, so it still precedes the
        refusals below rather than following them: there is no argument-level
        validator for these options -- the name verdicts come out of the
        computing call itself -- and ordering one of them first would mean
        computing and discarding a full descriptor table on a call about to be
        refused.
    :raises RuntimeError: If the descriptor layer refuses the request. Among
        the reasons: an unknown source, column, or group name; and fewer than
        two molecules, which is too few to fit a variance.
    """
    options = _oecluster.DescriptorStatisticsOptions()
    # Convert every sequence option before judging any of them, so that an
    # early empty one cannot stop a later one from being converted at all. The
    # emptiness rule is this layer's own advisory rule and a conversion failure
    # is a verdict on what was passed, so the second has to be able to win.
    # ``_descriptor_options_and_empties`` gives the comparison paths the same
    # shape; only the vector conversion and the message are shared with it,
    # because it fills a ``DescriptorOptions`` from a kwargs dict of five
    # sequence options and these are three explicit keyword parameters on a
    # ``DescriptorStatisticsOptions``.
    empty = []
    for name, value in (('sources', sources), ('columns', columns),
                        ('groups', groups)):
        if value is None:
            continue
        vector = _comparisons._string_vector(value)
        setattr(options, name, vector)
        if len(vector) == 0:
            empty.append(name)
    options.inverse_covariance = _flag(inverse_covariance, "inverse_covariance")

    # Still ahead of every name verdict, which only the computing call below
    # can give: ordering one of those first would compute and discard a full
    # descriptor table on a call about to be refused.
    if empty:
        raise _comparisons._empty_option_error(empty[0])

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

    def __new__(cls, mols, *, fp_type=None, storage=None, metric=None,
                numbits=None, radius=None, min_distance=None,
                max_distance=None, torsion_atom_count=None,
                use_chirality=None, p=None, tversky_alpha=None,
                tversky_beta=None, similarity=False):
        """
        Construct a FingerprintComparison.

        :param mols: List of OEMolBase molecules.
        :param fp_type: Fingerprint family: "morgan" (the default),
            "atom_pair", "topological_atom_pair", or "topological_torsions".
            "topological_atom_pair" selects the same generator as "atom_pair".
        :param storage: Fingerprint storage: "binary" (the default), "count",
            "sparse", or "sparse_count". A counted storage needs a metric
            defined on counts, such as "bray_curtis" or "manhattan".
        :param metric: OEFP scalar metric name. Defaults to "tanimoto". Every
            metric but "tversky" has a distance form; "tversky" has only a
            similarity form and requires ``similarity=True``.
        :param numbits: Fingerprint size in bits. Defaults to 2048. A sparse
            storage keeps its family's own domain rather than folding to a
            width, and naming this alongside one raises.
        :param radius: Morgan radius. Defaults to 2.
        :param min_distance: Minimum atom-pair graph distance. Defaults to 1.
        :param max_distance: Maximum atom-pair graph distance. Defaults to 30.
            This no longer sets the Morgan radius; naming it with
            ``fp_type="morgan"`` raises. Use ``radius``.
        :param torsion_atom_count: Torsion path length. Defaults to 4.
        :param use_chirality: Distinguish stereocenters. Defaults to False.
        :param p: Minkowski order. Defaults to 2.0.
        :param tversky_alpha: Tversky reference weight. Defaults to 0.5.
        :param tversky_beta: Tversky fit weight. Defaults to 0.5. Unequal
            ``tversky_alpha`` and ``tversky_beta`` build a valid asymmetric
            comparison, but the object is good only for direct ``Compare(i, j)``
            calls: ``pdist()`` then refuses it for the asymmetry, and
            ``cdist()`` takes no prebuilt comparison object at all. Ask
            ``pdist()``/``cdist()`` for the same configuration by keyword
            instead if you want the diagnostic that names the two weights and
            points at ``cdist()``.
        :param similarity: Return similarity instead of distance.
        :returns: C++ FingerprintComparison object.
        :raises TypeError: If a named option does not apply to the selected
            family, storage, or metric.
        :raises RuntimeError: If the C++ layer refuses the request. Among the
            reasons: an unknown family, storage, or metric; a metric with no
            similarity form under ``similarity=True``, or none with a distance
            form under ``similarity=False``; the one family and storage
            combination OEFP provides no batch type for; and a bit-set metric
            on a counted storage.
        """
        # Delegated to the keyword surface's own builder rather than
        # reimplemented, so this path cannot become the way around the
        # explicitness rules. The ordering is part of what has to match: the
        # builder constructs the C++ object before applying any advisory rule,
        # which is what keeps an advisory refusal from pre-empting an
        # authoritative one.
        #
        # symmetric=False because this builds a comparison, not a pdist. An
        # asymmetric Tversky is a legal comparison -- Compare(i, j) and
        # Compare(j, i) differ, which is the point of it -- and C++ refuses it
        # only at FingerprintComparison::TryPDist, where pdist meets it.
        # Refusing here instead would deny the direct-Compare caller an object
        # the C++ constructor accepts.
        comparison, _ = _comparisons._build_fingerprint(
            mols, similarity,
            {
                'fp_type': fp_type,
                'storage': storage,
                'metric': metric,
                'numbits': numbits,
                'radius': radius,
                'min_distance': min_distance,
                'max_distance': max_distance,
                'torsion_atom_count': torsion_atom_count,
                'use_chirality': use_chirality,
                'p': p,
                'tversky_alpha': tversky_alpha,
                'tversky_beta': tversky_beta,
            },
            symmetric=False)
        return comparison


class ROCSComparison:
    """ROCS-style shape overlay comparison."""

    def __new__(cls, mols, *, similarity=False):
        """
        Construct a ROCSComparison.

        :param mols: List of OEMolBase molecules with 3D coordinates. The
            first element decides how the whole list is read: if it is an
            ``OEGraphMol``, every later ``OEMol`` is silently reduced to its
            active conformer, and the reverse ordering raises ``TypeError``.
            Convert with ``oechem.OEMol(...)`` first when the list would
            otherwise be mixed.
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
        :raises ValueError: If ``sources``, ``columns``, ``groups``,
            ``variances`` or ``inverse_covariance`` is passed as an empty
            sequence, which C++ cannot tell apart from an omitted option.
            Reported after a rejected option value, which the caller has to fix
            whatever they do about the empty sequence.
        :raises RuntimeError: If the C++ layer refuses the request. Among the
            reasons: a rejected option value; a molecule with an absent or
            non-finite value for a selected descriptor under the
            "complete_case" policy; and too few molecules for a metric that
            fits its variances from the input.
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
        opts, empty = _comparisons._descriptor_options_and_empties(kwargs)
        if empty:
            # The C++ constructor gives the verdicts on option values itself,
            # but it is never reached on a call about to be refused. Running
            # the same check here is what keeps the emptiness refusal last,
            # so a caller who also misspelled the metric is told about the
            # metric first. Only on this path, so the common case does not pay
            # for a validation the constructor is about to repeat.
            _oecluster.validate_descriptor_options(opts)
            raise _comparisons._empty_option_error(empty[0])

        return _DescriptorComparison(mols, opts)


class RMSDComparison:
    """Coordinate RMSD between poses of one molecule."""

    def __new__(cls, mols, *, overlay=None, automorph=None, heavy_only=None):
        """
        Construct an RMSDComparison.

        The molecules are taken as given. Use ``pdist(mols, "rmsd")`` instead
        if the poses need conformer expansion first, since this factory
        receives the molecules exactly as passed.

        :param mols: List of OEMol molecules with coordinates.
        :param overlay: Superpose before measuring, removing rigid-body
            differences. Defaults to False (in-frame RMSD).
        :param automorph: Minimize over graph automorphisms. Defaults to True.
        :param heavy_only: Ignore hydrogens. Defaults to True.
        :returns: C++ RMSDComparison object.
        :raises RuntimeError: If the C++ layer refuses the request. Among the
            reasons: a molecule with no coordinates; a mix of 2D and 3D input;
            molecules that do not share one topology; with
            ``heavy_only=False``, a differing atom count once hydrogens are
            counted; and, with ``automorph=False``, a differing atom count,
            element order, or bond set.
        """
        kwargs = {
            'overlay': overlay,
            'automorph': automorph,
            'heavy_only': heavy_only,
        }
        return _RMSDComparison(mols, _comparisons.rmsd_options(kwargs))


class MCSComparison:
    """Maximum common substructure, scored as Tanimoto over matched bonds."""

    def __new__(cls, mols, *, search_mode=None, match_level=None,
                max_matches=None, similarity=False):
        """
        Construct an MCSComparison.

        Coordinates are never read, so molecules parsed from SMILES need no
        embedding step. Suppression folds a hydrogen into the implicit
        hydrogen count of the atom it hangs off, so isotopic hydrogens go too
        and a deuterated analogue scores against its parent as identical. A
        hydrogen whose fold has nowhere to go stays explicit, which in
        practice means a bridging hydrogen or one carrying a formal charge.

        Approximate search is asymmetric, so each pair is searched in both
        directions and the larger match wins. There is no metric guarantee: no
        triangle-inequality violation has been observed, but none is proven
        either, so the matrix reports ``triangle`` as ``"unknown"``.

        :param mols: List of OEMolBase molecules. The first element decides
            which mixed lists are accepted, as for ``ROCSComparison``, but the
            scores are unchanged because coordinates are never read.
        :param search_mode: ``"approximate"`` (default) or ``"exhaustive"``.
            Exhaustive is one to three orders of magnitude slower and is not
            reliably better: on a 53-bond against 54-bond macrolide pair it took
            16.2 s and matched 50 bonds where approximate took 8.5 ms and
            matched 51.
        :param match_level: ``"default"``, ``"exact"`` or ``"loose"``, setting
            how strictly atoms and bonds must correspond.
        :param max_matches: How many matches one directed search may enumerate.
            Defaults to 1024. Quality saturated by 256 on all 45 pairs
            measured; small values are destructive rather than merely faster.
            A bool is refused, because ``True`` would otherwise be read as a
            budget of 1.
        :param similarity: Return the bond Tanimoto rather than one minus it.
        :returns: C++ MCSComparison object.
        :raises RuntimeError: If the C++ layer refuses the request. Among the
            reasons: a molecule with no bonds after hydrogen suppression, such
            as methane, water or argon, for which bond Tanimoto has a zero
            denominator; and ``max_matches=0``, which the toolkit would read
            as a budget of zero matches.
        :raises TypeError: If ``max_matches`` is a bool, or if an option value
            is of a type the options struct will not take.
        :raises OverflowError: If ``max_matches`` is outside the range of
            an ``unsigned int``, such as ``-1`` or ``2**40``.
        :raises ValueError: If ``search_mode`` or ``match_level`` names a mode
            the comparison does not have.
        """
        kwargs = {
            'search_mode': search_mode,
            'match_level': match_level,
            'max_matches': max_matches,
        }
        opts = _comparisons.mcs_options(kwargs)
        # similarity is not an _MCS_KEYS member, so mcs_options leaves the
        # struct default of False. This path does not go through _build_mcs, so
        # it carries its own copy of the assignment; without it,
        # MCSComparison(mols, similarity=True) would return distances and report
        # is_distance = Yes with no error anywhere. Assigned raw, matching
        # ROCSComparison and SuperposeComparison, so the SWIG bool typemap
        # refuses a string the caller meant as false.
        opts.similarity = similarity
        return _MCSComparison(mols, opts)
