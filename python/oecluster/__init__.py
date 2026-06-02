"""
oecluster: High-performance pairwise distance computation for molecular datasets.

This package provides efficient computation of pairwise distance matrices for
molecular and protein structure datasets using OpenEye toolkits. It supports
multiple comparison methods including fingerprint similarity, ROCS shape
overlay, protein superposition, and binding site comparison.
"""

import json
import hashlib
import importlib.machinery
import importlib.util
import os
import re
import shutil
import sys
import warnings
import ctypes
from importlib import metadata
from pathlib import Path
from typing import Any

import numpy as np

__version__ = "4.0.0"
__version_info__ = (4, 0, 0)


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
    "ClusteringResult",
    "ButinaResult",
    "DBSCANResult",
    "HDBSCANResult",
    "AgglomerativeResult",
    "BitBirchResult",
    "RepresentativeMetrics",
    "ClusterRepresentative",
    "pdist",
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

from .oecluster import FingerprintComparison as _FingerprintComparison
from .oecluster import FingerprintOptions
from .oecluster import ROCSComparison as _ROCSComparison
from .oecluster import ROCSOptions
from .oecluster import SuperposeComparison as _SuperposeComparison
from .oecluster import SuperposeOptions


class DistanceMatrix:
    """
    A distance matrix computed from a set of items using a comparison method.

    Provides zero-copy access to condensed distance matrix data via numpy arrays,
    conversion to scipy sparse matrices, and serialization support.
    """

    def __init__(self, storage, comparison_name, labels=None, params=None):
        """
        Construct a DistanceMatrix wrapper.

        :param storage: C++ StorageBackend instance.
        :param comparison_name: Name of the comparison method used.
        :param labels: Optional list of labels for items.
        :param params: Optional dictionary of comparison parameters.
        """
        self._storage = storage
        self._comparison_name = comparison_name
        self._labels = labels if labels is not None else []
        self._params = params if params is not None else {}
        self._condensed_cache = None

    @property
    def storage(self):
        """Get the underlying storage backend."""
        return self._storage

    @property
    def comparison_name(self):
        """Get the name of the comparison method."""
        return self._comparison_name

    @property
    def params(self):
        """Get the comparison parameters dictionary."""
        return self._params

    @property
    def labels(self):
        """Get the item labels."""
        return self._labels

    @property
    def num_items(self):
        """Get the number of items in the dataset."""
        return self._storage.NumItems()

    @property
    def num_pairs(self):
        """Get the number of pairwise distances."""
        return self._storage.NumPairs()

    @property
    def shape(self):
        """Get the shape of the full distance matrix (N, N)."""
        n = self.num_items
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
            n = self.num_items
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

        # Dense or MMap storage: zero-copy access
        ptr = self._storage._data_ptr()
        num_pairs = self.num_pairs
        arr = np.ctypeslib.as_array(
            ctypes.cast(ptr, ctypes.POINTER(ctypes.c_double)),
            shape=(num_pairs,)
        )

        # Cache to keep storage alive
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

        n = self.num_items

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

        :param path: Output file path.
        """
        np.savez_compressed(
            path,
            condensed=self.condensed,
            comparison_name=np.array(self._comparison_name),
            params_json=np.array(json.dumps(self._params)),
            labels=np.array(self._labels),
            num_items=np.array(self.num_items),
        )

    @classmethod
    def from_file(cls, path):
        """
        Load a distance matrix from a .npz file.

        :param path: Input file path.
        :returns: DistanceMatrix instance.
        """
        data = np.load(path, allow_pickle=False)
        condensed = data['condensed']

        comparison_name = str(data['comparison_name'])
        labels = list(data.get('labels', np.array([])))

        try:
            params = json.loads(str(data['params_json']))
        except (KeyError, json.JSONDecodeError):
            params = {}

        num_items = int(data['num_items'])
        storage = DenseStorage(num_items)

        idx = 0
        for i in range(num_items):
            for j in range(i + 1, num_items):
                storage.Set(i, j, float(condensed[idx]))
                idx += 1

        return cls(storage, comparison_name, labels, params)

    def __array__(self):
        """Support numpy array interface."""
        return self.condensed

    def __len__(self):
        """Return the number of pairwise distances."""
        return self.num_pairs

    def __repr__(self):
        return (f"DistanceMatrix(comparison={self._comparison_name!r}, "
                f"num_items={self.num_items}, num_pairs={self.num_pairs})")


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
        :param native_owner: Optional native result object that owns borrowed
            storage (e.g. BitBirch centroids), kept alive by this reference.
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
    def num_clusters(self):
        """Number of clusters."""
        return len(self._clusters)

    @property
    def num_items(self):
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
                f"num_items={self.num_items})")


class ButinaResult(ClusteringResult):
    """Butina clustering result.

    The first member of each cluster is the highest-neighborhood representative.
    """


class DBSCANResult(ClusteringResult):
    """DBSCAN clustering result with core sample indices."""

    def __init__(self, labels, clusters, *, core_sample_indices=(),
                 native_owner=None):
        super().__init__(labels, clusters, native_owner=native_owner)
        self._core_sample_indices = tuple(int(i) for i in core_sample_indices)

    @property
    def core_sample_indices(self):
        """Indices of core samples."""
        return self._core_sample_indices


class HDBSCANResult(ClusteringResult):
    """HDBSCAN clustering result with membership probabilities."""

    def __init__(self, labels, clusters, *, probabilities=None,
                 native_owner=None):
        super().__init__(labels, clusters, native_owner=native_owner)
        self._probabilities = (
            None if probabilities is None
            else np.asarray(list(probabilities), dtype=np.float64)
        )

    @property
    def probabilities(self):
        """Per-item membership probabilities as a NumPy ``float64`` array."""
        return self._probabilities


class AgglomerativeResult(ClusteringResult):
    """Agglomerative clustering result with merge-tree metadata."""

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
        """Tuple of ``(left, right)`` merge child node indices."""
        return self._children

    @property
    def distances(self):
        """Merge distances as a NumPy ``float64`` array."""
        return self._distances

    @property
    def cluster_sizes(self):
        """Merged cluster size per merge."""
        return self._cluster_sizes


class BitBirchResult(ClusteringResult):
    """BitBirch clustering result with centroid fingerprints."""

    def __init__(self, labels, clusters, *, centroids=None, cluster_sizes=(),
                 native_owner=None):
        super().__init__(labels, clusters, native_owner=native_owner)
        self._centroids = centroids
        self._cluster_sizes = tuple(int(size) for size in cluster_sizes)

    @property
    def centroids(self):
        """Centroid fingerprints (an ``oefp.OEFPBatch``)."""
        return self._centroids

    @property
    def cluster_sizes(self):
        """Member count per cluster."""
        return self._cluster_sizes


def _extract_labels(items):
    """
    Extract labels from molecular items.

    :param items: List of molecules or design units.
    :returns: List of labels (molecule titles or indices).
    """
    labels = []

    try:
        # Try to import openeye to check types
        from openeye import oechem

        for idx, item in enumerate(items):
            if isinstance(item, oechem.OEMolBase):
                title = item.GetTitle()
                labels.append(title if title else f"mol_{idx}")
            else:
                labels.append(f"item_{idx}")
    except (ImportError, AttributeError):
        # If openeye not available or not a molecule, use indices
        labels = [f"item_{idx}" for idx in range(len(items))]

    return labels


def pdist(items,
          comparison,
          *,
          similarity=False,
          num_threads=0,
          chunk_size=256,
          cutoff=0.0,
          output=None,
          progress=None,
          **kwargs) -> DistanceMatrix:
    """
    Compute pairwise distances for a collection of items using a comparison.

    :param items: List of molecules, design units, or other items.
    :param comparison: Comparison method: "fingerprint", "rocs", "superpose",
                       "sitehopper", or a C++ comparison object.
    :param similarity: Return similarities instead of distances.
    :param num_threads: Number of threads (0 = auto).
    :param chunk_size: Pairs per work unit.
    :param cutoff: Distance cutoff for sparse storage (0 = store all).
    :param output: Optional file path for memory-mapped storage.
    :param progress: Optional callback(completed, total).
    :param kwargs: Comparison-specific options.
    :returns: DistanceMatrix with computed distances/similarities.
    :raises TypeError: If unknown kwargs are passed.
    """
    if isinstance(comparison, str):
        comparison_lower = comparison.lower()
        labels = _extract_labels(items)
        params = {'comparison_type': comparison_lower, 'similarity': similarity}
        comparison_obj: Any

        if comparison_lower == "fingerprint":
            fp_opts = FingerprintOptions()
            fp_opts.similarity = similarity
            if 'fp_type' in kwargs:
                fp_opts.fp_type = kwargs.pop('fp_type')
            if 'numbits' in kwargs:
                fp_opts.numbits = kwargs.pop('numbits')
            if 'min_distance' in kwargs:
                fp_opts.min_distance = kwargs.pop('min_distance')
            if 'max_distance' in kwargs:
                fp_opts.max_distance = kwargs.pop('max_distance')
            if 'metric' in kwargs:
                fp_opts.metric = kwargs.pop('metric')
            if kwargs:
                raise TypeError(
                    f"Unknown kwargs for fingerprint comparison: {list(kwargs)}")
            comparison_obj = _FingerprintComparison(items, fp_opts)
            comparison_name = "fingerprint"

        elif comparison_lower == "rocs":
            rocs_opts = ROCSOptions()
            rocs_opts.similarity = similarity
            if 'score_type' in kwargs:
                st = kwargs.pop('score_type')
                score_map = {
                    'combo_norm': _oecluster.ROCSScoreType_ComboNorm,
                    'combo': _oecluster.ROCSScoreType_Combo,
                    'shape': _oecluster.ROCSScoreType_Shape,
                    'color': _oecluster.ROCSScoreType_Color,
                }
                if st not in score_map:
                    raise ValueError(f"Unknown ROCS score type: {st}")
                rocs_opts.score_type = score_map[st]
            if 'color_ff_type' in kwargs:
                rocs_opts.color_ff_type = kwargs.pop('color_ff_type')
            if kwargs:
                raise TypeError(
                    f"Unknown kwargs for rocs comparison: {list(kwargs)}")
            comparison_obj = _ROCSComparison(items, rocs_opts)
            comparison_name = "rocs"

        elif comparison_lower in ("superpose", "sitehopper"):
            superpose_opts = SuperposeOptions()
            superpose_opts.similarity = similarity
            method = kwargs.pop('method', None)
            if comparison_lower == "sitehopper":
                method = method or "sitehopper"
            if method is not None:
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
                superpose_opts.method = method_map[method]
            if 'score_type' in kwargs:
                st = kwargs.pop('score_type')
                st_map = {
                    'auto': _oecluster.SuperposeScoreType_Auto,
                    'rmsd': _oecluster.SuperposeScoreType_RMSD,
                    'tanimoto': _oecluster.SuperposeScoreType_Tanimoto,
                    'patch_score': _oecluster.SuperposeScoreType_PatchScore,
                }
                if st not in st_map:
                    raise ValueError(f"Unknown superpose score type: {st}")
                superpose_opts.score_type = st_map[st]
            if 'predicate' in kwargs:
                superpose_opts.predicate = kwargs.pop('predicate')
            if 'ref_predicate' in kwargs:
                superpose_opts.ref_predicate = kwargs.pop('ref_predicate')
            if 'fit_predicate' in kwargs:
                superpose_opts.fit_predicate = kwargs.pop('fit_predicate')
            if kwargs:
                raise TypeError(
                    f"Unknown kwargs for superpose comparison: {list(kwargs)}")
            comparison_obj = _SuperposeComparison(items, superpose_opts)
            comparison_name = comparison_obj.ComparisonName()

        else:
            raise ValueError(
                f"Unknown comparison: {comparison!r}. "
                f"Valid options: 'fingerprint', 'rocs', 'superpose', 'sitehopper'")

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
    return DistanceMatrix(storage, comparison_name, labels, params)


def butina(distance_matrix, threshold, *, reordering=False,
           num_threads=0, chunk_size=4096):
    """
    Cluster a precomputed distance matrix using the Butina algorithm.

    :param distance_matrix: DistanceMatrix returned by :func:`pdist`.
    :param threshold: Maximum distance for two items to be neighbors.
    :param reordering: Recompute candidate neighbor counts after each cluster.
    :param num_threads: Thread count for threshold graph construction.
    :param chunk_size: Condensed-distance pairs per work unit.
    :returns: ButinaResult with per-item labels and grouped clusters. The
        first member of each cluster is the highest-neighborhood representative,
        and each member's label equals its cluster index.
    :raises TypeError: If distance_matrix is not a DistanceMatrix.
    :raises ValueError: If threshold is negative.
    """
    if threshold < 0.0:
        raise ValueError("Butina threshold must be non-negative")
    if not isinstance(distance_matrix, DistanceMatrix):
        raise TypeError("butina() expects a DistanceMatrix")

    options = ButinaOptions()
    options.distance_threshold = float(threshold)
    options.reordering = bool(reordering)
    options.num_threads = int(num_threads)
    options.chunk_size = int(chunk_size)

    result = _butina_cluster(distance_matrix.storage, options)
    return ButinaResult(result.Labels(), result.Members())


class RepresentativeMetrics:
    """Quality metrics for one cluster representative."""

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
    """A scored cluster representative and its quality metrics."""

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
    if converted and len(converted) != distance_matrix.num_items:
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
    if converted and len(converted) != distance_matrix.num_items:
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
    if not isinstance(distance_matrix, DistanceMatrix):
        raise TypeError("representative functions expect a DistanceMatrix")

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
    :param distance_matrix: DistanceMatrix used to compute representative scores.
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
    :raises TypeError: If distance_matrix is not a DistanceMatrix.
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
    :param distance_matrix: DistanceMatrix used to compute representative scores.
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
    :raises TypeError: If distance_matrix is not a DistanceMatrix.
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
    :param distance_matrix: DistanceMatrix used to compute representative scores.
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
    :raises TypeError: If distance_matrix is not a DistanceMatrix.
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


def dbscan(distance_matrix, eps, *, min_samples=5, num_threads=0, chunk_size=4096):
    """
    Cluster a precomputed distance matrix using DBSCAN.

    :param distance_matrix: DistanceMatrix returned by :func:`pdist`.
    :param eps: Maximum distance for two items to be neighbors.
    :param min_samples: Minimum self-inclusive neighbor count for a core sample.
    :param num_threads: Thread count for threshold graph construction.
    :param chunk_size: Condensed-distance pairs per work unit.
    :returns: DBSCANResult with labels, clusters, and core sample indices.
    :raises TypeError: If distance_matrix is not a DistanceMatrix.
    :raises ValueError: If eps or min_samples are invalid.
    """
    if eps < 0.0:
        raise ValueError("DBSCAN eps must be non-negative")
    if min_samples < 1:
        raise ValueError("DBSCAN min_samples must be at least one")
    if not isinstance(distance_matrix, DistanceMatrix):
        raise TypeError("dbscan() expects a DistanceMatrix")

    options = DBSCANOptions()
    options.eps = float(eps)
    options.min_samples = int(min_samples)
    options.num_threads = int(num_threads)
    options.chunk_size = int(chunk_size)

    result = _dbscan_cluster(distance_matrix.storage, options)
    return DBSCANResult(
        result.Labels(),
        result.Members(),
        core_sample_indices=result.CoreSampleIndices(),
    )


def hdbscan(distance_matrix, *, min_cluster_size=5, min_samples=None,
            cluster_selection_epsilon=0.0, max_cluster_size=None, alpha=1.0,
            cluster_selection_method="eom", allow_single_cluster=False,
            num_threads=0, chunk_size=4096):
    """
    Cluster a precomputed distance matrix using HDBSCAN.

    :param distance_matrix: DistanceMatrix returned by :func:`pdist`.
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
    :returns: HDBSCANResult with labels, clusters, and probabilities.
    :raises TypeError: If distance_matrix is not a DistanceMatrix.
    :raises ValueError: If options are invalid.
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
    if not isinstance(distance_matrix, DistanceMatrix):
        raise TypeError("hdbscan() expects a DistanceMatrix")

    method_map = {
        "eom": _oecluster.HDBSCANClusterSelectionMethod_EOM,
        "leaf": _oecluster.HDBSCANClusterSelectionMethod_Leaf,
    }
    method_key = str(cluster_selection_method).lower()
    if method_key not in method_map:
        raise ValueError(
            f"Unknown HDBSCAN cluster_selection_method: {cluster_selection_method!r}"
        )

    options = HDBSCANOptions()
    options.min_cluster_size = int(min_cluster_size)
    options.min_samples = 0 if min_samples is None else int(min_samples)
    options.cluster_selection_epsilon = float(cluster_selection_epsilon)
    options.max_cluster_size = 0 if max_cluster_size is None else int(max_cluster_size)
    options.alpha = float(alpha)
    options.cluster_selection_method = method_map[method_key]
    options.allow_single_cluster = bool(allow_single_cluster)
    options.num_threads = int(num_threads)
    options.chunk_size = int(chunk_size)

    result = _hdbscan_cluster(distance_matrix.storage, options)
    return HDBSCANResult(
        result.Labels(),
        result.Members(),
        probabilities=result.Probabilities(),
    )


def agglomerative(distance_matrix, *, n_clusters=2, distance_threshold=None,
                  linkage="average", compute_full_tree=True,
                  num_threads=0, chunk_size=4096):
    """
    Cluster a precomputed distance matrix using hierarchical agglomerative clustering.

    :param distance_matrix: DistanceMatrix returned by :func:`pdist`.
    :param n_clusters: Number of flat clusters when distance_threshold is omitted.
    :param distance_threshold: Optional merge-distance cutoff for flat clusters.
    :param linkage: Linkage method: "single", "complete", "average", or "weighted".
    :param compute_full_tree: Whether to request full-tree computation.
    :param num_threads: Thread count for initial distance materialization.
    :param chunk_size: Rows per work unit during distance materialization.
    :returns: AgglomerativeResult with labels, clusters, children, distances, and cluster sizes.
    :raises TypeError: If distance_matrix is not a DistanceMatrix.
    :raises ValueError: If options are invalid.
    """
    if not isinstance(distance_matrix, DistanceMatrix):
        raise TypeError("agglomerative() expects a DistanceMatrix")
    if distance_threshold is None and n_clusters < 1:
        raise ValueError("Agglomerative n_clusters must be at least one")
    if distance_threshold is not None and distance_threshold < 0.0:
        raise ValueError("Agglomerative distance_threshold must be non-negative")

    linkage_map = {
        "single": _oecluster.AgglomerativeLinkageMethod_Single,
        "complete": _oecluster.AgglomerativeLinkageMethod_Complete,
        "average": _oecluster.AgglomerativeLinkageMethod_Average,
        "weighted": _oecluster.AgglomerativeLinkageMethod_Weighted,
    }
    linkage_key = str(linkage).lower()
    if linkage_key not in linkage_map:
        raise ValueError(f"Unknown agglomerative linkage: {linkage!r}")

    options = AgglomerativeOptions()
    options.n_clusters = int(n_clusters)
    options.distance_threshold = (
        -1.0 if distance_threshold is None else float(distance_threshold)
    )
    options.linkage = linkage_map[linkage_key]
    options.compute_full_tree = bool(compute_full_tree)
    options.num_threads = int(num_threads)
    options.chunk_size = int(chunk_size)

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
    :param mode: "strict_parity" or "fast"; fast currently uses the
        strict-parity path and is reserved for future optimized behavior.
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
    :param mode: "strict_parity" or "fast"; fast currently uses the
        strict-parity path and is reserved for future optimized behavior.
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
    :param mode: "strict_parity" or "fast"; fast currently uses the
        strict-parity path and is reserved for future optimized behavior.
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
