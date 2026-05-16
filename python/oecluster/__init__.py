"""
oecluster: High-performance pairwise distance computation for molecular datasets.

This package provides efficient computation of pairwise distance matrices for
molecular and protein structure datasets using OpenEye toolkits. It supports
multiple distance metrics including fingerprint similarity, ROCS shape overlay,
protein superposition, and binding site comparison.
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

__version__ = "3.3.0"
__version_info__ = (3, 3, 0)


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
    "pdist",
    "FingerprintMetric",
    "ROCSMetric",
    "SuperposeMetric",
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
        pdist as _cpp_pdist,
    )
except ImportError as e:
    raise ImportError(
        f"Failed to import _oecluster C++ extension: {e}. "
        "The package may not be built correctly."
    ) from e

from . import oecluster as _oecluster

from .oecluster import FingerprintMetric as _FingerprintMetric
from .oecluster import FingerprintOptions
from .oecluster import ROCSMetric as _ROCSMetric
from .oecluster import ROCSOptions
from .oecluster import SuperposeMetric as _SuperposeMetric
from .oecluster import SuperposeOptions


class DistanceMatrix:
    """
    A distance matrix computed from a set of items using a specific metric.

    Provides zero-copy access to condensed distance matrix data via numpy arrays,
    conversion to scipy sparse matrices, and serialization support.
    """

    def __init__(self, storage, metric_name, labels=None, params=None):
        """
        Construct a DistanceMatrix wrapper.

        :param storage: C++ StorageBackend instance.
        :param metric_name: Name of the distance metric used.
        :param labels: Optional list of labels for items.
        :param params: Optional dictionary of metric parameters.
        """
        self._storage = storage
        self._metric_name = metric_name
        self._labels = labels if labels is not None else []
        self._params = params if params is not None else {}
        self._condensed_cache = None

    @property
    def storage(self):
        """Get the underlying storage backend."""
        return self._storage

    @property
    def metric_name(self):
        """Get the name of the distance metric."""
        return self._metric_name

    @property
    def params(self):
        """Get the metric parameters dictionary."""
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
            metric_name=np.array(self._metric_name),
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

        metric_name = str(data['metric_name'])
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

        return cls(storage, metric_name, labels, params)

    def __array__(self):
        """Support numpy array interface."""
        return self.condensed

    def __len__(self):
        """Return the number of pairwise distances."""
        return self.num_pairs

    def __repr__(self):
        return (f"DistanceMatrix(metric={self._metric_name!r}, "
                f"num_items={self.num_items}, num_pairs={self.num_pairs})")


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
          metric,
          *,
          similarity=False,
          num_threads=0,
          chunk_size=256,
          cutoff=0.0,
          output=None,
          progress=None,
          **kwargs) -> DistanceMatrix:
    """
    Compute pairwise distances for a collection of items using a specified metric.

    :param items: List of molecules, design units, or other items.
    :param metric: Distance metric: "fingerprint", "rocs", "superpose",
                   "sitehopper", or a C++ metric object.
    :param similarity: Return similarities instead of distances.
    :param num_threads: Number of threads (0 = auto).
    :param chunk_size: Pairs per work unit.
    :param cutoff: Distance cutoff for sparse storage (0 = store all).
    :param output: Optional file path for memory-mapped storage.
    :param progress: Optional callback(completed, total).
    :param kwargs: Metric-specific options.
    :returns: DistanceMatrix with computed distances/similarities.
    :raises TypeError: If unknown kwargs are passed.
    """
    if isinstance(metric, str):
        metric_lower = metric.lower()
        labels = _extract_labels(items)
        params = {'metric_type': metric_lower, 'similarity': similarity}
        metric_obj: Any

        if metric_lower == "fingerprint":
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
            if 'similarity_func' in kwargs:
                fp_opts.similarity_func = kwargs.pop('similarity_func')
            if kwargs:
                raise TypeError(
                    f"Unknown kwargs for fingerprint metric: {list(kwargs)}")
            metric_obj = _FingerprintMetric(items, fp_opts)
            metric_name = "fingerprint"

        elif metric_lower == "rocs":
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
                    f"Unknown kwargs for rocs metric: {list(kwargs)}")
            metric_obj = _ROCSMetric(items, rocs_opts)
            metric_name = "rocs"

        elif metric_lower in ("superpose", "sitehopper"):
            superpose_opts = SuperposeOptions()
            superpose_opts.similarity = similarity
            method = kwargs.pop('method', None)
            if metric_lower == "sitehopper":
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
                    f"Unknown kwargs for superpose metric: {list(kwargs)}")
            metric_obj = _SuperposeMetric(items, superpose_opts)
            metric_name = metric_obj.Name()

        else:
            raise ValueError(
                f"Unknown metric: {metric!r}. "
                f"Valid options: 'fingerprint', 'rocs', 'superpose', 'sitehopper'")

    else:
        metric_obj = metric
        metric_name = metric_obj.Name()
        labels = []
        params = {}

    n = metric_obj.Size()

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

    _cpp_pdist(metric_obj, storage, options)
    return DistanceMatrix(storage, metric_name, labels, params)


# Python wrapper classes for metric construction
class FingerprintMetric:
    """Fingerprint-based distance metric using Tanimoto similarity."""

    def __new__(cls, mols, *, fp_type=None, similarity=False):
        """
        Construct a FingerprintMetric.

        :param mols: List of OEMolBase molecules.
        :param fp_type: Fingerprint type (unsigned int).
        :param similarity: Return similarity instead of distance.
        :returns: C++ FingerprintMetric object.
        """
        opts = FingerprintOptions()
        opts.similarity = similarity
        if fp_type is not None:
            opts.fp_type = fp_type
        return _FingerprintMetric(mols, opts)


class ROCSMetric:
    """ROCS-style shape overlay distance metric."""

    def __new__(cls, mols, *, similarity=False):
        """
        Construct a ROCSMetric.

        :param mols: List of OEMol molecules with 3D coordinates.
        :param similarity: Return similarity instead of distance.
        :returns: C++ ROCSMetric object.
        """
        opts = ROCSOptions()
        opts.similarity = similarity
        return _ROCSMetric(mols, opts)


class SuperposeMetric:
    """Protein superposition distance metric using oespruce OESuperpose."""

    def __new__(cls, items, *, method="global_carbon_alpha", similarity=False,
                predicate=None, ref_predicate=None, fit_predicate=None):
        """
        Construct a SuperposeMetric.

        :param items: List of OEDesignUnit or OEMolBase objects.
        :param method: Superposition method name.
        :param similarity: Return similarity instead of distance.
        :param predicate: oeselect expression for both ref and fit.
        :param ref_predicate: Override predicate for ref.
        :param fit_predicate: Override predicate for fit.
        :returns: C++ SuperposeMetric object.
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
        return _SuperposeMetric(items, opts)
