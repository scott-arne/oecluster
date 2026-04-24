"""
oecluster: High-performance pairwise distance computation for molecular datasets.

This package provides efficient computation of pairwise distance matrices for
molecular and protein structure datasets using OpenEye toolkits. It supports
multiple distance metrics including fingerprint similarity, ROCS shape overlay,
protein superposition, and binding site comparison.
"""

import json
import os
import re
import warnings
import ctypes

import numpy as np

__version__ = "3.2.0"
__version_info__ = (3, 2, 0)

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


def _ensure_library_compat():
    """Create compatibility symlinks when OpenEye library versions differ from build time."""
    try:
        from . import _build_info
    except ImportError:
        return False

    if getattr(_build_info, 'OPENEYE_LIBRARY_TYPE', 'STATIC') != 'SHARED':
        return False

    expected_libs = getattr(_build_info, 'OPENEYE_EXPECTED_LIBS', [])
    if not expected_libs:
        return False

    try:
        from openeye import libs
        oe_lib_dir = libs.FindOpenEyeDLLSDirectory()
    except (ImportError, Exception):
        return False

    if not os.path.isdir(oe_lib_dir):
        return False

    pkg_dir = os.path.dirname(__file__)
    created_any = False

    for expected_name in expected_libs:
        if os.path.exists(os.path.join(oe_lib_dir, expected_name)):
            continue

        symlink_path = os.path.join(pkg_dir, expected_name)
        if os.path.islink(symlink_path):
            if os.path.exists(symlink_path):
                continue
            try:
                os.unlink(symlink_path)
            except OSError:
                continue
        elif os.path.exists(symlink_path):
            continue

        match = re.match(r'(lib\w+?)(-[\d.]+)?(\.[\d.]*\w+)$', expected_name)
        if not match:
            continue

        base_name = match.group(1)
        actual_path = None
        for f in os.listdir(oe_lib_dir):
            if f.startswith(base_name + '-') or f.startswith(base_name + '.'):
                actual_path = os.path.join(oe_lib_dir, f)
                break

        if actual_path:
            symlink_path = os.path.join(pkg_dir, expected_name)
            try:
                os.symlink(actual_path, symlink_path)
                created_any = True
            except OSError:
                pass

    return created_any


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

    expected_libs = getattr(_build_info, 'OPENEYE_EXPECTED_LIBS', [])
    if not expected_libs:
        return

    try:
        from openeye import libs
        oe_lib_dir = libs.FindOpenEyeDLLSDirectory()
    except (ImportError, Exception):
        return

    if not os.path.isdir(oe_lib_dir):
        return

    pkg_dir = os.path.dirname(__file__)
    for lib_name in expected_libs:
        oe_path = os.path.join(oe_lib_dir, lib_name)
        local_path = os.path.join(pkg_dir, lib_name)
        path = oe_path if os.path.exists(oe_path) else local_path
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
        from openeye import oechem
        runtime_version = oechem.OEToolkitsGetRelease()
        if runtime_version and build_version:
            build_parts = build_version.split('.')[:2]
            runtime_parts = runtime_version.split('.')[:2]
            if build_parts != runtime_parts:
                warnings.warn(
                    f"OpenEye version mismatch: oecluster was built with OpenEye Toolkits {build_version} "
                    f"but runtime has OpenEye Toolkits {runtime_version}. "
                    f"This may cause compatibility issues.",
                    RuntimeWarning
                )
    except ImportError:
        warnings.warn(
            "openeye-toolkits package not found. "
            "This package requires openeye-toolkits to be installed.",
            ImportWarning
        )


# Initialize compatibility checks
_ensure_library_compat()
_preload_shared_libs()
_preload_bundled_libs()
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
            return coo_matrix((data, (rows, cols)), shape=(n, n), dtype=np.float64)

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

        return coo_matrix((data, (rows, cols)), shape=(n, n), dtype=np.float64)

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

        if metric_lower == "fingerprint":
            opts = FingerprintOptions()
            opts.similarity = similarity
            if 'fp_type' in kwargs:
                opts.fp_type = kwargs.pop('fp_type')
            if 'numbits' in kwargs:
                opts.numbits = kwargs.pop('numbits')
            if 'min_distance' in kwargs:
                opts.min_distance = kwargs.pop('min_distance')
            if 'max_distance' in kwargs:
                opts.max_distance = kwargs.pop('max_distance')
            if 'atom_type_mask' in kwargs:
                opts.atom_type_mask = kwargs.pop('atom_type_mask')
            if 'bond_type_mask' in kwargs:
                opts.bond_type_mask = kwargs.pop('bond_type_mask')
            if 'similarity_func' in kwargs:
                opts.similarity_func = kwargs.pop('similarity_func')
            if kwargs:
                raise TypeError(
                    f"Unknown kwargs for fingerprint metric: {list(kwargs)}")
            metric_obj = _FingerprintMetric(items, opts)
            metric_name = "fingerprint"

        elif metric_lower == "rocs":
            opts = ROCSOptions()
            opts.similarity = similarity
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
                opts.score_type = score_map[st]
            if 'color_ff_type' in kwargs:
                opts.color_ff_type = kwargs.pop('color_ff_type')
            if kwargs:
                raise TypeError(
                    f"Unknown kwargs for rocs metric: {list(kwargs)}")
            metric_obj = _ROCSMetric(items, opts)
            metric_name = "rocs"

        elif metric_lower in ("superpose", "sitehopper"):
            opts = SuperposeOptions()
            opts.similarity = similarity
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
                opts.method = method_map[method]
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
                opts.score_type = st_map[st]
            if 'predicate' in kwargs:
                opts.predicate = kwargs.pop('predicate')
            if 'ref_predicate' in kwargs:
                opts.ref_predicate = kwargs.pop('ref_predicate')
            if 'fit_predicate' in kwargs:
                opts.fit_predicate = kwargs.pop('fit_predicate')
            if kwargs:
                raise TypeError(
                    f"Unknown kwargs for superpose metric: {list(kwargs)}")
            metric_obj = _SuperposeMetric(items, opts)
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
