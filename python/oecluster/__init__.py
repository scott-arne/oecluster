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

__version__ = "0.1.0"

__all__ = [
    "__version__",
    "DenseStorage",
    "MMapStorage",
    "SparseStorage",
    "PDistOptions",
    "DistanceMatrix",
    "pdist",
    "FingerprintMetric",
    "ROCSMetric",
    "SuperposeMetric",
    "SiteHopperMetric",
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
            build_parts = build_version.split('.')[:3]
            runtime_parts = runtime_version.split('.')[:3]
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
_check_openeye_version()


# Import C++ bindings from SWIG module
try:
    from ._oecluster import (
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


# Conditionally import metric classes
_HAVE_FINGERPRINT = False
_HAVE_ROCS = False
_HAVE_SUPERPOSE = False
_HAVE_SITEHOPPER = False

try:
    from ._oecluster import FingerprintMetric as _FingerprintMetric
    from ._oecluster import FingerprintOptions
    _HAVE_FINGERPRINT = True
except ImportError:
    pass

try:
    from ._oecluster import ROCSMetric as _ROCSMetric
    from ._oecluster import ROCSOptions
    _HAVE_ROCS = True
except ImportError:
    pass

try:
    from ._oecluster import SuperposeMetric as _SuperposeMetric
    from ._oecluster import SuperposeOptions
    _HAVE_SUPERPOSE = True
except ImportError:
    pass

try:
    from ._oecluster import SiteHopperMetric as _SiteHopperMetric
    from ._oecluster import SiteHopperOptions
    _HAVE_SITEHOPPER = True
except ImportError:
    pass


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
          num_threads=0,
          chunk_size=256,
          cutoff=0.0,
          output=None,
          progress=None) -> DistanceMatrix:
    """
    Compute pairwise distances for a collection of items using a specified metric.

    :param items: List of molecules, design units, or other items (type depends on metric).
    :param metric: Distance metric specification. Can be:
                   - A string: "fingerprint", "rocs", "superpose", "sitehopper"
                   - A C++ metric object (FingerprintMetric, ROCSMetric, etc.)
    :param num_threads: Number of threads to use (0 = auto).
    :param chunk_size: Size of work chunks for parallel processing.
    :param cutoff: Distance cutoff for sparse storage (only store distances <= cutoff, 0 = store all).
    :param output: Optional file path for memory-mapped storage.
    :param progress: Optional progress callback function(completed, total).
    :returns: DistanceMatrix with computed distances.
    :raises ValueError: If metric string is not recognized or not available.
    :raises TypeError: If items type doesn't match the metric.

    Example::

        from openeye import oechem
        import oecluster

        # Load molecules
        mols = []
        ifs = oechem.oemolistream("molecules.sdf")
        for mol in ifs.GetOEMols():
            mols.append(oechem.OEMol(mol))

        # Compute fingerprint distances
        dm = oecluster.pdist(mols, "fingerprint", num_threads=4)

        # Access condensed array
        distances = dm.condensed

        # Convert to square form
        square = dm.squareform()
    """
    # Handle string metric specifications
    if isinstance(metric, str):
        metric_lower = metric.lower()
        labels = _extract_labels(items)
        params = {'metric_type': metric_lower}

        if metric_lower == "fingerprint":
            if not _HAVE_FINGERPRINT:
                raise ValueError(
                    "FingerprintMetric is not available. "
                    "oecluster was not compiled with OEGraphSim support."
                )
            metric_obj = _FingerprintMetric(items)
            metric_name = "fingerprint"

        elif metric_lower == "rocs":
            if not _HAVE_ROCS:
                raise ValueError(
                    "ROCSMetric is not available. "
                    "oecluster was not compiled with OEShape support."
                )
            metric_obj = _ROCSMetric(items)
            metric_name = "rocs"

        elif metric_lower == "superpose":
            if not _HAVE_SUPERPOSE:
                raise ValueError(
                    "SuperposeMetric is not available. "
                    "oecluster was not compiled with OEBio support."
                )
            metric_obj = _SuperposeMetric(items)
            metric_name = "superpose"

        elif metric_lower == "sitehopper":
            if not _HAVE_SITEHOPPER:
                raise ValueError(
                    "SiteHopperMetric is not available. "
                    "oecluster was not compiled with OEBio support."
                )
            metric_obj = _SiteHopperMetric(items)
            metric_name = "sitehopper"

        else:
            raise ValueError(
                f"Unknown metric: {metric!r}. "
                f"Valid options are: 'fingerprint', 'rocs', 'superpose', 'sitehopper'"
            )

    else:
        # Assume metric is already a C++ metric object
        metric_obj = metric
        metric_name = metric_obj.Name()
        labels = []
        params = {}

    # Determine storage backend
    n = metric_obj.Size()

    if output is not None:
        # Memory-mapped storage
        storage = MMapStorage(output, n)
    elif cutoff > 0.0:
        # Sparse storage
        storage = SparseStorage(n, cutoff)
    else:
        # Dense storage
        storage = DenseStorage(n)

    # Build options
    options = PDistOptions()
    options.num_threads = num_threads
    options.chunk_size = chunk_size
    options.cutoff = cutoff
    if progress is not None:
        options.progress = progress

    # Call C++ pdist (which calls storage.Finalize() internally)
    _cpp_pdist(metric_obj, storage, options)

    # Wrap in DistanceMatrix
    return DistanceMatrix(storage, metric_name, labels, params)


# Python wrapper classes for metric construction
class FingerprintMetric:
    """
    Fingerprint-based distance metric using Tanimoto similarity.

    Computes distances as ``1 - Tanimoto(fp_i, fp_j)`` where fingerprints
    are generated using OpenEye's OEGraphSim toolkit.
    """

    def __new__(cls, mols, *, fp_type=None):
        """
        Construct a FingerprintMetric.

        :param mols: List of OEMolBase molecules.
        :param fp_type: Fingerprint type (unsigned int, default 105 = OEFPType::Tree).
        :returns: C++ FingerprintMetric object.
        :raises ValueError: If FingerprintMetric is not available.
        """
        if not _HAVE_FINGERPRINT:
            raise ValueError(
                "FingerprintMetric is not available. "
                "oecluster was not compiled with OEGraphSim support."
            )

        if fp_type is None:
            opts = FingerprintOptions()
        else:
            opts = FingerprintOptions()
            opts.fp_type = fp_type

        return _FingerprintMetric(mols, opts)


class ROCSMetric:
    """
    ROCS-style shape overlay distance metric.

    Computes distances as ``1 - score`` where score is the TanimotoCombo,
    color Tanimoto, or shape Tanimoto from OEShape overlay.
    """

    def __new__(cls, mols, *, color=False, combo=True):
        """
        Construct a ROCSMetric.

        :param mols: List of OEMol molecules with 3D coordinates.
        :param color: Use color Tanimoto only (default False).
        :param combo: Use TanimotoCombo (shape + color, default True).
        :returns: C++ ROCSMetric object.
        :raises ValueError: If ROCSMetric is not available.
        """
        if not _HAVE_ROCS:
            raise ValueError(
                "ROCSMetric is not available. "
                "oecluster was not compiled with OEShape support."
            )

        opts = ROCSOptions()
        opts.color_score = color
        opts.combo_score = combo

        return _ROCSMetric(mols, opts)


class SuperposeMetric:
    """
    Protein superposition distance metric using RMSD.

    Computes pairwise RMSD between protein structures after sequence
    alignment and optimal overlay using OEBio.
    """

    def __new__(cls, items, *, sequence_align=True, only_calpha=True):
        """
        Construct a SuperposeMetric.

        :param items: List of OEDesignUnit objects or OEMolBase molecules.
        :param sequence_align: Align by sequence first (default True).
        :param only_calpha: Use only C-alpha atoms for RMSD (default True).
        :returns: C++ SuperposeMetric object.
        :raises ValueError: If SuperposeMetric is not available.
        """
        if not _HAVE_SUPERPOSE:
            raise ValueError(
                "SuperposeMetric is not available. "
                "oecluster was not compiled with OEBio support."
            )

        opts = SuperposeOptions()
        opts.sequence_align = sequence_align
        opts.only_calpha = only_calpha

        return _SuperposeMetric(items, opts)


class SiteHopperMetric:
    """
    Binding site comparison distance metric using RMSD.

    Extracts binding site protein components from design units and computes
    pairwise RMSD using OEBio.
    """

    def __new__(cls, dus, *, only_calpha=True):
        """
        Construct a SiteHopperMetric.

        :param dus: List of OEDesignUnit objects.
        :param only_calpha: Use only C-alpha atoms for RMSD (default True).
        :returns: C++ SiteHopperMetric object.
        :raises ValueError: If SiteHopperMetric is not available.
        """
        if not _HAVE_SITEHOPPER:
            raise ValueError(
                "SiteHopperMetric is not available. "
                "oecluster was not compiled with OEBio support."
            )

        opts = SiteHopperOptions()
        opts.only_calpha = only_calpha

        return _SiteHopperMetric(dus, opts)
