"""Load a distance matrix for the command line.

Two paths, deliberately different. A ``.npz`` was written by the Python API
and already carries its facts, so it is used as loaded; rebuilding it would
discard the orientation evidence the gate depends on. The raw formats carry
nothing, so they are reconstructed from their JSON sidecar.
"""
import json
import math
import os
import sys
import zipfile

import numpy as np

from . import SymmetricDistanceMatrix, load_distance_matrix


class InputError(Exception):
    """A distance-matrix file that cannot be used, with the reason."""


def _sidecar_path(path):
    stem, _ = os.path.splitext(path)
    return stem + ".json"


def _items_from_pairs(count):
    """:returns: n where count == n(n-1)/2, or None if there is no such n."""
    n = (1 + math.isqrt(1 + 8 * count)) // 2
    return n if n * (n - 1) // 2 == count else None


def _read_sidecar(path, raw_path):
    if not os.path.isfile(path):
        raise InputError(
            f"{os.path.basename(raw_path)} needs its sidecar "
            f"{os.path.basename(path)}, which is missing; it carries the "
            "shape and provenance the raw file does not")
    try:
        with open(path, encoding="utf-8") as handle:
            return json.load(handle)
    except UnicodeDecodeError:
        raise InputError(
            f"{os.path.basename(path)} is not valid UTF-8; a molecule title "
            "with non-UTF-8 bytes is the usual cause") from None
    except OSError as error:
        raise InputError(
            f"cannot read {os.path.basename(path)}: {error}") from None
    except json.JSONDecodeError as error:
        # oepdist's JsonEscape handles only quotes and backslashes, so a
        # molecule title with a newline produces an unparseable sidecar.
        raise InputError(
            f"{os.path.basename(path)} is not valid JSON ({error.msg}); a "
            "molecule title containing a newline is the usual cause") from None


def _values(path):
    """:raises InputError: If the file cannot be read as float64 values."""
    try:
        if path.endswith(".npy"):
            return np.asarray(np.load(path, allow_pickle=False),
                              dtype=np.float64).ravel()
        with open(path, "rb") as handle:
            return np.frombuffer(handle.read(), dtype=np.float64)
    except (OSError, ValueError, EOFError, zipfile.BadZipFile) as error:
        raise InputError(
            f"cannot read {os.path.basename(path)}: {error}") from None


def load(path, *, warn=None):
    """Load ``path`` as a symmetric distance matrix.

    :param path: Input file, ``.npz``, ``.npy`` or ``.bin``.
    :param warn: Callable taking one message, used for unproven orientation.
    :returns: A :class:`SymmetricDistanceMatrix`.
    :raises InputError: For any file that cannot be used, with the reason.
    """
    warn = warn or (lambda message: print(message, file=sys.stderr))
    if not os.path.isfile(path):
        raise InputError(f"no such file: {path}")
    if path.endswith(".csv"):
        raise InputError(
            "csv is not accepted: oepdist writes titles unquoted and values "
            "at 8 significant digits, and records no provenance; re-run "
            "oepdist with a .npy output, or save a .npz from Python")
    if path.endswith(".npz"):
        # TypeError is in the net because the archive's metadata is decoded
        # JSON of any shape: a scalar facts_json reaches dict.update and
        # raises "'int' object is not iterable" from deep inside the library.
        try:
            matrix = load_distance_matrix(path)
        except (OSError, TypeError, ValueError, KeyError, EOFError,
                zipfile.BadZipFile) as error:
            raise InputError(
                f"cannot read {os.path.basename(path)}: {error}") from None
        if not isinstance(matrix, SymmetricDistanceMatrix):
            raise InputError(
                "this is a cross-distance matrix; clustering needs a "
                "symmetric one")
        # from_condensed enforces this, but the .npz path restores state
        # directly and does not. A short list would silently truncate every
        # per-item output to its length.
        count = len(matrix.labels)
        if count and count != matrix.num_samples:
            raise InputError(
                f"{os.path.basename(path)} carries {count} labels for "
                f"{matrix.num_samples} items")
        if matrix.is_distance == "unknown":
            warn(f"{os.path.basename(path)}: orientation unproven, treating "
                 "values as distances")
        return matrix
    if not path.endswith((".npy", ".bin")):
        raise InputError(f"unsupported input format: {path}")

    sidecar = _read_sidecar(_sidecar_path(path), path)
    if not isinstance(sidecar, dict):
        raise InputError(
            f"{os.path.basename(_sidecar_path(path))} is not a JSON object")
    if sidecar.get("mode") != "pdist":
        raise InputError(
            f"sidecar mode is {sidecar.get('mode')!r}, not 'pdist'; a "
            "cross-distance matrix cannot be clustered")
    rows, cols = sidecar.get("n_rows"), sidecar.get("n_cols")
    if rows != cols:
        raise InputError(f"sidecar is {rows}x{cols}, not square")
    values = _values(path)
    items = _items_from_pairs(values.size)
    if items is None:
        raise InputError(
            f"{values.size} values is not a condensed symmetric matrix")
    if items != rows:
        raise InputError(
            f"{values.size} values describe {items} items, but the sidecar "
            f"says {rows}")
    params = sidecar.get("params") or {}
    if not isinstance(params, dict):
        raise InputError(
            f"{os.path.basename(_sidecar_path(path))} has a non-object "
            "'params' field")
    similarity = params.get("similarity")
    if similarity is True:
        raise InputError(
            "this file holds similarities, not distances; re-run oepdist "
            "without --sim, or convert it with the Python API")
    if similarity is None:
        warn(f"{os.path.basename(path)}: orientation unproven (the sidecar "
             "records no similarity flag), treating values as distances")
    labels = sidecar.get("row_labels") or None
    # from_condensed runs the library's own validation and metric probe; its
    # refusals are the user's problem with this file, so they arrive as
    # InputError like every other unusable input rather than as a bare
    # ValueError that the caller would map to a usage error.
    try:
        return SymmetricDistanceMatrix.from_condensed(
            values, labels=labels,
            comparison_name=sidecar.get("comparison") or "precomputed",
            params=params or None)
    except (TypeError, ValueError) as error:
        raise InputError(
            f"{os.path.basename(path)} is not a usable distance matrix: "
            f"{error}") from None
