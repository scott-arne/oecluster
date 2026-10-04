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
from typing import NoReturn

import numpy as np

from . import SymmetricDistanceMatrix, load_distance_matrix


class InputError(Exception):
    """A distance-matrix file that cannot be used, with the reason."""


def _unreadable(path, error) -> NoReturn:
    """Report a failed library call as an unusable file.

    The catch around each library call is inverted rather than enumerated.
    Four rounds of review each found a new escaping type -- ``EOFError``,
    ``BadZipFile``, ``UnicodeDecodeError``, ``TypeError``, ``IndexError``
    from a 0-d ``condensed``, ``OverflowError`` from a negative
    ``num_samples`` reaching ``size_t`` -- because the metadata is decoded
    JSON of any shape and the arrays carry any dtype and rank. A list of
    types cannot win that race, so every failure of a library call on
    user-supplied bytes is an unusable file.

    :param path: The file being read, named in the message.
    :param error: The exception the library raised.
    :raises InputError: Always.
    """
    raise InputError(
        f"cannot read {os.path.basename(path)}: {error}") from None


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
    except json.JSONDecodeError as error:
        # oepdist's JsonEscape handles only quotes and backslashes, so a
        # molecule title with a newline produces an unparseable sidecar.
        raise InputError(
            f"{os.path.basename(path)} is not valid JSON ({error.msg}); a "
            "molecule title containing a newline is the usual cause") from None
    except MemoryError:
        raise
    except Exception as error:  # noqa: BLE001
        _unreadable(path, error)


def _values(path):
    """Read ``path`` as a flat float64 array.

    :returns: The file's values as one-dimensional float64.
    :raises InputError: If the file cannot be read, or holds anything but
        real numbers.
    """
    try:
        if path.endswith(".npy"):
            array = np.load(path, allow_pickle=False)
        else:
            with open(path, "rb") as handle:
                array = np.frombuffer(handle.read(), dtype=np.float64)
    except MemoryError:
        raise
    except Exception as error:  # noqa: BLE001
        _unreadable(path, error)
    # np.load reads the bytes, not the name: a .npz archive renamed .npy
    # comes back as a lazy NpzFile, and every check below is an ndarray's.
    if not isinstance(array, np.ndarray):
        array.close()
        raise InputError(
            f"{os.path.basename(path)} is a .npz archive under a .npy name; "
            "rename it to .npz so its own metadata is read")
    # Checked before the cast, which is where the evidence is lost:
    # float64(complex) keeps the real part behind a ComplexWarning the
    # caller may have filtered, float64(datetime64) yields epoch seconds,
    # and float64(bool) yields 0/1 -- an adjacency matrix read as distances
    # is a similarity, the one error this module exists to refuse.
    # from_condensed refuses complex for the same reason and never sees the
    # dtype, because the cast has already happened by the time it is called.
    if array.dtype.kind not in "fiu":
        raise InputError(
            f"{os.path.basename(path)} holds {array.dtype} values; a "
            "distance matrix must hold real numbers")
    # Left outside the guard above: on an array already known to be
    # real-numeric the only failure left is MemoryError, which must reach
    # the caller as itself.
    return np.asarray(array, dtype=np.float64).ravel()


def load(path, *, warn=None):
    """Load ``path`` as a symmetric distance matrix.

    :param path: Input file, ``.npz``, ``.npy`` or ``.bin``.
    :param warn: Callable taking one message, used for unproven orientation.
    :returns: A :class:`SymmetricDistanceMatrix`.
    :raises InputError: For any file that cannot be used, with the reason.
    """
    warn = warn or (lambda message: print(message, file=sys.stderr))
    # click hands back a pathlib.Path when the argument is declared with
    # click.Path(path_type=Path), and every format check below is a string
    # operation that would raise AttributeError on one. Decoded rather than
    # just fspath'd, because fspath passes bytes straight through and the
    # checks would fail the same way on those.
    path = os.fsdecode(path)
    if not os.path.isfile(path):
        raise InputError(f"no such file: {path}")
    if path.endswith(".csv"):
        raise InputError(
            "csv is not accepted: oepdist writes titles unquoted and values "
            "at 8 significant digits, and records no provenance; re-run "
            "oepdist with a .npy output, or save a .npz from Python")
    if path.endswith(".npz"):
        try:
            matrix = load_distance_matrix(path)
        except MemoryError:
            raise
        except Exception as error:  # noqa: BLE001
            _unreadable(path, error)
        if not isinstance(matrix, SymmetricDistanceMatrix):
            raise InputError(
                "this is a cross-distance matrix; clustering needs a "
                "symmetric one")
        # Unguarded: these three read state the call above already built,
        # and 34 corrupted archives produced none that loads and then fails
        # a read. A guard with no reachable catch is surface, not safety.
        labelled, items = len(matrix.labels), matrix.num_samples
        orientation = matrix.is_distance
        # from_condensed enforces this, but the .npz path restores state
        # directly and does not. A short list would silently truncate every
        # per-item output to its length.
        if labelled and labelled != items:
            raise InputError(
                f"{os.path.basename(path)} carries {labelled} labels for "
                f"{items} items")
        # Tested against the two booleans rather than for "unknown": facts
        # are arbitrary JSON here, the library's gate refuses only
        # `is_distance is False`, and a corrupted fact of 0 is neither a
        # refusal nor a proof. Anything that is not a boolean is unproven.
        if orientation is not True and orientation is not False:
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
    elif similarity is not False:
        # Only the JSON booleans mean anything here. Tested by identity, so
        # 1 and "true" fall through both branches above: neither refused nor
        # warned about, they would be clustered as proven distances on the
        # strength of a flag that in fact says the opposite.
        raise InputError(
            f"sidecar records similarity={similarity!r}, which is neither "
            "true nor false; orientation cannot be read from it")
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
    except MemoryError:
        raise
    except Exception as error:  # noqa: BLE001
        raise InputError(
            f"{os.path.basename(path)} is not a usable distance matrix: "
            f"{error}") from None
