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

from . import (
    SparseStorage,
    SymmetricDistanceMatrix,
    _cli_render,
    load_distance_matrix,
)


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


def _refuse_negative(matrix, path):
    """Refuse a matrix holding a negative distance.

    Shared by both paths so they cannot drift again. ``from_condensed``
    already refuses this on the raw path; a ``.npz`` is used as loaded and
    so never saw the check, and nothing downstream covers the difference --
    the clustering gate re-scans for non-finite values, not for negatives.
    A negative distance reaches an archive through the public API, not only
    by editing one: ``condensed`` is writeable and ``to_file`` writes the
    mutated values under the facts the matrix was built with.

    :param matrix: The loaded matrix.
    :param path: The file it came from, named in the message.
    :raises InputError: If any stored distance is below zero.
    """
    if isinstance(matrix.storage, SparseStorage):
        # The entry list, not ``condensed``: that property answers the
        # question by densifying a sparse matrix into n(n-1)/2 floats, and
        # the omitted pairs it invents are zeros that cannot be negative.
        values = np.array([entry[2] for entry in matrix.storage._entries()],
                          dtype=np.float64)
    else:
        values = matrix.condensed
    # ``min`` propagates NaN and ``nan < 0`` is False, so a non-finite value
    # cannot switch this check off; refusing those is the gate's job.
    if values.size and float(values.min()) < 0.0:
        raise InputError(
            f"{os.path.basename(path)} holds negative values, so it is not "
            "a distance matrix; a distance cannot be below zero")


def _refuse_colliding_labels(matrix, path):
    """Refuse labels the id column could not tell apart.

    Delegated to the renderer rather than reimplemented, so what is checked
    here is literally what will be exported: :func:`_cli_render.check_ids`
    is the one function that turns a label into an id, and the command
    layer calls it again on the way out. Checking it here is what makes the
    refusal cheap -- it lands before the clustering rather than after it.

    :param matrix: The loaded matrix.
    :param path: The file it came from, named in the message.
    :raises InputError: If two distinct labels render as one id, or if an
        id cannot be encoded as UTF-8. Both are facts about this file, so
        they arrive as an unusable input like every other.
    """
    labels = matrix.labels if matrix.labels is not None else []
    try:
        _cli_render.check_ids(labels)
    except ValueError as error:
        raise InputError(f"{os.path.basename(path)}: {error}") from None


def _refuse_unrenderable_labels(matrix, path):
    """Refuse labels that have no distinct text form.

    A ``.npz`` restores whatever labels the Python API was handed, and numpy
    brings a bytes label back as ``np.bytes_``. Bytes that are not UTF-8
    have no faithful text form, and the lossy decodings are worse than the
    refusal: ``errors="replace"`` maps every undecodable byte onto the one
    replacement character, so two items export the same id and the result
    can no longer be joined back to what the user clustered. The raw path
    already refuses a non-UTF-8 title in its sidecar (:func:`_read_sidecar`)
    for the same reason; refusing here is what keeps the two paths from
    disagreeing about the same bad title.

    :param matrix: The loaded matrix.
    :param path: The file it came from, named in the message.
    :raises InputError: If a bytes label is not valid UTF-8.
    """
    labels = matrix.labels if matrix.labels is not None else []
    for index, label in enumerate(labels):
        if not isinstance(label, bytes):
            continue
        try:
            label.decode("utf-8")
        except UnicodeDecodeError:
            raise InputError(
                f"{os.path.basename(path)} has a label at index {index} "
                "whose bytes are not valid UTF-8, so it has no text form "
                "that stays distinct from its neighbours; relabel the "
                "matrix through the Python API") from None


def _label_list(labels, sidecar):
    """Validate a sidecar's ``row_labels`` field.

    ``from_condensed`` runs ``list()`` over whatever it is handed, so a
    string becomes one label per character and an object becomes its keys.
    With a matching entry count both pass every length check and the labels
    reach the output as identities the user never wrote. A non-string entry
    is refused rather than coerced for the same reason: ``str(None)`` is an
    invented identity, and oepdist writes every title as a JSON string.

    :param labels: The raw ``row_labels`` value, or None if absent.
    :param sidecar: Sidecar path, named in the message.
    :returns: The labels, or None if there are none to use.
    :raises InputError: If the field is present but not a list of strings.
    """
    if labels is None:
        return None
    if not isinstance(labels, list):
        raise InputError(
            f"{os.path.basename(sidecar)} has a 'row_labels' field that is "
            f"a {type(labels).__name__}, not an array of strings")
    for index, label in enumerate(labels):
        if not isinstance(label, str):
            raise InputError(
                f"{os.path.basename(sidecar)} has a 'row_labels' entry at "
                f"index {index} that is not a string: {label!r}")
    return labels or None


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


def _values(path, suffix):
    """Read ``path`` as a flat float64 array.

    :param path: The raw input file.
    :param suffix: Its lowercased extension, ``.npy`` or ``.bin``.
    :returns: The file's values as one-dimensional float64.
    :raises InputError: If the file cannot be read, or holds anything but
        real numbers.
    """
    try:
        if suffix == ".npy":
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
    # Lowercased once and dispatched on throughout: oepdist lowercases the
    # extension before choosing a writer (tools/OutputWriter.cpp:19-24), so
    # `-o out.NPY` writes a real .NPY that a case-sensitive check refuses as
    # an unsupported format.
    suffix = os.path.splitext(path)[1].lower()
    if suffix == ".csv":
        raise InputError(
            "csv is not accepted: oepdist writes titles unquoted and values "
            "at 8 significant digits, and records no provenance; re-run "
            "oepdist with a .npy output, or save a .npz from Python")
    if suffix == ".npz":
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
        # Refused here rather than left to the library's clustering gate.
        # The gate raises ValueError, which the command layer maps to a
        # usage error -- exit 2, with Click's "try --help" advice, for a
        # fact recorded in the file that no change to the invocation can
        # answer. The sidecar path already refuses the byte-identical
        # evidence as an unusable input, so the same fact now gets the same
        # exit code whichever file carries it. The gate stays as the
        # library's own backstop for callers that do not come through here.
        if orientation is False:
            raise InputError(
                f"{os.path.basename(path)} holds similarities, not "
                "distances: its stored facts record is_distance false. "
                "Rebuild it from distances, or convert it with the Python "
                "API")
        # Tested by identity against True rather than for "unknown": facts
        # are arbitrary JSON here, so a corrupted fact of 0 -- falsy, and
        # plainly not a proven distance -- is neither a refusal nor a
        # proof. Anything that is not True is unproven.
        if orientation is not True:
            warn(f"{os.path.basename(path)}: orientation unproven, treating "
                 "values as distances")
        _refuse_unrenderable_labels(matrix, path)
        _refuse_colliding_labels(matrix, path)
        _refuse_negative(matrix, path)
        return matrix
    if suffix not in (".npy", ".bin"):
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
    values = _values(path, suffix)
    items = _items_from_pairs(values.size)
    if items is None:
        raise InputError(
            f"{values.size} values is not a condensed symmetric matrix")
    if items != rows:
        raise InputError(
            f"{values.size} values describe {items} items, but the sidecar "
            f"says {rows}")
    # Typed before the default is applied, not after: `or {}` turned every
    # falsy non-object -- [], 0, "" and False -- into "no params", so a
    # producer fault read as a sidecar that simply recorded nothing.
    params = sidecar.get("params")
    if params is None:
        params = {}
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
    labels = _label_list(sidecar.get("row_labels"), _sidecar_path(path))
    # from_condensed runs the library's own validation and metric probe; its
    # refusals are the user's problem with this file, so they arrive as
    # InputError like every other unusable input rather than as a bare
    # ValueError that the caller would map to a usage error.
    try:
        matrix = SymmetricDistanceMatrix.from_condensed(
            values, labels=labels,
            comparison_name=sidecar.get("comparison") or "precomputed",
            params=params or None)
    except MemoryError:
        raise
    except Exception as error:  # noqa: BLE001
        raise InputError(
            f"{os.path.basename(path)} is not a usable distance matrix: "
            f"{error}") from None
    # Run on this path too: _label_list proves the labels are strings, not
    # that they can be written. A JSON "\udcff" escape decodes to a lone
    # surrogate, which is a perfectly ordinary str until the writer tries
    # to encode it.
    _refuse_colliding_labels(matrix, path)
    _refuse_negative(matrix, path)
    return matrix
