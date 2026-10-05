"""Rich presentation for the command line.

A plain ``Console()`` is correct here: rich detects a non-TTY and drops colour
and markup by itself, so redirected output and the test runner both see clean
text without a flag.
"""
import csv
import io
import json
import math
import os

from rich.console import Console
from rich.panel import Panel
from rich.table import Table

SCHEMA_VERSION = 1

#: Extensions :func:`write_output` can write. Public so the command layer can
#: refuse a destination before doing the work rather than after it, and so the
#: two cannot drift into disagreeing about what is writable.
OUTPUT_SUFFIXES = (".csv", ".json")


def console():
    """:returns: The console every command prints through."""
    return Console()


def table(columns, rows, *, title=None):
    """Build a table from a ``columns``/``to_table()`` pair.

    `ParameterSelection`, `ClusterStability` and `ConsensusResult` all expose
    that pair, so one renderer serves three commands.

    :param columns: Column names.
    :param rows: Iterable of row tuples in ``columns`` order.
    :param title: Optional table title.
    :returns: A populated :class:`rich.table.Table`.
    """
    built = Table(*[str(name) for name in columns], title=title)
    for row in rows:
        built.add_row(*["" if value is None else str(value) for value in row])
    return built


def text(value):
    """Render one label as text.

    A ``.npz`` restores whatever labels the Python API was handed, so a label
    is not necessarily a string: numpy brings a bytes label back as
    ``np.bytes_``, and ``str`` on that writes ``b'x'`` -- an identity the user
    never wrote. The sidecar path refuses non-strings because oepdist always
    writes titles as JSON strings; the ``.npz`` path accepts them because they
    are the user's own API input, so the rendering happens here instead.

    Undecodable bytes cannot reach here from a loaded matrix, because
    ``_cli_input`` refuses them outright. The fallback is
    ``backslashreplace`` rather than ``replace`` so that if one ever did,
    two labels would still not collapse onto a single id: ``replace`` maps
    every undecodable byte onto the same character, which is the exact
    failure the loader's refusal exists to prevent.

    :param value: A label of any type.
    :returns: Its text form.
    """
    if isinstance(value, bytes):
        return value.decode("utf-8", "backslashreplace")
    return str(value)


def _identical(first, second):
    """:returns: True if two labels are provably one and the same label.

    Conservative by construction: anything this cannot settle counts as
    "not provably the same", which refuses the run. A refusal naming the
    two items beats an id column that silently merges them.

    ``==`` alone is not enough. A numpy array answers it with an
    elementwise array rather than a bool, and ``bool()`` on that raises for
    anything but a single element, so the result is reduced with ``all()``
    whenever the comparison offers one.
    """
    if first is second:
        return True
    try:
        equal = first == second
        reduce = getattr(equal, "all", None)
        return bool(equal if reduce is None else reduce())
    except Exception:  # noqa: BLE001
        return False


def check_ids(labels):
    """Render labels as the exported id column, refusing a lossy rendering.

    The id column has one job -- to say which item a row is about -- so the
    rendering has to keep distinct labels distinct. ``str`` does not:
    numpy prints ``array([1.000000001])`` as ``[1.]``, so three different
    labels exported one id and the result could no longer be joined back to
    what was clustered. The rule is stated over the class rather than over
    a type because the same defect arrived first through ``bytes`` and
    would arrive next through a third type.

    Duplicate labels stay legitimate: two items both called ``mol1`` share
    an id on purpose and still export two rows. What is refused is two
    *different* labels arriving at one id.

    An id that cannot be encoded as UTF-8 is refused here too. It is the
    same question -- whether a label can be an id at all -- and asking it
    here means the answer arrives before the run rather than from inside
    the writer.

    :param labels: The labels, in item order.
    :returns: The rendered id column, in the same order.
    :raises ValueError: If two distinct labels render alike, or if a
        rendered id cannot be encoded as UTF-8.
    """
    rendered = [text(label) for label in labels]
    first_seen = {}
    for index, name in enumerate(rendered):
        try:
            name.encode("utf-8")
        except UnicodeEncodeError:
            raise ValueError(
                f"the label of item {index} cannot be encoded as UTF-8, so "
                "it cannot be written as an id; a lone surrogate is the "
                "usual cause") from None
        earlier = first_seen.setdefault(name, index)
        if earlier != index and not _identical(labels[earlier], labels[index]):
            raise ValueError(
                f"items {earlier} and {index} have different labels that "
                f"both render as {name!r}, so the id column cannot tell "
                "them apart; relabel the matrix through the Python API")
    return rendered


def _named(spec):
    """:returns: One algorithm and its options as a single line of text."""
    # num_threads is routed by --threads rather than chosen per run, so
    # echoing it back as an option would be noise.
    shown = " ".join(f"{key}={value}"
                     for key, value in sorted(spec["options"].items())
                     if key != "num_threads")
    return f"{spec['algorithm']} {shown}".strip()


def panel(command, spec, path, num_items):
    """Build the header panel naming what ran.

    Ensemble members carry their options for the same reason a single spec
    does, and more urgently: names alone rendered two butina members at
    different thresholds as "butina, butina", erasing the one thing that
    distinguished them.

    :param spec: Either ``{"algorithm": name, "options": {...}}`` or
        ``{"members": [{"algorithm": ..., "options": {...}}, ...]}``.
    :returns: A :class:`rich.panel.Panel`.
    """
    if "members" in spec:
        what = ", ".join(_named(member) for member in spec["members"])
    else:
        what = _named(spec)
    return Panel(f"{what}\n{path}  ({num_items} items)", title=command)


def orientation(is_distance):
    """Spell a matrix's orientation fact for the JSON document.

    ``str()`` on the fact cannot do this job: the library reports a proven
    orientation as a Python ``bool`` and an unproven one as the string
    ``"unknown"``, so the repr put ``"True"`` and ``"unknown"`` in one
    ``schema_version`` 1 field -- two casing conventions a reader would have
    to know about. The field stays a string rather than becoming a JSON
    boolean beside a ``"unknown"`` string, because a union-typed field is
    the harder thing for a typed consumer to read.

    Anything that is not one of the two booleans is ``"unknown"``, matching
    the loader: ``facts_json`` is arbitrary JSON, so a corrupted fact of 0
    is unproven rather than false.

    :param is_distance: ``SymmetricDistanceMatrix.is_distance``.
    :returns: ``"true"``, ``"false"`` or ``"unknown"``.
    """
    if is_distance is True:
        return "true"
    if is_distance is False:
        return "false"
    return "unknown"


def encode(value):
    """Encode one statistic for JSON.

    NaN becomes ``null`` and the infinities become sentinel strings:
    collapsing them together would erase a real distinction, since NaN means
    undefined while an infinity is a defined extreme some criteria produce.

    :param value: Any statistic.
    :returns: A JSON-safe value.
    """
    if isinstance(value, float):
        if math.isnan(value):
            return None
        if math.isinf(value):
            return "inf" if value > 0 else "-inf"
    return value


def check_output(path):
    """Refuse a destination :func:`write_output` could not write.

    Separated from the write so a command can call it before the run: the
    extension was checked where the file is opened, which is after the
    clustering it was meant to save.

    :param path: Destination path.
    :raises ValueError: For any extension outside :data:`OUTPUT_SUFFIXES`.
    """
    if os.path.splitext(path)[1].lower() not in OUTPUT_SUFFIXES:
        # Spelled from the tuple so the message cannot outlive the list it
        # describes.
        raise ValueError(f"--output must end in {' or '.join(OUTPUT_SUFFIXES)}"
                         f", got {path!r}")


def _encoded(value):
    """:returns: ``value`` with every float encoded for strict JSON.

    The walk has to happen before :func:`json.dump`, not inside it: its
    ``default`` hook fires only for types json cannot serialize, and a NaN
    float is not one of them -- json writes it as the bare ``NaN`` literal
    that strict parsers reject.
    """
    if isinstance(value, dict):
        return {key: _encoded(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_encoded(item) for item in value]
    return encode(value)


def _csv_cell(value):
    # Spelled here rather than at the call site so one convention covers
    # every column: csv.writer falls back to str(), which writes Python's
    # repr "True", while select-parameter's own `winner` column and the
    # JSON document both spell a boolean lowercase. One select-parameter
    # row was carrying both conventions at once.
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, float):
        if math.isnan(value):
            return ""
        if math.isinf(value):
            return "inf" if value > 0 else "-inf"
    return "" if value is None else value


def write_output(path, payload):
    """Write results as CSV or JSON, chosen by the extension.

    The document is run through :func:`encode` here rather than at each call
    site, so a caller cannot leave a NaN in it by forgetting to; consensus
    and stability produce NaN routinely, and the file is read by tools that
    are stricter than Python's parser. ``allow_nan=False`` is the tripwire
    behind that: a value the walk somehow did not reach fails the write
    loudly instead of becoming a file no strict parser will load.

    The file is rendered and encoded in full before it is opened, so a
    failure part way through leaves no file rather than a truncated one: a
    header-only CSV on disk reads as a run that succeeded and found
    nothing, which is worse than no file at all. The bytes are written
    through a binary handle because a text handle encodes during the write,
    which is on the far side of the open.

    :param path: Destination; must end in ``.csv`` or ``.json``.
    :param payload: Mapping with ``header`` and ``rows`` for CSV and
        ``document`` for JSON.
    :raises ValueError: For any other extension, or for a non-finite value
        :func:`_encoded` did not reach.
    """
    check_output(path)
    if os.path.splitext(path)[1].lower() == ".csv":
        buffer = io.StringIO()
        writer = csv.writer(buffer)
        writer.writerow(payload["header"])
        for row in payload["rows"]:
            writer.writerow([_csv_cell(value) for value in row])
        body = buffer.getvalue()
    else:
        body = json.dumps(_encoded(payload["document"]), indent=2,
                          allow_nan=False) + "\n"
    data = body.encode("utf-8")
    with open(path, "wb") as handle:
        handle.write(data)
