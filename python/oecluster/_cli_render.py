"""Rich presentation for the command line.

A plain ``Console()`` is correct here: rich detects a non-TTY and drops colour
and markup by itself, so redirected output and the test runner both see clean
text without a flag.
"""
import csv
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


def panel(command, spec, path, num_items):
    """Build the header panel naming what ran.

    :param spec: Either ``{"algorithm": name, "options": {...}}`` or
        ``{"members": [{"algorithm": ..., "options": {...}}, ...]}``.
    :returns: A :class:`rich.panel.Panel`.
    """
    if "members" in spec:
        what = ", ".join(member["algorithm"] for member in spec["members"])
    else:
        # num_threads is routed by --threads rather than chosen per run, so
        # echoing it back as an option would be noise.
        shown = " ".join(f"{key}={value}"
                         for key, value in sorted(spec["options"].items())
                         if key != "num_threads")
        what = f"{spec['algorithm']} {shown}".strip()
    return Panel(f"{what}\n{path}  ({num_items} items)", title=command)


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
    if isinstance(value, float):
        if math.isnan(value):
            return ""
        if math.isinf(value):
            return "inf" if value > 0 else "-inf"
    # csv.writer falls back to str() for anything it is not given as text,
    # which turns a np.bytes_ label into b'x'.
    if isinstance(value, bytes):
        return text(value)
    return "" if value is None else value


def write_output(path, payload):
    """Write results as CSV or JSON, chosen by the extension.

    The document is run through :func:`encode` here rather than at each call
    site, so a caller cannot leave a NaN in it by forgetting to; consensus
    and stability produce NaN routinely, and the file is read by tools that
    are stricter than Python's parser. ``allow_nan=False`` is the tripwire
    behind that: a value the walk somehow did not reach fails the write
    loudly instead of becoming a file no strict parser will load.

    :param path: Destination; must end in ``.csv`` or ``.json``.
    :param payload: Mapping with ``header`` and ``rows`` for CSV and
        ``document`` for JSON.
    :raises ValueError: For any other extension, or for a non-finite value
        :func:`_encoded` did not reach.
    """
    check_output(path)
    if os.path.splitext(path)[1].lower() == ".csv":
        with open(path, "w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle)
            writer.writerow(payload["header"])
            for row in payload["rows"]:
                writer.writerow([_csv_cell(value) for value in row])
    else:
        with open(path, "w", encoding="utf-8") as handle:
            json.dump(_encoded(payload["document"]), handle, indent=2,
                      allow_nan=False)
            handle.write("\n")
