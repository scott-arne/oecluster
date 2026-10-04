"""Rich presentation for the command line.

A plain ``Console()`` is correct here: rich detects a non-TTY and drops colour
and markup by itself, so redirected output and the test runner both see clean
text without a flag.
"""
from rich.console import Console
from rich.table import Table


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
