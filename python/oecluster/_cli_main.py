"""The ``oecluster`` command line.

Commands are flat and self-contained: each takes a distance matrix and does a
whole job. The algorithm and its options arrive as ``--algorithm NAME`` plus
repeatable ``--set key=value``, validated against the registry before any file
is read.
"""
import rich_click as click

from . import __version__, _cli_registry, _cli_render


@click.group()
@click.version_option(__version__, prog_name="oecluster")
def cli():
    """Cluster a precomputed distance matrix."""


@cli.command()
@click.argument("name", required=False)
def algorithms(name):
    """List the clustering algorithms, or show one algorithm's options."""
    registry = _cli_registry.build()
    out = _cli_render.console()
    if name is None:
        rows = []
        for key in sorted(registry):
            entry = registry[key]
            rows.append((key, entry.kind,
                         "yes" if entry.eligible else "no", entry.summary))
        out.print(_cli_render.table(
            ("algorithm", "input", "usable", "summary"), rows))
        return
    if name not in registry:
        raise click.UsageError(f"unknown algorithm {name!r}")
    entry = registry[name]
    if not entry.eligible:
        # Listed rather than hidden: a user who knows the name should learn
        # why it cannot be used, not that it does not exist.
        out.print(f"{name} takes {entry.kind}, not a distance matrix")
        return
    rows = [(option, _cli_registry.option_type(name, option, entry),
             "required", "")
            for option in sorted(entry.required)]
    rows += [(option, _cli_registry.option_type(name, option, entry),
              "optional", repr(entry.optional[option]))
             for option in sorted(entry.optional)]
    out.print(_cli_render.table(("option", "type", "required", "default"),
                                rows, title=name))


def main():
    """Console-script entry point."""
    cli()
