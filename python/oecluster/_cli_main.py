"""The ``oecluster`` command line.

Commands are flat and self-contained: each takes a distance matrix and does a
whole job. The algorithm and its options arrive as ``--algorithm NAME`` plus
repeatable ``--set key=value``, validated against the registry before any file
is read.
"""
import functools
import os

import rich_click as click

from . import ClusteringSpec, __version__, _cli_input, _cli_registry, _cli_render
from ._stability import _thread_count


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


def _translate(function):
    """Turn library exceptions into CLI exits instead of tracebacks.

    Validation errors are the user's invocation, so they exit 2 as Click's own
    parse errors do; everything else is a failure of the run and exits 1.
    """
    @functools.wraps(function)
    def wrapper(*args, **kwargs):
        # --traceback is consumed here, so command bodies never see it.
        show = (kwargs.pop("traceback", False)
                or os.environ.get("OECLUSTER_CLI_TRACEBACK"))
        if show:
            try:
                return function(*args, **kwargs)
            except Exception as error:
                # _cli_input raises `from None` throughout, which does not
                # drop the original exception -- it only hides the context
                # already recorded on this one. The flag is not knowable at
                # those raise sites, so the suppression is lifted here
                # instead; otherwise --traceback shows the frame that
                # renamed the failure rather than the one that failed.
                error.__suppress_context__ = False
                raise
        try:
            return function(*args, **kwargs)
        except (TypeError, ValueError) as error:
            raise click.UsageError(str(error)) from None
        except (_cli_input.InputError, RuntimeError, OSError) as error:
            raise click.ClickException(str(error)) from None
        except MemoryError as error:
            raise click.ClickException(
                f"out of memory: {error}; for consensus, --mmap keeps the "
                "matrix on disk") from None
    return wrapper


def _common(function):
    """Attach the options every working command shares."""
    function = click.option("--quiet", is_flag=True,
                            help="Suppress the terminal summary.")(function)
    function = click.option("--output", default=None, metavar="PATH",
                            help="Write results to a .csv or .json file.")(function)
    function = click.option("--allow-nonmetric", is_flag=True,
                            help="Proceed on a matrix with triangle violations.")(function)
    function = click.option("--threads", default=0, show_default=True,
                            help="Worker threads; 0 auto-detects.")(function)
    function = click.option("--traceback", is_flag=True,
                            help="Show the full traceback on failure.")(function)
    return function


def _spec(algorithm, assignments, registry, threads, nonmetric, *, swept=None):
    """Validate an algorithm and its options into a ClusteringSpec.

    :raises ValueError: For an unknown or ineligible algorithm, or any option
        the registry refuses.
    """
    if algorithm not in registry:
        raise ValueError(
            f"unknown algorithm {algorithm!r}; run 'oecluster algorithms'")
    entry = registry[algorithm]
    if not entry.eligible:
        raise ValueError(
            f"{algorithm} takes {entry.kind}, not a distance matrix")
    options = _cli_registry.resolve(algorithm, assignments, registry,
                                    swept=swept)
    options["num_threads"] = _thread_count(threads)
    if nonmetric and entry.accepts_nonmetric:
        options["allow_nonmetric"] = True
    return ClusteringSpec(algorithm, **options)


def _labels_of(matrix):
    """:returns: The matrix's own labels as text, or None when it has none.

    An unlabelled matrix reports ``[]`` rather than ``None``, so testing for
    None would treat "no labels" as a real empty list. Each label is rendered
    through :func:`_cli_render.text` because a ``.npz`` restores whatever the
    Python API was handed, which need not be a string.
    """
    labels = list(matrix.labels) if matrix.labels is not None else []
    return [_cli_render.text(label) for label in labels] or None


def _ids(matrix):
    """:returns: Labels for the CSV id column, or indices when there are none."""
    return _labels_of(matrix) or list(range(matrix.num_samples))


def _check_destinations(matrix, *destinations):
    """Refuse a destination that would clobber the input or another output.

    With input ``d.npy`` the obvious ``--output d.json`` truncates exactly
    the sidecar the run depends on, and ``--mmap`` pointed at the input
    overwrites the matrix while it is mapped. Every command that writes
    calls this before loading anything.

    :raises ValueError: If a destination is the input or its sidecar, or if
        two destinations resolve to one file.
    """
    protected = {os.path.realpath(matrix): "the input matrix"}
    stem, extension = os.path.splitext(matrix)
    if extension in (".npy", ".bin"):
        protected[os.path.realpath(stem + ".json")] = "the input sidecar"
    seen = {}
    for path in [item for item in destinations if item]:
        resolved = os.path.realpath(path)
        if resolved in protected:
            raise ValueError(f"{path} is {protected[resolved]}; "
                             "choose different files")
        if resolved in seen:
            raise ValueError(f"{path} and {seen[resolved]} are the same "
                             "file; choose different files")
        seen[resolved] = path


@cli.command()
@click.argument("matrix")
@click.option("--algorithm", required=True, help="Roster name.")
@click.option("--set", "assignments", multiple=True, metavar="KEY=VALUE",
              help="Algorithm option; repeatable.")
@_common
@_translate
def cluster(matrix, algorithm, assignments, threads, allow_nonmetric, output,
            quiet):
    """Run one clustering algorithm over a distance matrix."""
    _check_destinations(matrix, output)
    registry = _cli_registry.build()
    spec = _spec(algorithm, assignments, registry, threads, allow_nonmetric)
    loaded = _cli_input.load(matrix)
    result = spec.run(loaded)
    labels = [int(value) for value in result.labels]
    noise = sum(1 for value in labels if value < 0)
    described = {"algorithm": algorithm, "options": dict(spec.options)}
    if not quiet:
        out = _cli_render.console()
        out.print(_cli_render.panel("cluster", described, matrix,
                                    loaded.num_samples))
        out.print(_cli_render.table(
            ("cluster", "size"),
            [(index, len(members))
             for index, members in enumerate(result.clusters)]))
        fraction = noise / loaded.num_samples if loaded.num_samples else 0.0
        out.print(f"clusters={result.num_clusters} noise={noise} "
                  f"({fraction:.1%})")
    if output:
        ids = _ids(loaded)
        named = _labels_of(loaded)
        _cli_render.write_output(output, {
            "header": ["id", "label"],
            "rows": list(zip(ids, labels)),
            "document": {
                "schema_version": _cli_render.SCHEMA_VERSION,
                "command": "cluster",
                "input": {"path": matrix, "num_items": loaded.num_samples,
                          "orientation": str(loaded.is_distance)},
                "spec": described,
                "result": {"labels": labels,
                           # Conditional per the spec: present only when the
                           # matrix carried labels of its own.
                           **({"ids": named} if named else {}),
                           "num_clusters": result.num_clusters,
                           "num_noise": noise},
            }})


def main():
    """Console-script entry point."""
    cli()
