"""The ``oecluster`` command line.

Commands are flat and self-contained: each takes a distance matrix and does a
whole job. The algorithm and its options arrive as ``--algorithm NAME`` plus
repeatable ``--set key=value``, validated against the registry before any file
is read.
"""
import functools
import os

import rich_click as click

import oecluster

from . import (
    _SIZE_T_MAX,
    ClusteringSpec,
    __version__,
    _cli_input,
    _cli_registry,
    _cli_render,
)
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
        # The variable is read against an off-list rather than for truth:
        # every non-empty value is truthy, so OECLUSTER_CLI_TRACEBACK=0
        # turned tracebacks on.
        show = (kwargs.pop("traceback", False)
                or os.environ.get("OECLUSTER_CLI_TRACEBACK", "")
                not in ("", "0"))
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
        except (click.exceptions.Exit, click.Abort):
            # Click's own control flow, and both subclass RuntimeError: the
            # arm below would turn ctx.exit(0) into exit 1 with the message
            # "0", and an aborted confirmation into a failed run.
            raise
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


def _entry(algorithm, registry):
    """Look up an algorithm that is allowed to run on a distance matrix.

    Both questions are asked here rather than at each call site, and in this
    order. ``select-parameter`` asked only the first before checking its
    swept name, so ``--algorithm murcko`` was answered "murcko has no option
    'threshold' to sweep" where ``cluster`` says "murcko takes mols, not a
    distance matrix" -- and only because ``bitbirch`` happens to have a
    ``threshold``, so the wrong message was not even consistent.

    :param algorithm: Roster name as the user spelled it.
    :param registry: Mapping from :func:`_cli_registry.build`.
    :returns: The :class:`_cli_registry.Entry` for ``algorithm``.
    :raises ValueError: For an unknown name, or one whose input is not a
        distance matrix.
    """
    if algorithm not in registry:
        raise ValueError(
            f"unknown algorithm {algorithm!r}; run 'oecluster algorithms'")
    entry = registry[algorithm]
    if not entry.eligible:
        raise ValueError(
            f"{algorithm} takes {entry.kind}, not a distance matrix")
    return entry


def _spec(algorithm, assignments, registry, threads, nonmetric, *, swept=None):
    """Validate an algorithm and its options into a ClusteringSpec.

    :raises ValueError: For an unknown or ineligible algorithm, or any option
        the registry refuses.
    """
    entry = _entry(algorithm, registry)
    options = _cli_registry.resolve(algorithm, assignments, registry,
                                    swept=swept)
    threads = _thread_count(threads)
    # _thread_count bounds below but not above, and only some entry points
    # bound it themselves (k_medoids does, butina does not), so an oversized
    # value reached a native size_t setter and escaped as an OverflowError
    # traceback after the matrix had already been loaded. The bound matches
    # the library's own, so it refuses nothing the parameter can hold. The
    # --set options need no equivalent: the registry caps their magnitude
    # and the library range-checks each one with a message naming it.
    if threads > _SIZE_T_MAX:
        raise ValueError(f"--threads exceeds size_t maximum, got {threads}")
    options["num_threads"] = threads
    if nonmetric and entry.accepts_nonmetric:
        options["allow_nonmetric"] = True
    return ClusteringSpec(algorithm, **options)


def _labels_of(matrix):
    """:returns: The matrix's own labels as text, or None when it has none.

    An unlabelled matrix reports ``[]`` rather than ``None``, so testing for
    None would treat "no labels" as a real empty list.

    :func:`_cli_render.check_ids` both renders the labels and refuses a
    rendering that two distinct labels would share. ``_cli_input`` calls it
    too, which is what makes the refusal cheap; calling it again here is
    what makes it airtight, because the function that produces the exported
    ids is then the function that checks them and no third party can
    introduce a collision in between.
    """
    labels = list(matrix.labels) if matrix.labels is not None else []
    return _cli_render.check_ids(labels) or None


def _ids(matrix):
    """:returns: Labels for the CSV id column, or indices when there are none."""
    return _labels_of(matrix) or list(range(matrix.num_samples))


def _same_file(first, second):
    """:returns: True if the two paths name one file.

    ``realpath`` resolves symlinks but then compares bytes, so on a
    case-folding filesystem it reports ``d.json`` and ``D.JSON`` as two
    files where there is only one. ``samefile`` asks the filesystem, which
    is the only authority on that, but it needs both paths to exist; a
    destination not yet written can only be compared by name, which is what
    the first test covers.
    """
    if os.path.realpath(first) == os.path.realpath(second):
        return True
    try:
        return os.path.samefile(first, second)
    except OSError:
        return False


def _check_destinations(matrix, output, *destinations):
    """Refuse a destination that would clobber the input or another output.

    With input ``d.npy`` the obvious ``--output d.json`` truncates exactly
    the sidecar the run depends on, and ``--mmap`` pointed at the input
    overwrites the matrix while it is mapped. Every command that writes
    calls this before loading anything.

    Destinations are tested against None, not for truth. ``--output ""``
    is a shell variable that did not expand, and read as falsy it meant "no
    output requested": the run did all its work, exited 0 and wrote nothing,
    silently under ``--quiet``.

    :param matrix: The input path.
    :param output: The ``--output`` path, or None if it was not supplied.
        Its extension is checked against the writer's own list as well, so
        a destination the writer could not write is refused before the run
        rather than after it.
    :param destinations: Further destinations, such as ``--mmap``. They are
        checked for collisions only: nothing writes them through
        :func:`_cli_render.write_output`, so they carry no extension
        contract.
    :raises ValueError: If a destination is empty, is the input or its
        sidecar, if two destinations name one file, or if ``output`` has an
        extension :func:`_cli_render.write_output` cannot write.
    """
    protected = [(matrix, "the input matrix")]
    stem, extension = os.path.splitext(matrix)
    # Lowercased to match _cli_input, which dispatches on the lowercased
    # suffix because `oepdist -o out.NPY` writes a real .NPY. Compared
    # case-sensitively, U.NPY looked like a format that has no sidecar, so
    # the run overwrote the one file that makes the input readable.
    if extension.lower() in (".npy", ".bin"):
        protected.append((stem + ".json", "the input sidecar"))
    seen = []
    for path in [item for item in (output, *destinations) if item is not None]:
        if not path:
            raise ValueError(
                "a destination path is empty; an unset shell variable is "
                "the usual cause")
        for other, what in protected:
            if _same_file(path, other):
                raise ValueError(f"{path} is {what}; "
                                 "choose different files")
        for other in seen:
            if _same_file(path, other):
                raise ValueError(f"{path} and {other} are the same "
                                 "file; choose different files")
        seen.append(path)
    # Last, so a destination that is both unwritable and a collision is
    # reported as the collision: that is the message naming the file at risk.
    if output is not None:
        _cli_render.check_output(output)


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
    if output is not None:
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


@cli.command("select-parameter")
@click.argument("matrix")
@click.option("--algorithm", required=True, help="Roster name.")
@click.option("--parameter", required=True, help="Option to sweep.")
@click.option("--values", required=True,
              help="Comma-separated values for the swept option.")
@click.option("--set", "assignments", multiple=True, metavar="KEY=VALUE",
              help="Fixed algorithm option; repeatable.")
@click.option("--criterion", default=None,
              help="Validity index to rank by; the scorer picks one if unset.")
@click.option("--max-noise-fraction", type=float, default=None,
              help="Reject a value whose noise fraction exceeds this.")
@click.option("--min-clusters", type=int, default=None,
              help="Reject a value that yields fewer clusters than this.")
@click.option("--max-clusters", type=int, default=None,
              help="Reject a value that yields more clusters than this.")
@_common
@_translate
def select_parameter(matrix, algorithm, parameter, values, assignments,
                     criterion, max_noise_fraction, min_clusters,
                     max_clusters, threads, allow_nonmetric, output, quiet):
    """Sweep one option over several values and pick a winner."""
    _check_destinations(matrix, output)
    registry = _cli_registry.build()
    entry = _entry(algorithm, registry)
    # Check the swept name BEFORE building the spec. Reversing these reports
    # "butina requires --set threshold=…" for a misspelled --parameter,
    # because the swept name excused the real option from requiredness.
    if parameter not in entry.known():
        raise ValueError(f"{algorithm} has no option {parameter!r} to sweep")
    # The swept option satisfies requiredness: select_parameter supplies it
    # per value, so demanding it in --set would reject the normal form.
    spec = _spec(algorithm, assignments, registry, threads, allow_nonmetric,
                 swept=parameter)
    # Coerced through the same function --set uses: the grid does not pass
    # through _cli_registry.resolve(), so without this it reached the
    # library as untyped text -- a string threshold, an infinity, or a
    # mistyped exponent int() expands into a billion digits.
    swept = [_cli_registry.coerce(algorithm, parameter, item.strip(), entry)
             for item in values.split(",") if item.strip()]
    if not swept:
        raise ValueError("--values needs at least one value")
    loaded = _cli_input.load(matrix)
    # cluster_report runs its own metric gate, and select_parameter forwards
    # report_options to it unchanged -- so the flag has to reach both places
    # or clustering succeeds and scoring then fails.
    report_options = {"allow_nonmetric": True} if allow_nonmetric else None
    selection = oecluster.select_parameter(
        spec, loaded, parameter, swept, criterion=criterion,
        max_noise_fraction=max_noise_fraction, min_clusters=min_clusters,
        max_clusters=max_clusters, report_options=report_options)
    rows = selection.to_table()
    described = {"algorithm": algorithm, "options": dict(spec.options)}
    if not quiet:
        out = _cli_render.console()
        out.print(_cli_render.panel("select-parameter", described, matrix,
                                    loaded.num_samples))
        out.print(_cli_render.table(selection.columns, rows))
        # winner is legitimately None when every row was ineligible or
        # unscorable; that is a valid outcome, not a crash.
        if selection.winner is None:
            out.print(f"no winner: no eligible {parameter} value "
                      f"({len(rows)} evaluated)")
        else:
            won = selection.winner.result
            noise = sum(1 for value in won.labels if value < 0)
            fraction = noise / loaded.num_samples if loaded.num_samples else 0.0
            out.print(f"winner {parameter}={selection.winner.value} "
                      f"({len(rows)} evaluated) "
                      f"clusters={won.num_clusters} noise={noise} "
                      f"({fraction:.1%})")
    if output is not None:
        _cli_render.write_output(output, {
            "header": ["winner", *selection.columns],
            "rows": [[str(index == selection.winner_index).lower(), *row]
                     for index, row in enumerate(rows)],
            "document": {
                "schema_version": _cli_render.SCHEMA_VERSION,
                "command": "select-parameter",
                "input": {"path": matrix, "num_items": loaded.num_samples,
                          "orientation": str(loaded.is_distance)},
                "spec": described,
                "result": {
                    "columns": list(selection.columns),
                    # Not encoded here: write_output walks the whole document
                    # itself, precisely so a call site cannot leave a NaN in
                    # it by forgetting to.
                    "rows": [list(row) for row in rows],
                    "winner_index": selection.winner_index},
            }})


@cli.command()
@click.argument("matrix")
@click.option("--algorithm", required=True, help="Roster name.")
@click.option("--set", "assignments", multiple=True, metavar="KEY=VALUE",
              help="Algorithm option; repeatable.")
@click.option("--resamples", default=100, show_default=True,
              help="Number of resampled runs to score against.")
@click.option("--fraction", default=0.5, show_default=True,
              help="Fraction of the items drawn for each resample.")
@click.option("--seed", default=0, show_default=True,
              help="Seed for the resampling draws, for a reproducible run.")
@click.option("--noise", default="singletons", show_default=True,
              type=click.Choice(["singletons", "grouped", "excluded"]),
              help="How noise points count towards mean_agreement.")
@_common
@_translate
def stability(matrix, algorithm, assignments, resamples, fraction, seed,
              noise, threads, allow_nonmetric, output, quiet):
    """Score how well each cluster survives resampling."""
    _check_destinations(matrix, output)
    registry = _cli_registry.build()
    spec = _spec(algorithm, assignments, registry, threads, allow_nonmetric)
    loaded = _cli_input.load(matrix)
    out = _cli_render.console()
    # keep_partitions=False: this command exports only the reference
    # partition, so retaining every resample is pure memory cost. The
    # consensus command passes True, because it consenses them.
    with out.status(f"running {resamples} resamples"):
        scored = oecluster.cluster_stability(
            spec, loaded, resamples=resamples, fraction=fraction, seed=seed,
            noise=noise, keep_partitions=False,
            num_threads=_thread_count(threads))
    rows = scored.to_table()
    reference = [int(value) for value in scored.reference.labels]
    described = {"algorithm": algorithm, "options": dict(spec.options)}
    if not quiet:
        out.print(_cli_render.panel("stability", described, matrix,
                                    loaded.num_samples))
        out.print(_cli_render.table(scored.columns, rows))
        noise_count = sum(1 for value in reference if value < 0)
        noise_fraction = (noise_count / loaded.num_samples
                          if loaded.num_samples else 0.0)
        # mean_agreement is reported because --noise governs nothing else:
        # the Jaccard matching never treats noise as a cluster, so with the
        # agreement left out the flag changed no number the user could see.
        out.print(f"resamples={resamples} clusters={len(rows)} "
                  f"noise={noise_count} ({noise_fraction:.1%}) "
                  f"mean_jaccard={scored.mean_jaccard:.4f} "
                  f"mean_agreement={scored.mean_agreement:.4f}")
    if output is not None:
        _cli_render.write_output(output, {
            "header": list(scored.columns),
            "rows": rows,
            "document": {
                "schema_version": _cli_render.SCHEMA_VERSION,
                "command": "stability",
                "input": {"path": matrix, "num_items": loaded.num_samples,
                          "orientation": str(loaded.is_distance)},
                "spec": described,
                "result": {
                    "columns": list(scored.columns),
                    # Left to write_output's own walk, as its contract says.
                    "rows": [list(row) for row in rows],
                    "reference_labels": reference},
            }})


def main():
    """Console-script entry point."""
    cli()
