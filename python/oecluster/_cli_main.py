"""The ``oecluster`` command line.

Commands are flat and self-contained: each takes a distance matrix and does a
whole job. The algorithm and its options arrive as ``--algorithm NAME`` plus
repeatable ``--set key=value``, validated against the registry before any file
is read.
"""
import contextlib
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
        except (click.exceptions.Exit, click.Abort, click.ClickException):
            # Click's own control flow and its own already-carried errors.
            # Exit and Abort subclass RuntimeError, so the arm below would
            # turn ctx.exit(0) into exit 1 with the message "0" and an
            # aborted confirmation into a failed run; a ClickException
            # already names its exit code, which the catch-all at the end
            # would flatten to 1.
            raise
        except (_cli_input.InputError, RuntimeError, OSError) as error:
            raise click.ClickException(str(error)) from None
        except MemoryError as error:
            raise click.ClickException(
                f"out of memory: {error}; for consensus, --mmap keeps the "
                "matrix on disk") from None
        except Exception as error:  # noqa: BLE001
            # No reachable trigger for this arm was found, which is the
            # reason to keep it: the promise is that nothing reaches the
            # user as a traceback unless they asked for one, and an
            # enumerated list of types cannot hold that against library
            # calls that may raise anything. The type is named because the
            # message alone from an unanticipated exception is often
            # unreadable without it. --traceback never arrives here -- it
            # takes the branch above, which re-raises.
            raise click.ClickException(
                f"{type(error).__name__}: {error}") from None
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


def _spec(algorithm, assignments, registry, threads, nonmetric, *, swept=None,
          remedy="--set"):
    """Validate an algorithm and its options into a ClusteringSpec.

    :param remedy: The syntax this caller accepts options in, forwarded to
        :func:`_cli_registry.resolve` so its advice is followable; a
        ``consensus`` member is not configured with ``--set``.
    :raises ValueError: For an unknown or ineligible algorithm, or any option
        the registry refuses.
    """
    entry = _entry(algorithm, registry)
    options = _cli_registry.resolve(algorithm, assignments, registry,
                                    swept=swept, remedy=remedy)
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
    is the only authority on that, but it needs both paths to exist --
    hence :func:`_reserve`, which makes a destination exist so that this
    can be answered rather than guessed.
    """
    if os.path.realpath(first) == os.path.realpath(second):
        return True
    try:
        return os.path.samefile(first, second)
    except OSError:
        return False


def _reserve(targets):
    """Create every write target that does not exist yet, and say which.

    Two targets that do not exist cannot be told apart by asking the
    filesystem, and their names do not settle it in either direction: byte
    comparison calls ``r.json`` and ``R.JSON`` two files, so on a
    case-folding volume ``--output r.json --mmap R.JSON`` built the
    co-association matrix, mapped it, and then let the result writer
    truncate it in place -- exiting 0 with the artifact the user asked for
    destroyed. Folding the case unconditionally is wrong the other way: on
    a case-sensitive volume those really are two files and the invocation
    is legitimate, so a blanket fold would refuse a working command there.

    Making each one exist moves the question to the filesystem, which is
    the only thing that knows which behaviour it has.

    :param targets: Resolved write targets, none of them empty. They are
        resolved rather than as typed because that is the file the write
        lands on; passing a symbolic link here would reserve nothing, since
        POSIX requires ``O_CREAT | O_EXCL`` to fail with ``EEXIST`` on a
        link whatever it points at.
    :returns: The subset this call created, for the caller to remove again.
    """
    created = []
    for target in targets:
        try:
            handle = os.open(target, os.O_CREAT | os.O_EXCL, 0o600)
        except OSError:
            # Already there, or not creatable at all. Each of those is
            # settled by a check the caller runs anyway.
            continue
        os.close(handle)
        created.append(target)
    return created


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
        sidecar, if two destinations name one file, if a destination's
        parent is not an existing directory, or if ``output`` has an
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
    paths = [item for item in (output, *destinations) if item is not None]
    for path in paths:
        if not path:
            raise ValueError(
                "a destination path is empty; an unset shell variable is "
                "the usual cause")
    # Every question below is about the file the write lands on, and a
    # symbolic link is not that file -- `open` follows it. Resolving first
    # is what lets one rule cover a destination that exists, one that does
    # not, a link to either, and a dangling link, rather than a branch per
    # variant: three of those four shipped as separate data-loss bugs.
    # realpath resolves a dangling link to its target without requiring the
    # target to exist, which is what makes this a simplification.
    targets = [os.path.realpath(path) for path in paths]
    created = _reserve(targets)
    try:
        seen = []
        for path, target in zip(paths, targets):
            for other, what in protected:
                if _same_file(target, other):
                    raise ValueError(f"{path} is {what}; "
                                     "choose different files")
            for name, other in seen:
                if _same_file(target, other):
                    raise ValueError(f"{path} and {name} are the same "
                                     "file; choose different files")
            seen.append((path, target))
        # After the collisions, for the same reason the extension check is:
        # a destination that is both a collision and unwritable is reported
        # as the collision, which is the message naming the file at risk.
        # The directory is checked at all because `open` only fails where
        # the file is written -- for consensus that is after every member
        # has clustered and the co-association matrix has been built.
        for path, target in seen:
            parent = os.path.dirname(target)
            if not os.path.isdir(parent):
                raise ValueError(f"cannot write {path}: {parent} is not an "
                                 "existing directory")
        # Last, so a destination that is both unwritable and a collision is
        # reported as the collision: the message naming the file at risk.
        # The user's own path, not the target: write_output dispatches on
        # the path it is handed, so that is the extension that governs.
        if output is not None:
            _cli_render.check_output(output)
    finally:
        for target in created:
            # By target, not by the path the user typed: unlinking the
            # typed path would remove a link and leave the placeholder it
            # resolved to. Nothing has been written at this point either
            # way, so a placeholder left behind is an empty results file
            # for a run that never happened -- which reads as a run that
            # succeeded and found nothing.
            with contextlib.suppress(OSError):
                os.unlink(target)


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
                          "orientation": _cli_render.orientation(
                              loaded.is_distance)},
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
    if not values.strip():
        raise ValueError(
            "--values is empty; an unset shell variable is the usual cause")
    # Every field is kept and checked by position. Dropping the empty ones
    # instead meant `--values 0.8,,1.0` quietly evaluated two points and
    # named a winner from a grid the user had not written -- a silence no
    # amount of care at the keyboard could detect, since the three-point
    # run and the two-point run look identical in the output.
    fields = [item.strip() for item in values.split(",")]
    for position, field in enumerate(fields, start=1):
        if not field:
            raise ValueError(
                f"--values field {position} of {len(fields)} is empty; "
                "remove the comma or fill the gap")
    # Coerced through the same function --set uses: the grid does not pass
    # through _cli_registry.resolve(), so without this it reached the
    # library as untyped text -- a string threshold, an infinity, or a
    # mistyped exponent int() expands into a billion digits.
    swept = [_cli_registry.coerce(algorithm, parameter, field, entry)
             for field in fields]
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
        # unscorable; that is a valid outcome, not a crash. The two are
        # named apart because they call for opposite responses: loosen the
        # bounds, or pick a criterion the input can actually support.
        if selection.winner is None:
            reason = ("no eligible" if not any(row.eligible
                                               for row in selection.rows)
                      else "no scorable")
            out.print(f"no winner: {reason} {parameter} value "
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
            # Left as a boolean for write_output to spell. Spelling it here
            # instead is how this column and the `eligible` column beside
            # it ended up on two conventions in one row.
            "rows": [[index == selection.winner_index, *row]
                     for index, row in enumerate(rows)],
            "document": {
                "schema_version": _cli_render.SCHEMA_VERSION,
                "command": "select-parameter",
                "input": {"path": matrix, "num_items": loaded.num_samples,
                          "orientation": _cli_render.orientation(
                              loaded.is_distance)},
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
    # Before the load, as consensus checks the same option. The library
    # refuses this too, but only after the matrix has been read -- so the
    # unproven-orientation warning was printed for a run that could never
    # have happened, against the rule that nothing is read before
    # validation finishes.
    if resamples < 1:
        raise ValueError(f"--resamples must be positive, got {resamples}")
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
            # The CSV is the per-cluster table and nothing else, so
            # mean_agreement and the noise mode are JSON-only. That is the
            # convention already: consensus keeps its own mean_agreement out
            # of a CSV that is one row per item. A trailing summary row
            # would make every column change meaning on the last line, and
            # a repeated column would restate one scalar on every row.
            "header": list(scored.columns),
            "rows": rows,
            "document": {
                "schema_version": _cli_render.SCHEMA_VERSION,
                "command": "stability",
                "input": {"path": matrix, "num_items": loaded.num_samples,
                          "orientation": _cli_render.orientation(
                              loaded.is_distance)},
                "spec": described,
                "result": {
                    "columns": list(scored.columns),
                    # Left to write_output's own walk, as its contract says.
                    "rows": [list(row) for row in rows],
                    "reference_labels": reference,
                    # mean_agreement is the only number --noise governs --
                    # the Jaccard matching never treats noise as a cluster
                    # -- so without it here the flag changed nothing a
                    # script could read, which is the surface --output
                    # exists for. The mode travels with it because it is
                    # the statistic's unit: the same partitions score
                    # differently under each, so the number cannot be
                    # compared across runs that chose differently.
                    "mean_agreement": scored.mean_agreement,
                    "noise": noise},
            }})


def _member(text, registry, threads, nonmetric):
    """Parse one ``--member 'name;k=v;...'`` into a ClusteringSpec.

    The grammar is deliberately minimal: no roster option takes a value
    containing ``;`` or ``=``, so inventing an escape syntax for a case that
    does not exist would be the wrong trade. A value containing either is
    refused with that reason.

    :param text: One ``--member`` argument as the user spelled it.
    :param registry: Mapping from :func:`_cli_registry.build`.
    :param threads: The ``--threads`` value, routed into the member's spec.
    :param nonmetric: Whether ``--allow-nonmetric`` was given.
    :returns: The member's :class:`ClusteringSpec`.
    :raises ValueError: For a malformed member or any option the registry
        refuses.
    """
    parts = text.split(";")
    if not parts or not parts[0].strip():
        raise ValueError("--member needs an algorithm name")
    name = parts[0].strip()
    options = []
    for part in parts[1:]:
        item = part.strip()
        if not item:
            # An empty segment is a typo (a doubled or trailing ';'), not an
            # option; silently dropping it would hide the mistake.
            raise ValueError(f"empty option in --member {text!r}")
        if "=" not in item:
            raise ValueError(f"member option needs key=value, got {item!r}")
        key, value = item.split("=", 1)
        if "=" in value:
            raise ValueError(
                f"member option values may not contain '=', got {item!r}")
        options.append(f"{key.strip()}={value}")
    # The remedy names this member rather than --set, which the same command
    # refuses in cross-algorithm mode: `--member butina` was answered
    # "butina requires --set threshold=…", advice that cannot be followed.
    return _spec(name, options, registry, threads, nonmetric,
                 remedy=f"--member '{name};…'")


@cli.command()
@click.argument("matrix")
@click.option("--algorithm", default=None, help="Roster name, bootstrap mode.")
@click.option("--set", "assignments", multiple=True, metavar="KEY=VALUE",
              help="Algorithm option for bootstrap mode; repeatable.")
@click.option("--resamples", type=int, default=None,
              help="Bootstrap mode: resample count.")
@click.option("--member", "members", multiple=True, metavar="'NAME;K=V;...'",
              help="Cross-algorithm mode: one ensemble member; repeatable.")
@click.option("--threshold", type=float, default=None,
              help="Co-association support for the default extraction.")
@click.option("--noise", default="singletons", show_default=True,
              type=click.Choice(["singletons", "grouped", "excluded"]),
              help="How noise points count towards mean_agreement.")
@click.option("--mmap", "mmap", default=None, metavar="PATH",
              help="Keep the co-association matrix on disk, for large N.")
@_common
@_translate
def consensus(matrix, algorithm, assignments, resamples, members, threshold,
              noise, mmap, threads, allow_nonmetric, output, quiet):
    """Combine an ensemble of partitions into one."""
    bootstrap = algorithm is not None or resamples is not None
    if bootstrap and members:
        raise ValueError(
            "--member is the cross-algorithm mode; drop --algorithm and "
            "--resamples, or drop --member")
    if not bootstrap and not members:
        raise ValueError(
            "consensus needs either --algorithm with --resamples, or one or "
            "more --member")
    # Every destination must differ from every other and from the input, and
    # from the input's sidecar: with `d.npy` the obvious `--output d.json`
    # truncates exactly the sidecar the run depends on, and `--mmap` pointed
    # at the input overwrites the matrix while it is mapped. --mmap is a
    # collision only: nothing writes it through write_output, so the
    # .csv/.json contract applies to --output alone.
    _check_destinations(matrix, output, mmap)

    if threshold is not None and not 0.0 <= threshold <= 1.0:
        raise ValueError(f"--threshold must lie in [0, 1], got {threshold}")

    # Everything above and below this line runs before the matrix is read:
    # an invalid member or threshold must not cost a resampling run first.
    registry = _cli_registry.build()
    if bootstrap:
        if algorithm is None or resamples is None:
            raise ValueError(
                "bootstrap mode needs both --algorithm and --resamples")
        if resamples < 1:
            raise ValueError(f"--resamples must be positive, got {resamples}")
        spec = _spec(algorithm, assignments, registry, threads,
                     allow_nonmetric)
        described = {"algorithm": algorithm, "options": dict(spec.options),
                     "resamples": resamples}
        specs = None
    else:
        if assignments:
            raise ValueError(
                "--set belongs to bootstrap mode; put each member's options "
                "inside its own --member 'name;key=value'")
        specs = [_member(text, registry, threads, allow_nonmetric)
                 for text in members]
        described = {"members": [{"algorithm": item.name,
                                  "options": dict(item.options)}
                                 for item in specs]}

    # Every member runs over this one loaded matrix. Building a member's
    # input any other way would bypass the loader's label-collision refusal,
    # which _labels_of would then raise as a usage error at export time.
    loaded = _cli_input.load(matrix)
    out = _cli_render.console()
    if specs is None:
        with out.status(f"running {resamples} resamples"):
            ensemble = oecluster.cluster_stability(
                spec, loaded, resamples=resamples, keep_partitions=True,
                num_threads=_thread_count(threads))
    else:
        with out.status(f"running {len(specs)} members"):
            ensemble = [item.run(loaded) for item in specs]

    agreed = oecluster.consensus(
        ensemble,
        num_items=None if bootstrap else loaded.num_samples,
        threshold=threshold, noise=noise,
        num_threads=_thread_count(threads), output=mmap)
    rows = agreed.to_table()
    labels = [int(value) for value in agreed.labels]
    # Named apart from the `noise` mode it would otherwise shadow: the mode
    # is still needed below, as the unit mean_agreement is reported in.
    noise_count = sum(1 for value in labels if value < 0)
    if not quiet:
        out.print(_cli_render.panel("consensus", described, matrix,
                                    loaded.num_samples))
        out.print(_cli_render.table(agreed.columns, rows))
        fraction = (noise_count / loaded.num_samples
                    if loaded.num_samples else 0.0)
        out.print(f"clusters={agreed.num_clusters} noise={noise_count} "
                  f"({fraction:.1%}) partitions={agreed.num_partitions} "
                  f"mean_agreement={agreed.mean_agreement:.4f}")
    if output is not None:
        ids = _ids(loaded)
        named = _labels_of(loaded)
        _cli_render.write_output(output, {
            # One row per item, so the scalars below are JSON-only, exactly
            # as stability keeps its own out of a per-cluster table.
            "header": ["id", "label"],
            "rows": list(zip(ids, labels)),
            "document": {
                "schema_version": _cli_render.SCHEMA_VERSION,
                "command": "consensus",
                "input": {"path": matrix, "num_items": loaded.num_samples,
                          "orientation": _cli_render.orientation(
                              loaded.is_distance)},
                "spec": described,
                "result": {
                    "labels": labels,
                    **({"ids": named} if named else {}),
                    "num_clusters": agreed.num_clusters,
                    "num_noise": noise_count,
                    "num_partitions": agreed.num_partitions,
                    "unobserved_pairs": agreed.unobserved_pairs,
                    "threshold": _cli_render.encode(agreed.threshold),
                    "mean_agreement": _cli_render.encode(agreed.mean_agreement),
                    # The mode travels beside the statistic because it is
                    # its unit: --noise governs the agreement call behind
                    # mean_agreement and nothing else, so the same
                    # partitions score differently under each and the
                    # number cannot be compared across runs that chose
                    # differently. stability records the pair on the same
                    # terms.
                    "noise": noise,
                    "agreement": [_cli_render.encode(v)
                                  for v in agreed.agreement],
                    "item_consensus": [_cli_render.encode(float(v))
                                       for v in agreed.item_consensus],
                    "records": [{"label": r.label, "size": r.size,
                                 "cluster_consensus":
                                     _cli_render.encode(r.cluster_consensus)}
                                for r in agreed.records]},
            }})


def main():
    """Console-script entry point."""
    cli()
