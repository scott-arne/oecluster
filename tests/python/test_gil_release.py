"""The registered entry points release the GIL in the built module.

The GIL assertions in test_native_bindings.py read swig/oecluster.i, which
describes the interface rather than what SWIG generated from it. A registration
naming a function SWIG never wraps passes them, and so does an extension built
before a registration existed. The tests here look at what a build produced.

The first scans the wrapper source in the build tree whose shadow module the
import uses. For each function that runs a registered call, it replays the GIL
state through the function's text, from held at the top through every release
and re-acquire the text names, and requires the call to be reached released. It
fails a function whose text it cannot follow: a release or re-acquire out of
turn, a second call, a preprocessor conditional, a release that does not begin a
statement, or a branch, label or block end between that release and the call. A
GIL change or jump the text does not spell out, such as a helper function that
re-acquires the GIL internally, is invisible to that scan.

The second times calls into the imported extension while a second thread tries
to run, which is the check for such a hidden re-acquire on the calls it covers.
"""

import logging
import pathlib
import re
import sys
import threading
import time

import pytest

_LOG = logging.getLogger(__name__)

_REPO = pathlib.Path(__file__).resolve().parents[2]
_INTERFACE = _REPO / "swig" / "oecluster.i"
_BUILD_ROOT = _REPO / "build"

# Both comment forms SWIG honours, stripped before registrations are read so
# that a commented-out registration does not count as one. Double-quoted
# strings are matched too and kept, so that a comment marker inside one, such
# as the "/*" of a glob pattern, does not start a comment.
_SWIG_COMMENT_OR_STRING = re.compile(r'"(?:\\.|[^"\\\n])*"|/\*.*?\*/|//[^\n]*', re.DOTALL)
_REGISTRATION = re.compile(r"\bOECLUSTER_GIL_(?:MEMORY_)?EXCEPTION\(OECluster::(\w+),")
# Looser than _REGISTRATION, to count the invocations it fails to parse.
_INVOCATION = re.compile(r"\bOECLUSTER_GIL_\w*EXCEPTION\s*\(")
_DEFINE_LINE = re.compile(r"^[ \t]*%define\b[^\n]*", re.MULTILINE)
# No brace inside a SWIG wrapper body sits alone at the start of a line, so the
# first line holding only "}" closes the function.
_WRAPPER_FUNCTION = re.compile(r"^SWIGINTERN PyObject \*(_wrap_\w+)\(.*?^\}$",
                               re.MULTILINE | re.DOTALL)
# The calls that change whether the running thread holds the GIL. SWIG's
# SWIG_PYTHON_THREAD_* macros expand to PyGILState_Ensure, PyGILState_Release,
# PyEval_SaveThread and PyEval_RestoreThread when SWIG's thread support is on,
# and to nothing otherwise; they are counted as the calls they can become.
_RELEASING = ("Py_BEGIN_ALLOW_THREADS", "Py_UNBLOCK_THREADS", "PyEval_SaveThread",
              "PyEval_ReleaseThread", "PyGILState_Release",
              "SWIG_PYTHON_THREAD_BEGIN_ALLOW", "SWIG_PYTHON_THREAD_END_BLOCK")
_ACQUIRING = ("Py_END_ALLOW_THREADS", "Py_BLOCK_THREADS", "PyEval_RestoreThread",
              "PyEval_AcquireThread", "PyGILState_Ensure",
              "SWIG_PYTHON_THREAD_END_ALLOW", "SWIG_PYTHON_THREAD_BEGIN_BLOCK")
_GIL_CALL = re.compile(rf"\b(?:{'|'.join(_RELEASING + _ACQUIRING)})\b")
_BRANCH = re.compile(r"\b(?:if|else|for|while|do|switch|case|default|goto)\b")
_LABEL = re.compile(r"(?:^|[;{}])[ \t]*(\w+)[ \t]*:(?!:)", re.MULTILINE)
_CONDITIONAL_LINE = re.compile(r"^[ \t]*#[ \t]*(?:if|ifdef|ifndef|elif|else|endif)\b",
                               re.MULTILINE)
_QUALIFIED_CALL = re.compile(r"\bOECluster::\w+\(")
# Under build/, this repository builds into build/ itself (the default CMake
# preset) and build/{wheel_tag} (pyproject's build-dir).
_SHADOW_MODULE_PATTERNS = ("swig/_oecluster.py", "*/swig/_oecluster.py")


def _registered_names(interface_text):
    """Return the function names registered for GIL release, in file order.

    :param interface_text: The text of swig/oecluster.i.
    :returns: One unqualified name per registration.
    :raises AssertionError: When an invocation of a GIL macro, outside
        comments and ``%define`` lines, does not parse as a registration, since
        the scan would otherwise leave it unchecked.
    """
    code = _SWIG_COMMENT_OR_STRING.sub(
        lambda match: match.group(0) if match.group(0).startswith('"') else "",
        interface_text)
    names = _REGISTRATION.findall(code)
    invocations = len(_INVOCATION.findall(_DEFINE_LINE.sub("", code)))
    assert invocations == len(names), (
        f"{invocations} GIL macro invocations, but only {len(names)} parse as "
        "OECLUSTER_GIL_[MEMORY_]EXCEPTION(OECluster::<name>, ...); the rest "
        "would go unchecked")
    return names


def _blank(match):
    """Blank a string literal to ``""`` and a comment to its own line breaks."""
    text = match.group(0)
    if text.startswith('"'):
        return '""'
    return "\n" * text.count("\n") or " "


def _release_problem(function, body, name):
    """Say why ``body`` does not call ``OECluster::<name>(`` released, or None.

    :param function: The wrapper function's name, for the message.
    :param body: The wrapper function's text.
    :param name: The registered function name.
    :returns: A description of the first violation found, or None.
    """
    code = _SWIG_COMMENT_OR_STRING.sub(_blank, body)

    def line(position):
        return code.count("\n", 0, position) + 1

    action = f"OECluster::{name}("
    calls = [match.start() for match in re.finditer(rf"\bOECluster::{name}\(", code)]
    if not calls:
        held = sorted({f"{call}...)" for call in _QUALIFIED_CALL.findall(code)})
        return (f"{function} never calls {action}...); it calls "
                f"{', '.join(held) or 'nothing in OECluster::'}")
    if len(calls) > 1:
        return (f"{function} calls {action}...) {len(calls)} times, and the replay "
                "vouches for one call only")
    call = calls[0]
    conditional = _CONDITIONAL_LINE.search(code, 0, call)
    if conditional:
        return (f"{function} has a preprocessor conditional on line "
                f"{line(conditional.start())}, before its call to {action}...)")

    held, last = True, None
    for match in _GIL_CALL.finditer(code, 0, call):
        releases = match.group(0) in _RELEASING
        if releases != held:
            state = "released" if releases else "held"
            return (f"{function} is malformed: {match.group(0)} on line "
                    f"{line(match.start())} finds the GIL already {state}")
        held, last = not releases, match
    if held:
        after = (f"after {last.group(0)} on line {line(last.start())}" if last
                 else "with no release before it")
        return f"{function} reaches {action}...) holding the GIL, {after}"

    # The replay reads the text top to bottom. That is the order that runs only
    # if the release it ended on always runs before the call: it must begin a
    # statement, so no condition governs it alone, and the text between it and
    # the call may hold no branch, no jump target, and no end of the block the
    # release sits in, which would also run the destructor of a scoped release.
    release = last.start()
    preceding = code[:release].rstrip()
    if preceding and preceding[-1] not in ";{}":
        return (f"{function}'s {last.group(0)} on line {line(release)} does not "
                "begin a statement, so a condition may skip it")
    branch = _BRANCH.search(code, release, call)
    if branch:
        return (f"{function} has `{branch.group(0)}` on line {line(branch.start())}, "
                f"between its release and its call to {action}...)")
    label = _LABEL.search(code, release, call)
    if label:
        return (f"{function} has the label {label.group(1)}: on line "
                f"{line(label.start(1))}, between its release and its call to "
                f"{action}...)")
    depth = 0
    for position in range(release, call):
        depth += {"{": 1, "}": -1}.get(code[position], 0)
        if depth < 0:
            return (f"{function} closes the block holding its release, on line "
                    f"{line(position)}, before its call to {action}...)")
    return None


def _unreleased(names, wrapper_text):
    """Describe each registered call the wrapper does not run with the GIL released.

    The functions that run a registered call are the leaves of an overloaded
    name, ``_wrap_<name>__SWIG_<n>``, or ``_wrap_<name>`` for a name with none.
    Comments and string literals blanked, each must call ``OECluster::<name>(``
    exactly once, with no preprocessor conditional before the call, and must
    pass a replay of its text from the top down to the call. The replay starts
    with the GIL held; each call in _RELEASING releases it and each in
    _ACQUIRING re-acquires it, and a release while released or a re-acquire
    while held fails the function as malformed. The call must be reached
    released. So that the text order is the order that runs, the release the
    replay ended on must begin a statement, and between it and the call there
    may be no if, else, for, while, do, switch, case, default or goto, no
    label, and no brace closing the block the release sits in.

    A GIL change or jump the text does not spell out, such as a helper function
    that re-acquires the GIL internally, is invisible to this scan; layer 2's
    runtime probes are the check for that, on the calls they cover.

    :param names: Registered function names.
    :param wrapper_text: The text of the generated oeclusterPYTHON_wrap.cxx.
    :returns: One line per registration with no wrapper function, and one per
        function that never calls ``OECluster::<name>(`` or breaks the rule
        above.
    """
    functions = {match.group(1): match.group(0)
                 for match in _WRAPPER_FUNCTION.finditer(wrapper_text)}
    problems = []
    for name in names:
        # An overloaded name gets one leaf per overload, which runs the call,
        # and a dispatcher under the plain name, which only picks a leaf and
        # carries no marker.
        leaf = re.compile(rf"_wrap_{name}__SWIG_\d+")
        checked = sorted(function for function in functions if leaf.fullmatch(function))
        if not checked and f"_wrap_{name}" in functions:
            checked = [f"_wrap_{name}"]
        if not checked:
            problems.append(f"{name} is registered, but the wrapper has no "
                            f"_wrap_{name} function")
        for function in checked:
            problem = _release_problem(function, functions[function], name)
            if problem is not None:
                problems.append(f"{name}: {problem}")
    return problems


def _bound_wrapper(build_root, module_file):
    """Return the wrapper source in the build tree the import uses, or skip.

    A build tree's swig/_oecluster.py is the shadow module a wheel installs,
    and the development build copies, as oecluster/oecluster.py, so a
    byte-identical copy identifies the build tree whose shadow module is
    imported. It does not identify the extension: a GIL registration is an
    ``%exception`` directive, which changes the C++ wrapper and never the
    shadow module, so two trees whose interfaces differ only in registrations
    match alike. The scan therefore assumes that the matching tree's wrapper is
    the one the imported extension was compiled from. That fails when the
    wrapper was regenerated but did not compile, or when a tree outside the
    search, such as build-debug, wrote the imported module and an older tree
    here matches it.

    :param build_root: Directory whose swig/ and */swig/ subdirectories are
        searched for _oecluster.py.
    :param module_file: The imported oecluster.oecluster module's file.
    :returns: Path of the matching tree's oeclusterPYTHON_wrap.cxx.
    """
    module_bytes = module_file.read_bytes()
    candidates = (sorted(path.parent for pattern in _SHADOW_MODULE_PATTERNS
                         for path in build_root.glob(pattern))
                  if build_root.is_dir() else [])
    matches = [directory for directory in candidates
               if (directory / "_oecluster.py").read_bytes() == module_bytes]
    if len(matches) == 1:
        return matches[0] / "oeclusterPYTHON_wrap.cxx"
    if matches:
        pytest.skip(f"{len(matches)} build trees hold a swig/_oecluster.py "
                    f"identical to {module_file}, so which wrapper to scan is "
                    f"ambiguous: {', '.join(map(str, matches))}")
    inspected = ", ".join(map(str, candidates)) or "none found"
    pytest.skip(f"neither {build_root}/swig/_oecluster.py nor any "
                f"{build_root}/*/swig/_oecluster.py is identical to "
                f"{module_file}; directories inspected: {inspected}")


def test_every_registered_wrapper_releases_the_gil():
    """Every function that runs a registered call releases the GIL.

    Scans the wrapper source in the build tree whose shadow module the import
    uses, assuming that wrapper is the one the imported extension was compiled
    from (see _bound_wrapper), and skips when no single build tree matches. A
    registration SWIG never applied fails it: one naming a function SWIG does
    not wrap has no wrapper function, and one placed after the declaration
    leaves the call outside any release.
    """
    from oecluster import oecluster as native

    wrapper = _bound_wrapper(_BUILD_ROOT, pathlib.Path(native.__file__))
    names = _registered_names(_INTERFACE.read_text(encoding="utf-8"))
    assert names, f"no GIL registrations found in {_INTERFACE}"
    problems = _unreleased(names, wrapper.read_text(encoding="utf-8"))
    assert not problems, "\n".join(problems)


# One registration among lines that name a GIL macro but are not registrations:
# both %define lines and two commented-out invocations.
_SYNTHETIC_INTERFACE = """\
%define OECLUSTER_GIL_EXCEPTION(QUALIFIED, NAME)
%enddef
%define OECLUSTER_GIL_MEMORY_EXCEPTION(QUALIFIED, NAME)
%enddef
// OECLUSTER_GIL_EXCEPTION(OECluster::commented, commented)
/* OECLUSTER_GIL_MEMORY_EXCEPTION(OECluster::blocked, blocked) */
OECLUSTER_GIL_EXCEPTION(OECluster::first, first)
"""


def test_the_registration_parser_skips_defines_and_comments():
    assert _registered_names(_SYNTHETIC_INTERFACE) == ["first"]


@pytest.mark.parametrize("registration", [
    "OECLUSTER_GIL_EXCEPTION( OECluster::spaced, spaced)",
    "OECLUSTER_GIL_MEMORY_EXCEPTION(OECluster::detail::nested, nested)",
], ids=["space-after-parenthesis", "nested-scope"])
def test_a_registration_that_does_not_parse_is_refused(registration):
    with pytest.raises(AssertionError, match="2 GIL macro invocations, but only 1 parse"):
        _registered_names(_SYNTHETIC_INTERFACE + registration + "\n")


def test_a_comment_marker_in_a_string_does_not_hide_a_registration():
    text = (_SYNTHETIC_INTERFACE
            + '%pythoncode %{\n_SOURCES = "src/*.py"\n%}\n'
            + "OECLUSTER_GIL_EXCEPTION(OECluster::second, second)\n"
            + "// a later */ that a comment opened inside the string would close on\n")
    assert _registered_names(text) == ["first", "second"]


def _synthetic_wrapper(body):
    """Wrap ``body`` in the shape SWIG gives an unoverloaded wrapper of demo."""
    return ("SWIGINTERN PyObject *_wrap_demo(PyObject *self, PyObject *args) {\n"
            "  PyObject *resultobj = 0;\n"
            f"{body}"
            "  return resultobj;\n"
            "fail:\n"
            "  return NULL;\n"
            "}\n")


def test_a_registered_call_inside_its_release_is_accepted():
    body = ("  {\n"
            "    Py_BEGIN_ALLOW_THREADS\n"
            "    try {\n"
            "      result = OECluster::demo(arg1);\n"
            "    } catch (...) {\n"
            "      Py_BLOCK_THREADS\n"
            '      SWIG_exception(SWIG_RuntimeError, "Unknown C++ exception in " "demo");\n'
            "    }\n"
            "    Py_END_ALLOW_THREADS\n"
            "  }\n")
    assert _unreleased(["demo"], _synthetic_wrapper(body)) == []


def test_a_registration_with_no_wrapper_function_is_refused():
    """A registration naming a function SWIG never wrapped, as the removed
    ``cdist`` one did, has no function to check, and is reported rather than
    passed."""
    problems = _unreleased(["cdist"], _synthetic_wrapper(""))
    assert problems == [
        "cdist is registered, but the wrapper has no _wrap_cdist function"]


@pytest.mark.parametrize("body, reason", [
    (("  result = OECluster::demo(arg1);\n"
      "  Py_BEGIN_ALLOW_THREADS\n"
      "  Py_END_ALLOW_THREADS\n"),
     "holding the GIL, with no release before it"),
    (("#if 0\n"
      "  Py_BEGIN_ALLOW_THREADS\n"
      "#endif\n"
      "  result = OECluster::demo(arg1);\n"
      "#if 0\n"
      "  Py_END_ALLOW_THREADS\n"
      "#endif\n"),
     "preprocessor conditional"),
    (("  Py_BEGIN_ALLOW_THREADS\n"
      "  prepare(arg1);\n"
      "  Py_END_ALLOW_THREADS\n"
      "  result = OECluster::demo(arg1);\n"),
     "holding the GIL, after Py_END_ALLOW_THREADS"),
    (("  Py_BEGIN_ALLOW_THREADS\n"
      "  result = OECluster::other(arg1);\n"
      "  Py_END_ALLOW_THREADS\n"),
     "never calls OECluster::demo(...); it calls OECluster::other("),
    (("  Py_BEGIN_ALLOW_THREADS\n"
      "  Py_BLOCK_THREADS\n"
      "  result = OECluster::demo(arg1);\n"
      "  Py_UNBLOCK_THREADS\n"
      "  Py_END_ALLOW_THREADS\n"),
     "holding the GIL, after Py_BLOCK_THREADS"),
    (("  Py_BEGIN_ALLOW_THREADS\n"
      "  PyGILState_STATE state = PyGILState_Ensure();\n"
      "  result = OECluster::demo(arg1);\n"
      "  PyGILState_Release(state);\n"
      "  Py_END_ALLOW_THREADS\n"),
     "holding the GIL, after PyGILState_Ensure"),
    (("  Py_BEGIN_ALLOW_THREADS\n"
      "  Py_UNBLOCK_THREADS\n"
      "  result = OECluster::demo(arg1);\n"
      "  Py_END_ALLOW_THREADS\n"),
     "malformed: Py_UNBLOCK_THREADS on line 4 finds the GIL already released"),
    (("  Py_BEGIN_ALLOW_THREADS\n"
      "  if (ready) {\n"
      "    result = OECluster::demo(arg1);\n"
      "  }\n"
      "  Py_END_ALLOW_THREADS\n"),
     "has `if` on line 4"),
    (("  if (ready) {\n"
      "    Py_BEGIN_ALLOW_THREADS\n"
      "  }\n"
      "  result = OECluster::demo(arg1);\n"
      "  Py_END_ALLOW_THREADS\n"),
     "closes the block"),
    (("  if (ready)\n"
      "    Py_BEGIN_ALLOW_THREADS\n"
      "  result = OECluster::demo(arg1);\n"
      "  Py_END_ALLOW_THREADS\n"),
     "does not begin a statement"),
    (("  goto released;\n"
      "  Py_BEGIN_ALLOW_THREADS\n"
      "released:\n"
      "  result = OECluster::demo(arg1);\n"
      "  Py_END_ALLOW_THREADS\n"),
     "label released: on line 5"),
], ids=["marker-after-the-call", "marker-only-under-if-0", "call-after-the-end",
        "no-call", "blocked-before-the-call", "gilstate-ensure-before-the-call",
        "double-release", "branch-after-the-release", "release-in-a-closed-block",
        "release-under-a-braceless-if", "label-after-the-release"])
def test_a_registered_call_outside_its_release_is_refused(body, reason):
    problems = _unreleased(["demo"], _synthetic_wrapper(body))
    assert len(problems) == 1, problems
    assert reason in problems[0], problems


# A call reads as released when the ticker's longest gap stays under this
# fraction of the call's duration.
_RELEASED_FRACTION = 0.25
# Shorter than this, the switch interval and scheduler jitter are too large a
# share of a released call for the fraction to separate it from a held one.
_MIN_DURATION = 0.050
_SWITCH_INTERVAL = 0.005
_ATTEMPTS = 3
_THREAD_TIMEOUT = 10.0

# _PROBE_ITEMS, with the multiples of it that murcko_scaffolds and the control
# take, was sized so that every released probe ran for at least twice
# _MIN_DURATION on the machine this test was written on.
_PROBE_ITEMS = 3000
_PROBE_SMILES = (
    "CC(=O)Oc1ccccc1C(=O)O", "c1ccc2c(c1)cccc2CCN", "O=C(Nc1ccc(F)cc1)c1ccncc1",
    "CN1CCN(CC1)c1ccc(cc1)C(=O)Nc1ccccc1", "c1ccc(cc1)C1CCCCC1",
    "Cc1ccc(S(=O)(=O)N)cc1", "COc1ccc2[nH]cc(CCN(C)C)c2c1",
    "OC(=O)C1CCN(Cc2ccccc2)CC1",
)


def _longest_gap(call):
    """Run ``call`` on this thread while another ticks, and time both.

    The ticker is pure Python, so it can only tick while it holds the GIL. A
    call that releases the GIL around its native work leaves it free to tick
    during that work; one that holds it stops the ticker for the whole native
    call, even when the call is so short that no tick lands inside it. The
    ticker starts before the call and stops after it, so the longest gap can
    only overstate the call's own, which can make a released call look held but
    never the reverse.

    :param call: Zero-argument callable to time.
    :returns: The longest gap between ticks, and the call's duration, in
        seconds.
    """
    stop = threading.Event()
    started = threading.Event()
    longest = 0.0

    def tick():
        nonlocal longest
        last = time.perf_counter()
        started.set()
        while True:
            # Reading the flag before the clock puts the last tick after stop
            # was set, and so after the call returned: the gap that spans the
            # end of the call is always counted.
            stopping = stop.is_set()
            now = time.perf_counter()
            longest = max(longest, now - last)
            last = now
            if stopping:
                return

    ticker = threading.Thread(target=tick, name="gil-probe-ticker", daemon=True)
    previous = sys.getswitchinterval()
    sys.setswitchinterval(_SWITCH_INTERVAL)
    try:
        ticker.start()
        assert started.wait(_THREAD_TIMEOUT), "the ticker thread never started"
        begin = time.perf_counter()
        call()
        duration = time.perf_counter() - begin
    finally:
        stop.set()
        if ticker.is_alive():
            ticker.join(_THREAD_TIMEOUT)
        sys.setswitchinterval(previous)
    assert not ticker.is_alive(), "the ticker thread outlived its probe"
    return longest, duration


@pytest.fixture(scope="module")
def probe_inputs():
    """Molecules, their fingerprint comparison, and its distance matrix."""
    from oecluster import oecluster as native
    from openeye import oechem

    mols = []
    for index in range(_PROBE_ITEMS):
        mol = oechem.OEGraphMol()
        assert oechem.OESmilesToMol(mol, _PROBE_SMILES[index % len(_PROBE_SMILES)])
        mols.append(mol)
    comparison = native.FingerprintComparison(mols, native.FingerprintOptions())
    storage = native.DenseStorage(len(mols))
    native.pdist(comparison, storage)
    return mols, comparison, storage


# Every probe asks for one worker thread, so on a runner with more than one
# core the call leaves a core free for the ticker.
def _butina_cluster_on_a_comparison(native, mols, comparison, storage):
    options = native.ButinaOptions()
    options.distance_threshold = 0.3
    options.num_threads = 1
    return lambda: native.butina_cluster(comparison, options)


def _k_medoids_cluster_on_a_matrix(native, mols, comparison, storage):
    options = native.KMedoidsOptions()
    options.n_clusters = 6
    options.num_threads = 1
    return lambda: native.k_medoids_cluster(storage, options)


def _murcko_scaffolds(native, mols, comparison, storage):
    options = native.MurckoOptions()
    options.num_threads = 1
    doubled = mols * 2
    return lambda: native.murcko_scaffolds(doubled, options)


@pytest.mark.parametrize("probe", [
    _butina_cluster_on_a_comparison,
    _k_medoids_cluster_on_a_matrix,
    _murcko_scaffolds,
], ids=["butina_cluster-comparison", "k_medoids_cluster-matrix", "murcko_scaffolds"])
def test_a_released_call_lets_another_thread_run(probe, probe_inputs):
    """A registered call leaves the GIL to a second thread while it runs.

    descriptor_statistics is not registered, so it holds the GIL, and it must
    read as held first: without that, a probe that read every call as released
    would pass vacuously.
    """
    from oecluster import oecluster as native

    mols, comparison, storage = probe_inputs
    name = probe.__name__.lstrip("_")

    # The control is not held to _MIN_DURATION: the ticker cannot run while
    # the call holds the GIL, so its longest gap covers the native call however
    # long that takes.
    control_mols = mols * 4
    control_options = native.DescriptorStatisticsOptions()
    gap, duration = _longest_gap(
        lambda: native.descriptor_statistics(control_mols, control_options))
    _LOG.info("held control: %.3f of %.1f ms", gap / duration, duration * 1e3)
    assert gap >= _RELEASED_FRACTION * duration, (
        f"the held control read as released (longest gap {gap * 1e3:.1f} ms of "
        f"a {duration * 1e3:.1f} ms call), so the probe cannot see a held GIL")

    # A released call passes on the first of up to _ATTEMPTS attempts that
    # reads as released. Retrying cannot pass a call that holds the GIL: each
    # probe spends nearly all of its duration in native work, and a held call
    # stops the ticker for that work on every attempt. What retrying absorbs is
    # the descheduling of the ticker itself, which on a loaded runner can make
    # a released attempt read as held. The control above runs once, because
    # passing it on any of several attempts would weaken the proof that the
    # probe can see a held GIL.
    call = probe(native, mols, comparison, storage)
    attempts = []
    for _ in range(_ATTEMPTS):
        gap, duration = _longest_gap(call)
        _LOG.info("%s: %.3f of %.1f ms", name, gap / duration, duration * 1e3)
        if duration < _MIN_DURATION:
            pytest.skip(f"{name} took {duration * 1e3:.1f} ms, under the "
                        f"{_MIN_DURATION * 1e3:.0f} ms a released call needs to "
                        "be told apart from a held one")
        attempts.append((gap, duration))
        if gap < _RELEASED_FRACTION * duration:
            return
    pytest.fail(f"{name} read as held on all {_ATTEMPTS} attempts; longest gap "
                "as a fraction of the call: " + ", ".join(
                    f"{gap / duration:.2f} ({gap * 1e3:.1f} of {duration * 1e3:.1f} ms)"
                    for gap, duration in attempts))
