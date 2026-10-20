"""The registered entry points release the GIL in the built module.

The GIL assertions in test_native_bindings.py read swig/oecluster.i, which
describes the interface rather than what SWIG generated from it. A registration
naming a function SWIG never wraps passes them, and so does an extension built
before a registration existed. The tests here look at what the build produced:
the first scans the wrapper source SWIG generated with the imported module, and
the second times calls into the imported extension while a second thread tries
to run.
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
# that a commented-out registration does not count as one.
_SWIG_COMMENT = re.compile(r"/\*.*?\*/|//[^\n]*", re.DOTALL)
_REGISTRATION = re.compile(r"\bOECLUSTER_GIL_(?:MEMORY_)?EXCEPTION\(OECluster::(\w+),")
# No brace inside a SWIG wrapper body sits alone at the start of a line, so the
# first line holding only "}" closes the function.
_WRAPPER_FUNCTION = re.compile(r"^SWIGINTERN PyObject \*(_wrap_\w+)\(.*?^\}$",
                               re.MULTILINE | re.DOTALL)
# The first line of both OECLUSTER_GIL_EXCEPTION and
# OECLUSTER_GIL_MEMORY_EXCEPTION.
_RELEASE_MARKER = "Py_BEGIN_ALLOW_THREADS"


def _registered_names(interface_text):
    """Return the function names registered for GIL release, in file order.

    :param interface_text: The text of swig/oecluster.i.
    :returns: One unqualified name per registration.
    """
    return _REGISTRATION.findall(_SWIG_COMMENT.sub("", interface_text))


def _unreleased(names, wrapper_text):
    """Describe each registered call the wrapper does not release the GIL for.

    :param names: Registered function names.
    :param wrapper_text: The text of the generated oeclusterPYTHON_wrap.cxx.
    :returns: One line per registration with no wrapper function, and one per
        function that runs a registered call without the release marker in its
        own body.
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
        problems.extend(f"{name}: {function} does not contain {_RELEASE_MARKER}"
                        for function in checked
                        if _RELEASE_MARKER not in functions[function])
    return problems


def _bound_wrapper(build_root, module_file):
    """Return the wrapper source generated with the imported module, or skip.

    A build directory's swig/_oecluster.py is the file a wheel installs, and
    the development build copies, as oecluster/oecluster.py, so a
    byte-identical copy ties the imported module to the wrapper SWIG generated
    beside it.

    :param build_root: Directory searched recursively for swig/_oecluster.py.
    :param module_file: The imported oecluster.oecluster module's file.
    :returns: Path of the matching directory's oeclusterPYTHON_wrap.cxx.
    """
    module_bytes = module_file.read_bytes()
    candidates = (sorted(path.parent for path in build_root.rglob("swig/_oecluster.py"))
                  if build_root.is_dir() else [])
    matches = [directory for directory in candidates
               if (directory / "_oecluster.py").read_bytes() == module_bytes]
    if len(matches) == 1:
        return matches[0] / "oeclusterPYTHON_wrap.cxx"
    if matches:
        pytest.skip(f"{len(matches)} build directories hold a swig/_oecluster.py "
                    f"identical to {module_file}, so the wrapper behind it is "
                    f"ambiguous: {', '.join(map(str, matches))}")
    inspected = ", ".join(map(str, candidates)) or "none found"
    pytest.skip(f"no swig/_oecluster.py under {build_root} is identical to "
                f"{module_file}; directories inspected: {inspected}")


def test_every_registered_wrapper_releases_the_gil():
    """Every function that runs a registered call releases the GIL.

    Reads the wrapper source SWIG generated alongside the imported module,
    and skips when no single build directory holds it. A registration SWIG
    never applied, whether to a function it does not wrap or because the
    directive sits after the declaration, leaves no marker to find.
    """
    from oecluster import oecluster as native

    wrapper = _bound_wrapper(_BUILD_ROOT, pathlib.Path(native.__file__))
    names = _registered_names(_INTERFACE.read_text(encoding="utf-8"))
    assert names, f"no GIL registrations found in {_INTERFACE}"
    problems = _unreleased(names, wrapper.read_text(encoding="utf-8"))
    assert not problems, "\n".join(problems)


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
