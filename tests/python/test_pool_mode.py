"""Under an OpenEye memory-pool mode that is not thread-safe, oecluster runs
OpenEye work one call and one thread at a time.

``OESetMemPoolMode`` is fatal on a second call, so every case runs in a fresh
child interpreter that sets the mode before it creates a molecule. Mode 8
(UnboundedCache) has neither the Mutexed nor the ThreadLocal bit and is
unsafe; mode 10 (Default) is the control. Before the gate, the comparison
cases crashed the child under mode 8 on at least one build type.
"""

import ctypes
import ctypes.util
import json
import os
import queue
import subprocess
import sys
import textwrap
import threading
import time

import pytest

pytest.importorskip("openeye.oechem")

SAFE = 10
UNSAFE = 8

_PRELUDE = textwrap.dedent(r'''
    import hashlib
    import json
    import math
    import sys
    import threading
    import time

    mode, case = int(sys.argv[1]), sys.argv[2]

    from openeye import oechem

    oechem.OESetMemPoolMode(mode)

    import numpy as np
    import oecluster
    from oecluster import oecluster as native

    SMILES = ["c1ccc2ccccc2c1CCN", "C1CCC(CC1)c1ccncc1",
              "O=C(Nc1ccccc1)c1ccco1", "c1ccc(cc1)C1CCNCC1",
              "CC(=O)Nc1ccc(O)cc1", "c1ccc2c(c1)[nH]c1ccccc12",
              "O=C1CCCN1c1ccccc1", "CCOC(=O)c1ccc(N)cc1",
              "c1ccc(-c2ccccn2)nc1", "CN1CCN(CC1)c1ccc2ncccc2c1"]


    def molecules(kind, count, one_graph=False, three_d=False):
        """Deterministic inputs: varied graphs, or one graph for RMSD, with
        synthetic non-planar coordinates where the comparison needs 3D."""
        made = []
        for index in range(count):
            mol = kind()
            smiles = SMILES[0] if one_graph else SMILES[index % len(SMILES)]
            oechem.OESmilesToMol(mol, smiles)
            if three_d:
                coords = []
                for atom_index, _ in enumerate(mol.GetAtoms()):
                    angle = 0.9 * atom_index + 0.01 * index
                    coords.extend([1.5 * math.cos(angle), 1.5 * math.sin(angle),
                                   0.4 * atom_index + 0.001 * index])
                mol.SetCoords(oechem.OEFloatArray(coords))
                mol.SetDimension(3)
            made.append(mol)
        return made


    def digest(values):
        array = np.round(np.asarray(values, dtype=np.float64), 6)
        return hashlib.sha256(array.tobytes()).hexdigest()[:16]


    def pdist_digest(mols, name, **kwargs):
        matrix = oecluster.pdist(mols, name, num_threads=8, chunk_size=16, **kwargs)
        return digest(matrix.condensed)


    def murcko_digest(mols):
        found = oecluster.murcko_scaffolds(mols, num_threads=8)
        return hashlib.sha256(repr(list(found)).encode()).hexdigest()[:16]


    def two_at_once(mols):
        """Two direct native pdist calls released together by a barrier, so
        they overlap unless the binding keeps them apart."""
        jobs = []
        for _ in range(2):
            options = native.PDistOptions()
            options.num_threads = 8
            options.chunk_size = 16
            jobs.append((native.RMSDComparison(mols),
                         native.DenseStorage(len(mols)), options))
        barrier = threading.Barrier(2)
        results = [None, None]

        def run(slot):
            comparison, storage, options = jobs[slot]
            barrier.wait()
            native.pdist(comparison, storage, options)
            results[slot] = digest([storage.Get(i, j)
                                    for i in range(len(mols))
                                    for j in range(i + 1, len(mols))])

        workers = [threading.Thread(target=run, args=(slot,)) for slot in (0, 1)]
        for worker in workers:
            worker.start()
        for worker in workers:
            worker.join()
        return results


    def with_progress(mols, name):
        calls = [0]

        def progress(done, total):
            calls[0] += 1

        return {"digest": pdist_digest(mols, name, progress=progress),
                "called": calls[0] > 0}


    def gil_fraction(mols):
        """The share of its free rate a busy Python thread keeps while a
        direct native pdist runs. With no Python callback, the only window in
        which it can run while the call keeps the GIL is the call's boundary,
        which a short switch interval keeps short."""
        sys.setswitchinterval(0.0005)
        ticks = [0]
        running = [True]

        def spin():
            while running[0]:
                ticks[0] += 1

        spinner = threading.Thread(target=spin, daemon=True)
        spinner.start()
        start = ticks[0]
        time.sleep(0.2)
        free_rate = (ticks[0] - start) / 0.2
        comparison = native.RMSDComparison(mols)
        storage = native.DenseStorage(len(mols))
        options = native.PDistOptions()
        options.num_threads = 2
        options.chunk_size = 64
        before = ticks[0]
        began = time.perf_counter()
        native.pdist(comparison, storage, options)
        seconds = time.perf_counter() - began
        after = ticks[0]
        running[0] = False
        return {"fraction": ((after - before) / seconds) / free_rate,
                "seconds": seconds}


    # case -> (prepare inputs, run on them)
    CASES = {
        "rmsd": (lambda: molecules(oechem.OEMol, 300, one_graph=True, three_d=True),
                 lambda mols: pdist_digest(mols, "rmsd")),
        "mcs-graph": (lambda: molecules(oechem.OEGraphMol, 300),
                      lambda mols: pdist_digest(mols, "mcs")),
        "mcs-mol": (lambda: molecules(oechem.OEMol, 300),
                    lambda mols: pdist_digest(mols, "mcs")),
        "rocs-graph": (lambda: molecules(oechem.OEGraphMol, 60, three_d=True),
                       lambda mols: pdist_digest(mols, "rocs")),
        "rocs-mol": (lambda: molecules(oechem.OEMol, 60, three_d=True),
                     lambda mols: pdist_digest(mols, "rocs")),
        "descriptor": (lambda: molecules(oechem.OEGraphMol, 2000),
                       lambda mols: pdist_digest(mols, "descriptor")),
        "descriptor-build": (lambda: molecules(oechem.OEGraphMol, 2000),
                             lambda mols: oecluster.descriptor_statistics(mols)),
        "murcko": (lambda: molecules(oechem.OEGraphMol, 5000), murcko_digest),
        "two-at-once": (lambda: molecules(oechem.OEMol, 150, one_graph=True, three_d=True),
                        two_at_once),
        "progress-fingerprint": (lambda: molecules(oechem.OEGraphMol, 200),
                                 lambda mols: with_progress(mols, "fingerprint")),
        "progress-rmsd": (lambda: molecules(oechem.OEMol, 120, one_graph=True, three_d=True),
                          lambda mols: with_progress(mols, "rmsd")),
        "gil": (lambda: molecules(oechem.OEMol, 400, one_graph=True, three_d=True),
                gil_fraction),
    }
    prepare, execute = CASES[case]
''')

_RUN = '''
inputs = prepare()
print(json.dumps(execute(inputs)))
'''

# The parent samples this child's thread count between its ready and done
# lines. A sampler inside the child could not run while an unsafe-mode call
# keeps the GIL, so it would pass without measuring anything.
_READY = "__oecluster_sampled_ready__"
_DONE = "__oecluster_sampled_done__"

_SAMPLED = f'''
inputs = prepare()
print("{_READY}", flush=True)
sys.stdin.readline()
execute(inputs)
print("{_DONE}", flush=True)
'''


def _run(mode, case, timeout=900):
    result = subprocess.run(
        [sys.executable, "-X", "faulthandler", "-c", _PRELUDE + _RUN,
         str(mode), case],
        capture_output=True, text=True, timeout=timeout, check=False)
    assert result.returncode == 0, (
        mode, case, result.returncode, result.stderr[-3000:])
    return json.loads(result.stdout.strip().splitlines()[-1])


def _thread_counter():
    """A function from a process id to its OS thread count, or None where
    this platform offers no simple way to count."""
    if sys.platform.startswith("linux"):
        return lambda pid: len(os.listdir(f"/proc/{pid}/task"))
    if sys.platform == "darwin":
        # libSystem (find_library("c")) re-exports proc_pidinfo from libproc;
        # the C1-8 probe called it this way.
        libc = ctypes.CDLL(ctypes.util.find_library("c"), use_errno=True)

        class TaskInfo(ctypes.Structure):
            _fields_ = [("virtual", ctypes.c_uint64), ("resident", ctypes.c_uint64),
                        ("tuser", ctypes.c_uint64), ("tsys", ctypes.c_uint64),
                        ("thuser", ctypes.c_uint64), ("thsys", ctypes.c_uint64),
                        ("policy", ctypes.c_int32), ("faults", ctypes.c_int32),
                        ("pageins", ctypes.c_int32), ("cow", ctypes.c_int32),
                        ("msent", ctypes.c_int32), ("mrecv", ctypes.c_int32),
                        ("smach", ctypes.c_int32), ("sunix", ctypes.c_int32),
                        ("csw", ctypes.c_int32), ("threadnum", ctypes.c_int32),
                        ("numrunning", ctypes.c_int32), ("priority", ctypes.c_int32)]

        info = TaskInfo()
        libc.proc_pidinfo.argtypes = [ctypes.c_int, ctypes.c_int, ctypes.c_uint64,
                                      ctypes.POINTER(TaskInfo), ctypes.c_int]
        libc.proc_pidinfo.restype = ctypes.c_int

        def count(pid):
            # PROC_PIDTASKINFO is 4.
            size = libc.proc_pidinfo(pid, 4, ctypes.c_uint64(0),
                                     ctypes.byref(info), ctypes.sizeof(info))
            if size != ctypes.sizeof(info):
                raise OSError(ctypes.get_errno(), "proc_pidinfo")
            return info.threadnum

        return count
    return None


_COUNT_THREADS = _thread_counter()


def _peak_workers(mode, case, deadline=600):
    """Threads above the ready-time baseline at the busiest moment of the
    child's call. Both of the child's streams are drained by their own
    threads and one absolute deadline covers the whole run, so neither a
    missing sentinel nor a full pipe can block the parent; diagnostics are
    read only after the child has gone."""
    process = subprocess.Popen(
        [sys.executable, "-X", "faulthandler", "-c", _PRELUDE + _SAMPLED,
         str(mode), case],
        stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        text=True)
    lines = queue.Queue()
    errors = []

    def forward_out():
        for line in process.stdout:  # pyright: ignore[reportOptionalIterable]
            lines.put(line.strip())
        lines.put(None)

    def drain_err():
        errors.append(process.stderr.read())  # pyright: ignore[reportOptionalMemberAccess]

    out_reader = threading.Thread(target=forward_out, daemon=True)
    err_reader = threading.Thread(target=drain_err, daemon=True)
    out_reader.start()
    err_reader.start()
    finish = time.monotonic() + deadline
    failure = None
    baseline = peak = 0
    try:
        while True:
            remaining = finish - time.monotonic()
            if remaining <= 0:
                failure = f"{case}: no {_READY} line within {deadline} s"
                break
            try:
                line = lines.get(timeout=remaining)
            except queue.Empty:
                continue
            if line is None:
                failure = f"{case}: the child exited before it was ready"
                break
            if line == _READY:
                break
        if failure is None:
            for _ in range(50):
                baseline = max(baseline, _COUNT_THREADS(process.pid))
                time.sleep(0.001)
            process.stdin.write("go\n")  # pyright: ignore[reportOptionalMemberAccess]
            process.stdin.flush()  # pyright: ignore[reportOptionalMemberAccess]
            peak = baseline
            done = False
            while not done and failure is None:
                if time.monotonic() > finish:
                    failure = f"{case}: no {_DONE} line within {deadline} s"
                    break
                try:
                    peak = max(peak, _COUNT_THREADS(process.pid))
                except OSError:
                    pass  # the child may be exiting; the queue decides
                while True:
                    try:
                        line = lines.get_nowait()
                    except queue.Empty:
                        break
                    if line == _DONE:
                        done = True
                    elif line is None:
                        failure = f"{case}: the child exited before {_DONE}"
                        break
                time.sleep(0.0002)
        if failure is None:
            # The done line can arrive before the interpreter has shut down:
            # wait for a graceful exit within what is left of the deadline,
            # never kill a child that is finishing normally.
            remaining = finish - time.monotonic()
            if remaining <= 0:
                failure = f"{case}: the child did not exit within {deadline} s"
            else:
                try:
                    process.wait(timeout=remaining)
                except subprocess.TimeoutExpired:
                    failure = f"{case}: the child did not exit within {deadline} s"
    finally:
        if process.poll() is None:
            process.kill()
        process.wait()
        out_reader.join(timeout=30)
        err_reader.join(timeout=30)
    diagnostics = "".join(errors)[-3000:]
    assert failure is None, (failure, diagnostics)
    assert process.returncode == 0, diagnostics
    return peak - baseline


@pytest.mark.parametrize("case", [
    "rmsd", "mcs-graph", "mcs-mol", "rocs-graph", "rocs-mol", "descriptor",
    "murcko"])
def test_an_unsafe_pool_gives_the_default_pool_s_answer(case):
    assert _run(UNSAFE, case) == _run(SAFE, case)


def test_two_calls_at_once_under_an_unsafe_pool():
    unsafe = _run(UNSAFE, "two-at-once")
    assert unsafe[0] == unsafe[1]
    assert unsafe == _run(SAFE, "two-at-once")


@pytest.mark.parametrize("case", ["progress-fingerprint", "progress-rmsd"])
def test_a_progress_callback_finishes_under_an_unsafe_pool(case):
    unsafe = _run(UNSAFE, case, timeout=300)
    assert unsafe["called"]
    assert unsafe == _run(SAFE, case)


def test_an_unsafe_pool_keeps_the_gil_for_the_call():
    unsafe = _run(UNSAFE, "gil")
    safe = _run(SAFE, "gil")
    assert safe["seconds"] >= 1.0, safe
    assert unsafe["fraction"] < 0.05, unsafe
    assert safe["fraction"] > 0.2, safe


@pytest.mark.skipif(_COUNT_THREADS is None,
                    reason="counts OS threads on Linux and macOS only")
def test_murcko_runs_one_worker_under_an_unsafe_pool():
    assert _peak_workers(SAFE, "murcko") >= 2
    assert _peak_workers(UNSAFE, "murcko") <= 1


@pytest.mark.skipif(_COUNT_THREADS is None,
                    reason="counts OS threads on Linux and macOS only")
def test_the_rocs_case_runs_in_parallel_under_the_default_pool():
    """Without this, the unsafe ROCS cases could pass by never running in
    parallel at all."""
    assert _peak_workers(SAFE, "rocs-graph") >= 2


@pytest.mark.skipif(_COUNT_THREADS is None,
                    reason="counts OS threads on Linux and macOS only")
def test_descriptor_calculation_follows_the_pool_mode():
    """The descriptor crash case cannot tell a working gate from one
    that is always serial; the thread count can."""
    assert _peak_workers(SAFE, "descriptor-build") >= 2
    assert _peak_workers(UNSAFE, "descriptor-build") <= 1


def test_runs_openeye_code_is_not_python_api():
    from oecluster import oecluster as native
    assert not hasattr(native.PairwiseComparison, "RunsOpenEyeCode")
