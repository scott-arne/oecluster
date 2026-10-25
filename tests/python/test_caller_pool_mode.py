"""The memory-pool mode oecluster reads from the caller's ``openeye.oechem``.

A statically linked extension has its own OpenEye copy, so the descriptor
entry points read the application's mode from ``openeye.oechem`` in
``sys.modules`` and keep OpenEye work on one thread when it is not
thread-safe. The read must fail closed: anything but a plain int within
``unsigned int``'s range answers 0, an unsafe mode, because a call wrongly
kept serial is slow while one wrongly run in parallel crashes. A bool would
otherwise read as 1 (Mutexed) and a value past 2**32 would truncate.

``_caller_pool_mode`` exposes the reader for these tests. Each
test swaps the ``sys.modules`` entry through ``monkeypatch``, which restores
the real module afterwards. A call that returns normally also shows that the
reader left no Python error set, since the wrapper would otherwise raise.

The last two tests raise from an entry point under each GIL macro while the
stub reports an unsafe mode, so the wrapper keeps the GIL and its catch
clause must not restore a thread state it never saved. The stub changes only
what oecluster reads; OpenEye's own pool keeps its default mode, so the calls
are safe in this process.
"""

import sys
import types

import pytest
from oecluster import oecluster as native

_KEY = "openeye.oechem"


class _Index:
    """Converts to 10 through ``__index__`` without being an int."""

    def __index__(self):
        return 10


class _Int(int):
    """An int subclass, which the reader must not take for a plain int."""


def _returning(value):
    module = types.ModuleType(_KEY)
    module.OEGetMemPoolMode = lambda: value  # pyright: ignore[reportAttributeAccessIssue]
    return module


def _raising():
    module = types.ModuleType(_KEY)

    def get():
        raise RuntimeError("OEGetMemPoolMode failed")

    module.OEGetMemPoolMode = get  # pyright: ignore[reportAttributeAccessIssue]
    return module


def _without_getter():
    return types.ModuleType(_KEY)


def test_the_real_module_s_mode_is_read():
    oechem = pytest.importorskip(_KEY)
    assert native._caller_pool_mode() == oechem.OEGetMemPoolMode()


def test_a_module_never_imported_gives_no_value(monkeypatch):
    monkeypatch.delitem(sys.modules, _KEY, raising=False)
    assert native._caller_pool_mode() is None


@pytest.mark.parametrize("mode", [8, 10, 2**32 - 1],
                         ids=["unbounded-cache", "thread-local-cache", "top-of-range"])
def test_a_plain_int_in_range_is_read_as_is(monkeypatch, mode):
    monkeypatch.setitem(sys.modules, _KEY, _returning(mode))
    assert native._caller_pool_mode() == mode


@pytest.mark.parametrize("make_module", [
    lambda: _returning(True),
    lambda: _returning(False),
    lambda: _returning(None),
    lambda: _returning(1.0),
    lambda: _returning("10"),
    lambda: _returning(_Index()),
    lambda: _returning(_Int(10)),
    lambda: _returning(-1),
    lambda: _returning(2**32),
    lambda: _returning(2**32 + 1),
    lambda: _returning(2**32 + 10),
    lambda: _returning(2**64 + 10),
    _raising,
    _without_getter,
], ids=["true", "false", "none", "float", "str", "index", "int-subclass",
        "negative", "2**32", "2**32+1", "2**32+10", "2**64+10", "raises",
        "no-getter"])
def test_a_failed_or_malformed_read_answers_the_unsafe_mode_0(monkeypatch, make_module):
    monkeypatch.setitem(sys.modules, _KEY, make_module())
    assert native._caller_pool_mode() == 0


def _keep_the_gil(monkeypatch):
    """Make the GIL macros read mode 8, UnboundedCache, and keep the GIL."""
    monkeypatch.setitem(sys.modules, _KEY, _returning(8))
    assert native._caller_pool_mode() == 8


def test_pdist_raises_and_recovers_under_an_unsafe_mode(monkeypatch):
    """OECLUSTER_GIL_EXCEPTION's OEClusterError clause, reached with the GIL
    kept, raises RuntimeError and leaves the interpreter usable."""
    oechem = pytest.importorskip(_KEY)
    mols = []
    for smiles in ("c1ccccc1", "Cc1ccccc1", "CCc1ccccc1", "C1CCCCC1"):
        mol = oechem.OEGraphMol()
        assert oechem.OESmilesToMol(mol, smiles)
        mols.append(mol)
    comparison = native.FingerprintComparison(mols, native.FingerprintOptions())
    expected = native.DenseStorage(len(mols))
    native.pdist(comparison, expected)

    _keep_the_gil(monkeypatch)
    with pytest.raises(RuntimeError, match="pdist storage size does not match"):
        native.pdist(comparison, native.DenseStorage(len(mols) - 1))

    storage = native.DenseStorage(len(mols))
    native.pdist(comparison, storage)
    pairs = [(i, j) for i in range(len(mols)) for j in range(i + 1, len(mols))]
    assert [storage.Get(i, j) for i, j in pairs] == [expected.Get(i, j) for i, j in pairs]


def test_butina_raises_and_recovers_under_an_unsafe_mode(monkeypatch):
    """OECLUSTER_GIL_MEMORY_EXCEPTION's std::exception clause, reached with the
    GIL kept, raises RuntimeError and leaves the interpreter usable."""
    storage = native.DenseStorage(4)
    for (i, j), distance in {(0, 1): 0.1, (0, 2): 0.9, (0, 3): 0.9,
                             (1, 2): 0.9, (1, 3): 0.9, (2, 3): 0.1}.items():
        storage.Set(i, j, distance)

    _keep_the_gil(monkeypatch)
    refused = native.ButinaOptions()
    refused.distance_threshold = -1.0
    with pytest.raises(RuntimeError, match="Butina distance threshold must be non-negative"):
        native.butina_cluster(storage, refused)

    options = native.ButinaOptions()
    options.distance_threshold = 0.5
    assert list(native.butina_cluster(storage, options).Labels()) == [1, 1, 0, 0]
