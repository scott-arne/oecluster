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
                         ids=["unbounded-cache", "mutexed-cache", "top-of-range"])
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
