"""Integer arguments are coerced with operator.index before any comparison.

Each contract names one integer parameter of a public entry point. A float or
a numeric string must raise TypeError, a numpy integer must behave as the
equal int, and an object usable only through ``__index__`` must be accepted,
which proves the coercion runs before any comparison or ``int()``.
"""

import numpy as np
import oecluster
import pytest


def _line_matrix():
    """Six items on a line, in two groups of three."""
    positions = np.array([0.0, 1.0, 2.0, 10.0, 11.0, 12.0])
    square = np.abs(positions[:, None] - positions[None, :])
    return oecluster.SymmetricDistanceMatrix.from_condensed(square)


_DM = _line_matrix()


class _IndexOnly:
    """An integer through ``__index__`` only; ``int()`` and comparisons fail."""

    def __init__(self, value):
        self._value = value

    def __index__(self):
        return self._value

    def __int__(self):
        raise AssertionError("int() reached a value meant for operator.index")

    def _compared(self, other):
        raise AssertionError(
            "a comparison reached a value meant for operator.index")

    __lt__ = __le__ = __gt__ = __ge__ = __eq__ = __ne__ = _compared


def _outcome(result):
    """A comparable summary of any entry point's result."""
    if hasattr(result, "labels"):
        return [int(label) for label in result.labels]
    if isinstance(result, tuple):
        return [representative.member for representative in result]
    return result


# name: (call taking the parameter's value, a legal value, a float at or above
# the minimum, a float below the minimum)
_CONTRACTS = {
    "dbscan.min_samples": (
        lambda v: oecluster.dbscan(_DM, 1.5, min_samples=v), 2, 2.5, 0.5),
    "hdbscan.min_cluster_size": (
        lambda v: oecluster.hdbscan(_DM, min_cluster_size=v), 2, 2.5, 1.5),
    "hdbscan.min_samples": (
        lambda v: oecluster.hdbscan(_DM, min_cluster_size=2, min_samples=v),
        1, 1.5, 0.5),
    "hdbscan.max_cluster_size": (
        lambda v: oecluster.hdbscan(_DM, min_cluster_size=2,
                                    max_cluster_size=v),
        6, 6.5, 0.5),
    "agglomerative.n_clusters": (
        lambda v: oecluster.agglomerative(_DM, n_clusters=v), 2, 2.5, 0.5),
    "representative.cluster": (
        lambda v: oecluster.representative([v, 1, 2], _DM), 0, 0.5, -0.5),
    "rank_representatives.cluster": (
        lambda v: oecluster.rank_representatives([v, 1, 2], _DM),
        0, 0.5, -0.5),
    "select_representatives.cluster": (
        lambda v: oecluster.select_representatives([v, 1, 2], _DM, k=1),
        0, 0.5, -0.5),
    "select_representatives.k": (
        lambda v: oecluster.select_representatives([0, 1, 2], _DM, k=v),
        1, 1.5, -0.5),
}


@pytest.mark.parametrize("name", sorted(_CONTRACTS))
def test_a_float_at_or_above_the_minimum_is_a_type_error(name):
    call, _, above, _ = _CONTRACTS[name]
    with pytest.raises(TypeError):
        call(above)


@pytest.mark.parametrize("name", sorted(_CONTRACTS))
def test_a_float_below_the_minimum_is_a_type_error(name):
    call, _, _, below = _CONTRACTS[name]
    with pytest.raises(TypeError):
        call(below)


@pytest.mark.parametrize("name", sorted(_CONTRACTS))
def test_a_numeric_string_is_a_type_error(name):
    call, legal, _, _ = _CONTRACTS[name]
    with pytest.raises(TypeError):
        call(str(legal))


@pytest.mark.parametrize("name", sorted(_CONTRACTS))
def test_a_numpy_integer_behaves_as_the_equal_int(name):
    call, legal, _, _ = _CONTRACTS[name]
    assert _outcome(call(np.int64(legal))) == _outcome(call(legal))


@pytest.mark.parametrize("name", sorted(_CONTRACTS))
def test_an_index_only_value_is_accepted(name):
    call, legal, _, _ = _CONTRACTS[name]
    assert _outcome(call(_IndexOnly(legal))) == _outcome(call(legal))


def test_hdbscan_still_takes_none_for_min_samples_and_max_cluster_size():
    oecluster.hdbscan(_DM, min_cluster_size=2, min_samples=None,
                      max_cluster_size=None)


def test_hdbscan_still_reports_min_cluster_size_before_max_cluster_size():
    with pytest.raises(ValueError, match="min_cluster_size must be at least two"):
        oecluster.hdbscan(_DM, min_cluster_size=1, max_cluster_size=2.5)


def test_dbscan_still_reports_eps_before_min_samples():
    with pytest.raises(ValueError, match="eps must be non-negative"):
        oecluster.dbscan(_DM, -1.0, min_samples=2.5)


@pytest.mark.parametrize("value", [2.5, "2"], ids=["float", "str"])
def test_agglomerative_refuses_a_non_integer_n_clusters_with_a_threshold(value):
    with pytest.raises(TypeError):
        oecluster.agglomerative(_DM, n_clusters=value, distance_threshold=3.0)


def test_agglomerative_accepts_an_index_only_n_clusters_with_a_threshold():
    expected = oecluster.agglomerative(_DM, n_clusters=2,
                                       distance_threshold=3.0)
    got = oecluster.agglomerative(_DM, n_clusters=_IndexOnly(2),
                                  distance_threshold=3.0)
    assert _outcome(got) == _outcome(expected)


def test_agglomerative_reports_a_bad_threshold_before_a_float_n_clusters():
    with pytest.raises(ValueError, match="distance_threshold must be non-negative"):
        oecluster.agglomerative(_DM, n_clusters=2.5, distance_threshold=-1.0)


def test_agglomerative_reports_a_bad_linkage_before_a_float_n_clusters():
    with pytest.raises(ValueError, match="Unknown agglomerative linkage"):
        oecluster.agglomerative(_DM, n_clusters=2.5, distance_threshold=3.0,
                                linkage="bogus")
