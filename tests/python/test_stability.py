"""Python surface of take and cluster_stability."""
import numpy as np
import oecluster
import oefp
import pytest

# --- fixtures ----------------------------------------------------------------

FP_ROWS = ([set(range(8)) | {16 + i} for i in range(5)]
           + [set(range(32, 40)) | {48 + i} for i in range(5)])


def _batch(rows):
    return oefp.OEFPBatch.from_fingerprints(
        [oefp.OEFP.from_on_bits(64, sorted(on)) for on in rows])


def _fps():
    """Two groups of five sharing an eight-bit block, one private bit each."""
    return _batch(FP_ROWS)


# --- native bindings ---------------------------------------------------------

def test_native_take_pairs_gathers_in_index_order():
    source = oecluster.DenseStorage(4)
    for i in range(4):
        for j in range(i + 1, 4):
            source.Set(i, j, 10 * i + j)
    destination = oecluster.DenseStorage(3)
    oecluster.oecluster.take_pairs(
        source, oecluster.oecluster.SizeTVector([3, 0, 2]), destination)
    assert destination.Get(0, 1) == source.Get(3, 0)
    assert destination.Get(0, 2) == source.Get(3, 2)
    assert destination.Get(1, 2) == source.Get(0, 2)


def test_native_take_pairs_refusal_is_a_runtime_error():
    source = oecluster.DenseStorage(3)
    with pytest.raises(RuntimeError, match="more than once"):
        oecluster.oecluster.take_pairs(
            source, oecluster.oecluster.SizeTVector([1, 1]),
            oecluster.DenseStorage(2))


def test_native_take_fingerprints_returns_an_oefp_native_batch():
    batch = _fps()
    native = oecluster.oecluster.take_fingerprints(
        batch, oecluster.oecluster.SizeTVector([7, 0]))
    assert type(native).__module__ == "oefp._native"
    subset = oefp.OEFPBatch._from_native(native)
    np.testing.assert_array_equal(subset.words, batch.words[[7, 0]])
    assert subset.spec == batch.spec
