"""Python surface of consensus clustering."""
import numpy as np
import oecluster
import pytest

# --- fixtures ----------------------------------------------------------------


def _blobs20():
    """Twenty points, two blobs of ten: intra-blob 0.1, inter-blob 0.9."""
    values = [0.1 if (i < 10) == (j < 10) else 0.9
              for i in range(20) for j in range(i + 1, 20)]
    return oecluster.SymmetricDistanceMatrix.from_condensed(np.array(values))


def _result(labels):
    """A ClusteringResult from labels alone; -1 is noise."""
    labels = np.asarray(labels)
    clusters = [tuple(int(p) for p in np.flatnonzero(labels == label))
                for label in sorted(set(labels.tolist())) if label >= 0]
    return oecluster.ClusteringResult(labels, clusters)


# --- native bindings ---------------------------------------------------------

def test_native_coassociation_distances_fills_a_destination():
    native = oecluster.oecluster
    destination = oecluster.DenseStorage(4)
    summary = native.coassociation_distances(
        4,
        native.SizeTVector([0, 4, 8]),
        native.SizeTVector([0, 1, 2, 3, 0, 1, 2, 3]),
        native.IntVector([0, 0, 1, 1, 0, 0, 0, 1]),
        destination)
    assert summary.num_partitions == 2
    assert summary.unobserved_pairs == 0
    assert destination.Get(0, 1) == 0.0
    assert destination.Get(0, 2) == 0.5
    assert destination.Get(0, 3) == 1.0


def test_native_consensus_components_and_strength_round_trip():
    native = oecluster.oecluster
    destination = oecluster.DenseStorage(4)
    native.coassociation_distances(
        4,
        native.SizeTVector([0, 4]),
        native.SizeTVector([0, 1, 2, 3]),
        native.IntVector([0, 0, 1, 1]),
        destination)
    components = native.consensus_components(destination, 0.5)
    assert isinstance(components, native.IntVector)
    labels = list(components)
    assert labels == [0, 0, 1, 1]
    strength = native.consensus_strength(
        destination, native.IntVector(labels))
    assert list(strength.item_consensus) == [1.0, 1.0, 1.0, 1.0]
    assert list(strength.cluster_consensus) == [1.0, 1.0]


def test_native_refusal_is_a_runtime_error():
    native = oecluster.oecluster
    with pytest.raises(RuntimeError, match="more than once"):
        native.coassociation_distances(
            4,
            native.SizeTVector([0, 2]),
            native.SizeTVector([1, 1]),
            native.IntVector([0, 0]),
            oecluster.DenseStorage(4))


def test_native_strength_members_survive_a_temporary():
    native = oecluster.oecluster
    destination = oecluster.DenseStorage(4)
    native.coassociation_distances(
        4,
        native.SizeTVector([0, 4]),
        native.SizeTVector([0, 1, 2, 3]),
        native.IntVector([0, 0, 1, 1]),
        destination)
    labels = native.IntVector([0, 0, 1, 1])
    assert list(native.consensus_strength(
        destination, labels).item_consensus) == [1.0, 1.0, 1.0, 1.0]
    assert list(native.consensus_strength(
        destination, labels).cluster_consensus) == [1.0, 1.0]
