"""Shared Python clustering result behavior."""

import numpy as np


def test_clustering_result_normalizes_labels_and_clusters():
    import oecluster

    result = oecluster.ClusteringResult(
        labels=[0, -1, 1, 0],
        clusters=((0, 3), (2,)),
    )

    assert isinstance(result.labels, np.ndarray)
    assert result.labels.dtype == np.intp
    assert result.labels.tolist() == [0, -1, 1, 0]
    assert result.clusters == ((0, 3), (2,))
    assert result.core_sample_indices == ()
    assert result.probabilities is None
