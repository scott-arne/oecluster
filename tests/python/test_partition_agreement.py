"""Partition-agreement metrics: parity, divergences and the Python surface."""

import numpy as np
import oecluster
import pytest


def test_partition_agreement_matches_reference_values():
    agreement = oecluster.partition_agreement(
        [0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2],
        [0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3],
    )
    assert agreement.num_samples == 12
    assert agreement.num_clusters_a == 3
    assert agreement.num_clusters_b == 4
    assert agreement.adjusted_rand_index == pytest.approx(0.40310077519379844)
    assert agreement.v_measure == agreement.normalized_mutual_information
    assert agreement.requested.adjusted_mutual_information is False
    assert isinstance(agreement.requested, oecluster.PartitionAgreementRequested)


def test_partition_agreement_rejects_float_labels():
    with pytest.raises(TypeError):
        oecluster.partition_agreement([0.1, 0.9, 0.1], [0, 1, 0])


def test_scaffold_agreement_rejects_non_string_scaffolds():
    with pytest.raises(TypeError):
        oecluster.scaffold_agreement([0, 0, 1], [None, None, "X"])


def test_partition_agreement_accepts_numpy_intp_labels():
    # ClusteringResult.labels is np.intp, which must keep working
    labels_a = np.array([0, 0, 1, 1], dtype=np.intp)
    labels_b = np.array([0, 1, 2, 2], dtype=np.intp)
    agreement = oecluster.partition_agreement(labels_a, labels_b)
    assert agreement.num_samples == 4
