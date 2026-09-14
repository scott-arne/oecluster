"""Partition-agreement metrics: parity, divergences and the Python surface."""

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
