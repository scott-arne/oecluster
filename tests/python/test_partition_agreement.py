"""Partition-agreement metrics: parity, divergences and the Python surface."""

import math

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


MAIN_A = [0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2]
MAIN_B = [0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3]
NOISE_A = [0, 0, 0, -1, 1, 1, 1, -1, 2, 2, -1, 2]
NOISE_B = [0, 0, 1, 1, 1, -1, 2, 2, 2, 3, 3, -1]


def test_live_sklearn_parity_on_non_degenerate_fixtures():
    from sklearn.metrics import (
        adjusted_mutual_info_score,
        adjusted_rand_score,
        fowlkes_mallows_score,
        homogeneity_completeness_v_measure,
        normalized_mutual_info_score,
    )

    fixtures = [
        (MAIN_A, MAIN_B),
        ([0, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 5],
         [0, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4]),
        ([0, 0, 1, 1, 2, 2, 3, 3], [0, 1, 1, 2, 2, 3, 3, 0]),
    ]
    for a, b in fixtures:
        agreement = oecluster.partition_agreement(
            a, b, adjusted_mutual_information=True)
        homogeneity, completeness, v_measure = \
            homogeneity_completeness_v_measure(a, b)
        assert agreement.adjusted_rand_index == pytest.approx(
            adjusted_rand_score(a, b))
        assert agreement.fowlkes_mallows == pytest.approx(
            fowlkes_mallows_score(a, b))
        assert agreement.normalized_mutual_information == pytest.approx(
            normalized_mutual_info_score(a, b))
        assert agreement.homogeneity == pytest.approx(homogeneity)
        assert agreement.completeness == pytest.approx(completeness)
        assert agreement.v_measure == pytest.approx(v_measure)
        assert agreement.adjusted_mutual_information == pytest.approx(
            adjusted_mutual_info_score(a, b))


def test_rule_one_diverges_from_sklearn_on_a_single_sample():
    from sklearn.metrics import (
        adjusted_mutual_info_score,
        adjusted_rand_score,
        fowlkes_mallows_score,
        homogeneity_completeness_v_measure,
        normalized_mutual_info_score,
    )

    # AMI is requested so the divergence is asserted over all seven fields:
    # rule 1 is a whole-struct rule, and an AMI that leaked a value here would
    # be the one metric this test failed to notice.
    agreement = oecluster.partition_agreement(
        [0], [0], adjusted_mutual_information=True)
    assert agreement.num_samples == 1
    for name in ("adjusted_rand_index", "fowlkes_mallows",
                 "normalized_mutual_information", "homogeneity",
                 "completeness", "v_measure", "adjusted_mutual_information"):
        assert math.isnan(getattr(agreement, name)), name

    # scikit-learn calls a one-sample labeling perfect agreement -- on every
    # metric but Fowlkes-Mallows, which it reports as 0.0. None of these raise.
    assert adjusted_rand_score([0], [0]) == 1.0
    assert fowlkes_mallows_score([0], [0]) == 0.0
    assert homogeneity_completeness_v_measure([0], [0]) == (1.0, 1.0, 1.0)
    assert normalized_mutual_info_score([0], [0]) == 1.0
    assert adjusted_mutual_info_score([0], [0]) == 1.0


def test_single_cluster_side_a_diverges_on_homogeneity():
    from sklearn.metrics import homogeneity_completeness_v_measure

    a, b = [0, 0, 0, 0], [0, 0, 1, 1]
    agreement = oecluster.partition_agreement(a, b)
    homogeneity, completeness, v_measure = \
        homogeneity_completeness_v_measure(a, b)

    assert math.isnan(agreement.homogeneity)
    assert homogeneity == 1.0
    # The divergence is confined to the component whose entropy is zero. Pinning
    # the transpose component, which stays defined and agrees with
    # scikit-learn, is what makes a swapped entropy guard detectable.
    assert agreement.completeness == 0.0
    assert agreement.completeness == pytest.approx(completeness)
    # The composite agrees even though the component does not.
    assert agreement.v_measure == pytest.approx(v_measure)
    assert agreement.v_measure == 0.0


def test_single_cluster_side_b_diverges_on_completeness():
    from sklearn.metrics import homogeneity_completeness_v_measure

    a, b = [0, 0, 1, 1], [0, 0, 0, 0]
    agreement = oecluster.partition_agreement(a, b)
    homogeneity, completeness, _ = homogeneity_completeness_v_measure(a, b)

    assert math.isnan(agreement.completeness)
    assert completeness == 1.0
    # The divergence is confined to the component whose entropy is zero. Pinning
    # the transpose component, which stays defined and agrees with
    # scikit-learn, is what makes a swapped entropy guard detectable.
    assert agreement.homogeneity == 0.0
    assert agreement.homogeneity == pytest.approx(homogeneity)


def test_all_singletons_one_side_diverges_on_fowlkes_mallows():
    from sklearn.metrics import fowlkes_mallows_score

    a, b = [0, 1, 2, 3], [0, 0, 1, 1]
    agreement = oecluster.partition_agreement(a, b)

    assert math.isnan(agreement.fowlkes_mallows)
    assert fowlkes_mallows_score(a, b) == 0.0


def test_both_sides_all_singletons_are_identical_not_zero():
    from sklearn.metrics import fowlkes_mallows_score

    a, b = [0, 1, 2, 3], [3, 2, 1, 0]
    agreement = oecluster.partition_agreement(a, b)

    assert agreement.fowlkes_mallows == 1.0
    assert fowlkes_mallows_score(a, b) == 0.0


def test_independent_partitions_agree_with_sklearn_everywhere():
    from sklearn.metrics import (
        adjusted_mutual_info_score,
        adjusted_rand_score,
        fowlkes_mallows_score,
        normalized_mutual_info_score,
    )

    a, b = [0, 0, 1, 1], [0, 1, 0, 1]
    agreement = oecluster.partition_agreement(
        a, b, adjusted_mutual_information=True)

    assert agreement.adjusted_rand_index == pytest.approx(
        adjusted_rand_score(a, b))
    assert agreement.fowlkes_mallows == pytest.approx(fowlkes_mallows_score(a, b))
    assert agreement.normalized_mutual_information == pytest.approx(
        normalized_mutual_info_score(a, b))
    assert agreement.adjusted_mutual_information == pytest.approx(
        adjusted_mutual_info_score(a, b))
    assert agreement.homogeneity == 0.0
    assert agreement.completeness == 0.0
    assert agreement.v_measure == 0.0


def test_the_three_input_forms_agree():
    # A real ClusteringResult, not a stand-in with a `.labels` attribute: the
    # spec asks for the clustering-result form, and `ClusteringResult` stores
    # its labels as an np.intp array, which a hand-rolled object would not
    # exercise.
    result_a = oecluster.ClusteringResult(
        MAIN_A, [[i for i, label in enumerate(MAIN_A) if label == k]
                 for k in (0, 1, 2)])
    result_b = oecluster.ClusteringResult(
        MAIN_B, [[i for i, label in enumerate(MAIN_B) if label == k]
                 for k in (0, 1, 2, 3)])

    def agreement(a, b):
        return oecluster.partition_agreement(
            a, b, adjusted_mutual_information=True)

    as_list = agreement(MAIN_A, MAIN_B)
    as_tuple = agreement(tuple(MAIN_A), tuple(MAIN_B))
    as_array = agreement(np.array(MAIN_A, dtype=np.int32),
                         np.array(MAIN_B, dtype=np.int64))
    as_result_a = agreement(result_a, MAIN_B)
    as_result_b = agreement(MAIN_A, result_b)

    # "Identical" is every field, not a sample of them, so compare the whole
    # table plus `requested`. AMI is requested so the table's one opt-in row
    # carries a value rather than None in all four forms. The rows cannot be
    # compared with a bare == because a NaN cell never equals itself.
    for other in (as_tuple, as_array, as_result_a, as_result_b):
        assert other.requested == as_list.requested
        for (name, value), (expected_name, expected) in zip(
                other.to_table(), as_list.to_table(), strict=True):
            assert name == expected_name
            if isinstance(expected, float) and math.isnan(expected):
                assert math.isnan(value), name
            else:
                assert value == expected, name


def test_noise_strings_reach_the_right_enum():
    singletons = oecluster.partition_agreement(NOISE_A, NOISE_B)
    grouped = oecluster.partition_agreement(NOISE_A, NOISE_B, noise="grouped")
    excluded = oecluster.partition_agreement(NOISE_A, NOISE_B, noise="excluded")

    assert (singletons.num_samples, singletons.num_clusters_a,
            singletons.num_clusters_b) == (12, 6, 6)
    assert (grouped.num_samples, grouped.num_clusters_a,
            grouped.num_clusters_b) == (12, 4, 5)
    assert (excluded.num_samples, excluded.num_clusters_a,
            excluded.num_clusters_b) == (7, 3, 4)

    assert singletons.adjusted_rand_index == pytest.approx(
        -0.012269938650306749)
    assert grouped.adjusted_rand_index == pytest.approx(-0.07179487179487179)
    assert excluded.adjusted_rand_index == pytest.approx(0.08695652173913043)


def test_unknown_noise_string_names_the_three_valid_values():
    with pytest.raises(ValueError) as excinfo:
        oecluster.partition_agreement(MAIN_A, MAIN_B, noise="drop")
    message = str(excinfo.value)
    assert "singletons" in message
    assert "grouped" in message
    assert "excluded" in message


def test_non_sequence_arguments_raise_type_error_naming_the_argument():
    with pytest.raises(TypeError) as excinfo:
        oecluster.partition_agreement(3, MAIN_B)
    assert str(excinfo.value).startswith("a ")

    with pytest.raises(TypeError) as excinfo:
        oecluster.partition_agreement(MAIN_A, object())
    assert str(excinfo.value).startswith("b ")

    with pytest.raises(TypeError) as excinfo:
        oecluster.scaffold_agreement([0, 0, 1], "abc")
    assert "scaffold_labels" in str(excinfo.value)


def test_out_of_range_label_raises_value_error_naming_the_argument():
    # 2**40 is a perfectly good int, so the rejection comes from the native
    # vector<int> rather than from the type check, and the raw OverflowError
    # would name an internal container type instead of the real constraint.
    with pytest.raises(ValueError) as excinfo:
        oecluster.partition_agreement([0, 2**40, 1], [0, 1, 1])
    message = str(excinfo.value)
    assert message.startswith("a ")
    assert "32-bit" in message

    with pytest.raises(ValueError) as excinfo:
        oecluster.scaffold_agreement([0, 2**40, 1], ["x", "y", "z"])
    message = str(excinfo.value)
    assert "result" in message
    assert "32-bit" in message


def test_unencodable_scaffold_raises_type_error_naming_the_argument():
    # A lone surrogate is a str, so it clears the isinstance check and is
    # rejected only when SWIG encodes it. surrogateescape decoding of a
    # mis-encoded scaffold file is how one reaches a caller in practice.
    with pytest.raises(TypeError) as excinfo:
        oecluster.scaffold_agreement([0, 1, 1], ["a", "\ud800", "c"])
    message = str(excinfo.value)
    assert "scaffold_labels" in message
    assert "StringVector" not in message


def test_mapping_arguments_raise_type_error_naming_the_argument():
    # A dict is iterable, so without the guard these score its keys and return
    # a plausible number for a labeling nobody passed. The keys here are valid
    # labels and the lengths line up, so nothing downstream would object.
    labels = dict(enumerate(MAIN_B))

    with pytest.raises(TypeError) as excinfo:
        oecluster.partition_agreement(labels, MAIN_B)
    message = str(excinfo.value)
    assert message.startswith("a ")
    assert "mapping" in message

    with pytest.raises(TypeError) as excinfo:
        oecluster.partition_agreement(MAIN_A, labels)
    assert str(excinfo.value).startswith("b ")

    with pytest.raises(TypeError) as excinfo:
        oecluster.scaffold_agreement([0, 1, 1], {"x": 1, "y": 2, "z": 3})
    message = str(excinfo.value)
    assert "scaffold_labels" in message
    assert "mapping" in message


def test_string_flag_is_refused_rather_than_read_as_true():
    # bool("false") is True, so the string a caller most plausibly means as
    # "off" would have switched AMI on and reported it as requested.
    # Pyright infers the parameter type from its False default, so the three
    # off-type arguments below are exactly what it flags -- and exactly what
    # this test exists to pass at runtime.
    with pytest.raises(TypeError) as excinfo:
        oecluster.partition_agreement(
            MAIN_A, MAIN_B,
            adjusted_mutual_information="false")  # pyright: ignore[reportArgumentType]
    message = str(excinfo.value)
    assert "adjusted_mutual_information" in message
    assert "str" in message

    # Ints and numpy bools stay admissible: this rejects strings only.
    assert oecluster.partition_agreement(
        MAIN_A, MAIN_B,
        adjusted_mutual_information=1  # pyright: ignore[reportArgumentType]
    ).requested.adjusted_mutual_information
    assert not oecluster.partition_agreement(
        MAIN_A, MAIN_B,
        adjusted_mutual_information=np.bool_(False)  # pyright: ignore[reportArgumentType]
    ).requested.adjusted_mutual_information


def test_validation_surfaces_as_value_error_not_runtime_error():
    # SWIG maps every std::exception to RuntimeError, so these only pass while
    # the Python-side checks run ahead of the native call.
    with pytest.raises(ValueError) as excinfo:
        oecluster.scaffold_agreement([0, 0, 1], ["a", "b"])
    message = str(excinfo.value)
    assert "2" in message and "3" in message and "scaffold_labels" in message

    with pytest.raises(ValueError) as excinfo:
        oecluster.partition_agreement([0, 0, 1], [0, 1])
    message = str(excinfo.value)
    assert "2" in message and "3" in message

    with pytest.raises(ValueError):
        oecluster.partition_agreement([], [])
    with pytest.raises(ValueError):
        oecluster.scaffold_agreement([], [])
    with pytest.raises(ValueError):
        oecluster.scaffold_agreement([0, 1], [])


def test_scaffold_agreement_treats_empty_strings_as_missing():
    labels = [0, 0, 0, 1, 1, 1, 2, 2, 2]
    scaffolds = ["ar", "ar", "ar", "pi", "", "al", "al", "al", ""]

    singletons = oecluster.scaffold_agreement(labels, scaffolds)
    grouped = oecluster.scaffold_agreement(labels, scaffolds, noise="grouped")
    excluded = oecluster.scaffold_agreement(labels, scaffolds, noise="excluded")

    assert (singletons.num_samples, singletons.num_clusters_b) == (9, 5)
    assert (grouped.num_samples, grouped.num_clusters_b) == (9, 4)
    assert (excluded.num_samples, excluded.num_clusters_b) == (7, 3)
    assert singletons.adjusted_rand_index == pytest.approx(0.4166666666666667)
    assert grouped.adjusted_rand_index == pytest.approx(0.36)
    assert excluded.adjusted_rand_index == pytest.approx(0.631578947368421)


def test_scaffold_agreement_matches_sklearn_on_interned_strings():
    from sklearn.metrics import (
        adjusted_mutual_info_score,
        adjusted_rand_score,
        completeness_score,
        fowlkes_mallows_score,
        homogeneity_score,
        normalized_mutual_info_score,
        v_measure_score,
    )

    labels = [0, 0, 0, 1, 1, 1, 2, 2, 2]
    scaffolds = ["ar", "ar", "ar", "pi", "pi", "al", "al", "al", "al"]

    # Default call: AMI not requested. Pin the full set of reported metrics.
    agreement = oecluster.scaffold_agreement(labels, scaffolds)
    assert agreement.adjusted_rand_index == pytest.approx(
        adjusted_rand_score(labels, scaffolds))
    assert agreement.fowlkes_mallows == pytest.approx(
        fowlkes_mallows_score(labels, scaffolds))
    assert agreement.normalized_mutual_information == pytest.approx(
        normalized_mutual_info_score(labels, scaffolds))
    assert agreement.homogeneity == pytest.approx(
        homogeneity_score(labels, scaffolds))
    assert agreement.completeness == pytest.approx(
        completeness_score(labels, scaffolds))
    assert agreement.v_measure == pytest.approx(
        v_measure_score(labels, scaffolds))
    assert agreement.requested.adjusted_mutual_information is False
    assert dict(agreement.to_table())["adjusted_mutual_information"] is None

    # AMI requested.
    with_ami = oecluster.scaffold_agreement(
        labels, scaffolds, adjusted_mutual_information=True)
    assert with_ami.requested.adjusted_mutual_information is True
    assert with_ami.adjusted_mutual_information == pytest.approx(
        adjusted_mutual_info_score(labels, scaffolds))

    # ClusteringResult as first argument.
    result = oecluster.ClusteringResult(
        labels, [[i for i, label in enumerate(labels) if label == k]
                 for k in (0, 1, 2)])
    as_result = oecluster.scaffold_agreement(result, scaffolds)
    assert as_result.requested == agreement.requested
    for (name, value), (expected_name, expected) in zip(
            as_result.to_table(), agreement.to_table(), strict=True):
        assert name == expected_name
        if isinstance(expected, float) and math.isnan(expected):
            assert math.isnan(value), name
        else:
            assert value == expected, name


def test_unasked_ami_renders_as_none_and_asked_undefined_as_nan():
    unasked = oecluster.partition_agreement(MAIN_A, MAIN_B)
    rows = dict(unasked.to_table())
    assert rows["adjusted_mutual_information"] is None
    repr_lines = repr(unasked).split('\n')
    ami_line = next(line for line in repr_lines if line.strip().startswith('adjusted_mutual_information'))
    assert "--" in ami_line
    # AMI is the only row rendering -- since every other metric has a real value.
    assert repr(unasked).count('--') == 1

    # Asked, but rule 1 leaves it undefined.
    asked = oecluster.partition_agreement(
        [0], [0], adjusted_mutual_information=True)
    rows = dict(asked.to_table())
    assert math.isnan(rows["adjusted_mutual_information"])
    repr_lines = repr(asked).split('\n')
    ami_line = next(line for line in repr_lines if line.strip().startswith('adjusted_mutual_information'))
    assert "nan" in ami_line


def test_to_table_covers_every_reported_field_in_struct_order():
    from sklearn.metrics import (
        adjusted_mutual_info_score,
        adjusted_rand_score,
        completeness_score,
        fowlkes_mallows_score,
        homogeneity_score,
        normalized_mutual_info_score,
        v_measure_score,
    )

    agreement = oecluster.partition_agreement(
        MAIN_A, MAIN_B, adjusted_mutual_information=True)

    # Struct-order check is load-bearing.
    assert [name for name, _ in agreement.to_table()] == [
        "num_samples",
        "num_clusters_a",
        "num_clusters_b",
        "adjusted_rand_index",
        "fowlkes_mallows",
        "normalized_mutual_information",
        "homogeneity",
        "completeness",
        "v_measure",
        "adjusted_mutual_information",
    ]

    # Pin the row values against live scikit-learn.
    rows = dict(agreement.to_table())
    assert rows["num_samples"] == 12
    assert rows["num_clusters_a"] == 3
    assert rows["num_clusters_b"] == 4
    assert rows["adjusted_rand_index"] == pytest.approx(
        adjusted_rand_score(MAIN_A, MAIN_B))
    assert rows["fowlkes_mallows"] == pytest.approx(
        fowlkes_mallows_score(MAIN_A, MAIN_B))
    assert rows["normalized_mutual_information"] == pytest.approx(
        normalized_mutual_info_score(MAIN_A, MAIN_B))
    assert rows["homogeneity"] == pytest.approx(
        homogeneity_score(MAIN_A, MAIN_B))
    assert rows["completeness"] == pytest.approx(
        completeness_score(MAIN_A, MAIN_B))
    assert rows["v_measure"] == pytest.approx(
        v_measure_score(MAIN_A, MAIN_B))
    assert rows["adjusted_mutual_information"] == pytest.approx(
        adjusted_mutual_info_score(MAIN_A, MAIN_B))

    # normalized_mutual_information and v_measure are identically equal by
    # construction: homogeneity is I/H(a), completeness is I/H(b), and their
    # harmonic mean reduces to 2I/(H(a)+H(b)), which is arithmetic-mean NMI.
    # No fixture can separate a swap between those two rows.
