"""Tests for the precomputed-distance ingress and the triangle probe."""

import numpy as np
import oecluster
import pytest
from oecluster import SymmetricDistanceMatrix, _gate


def _metric_condensed(n=6, seed=1):
    """A genuine Euclidean distance matrix, in condensed form."""
    rng = np.random.default_rng(seed)
    points = rng.normal(size=(n, 3))
    diff = points[:, None, :] - points[None, :, :]
    square = np.sqrt((diff ** 2).sum(axis=-1))
    return square[np.triu_indices(n, k=1)], square


def test_from_condensed_accepts_a_condensed_vector():
    condensed, _ = _metric_condensed()
    dm = SymmetricDistanceMatrix.from_condensed(condensed)
    assert dm.num_samples == 6
    assert dm.comparison_name == "precomputed"
    np.testing.assert_allclose(dm.condensed, condensed)


def test_from_condensed_accepts_a_square_matrix():
    condensed, square = _metric_condensed()
    dm = SymmetricDistanceMatrix.from_condensed(square)
    np.testing.assert_allclose(dm.condensed, condensed)


def test_from_condensed_honors_labels():
    condensed, _ = _metric_condensed(n=3)
    dm = SymmetricDistanceMatrix.from_condensed(
        condensed, labels=["a", "b", "c"], comparison_name="mine")
    assert dm.labels == ["a", "b", "c"]
    assert dm.comparison_name == "mine"


def test_from_condensed_rejects_a_non_triangular_length():
    with pytest.raises(ValueError, match="not a valid condensed length"):
        SymmetricDistanceMatrix.from_condensed(np.zeros(4))


def test_from_condensed_rejects_a_label_count_mismatch():
    condensed, _ = _metric_condensed(n=3)
    with pytest.raises(ValueError, match="labels"):
        SymmetricDistanceMatrix.from_condensed(condensed, labels=["a", "b"])


def test_from_condensed_rejects_non_finite_values():
    condensed, _ = _metric_condensed(n=3)
    condensed = condensed.copy()
    condensed[0] = np.nan
    with pytest.raises(ValueError, match="non-finite"):
        SymmetricDistanceMatrix.from_condensed(condensed)


def test_from_condensed_rejects_negative_values():
    condensed, _ = _metric_condensed(n=3)
    condensed = condensed.copy()
    condensed[1] = -0.5
    with pytest.raises(ValueError, match="negative"):
        SymmetricDistanceMatrix.from_condensed(condensed)


def test_from_condensed_rejects_an_asymmetric_square_matrix():
    _, square = _metric_condensed(n=4)
    square = square.copy()
    square[0, 1] += 1.0
    with pytest.raises(ValueError, match="symmetric"):
        SymmetricDistanceMatrix.from_condensed(square)


def test_from_condensed_rejects_a_non_zero_diagonal():
    _, square = _metric_condensed(n=4)
    square = square.copy()
    square[2, 2] = 0.5
    with pytest.raises(ValueError, match="diagonal"):
        SymmetricDistanceMatrix.from_condensed(square)


def test_from_condensed_rejects_a_non_square_matrix():
    with pytest.raises(ValueError, match="square"):
        SymmetricDistanceMatrix.from_condensed(np.zeros((3, 4)))


def test_from_condensed_rejects_a_three_dimensional_input():
    with pytest.raises(ValueError, match="1-D or 2-D"):
        SymmetricDistanceMatrix.from_condensed(np.zeros((2, 2, 2)))


def test_a_metric_matrix_probes_clean_and_clusters():
    """Both capabilities stay unknown; the probe carries the whole burden.

    There is no metric object to interrogate here, so stamping ``zero_self``
    True would assert something the caller never said. Unknown is permissive,
    so the matrix still clusters. ``is_distance`` is unknown for the same
    reason: the argument's name is not evidence about the numbers in it.
    """
    condensed, _ = _metric_condensed(n=12)
    dm = SymmetricDistanceMatrix.from_condensed(condensed)
    assert dm.metric_probe == "no_violations_found"
    assert dm.probe_violations == 0
    assert dm.probe_sampled > 0
    assert dm.is_distance == "unknown"
    assert dm.metric_capabilities == {'zero_self': "unknown",
                                      'triangle': "unknown"}
    assert dm.data_integrity == "complete"
    oecluster.butina(dm, 1.0)


def test_a_planted_violation_is_found_and_refused():
    condensed = np.full(6, 0.1)   # n = 4
    condensed[1] = 10.0           # d(0, 2) >> d(0, 1) + d(1, 2)
    dm = SymmetricDistanceMatrix.from_condensed(condensed, probe_triples=2000)
    assert dm.metric_probe == "violations_found"
    assert dm.probe_violations > 0
    with pytest.raises(ValueError, match="sampled triples"):
        oecluster.butina(dm, 0.5)
    oecluster.butina(dm, 0.5, allow_nonmetric=True)


def test_the_probe_can_be_disabled():
    condensed = np.full(6, 0.1)
    condensed[1] = 10.0
    dm = SymmetricDistanceMatrix.from_condensed(condensed, probe_triples=0)
    assert dm.metric_probe == "not_run"
    oecluster.butina(dm, 0.5)


def test_the_probe_is_deterministic():
    condensed = np.full(6, 0.1)
    condensed[1] = 10.0
    first = SymmetricDistanceMatrix.from_condensed(condensed,
                                                   probe_triples=500)
    second = SymmetricDistanceMatrix.from_condensed(condensed,
                                                    probe_triples=500)
    assert first.probe_violations == second.probe_violations
    assert first.probe_sampled == second.probe_sampled


def test_check_false_skips_validation_and_the_probe():
    """The user-asserted-safe escape hatch.

    A negative entry and a planted violation both pass, and the probe does not
    run. The matrix is then indistinguishable at the gate from a 4.x matrix
    with no recorded provenance.
    """
    condensed = np.full(6, 0.1)
    condensed[1] = 10.0
    condensed[2] = -1.0
    dm = SymmetricDistanceMatrix.from_condensed(condensed, check=False)
    assert dm.metric_probe == "not_run"
    assert dm.metric_capabilities == {'zero_self': "unknown",
                                      'triangle': "unknown"}
    assert dm.data_integrity == "complete"
    oecluster.butina(dm, 0.5)


def test_check_false_still_rejects_a_structurally_impossible_input():
    """Skipping the value checks does not skip the shape arithmetic.

    A length that maps to no item count cannot be stored at all, so it raises
    whatever ``check`` says.
    """
    with pytest.raises(ValueError, match="not a valid condensed length"):
        SymmetricDistanceMatrix.from_condensed(np.zeros(4), check=False)


def test_params_are_recorded_and_default_to_empty():
    condensed, _ = _metric_condensed(n=3)
    assert SymmetricDistanceMatrix.from_condensed(condensed).params == {}
    dm = SymmetricDistanceMatrix.from_condensed(
        condensed, params={"source": "scipy", "metric": "cosine"})
    assert dm.params == {"source": "scipy", "metric": "cosine"}


def test_params_are_copied_from_the_callers_dict():
    """The docstring promises a copy, so the promise is pinned here."""
    condensed, _ = _metric_condensed(n=3)
    supplied = {"source": "scipy"}
    dm = SymmetricDistanceMatrix.from_condensed(condensed, params=supplied)
    supplied["source"] = "somewhere else"
    assert dm.params == {"source": "scipy"}


def test_the_probe_tolerance_keeps_a_collinear_metric_clean():
    """Over-refusal is the failure mode the probe's tolerance exists to stop.

    Collinear points make the triangle inequality an equality, which is where
    floating-point rounding pushes a genuine metric a few ulps over the line.
    Measured on this 200-point set: 1256 violations in 97273 sampled triples
    with the ``1e-9 * scale`` term deleted, and 0 with it in place.
    """
    points = np.arange(200, dtype=np.float64)[:, None] * np.array([0.1, 0.0])
    diff = points[:, None, :] - points[None, :, :]
    square = np.sqrt((diff ** 2).sum(axis=-1))
    condensed = square[np.triu_indices(200, k=1)]
    assert _gate.probe_triangle(condensed, 200) == {
        'metric_probe': "no_violations_found",
        'probe_violations': 0,
        'probe_sampled': 97273,
    }


def test_the_probe_counts_distinct_triples_not_draws():
    """A draw count would overstate the evidence the refusal cites.

    A 4-item matrix admits exactly 12 inequalities of the form
    ``d(i, k) <= d(i, j) + d(j, k)``: ``C(4, 2)`` unordered end pairs times the
    2 remaining choices of ``j``. The sampler draws with replacement, so a
    count of surviving draws puts ``probe_sampled`` in the tens of thousands
    for a matrix that contains twelve tests.
    """
    condensed = np.full(6, 0.1)
    condensed[1] = 10.0
    result = _gate.probe_triangle(condensed, 4)
    assert result['probe_sampled'] == 12
    # d(0, 2) = 10.0 against 0.1 everywhere else, so the two inequalities that
    # route from 0 to 2 through 1 and through 3 are the violations that exist.
    assert result['probe_violations'] == 2


def test_deduplication_still_leaves_a_large_sample_on_a_large_set():
    """Deduplication must not quietly shrink the probe into uselessness.

    ``probe_sampled`` depends only on ``n``, the draw count and the seed, so
    this figure is a property of the sampler rather than of the distances.
    """
    condensed, _ = _metric_condensed(n=60)
    result = _gate.probe_triangle(condensed, 60)
    assert result['probe_sampled'] == 62129


def test_the_probe_is_skipped_below_three_items():
    dm = SymmetricDistanceMatrix.from_condensed(np.array([0.5]))
    assert dm.metric_probe == "not_run"
    assert dm.num_samples == 2


@pytest.mark.parametrize("metric", ["cosine", "correlation", "sqeuclidean"])
def test_the_probe_catches_real_non_metric_distances(metric):
    """scipy's non-metric distances are caught at the default sample rate.

    A planted violation proves the arithmetic; this proves the sampler finds
    violations at the density real non-metrics actually produce. Measured on
    this 60-point set: cosine 7708 violations in 62129 sampled triples.
    """
    scipy_distance = pytest.importorskip("scipy.spatial.distance")
    rng = np.random.default_rng(7)
    points = rng.normal(size=(60, 4))
    condensed = scipy_distance.pdist(points, metric=metric)
    dm = SymmetricDistanceMatrix.from_condensed(condensed)
    assert dm.metric_probe == "violations_found"
    assert dm.probe_violations > 0
    with pytest.raises(ValueError, match="sampled triples"):
        oecluster.butina(dm, 0.5)
    oecluster.butina(dm, 0.5, allow_nonmetric=True)


def test_the_probe_does_not_flag_scipy_euclidean():
    """The counterpart to the parametrized case: no false positives."""
    scipy_distance = pytest.importorskip("scipy.spatial.distance")
    rng = np.random.default_rng(7)
    points = rng.normal(size=(60, 4))
    condensed = scipy_distance.pdist(points, metric="euclidean")
    dm = SymmetricDistanceMatrix.from_condensed(condensed)
    assert dm.metric_probe == "no_violations_found"
    assert dm.probe_violations == 0


def test_pdist_round_trips_through_from_condensed():
    """A computed matrix survives export and re-ingress bit-identically."""
    from openeye import oechem

    smiles = ["c1ccccc1", "c1ccc(O)cc1", "CCCCCCCC", "c1ccncc1", "CCO"]
    mols = []
    for smi in smiles:
        mol = oechem.OEGraphMol()
        oechem.OESmilesToMol(mol, smi)
        mols.append(mol)

    computed = oecluster.pdist(mols, "fingerprint")
    reloaded = SymmetricDistanceMatrix.from_condensed(
        np.asarray(computed), labels=computed.labels)
    np.testing.assert_array_equal(np.asarray(reloaded), np.asarray(computed))
    assert reloaded.labels == computed.labels
    assert reloaded.data_integrity == "complete"
    # The numbers survive; the provenance does not. Tanimoto's proven zero
    # self-distance is not recoverable from the array it produced.
    assert computed.metric_capabilities['zero_self'] is True
    assert reloaded.metric_capabilities['zero_self'] == "unknown"


def test_condensed_lookup_matches_the_square_form():
    condensed, square = _metric_condensed(n=7)
    i = np.array([0, 3, 6, 2])
    j = np.array([5, 1, 0, 4])
    np.testing.assert_allclose(
        _gate.condensed_lookup(condensed, 7, i, j), square[i, j])


def test_from_file_reproduces_the_matrix(tmp_path):
    condensed, _ = _metric_condensed(n=40)
    path = tmp_path / "big.npz"
    SymmetricDistanceMatrix.from_condensed(
        condensed, labels=[f"m{i}" for i in range(40)]).to_file(str(path))
    loaded = oecluster.load_distance_matrix(str(path))
    assert loaded.num_samples == 40
    np.testing.assert_array_equal(loaded.condensed, condensed)
    assert loaded.labels[-1] == "m39"
