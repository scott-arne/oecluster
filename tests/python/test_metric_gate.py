import json

import numpy as np
import oecluster
import pytest
from oecluster import _gate
from openeye import oechem

SMILES = ["CCO", "CCC", "CCCC", "c1ccccc1", "CCN", "CCOC"]


def _mols(smiles_list=None):
    mols = []
    for idx, smi in enumerate(smiles_list or SMILES):
        mol = oechem.OEGraphMol()
        oechem.OESmilesToMol(mol, smi)
        mol.SetTitle(f"mol{idx}")
        mols.append(mol)
    return mols


def _conformer_mols(smiles_list=None):
    """Molecules embedded in 3D, so ROCS has something real to overlay.

    ``ROCSComparison`` refuses a molecule whose recomputed OEChem dimension
    attribute -- an axis count, not a geometric rank -- is below three, and a
    molecule straight from a SMILES parse carries no coordinates at all, so
    these fixtures have to come from Omega. The typemap needs ``OEMol``: an
    ``OEGraphMol`` fails earlier and differently, with a SWIG ``TypeError``
    rather than the dimension refusal.
    """
    pytest.importorskip("openeye.oeomega")
    from openeye import oeomega

    omega = oeomega.OEOmega()
    omega.SetMaxConfs(1)
    omega.SetStrictStereo(False)

    mols = []
    for idx, smi in enumerate(smiles_list or SMILES):
        mol = oechem.OEMol()
        oechem.OESmilesToMol(mol, smi)
        mol.SetTitle(f"mol{idx}")
        assert omega(mol)
        mols.append(mol)
    return mols


def test_default_fingerprint_path_is_stamped_metric():
    dist = oecluster.pdist(_mols(), "fingerprint")
    assert dist.metric_capabilities == {'zero_self': True, 'triangle': True}
    assert dist.data_integrity == "complete"
    assert dist.metric_probe == "not_run"


def test_default_fingerprint_path_clusters_without_an_override():
    dist = oecluster.pdist(_mols(), "fingerprint")
    result = oecluster.butina(dist, 0.5)
    assert len(result.labels) == 6


def test_similarity_matrix_is_refused_by_butina():
    dist = oecluster.pdist(_mols(), "fingerprint", similarity=True)
    # ``is_distance`` is the fact this refusal fires on; ``zero_self`` is
    # asserted alongside it only to record that both are False here.
    assert dist.is_distance is False
    assert dist.metric_capabilities['zero_self'] is False
    with pytest.raises(ValueError, match="similarity=False"):
        oecluster.butina(dist, 0.5)


def test_the_similarity_refusal_cannot_be_overridden():
    dist = oecluster.pdist(_mols(), "fingerprint", similarity=True)
    with pytest.raises(ValueError, match="cannot be overridden"):
        oecluster.butina(dist, 0.5, allow_nonmetric=True)


def test_dice_is_refused_for_violating_the_triangle_inequality():
    dist = oecluster.pdist(_mols(), "fingerprint", metric="dice")
    assert dist.metric_capabilities == {'zero_self': True, 'triangle': False}
    with pytest.raises(ValueError, match="triangle inequality"):
        oecluster.butina(dist, 0.5)


def test_dice_is_allowed_with_allow_nonmetric():
    dist = oecluster.pdist(_mols(), "fingerprint", metric="dice")
    result = oecluster.butina(dist, 0.5, allow_nonmetric=True)
    assert len(result.labels) == 6


@pytest.mark.parametrize("call", [
    lambda dm, **kw: oecluster.butina(dm, 0.5, **kw),
    lambda dm, **kw: oecluster.dbscan(dm, 0.5, **kw),
    lambda dm, **kw: oecluster.hdbscan(dm, min_cluster_size=2, **kw),
    lambda dm, **kw: oecluster.agglomerative(dm, n_clusters=2, **kw),
])
def test_every_clustering_entry_point_refuses_and_overrides(call):
    dist = oecluster.pdist(_mols(), "fingerprint", metric="dice")
    with pytest.raises(ValueError, match="triangle inequality"):
        call(dist)
    call(dist, allow_nonmetric=True)


def test_cluster_report_refuses_and_overrides():
    dist = oecluster.pdist(_mols(), "fingerprint", metric="dice")
    result = oecluster.butina(dist, 0.5, allow_nonmetric=True)
    with pytest.raises(ValueError, match="triangle inequality"):
        oecluster.cluster_report(result, dist)
    oecluster.cluster_report(result, dist, allow_nonmetric=True)


def test_a_prebuilt_comparison_object_is_still_stamped():
    opts = oecluster.FingerprintOptions()
    opts.metric = "dice"
    comparison = oecluster.oecluster.FingerprintComparison(_mols(), opts)
    dist = oecluster.pdist(_mols(), comparison)
    assert dist.metric_capabilities['triangle'] is False
    with pytest.raises(ValueError, match="triangle inequality"):
        oecluster.butina(dist, 0.5)


def test_rocs_similarity_is_refused_though_no_metric_is_involved():
    """The gate reads facts off the comparison, not off a metric name.

    ROCS never consults the metric table, so a metric-name-only gate would let
    this similarity matrix through. This is the hole ``GateFacts`` closes.

    ``is_distance`` is what makes this a similarity, and it is the fact the
    gate refuses on. ``zero_self`` reports the diagonal, which for these six
    molecules is 1.0 -- so here the two facts happen to agree. The test below
    is the one where they come apart.
    """
    dist = oecluster.pdist(_conformer_mols(), "rocs", similarity=True)
    assert dist.is_distance is False
    assert dist.metric_capabilities['zero_self'] is False
    with pytest.raises(ValueError, match="similarity=False"):
        oecluster.butina(dist, 0.5)
    with pytest.raises(ValueError, match="cannot be overridden"):
        oecluster.butina(dist, 0.5, allow_nonmetric=True)


def test_a_similarity_that_self_scores_zero_is_still_refused():
    """Why ``is_distance`` and ``zero_self`` have to be separate facts.

    ``GateFacts.h`` keeps them apart on the grounds that "a comparison can be
    a similarity whose self-value happens to be zero." Methane is that case:
    it carries no colour features, so its colour self-similarity is 0.0, and a
    set containing only methane stamps ``zero_self`` True on a *similarity*.
    A gate reading ``zero_self`` would admit this matrix and cluster a
    similarity as if it were a distance. Reading ``is_distance`` refuses it.

    This is the test that makes the two-fact split non-vacuous. If it is ever
    deleted, nothing else in the suite distinguishes the two facts.
    """
    dist = oecluster.pdist(_conformer_mols(["C", "C"]), "rocs",
                           score_type="color", similarity=True)
    assert dist.metric_capabilities['zero_self'] is True
    assert dist.is_distance is False
    with pytest.raises(ValueError, match="similarity=False"):
        oecluster.butina(dist, 0.5)
    with pytest.raises(ValueError, match="cannot be overridden"):
        oecluster.butina(dist, 0.5, allow_nonmetric=True)


def test_a_distance_whose_diagonal_does_not_vanish_is_refused():
    """The complement: a genuine distance that fails the other tier-1 axiom.

    ``combo_norm`` averages a shape Tanimoto with a colour Tanimoto. Methane
    carries no colour features, so its colour self-similarity is 0.0 and only
    the shape half of its self-score saturates: the diagonal sits at 0.5
    rather than 0.0. The matrix is oriented correctly -- larger still means
    further apart, so ``is_distance`` is True -- and it still fails the
    zero-self axiom, which is its own non-overridable refusal with its own
    message.

    Together with ``test_a_similarity_that_self_scores_zero_is_still_refused``
    this closes the square. That test is ``is_distance`` False with
    ``zero_self`` True; this one is ``is_distance`` True with ``zero_self``
    False. Neither fact implies the other in either direction, and each has a
    refusal the other cannot stand in for -- which is the whole justification
    for ``GateFacts`` keeping them apart. Delete either one and a tier-1
    refusal goes untested.
    """
    dist = oecluster.pdist(_conformer_mols(["C", "C"]), "rocs")
    assert dist.is_distance is True
    assert dist.metric_capabilities['zero_self'] is False
    with pytest.raises(ValueError, match="zero self-distance"):
        oecluster.butina(dist, 0.5)
    with pytest.raises(ValueError, match="cannot be overridden"):
        oecluster.butina(dist, 0.5, allow_nonmetric=True)


def test_rocs_shape_distance_is_admitted_with_no_override():
    """The counterpart: a ROCS configuration that is a true distance.

    On this six-molecule fixture all four score types are true distances now
    that colour atoms are prepared -- see the test below. ``shape`` is tested
    here because it is the one whose self-distance is zero for *every*
    molecule: on a featureless fragment such as methane the other three stamp
    ``zero_self`` No, because their colour term contributes nothing to the
    self-score. Only ``shape`` never depended on the colour repair.
    """
    dist = oecluster.pdist(_conformer_mols(), "rocs", score_type="shape")
    assert dist.is_distance is True
    assert dist.metric_capabilities == {'zero_self': True,
                                        'triangle': "unknown"}
    assert dist.data_integrity == "complete"
    result = oecluster.butina(dist, 0.5)
    assert len(result.labels) == 6


def test_rocs_combo_norm_is_admitted_now_that_colour_is_prepared():
    """The default ROCS distance vanishes on the diagonal for this fixture.

    Not universally: ``combo_norm`` averages a shape term with a colour term,
    so its self-distance is zero only for molecules that carry colour
    features. Every molecule here does. Methane does not, and
    ``test_a_distance_whose_diagonal_does_not_vanish_is_refused`` covers that
    side.

    This test used to be a refusal. Nothing in the repository prepared colour
    atoms, so ``GetColorTanimoto()`` returned 0.0 for every pair including
    self-pairs, and ``combo_norm``'s self-distance was 0.5. The stamp reported
    what the scorer did, and the earlier version of this test recorded that
    honestly while saying it should flip to an admission once the colour term
    was repaired. It has been repaired, so this is that flip.

    Keep it. It is the regression test for the colour preparation: if
    ``InitOverlay`` stops assigning colour atoms, the diagonal returns to 0.5
    and this admission fails.
    """
    dist = oecluster.pdist(_conformer_mols(), "rocs")
    assert dist.is_distance is True
    assert dist.metric_capabilities['zero_self'] is True
    result = oecluster.butina(dist, 0.5)
    assert len(result.labels) == 6


def test_a_prebuilt_rocs_object_is_stamped_too():
    """The prebuilt-object branch skips the builder; the facts still apply."""
    mols = _conformer_mols()
    comparison = oecluster.ROCSComparison(mols, similarity=True)
    dist = oecluster.pdist(mols, comparison)
    assert dist.is_distance is False
    with pytest.raises(ValueError, match="similarity=False"):
        oecluster.butina(dist, 0.5)


def test_cdist_stamps_facts_too():
    cross = oecluster.cdist(_mols()[:2], _mols()[2:], "fingerprint",
                            metric="dice")
    assert cross.metric_capabilities['triangle'] is False


def test_facts_round_trip_through_to_file(tmp_path):
    path = tmp_path / "dice.npz"
    oecluster.pdist(_mols(), "fingerprint", metric="dice").to_file(str(path))
    loaded = oecluster.load_distance_matrix(str(path))
    assert loaded.metric_capabilities == {'zero_self': True, 'triangle': False}
    assert loaded.data_integrity == "complete"


def test_cross_facts_round_trip_through_to_file(tmp_path):
    path = tmp_path / "cross.npz"
    oecluster.cdist(_mols()[:2], _mols()[2:], "fingerprint",
                    metric="dice").to_file(str(path))
    loaded = oecluster.load_distance_matrix(str(path))
    assert loaded.metric_capabilities == {'zero_self': True, 'triangle': False}


def test_a_file_without_facts_loads_as_unknown(tmp_path):
    """A matrix written before 5.0.0 must not be credited with a claim."""
    path = tmp_path / "legacy.npz"
    np.savez_compressed(
        str(path),
        condensed=np.zeros(3, dtype=np.float64),
        comparison_name=np.array("fingerprint"),
        params_json=np.array(json.dumps({})),
        labels=np.array(["a", "b", "c"]),
        num_samples=np.array(3),
    )
    loaded = oecluster.load_distance_matrix(str(path))
    assert loaded.metric_capabilities == {
        'zero_self': "unknown", 'triangle': "unknown"}
    assert loaded.data_integrity == "unknown"


def test_unknown_facts_never_refuse(tmp_path):
    path = tmp_path / "legacy.npz"
    np.savez_compressed(
        str(path),
        condensed=np.array([0.1, 0.2, 0.3], dtype=np.float64),
        comparison_name=np.array("fingerprint"),
        params_json=np.array(json.dumps({})),
        labels=np.array(["a", "b", "c"]),
        num_samples=np.array(3),
    )
    loaded = oecluster.load_distance_matrix(str(path))
    assert len(oecluster.butina(loaded, 0.5).labels) == 3


def test_require_metric_reports_probe_violations():
    dist = oecluster.pdist(_mols(), "fingerprint")
    dist._facts.update(metric_probe="violations_found",
                       probe_violations=668, probe_sampled=99926)
    with pytest.raises(ValueError) as excinfo:
        _gate.require_metric(dist, "butina")
    message = str(excinfo.value)
    assert "668 of 99926 sampled triples" in message
    assert "allow_nonmetric=True" in message


def test_require_metric_accepts_probe_violations_with_the_override():
    dist = oecluster.pdist(_mols(), "fingerprint")
    dist._facts.update(metric_probe="violations_found",
                       probe_violations=1, probe_sampled=10)
    _gate.require_metric(dist, "butina", allow_nonmetric=True)


def test_require_metric_refuses_nan_present_even_with_the_override():
    dist = oecluster.pdist(_mols(), "fingerprint")
    dist._facts['data_integrity'] = "nan_present"
    with pytest.raises(ValueError, match="cannot be overridden"):
        _gate.require_metric(dist, "butina", allow_nonmetric=True)


def test_require_metric_refuses_subset_scored_but_allows_the_override():
    dist = oecluster.pdist(_mols(), "fingerprint")
    dist._facts['data_integrity'] = "subset_scored"
    with pytest.raises(ValueError, match="subset"):
        _gate.require_metric(dist, "butina")
    _gate.require_metric(dist, "butina", allow_nonmetric=True)


def test_butina_coercion_error_before_gate():
    """Caller argument coercion errors come before the gate's advisory refusal.

    Discriminator: passing a list to int() raises TypeError, while the gate
    raises ValueError. If the gate ran first, the caller would see ValueError
    with the "allow_nonmetric=True" advisory -- a remedy that cannot rescue
    the invalid type.
    """
    dist = oecluster.pdist(_mols(), "fingerprint", metric="dice")
    assert dist.metric_capabilities['triangle'] is False
    with pytest.raises(TypeError):
        oecluster.butina(dist, 0.5, num_threads=["not_an_int"])


def test_agglomerative_coercion_error_before_gate():
    """Caller argument coercion errors come before the gate's advisory refusal.

    This covers the n_clusters coercion, which is only guarded when
    distance_threshold is None. When distance_threshold is provided,
    n_clusters validation is skipped above the gate, so the coercion must
    be hoisted to avoid the gate pre-empting the type error.

    Discriminator: passing a list to int() raises TypeError, while the gate
    raises ValueError.
    """
    dist = oecluster.pdist(_mols(), "fingerprint", metric="dice")
    assert dist.metric_capabilities['triangle'] is False
    with pytest.raises(TypeError):
        oecluster.agglomerative(dist, n_clusters=["not_an_int"],
                                distance_threshold=0.5)


def _nonmetric():
    """A dense six-item matrix whose facts arm the gate's tier-2 advisory."""
    dist = oecluster.pdist(_mols(), "fingerprint", metric="dice")
    assert dist.metric_capabilities['triangle'] is False
    return dist


def _nonmetric_sparse():
    """The same measure stored sparsely, so the gate is armed here too."""
    sparse = oecluster.pdist(_mols(), "fingerprint", metric="dice",
                             cutoff=0.2)
    assert sparse.metric_capabilities['triangle'] is False
    assert isinstance(sparse.storage, oecluster.SparseStorage)
    return sparse


def test_hdbscan_item_count_bound_before_gate():
    """An out-of-range min_samples outranks the gate's advisory.

    ``allow_nonmetric=True`` cannot make min_samples fit the item count, so
    naming it as the remedy sends the caller down a dead end. Matching on
    "at most the item count" is what distinguishes the authoritative message
    from the advisory, which shares the ValueError type.
    """
    dist = _nonmetric()
    with pytest.raises(ValueError, match="min_samples must be at most"):
        oecluster.hdbscan(dist, min_samples=99)


def test_hdbscan_item_count_bound_reads_the_effective_min_samples():
    """The bound applies to the value HDBSCAN actually uses.

    ``HDBSCAN.cpp`` substitutes min_cluster_size when min_samples is unset,
    then bounds that substituted value. A mirror that inspected only the raw
    keyword would miss this call, which the native code refuses.
    """
    dist = _nonmetric()
    with pytest.raises(ValueError, match="min_samples must be at most"):
        oecluster.hdbscan(dist, min_cluster_size=7)


def test_agglomerative_item_count_bound_before_gate():
    dist = _nonmetric()
    with pytest.raises(ValueError, match="n_clusters must be at most"):
        oecluster.agglomerative(dist, n_clusters=99)


def test_agglomerative_item_count_bound_only_without_a_threshold():
    """Not an over-refusal: the native bound is conditional.

    ``Agglomerative.cpp`` bounds n_clusters only when distance_threshold is
    negative, which is how the wrapper encodes "no threshold". With a
    threshold supplied, n_clusters is ignored, and an out-of-range value must
    still cluster.
    """
    dist = _nonmetric()
    result = oecluster.agglomerative(dist, n_clusters=99,
                                     distance_threshold=0.5,
                                     allow_nonmetric=True)
    assert len(result.labels) == 6


@pytest.mark.parametrize("call", [
    lambda dm, **kw: oecluster.butina(dm, 0.5, **kw),
    lambda dm, **kw: oecluster.dbscan(dm, 0.5, **kw),
    lambda dm, **kw: oecluster.hdbscan(dm, min_cluster_size=2, **kw),
    lambda dm, **kw: oecluster.agglomerative(dm, **kw),
])
def test_negative_num_threads_before_gate(call):
    """A negative thread count outranks the gate's advisory.

    The native options field is a size_t, so the assignment raises
    OverflowError once the gate lets the call through -- another failure
    ``allow_nonmetric=True`` cannot rescue.
    """
    dist = _nonmetric()
    with pytest.raises(ValueError, match="num_threads must be non-negative"):
        call(dist, num_threads=-1)


def test_negative_num_threads_before_gate_in_cluster_report():
    dist = _nonmetric()
    result = oecluster.butina(dist, 0.5, allow_nonmetric=True)
    with pytest.raises(ValueError, match="num_threads must be non-negative"):
        oecluster.cluster_report(result, dist, num_threads=-1)


@pytest.mark.parametrize("call", [
    lambda dm, **kw: oecluster.butina(dm, 0.5, **kw),
    lambda dm, **kw: oecluster.dbscan(dm, 0.5, **kw),
    lambda dm, **kw: oecluster.hdbscan(dm, min_cluster_size=2, **kw),
    lambda dm, **kw: oecluster.agglomerative(dm, **kw),
])
def test_negative_chunk_size_before_gate(call):
    dist = _nonmetric()
    with pytest.raises(ValueError, match="chunk_size must be non-negative"):
        call(dist, chunk_size=-1)


def test_negative_n_clusters_before_gate_with_a_threshold():
    """The n_clusters size_t assignment is unguarded when a threshold is set.

    With distance_threshold supplied the wrapper skips its "at least one"
    check, so a negative n_clusters reached the size_t setter and raised
    OverflowError below the gate.
    """
    dist = _nonmetric()
    with pytest.raises(ValueError, match="n_clusters must be non-negative"):
        oecluster.agglomerative(dist, n_clusters=-1, distance_threshold=0.5)


@pytest.mark.parametrize("call", [
    lambda dm, **kw: oecluster.hdbscan(dm, min_cluster_size=2, **kw),
    lambda dm, **kw: oecluster.agglomerative(dm, **kw),
])
def test_sparse_storage_rejected_before_gate(call):
    """The algorithms that need every pair say so before the advisory.

    These three are the callers of ``validate_complete_distance_storage``;
    butina and dbscan are not, and the test below holds them to that.
    """
    sparse = _nonmetric_sparse()
    with pytest.raises(ValueError, match="SparseStorage is not supported"):
        call(sparse)


def test_sparse_storage_rejected_before_gate_in_cluster_report():
    dist = _nonmetric()
    sparse = _nonmetric_sparse()
    result = oecluster.butina(dist, 0.5, allow_nonmetric=True)
    with pytest.raises(ValueError, match="SparseStorage is not supported"):
        oecluster.cluster_report(result, sparse)


def test_butina_and_dbscan_still_accept_sparse_storage():
    """Not an over-refusal: both build a threshold graph from sparse entries.

    A mirror stricter than the native rule breaks working code silently,
    because no native error contradicts it. This is the guard against that.
    """
    sparse = _nonmetric_sparse()
    assert len(oecluster.butina(sparse, 0.1, allow_nonmetric=True).labels) == 6
    assert len(oecluster.dbscan(sparse, 0.1, allow_nonmetric=True).labels) == 6


def test_zero_is_still_a_valid_size_t_argument():
    """Not an over-refusal: zero means "choose for me", and stays legal."""
    dist = _nonmetric()
    result = oecluster.butina(dist, 0.5, num_threads=0, chunk_size=0,
                              allow_nonmetric=True)
    assert len(result.labels) == 6
    assert len(oecluster.hdbscan(dist, min_cluster_size=2, num_threads=0,
                                 chunk_size=0,
                                 allow_nonmetric=True).labels) == 6
    report = oecluster.cluster_report(result, dist, num_threads=0,
                                      allow_nonmetric=True)
    assert report.num_samples == 6


@pytest.mark.parametrize("call", [
    lambda dm, **kw: oecluster.butina(dm, 0.5, **kw),
    lambda dm, **kw: oecluster.dbscan(dm, 0.5, **kw),
    lambda dm, **kw: oecluster.hdbscan(dm, min_cluster_size=2, **kw),
    lambda dm, **kw: oecluster.agglomerative(dm, **kw),
])
def test_a_non_bool_allow_nonmetric_is_refused(call):
    """A truthy string silently disabled the gate; now it is a TypeError.

    ``allow_nonmetric="False"`` reads as an override, so the safety check the
    caller believed they were keeping was switched off without a word.
    """
    dist = _nonmetric()
    with pytest.raises(TypeError, match="allow_nonmetric must be True or False"):
        call(dist, allow_nonmetric="False")


def test_a_non_bool_allow_nonmetric_is_refused_in_cluster_report():
    dist = _nonmetric()
    result = oecluster.butina(dist, 0.5, allow_nonmetric=True)
    with pytest.raises(TypeError, match="allow_nonmetric must be True or False"):
        oecluster.cluster_report(result, dist, allow_nonmetric="False")


def test_numpy_bool_true_overrides_the_gate():
    """numpy.bool_ is accepted and coerces faithfully."""
    dist = _nonmetric()
    result = oecluster.butina(dist, 0.5, allow_nonmetric=np.True_)
    assert len(result.labels) == 6


def test_numpy_bool_false_leaves_the_gate_enforced():
    """numpy.bool_ is accepted, and False still means the gate is on."""
    dist = _nonmetric()
    with pytest.raises(ValueError, match="triangle inequality"):
        oecluster.butina(dist, 0.5, allow_nonmetric=np.False_)


def test_a_string_false_still_raises_type_error():
    """Regression guard: the numpy.bool_ widening must not admit strings."""
    dist = _nonmetric()
    with pytest.raises(TypeError, match="allow_nonmetric must be True or False"):
        oecluster.butina(dist, 0.5, allow_nonmetric="False")


def test_agglomerative_chunk_size_zero_before_gate():
    """Zero chunk_size is rejected before the advisory for agglomerative.

    The override cannot make chunk_size valid, so the gate must not pre-empt
    this message. Match on "at least one" to distinguish from the advisory.
    """
    dist = _nonmetric()
    assert dist.metric_capabilities['triangle'] is False
    with pytest.raises(ValueError, match="chunk_size must be at least one"):
        oecluster.agglomerative(dist, chunk_size=0)


def test_agglomerative_nan_distance_threshold_before_gate():
    """NaN distance_threshold is rejected before the advisory.

    The override cannot make NaN valid, so the gate must not pre-empt this
    message. Match on "must not be NaN" to distinguish from the advisory.
    """
    dist = _nonmetric()
    assert dist.metric_capabilities['triangle'] is False
    with pytest.raises(ValueError, match="distance_threshold must not be NaN"):
        oecluster.agglomerative(dist, distance_threshold=float('nan'))


def test_butina_dbscan_hdbscan_still_accept_chunk_size_zero():
    """Zero chunk_size is legal for butina, dbscan, and hdbscan.

    Agglomerative rejects zero; the other three must not. A mirror stricter
    than the native rule breaks working code silently.
    """
    dist = oecluster.pdist(_mols(), "fingerprint")
    assert len(oecluster.butina(dist, 0.5, chunk_size=0).labels) == 6
    assert len(oecluster.dbscan(dist, 0.5, chunk_size=0).labels) == 6
    assert len(oecluster.hdbscan(dist, min_cluster_size=2, chunk_size=0).labels) == 6
