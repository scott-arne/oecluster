import inspect
import json
import math

import numpy as np
import oecluster
import pytest
from oecluster import _gate
from oecluster.oecluster import butina_cluster as _butina_cluster
from openeye import oechem

SMILES = ["CCO", "CCC", "CCCC", "c1ccccc1", "CCN", "CCOC"]


# The three refusals require_comparable inherits from require_metric and never
# waives. Each acceptance test replays all three against its waived fact, so a
# tier-1 check that gets reordered behind a waiver is caught rather than passing
# on whichever fault that test happened to pick.
_TIER1_FAULTS = [
    ("is_distance", lambda d: d._facts.__setitem__('is_distance', False),
     "requires distances"),
    ("zero_self", lambda d: d._facts.__setitem__('zero_self', False),
     "zero self-distance"),
    ("non_finite", lambda d: d.condensed.__setitem__(0, math.nan),
     "non-finite entries"),
]


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
    carries no colour features, so its colour self-similarity is 0.0, and on
    the conformer ``_conformer_mols`` builds its shape half does saturate: the
    diagonal measures exactly 0.5. That 0.5 belongs to this fixture rather than
    to methane -- embedded from explicit hydrogens instead, methane's shape half
    stops saturating too and the diagonal moves to 5.05e-01. A colourless
    molecule is therefore one way to reach this refusal and not the only one;
    ``test_rocs_shape_distance_is_admitted_with_no_override`` names the other.
    The matrix is oriented correctly -- larger still means further apart, so
    ``is_distance`` is True -- and it still fails the zero-self axiom, which is
    its own non-overridable refusal with its own message.

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
    here because it is the one that never depended on the colour repair: on a
    featureless fragment such as methane the other three stamp ``zero_self``
    No, because their colour term contributes nothing to the self-score.

    That is a reason to choose ``shape``, not a guarantee about it. What the
    assertion below establishes is that ``shape`` vanishes on *these six*
    molecules, and that does not generalise: ``_conformer_mols(["CO", "CO"])``
    -- methanol, which does carry colour features, its colour self-distance
    measuring 0.0 -- stamps ``zero_self`` False under ``score_type="shape"``,
    with a self-distance of 1.04e-02. A small compact molecule whose
    self-overlay does not quite saturate fails the diagonal whatever its
    colour term does.
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
    and needs both to vanish. Carrying colour features is therefore necessary
    but not sufficient -- methanol carries them and still self-scores
    5.18e-03, because its shape term does not saturate. Every molecule here
    clears both terms. Methane clears the shape term and not the colour one,
    and ``test_a_distance_whose_diagonal_does_not_vanish_is_refused`` covers
    that side.

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


def test_a_nan_written_through_condensed_is_refused():
    """The gate measures the data, not just the stamp it was handed.

    ``.condensed`` is the live zero-copy view over dense storage, so a caller
    can write a NaN into a matrix long after the comparison stamped it
    ``complete``. The stamp assertion comes first on purpose: it records that
    the stamp still says ``complete``, so the refusal can only have come from
    the measurement.
    """
    dist = oecluster.pdist(_mols(), "fingerprint")
    dist.condensed[0] = float('nan')
    assert dist.data_integrity == "complete"
    with pytest.raises(ValueError, match="non-finite entries"):
        oecluster.butina(dist, threshold=0.4)


def test_a_nan_written_through_storage_set_is_refused():
    """The same defect by the other write route.

    ``storage.Set`` reaches the buffer without going through the ``.condensed``
    property, and like the property it leaves the facts untouched, so a gate
    reading only the stamp admits this route too. What makes the scan catch it
    is that the two routes are not independent: for dense storage
    ``.condensed`` is a zero-copy view over the very array ``DenseStorage::Set``
    writes into, so scanning the view sees the ``Set``.
    """
    dist = oecluster.pdist(_mols(), "fingerprint")
    dist.storage.Set(0, 1, float('nan'))
    assert dist.data_integrity == "complete"
    with pytest.raises(ValueError, match="non-finite entries"):
        oecluster.dbscan(dist, eps=0.4, min_samples=2)


def test_the_nan_refusal_names_only_a_remedy_the_caller_can_take():
    """``missing`` is a descriptor-only keyword.

    Every route the scan above exists for -- ``.condensed``, ``storage.Set``,
    a trusted file -- can deliver a NaN into any comparison's matrix, and
    offering ``missing='complete_case'`` to the ones that never took the
    keyword costs a second refusal. Removing the offending items needs no
    option, so that is the remedy both branches keep.
    """
    mols = _mols()

    descriptor = oecluster.pdist(mols, "descriptor")
    descriptor.condensed[0] = float('nan')
    with pytest.raises(ValueError, match=r"missing='complete_case'"):
        oecluster.butina(descriptor, threshold=0.4)

    fingerprint = oecluster.pdist(mols, "fingerprint")
    fingerprint.condensed[0] = float('nan')
    with pytest.raises(ValueError, match=r"Remove the offending items\.") as excinfo:
        oecluster.butina(fingerprint, threshold=0.4)
    assert "missing=" not in str(excinfo.value)

    # The named remedy has to be a keyword its own path accepts, and the
    # withheld one has to be a keyword the other path refuses. Neither half
    # holds by inspection of the gate alone.
    oecluster.pdist(mols, "descriptor", missing="complete_case")
    with pytest.raises(TypeError, match=r"Unknown kwargs for fingerprint"):
        oecluster.pdist(mols, "fingerprint", missing="complete_case")


def test_a_nan_bearing_file_is_refused_though_its_stamp_says_complete(tmp_path):
    """The file route needs no deliberate mutation of the loaded matrix.

    ``from_file`` trusts ``facts_json`` outright and never checks it against
    the ``condensed`` array it just loaded, so an untrusted or stale file
    arrives with NaN data under a ``complete`` stamp. Deliberately left as a
    known residual: ``loaded.data_integrity`` still reports ``complete``. The
    gate is what stops it clustering.
    """
    path = tmp_path / "tampered.npz"
    dist = oecluster.pdist(_mols(), "fingerprint")
    dist.condensed[0] = float('nan')
    dist.to_file(str(path))

    loaded = oecluster.load_distance_matrix(str(path))
    assert loaded.data_integrity == "complete"
    assert np.isnan(loaded.condensed[0])
    with pytest.raises(ValueError, match="non-finite entries"):
        oecluster.butina(loaded, threshold=0.4)


@pytest.mark.parametrize("call", [
    lambda dm, **kw: oecluster.butina(dm, 0.4, **kw),
    lambda dm, **kw: oecluster.dbscan(dm, 0.4, **kw),
    lambda dm, **kw: oecluster.hdbscan(dm, min_cluster_size=2, **kw),
    lambda dm, **kw: oecluster.agglomerative(dm, n_clusters=2, **kw),
])
def test_measured_nonfinite_data_cannot_be_overridden(call):
    """Measured non-finite data is tier 1, so the override cannot rescue it.

    Matching on "cannot be overridden" distinguishes this from the tier-2
    advisory, which names ``allow_nonmetric=True`` as the remedy.
    """
    dist = oecluster.pdist(_mols(), "fingerprint")
    dist.condensed[0] = float('nan')
    with pytest.raises(ValueError, match="cannot be overridden"):
        call(dist, allow_nonmetric=True)


def test_measured_nonfinite_data_cannot_be_overridden_in_cluster_report():
    dist = oecluster.pdist(_mols(), "fingerprint")
    result = oecluster.butina(dist, 0.4)
    dist.condensed[0] = float('nan')
    with pytest.raises(ValueError, match="cannot be overridden"):
        oecluster.cluster_report(result, dist, allow_nonmetric=True)


def test_an_infinity_is_refused_too():
    """Non-finite, not NaN-only.

    This is not an over-refusal: the native predicates the stamp is built from
    are already non-finite (``!std::isfinite``), and the refusal message says
    "non-finite entries".
    """
    dist = oecluster.pdist(_mols(), "fingerprint")
    dist.condensed[1] = float('inf')
    assert dist.data_integrity == "complete"
    with pytest.raises(ValueError, match="non-finite entries"):
        oecluster.butina(dist, threshold=0.4)


def test_a_clean_matrix_still_clusters_through_every_entry_point():
    """Not an over-refusal: the scan must be silent on finite data.

    A measurement stricter than the data warrants breaks working code without
    a word, because no native error contradicts it.
    """
    dist = oecluster.pdist(_mols(), "fingerprint")
    assert _gate._has_nonfinite(dist) is False
    result = oecluster.butina(dist, 0.4)
    assert len(result.labels) == 6
    assert len(oecluster.dbscan(dist, 0.4, min_samples=2).labels) == 6
    assert len(oecluster.hdbscan(dist, min_cluster_size=2).labels) == 6
    assert len(oecluster.agglomerative(dist, n_clusters=2).labels) == 6
    assert oecluster.cluster_report(result, dist).num_samples == 6


def _sparse_matrix():
    """A clean six-item matrix held in ``SparseStorage``.

    ``cutoff`` alone is what selects the sparse backend
    (``__init__.py`` picks ``SparseStorage`` when ``cutoff > 0.0``); the
    unrelated ``storage=`` keyword is a fingerprint bit-vector option and says
    nothing about how the distances are stored.
    """
    sparse = oecluster.pdist(_mols(), "fingerprint", cutoff=0.9)
    assert isinstance(sparse.storage, oecluster.SparseStorage)
    return sparse


def _finalized_sparse(value):
    """A sparse matrix carrying ``value`` as a merged extra entry.

    ``SparseStorage::Set`` only appends to an unmerged per-thread buffer, so
    the write is invisible until ``Finalize`` -- public on every storage handle
    a caller already holds -- folds it into the entries the algorithms iterate.
    Neither ``-inf`` nor NaN compares greater than the cutoff, so neither is
    dropped by ``Set``'s ``value > cutoff_`` filter.
    """
    sparse = _sparse_matrix()
    before = len(sparse.storage._entries())
    sparse.storage.Set(3, 0, value)
    assert len(sparse.storage._entries()) == before
    sparse.storage.Finalize()
    assert len(sparse.storage._entries()) == before + 1
    assert not math.isfinite(sparse.storage.Get(3, 0))
    return sparse


def test_the_scan_reads_sparse_entries_and_leaves_a_clean_one_clusterable():
    """Not an over-refusal, and not a new cache either.

    The scan walks the merged entry list, which is exactly what
    ``ThresholdGraph`` iterates, so a clean sparse matrix must still reach
    ``butina`` and ``dbscan`` -- the only two entry points a sparse matrix can
    get to the gate through -- and cluster as before. The
    ``_condensed_cache`` assertion is the other half: the scan must not have
    densified this matrix behind the caller's back, because that cache would
    then be stale with respect to the storage.
    """
    sparse = _sparse_matrix()
    assert _gate._has_nonfinite(sparse) is False
    assert sparse._condensed_cache is None
    assert len(oecluster.butina(sparse, 0.4).labels) == 6
    assert len(oecluster.dbscan(sparse, 0.4, min_samples=2).labels) == 6
    assert sparse._condensed_cache is None


def test_a_sparse_matrix_survives_a_file_round_trip_as_sparse(tmp_path):
    """The round trip used to hand back a dense matrix full of fake zeros.

    ``to_file`` wrote ``.condensed``, which for sparse storage is a densified
    copy where every pair the cutoff omitted -- the *farthest* pairs -- reads
    as ``0.0``. Reloading therefore both silenced the ``SparseStorage``
    refusals and inverted the omitted distances.
    """
    path = tmp_path / "sparse.npz"
    sparse = _sparse_matrix()
    sparse.to_file(str(path))

    loaded = oecluster.load_distance_matrix(str(path))
    assert isinstance(loaded.storage, oecluster.SparseStorage)
    assert loaded.storage.Cutoff() == sparse.storage.Cutoff()
    # Verbatim, in order: ThresholdGraph iterates every tuple Entries() holds,
    # so an entry list that merely agrees pairwise would still cluster apart.
    assert loaded.storage._entries() == sparse.storage._entries()

    with pytest.raises(ValueError, match="SparseStorage is not supported"):
        oecluster.hdbscan(loaded, min_cluster_size=2)
    with pytest.raises(ValueError, match="SparseStorage is not supported"):
        oecluster.agglomerative(loaded, n_clusters=2)


def test_a_reloaded_sparse_matrix_clusters_identically(tmp_path):
    """The half of the defect that no gate could have caught.

    ``butina`` and ``dbscan`` legitimately accept sparse storage, so for them
    the round trip was not refused -- it was silently answered against
    distances where the most distant pairs had become the closest.
    """
    path = tmp_path / "sparse.npz"
    sparse = _sparse_matrix()
    sparse.to_file(str(path))
    loaded = oecluster.load_distance_matrix(str(path))

    assert (list(oecluster.butina(loaded, 0.6).labels)
            == list(oecluster.butina(sparse, 0.6).labels))
    assert (list(oecluster.dbscan(loaded, 0.6, min_samples=2).labels)
            == list(oecluster.dbscan(sparse, 0.6, min_samples=2).labels))


def _write_sparse_file(path, *, num_samples=3, cutoff=0.5, rows=(0,),
                       cols=(1,), values=(0.1,), index_dtype=np.int64):
    """Hand-build a sparse ``.npz``, bypassing ``to_file`` and its guards.

    The shapes these tests need cannot be produced by ``to_file`` any more,
    which is the point of M1 -- but ``from_file`` still reads untrusted input,
    so it has to keep refusing them. ``index_dtype`` and a sequence ``cutoff``
    vary three of the four sparse arrays ``to_file`` writes with a fixed dtype:
    the two index arrays by dtype, the cutoff by shape. The fourth,
    ``sparse_v``, this helper always writes as float64.
    """
    np.savez_compressed(
        str(path),
        comparison_name=np.array("fingerprint"),
        params_json=np.array(json.dumps({})),
        labels=np.array([]),
        num_samples=np.array(num_samples),
        storage_kind=np.array("sparse"),
        sparse_cutoff=np.array(cutoff, dtype=np.float64),
        sparse_i=np.array(rows, dtype=index_dtype),
        sparse_j=np.array(cols, dtype=index_dtype),
        sparse_v=np.array(values, dtype=np.float64),
    )


class _PokedSparseStorage(oecluster.SparseStorage):
    """A sparse storage whose entry list holds a pair ``Set`` will not store.

    The test below used to poke the diagonal in through ``Set`` itself, which
    accepted it: the ``i != j`` precondition was a bare ``assert``, compiled out
    of release builds. ``Set`` refuses the diagonal now, so no supported call
    sequence produces a malformed entry list and the only way left to reach
    ``to_file``'s guard is to hand it one.
    """

    def _entries(self):
        return [(2, 2, 0.11)]


def test_to_file_refuses_a_sparse_matrix_from_file_would_refuse(tmp_path):
    """A silent write followed by a hard read failure is the worst shape.

    A diagonal entry that reached ``Entries()`` used to be written happily --
    then refused by ``from_file`` on the next load, by which time the caller had
    thrown the matrix away. The refusal has to come at the write, and before the
    file exists. ``Set`` is now the first line of that defence and
    ``test_storage_set_refuses_an_out_of_range_index`` pins it there; this keeps
    the backstop honest for an entry list that arrives by any other route.
    """
    path = tmp_path / "poked.npz"
    storage = _PokedSparseStorage(6, 0.9)
    storage.Finalize()
    sparse = oecluster.SymmetricDistanceMatrix(storage, "fingerprint")
    with pytest.raises(ValueError,
                       match=r"not a pair of distinct indices below 6"):
        sparse.to_file(str(path))
    assert not path.exists()


def test_what_to_file_accepts_from_file_accepts(tmp_path):
    """The two sides run one shared rule, so they cannot drift apart.

    Round 7 put the check only on the load side, which made agreement a hope
    rather than a guarantee. Cutoffs are swept because the cutoff is the one
    file field the write side derives from the storage it is saving.
    """
    for cutoff in (0.1, 0.5, 0.9, 1.0):
        path = tmp_path / f"clean_{cutoff}.npz"
        sparse = oecluster.pdist(_mols(), "fingerprint", cutoff=cutoff)
        assert isinstance(sparse.storage, oecluster.SparseStorage)
        sparse.to_file(str(path))
        loaded = oecluster.load_distance_matrix(str(path))
        assert loaded.storage._entries() == sparse.storage._entries()


@pytest.mark.parametrize("value", [float('-inf'), float('nan')])
def test_poisoned_sparse_storage_round_trips_and_the_gate_still_refuses(
        tmp_path, value):
    """The write-side guard must not seize the refusal the metric gate owns.

    ``to_file`` validates before writing now, so this pins where that rule
    stops: at replay fidelity. A non-finite entry is a legal thing to store, so
    it is written, reloaded unchanged, and convicted by the gate -- with a
    message about the data -- rather than by the file format.
    """
    path = tmp_path / "poisoned.npz"
    sparse = _finalized_sparse(value)
    sparse.to_file(str(path))

    loaded = oecluster.load_distance_matrix(str(path))
    assert isinstance(loaded.storage, oecluster.SparseStorage)
    assert not math.isfinite(loaded.storage.Get(3, 0))
    assert len(loaded.storage._entries()) == len(sparse.storage._entries())

    with pytest.raises(ValueError, match="non-finite entries"):
        oecluster.butina(loaded, 0.5)
    with pytest.raises(ValueError, match="cannot be overridden"):
        oecluster.butina(loaded, 0.5, allow_nonmetric=True)


def test_one_nan_entry_does_not_switch_off_the_cutoff_guard(tmp_path):
    """``values.max()`` propagates NaN, and ``nan > cutoff`` is False.

    So a single NaN turned the whole guard off, and the over-cutoff value it
    was hiding was then dropped by ``Set``: the matrix loaded was not the
    matrix in the file, which is precisely what the guard exists to prevent.
    """
    path = tmp_path / "nan_and_over.npz"
    _write_sparse_file(path, num_samples=3, cutoff=0.5, rows=(0, 1),
                       cols=(1, 2), values=(float('nan'), 0.9))
    with pytest.raises(ValueError, match=r"sparse entry 1 has value 0\.9"):
        oecluster.load_distance_matrix(str(path))


def test_a_nan_below_the_cutoff_is_still_loaded(tmp_path):
    """Not an over-refusal: NaN stays a legal sparse entry value.

    Refusing it in the file format would be a second, earlier refusal for data
    the metric gate already handles, and would break the round trip the gate's
    own sparse tests depend on.
    """
    path = tmp_path / "nan_only.npz"
    _write_sparse_file(path, num_samples=3, cutoff=0.5, rows=(0, 1),
                       cols=(1, 2), values=(float('nan'), 0.2))
    loaded = oecluster.load_distance_matrix(str(path))
    entries = loaded.storage._entries()
    assert [(e[0], e[1]) for e in entries] == [(0, 1), (1, 2)]
    assert math.isnan(entries[0][2])
    with pytest.raises(ValueError, match="non-finite entries"):
        oecluster.butina(loaded, 0.5)


def test_a_negative_sparse_index_is_a_value_error(tmp_path):
    """``from_file`` documents ValueError; this shape reached SWIG instead.

    The old bound check tested only ``.max()``, so ``-1`` passed it, reached
    the replay loop and came back as an ``OverflowError`` about a ``size_t``
    argument of a method the caller never named.
    """
    path = tmp_path / "negative.npz"
    _write_sparse_file(path, num_samples=3, rows=(-1,), cols=(1,))
    with pytest.raises(ValueError, match=r"sparse entry 0 is \(-1, 1\)"):
        oecluster.load_distance_matrix(str(path))


def test_a_two_dimensional_sparse_index_array_is_a_value_error(tmp_path):
    """The other escape from the documented type, previously a TypeError."""
    path = tmp_path / "twod.npz"
    _write_sparse_file(path, num_samples=3, rows=((0,),), cols=((1,),),
                       values=((0.1,),))
    with pytest.raises(ValueError,
                       match="sparse_i must be a 1-D array, not 2-D"):
        oecluster.load_distance_matrix(str(path))


def test_a_float_sparse_index_array_is_a_value_error(tmp_path):
    """A wrong load with no exception at all, which is the worse failure.

    ``0.7`` is not negative, not at or beyond ``num_samples``, and not equal to
    ``1.9``, so every index rule passed and the replay loop's ``int(i)`` then
    silently made the pair ``(0, 1)``. The file said one matrix and the loaded
    object was another.
    """
    path = tmp_path / "float_index.npz"
    _write_sparse_file(path, num_samples=3, rows=(0.7,), cols=(1.9,),
                       index_dtype=np.float64)
    with pytest.raises(ValueError,
                       match="sparse_i must hold integer indices, not "
                             "dtype float64"):
        oecluster.load_distance_matrix(str(path))


def test_a_non_scalar_sparse_cutoff_is_a_value_error(tmp_path):
    """The cutoff escaped the documented type one line before the index arrays.

    ``float()`` on a two-element array raises ``TypeError``, so a caller
    catching the ``ValueError`` ``from_file`` documents saw the load abort
    through them instead.
    """
    path = tmp_path / "vector_cutoff.npz"
    _write_sparse_file(path, num_samples=3, cutoff=(0.5, 0.6))
    with pytest.raises(ValueError,
                       match="sparse_cutoff must be a scalar, not a 1-D array"):
        oecluster.load_distance_matrix(str(path))


def test_a_one_element_sparse_cutoff_array_is_a_value_error(tmp_path):
    """One element is still not the 0-d scalar ``to_file`` writes.

    ``float()`` raises ``TypeError`` on this shape too, so weakening the rule
    from ``ndim`` to ``size`` would reopen the same escape.
    """
    path = tmp_path / "one_element_cutoff.npz"
    _write_sparse_file(path, num_samples=3, cutoff=(0.5,))
    with pytest.raises(ValueError,
                       match="sparse_cutoff must be a scalar, not a 1-D array"):
        oecluster.load_distance_matrix(str(path))


@pytest.mark.parametrize("cutoff", [float('nan'), float('inf')])
def test_a_nan_or_positive_infinite_sparse_cutoff_still_loads(tmp_path, cutoff):
    """Not an over-refusal: the new cutoff rule is about shape, nothing else.

    These two replay into the same matrix -- ``0.1 > nan`` and ``0.1 > inf`` are
    both False, so nothing is dropped -- and refusing them would be a second
    gate on data the format has no quarrel with.
    """
    path = tmp_path / "odd_cutoff.npz"
    _write_sparse_file(path, num_samples=3, cutoff=cutoff)
    loaded = oecluster.load_distance_matrix(str(path))
    assert loaded.storage._entries() == [(0, 1, 0.1)]


def test_an_entry_above_a_negative_infinite_sparse_cutoff_is_refused(tmp_path):
    """An entry above a ``-inf`` cutoff is refused by the over-cutoff rule.

    ``0.1 > -inf`` holds, so ``Set`` would drop the entry and the file would
    replay as an empty matrix. The pre-existing over-cutoff rule catches that;
    asserting on its message records which rule answers, so relaxing the shape
    rule cannot be mistaken for relaxing this one.
    """
    path = tmp_path / "negative_inf_cutoff.npz"
    _write_sparse_file(path, num_samples=3, cutoff=float('-inf'))
    with pytest.raises(ValueError,
                       match=r"sparse entry 0 has value 0\.1, above the cutoff "
                             r"-inf"):
        oecluster.load_distance_matrix(str(path))


def test_a_dense_file_without_storage_kind_still_loads_as_dense(tmp_path):
    """Absence of the new keys means dense, so 4.x files load unchanged."""
    path = tmp_path / "legacy_dense.npz"
    np.savez_compressed(
        str(path),
        condensed=np.array([0.1, 0.2, 0.3], dtype=np.float64),
        comparison_name=np.array("fingerprint"),
        params_json=np.array(json.dumps({})),
        labels=np.array(["a", "b", "c"]),
        num_samples=np.array(3),
    )
    loaded = oecluster.load_distance_matrix(str(path))
    assert isinstance(loaded.storage, oecluster.DenseStorage)
    assert list(loaded.condensed) == [0.1, 0.2, 0.3]


def test_an_unknown_storage_kind_is_refused(tmp_path):
    """An unrecognised kind must not fall through to the dense path."""
    path = tmp_path / "future.npz"
    np.savez_compressed(
        str(path),
        condensed=np.array([0.1, 0.2, 0.3], dtype=np.float64),
        comparison_name=np.array("fingerprint"),
        num_samples=np.array(3),
        storage_kind=np.array("mmap"),
    )
    with pytest.raises(ValueError, match="Unknown storage_kind 'mmap'"):
        oecluster.load_distance_matrix(str(path))


def test_an_infinity_finalized_into_sparse_storage_is_refused_by_butina():
    """The write route a ``Set``-only measurement misses.

    ``Finalize`` promotes the unmerged write into ``Entries()``, and
    ``ThresholdGraph`` admits ``-inf`` as an edge because ``-inf <= threshold``
    holds -- so this poisons the clusters rather than being ignored. The stamp
    assertion comes first on purpose: it records that the stamp still says
    ``complete``, so the refusal can only have come from the measurement. The
    cache assertion records that measuring sparse entries does not densify.
    """
    sparse = _finalized_sparse(float('-inf'))
    assert sparse.data_integrity == "complete"
    with pytest.raises(ValueError, match="non-finite entries"):
        oecluster.butina(sparse, 0.5)
    assert sparse._condensed_cache is None


def test_an_infinity_finalized_into_sparse_storage_is_refused_by_dbscan():
    sparse = _finalized_sparse(float('-inf'))
    assert sparse.data_integrity == "complete"
    with pytest.raises(ValueError, match="non-finite entries"):
        oecluster.dbscan(sparse, 0.5, min_samples=2)


def test_a_nan_finalized_into_sparse_storage_is_refused():
    """NaN takes the same route as ``-inf`` but leaves no trace in the labels.

    A NaN entry is inert in ``ThresholdGraph``'s ``value <= threshold`` test,
    so unlike ``-inf`` it fabricates no edge: driving the native clusterer
    directly, past the gate, returns the same labels for the poisoned storage
    as for a clean one. There is consequently no wrong answer downstream for a
    caller to notice, which is why this asserts the refusal rather than a
    change in labels.
    """
    options = oecluster.ButinaOptions()
    options.distance_threshold = 0.5
    clean = _sparse_matrix()
    sparse = _finalized_sparse(float('nan'))
    assert (list(_butina_cluster(sparse.storage, options).Labels())
            == list(_butina_cluster(clean.storage, options).Labels()))

    assert sparse.data_integrity == "complete"
    with pytest.raises(ValueError, match="non-finite entries"):
        oecluster.butina(sparse, 0.5)
    with pytest.raises(ValueError, match="non-finite entries"):
        oecluster.dbscan(sparse, 0.5, min_samples=2)


@pytest.mark.parametrize("value", [float('-inf'), float('nan')])
@pytest.mark.parametrize("call", [
    lambda dm, **kw: oecluster.butina(dm, 0.5, **kw),
    lambda dm, **kw: oecluster.dbscan(dm, 0.5, min_samples=2, **kw),
])
def test_finalized_sparse_nonfinite_data_cannot_be_overridden(call, value):
    """Sparse non-finite data is tier 1, like every other non-finite route.

    Matching on "cannot be overridden" distinguishes this from the tier-2
    advisory, which names ``allow_nonmetric=True`` as the remedy.
    """
    sparse = _finalized_sparse(value)
    with pytest.raises(ValueError, match="cannot be overridden"):
        call(sparse, allow_nonmetric=True)


def test_an_unfinalized_sparse_set_still_clusters():
    """Not an over-refusal: the gate measures what the algorithms read.

    Without ``Finalize`` the written value stays in an unmerged per-thread
    buffer, invisible to ``Entries()``, to ``Get`` and therefore to
    ``ThresholdGraph``. Refusing on it would break a call the native layer runs
    correctly, so the labels must match the untouched matrix exactly.
    """
    expected = list(oecluster.butina(_sparse_matrix(), 0.5).labels)

    sparse = _sparse_matrix()
    sparse.storage.Set(3, 0, float('nan'))
    assert sparse.storage.Get(3, 0) == 0.0
    assert _gate._has_nonfinite(sparse) is False
    assert list(oecluster.butina(sparse, 0.5).labels) == expected
    assert len(oecluster.dbscan(sparse, 0.5, min_samples=2).labels) == 6
    assert sparse._condensed_cache is None


def test_all_unknown_facts_with_clean_data_still_cluster(tmp_path):
    """Not an over-refusal: the scan must not convict "unknown" on its own.

    A pre-5.0.0 file records no claim at all. The measurement answers only the
    finiteness question, so an all-unknown matrix holding finite values still
    clusters.
    """
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
    assert loaded.data_integrity == "unknown"
    assert loaded.metric_capabilities == {
        'zero_self': "unknown", 'triangle': "unknown"}
    assert _gate._has_nonfinite(loaded) is False
    assert len(oecluster.butina(loaded, 0.5).labels) == 3


def test_an_all_unknown_matrix_with_nonfinite_data_is_still_refused(tmp_path):
    """The complement: "unknown" buys no pass on data that is measurably bad."""
    path = tmp_path / "legacy_nan.npz"
    np.savez_compressed(
        str(path),
        condensed=np.array([0.1, float('nan'), 0.3], dtype=np.float64),
        comparison_name=np.array("fingerprint"),
        params_json=np.array(json.dumps({})),
        labels=np.array(["a", "b", "c"]),
        num_samples=np.array(3),
    )
    loaded = oecluster.load_distance_matrix(str(path))
    assert loaded.data_integrity == "unknown"
    with pytest.raises(ValueError, match="non-finite entries"):
        oecluster.butina(loaded, 0.5, allow_nonmetric=True)


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


def test_storage_get_refuses_an_out_of_range_index():
    """An unchecked Get answered with a wrong number instead of refusing.

    Which wrong number depended on the backend and the index. The three
    off-diagonal dense probes all read memory the one-element buffer does not
    own, though not all of it past the end: the condensed formula runs in
    ``size_t``, so ``(0, 2)`` and ``(1, 1399)`` index 1 and 1398 elements in,
    while ``(500, 1399)`` underflows and lands about 964 KiB *before* the base.
    Each returned whatever occupied that memory -- ``0.0`` and an uninitialised
    denormal have both been seen, so neither is the behaviour.
    ``SparseStorage::Get`` never reads out of bounds at all: it looks the pair
    up, and ``(0, 6)`` computes index 5, which is the stored ``(1, 2)`` pair,
    so it answered ``0.5`` -- a real distance between the wrong two items. The
    ``(1000, 1000)`` probes reach neither: the ``i == j`` shortcut answered
    ``0.0`` without computing an index, for a pair naming no item at all.
    """
    dense = oecluster.pdist(_mols()[:2], "fingerprint")
    assert dense.num_samples == 2
    assert dense.storage.Get(0, 1) == dense.condensed[0]
    for i, j in [(0, 2), (1, 1399), (500, 1399), (1000, 1000)]:
        with pytest.raises(RuntimeError, match="outside the storage range"):
            dense.storage.Get(i, j)

    sparse = _sparse_matrix()
    for i, j in [(0, 6), (1000, 1000)]:
        with pytest.raises(RuntimeError, match="outside the storage range"):
            sparse.storage.Get(i, j)


def test_storage_set_refuses_an_out_of_range_index():
    """The write half of the same contract, and the worse half of it.

    ``Get`` answered about a pair that does not exist; ``Set`` overwrote one
    that does, then returned as though it had stored what the caller asked for.
    On six samples the condensed formula sends ``(2, 6)`` to offset 12, which
    is the ``(3, 4)`` pair, and the diagonal ``(2, 2)`` to offset 8, which is
    ``(1, 5)`` -- the diagonal is in range on both indices yet owns no slot,
    because ``Get`` answers it from a shortcut rather than from memory. Further
    out there is no pair left to collide with: ``(0, 5000)`` wrote roughly
    40 KiB past a fifteen-element buffer without faulting.

    The sentinel is negative so that a surviving write cannot be mistaken for a
    distance: no metric in the library produces one.
    """
    dense = oecluster.pdist(_mols(), "fingerprint")
    assert dense.num_samples == 6
    aliased = [(3, 4), (1, 5)]
    before = [dense.storage.Get(i, j) for i, j in aliased]

    for i, j in [(2, 6), (500, 1399), (0, 5000), (1000, 1000)]:
        with pytest.raises(RuntimeError, match="outside the storage range"):
            dense.storage.Set(i, j, -1.0)
    with pytest.raises(RuntimeError, match="cannot store the diagonal pair"):
        dense.storage.Set(2, 2, -1.0)

    assert [dense.storage.Get(i, j) for i, j in aliased] == before

    # SparseStorage never wrote out of bounds -- it keeps (i, j) tuples -- but
    # its lookup collides the same way, so a stored bad pair would come back
    # later as a real one's distance.
    sparse = _sparse_matrix()
    with pytest.raises(RuntimeError, match="outside the storage range"):
        sparse.storage.Set(0, 6, -1.0)
    with pytest.raises(RuntimeError, match="cannot store the diagonal pair"):
        sparse.storage.Set(2, 2, -1.0)


def test_cluster_report_refuses_a_result_from_a_smaller_matrix():
    """The under-range direction, which every range check misses.

    A three-sample result scored against a six-sample matrix keeps every
    member index in range, so nothing native objects, and the scorecard comes
    back over distances the result was never computed from.
    """
    big = oecluster.pdist(_mols(), "fingerprint")
    small = oecluster.pdist(_mols()[:3], "fingerprint")
    result = oecluster.butina(small, 0.5)
    with pytest.raises(ValueError, match="result covers 3 samples"):
        oecluster.cluster_report(result, big)


def test_cluster_report_refuses_a_result_from_a_larger_matrix():
    """The over-range direction, refused before the native reader is reached.

    ``validate_cluster_members`` already caught this, but only after
    ``ClusterReport.cpp`` had read the out-of-range pairs.
    """
    big = oecluster.pdist(_mols(), "fingerprint")
    small = oecluster.pdist(_mols()[:3], "fingerprint")
    result = oecluster.butina(big, 0.5)
    with pytest.raises(ValueError, match="the matrix 3"):
        oecluster.cluster_report(result, small)


def test_cluster_report_names_the_bad_member_not_the_storage_class():
    """The bounds check on ``Get`` pre-empted the diagnosis the caller can act on.

    Sample counts agree here, so the wrapper's mismatch guard passes and member
    99 reached ``storage.Get`` before anything validated it -- answering with a
    ``DenseStorage`` complaint about a backend the caller never chose. Matching
    on "Cluster member index" is the point of the test: both messages carry
    "outside the storage range", which is why the regression went unnoticed.
    """
    d3 = oecluster.pdist(_mols()[:3], "fingerprint")
    bad = oecluster.ClusteringResult([0, 0, 0], [(0, 1, 99)])
    with pytest.raises(RuntimeError, match="Cluster member index") as excinfo:
        oecluster.cluster_report(bad, d3)
    assert "DenseStorage" not in str(excinfo.value)
    # The representative entry points always validated first; they are asserted
    # alongside to pin that both paths now give the same domain message.
    with pytest.raises(RuntimeError, match="Cluster member index"):
        oecluster.representative((0, 1, 99), d3)


def test_cluster_report_sample_mismatch_outranks_the_storage_refusal():
    """A mismatched pairing is the error the caller has to fix first.

    Whether the matrix is sparse, or non-metric, is not the caller's problem
    when it is the wrong matrix -- so the mismatch must not be pre-empted by
    either of the refusals that follow it.
    """
    small = oecluster.pdist(_mols()[:3], "fingerprint", metric="dice")
    result = oecluster.butina(small, 0.5, allow_nonmetric=True)
    with pytest.raises(ValueError, match="same items"):
        oecluster.cluster_report(result, _nonmetric_sparse())
    with pytest.raises(ValueError, match="same items"):
        oecluster.cluster_report(result, _nonmetric())


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


def test_require_comparable_refuses_a_similarity():
    dist = oecluster.pdist(_mols(), "fingerprint")
    dist._facts['is_distance'] = False

    with pytest.raises(ValueError, match="requires distances"):
        _gate.require_comparable(dist, "activity_landscape")


def test_require_comparable_refuses_a_nonzero_self_distance():
    dist = oecluster.pdist(_mols(), "fingerprint")
    dist._facts['zero_self'] = False

    with pytest.raises(ValueError, match="zero self-distance"):
        _gate.require_comparable(dist, "activity_landscape")


def test_require_comparable_refuses_a_non_finite_entry():
    dist = oecluster.pdist(_mols(), "fingerprint")
    dist.condensed[0] = math.nan

    with pytest.raises(ValueError, match="non-finite entries"):
        _gate.require_comparable(dist, "modelability")


def test_require_comparable_refuses_subset_scored_without_an_override():
    """The one tier-2 check that still bites. Ranking a nearest neighbour or
    thresholding a cliff compares two distances against each other, and under
    missing='ignore' the two answer different questions."""
    dist = oecluster.pdist(_mols(), "fingerprint")
    dist._facts['data_integrity'] = "subset_scored"

    with pytest.raises(ValueError) as excinfo:
        _gate.require_comparable(dist, "activity_landscape")

    message = str(excinfo.value)
    assert "not mutually comparable" in message
    assert "missing='complete_case'" in message
    assert "cannot be overridden" in message
    assert "allow_nonmetric" not in message


def test_require_comparable_has_no_override_keyword():
    """The refusal is non-overridable as behaviour, not merely as wording.

    ``test_..._refuses_subset_scored_without_an_override`` asserts the message
    does not mention ``allow_nonmetric``; this asserts the function really has
    no such parameter, so the refusal cannot be switched off by a caller who
    guesses the keyword from ``require_metric``.
    """
    parameters = inspect.signature(_gate.require_comparable).parameters
    assert "allow_nonmetric" not in parameters

    dist = oecluster.pdist(_mols(), "fingerprint")
    dist._facts['data_integrity'] = "subset_scored"
    # The call below is intentionally invalid to verify the function rejects
    # an override parameter rather than silently accepting it.
    with pytest.raises(TypeError):
        _gate.require_comparable(dist, "activity_landscape",
                                 allow_nonmetric=True)  # pyright: ignore[reportCallIssue]


def test_require_comparable_accepts_a_triangle_violation():
    """These metrics never assume a metric: they rank and threshold distances,
    and a triangle-inequality violation leaves both operations meaningful.

    Asserts that the gate accepts the violation, then proves the gate is
    actually enforcing all three tier-1 checks by injecting each fault into
    a fresh matrix with the waived fact present."""
    dist = oecluster.pdist(_mols(), "fingerprint")
    dist._facts['triangle'] = False

    assert _gate.require_comparable(dist, "activity_landscape") is None

    for name, apply_fault, message in _TIER1_FAULTS:
        dist = oecluster.pdist(_mols(), "fingerprint")
        dist._facts['triangle'] = False
        apply_fault(dist)
        with pytest.raises(ValueError, match=message):
            _gate.require_comparable(dist, "activity_landscape")


def test_require_comparable_accepts_a_proven_probe_violation():
    """These metrics never assume a metric: they rank and threshold distances,
    and a proven triangle inequality violation (via probe sampling) leaves both
    operations meaningful.

    Asserts that the gate accepts the violation, then proves the gate is
    actually enforcing all three tier-1 checks by injecting each fault into
    a fresh matrix with the waived fact present."""
    dist = oecluster.pdist(_mols(), "fingerprint")
    dist._facts['metric_probe'] = "violations_found"
    dist._facts['probe_violations'] = 3
    dist._facts['probe_sampled'] = 100

    assert _gate.require_comparable(dist, "modelability") is None

    for name, apply_fault, message in _TIER1_FAULTS:
        dist = oecluster.pdist(_mols(), "fingerprint")
        dist._facts['metric_probe'] = "violations_found"
        dist._facts['probe_violations'] = 3
        dist._facts['probe_sampled'] = 100
        apply_fault(dist)
        with pytest.raises(ValueError, match=message):
            _gate.require_comparable(dist, "modelability")
