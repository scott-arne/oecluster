import numpy as np
import oecluster
import pytest
from oecluster import _comparisons
from oecluster import oecluster as _native
from openeye import oechem

SMILES = ["CCO", "CCC", "CCCC", "c1ccccc1", "CCN", "CCOC", "CC(=O)O", "CCCCCC"]

# Spec section 2.4's measured missingness set: 19 molecules, of which OpenEye
# cannot assign XLogP types to the last 6.
MISSINGNESS = [
    ("benzene", "c1ccccc1"),
    ("phenol", "c1ccc(O)cc1"),
    ("octane", "CCCCCCCC"),
    ("aspirin", "CC(=O)Oc1ccccc1C(=O)O"),
    ("procaine", "CCN(CC)CCOC(=O)c1ccccc1"),
    ("pyridine", "c1ccncc1"),
    ("ethanol", "CCO"),
    ("toluene", "Cc1ccccc1"),
    ("aniline", "Nc1ccccc1"),
    ("caffeine", "Cn1cnc2c1c(=O)n(C)c(=O)n2C"),
    ("acetaminophen", "CC(=O)Nc1ccc(O)cc1"),
    ("ibuprofen", "CC(C)Cc1ccc(cc1)C(C)C(=O)O"),
    ("naphthalene", "c1ccc2ccccc2c1"),
    ("sodium", "[Na+]"),
    ("iron", "[Fe]"),
    ("helium", "[He]"),
    ("platinum", "[Pt]"),
    ("uranium", "[U]"),
    ("silicon", "[Si]"),
]

UNTYPED = ["sodium", "iron", "helium", "platinum", "uranium", "silicon"]


def _mols(smiles_list=None):
    mols = []
    for idx, smi in enumerate(smiles_list or SMILES):
        mol = oechem.OEGraphMol()
        oechem.OESmilesToMol(mol, smi)
        mol.SetTitle(f"mol{idx}")
        mols.append(mol)
    return mols


def _missingness_mols():
    mols = []
    for title, smi in MISSINGNESS:
        mol = oechem.OEGraphMol()
        oechem.OESmilesToMol(mol, smi)
        mol.SetTitle(title)
        mols.append(mol)
    return mols


def test_descriptor_statistics_defaults_to_the_openeye_source():
    stats = oecluster.descriptor_statistics(_mols())
    assert stats['num_rows'] == 8
    assert stats['columns']
    assert len(stats['mean']) == len(stats['columns'])
    assert len(stats['variance']) == len(stats['columns'])
    assert len(stats['present_count']) == len(stats['columns'])


def test_descriptor_statistics_drops_zero_variance_columns():
    """Adding water drops a column the clean set keeps.

    The clean set drops nothing, so the drop report only becomes an oracle
    once an input that provokes a drop is used.
    """
    clean = oecluster.descriptor_statistics(_mols())
    stats = oecluster.descriptor_statistics(_mols(["O"] + SMILES))

    assert clean['dropped'] == []
    assert "FractionCsp3" in clean['columns']

    assert stats['dropped'] == [("FractionCsp3", "zero-variance")]
    assert "FractionCsp3" not in stats['columns']
    assert all(v > 0.0 for v in stats['variance'])


def test_descriptor_statistics_can_return_an_inverse_covariance():
    stats = oecluster.descriptor_statistics(_mols(), inverse_covariance=True)
    k = len(stats['columns'])
    assert stats['inverse_covariance'].shape == (k, k)
    # Pinned rather than bounded by k: the rank of a k x k matrix is at most k
    # by construction, so ``rank <= k`` holds whatever the wrapper reports.
    assert k == 11
    assert stats['inverse_covariance_rows'] == 8
    assert stats['inverse_covariance_rank'] == 7


def test_descriptor_statistics_reports_the_covariance_row_count():
    """Listwise deletion fits the covariance over fewer rows than the columns.

    OpenEye assigns sodium no XLogP value, and XLogP survives the
    zero-variance drop here, so the sodium row reaches the per-column
    statistics but not the covariance. A caller reading ``num_rows`` beside
    the matrix would overstate what produced it.
    """
    mols = _mols(["[Na+]"] + SMILES)
    stats = oecluster.descriptor_statistics(mols, inverse_covariance=True)

    assert stats['num_rows'] == 9
    assert stats['inverse_covariance_rows'] == 8

    assert "XLogP" in stats['columns']
    xlogp_present = dict(zip(stats['columns'], stats['present_count']))["XLogP"]
    assert xlogp_present == 8


def test_descriptor_statistics_skips_inverse_covariance_by_default():
    stats = oecluster.descriptor_statistics(_mols())
    assert stats['inverse_covariance'] is None
    assert stats['inverse_covariance_rank'] == 0
    assert stats['inverse_covariance_rows'] == 0


def test_descriptor_statistics_needs_two_molecules_to_fit_a_variance():
    with pytest.raises(RuntimeError, match="at least two molecules"):
        oecluster.descriptor_statistics(_mols(["CCO"]))


def test_descriptor_statistics_rejects_an_unknown_source():
    # SWIG maps every OEClusterError, including ComparisonError, to
    # RuntimeError, so C++-side rejections surface as RuntimeError.
    with pytest.raises(RuntimeError, match="nope"):
        oecluster.descriptor_statistics(_mols(), sources=["nope"])


def test_descriptor_statistics_rejects_an_unknown_column():
    with pytest.raises(RuntimeError):
        oecluster.descriptor_statistics(_mols(), columns=["not_a_column"])


def test_pdist_descriptor_uses_standardized_euclidean_by_default():
    dist = oecluster.pdist(_mols(), "descriptor")
    assert dist.num_samples == 8
    assert dist.comparison_name == "descriptor"
    assert np.all(np.isfinite(dist.condensed))
    assert np.any(dist.condensed > 0.0)


def test_pdist_descriptor_is_stamped_as_a_metric():
    dist = oecluster.pdist(_mols(), "descriptor")
    assert dist.metric_capabilities == {'zero_self': True, 'triangle': True}
    assert dist.data_integrity == "complete"
    oecluster.butina(dist, 1.0)


def test_pdist_descriptor_accepts_mahalanobis():
    dist = oecluster.pdist(_mols(), "descriptor", metric="mahalanobis")
    assert np.all(np.isfinite(dist.condensed))


def test_pdist_descriptor_accepts_minkowski_with_p():
    dist = oecluster.pdist(_mols(), "descriptor", metric="minkowski", p=3.0)
    assert np.all(np.isfinite(dist.condensed))


def test_explicit_variances_bypass_the_pooled_fit():
    fitted = oecluster.pdist(_mols(), "descriptor")
    stats = oecluster.descriptor_statistics(_mols())
    override = oecluster.pdist(_mols(), "descriptor",
                               columns=stats['columns'],
                               variances=stats['variance'])
    np.testing.assert_allclose(override.condensed, fitted.condensed, rtol=1e-9)


def test_propagate_is_stamped_nan_present_and_refused():
    dist = oecluster.pdist(_mols(), "descriptor", missing="propagate")
    assert dist.data_integrity == "nan_present"
    with pytest.raises(ValueError, match="cannot be overridden"):
        oecluster.butina(dist, 1.0, allow_nonmetric=True)


def test_ignore_is_stamped_subset_scored_and_overridable():
    dist = oecluster.pdist(_mols(), "descriptor", metric="euclidean",
                           missing="ignore")
    assert dist.data_integrity == "subset_scored"
    with pytest.raises(ValueError, match="subset"):
        oecluster.butina(dist, 1.0)
    oecluster.butina(dist, 1.0, allow_nonmetric=True)


def test_an_invalid_missing_policy_is_reported_as_invalid():
    """An explicit policy is never folded onto the default.

    Both inputs get the same verdict on the same bad argument. Were the
    normalizer to read "" as unspecified, the first call would filter under
    complete-case, empty its own input, and be refused for a shape problem the
    caller does not have.
    """
    untyped = _mols(["[Na+]", "[Fe]", "[He]"])
    with pytest.raises(RuntimeError, match="Unknown missing-value policy"):
        oecluster.pdist(untyped, "descriptor", missing="")
    with pytest.raises(RuntimeError, match="Unknown missing-value policy"):
        oecluster.pdist(_mols(), "descriptor", missing="")


def test_similarity_is_rejected_for_descriptors():
    with pytest.raises(ValueError, match="no similarity form"):
        oecluster.pdist(_mols(), "descriptor", similarity=True)


def test_unknown_descriptor_kwargs_are_rejected():
    with pytest.raises(TypeError, match="Unknown kwargs for descriptor"):
        oecluster.pdist(_mols(), "descriptor", bogus=1)


def test_complete_data_excludes_nothing():
    dist = oecluster.pdist(_mols(), "descriptor")
    assert 'excluded_items' not in dist.params
    assert len(dist.labels) == 8


def test_complete_case_excludes_the_untyped_species_end_to_end():
    """The real oracle for spec section 2.4, with nothing monkeypatched.

    19 molecules in, 6 of which OpenEye cannot assign XLogP types to, so a
    13x13 matrix comes out and it is stamped ``complete``: the exclusions
    happened before scoring, so nothing that survived is missing anything.
    Asserted by title, because an assertion on the count alone still passes
    when the mask lands on the wrong rows.
    """
    mols = _missingness_mols()
    dist = oecluster.pdist(mols, "descriptor")

    assert dist.num_samples == 13
    assert len(dist.condensed) == 13 * 12 // 2
    assert dist.labels == [title for title, _ in MISSINGNESS
                           if title not in UNTYPED]
    excluded = [index for index, _ in dist.params['excluded_items']]
    assert [MISSINGNESS[index][0] for index in excluded] == UNTYPED
    assert all(reason == "missing-descriptor"
               for _, reason in dist.params['excluded_items'])

    assert dist.data_integrity == "complete"
    assert np.all(np.isfinite(dist.condensed))
    oecluster.butina(dist, 1.0)


def test_a_present_but_non_finite_value_is_still_excluded():
    """Validity bit set, value NaN: the case a mask-only filter would admit.

    OpenEye reports ``FractionCsp3`` for water as present and NaN -- there is
    no carbon to be sp3 or otherwise -- while every validity bit on the row is
    set. Only the finiteness half of the check catches it, and admitting the
    row would put NaN into the distance matrix.
    """
    mols = _mols(["O"] + SMILES)
    excluded = list(_native.descriptor_excluded_indices(
        mols, _native.DescriptorOptions()))
    assert excluded == [0]

    dist = oecluster.pdist(mols, "descriptor")
    assert dist.num_samples == 8
    assert dist.labels == [f"mol{i}" for i in range(1, 9)]
    assert dist.params['excluded_items'] == [[0, "missing-descriptor"]]
    assert np.all(np.isfinite(dist.condensed))


def test_the_rdkit_bcut_columns_are_present_and_nan():
    """The same failure on the source the spec names, where it is total.

    RDKit's eight ``BCUT2D_*`` columns come back present-and-NaN for a
    single-atom species, and the whole 213-column RDKit matrix reports no
    absent values at all for this input. A filter reading the validity mask
    alone would find nothing to exclude here.
    """
    mols = _mols(["[Na+]"] + SMILES)
    excluded = list(_native.descriptor_excluded_indices(
        mols, _comparisons.descriptor_options({'sources': ["rdkit"]})))
    assert excluded == [0]

    dist = oecluster.pdist(mols, "descriptor", sources=["rdkit"])
    assert dist.num_samples == 8
    assert dist.params['excluded_items'] == [[0, "missing-descriptor"]]
    assert dist.data_integrity == "complete"
    assert np.all(np.isfinite(dist.condensed))


def test_complete_case_filtering_drops_and_records(monkeypatch):
    """The normalizer must remove the masked rows before labels are taken."""
    monkeypatch.setattr(_native, "descriptor_excluded_indices",
                        lambda items, opts: [1, 3])
    dist = oecluster.pdist(_mols(), "descriptor")
    assert dist.num_samples == 6
    assert dist.labels == ["mol0", "mol2", "mol4", "mol5", "mol6", "mol7"]
    assert dist.params['excluded_items'] == [[1, "missing-descriptor"],
                                             [3, "missing-descriptor"]]


def test_complete_case_filtering_is_skipped_for_other_policies(monkeypatch):
    monkeypatch.setattr(_native, "descriptor_excluded_indices",
                        lambda items, opts: [1, 3])
    dist = oecluster.pdist(_mols(), "descriptor", missing="propagate")
    assert dist.num_samples == 8


def test_pdist_raises_when_filtering_empties_the_input(monkeypatch):
    """Filtering everything out must blame the filtering, not the input.

    The caller passed eight molecules; the error has to say normalization
    emptied the list rather than report it as arriving empty.
    """
    monkeypatch.setattr(_native, "descriptor_excluded_indices",
                        lambda items, opts: list(range(len(items))))
    with pytest.raises(ValueError, match="normalizing the inputs"):
        oecluster.pdist(_mols(), "descriptor")


def test_pdist_still_returns_a_zero_sample_matrix_for_an_empty_input():
    """An input that arrived empty is not the emptied-by-filtering case.

    Nothing was excluded here, so the new guard must stand aside and leave
    this call returning the 0-sample matrix it has always returned.
    """
    dist = oecluster.pdist([], "fingerprint")
    assert dist.num_samples == 0
    assert len(dist.condensed) == 0


def test_an_empty_descriptor_input_is_not_blamed_on_the_filtering():
    """The same standing-aside, on the comparison that does normalize.

    The normalizer runs and excludes nothing, so the refusal comes from the
    C++ fit and describes what it needs instead of naming normalization.
    """
    with pytest.raises(RuntimeError, match="at least two molecules"):
        oecluster.pdist([], "descriptor")


def test_cdist_descriptor_filters_a_real_mask_end_to_end():
    """cdist with the computed mask, nothing monkeypatched."""
    a = _mols(["[Na+]", "CCO", "CCC", "CCCC"])
    b = _mols(SMILES[3:])
    cross = oecluster.cdist(a, b, "descriptor")

    assert cross.shape == (3, 5)
    assert cross.params['excluded_items_a'] == [[0, "missing-descriptor"]]
    assert 'excluded_items_b' not in cross.params
    assert cross.labels_a == ["mol1", "mol2", "mol3"]
    assert np.isfinite(np.asarray(cross)).all()


def test_cdist_filters_each_side_independently(monkeypatch):
    calls = []

    def fake_mask(items, opts):
        calls.append(len(items))
        return [0] if len(items) == 3 else []

    monkeypatch.setattr(_native, "descriptor_excluded_indices", fake_mask)
    cross = oecluster.cdist(_mols()[:3], _mols()[3:], "descriptor")
    assert calls == [3, 5]
    assert cross.shape == (2, 5)
    assert cross.params['excluded_items_a'] == [[0, "missing-descriptor"]]
    assert 'excluded_items_b' not in cross.params


def test_cdist_filters_the_b_side_only(monkeypatch):
    """Exclusions on B are indexed against B's own original list.

    B here is ``_mols()[3:]``, so its index 2 is "mol5". A mask read against
    the combined A+B list would have dropped "mol2" instead, which is not even
    in B.
    """
    monkeypatch.setattr(
        _native, "descriptor_excluded_indices",
        lambda items, opts: [2] if len(items) == 5 else [])
    cross = oecluster.cdist(_mols()[:3], _mols()[3:], "descriptor")
    assert cross.shape == (3, 4)
    assert 'excluded_items_a' not in cross.params
    assert cross.params['excluded_items_b'] == [[2, "missing-descriptor"]]
    assert cross.labels_b == ["mol3", "mol4", "mol6", "mol7"]


def test_cdist_filters_both_sides(monkeypatch):
    """Both masks apply and the Size() == n_a + n_b guard still passes."""
    monkeypatch.setattr(
        _native, "descriptor_excluded_indices",
        lambda items, opts: [0] if len(items) == 3 else [1, 4])
    cross = oecluster.cdist(_mols()[:3], _mols()[3:], "descriptor")
    assert cross.shape == (2, 3)
    assert cross.params['excluded_items_a'] == [[0, "missing-descriptor"]]
    assert cross.params['excluded_items_b'] == [[1, "missing-descriptor"],
                                                [4, "missing-descriptor"]]
    assert np.isfinite(np.asarray(cross)).all()


def test_cdist_raises_when_filtering_empties_a_side(monkeypatch):
    """Filtering a side down to nothing must blame normalization, not the input.

    The caller passed three molecules; the error has to say the filtering
    emptied set A, not that set A arrived empty.
    """
    monkeypatch.setattr(
        _native, "descriptor_excluded_indices",
        lambda items, opts: list(range(len(items))) if len(items) == 3 else [])
    with pytest.raises(ValueError, match="normalizing the inputs"):
        oecluster.cdist(_mols()[:3], _mols()[3:], "descriptor")


def test_overflow_during_scoring_downgrades_the_stamp():
    """Facts are read after the computation, never from the builder.

    Every descriptor here is present and finite; only the Minkowski
    accumulator overflows while pairs are scored. A gate that read facts at
    build time would stamp ``complete`` and admit an all-infinite matrix.
    """
    dist = oecluster.pdist(_mols(), "descriptor", metric="minkowski", p=400.0)
    assert dist.data_integrity == "nan_present"
    with pytest.raises(ValueError, match="cannot be overridden"):
        oecluster.butina(dist, 1.0)
    with pytest.raises(ValueError, match="cannot be overridden"):
        oecluster.butina(dist, 1.0, allow_nonmetric=True)


def test_descriptor_options_does_not_consume_kwargs():
    kwargs = {'metric': "euclidean", 'missing': "propagate"}
    _comparisons.descriptor_options(kwargs)
    assert kwargs == {'metric': "euclidean", 'missing': "propagate"}


def test_the_factory_class_builds_a_usable_comparison():
    comparison = oecluster.DescriptorComparison(_mols(), metric="euclidean")
    dist = oecluster.pdist(_mols(), comparison)
    assert dist.num_samples == 8
    assert dist.metric_capabilities['triangle'] is True


def test_the_factory_class_refuses_what_pdist_would_have_filtered():
    """Applying no filtering means refusing, not admitting.

    ``pdist(mols, "descriptor")`` drops the water row and scores the other
    eight. The factory keeps the caller's list intact, so the same input meets
    the default complete-case policy head-on.
    """
    mols = _mols(["O"] + SMILES)
    assert oecluster.pdist(mols, "descriptor").num_samples == 8

    with pytest.raises(RuntimeError, match="absent or non-finite"):
        oecluster.DescriptorComparison(mols)

    comparison = oecluster.DescriptorComparison(mols, missing="propagate")
    assert comparison.Size() == 9
