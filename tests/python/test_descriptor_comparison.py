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


def test_a_constant_column_is_dropped_as_zero_variance():
    """Six saturated non-aromatics make two columns genuinely constant."""
    saturated = ["CCO", "CCC", "CCCC", "CCN", "CCOC", "CCCCCC"]
    stats = oecluster.descriptor_statistics(_mols(saturated))

    assert stats['dropped'] == [("FractionCsp3", "zero-variance"),
                                ("AromaticRingCount", "zero-variance")]
    assert "FractionCsp3" not in stats['columns']
    assert "AromaticRingCount" not in stats['columns']
    assert all(v > 0.0 for v in stats['variance'])

    # A dropped column is absent from the reported arrays, so the report alone
    # cannot show constancy. Spread over the column by itself can.
    for name in ("FractionCsp3", "AromaticRingCount"):
        spread = oecluster.pdist(_mols(saturated), "descriptor",
                                 columns=[name], metric="euclidean")
        assert np.max(np.abs(spread.condensed)) == 0.0


def test_a_present_but_non_finite_value_drops_its_whole_column():
    """One NaN value costs the column, under the constant-column label.

    OpenEye reports water's FractionCsp3 as present and NaN, and a non-finite
    variance is dropped by the same branch as a zero one. The reason string is
    therefore the one ``test_a_constant_column_is_dropped_as_zero_variance``
    gets, even though this column varies perfectly well over the eight
    organics on its own. A caller reading the drop report cannot tell the two
    causes apart.
    """
    clean = oecluster.descriptor_statistics(_mols())
    clean_variance = dict(zip(clean['columns'], clean['variance']))
    stats = oecluster.descriptor_statistics(_mols(["O"] + SMILES))

    assert clean['dropped'] == []
    assert clean_variance["FractionCsp3"] > 0.0

    assert stats['dropped'] == [("FractionCsp3", "zero-variance")]
    assert "FractionCsp3" not in stats['columns']
    assert all(v > 0.0 for v in stats['variance'])


def test_descriptor_statistics_can_return_an_inverse_covariance():
    stats = oecluster.descriptor_statistics(_mols(), inverse_covariance=True)
    k = len(stats['columns'])
    assert stats['inverse_covariance'].shape == (k, k)
    # A ``rank <= k`` bound admitted every value from 0 to k, including the 0
    # the default path reports when nothing was fitted at all. Pin instead.
    assert k == 11
    assert stats['inverse_covariance_rows'] == 8
    assert stats['inverse_covariance_rank'] == 7


def test_descriptor_statistics_reports_the_covariance_row_count():
    """An absent value costs its row in the covariance but not elsewhere.

    OpenEye assigns sodium no XLogP value, and XLogP survives the
    zero-variance drop here, so the sodium row reaches the per-column
    statistics but is deleted listwise from the covariance. A caller reading
    ``num_rows`` beside the matrix would overstate what produced it.

    This is the absent-value branch only. Water is also missing a selected
    descriptor, but as a present NaN, which costs the whole column before the
    covariance is fitted and so leaves no row to delete. The two shapes are
    asserted side by side because the row count diverges between them.
    """
    absent = oecluster.descriptor_statistics(
        _mols(["[Na+]"] + SMILES), inverse_covariance=True)
    non_finite = oecluster.descriptor_statistics(
        _mols(["O"] + SMILES), inverse_covariance=True)

    assert absent['num_rows'] == 9
    assert absent['inverse_covariance_rows'] == 8

    assert non_finite['num_rows'] == 9
    assert non_finite['inverse_covariance_rows'] == 9

    assert "XLogP" in absent['columns']
    xlogp_present = dict(zip(absent['columns'], absent['present_count']))["XLogP"]
    assert xlogp_present == 8


def test_reusing_variances_needs_the_columns_they_were_fitted_over():
    """A dropped column desynchronizes ``variances=`` from the default selection.

    ``descriptor_statistics`` reports over surviving columns only, so once
    anything is dropped its ``variance`` list is shorter than the selection
    ``pdist`` makes on its own, and the two are matched by position. Dropping
    is by design, so this is the ordinary case, not an exotic one.
    """
    mols = _mols(["O"] + SMILES)
    stats = oecluster.descriptor_statistics(mols)

    assert stats['dropped'] == [("FractionCsp3", "zero-variance")]
    assert len(stats['columns']) == 10

    with pytest.raises(RuntimeError, match="11 columns are selected"):
        oecluster.pdist(mols, "descriptor", variances=stats['variance'])

    paired = oecluster.pdist(mols, "descriptor",
                             columns=stats['columns'],
                             variances=stats['variance'])
    assert paired.num_samples == 9


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
    """Both inputs get the same verdict on the same bad argument.

    The second input is one the complete-case filter empties. No input can
    make ``similarity=True`` valid here, so a complaint about the emptied list
    would describe a filter the caller never asked for instead of the argument
    they have to change.
    """
    for mols in (_mols(), _mols(["[Na+]", "[Fe]", "[He]"])):
        with pytest.raises(ValueError, match="no similarity form"):
            oecluster.pdist(mols, "descriptor", similarity=True)


def test_unknown_descriptor_kwargs_are_rejected():
    """Both inputs get the same verdict, message included."""
    for mols in (_mols(), _mols(["[Na+]", "[Fe]", "[He]"])):
        with pytest.raises(
                TypeError,
                match=r"Unknown kwargs for descriptor comparison: \['bogus'\]"):
            oecluster.pdist(mols, "descriptor", bogus=1)


def test_cdist_validates_arguments_before_filtering_either_side():
    """The precedence holds on the rectangular path, whichever side empties."""
    untyped = _mols(["[Na+]", "[Fe]", "[He]"])

    for a, b in ((untyped, _mols()), (_mols(), untyped)):
        with pytest.raises(ValueError, match="no similarity form"):
            oecluster.cdist(a, b, "descriptor", similarity=True)
        with pytest.raises(
                TypeError,
                match=r"Unknown kwargs for descriptor comparison: \['bogus'\]"):
            oecluster.cdist(a, b, "descriptor", bogus=1)


# Every option mistake no molecule set can rescue, paired with the message it
# has to produce. The last six need the descriptor schema, which is resolved
# from ``sources`` alone -- that is why the boundary is "needs no molecules"
# rather than "needs no schema", and it is what puts them ahead of a filter
# that can empty the item list. ``columns`` is named explicitly wherever an
# entry's value is what is wrong, so the case does not depend on how many
# columns the default selection happens to hold.
UNRESCUABLE = [
    ({'metric': "bogus"}, "Unknown metric 'bogus'"),
    ({'metric': "euclidean", 'variances': [1.0]},
     "variances applies only to metric='standardized_euclidean'"),
    ({'variances': [1.0, 2.0]},
     r"variances has 2 entries but \d+ columns are selected"),
    ({'metric': "mahalanobis", 'inverse_covariance': [1.0, 2.0, 3.0]},
     r"inverse_covariance has 3 entries but \d+ columns are selected"),
    ({'columns': ["MolecularWeight", "XLogP"], 'variances': [0.0, 1.0]},
     (r"variances\[0\] for column 'MolecularWeight' must be finite and "
      r"strictly positive")),
    ({'columns': ["MolecularWeight", "XLogP"],
      'variances': [float("nan"), 1.0]},
     r"variances\[0\] for column 'MolecularWeight' .* got nan"),
    ({'metric': "mahalanobis", 'columns': ["MolecularWeight", "XLogP"],
      'inverse_covariance': [float("nan"), 0.0, 0.0, 1.0]},
     r"inverse_covariance\[0\] must be finite, got nan"),
    ({'columns': ["XLogP", "MolecularWeight"], 'variances': [1.5, 2.5]},
     "columns must be in ascending schema order"),
]


@pytest.mark.parametrize(("kwargs", "expected"), UNRESCUABLE,
                         ids=[str(case[0]) for case in UNRESCUABLE])
def test_an_unusable_option_value_outranks_the_emptied_input(kwargs, expected):
    """A value no input can rescue is refused before the filter runs.

    Argument *names* were already checked ahead of the filter; values were not,
    so an emptied input answered for ``metric='bogus'`` with a complaint about
    a filter the caller never asked for. The clean input is what shows the
    message is unchanged, and the emptied one is what shows it now wins.

    The ``variances has N entries`` case is the one whose remedy is the
    documented ``descriptor_statistics`` recipe, so suppressing it withheld the
    fix from the caller who was part-way through following the documentation.
    """
    untyped = _mols(["[Na+]", "[Fe]", "[He]"])
    clean = _mols()

    for mols in (clean, untyped):
        with pytest.raises(RuntimeError, match=expected):
            oecluster.pdist(mols, "descriptor", **kwargs)

    for a, b in ((untyped, clean), (clean, untyped)):
        with pytest.raises(RuntimeError, match=expected):
            oecluster.cdist(a, b, "descriptor", **kwargs)


def test_the_option_check_reaches_python_without_molecules():
    """The pre-filter check is the C++ one, reached through SWIG.

    Python cannot reproduce the metric table -- ``resolve_metric`` lives in a
    private header -- so the validator has to call into C++ with the options
    alone. Calling it here with no molecules at all is what shows the check is
    genuinely molecule-independent rather than merely early.

    The override case is the one that fixes the boundary: matching
    ``variances`` to the selection needs the descriptor schema, and the schema
    comes from ``sources``, so it too is decided here with nothing to score.
    """
    with pytest.raises(RuntimeError, match="Unknown metric 'bogus'"):
        _native.validate_descriptor_options(
            _comparisons.descriptor_options({'metric': "bogus"}))
    with pytest.raises(RuntimeError, match="columns are selected"):
        _native.validate_descriptor_options(
            _comparisons.descriptor_options({'variances': [1.0, 2.0]}))
    _native.validate_descriptor_options(_comparisons.descriptor_options({}))


def test_pdist_and_cdist_agree_on_an_argument_no_input_can_rescue():
    """One bad argument, one answer, whichever entry point is used.

    ``cdist`` used to raise its own arrived-empty and cutoff refusals ahead of
    the argument check, so these three calls gave three different messages.
    The cutoff one was the worst: it named a remedy -- drop the cutoff, keep
    ``similarity=True`` -- that cannot make the call valid.
    """
    clean = _mols()
    untyped = _mols(["[Na+]", "[Fe]", "[He]"])

    with pytest.raises(ValueError, match="no similarity form"):
        oecluster.pdist([], "descriptor", similarity=True)
    with pytest.raises(ValueError, match="no similarity form"):
        oecluster.cdist([], clean, "descriptor", similarity=True)
    with pytest.raises(ValueError, match="no similarity form"):
        oecluster.cdist(untyped, clean, "descriptor", similarity=True,
                        cutoff=0.5)


def test_an_emptied_input_is_still_refused_when_the_arguments_are_valid():
    """Validating first must not cost the filter its own refusal.

    With nothing wrong in the arguments there is no authoritative error to
    outrank, so the emptied-list message is the right one and must survive.
    Over-refusing is the same size of defect as under-refusing, so two
    well-formed overrides are checked too. Between them they satisfy all five
    of the hoisted rules -- ascending columns, matching length, finite and
    positive entries, on both the ``variances`` and the ``inverse_covariance``
    side -- and must still lose to the emptied input.
    """
    untyped = _mols(["[Na+]", "[Fe]", "[He]"])
    stats = oecluster.descriptor_statistics(_mols(), inverse_covariance=True)
    good_variances = {'columns': list(stats['columns']),
                      'variances': list(stats['variance'])}
    good_inverse = {'metric': "mahalanobis",
                    'columns': list(stats['columns']),
                    'inverse_covariance': stats['inverse_covariance']}

    for kwargs in ({}, good_variances, good_inverse):
        with pytest.raises(ValueError, match="left 0 items"):
            oecluster.pdist(untyped, "descriptor", **kwargs)
        with pytest.raises(ValueError, match="left 0 item"):
            oecluster.cdist(untyped, _mols(), "descriptor", **kwargs)


def test_the_filter_refuses_the_options_it_is_handed():
    """The documented pre-filtering route validates what it is given.

    The header points a C++ caller at ``descriptor_excluded_indices`` to filter
    before constructing, so the options it gets are the ones the constructor
    sees next. It used to return a mask for options that could never build --
    a bogus metric, an unknown policy, both overrides at once -- and the caller
    learned nothing until the construction it was preparing for.
    """
    clean = _mols()
    for kwargs, expected in (({'metric': "bogus"}, "Unknown metric 'bogus'"),
                             ({'missing': "drop"},
                              "Unknown missing-value policy"),
                             ({'variances': [1.0],
                               'inverse_covariance': [1.0]},
                              "mutually exclusive")):
        with pytest.raises(RuntimeError, match=expected):
            _native.descriptor_excluded_indices(
                clean, _comparisons.descriptor_options(kwargs))

    # Still a filter, not only a gate: valid options get their mask.
    assert list(_native.descriptor_excluded_indices(
        clean, _native.DescriptorOptions())) == []


def test_a_selection_with_no_spread_is_refused_by_the_constructor():
    """A fitted metric needs variance, and only the molecules supply it.

    This is the third of the constructor's input-reading rules -- alongside
    the minimum input size and the complete-case row check -- and the one the
    header's prose once omitted. (The constructor keeps name checks as well,
    which read no input.) It cannot move upstream: the options here are
    perfectly valid, and identical inputs are what make them unusable.
    """
    # First, the placement pin: the validator must accept these options, so a
    # build that moved the rule upstream fails here rather than below.
    _native.validate_descriptor_options(
        _comparisons.descriptor_options({'columns': ["MolecularWeight"]}))

    identical = _mols(["CCO", "CCO"])
    with pytest.raises(RuntimeError,
                       match="Every selected descriptor column has zero "
                             "variance"):
        oecluster.pdist(identical, "descriptor", columns=["MolecularWeight"])

    # The options alone are never enough to know: the same request over
    # molecules that differ is fine, and an unfitted metric never asks. Each
    # control changes exactly one thing against the refused call above.
    assert oecluster.pdist(_mols(["CCO", "CCCC"]), "descriptor",
                           columns=["MolecularWeight"]).num_samples == 2
    assert oecluster.pdist(identical, "descriptor", metric="euclidean",
                           columns=["MolecularWeight"]).num_samples == 2


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


def test_validate_request_does_not_consume_kwargs():
    """The validator runs first, so consuming here would starve the builder.

    The failure is silent, not loud: ``build_comparison`` complains only about
    keys left *over*, so a popped one is hidden rather than reported. Popping
    ``missing`` from ``pdist(_mols(["O"] + SMILES[:4]), "descriptor",
    missing="propagate")`` restores the default complete-case policy and
    returns four molecules instead of five -- water dropped under a policy the
    caller explicitly declined.
    """
    kwargs = {'metric': "euclidean", 'missing': "propagate", 'p': 3.0}
    _comparisons.validate_request("descriptor", False, kwargs)
    assert kwargs == {'metric': "euclidean", 'missing': "propagate", 'p': 3.0}


def test_the_factory_class_builds_a_usable_comparison():
    comparison = oecluster.DescriptorComparison(_mols(), metric="euclidean")
    dist = oecluster.pdist(_mols(), comparison)
    assert dist.num_samples == 8
    assert dist.metric_capabilities['triangle'] is True


def test_the_factory_class_refuses_what_pdist_would_have_filtered():
    """The factory filters nothing, so the policy alone decides the outcome.

    ``pdist(mols, "descriptor")`` drops the water row and scores the other
    eight. The factory is handed all nine: under the default complete-case
    policy that is a refusal, and under "propagate" it is a nine-item
    comparison. Neither outcome is eight, which is what filtering would give.
    """
    mols = _mols(["O"] + SMILES)
    assert oecluster.pdist(mols, "descriptor").num_samples == 8

    with pytest.raises(RuntimeError, match="absent or non-finite"):
        oecluster.DescriptorComparison(mols)

    comparison = oecluster.DescriptorComparison(mols, missing="propagate")
    assert comparison.Size() == 9


def test_one_molecule_is_refused_only_by_a_metric_that_fits_variances():
    """The input-size floor belongs to the metric, not to the comparison.

    "standardized_euclidean" fits its scale from the input, so a single
    molecule gives it nothing to fit. "euclidean" fits nothing and accepts the
    same input, which is why the refusal cannot be stated unconditionally.
    """
    with pytest.raises(RuntimeError, match="at least two molecules"):
        oecluster.DescriptorComparison(_mols(["CCO"]))

    assert oecluster.DescriptorComparison(
        _mols(["CCO"]), metric="euclidean").Size() == 1
