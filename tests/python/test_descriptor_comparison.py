import re
import warnings

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
    """One NaN value costs the column, under its own label.

    OpenEye reports water's FractionCsp3 as present and NaN, so the column is
    present in all nine rows and its variance is not finite. That is a
    different cause from the genuinely constant column
    ``test_a_constant_column_is_dropped_as_zero_variance`` covers -- this one
    varies perfectly well over the eight organics on its own -- and the drop
    report now says so.
    """
    clean = oecluster.descriptor_statistics(_mols())
    clean_variance = dict(zip(clean['columns'], clean['variance']))
    stats = oecluster.descriptor_statistics(_mols(["O"] + SMILES))

    assert clean['dropped'] == []
    assert clean_variance["FractionCsp3"] > 0.0

    assert stats['dropped'] == [("FractionCsp3", "non-finite-variance")]
    assert "FractionCsp3" not in stats['columns']
    assert all(v > 0.0 for v in stats['variance'])


def test_a_column_present_in_one_row_is_not_called_constant():
    """It had too few values to measure, which is a different problem.

    XLogP is computable for CCO and for none of the ions, so its variance
    was never computed -- OEFP leaves it at NaN. Reporting that as
    "zero-variance" told the caller the column carries no signal, when the
    remedy is more molecules where the descriptor is defined.
    """
    mols = _mols(["[Na+]", "[K+]", "[Li+]", "CCO"])
    stats = oecluster.descriptor_statistics(mols)
    assert ("XLogP", "insufficient-data") in stats['dropped']


def test_a_genuinely_constant_column_still_says_zero_variance():
    """Regression shield: the unchanged case keeps its unchanged reason."""
    stats = oecluster.descriptor_statistics(_mols(["CCO", "CCC"]))
    reasons = dict(stats['dropped'])
    assert reasons["HeavyAtomCount"] == "zero-variance"


def test_the_comparison_surface_names_the_same_causes_as_statistics():
    # The two surfaces are separate code paths and must not disagree. They
    # were independent copies of the rule, which is how they could.
    mols = _mols(["[Na+]", "[K+]", "[Li+]", "CCO"])
    comparison = oecluster.DescriptorComparison(
        mols, missing="propagate", metric="standardized_euclidean")
    dropped = dict(zip(list(comparison.DroppedColumns()),
                       list(comparison.DroppedReasons())))
    assert dropped["XLogP"] == "insufficient-data"

    constant = oecluster.DescriptorComparison(
        _mols(["CCO", "CCC"]), metric="standardized_euclidean")
    reasons = dict(zip(list(constant.DroppedColumns()),
                       list(constant.DroppedReasons())))
    assert reasons["HeavyAtomCount"] == "zero-variance"


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
    variance filter here, so the sodium row reaches the per-column
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

    assert stats['dropped'] == [("FractionCsp3", "non-finite-variance")]
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


def test_a_finite_propagate_matrix_is_stamped_complete_and_clusters():
    """The stamp describes the data, not the policy that permitted a NaN.

    ``missing='propagate'`` used to stamp ``nan_present`` from the policy name
    alone, so a matrix with no non-finite entry anywhere was refused with
    "contains non-finite entries" -- a claim that was not true of it.
    """
    dist = oecluster.pdist(_mols(), "descriptor", missing="propagate")
    assert np.isfinite(dist.condensed).all()
    assert dist.data_integrity == "complete"
    oecluster.butina(dist, threshold=0.5)


def test_a_propagate_matrix_that_really_holds_a_nan_is_still_refused():
    """The half that proves the fix did not simply disable the refusal."""
    dist = oecluster.pdist(_missingness_mols(), "descriptor",
                           missing="propagate")
    assert not np.isfinite(dist.condensed).all()
    assert dist.data_integrity == "nan_present"
    with pytest.raises(ValueError, match="non-finite entries"):
        oecluster.butina(dist, threshold=0.5)


def test_a_sparse_matrix_keeps_its_stamp_because_the_cutoff_hid_the_infinity():
    """SparseStorage discards ``+inf``, so the data is not a witness.

    ``Set`` is ``if (value > cutoff_) return;`` (src/StorageBackend.cpp:436)
    and ``+inf`` exceeds every finite cutoff, so an infinite distance is
    dropped rather than stored. The matrix then scans clean while the stamp is
    truthful, and downgrading it would turn a true stamp into a false one.
    """
    # Pin the premise: without the cutoff the same call really does overflow.
    # Were a fixture change to stop producing ``+inf``, the assertion below
    # would still pass while testing nothing at all.
    uncut = oecluster.pdist(_mols(), "descriptor", metric="minkowski",
                            p=200.0, missing="propagate")
    assert not np.isfinite(uncut.condensed).all()
    assert uncut.data_integrity == "nan_present"

    dist = oecluster.pdist(_mols(), "descriptor", metric="minkowski",
                           p=200.0, missing="propagate", cutoff=1.0)
    assert dist.data_integrity == "nan_present"


def test_an_mmap_matrix_with_a_cutoff_is_still_downgraded(tmp_path):
    """mmap storage ignores the cutoff, so it stays a faithful record.

    This is what distinguishes the correct condition from ``cutoff == 0.0``:
    PDistOptions.cutoff is not read by pdist (include/oecluster/PDist.h:22-26).
    """
    dist = oecluster.pdist(_mols(), "descriptor", missing="propagate",
                           cutoff=0.5, output=str(tmp_path / "m.dat"))
    assert dist.data_integrity == "complete"


def test_a_finite_propagate_cross_matrix_is_stamped_complete():
    cross = oecluster.cdist(_mols()[:4], _mols()[4:], "descriptor",
                            missing="propagate")
    assert np.isfinite(cross.matrix).all()
    assert cross.data_integrity == "complete"


def test_a_cross_matrix_keeps_its_stamp_because_the_cutoff_hid_the_infinity():
    """cdist's counterpart to the sparse shield, on its own condition.

    The cutoff path rewrites a value above the threshold rather than storing
    it (src/comparisons/DescriptorComparison.cpp:520), so the returned matrix
    scans clean while the stamp is truthful. Without the ``cutoff > 0.0`` test
    in cdist, this true ``nan_present`` would silently be downgraded.
    """
    a, b = _mols()[:4], _mols()[4:]
    # Pin the premise, as the sparse test does: uncut, this call overflows.
    uncut = oecluster.cdist(a, b, "descriptor", metric="minkowski", p=200.0,
                            missing="propagate")
    assert not np.isfinite(uncut.matrix).all()
    assert uncut.data_integrity == "nan_present"

    cross = oecluster.cdist(a, b, "descriptor", metric="minkowski", p=200.0,
                            missing="propagate", cutoff=1.0)
    assert np.isfinite(cross.matrix).all()
    assert cross.data_integrity == "nan_present"


def test_a_cross_cutoff_is_read_from_the_value_native_filtered_with():
    """The stamp cannot rest on a second reading of the caller's cutoff.

    cdist's counterpart to test_pdist_cutoff_positivity_is_decided_once. The
    guard and the stamping test used to compare the caller's object
    independently, so a float subclass answering positive to the first and
    zero to the second cleared both: native received 1.0 and zeroed every
    overflowed pair, and the all-zero result was then stamped ``complete``.
    Converting once and reading the double back ties the decision to the
    value native filtered with, so the caller's ``__gt__`` is never consulted.
    """
    class MirrorCutoff(float):
        """Compares as positive once, then as zero."""

        def __init__(self, _value):
            self._comparisons = 0

        def __gt__(self, other):
            self._comparisons += 1
            return self._comparisons == 1

    a, b = _mols()[:4], _mols()[4:]
    cutoff = MirrorCutoff(1.0)
    cross = oecluster.cdist(a, b, "descriptor", metric="minkowski", p=200.0,
                            missing="propagate", cutoff=cutoff)

    assert cutoff._comparisons == 0
    # Native filtered with 1.0, so the overflow is gone from the output and
    # only the stamp still records it.
    assert np.isfinite(cross.matrix).all()
    assert cross.data_integrity == "nan_present"


def test_an_infinite_cross_cutoff_downgrades_because_it_filters_nothing():
    """Only a positive *finite* cutoff can hide a value.

    Both filters drop ``distance > cutoff`` (src/CDist.cpp:74,
    src/comparisons/DescriptorComparison.cpp:520-526) and no double exceeds
    ``+inf``, so this matrix is a faithful record and the data decides. This
    is not the residual spec section 5.3 accepts: that one is pdist's, where
    an infinite cutoff selects sparse storage instead of feeding a filter.
    """
    cross = oecluster.cdist(_mols()[:4], _mols()[4:], "descriptor",
                            missing="propagate", cutoff=float("inf"))
    assert np.isfinite(cross.matrix).all()
    assert cross.data_integrity == "complete"


def test_a_subset_scored_cross_matrix_keeps_its_own_stamp():
    """The downgrade acts on ``nan_present`` and on nothing else.

    ``missing='ignore'`` scores each pair over the columns both items have,
    which is a claim about coverage rather than about finiteness, so a finite
    matrix is no reason to rewrite it.
    """
    cross = oecluster.cdist(_mols()[:4], _mols()[4:], "descriptor",
                            metric="euclidean", missing="ignore")
    assert np.isfinite(cross.matrix).all()
    assert cross.data_integrity == "subset_scored"


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


# Option mistakes no molecule set can rescue, paired with the message each has
# to produce. Not the whole set -- the C++ table in
# ``tests/cpp/test_descriptor_comparison.cpp`` is wider, and unknown
# source/column/group names are unrescuable too but are reported downstream on
# a request with no override. These are the ones pinned through the Python
# surface, chosen to cover both halves of the boundary.
#
# Every case below the first two is applied only once the descriptor schema is
# resolved, and that schema comes from ``sources`` alone -- which is why the
# boundary is "needs no molecules" rather than "needs no schema", and what puts
# them ahead of a filter that can empty the item list. Most read the schema; the
# ``inverse_covariance`` finiteness case does not, and is applied there only to
# stay beside the length rule it follows. Where the fault is in the override's
# values rather than its length, ``columns`` is named explicitly, so the case
# does not depend on how many columns the default selection happens to hold.
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
    # Indefinite with a positive diagonal, so a diagonal-only screen would wave
    # it through; separating it takes a verdict on the matrix as a whole.
    ({'metric': "mahalanobis", 'columns': ["MolecularWeight", "XLogP"],
      'inverse_covariance': [1.0, 2.0, 2.0, 1.0]},
     "positive semidefinite"),
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


def _mahalanobis(matrix, columns):
    """Score the fixture molecules under an explicit inverse covariance."""
    return oecluster.pdist(_mols(), "descriptor", metric="mahalanobis",
                           columns=columns,
                           inverse_covariance=matrix.ravel().tolist()).condensed


def _first_columns(count):
    stats = oecluster.descriptor_statistics(_mols())
    assert len(stats['columns']) >= count
    return list(stats['columns'][:count])


# The three matrices that are reported rather than refused, one per failure
# mode a refusal rule was tried against and lost to.
ASYMMETRIC = [
    # The item's own example: a transposed pair, asymmetric by order 1.
    (np.array([[2.0, 1.0], [0.0, 2.0]]), r"\[0\]\[1\] is 1\.0 but \[1\]\[0\] is 0\.0"),
    # A triangular fill, weakly correlated enough that every off-diagonal sits
    # far below the diagonal scaling it. The fill signal is what reaches the
    # rest of that family: below the band an entry is caught only by it, and
    # its own band narrows with the dimension, so a 1000-column fill is still
    # reported where this one would also trip the magnitude signal.
    (np.array([[1.0, 9e-4], [0.0, 1.0]]), r"\[0\]\[1\] is 0\.0009 but \[1\]\[0\] is 0\.0"),
    # And the mirror image, which an upper-triangle-only test would miss.
    (np.array([[1.0, 0.0], [9e-4, 1.0]]), r"\[0\]\[1\] is 0\.0 but \[1\]\[0\] is 0\.0009"),
]


def test_a_refusal_is_not_pre_empted_by_the_asymmetry_warning():
    """An advisory must not stand in front of a real fault.

    ``[[1, 3], [0, 1]]`` is asymmetric *and* its symmetric part is
    indefinite, so it earns both the warning and the semidefinite refusal.
    The refusal is the one the caller has to act on, and under ``-W error``
    an advisory emitted first would replace it entirely. The constructor
    substitutes the matrix up front but holds the warning until C++ has
    accepted the options, which is the ordering the named paths get for
    free by validating first.
    """
    bad = np.array([[1.0, 3.0], [0.0, 1.0]])
    columns = _first_columns(2)

    with (warnings.catch_warnings(),
          pytest.raises(RuntimeError, match="positive semidefinite")):
        warnings.simplefilter("error")
        oecluster.DescriptorComparison(
            _mols(), columns=columns, metric="mahalanobis",
            inverse_covariance=bad.ravel().tolist())

    # And with warnings merely shown, the refusal still wins.
    with (warnings.catch_warnings(),
          pytest.raises(RuntimeError, match="positive semidefinite")):
        warnings.simplefilter("always")
        oecluster.DescriptorComparison(
            _mols(), columns=columns, metric="mahalanobis",
            inverse_covariance=bad.ravel().tolist())


@pytest.mark.parametrize(("matrix", "pair"), ASYMMETRIC,
                         ids=["transposed-pair", "triangular-fill", "mirrored-fill"])
def test_an_asymmetric_inverse_covariance_warns_and_is_symmetrized(matrix, pair):
    """Reported, not refused, because no tolerance can tell the two apart.

    Five successive refusal rules each fell to a reproduced counterexample.
    What ended the search was that a sub-tolerance perturbation is unbounded
    in relative terms -- ``v.T @ A @ v`` can be arbitrarily small along a
    weakly constrained direction, so however tight a refusal band is set,
    a matrix under it can still be scored arbitrarily far from what the
    caller supplied. The defect was never the asymmetry, it was the silence,
    so the substitution is announced instead of guarded against.
    """
    columns = _first_columns(2)
    with pytest.warns(UserWarning, match="inverse_covariance is not symmetric") as caught:
        warned = _mahalanobis(matrix, columns)
    assert len(caught) == 1
    message = str(caught[0].message)
    # The offending pair and both of its values, so a caller holding a large
    # matrix knows which entry to check.
    assert re.search(pair, message), message
    assert "(A + A.T) / 2" in message

    # What was scored is the symmetrization, and nothing else: the matrix the
    # caller would have had to build themselves.
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        explicit = _mahalanobis((matrix + matrix.T) / 2.0, columns)
    np.testing.assert_array_equal(warned, explicit)


def test_the_case_no_refusal_band_could_have_caught_still_warns():
    """The counterexample that ended the search for a refusal rule.

    One triangle is filled: the caller wrote 3.33e-4 at ``[0][1]`` and left
    ``[1][0]`` at zero, so they meant the mirrored matrix. What is scored is
    ``(A + A.T) / 2``, which halves that entry instead. Along the direction
    ``[0.499979, 0.499979, -0.707136]``, near the matrix's weakly constrained
    subspace, the intended distance is 9.12e-3 and the scored one 2.90e-6 --
    a change of 99.97% from a perturbation of 1.7e-4. ``v.T @ A @ v`` can be
    arbitrarily small there, so the relative effect of a perturbation is
    unbounded however small the perturbation is, and no refusal threshold
    could have been set low enough. That is why this is a warning.

    It is also what moved the band from 1e-3 to 1e-6: it scores 3.33e-4,
    which the band the refusal era required let through in silence -- the
    original defect, reproduced by the very case that caused the redesign.
    """
    a = 0.70716564537
    matrix = np.array([[1.0, 3.33e-4, a],
                       [0.0, 1.0, a],
                       [a, a, 1.0]])
    # Assert the consequence rather than assuming it, so the fixture cannot
    # quietly stop being the case this test is named for.
    intended = np.where(matrix != 0.0, matrix, matrix.T)
    symmetrized = (matrix + matrix.T) / 2.0
    v = np.array([0.499979, 0.499979, -0.707136])
    np.testing.assert_allclose(np.sqrt(v @ intended @ v), 9.12e-3, rtol=1e-2)
    np.testing.assert_allclose(np.sqrt(v @ symmetrized @ v), 2.90e-6, rtol=1e-2)

    columns = _first_columns(3)
    with pytest.warns(UserWarning, match="inverse_covariance is not symmetric"):
        oecluster.DescriptorComparison(
            _mols(), columns=columns, metric="mahalanobis",
            inverse_covariance=matrix.ravel().tolist())


def test_the_warning_reaches_the_prebuilt_constructor_too():
    """Both descriptor surfaces warn, because both share one helper.

    ``p`` diverged between these two before, which is why the rule lives in
    one place. The constructor builds its options before the substitution, so
    this also pins that it rebuilds them -- a comparison holding the original
    matrix would warn and then score it anyway.
    """
    matrix = np.array([[2.0, 1.0], [0.0, 2.0]])
    columns = _first_columns(2)
    with pytest.warns(UserWarning, match="inverse_covariance is not symmetric"):
        warned = oecluster.DescriptorComparison(
            _mols(), columns=columns, metric="mahalanobis",
            inverse_covariance=matrix.ravel().tolist())

    with warnings.catch_warnings():
        warnings.simplefilter("error")
        explicit = oecluster.DescriptorComparison(
            _mols(), columns=columns, metric="mahalanobis",
            inverse_covariance=((matrix + matrix.T) / 2.0).ravel().tolist())
    np.testing.assert_array_equal(oecluster.pdist(_mols(), warned).condensed,
                                  oecluster.pdist(_mols(), explicit).condensed)


# Every surface that takes descriptor options directly. The asymmetry
# advisory has to behave identically on all three, which is the property that
# broke twice: once on the constructor against native validation, and once on
# the named paths against the emptiness refusal.
ASYMMETRY_SURFACES = [
    ("pdist", lambda kwargs: oecluster.pdist(_mols(), "descriptor", **kwargs)),
    ("cdist", lambda kwargs: oecluster.cdist(_mols(), _mols(), "descriptor", **kwargs)),
    ("constructor", lambda kwargs: oecluster.DescriptorComparison(_mols(), **kwargs)),
]


@pytest.mark.parametrize(("surface", "call"), ASYMMETRY_SURFACES,
                         ids=[case[0] for case in ASYMMETRY_SURFACES])
def test_an_asymmetry_advisory_never_pre_empts_an_emptiness_refusal(surface, call):
    """An advisory in front of a refusal is the same defect, a third time.

    ``groups=[]`` is refused because C++ cannot tell it from an omitted
    option, and that refusal is last in the argument band -- after the native
    verdicts the warning used to wait for. Warning from there put the
    advisory first, and under ``-W error`` the ``UserWarning`` replaced the
    ``ValueError``, hiding the fault the caller actually has to fix. Measured
    before the move on ``pdist`` and ``cdist``; the constructor already
    raised first.
    """
    kwargs = {'metric': "mahalanobis", 'columns': _first_columns(2),
              'inverse_covariance': [2.0, 1.0, 0.0, 2.0], 'groups': []}

    # Under the default filters: the refusal, and not a warning beside it.
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        with pytest.raises(ValueError, match="groups= was given as an empty"):
            call(kwargs)
    assert [str(w.message) for w in caught] == []

    # And under -W error, where an advisory raised first would be the only
    # thing the caller ever sees.
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        with pytest.raises(ValueError, match="groups= was given as an empty"):
            call(kwargs)


@pytest.mark.parametrize(("surface", "call"), ASYMMETRY_SURFACES,
                         ids=[case[0] for case in ASYMMETRY_SURFACES])
def test_a_non_string_metric_is_refused_rather_than_read_by_the_advisory(surface, call):
    """The helper must not read a metric only C++ can rule on.

    It gated on ``metric`` to mirror the native rule that the override
    belongs to ``mahalanobis`` alone, and called ``.lower()`` to do it. On
    ``metric=5`` that raised ``AttributeError: 'int' object has no attribute
    'lower'`` from the constructor, which validates after substituting, while
    ``pdist`` gave ``RuntimeError`` for the same call -- a regression against
    5.22.0 and a disagreement between two surfaces of one library. The helper
    now stands aside for a metric it cannot read.
    """
    kwargs = {'metric': 5, 'inverse_covariance': [2.0, 1.0, 0.0, 2.0],
              'columns': _first_columns(2)}
    for filt in ("always", "error"):
        with warnings.catch_warnings():
            warnings.simplefilter(filt)
            with pytest.raises(RuntimeError,
                               match="inverse_covariance applies only to"):
                call(kwargs)


def test_an_ill_conditioned_numpy_inverse_covariance_is_silent():
    """A real inv(cov(X)) is never exactly symmetric, and must not be reported.

    At condition number 5.2e9 this matrix has entries -24294.14521055912 and
    -24294.14521062444, a per-element relative asymmetry of 2.7e-12, which a
    per-element tolerance tight enough to catch a real mistake would flag.
    The detector judges each pair against sqrt(|Aii|) * sqrt(|Ajj|) instead,
    under which this matrix scores 8.5e-17.
    """
    rng = np.random.default_rng(3238)
    q = np.linalg.qr(rng.normal(size=(5, 5)))[0]
    x = rng.normal(size=(15, 5)) @ q @ np.diag(np.geomspace(1, 1e-5, 5))
    a = np.linalg.inv(np.cov(x, rowvar=False))
    assert not np.array_equal(a, a.T)
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        oecluster.DescriptorComparison(
            _mols(), columns=_first_columns(5), metric="mahalanobis",
            inverse_covariance=a.ravel().tolist())


def test_an_inverse_whose_roundoff_lands_on_exact_zero_is_silent():
    """Round-off does produce a one-sided exact zero, so the fill signal bands it.

    This covariance is symmetric, positive definite and well conditioned
    (13.9), yet ``np.linalg.inv`` returns -8.88e-16 against an exact 0.0 at
    one mirrored pair. An unbanded "exactly one side zero" rule refused it --
    the first version of the triangular-fill check did exactly that.
    """
    cov = np.array([[0.25, 0.125, 0.125],
                    [0.125, 0.25, 0.125],
                    [0.125, 0.125, 0.125]])
    a = np.linalg.inv(cov)
    assert (a[0, 1] == 0.0) != (a[1, 0] == 0.0)
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        oecluster.DescriptorComparison(
            _mols(), columns=_first_columns(3), metric="mahalanobis",
            inverse_covariance=a.ravel().tolist())


def test_the_librarys_own_inverse_covariance_round_trips_silently():
    """A library must not warn about what it just handed the caller.

    ``descriptor_statistics(inverse_covariance=True)`` is the documented
    source for this override, so its own output passed straight back in is
    the one input that can never be a caller's mistake.
    """
    stats = oecluster.descriptor_statistics(_mols(), inverse_covariance=True)
    a = stats['inverse_covariance']
    assert a is not None
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        oecluster.pdist(_mols(), "descriptor", metric="mahalanobis",
                        columns=list(stats['columns']),
                        inverse_covariance=a)


def test_a_dense_correlated_ill_conditioned_inverse_is_silent():
    """The upper end of the regime the band promises to stay quiet over.

    An earlier stress test scaled singular values along a fixed orthonormal
    basis, which does not mix eigenvectors the way real correlated data does,
    and so understated the round-off in an inverse by several orders. With
    dense mixing a legitimate inv(cov(X)) scores about cond * 1e-17, which an
    earlier 1e-8 band would have reported. Measured over 30 draws per decade,
    1e-6 is the last band that stays silent on real data: the worst
    legitimate asymmetry is 2.7e-7 at condition number 1e10 and 1.8e-6 at
    1e11, where warnings begin.
    """
    rng = np.random.default_rng(32)
    q = np.linalg.qr(rng.normal(size=(5, 5)))[0]
    sv = np.geomspace(1, 1e-5, 5) * rng.uniform(0.5, 2.0, 5)
    x = rng.normal(size=(30, 5)) @ (q * sv) @ q.T
    cov = np.cov(x, rowvar=False)
    a = np.linalg.inv(cov)
    assert np.linalg.cond(cov) > 1e10
    assert not np.array_equal(a, a.T)
    # Bracket the fixture rather than pinning a magnitude. It has to be
    # asymmetric by far more than pure round-off, or the silence proves
    # nothing, and below the band, or the test is asserting the wrong
    # outcome. The exact value is not ours to fix: it depends on the
    # platform's LAPACK, and an earlier `> 1e-8` passed on macOS at 4.8e-8
    # while failing on Linux at 8.1e-9 and on Windows -- the same inverse,
    # a different library.
    scale = np.sqrt(np.abs(np.outer(np.diag(a), np.diag(a))))
    skew = np.abs(a - a.T) / np.where(scale > 0, scale, 1.0)
    np.fill_diagonal(skew, 0.0)
    assert 1e-12 < skew.max() < 1e-6
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        oecluster.DescriptorComparison(
            _mols(), columns=_first_columns(5), metric="mahalanobis",
            inverse_covariance=a.ravel().tolist())


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
    well-formed overrides are checked too. Between them they satisfy every
    override rule that waits on the resolved selection -- ascending columns,
    matching length, finite and positive entries, on both the ``variances``
    and the ``inverse_covariance`` side, and the semidefinite verdict on the
    latter -- and must still lose to the emptied input.
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

    This is one of the rules the constructor keeps because it has to read the
    input to decide -- the null-pointer scan and the complete-case row check
    are two others -- and the one the header's prose once omitted. It cannot
    move upstream: the options here are perfectly valid, and identical inputs
    are what make them unusable.
    """
    # First, the placement pin: the validator must accept these options, so a
    # build that moved the rule upstream fails here rather than below.
    _native.validate_descriptor_options(
        _comparisons.descriptor_options({'columns': ["MolecularWeight"]}))

    identical = _mols(["CCO", "CCO"])
    with pytest.raises(RuntimeError,
                       match=re.escape(
                           "Every selected descriptor column was dropped over "
                           "these molecules (zero-variance)")):
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

    The metric is "minkowski" because ``p`` now requires it; the subject here
    is keyword consumption, not applicability.
    """
    kwargs = {'metric': "minkowski", 'missing': "propagate", 'p': 3.0}
    _comparisons.validate_request("descriptor", False, kwargs)
    assert kwargs == {'metric': "minkowski", 'missing': "propagate", 'p': 3.0}


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


SEQUENCE_OPTIONS = ("sources", "columns", "groups", "variances",
                    "inverse_covariance")


@pytest.mark.parametrize("option", SEQUENCE_OPTIONS)
def test_an_empty_sequence_option_is_not_the_same_as_omitting_it(option):
    """Only ``None`` is unspecified, on every entry point that builds options.

    Emptiness is the native layer's own encoding of "no override", so an empty
    sequence used to reach C++ indistinguishable from an omitted option and
    silently produce the default result: over ``_mols()`` every one of these
    five scored ``3.4842533260889486`` and dropped no column, exactly as
    passing nothing does.

    The keyword path and the prebuilt path are separate code in the package --
    ``pdist`` and ``cdist`` are refused by ``_validate_descriptor``, while
    ``DescriptorComparison`` converts and refuses in its own ``__new__`` -- so
    both are checked. The message has to name the option, since a caller who
    passed several has no other way to tell which one was empty.
    """
    match = rf"{option}= was given as an empty sequence"

    with pytest.raises(ValueError, match=match):
        oecluster.pdist(_mols(), "descriptor", **{option: []})
    with pytest.raises(ValueError, match=match):
        oecluster.cdist(_mols(), _mols(), "descriptor", **{option: []})
    with pytest.raises(ValueError, match=match):
        oecluster.DescriptorComparison(_mols(), **{option: []})


def test_an_empty_inverse_covariance_array_is_refused_once_flattened():
    """The matrix form is raveled first, so a ``(0, 0)`` array is empty too.

    ``inverse_covariance`` is the one option accepted as a 2-D array, and the
    guard has to see what C++ will see rather than the shape handed in.
    """
    with pytest.raises(ValueError,
                       match="inverse_covariance= was given as an empty"):
        oecluster.pdist(_mols(), "descriptor", metric="mahalanobis",
                        inverse_covariance=np.zeros((0, 0)))


def test_an_empty_sequence_loses_to_the_two_rules_settled_without_a_value():
    """``similarity=True`` and an unknown keyword open the argument-level band.

    The guard closes that band -- the rules that need no molecules -- and these
    two open it: ``_validate_descriptor`` decides both before it converts
    anything, neither depending on a value being readable at all. Between them
    sit the C++ verdicts on the option values, which outrank the guard as well.

    What the guard does precede is everything past the band -- but which
    reports lie past it depends on the request. The source, column and group
    names are among them only while no non-empty override is passed, because
    that is when nothing in the band resolves a schema to check a name
    against. Pass one and ``validate_descriptor_options`` resolves the schema
    and rules on all three from inside the band, ahead of the guard. Both
    halves are asserted below.
    """
    with pytest.raises(ValueError, match="no similarity form"):
        oecluster.pdist(_mols(), "descriptor", similarity=True, variances=[])
    with pytest.raises(
            TypeError,
            match=r"Unknown kwargs for descriptor comparison: \['bogus'\]"):
        oecluster.pdist(_mols(), "descriptor", bogus=1, variances=[])

    # A non-empty override is what makes the validator resolve a schema, which
    # moves all three name verdicts into the band and ahead of the guard.
    for names, match in (
            ({"sources": ["nope"], "columns": []},
             "Unknown descriptor source: nope"),
            ({"columns": ["nope"], "groups": []},
             "Unknown descriptor column: nope"),
            ({"groups": ["nope"], "columns": []},
             "Unknown or empty descriptor group: nope")):
        with pytest.raises(RuntimeError, match=match):
            oecluster.pdist(_mols(), "descriptor", variances=[1.0], **names)

    # The same request without the override: no schema is resolved, so the bad
    # source falls past the band and the guard reports first.
    with pytest.raises(ValueError, match="columns= was given as an empty"):
        oecluster.pdist(_mols(), "descriptor", sources=["nope"], columns=[])


@pytest.mark.parametrize("option", SEQUENCE_OPTIONS)
def test_an_empty_sequence_loses_to_a_value_verdict_from_cpp(option):
    """``validate_descriptor_options`` rules on the values before the guard.

    An unknown metric and an unknown missing policy are verdicts on values the
    caller has to fix whatever they do about the empty sequence, so they
    outrank it. Both entry points reach the validator first: ``pdist`` through
    ``_validate_descriptor``, and the prebuilt constructor through a validation
    call of its own, because the C++ constructor that would otherwise give the
    verdict is never reached on a call about to be refused.
    """
    for fault, match in (({"metric": "bogus"}, "Unknown metric 'bogus'"),
                         ({"missing": "bogus"},
                          "Unknown missing-value policy: bogus")):
        with pytest.raises(RuntimeError, match=match):
            oecluster.pdist(_mols(), "descriptor", **{option: []}, **fault)
        with pytest.raises(RuntimeError, match=match):
            oecluster.DescriptorComparison(_mols(), **{option: []}, **fault)


def test_cdist_orders_the_guard_through_the_same_validator_as_pdist():
    """One case covers the ordering, because that much is not ``cdist``'s code.

    ``pdist`` and ``cdist`` both reach the descriptor rules through
    ``validate_request``, which dispatches to the single registered validator;
    the identity assertion is what makes one case representative. Where the
    band as a whole lands among ``cdist``'s own guards *is* its own code, and
    is covered by the test below.
    """
    assert (_comparisons._VALIDATORS["descriptor"]
            is _comparisons._validate_descriptor)

    with pytest.raises(RuntimeError, match="Unknown metric 'bogus'"):
        oecluster.cdist(_mols(), _mols(), "descriptor", metric="bogus",
                        variances=[])


def test_the_validator_holds_a_copy_of_the_guard_and_cdist_needs_it():
    """The refusal in ``_validate_descriptor`` is load-bearing, not a repeat.

    ``descriptor_options`` refuses inside the builder, and the descriptor
    normalizer reaches that same copy before the builder does, which answers
    almost every route to an empty option whether or not the validator repeats
    the check. ``cdist`` is the exception: it refuses an input set that
    arrived empty in between, so with only the builder's copy
    ``cdist([], mols, "descriptor", columns=[])`` reports the shape and names
    a remedy that would leave the option still empty. One case takes both
    sides, because one ``if not a or not b`` answers for both.

    ``pdist`` cannot stand in for it. It leaves an input that arrived empty
    for the comparison to answer -- refusing only when normalization is what
    emptied the list -- so nothing of its own sits between ``validate_request``
    and the builder's copy. Its case below is the contrast that shows why,
    not coverage: it holds with the validator's copy deleted.
    """
    match = "columns= was given as an empty sequence"

    with pytest.raises(ValueError, match=match):
        _comparisons.validate_request("descriptor", False, {"columns": []})

    for a, b in (([], _mols()), (_mols(), [])):
        with pytest.raises(ValueError, match=match):
            oecluster.cdist(a, b, "descriptor", columns=[])

    with pytest.raises(ValueError, match=match):
        oecluster.pdist([], "descriptor", columns=[])


@pytest.mark.parametrize("metric,option", [
    ("mahalanobis", "variances"),
    ("standardized_euclidean", "inverse_covariance"),
])
def test_an_empty_override_is_refused_rather_than_paired_with_its_metric(
        metric, option):
    """An empty override is no override, so no pairing rule fires against it.

    Supplied non-empty to the metric it does not belong to, each of these is
    refused by ``validate_descriptor_options``. Empty, neither sets
    ``has_override``, so the validator returns before it resolves a schema and
    the caller sees the emptiness refusal rather than a mismatch against a
    metric they never overrode.
    """
    match = rf"{option}= was given as an empty sequence"

    with pytest.raises(ValueError, match=match):
        oecluster.pdist(_mols(), "descriptor", metric=metric, **{option: []})
    with pytest.raises(ValueError, match=match):
        oecluster.DescriptorComparison(_mols(), metric=metric,
                                       **{option: []})


def test_an_empty_sequence_still_wins_against_a_name_resolved_downstream():
    """A source name is not an argument-level verdict, so the guard precedes it.

    ``validate_descriptor_options`` builds a schema only when an override is
    supplied, and so resolves ``sources`` on no other request; the unknown name
    is reported further down instead. Not for want of molecules -- an empty
    item list still draws the name verdict, ahead of the refusal it draws with
    a valid source -- but because nothing in the argument-level band resolves a
    schema to check the name against. The guard sits at the end of that band,
    which puts it ahead of the report -- and the same bad source with nothing
    empty alongside still reaches it.
    """
    with pytest.raises(ValueError, match="columns= was given as an empty"):
        oecluster.pdist(_mols(), "descriptor", sources=["nope"], columns=[])

    with pytest.raises(RuntimeError, match="Unknown descriptor source: nope"):
        oecluster.pdist(_mols(), "descriptor", sources=["nope"])

    with pytest.raises(RuntimeError, match="Unknown descriptor source: nope"):
        oecluster.pdist([], "descriptor", sources=["nope"])

    with pytest.raises(RuntimeError, match="fits variances from"):
        oecluster.pdist([], "descriptor")


def test_descriptor_statistics_refuses_an_empty_sequence_before_it_computes():
    """Its only name verdict comes out of the computing call, so the guard leads.

    ``descriptor_statistics`` has no argument-level validator to run first:
    ``validate_descriptor_options`` takes a ``DescriptorOptions`` and not the
    statistics options type, and ``Unknown descriptor source`` is raised by the
    call that reads the molecules. Ordering that verdict ahead of the guard
    would mean computing and discarding a full descriptor table on a call about
    to be refused, so the guard keeps the end of the argument-level band --
    which is the same place, and the same outcome, as ``pdist`` given a bad
    source and an empty option.
    """
    with pytest.raises(ValueError, match="columns= was given as an empty"):
        oecluster.descriptor_statistics(_mols(), sources=["nope"], columns=[])

    with pytest.raises(RuntimeError, match="Unknown descriptor source: nope"):
        oecluster.descriptor_statistics(_mols(), sources=["nope"])


def test_a_conversion_failure_outranks_the_emptiness_guard():
    """Every sequence option is converted before any of them is judged empty.

    An empty sequence is refused by an advisory rule this layer owns alone; a
    value that cannot be converted at all is a verdict about what was passed,
    and it outranks the advisory one. Converting the whole set first is what
    makes the ranking hold whichever option carries which fault, rather than
    only when the unconvertible one happens to come first in the signature.

    ``descriptor_statistics`` converted and refused one option at a time, so
    ``sources=[]`` stopped ``columns=1`` from ever reaching the conversion and
    the empty sequence answered for it -- disagreeing with the two comparison
    paths, which have converted the whole set since the guard was reordered.
    """
    for faults in ({'sources': [], 'columns': 1},
                   {'columns': [], 'groups': 1}):
        for call in (oecluster.descriptor_statistics,
                     oecluster.DescriptorComparison):
            with pytest.raises(TypeError, match="'int' object is not iterable"):
                call(_mols(), **faults)
        with pytest.raises(TypeError, match="'int' object is not iterable"):
            oecluster.pdist(_mols(), "descriptor", **faults)


@pytest.mark.parametrize("option", ("sources", "columns", "groups"))
def test_descriptor_statistics_refuses_an_empty_sequence_option(option):
    """It builds its own options object, so it needs the guard of its own.

    ``descriptor_statistics`` shares the vector conversion with
    ``descriptor_options`` but not the option building, and before the guard
    each of these three returned the same eleven default columns as passing
    nothing at all.
    """
    assert len(oecluster.descriptor_statistics(_mols())['columns']) == 11

    with pytest.raises(ValueError,
                       match=rf"{option}= was given as an empty sequence"):
        oecluster.descriptor_statistics(_mols(), **{option: []})


# The descriptor surfaces the Python layer exposes. The rule lives in one
# shared helper but is called from two places, and only a test over all three
# catches a call site that was missed. These are not every route to the
# descriptor options: the documented ``oecluster.oecluster`` bindings reach
# them underneath this layer, and
# ``test_the_raw_native_bindings_still_accept_an_inert_p`` pins that boundary.
_DESCRIPTOR_SURFACES = [
    pytest.param(lambda mols, **kw: oecluster.pdist(mols, "descriptor", **kw),
                 id="pdist"),
    pytest.param(lambda mols, **kw: oecluster.cdist(mols[:4], mols[4:],
                                                    "descriptor", **kw),
                 id="cdist"),
    pytest.param(lambda mols, **kw: oecluster.DescriptorComparison(mols, **kw),
                 id="constructor"),
]


@pytest.mark.parametrize("build", _DESCRIPTOR_SURFACES)
@pytest.mark.parametrize("metric", ["euclidean", "standardized_euclidean",
                                    "seuclidean"])
def test_p_is_refused_on_a_metric_it_does_not_apply_to(build, metric):
    """``p`` reached the descriptor path and was silently inert.

    Measured before the fix: the condensed values were bit-identical with and
    without it. ``seuclidean`` is included because it is an alias, and an
    alias must not route around the rule.
    """
    with pytest.raises(TypeError, match=r"p does not apply to metric="):
        build(_mols(), metric=metric, p=3.0)


@pytest.mark.parametrize("build", _DESCRIPTOR_SURFACES)
def test_a_metric_is_read_once_however_unstable_its_str_is(build):
    """The rules and C++ must judge the same text.

    ``metric`` is coerced with ``str`` when the options are built, so a
    selector whose ``__str__`` answers differently on successive calls
    could give C++ one metric and the Python rules another -- reading
    "mahalanobis" into the options while the advisory saw "euclidean"
    and stayed silent on an asymmetric matrix. It is coerced once and
    written back, as the cdist cutoff is.
    """
    class Shifty:
        def __init__(self):
            self.calls = 0

        def __str__(self):
            self.calls += 1
            return "mahalanobis" if self.calls == 1 else "euclidean"

    asymmetric = np.array([[2.0, 1.0], [0.0, 2.0]])
    with pytest.warns(UserWarning, match="inverse_covariance is not symmetric"):
        build(_mols(), columns=_first_columns(2), metric=Shifty(),
              inverse_covariance=asymmetric.ravel().tolist())

    # The sharper case: ``str()`` returns whatever ``__str__`` gives, and
    # that may itself be a ``str`` *subclass*. SWIG stores its payload --
    # "mahalanobis" -- while a later Python reader calling ``str`` on it
    # gets "euclidean", so C++ and the rules judged different metrics.
    # Taking the value back off the options settles it: the getter hands
    # over the plain string C++ holds, so this behaves as "mahalanobis"
    # and the asymmetry is reported rather than slipping past.
    class Masked(str):
        def __str__(self):
            return "euclidean"

    class Indirect(str):
        def __str__(self):
            return Masked("mahalanobis")

    with pytest.warns(UserWarning, match="inverse_covariance is not symmetric"):
        build(_mols(), columns=_first_columns(2),
              metric=Indirect("unused"),
              inverse_covariance=asymmetric.ravel().tolist())


@pytest.mark.parametrize("build", _DESCRIPTOR_SURFACES)
def test_p_is_refused_by_default_because_the_default_is_not_minkowski(build):
    with pytest.raises(TypeError, match=r"p does not apply to metric="):
        build(_mols(), p=3.0)


@pytest.mark.parametrize("build", _DESCRIPTOR_SURFACES)
@pytest.mark.parametrize("metric", ["minkowski", "MINKOWSKI", "Minkowski"])
def test_p_is_accepted_on_minkowski_whatever_its_case(build, metric):
    """Native metric resolution lower-cases (MetricTable.cpp:18), so a rule
    that compares the raw spelling would refuse a valid call."""
    build(_mols(), metric=metric, p=3.0)


@pytest.mark.parametrize("build", _DESCRIPTOR_SURFACES)
@pytest.mark.parametrize("metric,needle", [
    ("bogus", r"Unknown metric 'bogus'"),
    ("tanimoto", r"is a fingerprint bit-set metric"),
    ("haversine", r"great-circle distance over exactly two columns"),
])
def test_a_bad_metric_is_named_before_p_is(build, metric, needle):
    """C++ speaks first on each of these surfaces.

    A message about ``p`` not applying to ``metric='bogus'`` would name a
    metric that does not exist, and ``tanimoto``/``haversine`` are recognised
    natively but unavailable here, each with its own message. Python cannot
    tell these apart -- the metric table is private -- so it must not try.
    """
    with pytest.raises(RuntimeError, match=needle):
        build(_mols(), metric=metric, p=3.0)


@pytest.mark.parametrize("build", _DESCRIPTOR_SURFACES)
def test_any_other_rejected_option_is_reported_before_p(build):
    """A valid metric is not enough to guarantee the ``p`` TypeError.

    The ``p`` rule runs after option conversion and after the native
    validator, so a request that is also wrong in some unrelated way is
    answered by whichever layer refuses it first. Here the metric is valid
    and not minkowski -- the exact case that otherwise raises TypeError --
    but ``missing`` is nonsense, so C++ speaks and the caller gets
    RuntimeError. Pinned because the docstrings promise this precedence.
    """
    with pytest.raises(RuntimeError, match="missing-value policy"):
        build(_mols(), metric="euclidean", p=3.0, missing="bogus")


@pytest.mark.parametrize("build", _DESCRIPTOR_SURFACES)
def test_which_of_two_faults_is_reported_depends_on_the_override(build):
    """There is no fixed precedence here, and the docstrings promise none.

    A source, column or group name resolves only once a schema exists, so
    where that check runs relative to the ``p`` rule is data-dependent:
    supplying ``variances`` resolves the schema early and moves it ahead.
    Both directions are pinned because four attempts to document a fixed
    ordering were each falsified by a case the previous one missed --
    the honest contract is that fixing one fault may expose another.
    """
    # No override: the schema is not resolved yet, so p is reported.
    with pytest.raises(TypeError, match=r"p does not apply to metric="):
        build(_mols(), metric="euclidean", p=3.0, sources=["nope"])

    # An override resolves the schema early, so the source is reported.
    with pytest.raises(RuntimeError, match="[Ss]ource"):
        build(_mols(), metric="standardized_euclidean", p=3.0,
              sources=["nope"], variances=[1.0, 2.0])

    # And with no p at all, the source is reported either way.
    with pytest.raises(RuntimeError, match="[Ss]ource"):
        build(_mols(), metric="euclidean", sources=["nope"])


def test_the_raw_native_bindings_still_accept_an_inert_p():
    """The documented escape hatch sits below this layer's contract.

    ``docs/python-api.md`` publishes ``oecluster.oecluster`` for direct access
    to the C++ options and classes, and a ``DescriptorOptions`` built there
    reaches the comparison without passing any Python rule. ``p`` is inert on
    that route rather than refused, which is measured here: under
    ``metric='euclidean'`` the condensed values are bit-identical at ``p=2.0``
    and ``p=3.0``.

    C++ cannot apply the rule in Python's place. ``DescriptorOptions::p`` is a
    plain ``double`` with a default, so the validator cannot tell "not
    supplied" from "supplied as the default" without a presence flag on a
    public struct -- a design change, and not the one this layer's fix is.
    Pinned rather than closed, so the boundary stays a recorded fact instead
    of an unstated gap behind a claim that the refusal is universal.
    """
    def condensed(p):
        opts = _native.DescriptorOptions()
        opts.metric = "euclidean"
        opts.p = p
        return oecluster.pdist(
            _mols(), _native.DescriptorComparison(_mols(), opts)).condensed

    assert np.array_equal(condensed(2.0), condensed(3.0))
