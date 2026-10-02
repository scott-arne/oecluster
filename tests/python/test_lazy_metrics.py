"""Lazy paths of cluster_report, activity_landscape and modelability.

A lazy path reads its distances through a comparison instead of a stored
matrix. Exactness is defined against a matrix filled through the same
comparison's Compare(min, max), so every parity test below builds that matrix
rather than calling pdist, whose batch kernels may differ in the last bit.
"""

import math
import re

import oecluster
import pytest
from openeye import oechem

FP_SMILES = ["CCO", "CCCO", "CCCCO", "c1ccccc1", "Cc1ccccc1", "CCc1ccccc1",
             "CC(=O)O", "CC(=O)OC", "CCN", "CCCN", "C1CCCCC1", "c1ccncc1",
             "CCCCCO", "c1ccc2ccccc2c1", "CC(C)O", "OCCO"]

DESCRIPTOR_SMILES = ["CCO", "CCC", "CCCC", "c1ccccc1", "CCN", "CCOC",
                     "CC(=O)O", "CCCCCC"]

# Thresholds of the form k / 10000 with k coprime to 10. A Tanimoto distance
# on a 2048-bit fingerprint is a ratio whose denominator is at most 2048, so
# it can never equal one of these and always sits at least 1 / (10000 * 2048)
# away -- far outside the 1e-12 by which pdist and Compare may disagree. A
# count therefore cannot flip between the two paths.
_COVERAGE = [0.2537, 0.3561, 0.4519]
_BOUNDARY = 0.3109

_PAIR_RANK_MESSAGE = (
    "cluster_report cannot compute pair-rank indices from a comparison; pass "
    "a SymmetricDistanceMatrix or set compute_pair_rank_indices=False")


def _mols(smiles_list):
    mols = []
    for idx, smi in enumerate(smiles_list):
        mol = oechem.OEGraphMol()
        oechem.OESmilesToMol(mol, smi)
        mol.SetTitle(f"mol{idx}")
        mols.append(mol)
    return mols


def _compare_matrix(comparison):
    """The matrix a lazy path must reproduce: every pair read via Compare."""
    n = comparison.Size()
    storage = oecluster.DenseStorage(n)
    for i in range(n):
        for j in range(i + 1, n):
            storage.Set(i, j, comparison.Compare(i, j))
    return oecluster.SymmetricDistanceMatrix(
        storage, "test", [f"item_{i}" for i in range(n)], {})


def _same(a, b):
    """Equal, NaN matching NaN, recursing into tuples and NamedTuples."""
    if isinstance(a, tuple):
        return (isinstance(b, tuple) and len(a) == len(b)
                and all(_same(x, y) for x, y in zip(a, b, strict=True)))
    if isinstance(a, float) and math.isnan(a):
        return isinstance(b, float) and math.isnan(b)
    return a == b


def _assert_close(a, b):
    """Floats to pdist's rounding, everything else exactly."""
    if isinstance(a, tuple):
        assert isinstance(b, tuple)
        assert len(a) == len(b)
        for x, y in zip(a, b, strict=True):
            _assert_close(x, y)
    elif isinstance(a, float):
        assert a == pytest.approx(b, rel=1e-12, abs=1e-15, nan_ok=True)
    else:
        assert a == b


def _report_fields(report):
    scalars = tuple(getattr(report, name)
                    for name in oecluster.ClusterReport._SCALAR_FIELDS)
    return (scalars, report.coverage_thresholds, report.coverage_at,
            report.noise_coverage_at, report.records, tuple(report.requested),
            report.method)


def _results(matrix):
    # Butina leaves no noise; DBSCAN with min_samples=3 can, which brings the
    # noise accounting and the noise coverage into the comparison.
    return [oecluster.butina(matrix, 0.5),
            oecluster.dbscan(matrix, eps=0.45, min_samples=3)]


@pytest.mark.parametrize("chunk_size", [1, 7, 4096])
@pytest.mark.parametrize("num_threads", [1, 4])
def test_report_three_paths_match_a_compare_filled_matrix(chunk_size,
                                                          num_threads):
    mols = _mols(FP_SMILES)
    prebuilt = oecluster.FingerprintComparison(mols)
    matrix = _compare_matrix(prebuilt)
    for result in _results(matrix):
        expected = _report_fields(oecluster.cluster_report(
            result, matrix, compute_per_cluster_records=True))
        for items, extra in ((prebuilt, {}),
                             (mols, {"comparison": "fingerprint"})):
            report = oecluster.cluster_report(
                result, items, compute_per_cluster_records=True,
                num_threads=num_threads, chunk_size=chunk_size, **extra)
            assert _same(_report_fields(report), expected)


def test_report_alias_matches_the_positional_matrix():
    matrix = oecluster.pdist(_mols(FP_SMILES), "fingerprint")
    result = oecluster.butina(matrix, 0.5)
    assert _same(
        _report_fields(oecluster.cluster_report(result,
                                                distance_matrix=matrix)),
        _report_fields(oecluster.cluster_report(result, matrix)))


@pytest.mark.parametrize("call, message", [
    (lambda r, m, c: oecluster.cluster_report(r, m, distance_matrix=m),
     "cluster_report() got both items and distance_matrix; pass one"),
    (lambda r, m, c: oecluster.cluster_report(r),
     "cluster_report() missing required argument: 'items'"),
    (lambda r, m, c: oecluster.cluster_report(r, distance_matrix=c),
     "cluster_report() expects a SymmetricDistanceMatrix"),
    (lambda r, m, c: oecluster.cluster_report(r, m, comparison="fingerprint"),
     "cluster_report() takes no comparison"),
    (lambda r, m, c: oecluster.cluster_report(r, c, comparison="fingerprint"),
     "cluster_report() takes no comparison"),
    (lambda r, m, c: oecluster.cluster_report(r, m, preseet="tight"),
     "preseet"),
    (lambda r, m, c: oecluster.cluster_report(r, c, metric="dice"),
     "metric"),
    (lambda r, m, c: oecluster.cluster_report(r, _mols(FP_SMILES)),
     "a sequence of items with comparison="),
])
def test_report_arguments_that_fit_no_path_are_type_errors(call, message):
    mols = _mols(FP_SMILES)
    matrix = oecluster.pdist(mols, "fingerprint")
    result = oecluster.butina(matrix, 0.5)
    prebuilt = oecluster.FingerprintComparison(mols)
    with pytest.raises(TypeError, match=re.escape(message)):
        call(result, matrix, prebuilt)


def test_report_refuses_pair_rank_indices_on_a_lazy_path():
    mols = _mols(FP_SMILES)
    matrix = oecluster.pdist(mols, "fingerprint")
    result = oecluster.butina(matrix, 0.5)
    for items, extra in ((oecluster.FingerprintComparison(mols), {}),
                         (mols, {"comparison": "fingerprint"})):
        with pytest.raises(ValueError, match=re.escape(_PAIR_RANK_MESSAGE)):
            oecluster.cluster_report(result, items,
                                     compute_pair_rank_indices=True, **extra)
    report = oecluster.cluster_report(result, matrix,
                                      compute_pair_rank_indices=True)
    assert report.requested.pair_rank_indices is True


def test_a_nonmetric_comparison_needs_allow_nonmetric():
    mols = _mols(FP_SMILES)
    prebuilt = oecluster.FingerprintComparison(mols, metric="dice")
    matrix = _compare_matrix(prebuilt)
    result = oecluster.butina(matrix, 0.5)
    refusal = ("violates the triangle inequality; cluster_report assumes a "
               "metric. Pass allow_nonmetric=True to proceed anyway.")
    for items, extra in ((prebuilt, {}),
                         (mols, {"comparison": "fingerprint",
                                 "metric": "dice"})):
        with pytest.raises(ValueError, match=re.escape(refusal)):
            oecluster.cluster_report(result, items, **extra)
        report = oecluster.cluster_report(result, items,
                                          allow_nonmetric=True, **extra)
        assert _same(_report_fields(report), _report_fields(
            oecluster.cluster_report(result, matrix, allow_nonmetric=True)))


def test_a_subset_scored_comparison_needs_allow_nonmetric():
    mols = _mols(DESCRIPTOR_SMILES)
    options = {"comparison": "descriptor", "metric": "euclidean",
               "missing": "ignore"}
    result = oecluster.butina(
        oecluster.pdist(mols, "descriptor", metric="euclidean",
                        missing="ignore"),
        0.5, allow_nonmetric=True)
    with pytest.raises(ValueError, match="per-pair subset of features"):
        oecluster.cluster_report(result, mols, **options)
    report = oecluster.cluster_report(result, mols, allow_nonmetric=True,
                                      **options)
    assert report.num_samples == len(mols)


def test_a_similarity_comparison_is_refused_even_with_allow_nonmetric():
    mols = _mols(FP_SMILES)
    result = oecluster.butina(oecluster.pdist(mols, "fingerprint"), 0.5)
    prebuilt = oecluster.FingerprintComparison(mols, similarity=True)
    with pytest.raises(ValueError, match="requires distances"):
        oecluster.cluster_report(result, prebuilt, allow_nonmetric=True)


def test_a_non_bool_allow_nonmetric_is_refused_on_a_lazy_path():
    mols = _mols(FP_SMILES)
    result = oecluster.butina(oecluster.pdist(mols, "fingerprint"), 0.5)
    with pytest.raises(TypeError, match="allow_nonmetric must be True or False"):
        oecluster.cluster_report(result, mols, comparison="fingerprint",
                                 allow_nonmetric="yes")


@pytest.mark.parametrize("chunk_size, message", [
    (0, "chunk_size must be at least 1"),
    (2.5, "chunk_size must be an integer"),
])
def test_report_chunk_size_is_validated_on_every_path(chunk_size, message):
    mols = _mols(FP_SMILES)
    matrix = oecluster.pdist(mols, "fingerprint")
    result = oecluster.butina(matrix, 0.5)
    for items, extra in ((matrix, {}),
                         (oecluster.FingerprintComparison(mols), {}),
                         (mols, {"comparison": "fingerprint"})):
        with pytest.raises(ValueError, match=message):
            oecluster.cluster_report(result, items, chunk_size=chunk_size,
                                     **extra)


def test_report_refuses_an_item_normalization_dropped():
    mols = _mols(["O", *DESCRIPTOR_SMILES])
    result = oecluster.butina(
        oecluster.pdist(_mols(FP_SMILES[:9]), "fingerprint"), 0.5)
    with pytest.raises(ValueError,
                       match=r"dropped item 0 \(missing-descriptor\)"):
        # Euclidean is named so that no triangle fact can refuse first.
        oecluster.cluster_report(result, mols, comparison="descriptor",
                                 metric="euclidean")


def test_report_refuses_a_result_over_different_items():
    mols = _mols(FP_SMILES)
    result = oecluster.butina(oecluster.pdist(mols, "fingerprint"), 0.5)
    with pytest.raises(ValueError,
                       match="a result and a comparison over the same items"):
        oecluster.cluster_report(result,
                                 oecluster.FingerprintComparison(mols[:12]))


def test_report_agrees_with_pdist_to_rounding():
    mols = _mols(FP_SMILES)
    matrix = oecluster.pdist(mols, "fingerprint")
    options = {"coverage_thresholds": _COVERAGE,
               "boundary_threshold": _BOUNDARY,
               "compute_per_cluster_records": True}
    for result in _results(matrix):
        expected = oecluster.cluster_report(result, matrix, **options)
        lazy = oecluster.cluster_report(result, mols,
                                        comparison="fingerprint", **options)
        _assert_close(_report_fields(lazy), _report_fields(expected))


def _activity(n):
    # Irregular values with one missing measurement, so cliffs, the RMODI band
    # and a dropped sample all reach the comparison.
    return [math.nan if i == 5 else 4.0 + 0.37 * ((7 * i) % 11)
            for i in range(n)]


def _classes(n):
    names = ("active", "inactive", "unknown")
    return ["" if i % 9 == 4 else names[(i * 7 + i // 5) % 3]
            for i in range(n)]


def _fields(scorecard):
    return tuple(getattr(scorecard, name)
                 for name in type(scorecard).__slots__)


def _landscape(items, **options):
    return oecluster.activity_landscape(
        items, _activity(len(FP_SMILES)), distance_threshold=0.45, **options)


def _modelability(items, **options):
    return oecluster.modelability(items, _classes(len(FP_SMILES)), **options)


@pytest.mark.parametrize("score", [_landscape, _modelability])
@pytest.mark.parametrize("chunk_size", [1, 7, 4096])
@pytest.mark.parametrize("num_threads", [1, 4])
def test_sar_three_paths_match_a_compare_filled_matrix(score, chunk_size,
                                                       num_threads):
    mols = _mols(FP_SMILES)
    prebuilt = oecluster.FingerprintComparison(mols)
    expected = _fields(score(_compare_matrix(prebuilt)))
    for items, extra in ((prebuilt, {}),
                         (mols, {"comparison": "fingerprint"})):
        result = score(items, num_threads=num_threads,
                       chunk_size=chunk_size, **extra)
        assert _same(_fields(result), expected)


@pytest.mark.parametrize("name, second", [
    ("activity_landscape", "activity"),
    ("modelability", "activity_classes"),
])
def test_sar_alias_rules_are_type_errors(name, second):
    function = getattr(oecluster, name)
    mols = _mols(FP_SMILES)
    matrix = oecluster.pdist(mols, "fingerprint")
    prebuilt = oecluster.FingerprintComparison(mols)
    values = (_activity(len(mols)) if second == "activity"
              else _classes(len(mols)))
    assert _same(
        _fields(function(distance_matrix=matrix, **{second: values})),
        _fields(function(matrix, values)))
    cases = [
        (lambda: function(matrix, values, distance_matrix=matrix),
         f"{name}() got both items and distance_matrix; pass one"),
        (lambda: function(**{second: values}),
         f"{name}() missing required argument: 'items'"),
        (lambda: function(distance_matrix=matrix),
         f"{name}() missing required argument: '{second}'"),
        (lambda: function(distance_matrix=prebuilt, **{second: values}),
         f"{name}() expects a SymmetricDistanceMatrix"),
        (lambda: function(prebuilt, values, comparison="fingerprint"),
         f"{name}() takes no comparison"),
        (lambda: function(prebuilt, values, num_theads=2),
         "num_theads"),
        (lambda: function(mols, values),
         "a sequence of items with comparison="),
    ]
    for call, message in cases:
        with pytest.raises(TypeError, match=re.escape(message)):
            call()


@pytest.mark.parametrize("score", [_landscape, _modelability])
def test_sar_refuses_similarity_and_subset_scored_comparisons(score):
    with pytest.raises(ValueError, match="requires distances"):
        score(oecluster.FingerprintComparison(_mols(FP_SMILES),
                                              similarity=True))
    mols = _mols(DESCRIPTOR_SMILES * 2)
    with pytest.raises(ValueError, match="per-pair feature subsets"):
        score(mols, comparison="descriptor", metric="euclidean",
              missing="ignore")


@pytest.mark.parametrize("score", [_landscape, _modelability])
def test_sar_accepts_a_nonmetric_comparison(score):
    prebuilt = oecluster.FingerprintComparison(_mols(FP_SMILES),
                                               metric="dice")
    assert _same(_fields(score(prebuilt)),
                 _fields(score(_compare_matrix(prebuilt))))


@pytest.mark.parametrize("score", [_landscape, _modelability])
@pytest.mark.parametrize("chunk_size, message", [
    (0, "chunk_size must be at least 1"),
    (2.5, "chunk_size must be an integer"),
])
def test_sar_chunk_size_is_validated_on_every_path(score, chunk_size,
                                                   message):
    mols = _mols(FP_SMILES)
    for items, extra in ((oecluster.pdist(mols, "fingerprint"), {}),
                         (oecluster.FingerprintComparison(mols), {}),
                         (mols, {"comparison": "fingerprint"})):
        with pytest.raises(ValueError, match=message):
            score(items, chunk_size=chunk_size, **extra)


def test_sar_refuses_an_item_normalization_dropped():
    # Nine descriptor items -- "O" and the eight that score -- against
    # sixteen activities: the drop is refused before the length check.
    mols = _mols(["O", *DESCRIPTOR_SMILES])
    for score in (_landscape, _modelability):
        with pytest.raises(ValueError,
                           match=r"dropped item 0 \(missing-descriptor\)"):
            score(mols, comparison="descriptor", metric="euclidean")


def _multiconformer(smiles, shifts, title):
    """Build one OEMol carrying a conformer per shift."""
    mol = oechem.OEMol()
    oechem.OESmilesToMol(mol, smiles)
    oechem.OEGenerate2DCoordinates(mol)
    base = oechem.OEFloatArray(3 * mol.NumAtoms())
    mol.GetCoords(base)
    for shift in shifts[1:]:
        moved = oechem.OEFloatArray(list(base))
        for idx in range(0, len(moved), 3):
            moved[idx] += shift
        mol.NewConf(moved)
    mol.SetTitle(title)
    return mol


def test_every_lazy_path_refuses_conformer_expansion():
    # _diversity_source refuses the expansion itself; this pins that all three
    # functions route through it before any length or size check.
    mols = [_multiconformer("CCCO", [0.0, 1.0], "a"),
            _multiconformer("CCCO", [0.0, 2.0], "b")]
    result = oecluster.butina(_compare_matrix(
        oecluster.FingerprintComparison(_mols(FP_SMILES))), 0.5)
    calls = [
        lambda: oecluster.cluster_report(result, mols, comparison="rmsd"),
        lambda: oecluster.activity_landscape(mols, [1.0, 2.0],
                                             comparison="rmsd"),
        lambda: oecluster.modelability(mols, ["a", "b"], comparison="rmsd"),
    ]
    for call in calls:
        with pytest.raises(ValueError, match="expand_conformers=False"):
            call()


def test_sar_length_refusals_name_the_comparison():
    prebuilt = oecluster.FingerprintComparison(_mols(FP_SMILES[:12]))
    with pytest.raises(ValueError, match=re.escape(
            "activity has 16 entries but the comparison covers 12 samples")):
        _landscape(prebuilt)
    with pytest.raises(ValueError, match=re.escape(
            "activity_classes has 16 entries but the comparison covers 12 "
            "samples")):
        _modelability(prebuilt)


def test_sar_agrees_with_pdist_to_rounding():
    mols = _mols(FP_SMILES)
    matrix = oecluster.pdist(mols, "fingerprint")
    activity = _activity(len(mols))
    _assert_close(
        _fields(oecluster.activity_landscape(
            mols, activity, comparison="fingerprint",
            distance_threshold=_BOUNDARY)),
        _fields(oecluster.activity_landscape(
            matrix, activity, distance_threshold=_BOUNDARY)))
    classes = _classes(len(mols))
    _assert_close(
        _fields(oecluster.modelability(mols, classes,
                                       comparison="fingerprint")),
        _fields(oecluster.modelability(matrix, classes)))


def test_the_matrix_gate_refuses_a_non_finite_pair_no_metric_reads():
    """The matrix half of the unread-pair difference.

    Sample 5 has no activity, so no landscape metric reads pair (5, 9), yet
    the matrix gate scans every stored pair and refuses. The lazy half -- a
    comparison whose only non-finite value sits on such a pair succeeds --
    cannot be staged from Python with a real comparison, and is pinned in C++
    by SARComparisonTest.ANonFiniteValueOnAnUnreadPairIsNeverSeen.
    """
    matrix = _compare_matrix(
        oecluster.FingerprintComparison(_mols(FP_SMILES)))
    matrix.storage.Set(5, 9, math.nan)
    with pytest.raises(ValueError, match="non-finite entries"):
        _landscape(matrix)


def _lazy_entry_points(n):
    """Call each lazy entry point on n named-path items, sized consistently."""
    result = oecluster.butina(
        oecluster.pdist(_mols(FP_SMILES[:n]), "fingerprint"), 0.5)
    return [
        lambda mols, **kw: oecluster.cluster_report(result, mols, **kw),
        lambda mols, **kw: oecluster.activity_landscape(
            mols, [4.0 + i for i in range(n)], **kw),
        lambda mols, **kw: oecluster.modelability(
            mols, ["active", "inactive"] * (n // 2) + ["active"] * (n % 2),
            **kw),
    ]


@pytest.mark.parametrize("entry", range(3))
def test_a_drop_is_refused_before_the_comparison_is_built(entry):
    # One survivor cannot fit seuclidean's variances, so the constructor would
    # raise ComparisonError (a RuntimeError) if the drop were not refused
    # first.
    call = _lazy_entry_points(2)[entry]
    with pytest.raises(ValueError,
                       match=r"dropped item 0 \(missing-descriptor\)"):
        call(_mols(["O", "CCO"]), comparison="descriptor",
             metric="seuclidean")


@pytest.mark.parametrize("entry", range(3))
def test_dropping_every_item_names_the_first_dropped(entry):
    call = _lazy_entry_points(2)[entry]
    with pytest.raises(ValueError,
                       match=r"dropped item 0 \(missing-descriptor\)"):
        call(_mols(["O", "O"]), comparison="descriptor", metric="euclidean")


@pytest.mark.parametrize("entry", range(3))
def test_a_declared_non_finite_comparison_is_refused(entry):
    call = _lazy_entry_points(len(DESCRIPTOR_SMILES))[entry]
    with pytest.raises(ValueError, match=re.escape(
            "cannot rank distances the comparison declares may be non-finite "
            "(missing='propagate'); use missing='complete_case'")):
        call(_mols(DESCRIPTOR_SMILES), comparison="descriptor",
             metric="euclidean", missing="propagate")


def _linear(smiles, title, rise):
    """A linear molecule in 3D, whose ROCS self-distance is not zero."""
    mol = oechem.OEMol()
    oechem.OESmilesToMol(mol, smiles)
    coords = []
    for idx in range(mol.NumAtoms()):
        coords += [1.16 * idx, 0.1 * idx, rise * idx]
    mol.SetCoords(oechem.OEFloatArray(coords))
    mol.SetDimension(3)
    mol.SetTitle(title)
    return mol


@pytest.mark.parametrize("entry", range(3))
def test_a_nonzero_self_distance_is_refused(entry):
    call = _lazy_entry_points(2)[entry]
    mols = [_linear("O=C=O", "a", 0.3), _linear("N#N", "b", 0.2)]
    with pytest.raises(ValueError, match="requires a zero self-distance"):
        call(mols, comparison="rocs")
