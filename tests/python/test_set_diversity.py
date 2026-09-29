"""Python surface of vendi_score and logdet_diversity: scores, dispatch, refusals."""

import math

import numpy as np
import oecluster
import pytest
from openeye import oechem

FP_SMILES = ["CCO", "CCCO", "CCCCO", "c1ccccc1", "Cc1ccccc1", "CCc1ccccc1",
             "CC(=O)O", "CC(=O)OC", "CCN", "CCCN", "C1CCCCC1", "c1ccncc1"]

WIDE_SMILES = [*FP_SMILES, "CCOC(=O)C", "c1ccc2ccccc2c1", "OCC(O)CO",
               "CC(C)O", "NCCO", "c1ccoc1", "CCS", "ClCCl"]

# The descriptor set test_selection.py uses; its default distances exceed 1,
# so the complement kernel refuses them.
DESCRIPTOR_SMILES = ["CCO", "CCC", "CCCC", "c1ccccc1", "CCN", "CCOC",
                     "CC(=O)O", "CCCCCC"]


def _mols(smiles_list):
    mols = []
    for idx, smi in enumerate(smiles_list):
        mol = oechem.OEGraphMol()
        oechem.OESmilesToMol(mol, smi)
        mol.SetTitle(f"mol{idx}")
        mols.append(mol)
    return mols


def _dense_distance_matrix(square):
    """Build a dense SymmetricDistanceMatrix from a square distance list."""
    storage = oecluster.DenseStorage(len(square))
    for i in range(len(square)):
        for j in range(i + 1, len(square)):
            storage.Set(i, j, float(square[i][j]))
    return oecluster.SymmetricDistanceMatrix(
        storage, "test", [f"item_{i}" for i in range(len(square))], {})


def _constant_matrix(n, distance):
    return _dense_distance_matrix(
        [[0.0 if i == j else distance for j in range(n)] for i in range(n)])


def _square(matrix):
    """The full distance matrix, built without scipy."""
    n = matrix.num_samples
    square = np.zeros((n, n))
    rows, columns = np.triu_indices(n, 1)
    square[rows, columns] = matrix.condensed
    return square + square.T


def _reference_vendi(kernel, q):
    """score_K from vertaix/Vendi-Score, without its p and normalize options."""
    n = kernel.shape[0]
    w = np.linalg.eigvalsh(kernel / n)
    p = w[w > 0]
    if q == 1:
        return float(np.exp(-(p * np.log(p)).sum()))
    return float(np.exp(np.log((p ** q).sum()) / (1 - q)))


@pytest.mark.parametrize("smiles", [FP_SMILES, WIDE_SMILES])
@pytest.mark.parametrize("order", [1, 2])
def test_vendi_matches_the_reference_on_fingerprints(smiles, order):
    matrix = oecluster.pdist(_mols(smiles), "fingerprint")
    kernel = 1.0 - _square(matrix)
    result = oecluster.vendi_score(matrix, order=order)
    assert result.score == pytest.approx(_reference_vendi(kernel, order),
                                         rel=1e-9, abs=1e-9)
    assert result.order == order
    assert result.size == len(smiles)
    assert result.kernel == "complement"


def test_vendi_matches_the_reference_under_the_laplacian_kernel():
    matrix = oecluster.pdist(_mols(DESCRIPTOR_SMILES), "descriptor")
    kernel = np.exp(-_square(matrix) / 2.0)
    result = oecluster.vendi_score(matrix, kernel="laplacian", bandwidth=2.0)
    assert result.score == pytest.approx(_reference_vendi(kernel, 1),
                                         rel=1e-9, abs=1e-9)
    assert result.kernel == "laplacian"
    assert result.min_eigenvalue == pytest.approx(
        np.linalg.eigvalsh(kernel)[0], abs=1e-9)


def test_logdet_matches_numpy_slogdet():
    matrix = oecluster.pdist(_mols(FP_SMILES), "fingerprint")
    kernel = 1.0 - _square(matrix)
    sign, expected = np.linalg.slogdet(kernel + 0.1 * np.eye(len(FP_SMILES)))
    result = oecluster.logdet_diversity(matrix, ridge=0.1)
    assert sign == 1.0
    assert result.score == pytest.approx(expected, rel=1e-9, abs=1e-9)
    assert result.nonpositive_count == 0
    assert result.ridge == 0.1
    assert result.size == len(FP_SMILES)
    assert result.kernel == "complement"
    assert result.min_eigenvalue == pytest.approx(
        np.linalg.eigvalsh(kernel)[0], abs=1e-9)
    assert result.excluded == []


@pytest.mark.parametrize(("function", "options"), [
    ("vendi_score", {"order": 1}),
    ("vendi_score", {"order": 2}),
    ("logdet_diversity", {"ridge": 0.1}),
])
def test_the_three_call_styles_agree(function, options):
    score = getattr(oecluster, function)
    mols = _mols(FP_SMILES)
    by_matrix = score(oecluster.pdist(mols, "fingerprint"), **options)
    by_object = score(oecluster.FingerprintComparison(mols), **options)
    by_name = score(mols, comparison="fingerprint", **options)
    for other in (by_object, by_name):
        assert other.score == pytest.approx(by_matrix.score, rel=1e-12)
        assert other.size == by_matrix.size
        assert other.excluded == []


def test_normalization_drops_are_reported_in_excluded():
    mols = _mols(["O", *DESCRIPTOR_SMILES])
    by_name = oecluster.vendi_score(mols, comparison="descriptor",
                                    kernel="laplacian", bandwidth=2.0)
    by_matrix = oecluster.vendi_score(oecluster.pdist(mols, "descriptor"),
                                      kernel="laplacian", bandwidth=2.0)
    assert by_name.excluded == [[0, "missing-descriptor"]]
    assert by_name.size == 8
    assert by_name.score == pytest.approx(by_matrix.score, rel=1e-9)


def test_identical_items_score_one_and_a_singular_logdet():
    matrix = _constant_matrix(4, 0.0)
    assert oecluster.vendi_score(matrix).score == pytest.approx(1.0)
    assert oecluster.vendi_score(matrix, order=2).score == 1.0
    singular = oecluster.logdet_diversity(matrix)
    assert singular.score == -math.inf
    assert singular.nonpositive_count == 3
    ridged = oecluster.logdet_diversity(matrix, ridge=0.5)
    assert ridged.score == pytest.approx(3 * math.log(0.5) + math.log(4.5))
    assert ridged.nonpositive_count == 0


def test_the_exact_ceiling_names_order_2():
    mols = _mols(FP_SMILES)
    with pytest.raises(ValueError, match=r"order=2"):
        oecluster.vendi_score(mols, comparison="fingerprint", max_exact=5)
    with pytest.raises(ValueError, match=r"order=2"):
        oecluster.logdet_diversity(mols, comparison="fingerprint",
                                   max_exact=5)
    assert oecluster.vendi_score(mols, comparison="fingerprint",
                                 max_exact=12).size == 12
    assert oecluster.vendi_score(mols, comparison="fingerprint", order=2,
                                 max_exact=5).size == 12


@pytest.mark.parametrize(("function", "kwargs", "error", "match"), [
    ("vendi_score", {"order": 0}, ValueError, "order must be 1 or 2"),
    ("vendi_score", {"order": 3}, ValueError, "order must be 1 or 2"),
    ("vendi_score", {"order": True}, ValueError, "order must be 1 or 2"),
    ("vendi_score", {"order": 1.0}, ValueError, "order must be 1 or 2"),
    ("vendi_score", {"order": "1"}, ValueError, "order must be 1 or 2"),
    ("vendi_score", {"kernel": "gaussian"}, ValueError, "Unknown kernel"),
    ("vendi_score", {"kernel": None}, ValueError, "Unknown kernel"),
    ("vendi_score", {"bandwidth": 1.0}, ValueError,
     "bandwidth applies only to kernel='laplacian'"),
    ("vendi_score", {"kernel": "laplacian"}, ValueError,
     "kernel='laplacian' requires a bandwidth"),
    ("vendi_score", {"kernel": "laplacian", "bandwidth": 0.0}, ValueError,
     "bandwidth must be positive and finite"),
    ("vendi_score", {"kernel": "laplacian", "bandwidth": -1.0}, ValueError,
     "bandwidth must be positive and finite"),
    ("vendi_score", {"kernel": "laplacian", "bandwidth": math.inf},
     ValueError, "bandwidth must be positive and finite"),
    ("vendi_score", {"kernel": "laplacian", "bandwidth": math.nan},
     ValueError, "bandwidth must be positive and finite"),
    ("vendi_score", {"kernel": "laplacian", "bandwidth": "1"}, TypeError,
     "bandwidth must be a number"),
    ("vendi_score", {"kernel": "laplacian", "bandwidth": True}, TypeError,
     "bandwidth must be a number"),
    ("vendi_score", {"max_exact": 0}, ValueError, "max_exact must be at least 1"),
    ("vendi_score", {"num_threads": -1}, ValueError,
     "num_threads must be at least 0"),
    ("vendi_score", {"chunk_size": 0}, ValueError,
     "chunk_size must be at least 1"),
    ("vendi_score", {"similarity": True}, ValueError, "similarity=True"),
    ("logdet_diversity", {"ridge": -1.0}, ValueError,
     "ridge must be finite and non-negative"),
    ("logdet_diversity", {"ridge": math.inf}, ValueError,
     "ridge must be finite and non-negative"),
    ("logdet_diversity", {"ridge": math.nan}, ValueError,
     "ridge must be finite and non-negative"),
    ("logdet_diversity", {"ridge": "0"}, TypeError, "ridge must be a number"),
    ("logdet_diversity", {"kernel": "gaussian"}, ValueError, "Unknown kernel"),
    ("logdet_diversity", {"max_exact": 0}, ValueError,
     "max_exact must be at least 1"),
    ("logdet_diversity", {"similarity": True}, ValueError, "similarity=True"),
])
def test_invalid_arguments_are_refused_before_native_code(function, kwargs,
                                                          error, match):
    with pytest.raises(error, match=match):
        getattr(oecluster, function)(_constant_matrix(3, 0.5), **kwargs)


def test_arguments_that_fit_no_path_are_type_errors():
    mols = _mols(FP_SMILES)
    with pytest.raises(TypeError, match="CrossDistanceMatrix"):
        oecluster.vendi_score(oecluster.cdist(mols[:2], mols[2:4],
                                              "fingerprint"))
    with pytest.raises(TypeError, match="takes no comparison"):
        oecluster.logdet_diversity(_constant_matrix(3, 0.5),
                                   comparison="fingerprint")
    with pytest.raises(TypeError, match="requires comparison="):
        oecluster.vendi_score(mols)


@pytest.mark.parametrize("function", ["vendi_score", "logdet_diversity"])
def test_an_empty_matrix_is_refused(function):
    with pytest.raises(ValueError, match="requires at least one item"):
        getattr(oecluster, function)(_dense_distance_matrix([]))


def test_the_matrix_path_checks_the_kernel_range_up_front():
    too_far = _dense_distance_matrix([[0, 0.2, 1.5], [0.2, 0, 0.3],
                                      [1.5, 0.3, 0]])
    negative = _dense_distance_matrix([[0, 0.2, -0.25], [0.2, 0, 0.3],
                                       [-0.25, 0.3, 0]])
    with pytest.raises(ValueError,
                       match=r"requires distances in \[0, 1\], but "
                             r"d\(0, 2\) = 1\.5"):
        oecluster.vendi_score(too_far)
    with pytest.raises(ValueError, match=r"d\(0, 2\) = -0\.25"):
        oecluster.logdet_diversity(negative)
    with pytest.raises(ValueError,
                       match=r"requires non-negative distances, but "
                             r"d\(0, 2\) = -0\.25"):
        oecluster.vendi_score(negative, kernel="laplacian", bandwidth=1.0)
    all_bad = _constant_matrix(4, 1.5)
    with pytest.raises(ValueError, match=r"d\(0, 1\) = 1\.5"):
        oecluster.vendi_score(all_bad)
    # The Laplacian kernel takes distances above 1.
    assert oecluster.vendi_score(too_far, kernel="laplacian",
                                 bandwidth=1.0).size == 3


def test_the_comparison_paths_refuse_out_of_range_distances_at_run_time():
    mols = _mols(DESCRIPTOR_SMILES)
    with pytest.raises(RuntimeError, match=r"requires distances in \[0, 1\]"):
        oecluster.vendi_score(mols, comparison="descriptor")
    with pytest.raises(RuntimeError, match=r"requires distances in \[0, 1\]"):
        oecluster.logdet_diversity(oecluster.DescriptorComparison(mols))


def test_a_similarity_comparison_is_refused():
    prebuilt = oecluster.FingerprintComparison(_mols(FP_SMILES),
                                               similarity=True)
    with pytest.raises(ValueError, match="requires distances"):
        oecluster.vendi_score(prebuilt)


def test_threading_options_are_forwarded(monkeypatch):
    native = oecluster.oecluster
    real = native.vendi_score
    seen = []

    def spy(target, options):
        seen.append((options.num_threads, options.chunk_size,
                     options.max_exact))
        return real(target, options)

    monkeypatch.setattr(native, "vendi_score", spy)
    result = oecluster.vendi_score(_mols(FP_SMILES), comparison="fingerprint",
                                   num_threads=3, chunk_size=2, max_exact=40)
    assert seen == [(3, 2, 40)]
    assert result.size == 12


def test_order_two_leaves_the_spectrum_fields_empty():
    result = oecluster.vendi_score(_constant_matrix(3, 1.0), order=2)
    assert result.score == 3.0
    assert result.min_eigenvalue is None
    assert result.negative_mass is None


def test_the_result_reprs_name_their_fields():
    matrix = _constant_matrix(3, 1.0)
    assert repr(oecluster.vendi_score(matrix, order=2)) == (
        "VendiResult(score=3.0, order=2, size=3, kernel='complement', "
        "excluded=0)")
    assert repr(oecluster.logdet_diversity(matrix)) == (
        "LogDetResult(score=0.0, ridge=0.0, size=3, nonpositive_count=0, "
        "excluded=0)")


def test_the_set_diversity_surface_is_exported():
    """Attribute access does not consult __all__, so every test above passes
    with the names missing from it; star-import and API discovery do not."""
    exported = ("vendi_score", "logdet_diversity", "VendiResult",
                "LogDetResult")
    missing = [name for name in exported if name not in oecluster.__all__]
    assert missing == []
    assert all(hasattr(oecluster, name) for name in exported)
