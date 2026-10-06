"""Butina and DBSCAN from a comparison, holding no matrix.

Exactness is defined against a matrix filled through the same comparison's
Compare(min, max), so every parity test builds that matrix rather than
calling pdist, whose batched kernels may differ in the last bit.
"""

import inspect

import oecluster
import pytest
from openeye import oechem

FP_SMILES = ["CCO", "CCCO", "CCCCO", "c1ccccc1", "Cc1ccccc1", "CCc1ccccc1",
             "CC(=O)O", "CC(=O)OC", "CCN", "CCCN", "C1CCCCC1", "c1ccncc1",
             "CCCCCO", "c1ccc2ccccc2c1", "CC(C)O", "OCCO"]

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


def _compare_matrix(comparison):
    """The matrix a comparison form must reproduce: every pair via Compare."""
    n = comparison.Size()
    storage = oecluster.DenseStorage(n)
    for i in range(n):
        for j in range(i + 1, n):
            storage.Set(i, j, comparison.Compare(i, j))
    return oecluster.SymmetricDistanceMatrix(
        storage, "test", [f"item_{i}" for i in range(n)], {})


def _butina_fields(result):
    return result.labels.tolist(), result.clusters


def _dbscan_fields(result):
    return (result.labels.tolist(), result.clusters,
            result.core_sample_indices)


@pytest.mark.parametrize("chunk_size", [0, 1, 7, 4096])
@pytest.mark.parametrize("num_threads", [1, 4])
@pytest.mark.parametrize("reordering", [False, True])
def test_butina_three_inputs_match_a_compare_filled_matrix(
        reordering, num_threads, chunk_size):
    mols = _mols(FP_SMILES)
    prebuilt = oecluster.FingerprintComparison(mols)
    matrix = _compare_matrix(prebuilt)
    for threshold in (0.3, 0.55, 0.8):
        expected = _butina_fields(
            oecluster.butina(matrix, threshold, reordering=reordering))
        for items, extra in ((prebuilt, {}),
                             (mols, {"comparison": "fingerprint"})):
            result = oecluster.butina(
                items, threshold, reordering=reordering,
                num_threads=num_threads, chunk_size=chunk_size, **extra)
            assert _butina_fields(result) == expected


@pytest.mark.parametrize("chunk_size", [0, 1, 7, 4096])
@pytest.mark.parametrize("num_threads", [1, 4])
def test_dbscan_three_inputs_match_a_compare_filled_matrix(num_threads,
                                                           chunk_size):
    mols = _mols(FP_SMILES)
    prebuilt = oecluster.FingerprintComparison(mols)
    matrix = _compare_matrix(prebuilt)
    for eps, min_samples in ((0.3, 2), (0.55, 3), (0.8, 5), (0.55, 17)):
        expected = _dbscan_fields(
            oecluster.dbscan(matrix, eps, min_samples=min_samples))
        for items, extra in ((prebuilt, {}),
                             (mols, {"comparison": "fingerprint"})):
            result = oecluster.dbscan(
                items, eps, min_samples=min_samples, num_threads=num_threads,
                chunk_size=chunk_size, **extra)
            assert _dbscan_fields(result) == expected


def test_the_alias_still_takes_a_matrix():
    matrix = oecluster.pdist(_mols(FP_SMILES), "fingerprint")
    assert (_butina_fields(oecluster.butina(distance_matrix=matrix,
                                            threshold=0.55))
            == _butina_fields(oecluster.butina(matrix, 0.55)))
    assert (_butina_fields(oecluster.butina(matrix, threshold=0.55))
            == _butina_fields(oecluster.butina(matrix, 0.55)))
    assert (_dbscan_fields(oecluster.dbscan(distance_matrix=matrix, eps=0.55))
            == _dbscan_fields(oecluster.dbscan(matrix, 0.55)))


@pytest.mark.parametrize(("entry", "second"), [("butina", "threshold"),
                                               ("dbscan", "eps")])
def test_the_input_and_the_distance_are_required(entry, second):
    function = getattr(oecluster, entry)
    matrix = oecluster.pdist(_mols(FP_SMILES[:4]), "fingerprint")
    with pytest.raises(TypeError, match="got both items and distance_matrix"):
        function(matrix, 0.5, distance_matrix=matrix)
    with pytest.raises(TypeError, match="missing required argument: 'items'"):
        function(**{second: 0.5})
    with pytest.raises(TypeError,
                       match=f"missing required argument: '{second}'"):
        function(matrix)
    with pytest.raises(TypeError, match="expects a SymmetricDistanceMatrix"):
        function(distance_matrix=_mols(FP_SMILES[:4]), **{second: 0.5})


def test_there_is_one_sentinel_shared_by_every_entry_point():
    # A second _MISSING = object() would leave the older entry points
    # comparing against a different object and silently treating "not given"
    # as a real argument.
    sentinel = oecluster._MISSING
    for function, names in ((oecluster.butina, ("items", "threshold",
                                                "distance_matrix")),
                            (oecluster.dbscan, ("items", "eps",
                                                "distance_matrix")),
                            (oecluster.cluster_report, ("items",
                                                        "distance_matrix")),
                            (oecluster.activity_landscape, ("items",
                                                            "activity")),
                            (oecluster.modelability, ("items",
                                                      "activity_classes"))):
        parameters = inspect.signature(function).parameters
        for name in names:
            assert parameters[name].default is sentinel, (function, name)

    matrix = oecluster.pdist(_mols(FP_SMILES[:4]), "fingerprint")
    result = oecluster.butina(matrix, 0.5)
    with pytest.raises(TypeError, match="missing required argument: 'items'"):
        oecluster.cluster_report(result)
    with pytest.raises(TypeError,
                       match="missing required argument: 'activity'"):
        oecluster.activity_landscape(matrix)
    with pytest.raises(TypeError,
                       match="missing required argument: 'activity_classes'"):
        oecluster.modelability(matrix)


def test_a_dropped_item_is_refused():
    mols = _mols(["O", *DESCRIPTOR_SMILES])
    for function, second in ((oecluster.butina, 0.5),
                             (oecluster.dbscan, 0.5)):
        with pytest.raises(ValueError,
                           match=r"dropped item 0 \(missing-descriptor\)"):
            # Euclidean is named so that no triangle fact can refuse first.
            function(mols, second, comparison="descriptor", metric="euclidean")


def test_a_truthy_allow_nonmetric_is_refused_before_dispatch():
    # "nonexistent" names no comparison, so had dispatch run first the error
    # would be about the name, not about allow_nonmetric.
    mols = _mols(FP_SMILES[:4])
    for function in (oecluster.butina, oecluster.dbscan):
        with pytest.raises(TypeError, match="allow_nonmetric must be True or"):
            function(mols, 0.5, comparison="nonexistent",
                     allow_nonmetric="False")


def test_a_nonmetric_comparison_needs_allow_nonmetric():
    mols = _mols(FP_SMILES)
    for function in (oecluster.butina, oecluster.dbscan):
        with pytest.raises(ValueError, match="violates the triangle inequality"):
            function(mols, 0.5, comparison="fingerprint", metric="dice")
        result = function(mols, 0.5, comparison="fingerprint", metric="dice",
                          allow_nonmetric=True)
        assert len(result.labels) == len(mols)


def test_similarity_is_refused():
    matrix = oecluster.pdist(_mols(FP_SMILES[:4]), "fingerprint")
    for function in (oecluster.butina, oecluster.dbscan):
        with pytest.raises(ValueError, match="similarity=True is not supported"):
            function(matrix, 0.5, similarity=True)


def test_a_budget_is_accepted_on_the_comparison_forms():
    mols = _mols(FP_SMILES)
    expected = _butina_fields(oecluster.butina(mols, 0.55,
                                               comparison="fingerprint"))
    for budget in (None, 1 << 40):
        assert _butina_fields(oecluster.butina(
            mols, 0.55, comparison="fingerprint",
            max_graph_bytes=budget)) == expected
    assert len(oecluster.dbscan(oecluster.FingerprintComparison(mols), 0.55,
                                max_graph_bytes=1 << 40).labels) == len(mols)


@pytest.mark.parametrize(("value", "error", "message"), [
    (True, TypeError, "not a bool"),
    (False, TypeError, "not a bool"),
    (1.5, TypeError, "not float"),
    ("1024", TypeError, "not str"),
    (0, ValueError, "must be positive, got 0"),
    (-1, ValueError, "must be positive, got -1"),
    (1 << 64, ValueError, "exceeds size_t maximum"),
])
def test_a_malformed_budget_is_refused(value, error, message):
    mols = _mols(FP_SMILES[:4])
    for function in (oecluster.butina, oecluster.dbscan):
        with pytest.raises(error, match=message):
            function(mols, 0.5, comparison="fingerprint",
                     max_graph_bytes=value)


@pytest.mark.parametrize("value", [True, 0, -1, 1.5, 1024])
def test_any_budget_beside_a_matrix_is_refused_first(value):
    # The matrix path is deliberately untouched, so any value is refused,
    # ahead of the bool and range checks it would otherwise fail.
    matrix = oecluster.pdist(_mols(FP_SMILES[:4]), "fingerprint")
    for function in (oecluster.butina, oecluster.dbscan):
        with pytest.raises(TypeError, match="no graph budget"):
            function(matrix, 0.5, max_graph_bytes=value)


def test_rocs_is_refused_before_it_is_built():
    # ROCS fails the repeatability the two graph passes need
    # (tests/cpp/test_comparison_repeatability.cpp), so both the named and the
    # prebuilt forms are refused; the matrix path still takes ROCS distances.
    pytest.importorskip("openeye.oeomega")
    from openeye import oeomega

    omega = oeomega.OEOmega()
    omega.SetMaxConfs(1)
    omega.SetStrictStereo(False)
    mols = []
    for smi in ("c1ccccc1", "Cc1ccccc1", "c1ccc(O)cc1"):
        mol = oechem.OEMol()
        oechem.OESmilesToMol(mol, smi)
        assert omega(mol)
        mols.append(mol)
    prebuilt = oecluster.ROCSComparison(mols)
    for function in (oecluster.butina, oecluster.dbscan):
        for items, extra in ((mols, {"comparison": "ROCS"}), (prebuilt, {})):
            with pytest.raises(ValueError,
                               match="cannot cluster a ROCS comparison without "
                                     "a matrix"):
                function(items, 0.5, **extra)
        matrix = oecluster.pdist(mols, "rocs")
        assert len(function(matrix, 0.5, allow_nonmetric=True).labels) == 3


def test_a_named_rocs_comparison_is_refused_before_any_item_is_read():
    # Not molecules at all: had dispatch run first, normalizing these for
    # ROCS would fail with an unrelated error, so only a refusal that
    # precedes dispatch can produce this one.
    for function in (oecluster.butina, oecluster.dbscan):
        with pytest.raises(ValueError,
                           match="cannot cluster a ROCS comparison without "
                                 "a matrix"):
            function(["not", "molecules"], 0.5, comparison="rocs")


def test_an_oversized_graph_is_a_memory_error():
    mols = _mols(FP_SMILES)
    for function in (oecluster.butina, oecluster.dbscan):
        with pytest.raises(MemoryError,
                           match=r"for 16 items and \d+ edges, above its "
                                 r"max_graph_bytes limit of 64 bytes"):
            function(mols, 0.9, comparison="fingerprint", max_graph_bytes=64)
