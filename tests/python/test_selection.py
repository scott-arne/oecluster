"""Python surface of maxmin_select and circles: dispatch, coordinates, refusals."""

import math

import oecluster
import pytest
from openeye import oechem

# Default Morgan/Tanimoto distances on this set include exact 1.0 ties (no
# shared bits), so the smaller-index tie rule decides several picks, and the
# #Circles packing at 0.75 keeps 5 of the 12 items rather than all or one.
FP_SMILES = ["CCO", "CCCO", "CCCCO", "c1ccccc1", "Cc1ccccc1", "CCc1ccccc1",
             "CC(=O)O", "CC(=O)OC", "CCN", "CCCN", "C1CCCCC1", "c1ccncc1"]

# The descriptor set test_descriptor_comparison.py uses; prefixing water makes
# the complete-case mask drop position 0.
DESCRIPTOR_SMILES = ["CCO", "CCC", "CCCC", "c1ccccc1", "CCN", "CCOC",
                     "CC(=O)O", "CCCCCC"]

# Five points on a line at these coordinates; the distance is the gap.
_LINE = (0.0, 1.0, 3.0, 7.0, 8.0)


def _mols(smiles_list):
    mols = []
    for idx, smi in enumerate(smiles_list):
        mol = oechem.OEGraphMol()
        oechem.OESmilesToMol(mol, smi)
        mol.SetTitle(f"mol{idx}")
        mols.append(mol)
    return mols


def _multiconformer(smiles, shifts, title="mol"):
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


def _dense_distance_matrix(square):
    """Build a dense SymmetricDistanceMatrix from a square distance list."""
    storage = oecluster.DenseStorage(len(square))
    for i in range(len(square)):
        for j in range(i + 1, len(square)):
            storage.Set(i, j, float(square[i][j]))
    return oecluster.SymmetricDistanceMatrix(
        storage, "test", [f"item_{i}" for i in range(len(square))], {})


def _line_matrix():
    return _dense_distance_matrix(
        [[abs(a - b) for b in _LINE] for a in _LINE])


def test_the_three_call_styles_agree():
    mols = _mols(FP_SMILES)
    by_matrix = oecluster.maxmin_select(oecluster.pdist(mols, "fingerprint"),
                                        count=4)
    by_object = oecluster.maxmin_select(oecluster.FingerprintComparison(mols),
                                        count=4)
    by_name = oecluster.maxmin_select(mols, comparison="fingerprint", count=4)

    assert by_matrix.indices == [0, 3, 10, 7]
    assert by_matrix.pick_distances == pytest.approx(
        [math.nan, 1.0, 1.0, 6 / 7], nan_ok=True)
    for other in (by_object, by_name):
        assert other.indices == by_matrix.indices
        assert other.pick_distances == pytest.approx(
            by_matrix.pick_distances, nan_ok=True)
        assert other.stop == "count"
        assert other.excluded == []


def test_a_count_stops_the_selection():
    selection = oecluster.maxmin_select(_line_matrix(), count=3)
    assert selection.indices == [0, 4, 2]
    assert selection.pick_distances == pytest.approx([math.nan, 8.0, 3.0],
                                                     nan_ok=True)
    assert selection.stop == "count"
    assert selection.excluded == []


def test_a_threshold_equal_to_the_best_distance_stops_the_selection():
    """The next candidate is exactly 3 from the selection; at or within the
    threshold is not added."""
    selection = oecluster.maxmin_select(_line_matrix(), threshold=3)
    assert selection.indices == [0, 4]
    assert selection.stop == "threshold"


def test_a_small_threshold_runs_to_exhaustion():
    selection = oecluster.maxmin_select(_line_matrix(), threshold=0.5)
    assert selection.indices == [0, 4, 2, 1, 3]
    assert selection.pick_distances == pytest.approx(
        [math.nan, 8.0, 3.0, 1.0, 1.0], nan_ok=True)
    assert selection.stop == "exhausted"


def test_an_explicit_seed_starts_the_selection():
    selection = oecluster.maxmin_select(_line_matrix(), count=2, seed=3)
    assert selection.indices == [3, 0]
    assert selection.pick_distances == pytest.approx([math.nan, 7.0],
                                                     nan_ok=True)


def test_the_farthest_seed_starts_from_the_item_farthest_from_item_zero():
    selection = oecluster.maxmin_select(_line_matrix(), count=3,
                                        seed="farthest")
    assert selection.indices == [4, 0, 2]


def test_the_medoid_seed_starts_from_the_smallest_distance_sum():
    """Row sums on the line are 19, 16, 14, 18 and 21, so item 2 seeds."""
    selection = oecluster.maxmin_select(_line_matrix(), count=3,
                                        seed="MEDOID")
    assert selection.indices == [2, 4, 0]
    assert selection.pick_distances == pytest.approx([math.nan, 5.0, 3.0],
                                                     nan_ok=True)


def test_an_initial_selection_is_extended():
    selection = oecluster.maxmin_select(_line_matrix(), count=3,
                                        initial=[1, 3])
    assert selection.indices == [1, 3, 2]
    assert selection.pick_distances == pytest.approx(
        [math.nan, math.nan, 2.0], nan_ok=True)


def test_an_initial_selection_of_length_count_is_returned_as_is():
    selection = oecluster.maxmin_select(_line_matrix(), count=2,
                                        initial=[1, 3])
    assert selection.indices == [1, 3]
    assert selection.stop == "count"


def test_an_explicit_seed_with_initial_is_a_type_error():
    """seed=0 is the default's value, but spelling it out is still a
    conflict with initial rather than a no-op."""
    with pytest.raises(TypeError, match="initial or seed, not both"):
        oecluster.maxmin_select(_line_matrix(), count=3, seed=0,
                                initial=[1, 3])


def test_positions_refer_to_the_callers_items_after_normalization_drops_one():
    mols = _mols(["O"] + DESCRIPTOR_SMILES)
    by_name = oecluster.maxmin_select(mols, comparison="descriptor", count=3)
    by_matrix = oecluster.maxmin_select(oecluster.pdist(mols, "descriptor"),
                                        count=3)

    assert by_name.excluded == [[0, "missing-descriptor"]]
    assert by_name.indices == [index + 1 for index in by_matrix.indices]
    # The default seed is the first item that survived normalization.
    assert by_name.indices == [1, 8, 4]


def test_a_dropped_position_cannot_seed_or_start_the_selection():
    mols = _mols(["O"] + DESCRIPTOR_SMILES)
    with pytest.raises(ValueError,
                       match=r"seed 0 names an item that normalization "
                             r"dropped \(missing-descriptor\)"):
        oecluster.maxmin_select(mols, comparison="descriptor", count=3,
                                seed=0)
    with pytest.raises(ValueError, match=r"initial entry 0 names an item"):
        oecluster.maxmin_select(mols, comparison="descriptor", count=3,
                                initial=[0])


def test_a_seed_outside_the_callers_items_is_refused():
    mols = _mols(["O"] + DESCRIPTOR_SMILES)
    with pytest.raises(ValueError, match=r"outside the item range \(0 to 8\)"):
        oecluster.maxmin_select(mols, comparison="descriptor", count=3,
                                seed=9)


def test_conformer_expansion_is_refused():
    mols = [_multiconformer("CCCO", [0.0, 1.0], "a"),
            _multiconformer("CCCO", [0.0, 2.0], "b")]
    with pytest.raises(ValueError, match="expand_conformers=False"):
        oecluster.maxmin_select(mols, comparison="rmsd", count=1)

    selection = oecluster.maxmin_select(mols, comparison="rmsd", count=2,
                                        expand_conformers=False)
    assert selection.indices == [0, 1]


def test_a_cross_distance_matrix_is_a_type_error():
    mols = _mols(FP_SMILES)
    cross = oecluster.cdist(mols[:2], mols[2:4], "fingerprint")
    with pytest.raises(TypeError, match="CrossDistanceMatrix"):
        oecluster.maxmin_select(cross, count=1)


@pytest.mark.parametrize(("missing", "match"), [
    ("propagate", "non-finite"),
    ("ignore", "per-pair feature subsets"),
])
def test_descriptor_missingness_that_cannot_be_ranked_is_refused(missing,
                                                                 match):
    """count=1 never reads a pair, so only the declared facts can refuse."""
    options = {"missing": missing}
    if missing == "ignore":
        options["metric"] = "euclidean"
    mols = _mols(DESCRIPTOR_SMILES)
    with pytest.raises(ValueError, match=match):
        oecluster.maxmin_select(mols, comparison="descriptor", count=1,
                                **options)
    with pytest.raises(ValueError, match=match):
        oecluster.maxmin_select(
            oecluster.DescriptorComparison(mols, **options), count=1)


def test_a_similarity_comparison_is_refused():
    prebuilt = oecluster.FingerprintComparison(_mols(FP_SMILES),
                                               similarity=True)
    with pytest.raises(ValueError, match="requires distances"):
        oecluster.maxmin_select(prebuilt, count=1)


def test_a_similarity_matrix_is_refused_by_the_gate():
    matrix = oecluster.pdist(_mols(FP_SMILES), "fingerprint", similarity=True)
    with pytest.raises(ValueError, match="holds a similarity"):
        oecluster.maxmin_select(matrix, count=2)


def test_a_nonzero_self_distance_matrix_is_refused_by_the_gate():
    matrix = oecluster.pdist(_mols(FP_SMILES), "fingerprint")
    matrix._facts['zero_self'] = False
    with pytest.raises(ValueError, match="zero self-distance"):
        oecluster.maxmin_select(matrix, count=2)


def test_a_nonzero_self_distance_comparison_is_refused():
    """ROCS measures methane's diagonal at 0.5: its colour term is empty.

    count=1 never reads a pair, so only the declared facts can refuse.
    """
    pytest.importorskip("openeye.oeomega")
    from openeye import oeomega

    omega = oeomega.OEOmega()
    omega.SetMaxConfs(1)
    omega.SetStrictStereo(False)
    mols = []
    for idx in range(2):
        mol = oechem.OEMol()
        oechem.OESmilesToMol(mol, "C")
        mol.SetTitle(f"mol{idx}")
        assert omega(mol)
        mols.append(mol)
    prebuilt = oecluster.ROCSComparison(mols)
    with pytest.raises(ValueError, match="zero self-distance"):
        oecluster.maxmin_select(prebuilt, count=1)


def test_a_subset_scored_matrix_is_refused_by_the_gate():
    matrix = oecluster.pdist(_mols(DESCRIPTOR_SMILES), "descriptor",
                             metric="euclidean", missing="ignore")
    with pytest.raises(ValueError, match="per-pair subset"):
        oecluster.maxmin_select(matrix, count=2)


def test_a_non_finite_matrix_entry_is_refused_by_the_gate():
    """Refused up front even though a count of 2 from item 0 never reads the
    (1, 2) entry, and whatever the seed."""
    square = [[0.0, 1.0, 2.0], [1.0, 0.0, math.nan], [2.0, math.nan, 0.0]]
    for seed in (0, "medoid"):
        with pytest.raises(ValueError, match="non-finite"):
            oecluster.maxmin_select(_dense_distance_matrix(square), count=2,
                                    seed=seed)


def test_a_medoid_row_sum_that_overflows_is_a_runtime_error():
    """Every entry is finite, so the gate passes; the native pre-scan is what
    catches a sum that reaches infinity."""
    square = [[0.0, 1e308, 1e308], [1e308, 0.0, 1e308],
              [1e308, 1e308, 0.0]]
    with pytest.raises(RuntimeError, match="overflows"):
        oecluster.maxmin_select(_dense_distance_matrix(square), count=2,
                                seed="medoid")


def test_sparse_storage_is_refused():
    storage = oecluster.SparseStorage(4, 0.5)
    matrix = oecluster.SymmetricDistanceMatrix(storage, "test", list("abcd"),
                                               {})
    with pytest.raises(ValueError, match="SparseStorage is not supported"):
        oecluster.maxmin_select(matrix, count=2)


@pytest.mark.parametrize(("kwargs", "match"), [
    ({}, "requires count, threshold, or both"),
    ({"count": 0}, "count must be at least 1"),
    ({"count": -1}, "count must be at least 1"),
    ({"count": 2.5}, "count must be an integer"),
    ({"count": "3"}, "count must be an integer"),
    ({"count": 6}, r"count must be at most the item count \(5\)"),
    ({"threshold": math.nan}, "not NaN"),
    ({"threshold": math.inf}, "finite"),
    ({"threshold": 10**400}, "finite"),
    ({"threshold": -1.0}, "non-negative"),
    ({"count": 2, "chunk_size": 0}, "chunk_size must be at least 1"),
    ({"count": 2, "chunk_size": 1.5}, "chunk_size must be an integer"),
    ({"count": 2, "num_threads": -1}, "num_threads must be at least 0"),
    ({"count": 2, "num_threads": None}, "num_threads must be an integer"),
    ({"count": 2, "seed": 5}, r"seed 5 is outside the item range \(0 to 4\)"),
    ({"count": 2, "seed": -1}, "seed must be at least 0"),
    ({"count": 2, "seed": "centroid"}, "Unknown seed"),
    ({"count": 2, "seed": 1.0}, "seed must be an integer"),
    ({"count": True}, "count must be an integer"),
    ({"count": 2, "chunk_size": True}, "chunk_size must be an integer"),
    ({"count": 2, "num_threads": False}, "num_threads must be an integer"),
    ({"count": 2, "seed": True}, "seed must be an integer"),
    ({"count": 3, "initial": [1, 1]}, "initial entries must be unique"),
    ({"count": 3, "initial": [1, 5]}, "initial entry 5 is outside"),
    ({"count": 3, "initial": [-1]}, "initial entries must be non-negative"),
    ({"count": 1, "initial": [1, 3]}, "initial holds more entries than count"),
    ({"count": 2, "similarity": True}, "similarity=True is not supported"),
])
def test_invalid_maxmin_arguments_are_value_errors(kwargs, match):
    with pytest.raises(ValueError, match=match):
        oecluster.maxmin_select(_line_matrix(), **kwargs)


def test_a_threshold_only_selection_with_initial_is_accepted():
    selection = oecluster.maxmin_select(_line_matrix(), threshold=2.5,
                                        initial=[1, 3])
    assert selection.indices == [1, 3]
    assert selection.stop == "threshold"


@pytest.mark.parametrize("entry", [1.0, True])
def test_a_non_int_initial_entry_is_a_type_error(entry):
    with pytest.raises(TypeError, match="initial must be a sequence of ints"):
        oecluster.maxmin_select(_line_matrix(), count=2, initial=[entry])


def test_a_one_shot_initial_iterable_is_consumed_exactly_once():
    """A generator or iterator is materialized once; valid entries work."""
    by_iterator = oecluster.maxmin_select(_line_matrix(), count=3,
                                          initial=iter([1, 3]))
    by_list = oecluster.maxmin_select(_line_matrix(), count=3,
                                      initial=[1, 3])
    assert by_iterator.indices == by_list.indices == [1, 3, 2]


def test_an_invalid_one_shot_initial_iterable_is_refused():
    """A generator with invalid entries still raises TypeError."""
    with pytest.raises(TypeError, match="initial must be a sequence of ints"):
        oecluster.maxmin_select(_line_matrix(), count=2, initial=iter([1.5]))
    with pytest.raises(TypeError, match="initial must be a sequence of ints"):
        oecluster.maxmin_select(_line_matrix(), count=2, initial=iter([True]))


def test_the_medoid_seed_needs_a_matrix():
    mols = _mols(FP_SMILES)
    with pytest.raises(ValueError, match="requires a distance matrix"):
        oecluster.maxmin_select(mols, comparison="fingerprint", count=2,
                                seed="medoid")
    with pytest.raises(ValueError, match="requires a distance matrix"):
        oecluster.maxmin_select(oecluster.FingerprintComparison(mols),
                                count=2, seed="medoid")


@pytest.mark.parametrize("call", [
    lambda mols: oecluster.maxmin_select(
        oecluster.pdist(mols, "fingerprint"), count=2,
        comparison="fingerprint"),
    lambda mols: oecluster.maxmin_select(
        oecluster.pdist(mols, "fingerprint"), count=2, radius=1),
    lambda mols: oecluster.maxmin_select(
        oecluster.FingerprintComparison(mols), count=2,
        comparison="fingerprint"),
    lambda mols: oecluster.maxmin_select(mols, count=2),
    lambda mols: oecluster.maxmin_select(mols, count=2, comparison=3),
    lambda mols: oecluster.maxmin_select(mols, count=2,
                                         comparison="fingerprint", bogus=1),
])
def test_arguments_that_fit_no_path_are_type_errors(call):
    with pytest.raises(TypeError):
        call(_mols(FP_SMILES))


def test_an_input_that_normalization_empties_is_refused():
    with pytest.raises(ValueError, match="requires at least one item"):
        oecluster.maxmin_select(_mols(["O"]), comparison="descriptor",
                                count=1)
    with pytest.raises(ValueError, match="requires at least one item"):
        oecluster.maxmin_select([], comparison="fingerprint", threshold=0.5)


def test_comparison_options_are_forwarded():
    """radius=1 changes the distances enough to end the selection one pick
    earlier than the default radius does at the same threshold."""
    mols = _mols(FP_SMILES)
    by_name = oecluster.maxmin_select(mols, comparison="fingerprint",
                                      threshold=0.75, radius=1)
    by_matrix = oecluster.maxmin_select(
        oecluster.pdist(mols, "fingerprint", radius=1), threshold=0.75)
    default = oecluster.maxmin_select(mols, comparison="fingerprint",
                                      threshold=0.75)

    assert by_name.indices == by_matrix.indices == [0, 3, 10, 7]
    assert default.indices == [0, 3, 10, 7, 5]


def test_threading_options_are_forwarded(monkeypatch):
    native = oecluster.oecluster
    real = native.maxmin_select
    seen = []

    def spy(target, options):
        seen.append((options.num_threads, options.chunk_size))
        return real(target, options)

    monkeypatch.setattr(native, "maxmin_select", spy)
    selection = oecluster.maxmin_select(_mols(FP_SMILES),
                                        comparison="fingerprint", count=4,
                                        num_threads=3, chunk_size=2)
    assert seen == [(3, 2)]
    assert selection.indices == [0, 3, 10, 7]


def test_the_selection_repr_names_its_fields():
    selection = oecluster.maxmin_select(_line_matrix(), count=2)
    assert repr(selection) == (
        "MaxMinSelection(indices=[0, 4], stop='count', excluded=0)")


def test_the_diversity_surface_is_exported():
    """Attribute access does not consult __all__, so every test above passes
    with the names missing from it; star-import and API discovery do not."""
    exported = ("maxmin_select", "MaxMinSelection")

    missing = [name for name in exported if name not in oecluster.__all__]
    assert missing == []
    assert all(hasattr(oecluster, name) for name in exported)
