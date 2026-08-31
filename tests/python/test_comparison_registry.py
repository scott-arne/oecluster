import oecluster
import pytest
from oecluster import _comparisons
from openeye import oechem


def _mols(smiles_list):
    mols = []
    for idx, smi in enumerate(smiles_list):
        mol = oechem.OEGraphMol()
        oechem.OESmilesToMol(mol, smi)
        mol.SetTitle(f"mol{idx}")
        mols.append(mol)
    return mols


def test_extract_labels_uses_titles():
    mols = _mols(["CCO", "CCC"])
    assert _comparisons.extract_labels(mols) == ["mol0", "mol1"]


def test_extract_labels_falls_back_to_index_for_untitled():
    mols = _mols(["CCO", "CCC"])
    mols[1].SetTitle("")
    assert _comparisons.extract_labels(mols) == ["mol0", "mol_1"]


def test_supported_comparisons_lists_the_builtin_names():
    names = _comparisons.supported_comparisons()
    for expected in ("fingerprint", "rocs", "sitehopper", "superpose"):
        assert expected in names


def test_build_comparison_returns_object_name_and_params():
    mols = _mols(["CCO", "CCC", "CCCC"])
    obj, name, params = _comparisons.build_comparison(
        mols, "fingerprint", False, {}, symmetric=True)
    assert name == "fingerprint"
    assert obj.Size() == 3
    assert params == {"comparison_type": "fingerprint", "similarity": False}


def test_build_comparison_is_case_insensitive():
    mols = _mols(["CCO", "CCC"])
    _, name, params = _comparisons.build_comparison(
        mols, "FingerPrint", False, {}, symmetric=True)
    assert name == "fingerprint"
    assert params["comparison_type"] == "fingerprint"


def test_build_comparison_consumes_its_kwargs():
    mols = _mols(["CCO", "CCC"])
    kwargs = {"numbits": 512, "metric": "dice"}
    _comparisons.build_comparison(mols, "fingerprint", False, kwargs,
                                  symmetric=True)
    assert kwargs == {}


def test_build_comparison_rejects_unknown_kwargs():
    mols = _mols(["CCO", "CCC"])
    with pytest.raises(TypeError, match="Unknown kwargs for fingerprint"):
        _comparisons.build_comparison(mols, "fingerprint", False,
                                      {"bogus": 1}, symmetric=True)


def test_build_comparison_rejects_an_unknown_name():
    mols = _mols(["CCO", "CCC"])
    with pytest.raises(ValueError, match="Unknown comparison"):
        _comparisons.build_comparison(mols, "nope", False, {}, symmetric=True)


def test_asymmetric_tversky_is_rejected_for_pdist():
    mols = _mols(["CCO", "CCC"])
    with pytest.raises(ValueError, match="tversky_alpha"):
        _comparisons.build_comparison(
            mols, "fingerprint", True,
            {"metric": "tversky", "tversky_alpha": 0.9, "tversky_beta": 0.1},
            symmetric=True)


def test_asymmetric_tversky_is_allowed_for_cdist():
    mols = _mols(["CCO", "CCC"])
    obj, _, _ = _comparisons.build_comparison(
        mols, "fingerprint", True,
        {"metric": "tversky", "tversky_alpha": 0.9, "tversky_beta": 0.1},
        symmetric=False)
    assert obj.Size() == 2


def test_symmetric_tversky_is_allowed_for_pdist():
    mols = _mols(["CCO", "CCC"])
    obj, _, _ = _comparisons.build_comparison(
        mols, "fingerprint", True,
        {"metric": "tversky", "tversky_alpha": 0.4, "tversky_beta": 0.4},
        symmetric=True)
    assert obj.Size() == 2


def test_numbits_is_rejected_for_a_sparse_storage():
    """A sparse fingerprint has no folding width to set."""
    mols = _mols(["CCO", "CCC"])
    with pytest.raises(TypeError, match="numbits does not apply"):
        _comparisons.build_comparison(
            mols, "fingerprint", False,
            {"storage": "sparse", "numbits": 4096}, symmetric=True)


def test_numbits_is_accepted_for_a_folded_storage():
    mols = _mols(["CCO", "CCC"])
    obj, _, _ = _comparisons.build_comparison(
        mols, "fingerprint", False,
        {"storage": "count", "numbits": 4096, "metric": "bray_curtis"},
        symmetric=True)
    assert obj.Size() == 2


def test_max_distance_is_rejected_for_morgan_naming_radius():
    """The 5.0.0 break: max_distance no longer aliases the Morgan radius."""
    mols = _mols(["CCO", "CCC"])
    with pytest.raises(TypeError, match="radius"):
        _comparisons.build_comparison(
            mols, "fingerprint", False, {"max_distance": 2}, symmetric=True)


def test_radius_is_rejected_for_atom_pair():
    mols = _mols(["CCO", "CCC"])
    with pytest.raises(TypeError, match="min_distance/max_distance"):
        _comparisons.build_comparison(
            mols, "fingerprint", False,
            {"fp_type": "atom_pair", "radius": 3}, symmetric=True)


def test_torsion_atom_count_is_rejected_for_morgan():
    mols = _mols(["CCO", "CCC"])
    with pytest.raises(TypeError, match="torsion_atom_count does not apply"):
        _comparisons.build_comparison(
            mols, "fingerprint", False, {"torsion_atom_count": 5},
            symmetric=True)


def test_use_chirality_is_never_rejected():
    """It applies to all four families, so no family rule may claim it."""
    mols = _mols(["CCO", "CCC"])
    for family in ("morgan", "atom_pair", "topological_atom_pair",
                   "topological_torsions"):
        obj, _, _ = _comparisons.build_comparison(
            mols, "fingerprint", False,
            {"fp_type": family, "use_chirality": True}, symmetric=True)
        assert obj.Size() == 2


def test_p_is_rejected_for_a_non_minkowski_metric():
    mols = _mols(["CCO", "CCC"])
    with pytest.raises(TypeError, match="p does not apply"):
        _comparisons.build_comparison(
            mols, "fingerprint", False, {"p": 3.0}, symmetric=True)


def test_p_is_accepted_for_minkowski():
    mols = _mols(["CCO", "CCC"])
    obj, _, _ = _comparisons.build_comparison(
        mols, "fingerprint", False, {"metric": "minkowski", "p": 3.0},
        symmetric=True)
    assert obj.Size() == 2


def test_tversky_weights_are_rejected_for_a_non_tversky_metric():
    mols = _mols(["CCO", "CCC"])
    with pytest.raises(TypeError, match="tversky_alpha does not apply"):
        _comparisons.build_comparison(
            mols, "fingerprint", False, {"tversky_alpha": 0.3},
            symmetric=True)


def test_an_unknown_family_falls_through_to_the_cpp_error():
    """The explicitness rules must not pre-empt the family validation.

    ``maccs`` is not a family, so ``radius`` has no family to be compared
    against. The C++ constructor owns that message.
    """
    mols = _mols(["CCO", "CCC"])
    with pytest.raises(RuntimeError):
        _comparisons.build_comparison(
            mols, "fingerprint", False, {"fp_type": "maccs", "radius": 3},
            symmetric=True)


def test_an_unknown_family_outranks_the_metric_and_storage_rules():
    """The family error outranks every option rule, not just the family ones.

    ``p``, ``tversky_alpha`` and the sparse ``numbits`` rule are keyed on
    metric and storage rather than on family, so they used to fire before the
    C++ constructor ever saw an unsupported ``fp_type``. Each told the caller
    to pick a different metric or storage -- advice that cannot make such a
    call valid, and that hides the one thing they have to change.
    """
    mols = _mols(["CCO", "CCC"])
    unsupported = (
        {"fp_type": "nonsense", "p": 3.0},
        {"fp_type": "nonsense", "tversky_alpha": 0.5},
        {"fp_type": "maccs", "storage": "sparse", "numbits": 4096},
        {"fp_type": "distance_atom_pair", "p": 3.0},
    )
    for kwargs in unsupported:
        with pytest.raises(RuntimeError):
            _comparisons.build_comparison(
                mols, "fingerprint", False, dict(kwargs), symmetric=True)


def test_an_unknown_family_does_not_swallow_an_unknown_kwarg():
    """Falling through on the family must not fall through on a typo.

    The unknown-kwarg check lives in the builder rather than in the
    explicitness rules, so it still fires for a family the rules skip.
    """
    mols = _mols(["CCO", "CCC"])
    with pytest.raises(TypeError, match=r"\['bogus_kwarg'\]"):
        _comparisons.build_comparison(
            mols, "fingerprint", False,
            {"fp_type": "nonsense", "bogus_kwarg": 1}, symmetric=True)


def test_normalize_items_defaults_to_passthrough():
    mols = _mols(["CCO", "CCC"])
    kept, excluded = _comparisons.normalize_items("fingerprint", mols, {})
    assert kept == mols
    assert excluded == []


def test_normalize_items_rejects_an_unknown_name():
    with pytest.raises(ValueError, match="Unknown comparison"):
        _comparisons.normalize_items("nope", [], {})


def test_normalize_items_runs_a_registered_normalizer():
    seen = {}

    def builder(items, similarity, kwargs, symmetric):
        kwargs.clear()
        return object(), "toy"

    def normalizer(items, kwargs):
        seen["mode"] = kwargs.get("mode", "default")
        return items[:1], [[1, "dropped-by-test"]]

    _comparisons.register_comparison("toy", builder, normalizer)
    try:
        kept, excluded = _comparisons.normalize_items(
            "toy", ["a", "b"], {"mode": "loud"})
    finally:
        _comparisons._BUILDERS.pop("toy", None)
        _comparisons._NORMALIZERS.pop("toy", None)

    assert kept == ["a"]
    assert excluded == [[1, "dropped-by-test"]]
    assert seen["mode"] == "loud"


def test_rocs_builder_rejects_unknown_kwargs():
    with pytest.raises(TypeError, match="Unknown kwargs for rocs"):
        _comparisons.build_comparison(
            ["item"], "rocs", False, {"bogus": 1}, symmetric=True)


def test_superpose_builder_rejects_unknown_kwargs():
    with pytest.raises(TypeError, match="Unknown kwargs for superpose"):
        _comparisons.build_comparison(
            ["item"], "superpose", False, {"bogus": 1}, symmetric=True)


def test_backstop_catches_a_broken_builder_that_leaks_kwargs():
    """A builder that forgets to drain kwargs hits the backstop."""
    def broken_builder(items, similarity, kwargs, symmetric):
        # Deliberately does not pop or check kwargs
        return object(), "broken"

    _comparisons.register_comparison("broken", broken_builder)
    try:
        with pytest.raises(RuntimeError, match=r"'broken' builder.*\['leftover'\]"):
            _comparisons.build_comparison(
                ["item"], "broken", False, {"leftover": 42}, symmetric=True)
    finally:
        _comparisons._BUILDERS.pop("broken", None)


def test_family_aliases_get_the_same_explicitness_rules():
    """An alias spelling must not be a way around the family rules.

    The C++ constructor folds several spellings onto each generator
    (normalize_family, FingerprintComparison.cpp:124). A rule keyed on the
    raw spelling would apply to the canonical name and skip the aliases.
    """
    mols = _mols(["CCO", "CCC"])
    atom_pair_spellings = ("atom_pair", "atompair", "topological_atom_pair")
    torsion_spellings = ("topological_torsions", "topological_torsion")

    for family in atom_pair_spellings:
        with pytest.raises(TypeError, match="radius does not apply"):
            _comparisons.build_comparison(
                mols, "fingerprint", False,
                {"fp_type": family, "radius": 3}, symmetric=True)
        obj, _, _ = _comparisons.build_comparison(
            mols, "fingerprint", False,
            {"fp_type": family, "min_distance": 1, "max_distance": 5},
            symmetric=True)
        assert obj.Size() == 2

    for family in torsion_spellings:
        with pytest.raises(TypeError, match="max_distance does not apply"):
            _comparisons.build_comparison(
                mols, "fingerprint", False,
                {"fp_type": family, "max_distance": 4}, symmetric=True)
        obj, _, _ = _comparisons.build_comparison(
            mols, "fingerprint", False,
            {"fp_type": family, "torsion_atom_count": 4}, symmetric=True)
        assert obj.Size() == 2


def test_the_rejection_message_echoes_the_spelling_the_caller_used():
    mols = _mols(["CCO", "CCC"])
    with pytest.raises(TypeError, match=r"fp_type='atompair'"):
        _comparisons.build_comparison(
            mols, "fingerprint", False,
            {"fp_type": "atompair", "radius": 3}, symmetric=True)


def test_pdist_still_works_through_the_registry():
    mols = _mols(["CCO", "CCC", "CCCC"])
    dist = oecluster.pdist(mols, "fingerprint")
    assert dist.num_samples == 3
    assert dist.labels == ["mol0", "mol1", "mol2"]
    assert dist.params["comparison_type"] == "fingerprint"


def test_cdist_still_works_through_the_registry():
    a = _mols(["CCO", "CCC"])
    b = _mols(["CCCC", "CCCCC", "c1ccccc1"])
    cross = oecluster.cdist(a, b, "fingerprint")
    assert cross.shape == (2, 3)
