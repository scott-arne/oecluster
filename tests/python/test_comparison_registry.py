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


def _multiconf_mols(smiles_list):
    """ROCS and superpose need ``OEMol``; their typemaps reject ``OEGraphMol``."""
    mols = []
    for idx, smi in enumerate(smiles_list):
        mol = oechem.OEMol()
        oechem.OESmilesToMol(mol, smi)
        mol.SetTitle(f"mol{idx}")
        mols.append(mol)
    return mols


def _conformer_mols(smiles_list):
    """The same, embedded in 3D, for the builders that require coordinates.

    ``ROCSComparison`` refuses a molecule whose recomputed OEChem dimension
    attribute -- an axis count, not a geometric rank -- is below three, and a
    molecule straight from a SMILES parse carries no coordinates at all, so
    these fixtures have to come from Omega.
    """
    pytest.importorskip("openeye.oeomega")
    from openeye import oeomega

    omega = oeomega.OEOmega()
    omega.SetMaxConfs(1)
    omega.SetStrictStereo(False)

    mols = _multiconf_mols(smiles_list)
    for mol in mols:
        assert omega(mol)
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


def test_an_unknown_family_outranks_the_pdist_tversky_guard():
    """The last Python guard before construction has to stand aside too.

    The asymmetric-Tversky check lives in the builder rather than in the
    explicitness rules, so the early return added for the family rules did not
    cover it. Its advice -- equal weights, or cdist -- cannot make an
    unsupported family valid.
    """
    mols = _mols(["CCO", "CCC"])
    with pytest.raises(RuntimeError):
        oecluster.pdist(mols, "fingerprint", fp_type="nonsense",
                        metric="tversky", tversky_alpha=0.9, tversky_beta=0.1)


def test_a_recognized_family_still_gets_the_pdist_tversky_guard():
    """The control for the test above, through the public entry point.

    ``test_asymmetric_tversky_is_rejected_for_pdist`` asserts this through
    ``build_comparison``; skipping the guard for unrecognized families must not
    skip it for real ones.
    """
    mols = _mols(["CCO", "CCC"])
    with pytest.raises(ValueError, match="pdist requires a symmetric metric"):
        oecluster.pdist(mols, "fingerprint", metric="tversky",
                        tversky_alpha=0.9, tversky_beta=0.1)


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


def test_none_means_unspecified_on_the_public_entry_points():
    """``None`` is the sentinel the public constructors already use.

    ``FingerprintComparison.__new__`` (``__init__.py:2409``) skips the
    assignment when an option is ``None``. The registry builders have to agree,
    or a config-driven caller passing ``fp_type=None`` gets a SWIG internals
    message instead of the default fingerprint.
    """
    mols = _mols(["CCO", "CCC", "c1ccccc1"])
    # Both results are bound to a name rather than read off a temporary:
    # ``condensed`` is a non-owning view over the storage buffer, so reading it
    # from a result that has already been collected yields freed memory.
    baseline_result = oecluster.pdist(mols, "fingerprint")
    explicit_result = oecluster.pdist(
        mols, "fingerprint", fp_type=None, storage=None, metric=None,
        numbits=None, radius=None, use_chirality=None)
    assert list(explicit_result.condensed) == list(baseline_result.condensed)


def test_a_none_valued_option_is_not_a_named_option():
    """A rule must not fire on an option the caller declined to set.

    ``numbits=None`` means "no numbits", so the sparse-storage rule has nothing
    to reject, and ``radius=None`` is not a Morgan option intruding on
    ``atom_pair``.
    """
    mols = _mols(["CCO", "CCC"])
    assert oecluster.pdist(
        mols, "fingerprint", storage="sparse", numbits=None).num_samples == 2
    assert oecluster.pdist(
        mols, "fingerprint", fp_type="atom_pair", radius=None).num_samples == 2


def test_none_is_unspecified_for_rocs_too():
    """The same convention, in the builder next door."""
    # Embedded, because ROCS refuses input without 3D coordinates. This test is
    # about kwarg plumbing and never reads a score, but it still has to hand the
    # builder molecules the builder will accept. Kept separate from the superpose
    # case below so that the Omega dependency this brings in cannot skip a test
    # that does not need it.
    mols = _conformer_mols(["CCO", "CCC"])
    obj, _, _ = _comparisons.build_comparison(
        mols, "rocs", False, {"score_type": None, "color_ff_type": None},
        symmetric=True)
    assert obj.Size() == 2


def test_rocs_refuses_an_extreme_coordinate_extent():
    """Through the public ``pdist``, because the alternative was a crash.

    Displacing one atom to x = 1e10 leaves the OEChem dimension attribute at 3,
    so the guard above this one admits it, and the constructor-time self-overlay
    then took the interpreter down: this exact call exited 139 (SIGSEGV) before
    the extent guard existed, which is not a failure a caller can catch.
    """
    mols = _conformer_mols(["c1ccc(O)cc1", "c1ccccc1"])
    atom = next(iter(mols[0].GetAtoms()))
    coords = list(mols[0].GetCoords(atom))
    coords[0] = 1e10
    assert mols[0].SetCoords(atom, coords)

    with pytest.raises(RuntimeError, match="extent"):
        oecluster.pdist(mols, "rocs")


def test_none_is_unspecified_for_superpose_too():
    """And in the one next to that, which needs no coordinates to construct."""
    mols = _multiconf_mols(["CCO", "CCC"])
    obj, _, _ = _comparisons.build_comparison(
        mols, "superpose", False, {"score_type": None, "predicate": None},
        symmetric=True)
    assert obj.Size() == 2


def test_an_unknown_kwarg_is_still_unknown_when_its_value_is_none():
    """Treating ``None`` as unspecified must not turn a typo into a default."""
    mols = _mols(["CCO", "CCC"])
    with pytest.raises(TypeError, match=r"\['bogus'\]"):
        oecluster.pdist(mols, "fingerprint", bogus=None)


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


def test_an_explicit_empty_method_is_rejected_not_defaulted():
    """An explicit value is used or rejected by name, never discarded."""
    mols = _multiconf_mols(["CCO", "CCC"])
    with pytest.raises(ValueError, match="Unknown superpose method"):
        oecluster.pdist(mols, "superpose", method="")


def test_the_sitehopper_alias_keeps_its_method_fallback():
    """The alias implies its own method, so an unusable value falls back to it.

    This is the base chain's behavior and is deliberately preserved: unlike
    plain superpose, naming the sitehopper comparison names the method.
    """
    mols = _multiconf_mols(["CCO", "CCC"])
    with pytest.raises(Exception) as excinfo:
        oecluster.pdist(mols, "sitehopper", method="")
    # Superposition of two small organics fails downstream; what matters is
    # that it got past method resolution rather than raising Unknown method.
    assert "Unknown superpose method" not in str(excinfo.value)


def test_an_empty_fp_type_does_not_masquerade_as_morgan():
    mols = _mols(["CCO", "CCC", "c1ccccc1"])
    with pytest.raises(RuntimeError, match="Unknown OEFP fingerprint type"):
        oecluster.pdist(mols, "fingerprint", fp_type="", min_distance=1)


def test_an_empty_metric_does_not_masquerade_as_tanimoto():
    mols = _mols(["CCO", "CCC", "c1ccccc1"])
    with pytest.raises(RuntimeError, match="Unknown metric"):
        oecluster.pdist(mols, "fingerprint", metric="", p=3.0)


def test_an_unknown_metric_outranks_the_metric_only_rules():
    """A typo'd metric is C++'s error, not an occasion to advise about p."""
    mols = _mols(["CCO", "CCC", "c1ccccc1"])
    with pytest.raises(RuntimeError, match="Unknown metric"):
        oecluster.pdist(mols, "fingerprint", metric="tanimotoo", p=3.0)
    with pytest.raises(RuntimeError, match="Unknown metric"):
        oecluster.pdist(
            mols, "fingerprint", metric="nonsense", tversky_alpha=0.9)


def test_a_recognized_metric_still_gets_the_metric_only_rules():
    """The control: deferring to C++ must not cost real advice."""
    mols = _mols(["CCO", "CCC", "c1ccccc1"])
    with pytest.raises(TypeError, match="p does not apply to metric='tanimoto'"):
        oecluster.pdist(mols, "fingerprint", metric="tanimoto", p=3.0)


def test_an_unknown_storage_outranks_the_pdist_tversky_guard():
    """The Tversky guard runs after the rejector and needs its own guard."""
    mols = _mols(["CCO", "CCC", "c1ccccc1"])
    with pytest.raises(RuntimeError, match="Unknown fingerprint storage"):
        oecluster.pdist(
            mols, "fingerprint", storage="nonsense", metric="tversky",
            tversky_alpha=0.9, tversky_beta=0.1)


def test_a_recognized_storage_still_gets_the_pdist_tversky_guard():
    """The control for the test above."""
    mols = _mols(["CCO", "CCC", "c1ccccc1"])
    with pytest.raises(ValueError, match="pdist requires a symmetric metric"):
        oecluster.pdist(
            mols, "fingerprint", storage="binary", metric="tversky",
            tversky_alpha=0.9, tversky_beta=0.1)


def test_a_descriptor_only_metric_is_unrecognized_on_the_fingerprint_surface():
    """seuclidean exists in C++ but not on this surface, so Python defers.

    The authoritative message is the surface rejection rather than ``Unknown
    metric``: ``resolve_metric`` finds the row and then rejects it for the
    fingerprint surface (``MetricTable.cpp:154``). Either way it is C++'s to
    report, and the point of the test is that no advice about ``p`` pre-empts
    it.
    """
    mols = _mols(["CCO", "CCC", "c1ccccc1"])
    with pytest.raises(RuntimeError, match="descriptor-space metric"):
        oecluster.pdist(mols, "fingerprint", metric="seuclidean", p=3.0)


def test_an_out_of_range_tversky_weight_outranks_the_symmetry_guard():
    """Both remedies the symmetry message names fail when weights are invalid."""
    mols = _mols(["CCO", "CCC", "c1ccccc1"])
    with pytest.raises(RuntimeError, match=r"Tversky alpha and beta must be in"):
        oecluster.pdist(
            mols, "fingerprint", metric="tversky",
            tversky_alpha=5.0, tversky_beta=0.1)


def test_an_in_range_asymmetric_tversky_still_hits_the_symmetry_guard():
    """The control: standing aside on bad weights must not cost real advice."""
    mols = _mols(["CCO", "CCC", "c1ccccc1"])
    with pytest.raises(ValueError, match="pdist requires a symmetric metric"):
        oecluster.pdist(
            mols, "fingerprint", metric="tversky",
            tversky_alpha=0.9, tversky_beta=0.1)


def test_an_invalid_minkowski_p_outranks_the_numbits_rule():
    """The numbits remedies -- drop it, binary, count -- all fail on p=0.0."""
    mols = _mols(["CCO", "CCC", "c1ccccc1"])
    with pytest.raises(RuntimeError, match="Minkowski exponent p"):
        oecluster.pdist(mols, "fingerprint", metric="minkowski", p=0.0,
                        storage="sparse", numbits=4096)


def test_an_invalid_minkowski_p_outranks_the_family_only_rules():
    """Same for the family-only rules: no named remedy makes p=0.0 valid."""
    mols = _mols(["CCO", "CCC", "c1ccccc1"])
    with pytest.raises(RuntimeError, match="Minkowski exponent p"):
        oecluster.pdist(mols, "fingerprint", metric="minkowski", p=0.0,
                        fp_type="atom_pair", radius=3)


def test_an_invalid_minkowski_p_outranks_the_metric_only_rules():
    """And for the metric-only rules, where the remedy bounced to a second
    advisory rather than to the truth."""
    mols = _mols(["CCO", "CCC", "c1ccccc1"])
    with pytest.raises(RuntimeError, match="Minkowski exponent p"):
        oecluster.pdist(mols, "fingerprint", metric="minkowski", p=0.0,
                        tversky_alpha=0.3)


def test_an_out_of_range_tversky_weight_outranks_the_numbits_rule():
    """Round 7 covered the symmetry guard; the rejector needed it too."""
    mols = _mols(["CCO", "CCC", "c1ccccc1"])
    with pytest.raises(RuntimeError, match="Tversky alpha and beta must be in"):
        oecluster.pdist(mols, "fingerprint", metric="tversky",
                        tversky_alpha=5.0, tversky_beta=0.5,
                        storage="sparse", numbits=4096)


def test_a_missing_similarity_form_outranks_the_metric_only_rules():
    """similarity=True on a distance-only metric is authoritative too."""
    mols = _mols(["CCO", "CCC", "c1ccccc1"])
    with pytest.raises(RuntimeError, match="has no similarity form"):
        oecluster.pdist(mols, "fingerprint", similarity=True,
                        metric="bray_curtis", p=3.0)


def test_an_unsupported_storage_for_a_family_outranks_the_numbits_rule():
    """Deferred in round 7 as 'names a working remedy'; the authoritative
    message names a better one, so ordering by authority wins outright."""
    mols = _mols(["CCO", "CCC", "c1ccccc1"])
    with pytest.raises(RuntimeError, match="does not support storage="):
        oecluster.pdist(mols, "fingerprint", fp_type="topological_torsions",
                        storage="sparse_count", numbits=4096)


def test_a_count_discarding_metric_outranks_the_metric_only_rules():
    """The other round 7 deferral, closed the same way."""
    mols = _mols(["CCO", "CCC", "c1ccccc1"])
    with pytest.raises(RuntimeError, match="discards the counts"):
        oecluster.pdist(mols, "fingerprint", storage="count",
                        metric="tanimoto", p=3.0)


def test_cdist_orders_errors_the_same_way_as_pdist():
    """The rejector is shared, so the fix has to reach both entry points."""
    mols = _mols(["CCO", "CCC", "c1ccccc1"])
    with pytest.raises(RuntimeError, match="Minkowski exponent p"):
        oecluster.cdist(mols, mols, "fingerprint", metric="minkowski", p=0.0,
                        storage="sparse", numbits=4096)


def test_an_unknown_family_stands_aside_from_the_family_only_rules():
    """A family C++ accepts but the alias table lacks must lose its advisory
    rules, not misfire. Called directly: reaching this through the public API
    would need a family C++ knows and this module does not, which by
    construction does not exist today."""
    _comparisons.reject_inapplicable_fingerprint_kwargs(
        {"radius"}, fp_type="a_family_cpp_knows", storage=None, metric=None)


def test_the_advisory_rules_still_fire_when_cpp_accepts():
    """The control for this whole round: deferring to C++ must not cost the
    advice that C++ cannot give, since it accepts these silently."""
    mols = _mols(["CCO", "CCC", "c1ccccc1"])
    with pytest.raises(TypeError, match="numbits does not apply"):
        oecluster.pdist(mols, "fingerprint", storage="sparse", numbits=4096)
    with pytest.raises(TypeError, match="radius does not apply"):
        oecluster.pdist(mols, "fingerprint", fp_type="atom_pair", radius=3)
    with pytest.raises(TypeError, match="tversky_alpha does not apply"):
        oecluster.pdist(mols, "fingerprint", metric="tanimoto",
                        tversky_alpha=0.3)


def test_a_sparse_numbits_of_zero_reports_the_inapplicable_kwarg():
    """Morgan's sparse generator must not judge a numbits it never reads.

    Its own message names the only remedy that fails: a positive numbits is
    rejected here, and dropping numbits is what actually works.
    """
    mols = _mols(["CCO", "CCC", "c1ccccc1"])
    with pytest.raises(TypeError, match=r"numbits does not apply to storage"):
        oecluster.pdist(mols, "fingerprint", storage="sparse", numbits=0)


def test_every_remedy_the_sparse_numbits_message_names_is_reachable():
    """The message offers three ways out, and one of them is conditional.

    ``storage='count'`` used to be offered flat, which sent the caller into a
    second refusal: counted storage rejects the default bit-set metric. The
    message now says so, and this asserts both halves of that -- the bare
    switch still fails, and the switch the message describes works.
    """
    mols = _mols(["CCO", "CCC", "c1ccccc1"])
    with pytest.raises(TypeError, match=r"numeric metric such as 'manhattan'"):
        oecluster.pdist(mols, "fingerprint", storage="sparse", numbits=4096)

    oecluster.pdist(mols, "fingerprint", storage="sparse")
    oecluster.pdist(mols, "fingerprint", storage="binary", numbits=4096)
    with pytest.raises(RuntimeError, match=r"is a bit-set metric"):
        oecluster.pdist(mols, "fingerprint", storage="count", numbits=4096)
    oecluster.pdist(mols, "fingerprint", storage="count", numbits=4096,
                    metric="manhattan")


def test_a_sparse_count_numbits_of_zero_reports_the_inapplicable_kwarg():
    """The same hole exists under sparse_count storage."""
    mols = _mols(["CCO", "CCC", "c1ccccc1"])
    with pytest.raises(TypeError, match=r"numbits does not apply to storage"):
        oecluster.pdist(
            mols, "fingerprint", storage="sparse_count",
            metric="bray_curtis", numbits=0)


def test_cdist_reports_the_sparse_numbits_kwarg_the_same_way():
    """Both public surfaces build through the same registry."""
    mols = _mols(["CCO", "CCC", "c1ccccc1"])
    with pytest.raises(TypeError, match=r"numbits does not apply to storage"):
        oecluster.cdist(mols, mols, "fingerprint", storage="sparse", numbits=0)


def test_a_binary_numbits_of_zero_still_reports_the_cpp_bound():
    """The control: where numbits IS read, the C++ bound is authoritative.

    Resetting the field must be confined to the storages that ignore it.
    """
    mols = _mols(["CCO", "CCC", "c1ccccc1"])
    with pytest.raises(RuntimeError, match=r"num_bits must be greater than"):
        oecluster.pdist(mols, "fingerprint", storage="binary", numbits=0)


def test_an_unknown_metric_outranks_the_sparse_numbits_rule():
    """Resetting numbits must not cost the constructor its first word.

    This and the three tests below hold at HEAD already. They exist because the
    obvious alternative fix -- raising the numbits advisory before construction
    -- breaks every one of them, and nothing else in the suite would notice.
    """
    mols = _mols(["CCO", "CCC", "c1ccccc1"])
    with pytest.raises(RuntimeError, match=r"Unknown metric"):
        oecluster.pdist(
            mols, "fingerprint", storage="sparse", numbits=4096,
            metric="not_a_real_metric")


def test_an_unknown_family_outranks_the_sparse_numbits_rule():
    mols = _mols(["CCO", "CCC", "c1ccccc1"])
    with pytest.raises(RuntimeError, match=r"Unknown OEFP fingerprint type"):
        oecluster.pdist(
            mols, "fingerprint", storage="sparse", numbits=4096,
            fp_type="not_a_real_family")


def test_an_unknown_storage_outranks_the_sparse_numbits_rule():
    mols = _mols(["CCO", "CCC", "c1ccccc1"])
    with pytest.raises(RuntimeError, match=r"Unknown fingerprint storage"):
        oecluster.pdist(
            mols, "fingerprint", storage="not_a_real_storage", numbits=4096)


def test_an_invalid_minkowski_p_outranks_the_sparse_numbits_rule():
    mols = _mols(["CCO", "CCC", "c1ccccc1"])
    with pytest.raises(RuntimeError, match=r"Minkowski exponent p must be"):
        oecluster.pdist(
            mols, "fingerprint", storage="sparse", numbits=4096,
            metric="minkowski", p=0.0)


def test_a_reset_numbits_never_reaches_a_successful_call():
    """Invariant 2, asserted behaviorally rather than by reading the code.

    ``numbits_is_inapplicable`` guards both the reset and the rejection, so no
    value it discards can survive into a call that succeeds. If those two sites
    ever drift apart, some named numbits below will silently succeed.
    """
    mols = _mols(["CCO", "CCC", "c1ccccc1"])
    # sparse_count rejects the default bit-set metric whatever numbits holds,
    # which is an unrelated rule. bray_curtis satisfies both storages and keeps
    # this test measuring only the numbits invariant.
    for storage in _comparisons._SPARSE_STORAGES:
        for numbits in (0, 1, 2048, 4096):
            with pytest.raises(TypeError, match=r"numbits does not apply"):
                oecluster.pdist(
                    mols, "fingerprint", storage=storage, numbits=numbits,
                    metric="bray_curtis")
        # Dropping numbits is the remedy the message names, and it must work.
        oecluster.pdist(
            mols, "fingerprint", storage=storage, metric="bray_curtis")


def test_a_non_string_selector_names_the_argument():
    """A bare AttributeError blames ``.lower()``, which no caller wrote.

    Reached only by calling the rule functions directly: on the pdist and
    cdist paths the SWIG setters refuse a non-string first.
    """
    with pytest.raises(TypeError, match=r"fp_type must be a string or None"):
        _comparisons.canonical_fingerprint_family(2)
    with pytest.raises(TypeError, match=r"storage must be a string or None"):
        _comparisons.numbits_is_inapplicable({"numbits"}, 3)
    with pytest.raises(TypeError, match=r"metric must be a string or None"):
        _comparisons.reject_inapplicable_fingerprint_kwargs(
            {"p"}, fp_type=None, storage=None, metric=["minkowski"])


def test_the_selector_guard_still_lets_every_string_through():
    """Refusing a non-string must not start refusing a chosen string.

    The empty string is the case that decides it: it is not ``None``, so it is
    a selection, and it has to reach the rules folded rather than be turned
    back into the default or refused outright.
    """
    assert _comparisons._default_selector(None, "binary", "storage") == "binary"
    assert _comparisons._default_selector("", "binary", "storage") == ""
    assert _comparisons._default_selector("SPARSE", "binary", "storage") == "sparse"

    # A selector this module cannot name loses its advisory rules rather than
    # raising, and the empty string is one such selector.
    assert _comparisons.canonical_fingerprint_family("") is None
    assert not _comparisons.numbits_is_inapplicable({"numbits"}, "")
