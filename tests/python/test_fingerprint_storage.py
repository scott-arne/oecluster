import numpy as np
import oecluster
import pytest
from openeye import oechem

SMILES = ["CCO", "CCC", "CCCC", "c1ccccc1", "CCN", "CCCCCCCCCCCCCCCC"]

FAMILIES = ["morgan", "atom_pair", "topological_atom_pair",
            "topological_torsions"]
STORAGES = ["binary", "count", "sparse", "sparse_count"]

# topological_torsions has no sparse-count batch type in OEFP.
UNSUPPORTED = {("topological_torsions", "sparse_count")}


def _mols():
    mols = []
    for idx, smi in enumerate(SMILES):
        mol = oechem.OEGraphMol()
        oechem.OESmilesToMol(mol, smi)
        mol.SetTitle(f"mol{idx}")
        mols.append(mol)
    return mols


@pytest.mark.parametrize("family", FAMILIES)
@pytest.mark.parametrize("storage", STORAGES)
def test_the_family_by_storage_grid(family, storage):
    metric = "tanimoto" if storage in ("binary", "sparse") else "bray_curtis"
    if (family, storage) in UNSUPPORTED:
        with pytest.raises(RuntimeError):
            oecluster.pdist(_mols(), "fingerprint", fp_type=family,
                            storage=storage, metric=metric)
        return
    dist = oecluster.pdist(_mols(), "fingerprint", fp_type=family,
                           storage=storage, metric=metric)
    assert dist.num_samples == 6
    assert np.all(np.isfinite(dist.condensed))


def test_a_boolean_metric_is_rejected_on_count_storage():
    with pytest.raises(RuntimeError, match="count"):
        oecluster.pdist(_mols(), "fingerprint", storage="count",
                        metric="tanimoto")


def test_count_storage_separates_homologs_that_binary_cannot():
    """C16 vs C30: binary Morgan saturates, count fingerprints do not."""
    mols = []
    for idx, smi in enumerate(["C" * 16, "C" * 30]):
        mol = oechem.OEGraphMol()
        oechem.OESmilesToMol(mol, smi)
        mol.SetTitle(f"c{idx}")
        mols.append(mol)

    binary = oecluster.pdist(mols, "fingerprint", storage="binary")
    counted = oecluster.pdist(mols, "fingerprint", storage="count",
                              metric="bray_curtis")
    assert binary.condensed[0] == pytest.approx(0.0, abs=1e-12)
    assert counted.condensed[0] > 0.1


def test_the_declared_atom_pair_window_is_1_to_30_and_output_honours_it():
    """The declared window is 1-30, and the output vector honours it.

    Two levels, because neither pins what the other does.
    ``FingerprintOptions()`` carries the declared defaults, so asserting on it
    pins both bounds exactly -- and it is the struct the keyword path uses:
    ``_build_fingerprint`` constructs a default one and overwrites only the
    fields the caller named. The condensed vector then pins that the pipeline
    honours the struct rather than reaching some other window, which the
    struct assertion on its own does not, and that the window in force is not
    the 0-2 the 4.x default truncated to.

    A C31 chain is added to the fixture for this test alone, and it is what
    makes the upper bound observable in output: the six molecules the rest of
    the file uses top out at a graph distance of 15, so every window from 1-15
    to 1-30 -- and the old default's replacement, whatever it were --
    reproduces the same vector over them. With the chain in the vector has
    twenty-one entries, and 1-29 moves five of them.

    What the output vector resolves, measured by moving each bound in turn:

    * ``max_distance`` is pinned exactly at 30. 29 moves five entries, and 31
      is refused outright by OEFP ("max_distance must be smaller than 31"), so
      there is no larger value the default could be.
    * ``min_distance`` is pinned below 2 -- 2-30 moves fifteen of the
      twenty-one -- but 0 and 1 are indistinguishable here and, as far as this
      fixture can tell, everywhere: an atom pair at graph distance 0 is an atom
      with itself, which OEFP does not enumerate. 0-30 reproduces the default
      vector exactly, so a default of 0 would pass the vector assertions
      unnoticed. Only the declared value separates the two, which is what the
      ``FingerprintOptions`` assertions are for.
    """
    declared = oecluster.FingerprintOptions()
    assert declared.min_distance == 1
    assert declared.max_distance == 30

    mols = _mols()
    chain = oechem.OEGraphMol()
    oechem.OESmilesToMol(chain, "C" * 31)
    chain.SetTitle("c31")
    mols.append(chain)

    def window(**kwargs):
        return oecluster.pdist(mols, "fingerprint", fp_type="atom_pair",
                               **kwargs).condensed

    default = window()
    np.testing.assert_array_equal(default,
                                  window(min_distance=1, max_distance=30))
    assert not np.allclose(default, window(min_distance=1, max_distance=29))
    assert not np.allclose(default, window(min_distance=2, max_distance=30))
    assert not np.allclose(default, window(min_distance=0, max_distance=2))
    with pytest.raises(RuntimeError, match="smaller than 31"):
        window(min_distance=1, max_distance=31)


def test_topological_atom_pair_is_an_alias_of_atom_pair():
    alias = oecluster.pdist(_mols(), "fingerprint",
                            fp_type="topological_atom_pair")
    base = oecluster.pdist(_mols(), "fingerprint", fp_type="atom_pair")
    np.testing.assert_allclose(alias.condensed, base.condensed)


def test_distance_atom_pair_is_rejected_as_unimplemented():
    with pytest.raises(RuntimeError, match="distance_atom_pair"):
        oecluster.pdist(_mols(), "fingerprint", fp_type="distance_atom_pair")


def test_an_unknown_family_is_rejected():
    with pytest.raises(RuntimeError):
        oecluster.pdist(_mols(), "fingerprint", fp_type="maccs")


def test_haversine_is_rejected_on_fingerprints():
    with pytest.raises(RuntimeError, match="haversine"):
        oecluster.pdist(_mols(), "fingerprint", metric="haversine")


def test_radius_controls_the_morgan_radius():
    small = oecluster.pdist(_mols(), "fingerprint", radius=1)
    large = oecluster.pdist(_mols(), "fingerprint", radius=4)
    assert not np.allclose(small.condensed, large.condensed)


def test_max_distance_no_longer_aliases_the_morgan_radius():
    """The 5.0.0 break: max_distance is the atom-pair window only.

    4.x silently routed ``max_distance`` into Morgan's radius. Ignoring it
    instead would be just as silent, so the keyword surface refuses it and
    names the replacement.
    """
    with pytest.raises(TypeError, match="radius"):
        oecluster.pdist(_mols(), "fingerprint", max_distance=2)


def test_the_factory_class_accepts_the_new_options():
    """Every option named here reaches the comparison the factory builds.

    ``dist.num_samples == 6`` was the old assertion, and it holds just as well
    for ``FingerprintComparison(_mols())`` with no options at all, so it could
    not tell a forwarded option from a discarded one. Equality with the keyword
    path is what pins the forwarding.

    ``torsion_atom_count=3`` rather than the default 4 is deliberate. At 4 the
    option is indistinguishable from leaving it off -- the condensed vector is
    the same either way -- so the equality above would survive the factory
    dropping it. At 3 it moves four entries. The probes below establish that
    the other three named options are live on this fixture too, so the equality
    is over a configuration none of them is inert in.
    """
    options = {'fp_type': "topological_torsions", 'storage': "count",
               'metric': "bray_curtis", 'torsion_atom_count': 3}
    comparison = oecluster.FingerprintComparison(_mols(), **options)
    factory = oecluster.pdist(_mols(), comparison).condensed
    np.testing.assert_array_equal(
        factory, oecluster.pdist(_mols(), "fingerprint", **options).condensed)

    # One probe per named option, dropping or moving only that one.
    def without(key):
        return {k: v for k, v in options.items() if k != key}

    with pytest.raises(TypeError, match="torsion_atom_count does not apply"):
        oecluster.pdist(_mols(), "fingerprint", **without("fp_type"))
    with pytest.raises(RuntimeError, match="bit-set metric"):
        oecluster.pdist(_mols(), "fingerprint", **without("metric"))
    assert not np.allclose(factory, oecluster.pdist(
        _mols(), "fingerprint", **{**options, 'storage': "binary"}).condensed)
    assert not np.allclose(factory, oecluster.pdist(
        _mols(), "fingerprint",
        **{**options, 'torsion_atom_count': 4}).condensed)


def test_the_factory_class_enforces_the_same_explicitness_rules():
    """The prebuilt-object path must not be the way around the rules.

    All eight advisory rules, so that no single one of them can be lost to the
    factory path unnoticed: ``numbits`` against sparse storage, the four
    family-only options against a family that ignores them, and the three
    metric-only options against a metric that ignores them.
    """
    with pytest.raises(TypeError,
                       match="max_distance does not apply.*Use radius instead"):
        oecluster.FingerprintComparison(_mols(), max_distance=2)
    with pytest.raises(TypeError, match="min_distance does not apply"):
        oecluster.FingerprintComparison(_mols(), min_distance=2)
    with pytest.raises(TypeError, match="radius does not apply"):
        oecluster.FingerprintComparison(_mols(), fp_type="atom_pair", radius=3)
    with pytest.raises(TypeError, match="torsion_atom_count does not apply"):
        oecluster.FingerprintComparison(_mols(), torsion_atom_count=4)
    with pytest.raises(TypeError, match="numbits does not apply"):
        oecluster.FingerprintComparison(_mols(), storage="sparse",
                                        numbits=4096)
    with pytest.raises(TypeError, match="p does not apply"):
        oecluster.FingerprintComparison(_mols(), p=3.0)
    with pytest.raises(TypeError, match="tversky_alpha does not apply"):
        oecluster.FingerprintComparison(_mols(), tversky_alpha=0.5)
    with pytest.raises(TypeError, match="tversky_beta does not apply"):
        oecluster.FingerprintComparison(_mols(), tversky_beta=0.5)


@pytest.mark.parametrize("family", FAMILIES)
def test_use_chirality_is_accepted_by_every_family(family):
    """Not an over-refusal: ``use_chirality`` is the one option every family reads.

    A rule table that swept it in with the family-only options would refuse a
    call the C++ layer honours, which the explicitness rules exist to prevent
    rather than to cause.
    """
    oecluster.FingerprintComparison(_mols(), fp_type=family,
                                    use_chirality=True)
