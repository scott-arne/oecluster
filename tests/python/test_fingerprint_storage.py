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


def test_the_atom_pair_window_defaults_to_the_oefp_values():
    """The old default truncated the window to 0-2 bonds."""
    dist = oecluster.pdist(_mols(), "fingerprint", fp_type="atom_pair")
    narrow = oecluster.pdist(_mols(), "fingerprint", fp_type="atom_pair",
                             min_distance=0, max_distance=2)
    assert not np.allclose(dist.condensed, narrow.condensed)


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
    comparison = oecluster.FingerprintComparison(
        _mols(), fp_type="topological_torsions", storage="count",
        metric="bray_curtis", torsion_atom_count=4)
    dist = oecluster.pdist(_mols(), comparison)
    assert dist.num_samples == 6


def test_the_factory_class_enforces_the_same_explicitness_rules():
    """The prebuilt-object path must not be the way around the rules."""
    with pytest.raises(TypeError, match="radius"):
        oecluster.FingerprintComparison(_mols(), max_distance=2)
    with pytest.raises(TypeError, match="numbits does not apply"):
        oecluster.FingerprintComparison(_mols(), storage="sparse",
                                        numbits=4096)
    with pytest.raises(TypeError, match="p does not apply"):
        oecluster.FingerprintComparison(_mols(), p=3.0)
