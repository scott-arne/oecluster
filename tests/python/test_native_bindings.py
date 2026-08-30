"""The SWIG surface the Python layer is built on.

These assertions are deliberately about names and shapes rather than behavior:
they fail loudly when an interface-file edit silently drops a symbol, which is
otherwise only visible as an AttributeError deep inside a builder.
"""

import pytest


@pytest.fixture
def native():
    from oecluster import oecluster as _oecluster

    return _oecluster


def test_capability_enum_is_exposed(native):
    assert native.Capability_Unknown != native.Capability_No
    assert native.Capability_No != native.Capability_Yes
    assert native.Capability_Yes != native.Capability_Unknown


def test_data_integrity_enum_is_exposed(native):
    values = {
        native.DataIntegrity_Complete,
        native.DataIntegrity_NaNPresent,
        native.DataIntegrity_SubsetScored,
    }
    assert len(values) == 3


def test_gate_facts_defaults(native):
    facts = native.GateFacts()
    assert facts.is_distance == native.Capability_Unknown
    assert facts.zero_self == native.Capability_Unknown
    assert facts.triangle == native.Capability_Unknown
    assert facts.data_integrity == native.DataIntegrity_Complete


def _benzene_series():
    from openeye import oechem

    mols = []
    for smi in ("c1ccccc1", "c1ccc(O)cc1", "CCCCCCCC", "c1ccncc1"):
        mol = oechem.OEGraphMol()
        oechem.OESmilesToMol(mol, smi)
        mols.append(mol)
    return mols


def test_fingerprint_comparison_reports_facts(native):
    comparison = native.FingerprintComparison(_benzene_series(), native.FingerprintOptions())
    facts = comparison.Facts()
    assert facts.is_distance == native.Capability_Yes
    assert facts.zero_self == native.Capability_Yes
    assert facts.triangle == native.Capability_Yes
    assert facts.data_integrity == native.DataIntegrity_Complete


def test_descriptor_comparison_is_exposed(native):
    comparison = native.DescriptorComparison(_benzene_series(), native.DescriptorOptions())
    assert comparison.ComparisonName() == "descriptor"
    assert len(list(comparison.Columns())) > 0
    assert len(list(comparison.Variances())) == len(list(comparison.Columns()))
    assert comparison.Facts().data_integrity == native.DataIntegrity_Complete


def test_descriptor_excluded_indices_is_exposed(native):
    excluded = native.descriptor_excluded_indices(_benzene_series(), native.DescriptorOptions())
    assert list(excluded) == []


def test_descriptor_statistics_is_exposed(native):
    stats = native.descriptor_statistics(_benzene_series(), native.DescriptorStatisticsOptions())
    assert stats.num_rows == 4
    assert len(list(stats.columns)) == len(list(stats.variance))
    assert len(list(stats.dropped_columns)) == len(list(stats.dropped_reasons))
    assert list(stats.inverse_covariance) == []


def test_descriptor_statistics_options_accept_lists(native):
    options = native.DescriptorStatisticsOptions()
    options.sources = native.StringVector(["openeye"])
    options.inverse_covariance = True
    stats = native.descriptor_statistics(_benzene_series(), options)
    k = len(list(stats.columns))
    assert len(list(stats.inverse_covariance)) == k * k


def test_rmsd_comparison_is_exposed(native):
    from openeye import oechem

    mols = []
    for shift in range(3):
        mol = oechem.OEMol()
        oechem.OESmilesToMol(mol, "c1ccccc1")
        oechem.OEGenerate2DCoordinates(mol)
        for atom in mol.GetAtoms():
            x, y, z = mol.GetCoords(atom)
            mol.SetCoords(atom, (x + shift, y, z))
        mols.append(mol)

    comparison = native.RMSDComparison(mols, native.RMSDOptions())
    assert comparison.ComparisonName() == "rmsd"
    assert comparison.Size() == 3
    assert comparison.Compare(0, 0) == pytest.approx(0.0, abs=1e-6)
    assert comparison.Facts().zero_self == native.Capability_Yes
