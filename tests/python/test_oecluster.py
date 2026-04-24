"""Template-aligned tests for oecluster Python bindings."""

import pytest


class TestVersion:
    """Test version information is exposed correctly."""

    def test_version_available(self):
        """Verify __version__ and __version_info__ are accessible."""
        import oecluster
        assert hasattr(oecluster, '__version__')
        assert hasattr(oecluster, '__version_info__')
        assert oecluster.__version__ == "3.1.7"
        assert oecluster.__version_info__ == (3, 1, 7)


class TestPackageImports:
    """Test that core package symbols import correctly."""

    def test_core_symbols(self):
        """Verify core classes and functions are available."""
        import oecluster
        assert hasattr(oecluster, 'pdist')
        assert hasattr(oecluster, 'DenseStorage')
        assert hasattr(oecluster, 'DistanceMatrix')
        assert hasattr(oecluster, 'FingerprintMetric')
        assert hasattr(oecluster, 'ROCSMetric')
        assert hasattr(oecluster, 'SuperposeMetric')


class TestNativeMoleculePassing:
    """Verify cross-SWIG-runtime typemaps work for molecule arguments."""

    def test_fingerprint_metric_accepts_oemolbase(self, aspirin_mol, ethanol_mol):
        """Verify FingerprintMetric accepts OEGraphMol from openeye.oechem.

        This test confirms the cross-SWIG-runtime typemap works: molecules
        created by openeye.oechem (SWIG v4) are accepted by our module
        (SWIG v5) without needing SMILES serialization.
        """
        import oecluster
        dist = oecluster.pdist([aspirin_mol, ethanol_mol], "fingerprint")
        assert dist.num_items == 2
        assert dist.num_pairs == 1
