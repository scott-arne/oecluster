"""Tests for Murcko scaffold assignment and clustering."""

import oecluster
import pytest

oechem = pytest.importorskip("openeye.oechem")


def mols(*smiles):
    """Build a list of OEGraphMol from SMILES strings."""
    built = []
    for value in smiles:
        mol = oechem.OEGraphMol()
        oechem.OESmilesToMol(mol, value)
        built.append(mol)
    return built


BENZENE = "c1ccccc1"
TOLUENE = "Cc1ccccc1"
PYRIDINE = "c1ccncc1"
HEXANE = "CCCCCC"
DIPHENYL_SULFONE = "O=S(=O)(c1ccccc1)c1ccccc1"
DIPHENYL_SULFOXIDE = "O=S(c1ccccc1)c1ccccc1"
DIPHENYL_SULFIDE = "S(c1ccccc1)c1ccccc1"
AZONIASPIRO = "C1CC[N+]2(CCCCC2)CC1"


class TestScaffoldKeyword:
    def test_framework_is_the_default(self):
        assert oecluster.murcko_scaffolds(mols(BENZENE)) == ["c1ccccc1"]

    def test_generic_reduces_to_topology(self):
        scaffolds = oecluster.murcko_scaffolds(
            mols(BENZENE, PYRIDINE), scaffold="generic")
        assert scaffolds == ["C1CCCCC1", "C1CCCCC1"]

    def test_the_keyword_is_case_insensitive(self):
        assert (oecluster.murcko_scaffolds(mols(PYRIDINE), scaffold="GENERIC")
                == oecluster.murcko_scaffolds(mols(PYRIDINE), scaffold="generic"))

    def test_an_unknown_name_names_both_accepted_values(self):
        with pytest.raises(ValueError, match="framework.*generic"):
            oecluster.murcko_scaffolds(mols(BENZENE), scaffold="bogus")

    def test_a_non_string_is_a_type_error(self):
        with pytest.raises(TypeError):
            oecluster.murcko_scaffolds(mols(BENZENE), scaffold=1)


class TestScaffoldKeywordReachesTheClusterer:
    def test_generic_merges_what_framework_separates(self):
        # The contrast is the point: under a defect that drops the caller's
        # level on the way to the native call, both runs come back Framework
        # and every other test stays green.
        molecules = mols(BENZENE, PYRIDINE)
        framework = oecluster.murcko(molecules)
        assert framework.cluster_scaffolds == ("c1ccccc1", "c1ccncc1")
        generic = oecluster.murcko(molecules, scaffold="generic")
        assert generic.cluster_scaffolds == ("C1CCCCC1",)
        assert list(generic.labels) == [0, 0]


class TestHeteroatomLinkerNormalization:
    # The Python mirror of the C++ fixtures: the normalization lives in the
    # kernel, so what these pin is that the binding surfaces it unchanged.

    def test_the_sulfur_oxidation_states_share_one_scaffold(self):
        molecules = mols(DIPHENYL_SULFONE, DIPHENYL_SULFOXIDE, DIPHENYL_SULFIDE)
        assert (oecluster.murcko_scaffolds(molecules)
                == ["c1ccc(cc1)Sc2ccccc2"] * 3)

    def test_the_oxidation_states_land_in_one_cluster(self):
        result = oecluster.murcko(
            mols(DIPHENYL_SULFONE, DIPHENYL_SULFOXIDE, DIPHENYL_SULFIDE,
                 BENZENE))
        assert result.cluster_scaffolds == ("c1ccc(cc1)Sc2ccccc2", "c1ccccc1")
        assert list(result.labels) == [0, 0, 0, 1]

    def test_an_n_methylated_ring_reduces_to_the_parent(self):
        assert (oecluster.murcko_scaffolds(mols("C[n+]1ccccc1"))
                == oecluster.murcko_scaffolds(mols(PYRIDINE)))

    def test_a_ring_charge_the_cut_never_touched_survives(self):
        # The scope limit: normalizing this one would produce a five-bonded
        # neutral nitrogen.
        assert (oecluster.murcko_scaffolds(mols(AZONIASPIRO))
                == ["C1CC[N+]2(CC1)CCCCC2"])


class TestValidationMirror:
    def test_an_empty_list_is_a_value_error(self):
        # Pins the caller name as well as the condition. murcko() raises the
        # same sentence, so matching only the condition leaves the two callers
        # interchangeable -- which is exactly the mutation that survived.
        with pytest.raises(
                ValueError,
                match=r"murcko_scaffolds\(\) requires at least one molecule"):
            oecluster.murcko_scaffolds([])

    def test_an_empty_tuple_is_a_value_error_not_a_type_error(self):
        # The mirror's length check fires before the typemap would have had an
        # element to object to; "empty container of the wrong type" is
        # otherwise genuinely ambiguous between the two.
        with pytest.raises(ValueError):
            oecluster.murcko_scaffolds(())

    def test_a_non_integer_thread_count_is_a_type_error(self):
        with pytest.raises(TypeError):
            oecluster.murcko_scaffolds(mols(BENZENE), num_threads=1.5)

    def test_a_negative_thread_count_is_a_value_error(self):
        with pytest.raises(ValueError):
            oecluster.murcko_scaffolds(mols(BENZENE), num_threads=-1)

    def test_an_oversized_thread_count_is_a_value_error_not_an_overflow(self):
        # Without the mirror this reaches SWIG's size_t conversion and escapes
        # as OverflowError, which says nothing about which argument was wrong.
        with pytest.raises(ValueError, match="size_t"):
            oecluster.murcko_scaffolds(mols(BENZENE), num_threads=1 << 128)

    def test_a_non_list_is_a_type_error(self):
        with pytest.raises(TypeError):
            oecluster.murcko_scaffolds("not molecules")

    def test_an_unsized_input_names_the_function(self):
        # Without the guard this escapes as "object of type 'generator' has no
        # len()" -- the right exception type, carrying neither the function name
        # nor what it wanted. A string reaches the typemap and gets a good
        # message; a generator did not.
        with pytest.raises(
                TypeError,
                match=r"murcko_scaffolds\(\) requires a list of molecules"):
            oecluster.murcko_scaffolds(mol for mol in mols(BENZENE))

    def test_a_non_molecule_element_is_a_type_error(self):
        with pytest.raises(TypeError):
            oecluster.murcko_scaffolds([1, 2, 3])

    def test_the_scaffold_verdict_precedes_the_empty_input_verdict(self):
        with pytest.raises(ValueError, match="framework.*generic"):
            oecluster.murcko_scaffolds([], scaffold="bogus")

    def test_the_clusterer_names_itself_in_the_message(self):
        # murcko_scaffolds() does not contain the substring "murcko()", so this
        # regex fails if the wrong caller name is threaded through.
        with pytest.raises(ValueError,
                           match=r"murcko\(\) requires at least one molecule"):
            oecluster.murcko([])


class TestMurckoResult:
    def test_exposes_scaffolds_and_clusters(self):
        result = oecluster.murcko(mols(BENZENE, TOLUENE, PYRIDINE, HEXANE))

        assert result.method == "murcko"
        assert result.scaffolds == ("c1ccccc1", "c1ccccc1", "c1ccncc1", "")
        assert result.cluster_scaffolds == ("c1ccccc1", "c1ccncc1")
        assert result.num_samples == 4
        assert result.num_clusters == 2
        assert result.clusters == ((0, 1), (2,))
        assert len(result) == 2
        assert result[0] == (0, 1)
        assert [tuple(cluster) for cluster in result] == [(0, 1), (2,)]
        assert list(result.labels) == [0, 0, 1, -1]

    def test_is_the_pure_python_class_not_the_swig_type(self):
        assert oecluster.MurckoResult is not oecluster.oecluster.MurckoResult
        result = oecluster.murcko(mols(BENZENE))
        assert isinstance(result, oecluster.MurckoResult)
        assert isinstance(result, oecluster.ClusteringResult)


class TestPublication:
    def test_star_import_resolves_every_exported_name(self):
        # The assertion that would have caught exporting ScaffoldType, which
        # SWIG never binds as a class.
        for name in oecluster.__all__:
            assert hasattr(oecluster, name), name

    def test_no_scaffold_type_class_is_published(self):
        assert not hasattr(oecluster, "ScaffoldType")

    def test_the_new_names_are_exported(self):
        for name in ("murcko", "murcko_scaffolds", "MurckoResult"):
            assert name in oecluster.__all__


class TestComposition:
    def test_the_two_entry_points_agree(self):
        molecules = mols(BENZENE, TOLUENE, PYRIDINE, HEXANE)
        assert (list(oecluster.murcko(molecules).scaffolds)
                == oecluster.murcko_scaffolds(molecules))

    def test_scaffold_agreement_consumes_the_output(self):
        molecules = mols(BENZENE, TOLUENE, PYRIDINE, HEXANE)
        scaffolds = oecluster.murcko_scaffolds(molecules)
        result = oecluster.murcko(molecules)
        agreement = oecluster.scaffold_agreement(result.labels, scaffolds)
        assert agreement.completeness == pytest.approx(1.0)

    def test_an_all_acyclic_annotation_routes_through_noise_handling(self):
        molecules = mols(HEXANE, "CCO", "C")
        scaffolds = oecluster.murcko_scaffolds(molecules)
        assert scaffolds == ["", "", ""]
        result = oecluster.murcko(molecules)
        assert list(result.labels) == [-1, -1, -1]
        assert result.num_clusters == 0
        assert result.cluster_scaffolds == ()
        # Producer and consumer agree that "" is missing data, not a category.
        oecluster.scaffold_agreement(result.labels, scaffolds, noise="singletons")
