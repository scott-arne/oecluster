/**
 * @file test_murcko_scaffold.cpp
 * @brief Tests for Bemis-Murcko scaffold extraction and clustering.
 */

#include <cstddef>
#include <optional>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include <oechem.h>

#include "oecluster/clustering/MurckoScaffold.h"

#include "../../src/clustering/MurckoKernels.h"

namespace {

using OECluster::ScaffoldType;
using OECluster::detail::scaffold_of;

/// Extract a scaffold from a SMILES string, reporting a failure as a string
/// the assertion messages can show rather than as an empty optional.
std::string scaffold_from_smiles(const char* smiles, ScaffoldType type) {
    OEChem::OEGraphMol mol;
    OEChem::OESmilesToMol(mol, smiles);
    return scaffold_of(mol, type).value_or(std::string("<extraction failed>"));
}

std::string framework(const char* smiles) {
    return scaffold_from_smiles(smiles, ScaffoldType::Framework);
}

std::string generic(const char* smiles) {
    return scaffold_from_smiles(smiles, ScaffoldType::Generic);
}

/// Heavy-atom count of a SMILES string. Re-parsing an aromatic framework can
/// emit a kekulization warning on stderr; the atom count is still correct and
/// the warning is not a failure.
size_t atom_count(const std::string& smiles) {
    OEChem::OEGraphMol mol;
    OEChem::OESmilesToMol(mol, smiles.c_str());
    return mol.NumAtoms();
}

TEST(MurckoExtractionTest, RemovesSidechains) {
    EXPECT_EQ(framework("Cc1ccccc1"), framework("c1ccccc1"));
    EXPECT_EQ(framework("CC(=O)Oc1ccccc1C(=O)O"), framework("c1ccccc1"));
}

TEST(MurckoExtractionTest, PinsTheBenzeneFrameworkLiteral) {
    // One exact literal, so a change in OEChem's canonicalization is caught
    // rather than silently absorbed by the relational assertions around it.
    EXPECT_EQ(framework("c1ccccc1"), "c1ccccc1");
}

TEST(MurckoExtractionTest, KeepsTheLinkerBetweenRings) {
    EXPECT_NE(framework("c1ccc(cc1)Cc1ccccc1"), framework("c1ccc(cc1)-c1ccccc1"));
    EXPECT_EQ(framework("c1ccc(cc1)Cc1ccccc1"), "c1ccc(cc1)Cc2ccccc2");
    EXPECT_EQ(framework("c1ccc(cc1)-c1ccccc1"), "c1ccc(cc1)c2ccccc2");
}

TEST(MurckoExtractionTest, TheDiphenylmethaneLinkerIsExactlyOneAtom) {
    // The structural form of the assertion above: twelve ring atoms plus the
    // single methylene that joins them, which is what "the linker is kept"
    // means without appealing to the canonical string.
    OEChem::OEGraphMol scaffold;
    OEChem::OESmilesToMol(scaffold, framework("c1ccc(cc1)Cc1ccccc1").c_str());
    OEChem::OEFindRingAtomsAndBonds(scaffold);

    EXPECT_EQ(scaffold.NumAtoms(), 13u);
    EXPECT_EQ(scaffold.NumBonds(), 14u);
    size_t non_ring_atoms = 0;
    for (OESystem::OEIter<OEChem::OEAtomBase> atom = scaffold.GetAtoms(); atom;
         ++atom) {
        if (!atom->IsInRing()) {
            ++non_ring_atoms;
        }
    }
    EXPECT_EQ(non_ring_atoms, 1u);
}

TEST(MurckoExtractionTest, DropsNMethylsFromCaffeine) {
    const char* caffeine = "Cn1cnc2c1c(=O)n(C)c(=O)n2C";
    EXPECT_EQ(framework(caffeine), "c1c2c([nH]c[nH]1)nc[nH]2");
    EXPECT_LT(atom_count(framework(caffeine)), atom_count(caffeine));
}

TEST(MurckoExtractionTest, PreservesElementIdentityAtFramework) {
    EXPECT_NE(framework("c1ccccc1"), framework("c1ccncc1"));
    EXPECT_EQ(framework("c1ccncc1"), "c1ccncc1");
}

TEST(MurckoExtractionTest, GenericReductionConverges) {
    // The binding contract of the Generic level: benzene and pyridine must
    // become one string. A hand-rolled atom-type rewrite does not achieve this.
    EXPECT_EQ(generic("c1ccccc1"), generic("c1ccncc1"));
    EXPECT_EQ(generic("c1ccccc1"), "C1CCCCC1");
}

TEST(MurckoExtractionTest, GenericErasesBondOrder) {
    // The lactam pair from the spec: its exocyclic carbonyl is a sidechain, so
    // what converges here is the ring heteroatom.
    EXPECT_EQ(generic("O=C1CCCCN1"), generic("C1CCCCC1"));
    // A ring double bond, which the pair above does not actually exercise.
    EXPECT_NE(framework("FC1=CCCCC1"), framework("FC1CCCCC1"));
    EXPECT_EQ(generic("FC1=CCCCC1"), generic("FC1CCCCC1"));
    EXPECT_EQ(generic("FC1=CCCCC1"), "C1CCCCC1");
}

TEST(MurckoExtractionTest, StripsStereochemistry) {
    // Decalin keeps both stereocenters inside the framework, so this pair
    // actually tests stereo stripping rather than sidechain removal.
    EXPECT_EQ(framework("C1CC[C@H]2CCCC[C@@H]2C1"),
              framework("C1CC[C@H]2CCCC[C@H]2C1"));
    EXPECT_EQ(framework("C1CC[C@H]2CCCC[C@@H]2C1"), "C1CCC2CCCCC2C1");
    // A plain enantiomer pair, whose stereocenter is in the sidechain.
    EXPECT_EQ(framework("C[C@H](N)c1ccccc1"), framework("C[C@@H](N)c1ccccc1"));
}

TEST(MurckoExtractionTest, IgnoresExplicitHydrogens) {
    OEChem::OEGraphMol implicit_h;
    OEChem::OESmilesToMol(implicit_h, "Cc1ccccc1");
    OEChem::OEGraphMol explicit_h;
    OEChem::OESmilesToMol(explicit_h, "Cc1ccccc1");
    OEChem::OEAddExplicitHydrogens(explicit_h);

    EXPECT_EQ(scaffold_of(implicit_h, ScaffoldType::Framework),
              scaffold_of(explicit_h, ScaffoldType::Framework));
    EXPECT_EQ(scaffold_of(implicit_h, ScaffoldType::Generic),
              scaffold_of(explicit_h, ScaffoldType::Generic));
}

TEST(MurckoExtractionTest, JoinsDisconnectedComponentsIntoOneString) {
    EXPECT_EQ(framework("c1ccccc1C(=O)O.c1ccncc1"), "c1ccccc1.c1ccncc1");
    EXPECT_EQ(generic("c1ccccc1C(=O)O.c1ccncc1"), "C1CCCCC1.C1CCCCC1");
}

TEST(MurckoExtractionTest, AcyclicMoleculesYieldTheEmptyString) {
    EXPECT_EQ(framework("CCCCCC"), "");
    EXPECT_EQ(framework("C"), "");
    EXPECT_EQ(framework(""), "");
    // The acyclic rule survives the Generic reduction unchanged.
    EXPECT_EQ(generic("CCCCCC"), "");
}

TEST(MurckoExtractionTest, LeavesTheInputMoleculeUnmodified) {
    // Ring-membership flags are checked alongside the SMILES because
    // OEFindRingAtomsAndBonds run on the caller's molecule instead of the copy
    // would change flags without changing either the SMILES or the counts.
    OEChem::OEGraphMol mol;
    OEChem::OESmilesToMol(mol, "CC(=O)Oc1ccccc1C(=O)O");

    std::string before_smiles;
    OEChem::OECreateCanSmiString(before_smiles, mol);
    const unsigned int before_atoms = mol.NumAtoms();
    const unsigned int before_bonds = mol.NumBonds();
    std::vector<bool> before_ring_flags;
    for (OESystem::OEIter<OEChem::OEBondBase> bond = mol.GetBonds(); bond; ++bond) {
        before_ring_flags.push_back(bond->IsInRing());
    }

    ASSERT_TRUE(scaffold_of(mol, ScaffoldType::Framework).has_value());
    ASSERT_TRUE(scaffold_of(mol, ScaffoldType::Generic).has_value());

    std::string after_smiles;
    OEChem::OECreateCanSmiString(after_smiles, mol);
    EXPECT_EQ(after_smiles, before_smiles);
    EXPECT_EQ(mol.NumAtoms(), before_atoms);
    EXPECT_EQ(mol.NumBonds(), before_bonds);
    std::vector<bool> after_ring_flags;
    for (OESystem::OEIter<OEChem::OEBondBase> bond = mol.GetBonds(); bond; ++bond) {
        after_ring_flags.push_back(bond->IsInRing());
    }
    EXPECT_EQ(after_ring_flags, before_ring_flags);
}

}  // namespace
