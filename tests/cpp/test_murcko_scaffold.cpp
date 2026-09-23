/**
 * @file test_murcko_scaffold.cpp
 * @brief Tests for Bemis-Murcko scaffold extraction and clustering.
 */

#include <algorithm>
#include <cstddef>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#include <gtest/gtest.h>

#include <oechem.h>
#include <oesystem.h>

#include "oecluster/Error.h"
#include "oecluster/clustering/MurckoScaffold.h"

#include "../../src/clustering/MurckoKernels.h"

namespace {

using OECluster::ScaffoldType;
using OECluster::detail::dispatch_thread_count;
using OECluster::detail::effective_thread_count;
using OECluster::detail::extract_all;
using OECluster::detail::first_failure;
using OECluster::detail::pool_is_thread_safe;
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

/// Heavy-atom count of a SMILES string, asserting that the string parses. A
/// scaffold OEChem cannot read back still yields a plausible count from the
/// partial parse, so an unchecked status hides exactly the bug this helper
/// would otherwise be measuring.
size_t atom_count(const std::string& smiles) {
    OEChem::OEGraphMol mol;
    EXPECT_TRUE(OEChem::OESmilesToMol(mol, smiles.c_str())) << smiles;
    return mol.NumAtoms();
}

/// Whether a scaffold string reparses. A scaffold the SDK that produced it
/// cannot read back is not a usable canonical SMILES, however plausible it
/// looks.
bool round_trips(const std::string& smiles) {
    OEChem::OEGraphMol mol;
    return OEChem::OESmilesToMol(mol, smiles.c_str());
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
    EXPECT_EQ(framework(caffeine), "c1[nH]c2c(n1)NCNC2");
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

TEST(MurckoExtractionTest, ScaffoldsRoundTripThroughOEChem) {
    // Aromatic flags inherited from the parent once produced a caffeine
    // framework that OEChem itself could not kekulize, and the suite passed
    // regardless because nothing checked a parse status.
    const char* inputs[] = {
        "c1ccccc1",
        "c1ccncc1",
        "c1ccc(cc1)Cc1ccccc1",
        "c1ccc(cc1)-c1ccccc1",
        "Cn1cnc2c1c(=O)n(C)c(=O)n2C",
        "C1CC[C@H]2CCCC[C@@H]2C1",
        "FC1=CCCCC1",
        "O=C1CCCCN1",
        "c1ccccc1C(=O)O.c1ccncc1",
    };
    for (const char* smiles : inputs) {
        EXPECT_TRUE(round_trips(framework(smiles))) << "framework of " << smiles;
        EXPECT_TRUE(round_trips(generic(smiles))) << "generic of " << smiles;
    }
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

/// Owns a set of molecules parsed from SMILES and hands out the pointer vector
/// the public entry points take. Pointers are collected only after every
/// molecule is allocated, so no reallocation can invalidate them.
class MolSet {
public:
    explicit MolSet(const std::vector<std::string>& smiles) {
        for (const std::string& value : smiles) {
            auto mol = std::make_unique<OEChem::OEGraphMol>();
            OEChem::OESmilesToMol(*mol, value.c_str());
            owned_.push_back(std::move(mol));
        }
        // OEGraphMol is not derived from OEMolBase -- it owns one and converts
        // to it -- so the cast is the repository's established idiom rather
        // than a pointer upcast.
        for (const auto& mol : owned_) {
            pointers_.push_back(&static_cast<OEChem::OEMolBase&>(*mol));
        }
    }

    const std::vector<OEChem::OEMolBase*>& Pointers() const { return pointers_; }

private:
    std::vector<std::unique_ptr<OEChem::OEGraphMol>> owned_;
    std::vector<OEChem::OEMolBase*> pointers_;
};

/// 200 molecules over 10 distinct frameworks: each core carries an alkyl chain
/// of length 0 to 19. Large enough to span several chunks at the internal chunk
/// size of 64, which a fixture of a dozen molecules would not.
std::vector<std::string> threading_fixture_smiles() {
    static const char* const CORES[] = {
        "c1ccccc1", "c1ccncc1", "c1ccc2ccccc2c1", "C1CCCCC1", "c1ccsc1",
        "c1cc[nH]c1", "C1CCNCC1", "c1ccc(cc1)-c1ccccc1", "C1CCOCC1", "c1cnccn1",
    };
    std::vector<std::string> smiles;
    for (const char* core : CORES) {
        for (size_t length = 0; length < 20; ++length) {
            smiles.push_back(std::string(length, 'C') + core);
        }
    }
    return smiles;
}

TEST(MurckoFailureScanTest, ReturnsTheFirstFailureInInputOrder) {
    std::vector<std::optional<std::string>> raw(10, std::string("c1ccccc1"));
    raw[7] = std::nullopt;
    raw[3] = std::nullopt;
    const std::optional<size_t> bad = first_failure(raw);
    ASSERT_TRUE(bad.has_value());
    EXPECT_EQ(*bad, 3u);
}

TEST(MurckoFailureScanTest, ReturnsNulloptWhenEveryExtractionSucceeded) {
    const std::vector<std::optional<std::string>> raw(10, std::string("c1ccccc1"));
    EXPECT_FALSE(first_failure(raw).has_value());
}

TEST(MurckoFailureScanTest, ReturnsNulloptForAnEmptyBuffer) {
    EXPECT_FALSE(first_failure({}).has_value());
}

TEST(MurckoExtractAllTest, FillsEverySlotAcrossChunks) {
    const std::vector<std::optional<std::string>> raw = extract_all(
        200u, 4u, 64u,
        [](size_t i) -> std::optional<std::string> { return std::to_string(i); });
    ASSERT_EQ(raw.size(), 200u);
    for (size_t i = 0; i < raw.size(); ++i) {
        ASSERT_TRUE(raw[i].has_value()) << "index " << i;
        EXPECT_EQ(*raw[i], std::to_string(i));
    }
}

TEST(MurckoExtractAllTest, PreservesFailureOrderAtEveryThreadCount) {
    // The regression guard for reporting a race winner instead of the first
    // failing index. The two indices must straddle a chunk boundary: at chunk
    // size 64, index 60 is in chunk 0 and index 64 in chunk 1, so different
    // workers can reach them in either time order. A pair inside one chunk --
    // 3 and 7, say -- is always visited in ascending order by a single worker,
    // so it cannot fail against the race this test exists for.
    for (size_t threads : {size_t(1), size_t(8)}) {
        const std::vector<std::optional<std::string>> raw = extract_all(
            200u, threads, 64u, [](size_t i) -> std::optional<std::string> {
                if (i == 64u || i == 60u) {
                    return std::nullopt;
                }
                return std::string("c1ccccc1");
            });
        ASSERT_EQ(raw.size(), 200u) << "threads=" << threads;
        const std::optional<size_t> bad = first_failure(raw);
        ASSERT_TRUE(bad.has_value()) << "threads=" << threads;
        EXPECT_EQ(*bad, 60u) << "threads=" << threads;
    }
}

TEST(MurckoThreadClampTest, ZeroMeansAutoDetectAndSurvivesTheClamp) {
    EXPECT_EQ(effective_thread_count(0u, 0u), 0u);
    EXPECT_EQ(effective_thread_count(0u, 1u), 0u);
    EXPECT_EQ(effective_thread_count(0u, 1000000u), 0u);
}

TEST(MurckoThreadClampTest, ClampsByItemCountAndConcurrencyCeiling) {
    const size_t hw = std::max<size_t>(1u, std::thread::hardware_concurrency());
    const size_t ceiling = 4u * hw;

    EXPECT_EQ(effective_thread_count(1u, 1u), 1u);
    EXPECT_EQ(effective_thread_count(1u, 1000u), 1u);
    // Below both bounds: returned unchanged.
    EXPECT_EQ(effective_thread_count(2u, 1000u), 2u);
    // Above the item count: clamped to the item count.
    EXPECT_EQ(effective_thread_count(1000u, 3u), 3u);
    // Above the concurrency ceiling on a large n. This is the case a clamp by
    // item count alone cannot reach, and the one that closes the ParallelFor
    // reserve-then-fill terminate path.
    EXPECT_EQ(effective_thread_count(1000000u, 1000000u), ceiling);
    // An empty item count must not produce the auto-detect sentinel.
    EXPECT_EQ(effective_thread_count(4u, 0u), 1u);
}

TEST(MurckoMemPoolGateTest, RecognizesOnlyTheThreadSafeModes) {
    // Every single bit OEMemPoolMode defines as of 2026.1.0. The members are
    // namespace-scoped const unsigned int rather than enumerators, so a vendor
    // addition compiles silently and simply goes unlisted -- this table is a
    // record of what was checked, not a compile-time guard.
    EXPECT_TRUE(pool_is_thread_safe(OESystem::OEMemPoolMode::Mutexed));
    EXPECT_TRUE(pool_is_thread_safe(OESystem::OEMemPoolMode::ThreadLocal));
    EXPECT_TRUE(pool_is_thread_safe(OESystem::OEMemPoolMode::Default));

    // SingleThreaded is 0, so it is the absence of a safe bit rather than a
    // flag: asserting it is unsafe is asserting the predicate's default.
    EXPECT_FALSE(pool_is_thread_safe(OESystem::OEMemPoolMode::SingleThreaded));
    // The cache and allocator-strategy bits say nothing about thread safety on
    // their own. Spinlocked is treated as unsafe deliberately: the SDK header
    // describes it as "Spinlock - for Java threading", which names an intended
    // consumer rather than a guarantee for this caller, and the cost of being
    // wrong in this direction is lost parallelism rather than corruption.
    EXPECT_FALSE(pool_is_thread_safe(OESystem::OEMemPoolMode::BoundedCache));
    EXPECT_FALSE(pool_is_thread_safe(OESystem::OEMemPoolMode::UnboundedCache));
    EXPECT_FALSE(pool_is_thread_safe(OESystem::OEMemPoolMode::System));
    EXPECT_FALSE(pool_is_thread_safe(OESystem::OEMemPoolMode::Spinlocked));

    // A safe bit stays safe when a cache bit joins it, which is the only
    // combination Default itself relies on.
    EXPECT_TRUE(pool_is_thread_safe(OESystem::OEMemPoolMode::Mutexed |
                                    OESystem::OEMemPoolMode::BoundedCache));
    EXPECT_FALSE(pool_is_thread_safe(OESystem::OEMemPoolMode::System |
                                     OESystem::OEMemPoolMode::UnboundedCache));
}

TEST(MurckoMemPoolGateTest, ForcesSerialExtractionWhenThePoolIsUnsafe) {
    // The line murcko_scaffolds actually executes, as a pure function. An
    // unsafe mode collapses every thread request to 1, whatever the caller
    // asked for and however many molecules there are. This is the branch whose
    // end-to-end form cannot be written in this binary.
    const unsigned int unsafe = OESystem::OEMemPoolMode::SingleThreaded |
                                OESystem::OEMemPoolMode::UnboundedCache;
    EXPECT_EQ(dispatch_thread_count(unsafe, 0u, 1000u), 1u);
    EXPECT_EQ(dispatch_thread_count(unsafe, 8u, 1000u), 1u);
    EXPECT_EQ(dispatch_thread_count(unsafe, 1000000u, 1000u), 1u);
}

TEST(MurckoMemPoolGateTest, DefersToTheClampWhenThePoolIsSafe) {
    const unsigned int safe = OESystem::OEMemPoolMode::Default;
    // The auto-detect sentinel must survive the gate rather than become 1.
    EXPECT_EQ(dispatch_thread_count(safe, 0u, 1000u), 0u);
    EXPECT_EQ(dispatch_thread_count(safe, 2u, 1000u), 2u);
    EXPECT_EQ(dispatch_thread_count(safe, 1000000u, 1000u),
              effective_thread_count(1000000u, 1000u));
}

TEST(MurckoScaffoldsTest, RefusesAnUnknownScaffoldType) {
    const MolSet mols({"c1ccccc1"});
    OECluster::MurckoOptions options;
    options.scaffold = static_cast<ScaffoldType>(42);
    EXPECT_THROW(OECluster::murcko_scaffolds(mols.Pointers(), options),
                 std::invalid_argument);
}

TEST(MurckoScaffoldsTest, RefusesAnEmptyInput) {
    const std::vector<OEChem::OEMolBase*> mols;
    EXPECT_THROW(OECluster::murcko_scaffolds(mols), OECluster::ComparisonError);
}

TEST(MurckoScaffoldsTest, RefusesANullMolecule) {
    const MolSet mols({"c1ccccc1"});
    std::vector<OEChem::OEMolBase*> with_null = mols.Pointers();
    with_null.push_back(nullptr);
    EXPECT_THROW(OECluster::murcko_scaffolds(with_null),
                 OECluster::ComparisonError);
}

TEST(MurckoScaffoldsTest, JudgesTheOptionsBeforeTheInput) {
    // A call wrong in two ways reports the enumerator -- the cheaper fix, and
    // the one that does not depend on data.
    const std::vector<OEChem::OEMolBase*> mols;
    OECluster::MurckoOptions options;
    options.scaffold = static_cast<ScaffoldType>(42);
    EXPECT_THROW(OECluster::murcko_scaffolds(mols, options),
                 std::invalid_argument);
}

TEST(MurckoScaffoldsTest, MatchesTheKernelForEveryFixture) {
    // The end-to-end form of the one-extraction guarantee: the public entry
    // point returns exactly what the kernel returns, so the extraction table in
    // MurckoExtractionTest covers both.
    const std::vector<std::string> smiles = {
        "c1ccccc1", "Cc1ccccc1", "CC(=O)Oc1ccccc1C(=O)O", "c1ccncc1",
        "c1ccc(cc1)-c1ccccc1", "c1ccc(cc1)Cc1ccccc1", "Cn1cnc2c1c(=O)n(C)c(=O)n2C",
        "O=C1CCCCN1", "C1CCCCC1", "FC1=CCCCC1", "FC1CCCCC1",
        "c1ccccc1C(=O)O.c1ccncc1", "CCCCCC", "C", "",
        // The stereo fixtures travel the driver too: the kernel strips
        // stereochemistry, and nothing on the threaded path may reintroduce it.
        "C1CC[C@H]2CCCC[C@@H]2C1", "C1CC[C@H]2CCCC[C@H]2C1",
        "C[C@H](N)c1ccccc1", "C[C@@H](N)c1ccccc1",
    };
    const MolSet mols(smiles);
    for (ScaffoldType type : {ScaffoldType::Framework, ScaffoldType::Generic}) {
        OECluster::MurckoOptions options;
        options.scaffold = type;
        const std::vector<std::string> actual =
            OECluster::murcko_scaffolds(mols.Pointers(), options);
        ASSERT_EQ(actual.size(), smiles.size());
        for (size_t i = 0; i < smiles.size(); ++i) {
            const std::optional<std::string> expected =
                scaffold_of(*mols.Pointers()[i], type);
            ASSERT_TRUE(expected.has_value()) << "index " << i;
            EXPECT_EQ(actual[i], *expected) << "index " << i;
        }
    }
}

TEST(MurckoScaffoldsTest, RaisesNothingForOrdinaryMolecules) {
    const MolSet mols(threading_fixture_smiles());
    EXPECT_NO_THROW(OECluster::murcko_scaffolds(mols.Pointers()));
}

TEST(MurckoScaffoldsTest, IsDeterministicAcrossThreadCounts) {
    // A regression smoke test, not a proof of thread safety: the safety
    // argument is the checked mem-pool precondition.
    const MolSet mols(threading_fixture_smiles());
    OECluster::MurckoOptions serial;
    serial.num_threads = 1;
    const std::vector<std::string> expected =
        OECluster::murcko_scaffolds(mols.Pointers(), serial);
    ASSERT_EQ(expected.size(), 200u);

    for (size_t threads : {size_t(2), size_t(8)}) {
        OECluster::MurckoOptions options;
        options.num_threads = threads;
        EXPECT_EQ(OECluster::murcko_scaffolds(mols.Pointers(), options), expected)
            << "threads=" << threads;
    }
}

TEST(MurckoScaffoldsTest, ClampsAnAbsurdThreadRequest) {
    const MolSet mols({"c1ccccc1", "Cc1ccccc1", "CCCCCC"});
    OECluster::MurckoOptions serial;
    serial.num_threads = 1;
    OECluster::MurckoOptions absurd;
    absurd.num_threads = 1000000;
    EXPECT_EQ(OECluster::murcko_scaffolds(mols.Pointers(), absurd),
              OECluster::murcko_scaffolds(mols.Pointers(), serial));
}

// The end-to-end serial-fallback case cannot be written in this binary:
// OESetMemPoolMode is fatal when called a second time in a process ("Fatal:
// OESetMemPoolMode called twice!"), so a test cannot restore the mode it
// changed, and leaving it changed would silently reconfigure every later test.
// MurckoMemPoolGateTest covers the decision instead, including the wiring line
// itself via dispatch_thread_count. What remains untested is only the call to
// OEGetMemPoolMode that feeds it.

}  // namespace
