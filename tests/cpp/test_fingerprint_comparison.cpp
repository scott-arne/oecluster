#include <gtest/gtest.h>
#include "oecluster/oecluster.h"
#include "oecluster/comparisons/FingerprintComparison.h"
#include "oecluster/GateFacts.h"
#include "oecluster/StorageBackend.h"
#include "oecluster/PDist.h"
#include "oecluster/Error.h"
#include <oechem.h>
#include <cmath>

using namespace OECluster;

class FingerprintComparisonTest : public ::testing::Test {
protected:
    void SetUp() override {
        auto add_mol = [this](const char* smi) {
            graph_mols_.emplace_back();
            OEChem::OESmilesToMol(graph_mols_.back(), smi);
        };
        add_mol("c1ccccc1");      // benzene
        add_mol("c1ccc(O)cc1");   // phenol
        add_mol("CCCCCCCC");      // octane

        for (auto& gm : graph_mols_) {
            mols_.push_back(&static_cast<OEChem::OEMolBase&>(gm));
        }
    }

    std::vector<OEChem::OEGraphMol> graph_mols_;
    std::vector<OEChem::OEMolBase*> mols_;
};

TEST_F(FingerprintComparisonTest, ConstructAndSize) {
    FingerprintOptions opts;
    EXPECT_EQ(opts.fp_type, "morgan");

    FingerprintComparison comparison(mols_);
    EXPECT_EQ(comparison.Size(), 3);
    EXPECT_EQ(comparison.ComparisonName(), "fingerprint");
}

TEST_F(FingerprintComparisonTest, SimilarMoleculesCloser) {
    FingerprintComparison comparison(mols_);
    double d_benzene_phenol = comparison.Compare(0, 1);
    double d_benzene_octane = comparison.Compare(0, 2);
    EXPECT_LT(d_benzene_phenol, d_benzene_octane);
}

TEST_F(FingerprintComparisonTest, DistanceRange) {
    FingerprintComparison comparison(mols_);
    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = i + 1; j < 3; ++j) {
            double d = comparison.Compare(i, j);
            EXPECT_GE(d, 0.0);
            EXPECT_LE(d, 1.0);
        }
    }
}

TEST_F(FingerprintComparisonTest, CloneSharesData) {
    FingerprintComparison comparison(mols_);
    auto clone = comparison.Clone();
    EXPECT_EQ(clone->Size(), 3);
    EXPECT_DOUBLE_EQ(clone->Compare(0, 1), comparison.Compare(0, 1));
}

TEST_F(FingerprintComparisonTest, IntegrationWithPDist) {
    FingerprintComparison comparison(mols_);
    DenseStorage storage(3);
    pdist(comparison, storage);
    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = i + 1; j < 3; ++j) {
            EXPECT_GT(storage.Get(i, j), 0.0);
        }
    }
}

// SparseStorage is the only backend whose Data() is null, so it is the only one
// that reaches TryPDist's fallback branch through the shared condensed_to_pair.
TEST_F(FingerprintComparisonTest, IntegrationWithPDistSparseStorage) {
    FingerprintComparison comparison(mols_);

    DenseStorage dense(3);
    pdist(comparison, dense);

    SparseStorage sparse(3, 1.0);
    pdist(comparison, sparse);

    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = i + 1; j < 3; ++j) {
            EXPECT_GT(sparse.Get(i, j), 0.0);
            EXPECT_DOUBLE_EQ(sparse.Get(i, j), dense.Get(i, j));
        }
    }
}

TEST_F(FingerprintComparisonTest, MorganFingerprintType) {
    FingerprintOptions opts;
    opts.fp_type = "morgan";
    opts.numbits = 2048;
    opts.radius = 2;
    FingerprintComparison comparison(mols_, opts);
    EXPECT_EQ(comparison.Size(), 3);
    double d = comparison.Compare(0, 1);
    EXPECT_GE(d, 0.0);
    EXPECT_LE(d, 1.0);
}

TEST_F(FingerprintComparisonTest, AtomPairFingerprintType) {
    FingerprintOptions opts;
    opts.fp_type = "atom_pair";
    FingerprintComparison comparison(mols_, opts);
    double d = comparison.Compare(0, 1);
    EXPECT_GE(d, 0.0);
    EXPECT_LE(d, 1.0);
}

TEST_F(FingerprintComparisonTest, TanimotoSimilarityComplementsJaccardDistance) {
    FingerprintOptions opts;
    opts.similarity = false;
    FingerprintComparison distance_comparison(mols_, opts);

    opts.similarity = true;
    FingerprintComparison similarity_comparison(mols_, opts);

    EXPECT_NEAR(
        distance_comparison.Compare(0, 1) + similarity_comparison.Compare(0, 1),
        1.0,
        1.0e-12);
}

TEST_F(FingerprintComparisonTest, DiceMetric) {
    FingerprintOptions opts;
    opts.metric = "dice";
    FingerprintComparison comparison(mols_, opts);
    double d = comparison.Compare(0, 1);
    EXPECT_GE(d, 0.0);
    EXPECT_LE(d, 1.0);
}

TEST_F(FingerprintComparisonTest, RemovedOpenEyeFingerprintTypesThrow) {
    FingerprintOptions opts;

    for (const auto* fp_type : {"circular", "tree", "path", "maccs", "lingo"}) {
        opts.fp_type = fp_type;
        EXPECT_THROW(FingerprintComparison(mols_, opts), ComparisonError)
            << "fp_type=" << fp_type;
    }
}

TEST_F(FingerprintComparisonTest, EuclideanMetric) {
    FingerprintOptions opts;
    opts.metric = "euclidean";
    FingerprintComparison comparison(mols_, opts);
    // Euclidean over a bit-set is sqrt(popcount of the symmetric difference),
    // which is above 1.0 here — no bit-set overlap coefficient can return this,
    // so the assertion fails if the table row stops resolving to Euclidean.
    EXPECT_DOUBLE_EQ(comparison.Compare(0, 1), std::sqrt(8.0));
}

TEST_F(FingerprintComparisonTest, InvalidFpTypeThrows) {
    FingerprintOptions opts;
    opts.fp_type = "invalid";
    EXPECT_THROW(FingerprintComparison(mols_, opts), ComparisonError);
}

TEST_F(FingerprintComparisonTest, DefaultOptionsMatchTheSpecifiedDefaults) {
    FingerprintOptions opts;
    EXPECT_EQ(opts.fp_type, "morgan");
    EXPECT_EQ(opts.storage, "binary");
    EXPECT_EQ(opts.numbits, 2048u);
    EXPECT_EQ(opts.metric, "tanimoto");
    EXPECT_FALSE(opts.similarity);
    EXPECT_EQ(opts.radius, 2u);
    EXPECT_EQ(opts.min_distance, 1u);
    EXPECT_EQ(opts.max_distance, 30u);
    EXPECT_EQ(opts.torsion_atom_count, 4u);
    EXPECT_FALSE(opts.use_chirality);
    EXPECT_DOUBLE_EQ(opts.p, 2.0);
    EXPECT_DOUBLE_EQ(opts.tversky_alpha, 0.5);
    EXPECT_DOUBLE_EQ(opts.tversky_beta, 0.5);
}

TEST_F(FingerprintComparisonTest, DefaultFactsPassBothGateTiers) {
    FingerprintComparison comparison(mols_);
    const GateFacts facts = comparison.Facts();
    EXPECT_EQ(facts.is_distance, Capability::Yes);
    EXPECT_EQ(facts.zero_self, Capability::Yes);
    EXPECT_EQ(facts.triangle, Capability::Yes);
    EXPECT_EQ(facts.data_integrity, DataIntegrity::Complete);
}

TEST_F(FingerprintComparisonTest, SimilarityFactsFailTierOne) {
    FingerprintOptions opts;
    opts.similarity = true;
    FingerprintComparison comparison(mols_, opts);
    const GateFacts facts = comparison.Facts();
    EXPECT_EQ(facts.is_distance, Capability::No);
    EXPECT_EQ(facts.zero_self, Capability::No);
}

TEST_F(FingerprintComparisonTest, DiceFactsFailTierTwoOnly) {
    FingerprintOptions opts;
    opts.metric = "dice";
    FingerprintComparison comparison(mols_, opts);
    const GateFacts facts = comparison.Facts();
    EXPECT_EQ(facts.is_distance, Capability::Yes);
    EXPECT_EQ(facts.zero_self, Capability::Yes);
    EXPECT_EQ(facts.triangle, Capability::No);
}

TEST_F(FingerprintComparisonTest, TverskyIsStampedASimilarityInBothDirections) {
    // The metric table builds OEFP's Tversky for tversky whatever the
    // similarity flag says (src/comparisons/MetricTable.cpp), because Tversky
    // has no distance form. The flag is therefore ignored here, and both
    // directions must report the same orientation rather than one of them
    // quietly claiming to be a distance.
    for (bool similarity : {false, true}) {
        FingerprintOptions opts;
        opts.metric = "tversky";
        opts.similarity = similarity;
        FingerprintComparison comparison(mols_, opts);
        const GateFacts facts = comparison.Facts();
        EXPECT_EQ(facts.is_distance, Capability::No) << similarity;
        EXPECT_EQ(facts.zero_self, Capability::No) << similarity;
    }
}

TEST_F(FingerprintComparisonTest, MorganRadiusIsSeparateFromAtomPairWindow) {
    // Exact values are pinned because an off-by-one in the radius forwarding
    // survives any inequality assertion.
    FingerprintOptions opts;
    opts.radius = 1;
    FingerprintComparison r1(mols_, opts);
    EXPECT_DOUBLE_EQ(r1.Compare(0, 1), 5.0 / 7.0);

    opts.radius = 2;
    FingerprintComparison r2(mols_, opts);
    EXPECT_DOUBLE_EQ(r2.Compare(0, 1), 8.0 / 11.0);

    opts.radius = 3;
    FingerprintComparison r3(mols_, opts);
    EXPECT_DOUBLE_EQ(r3.Compare(0, 1), 11.0 / 14.0);
}

TEST_F(FingerprintComparisonTest, MorganUseChiralityDistinguishesEnantiomers) {
    // The fixture's three molecules are achiral, so this test builds its own
    // pair: L- and D-alanine, identical except at the stereocenter.
    std::vector<OEChem::OEGraphMol> chiral_mols(2);
    OEChem::OESmilesToMol(chiral_mols[0], "C[C@H](N)C(=O)O");
    OEChem::OESmilesToMol(chiral_mols[1], "C[C@@H](N)C(=O)O");
    std::vector<OEChem::OEMolBase*> enantiomers;
    for (auto& gm : chiral_mols) {
        enantiomers.push_back(&static_cast<OEChem::OEMolBase&>(gm));
    }

    FingerprintOptions achiral;
    FingerprintComparison blind(enantiomers, achiral);
    EXPECT_DOUBLE_EQ(blind.Compare(0, 1), 0.0);

    FingerprintOptions chiral = achiral;
    chiral.use_chirality = true;
    FingerprintComparison aware(enantiomers, chiral);
    EXPECT_GT(aware.Compare(0, 1), 0.0);
}

TEST_F(FingerprintComparisonTest, AtomPairUseChiralityDistinguishesSpecifiedStereo) {
    // Atom-pair enantiomers do NOT distinguish with use_chirality at OEFP 0.3.0;
    // the flag widens the atom code but does not flip the bit pattern for
    // opposite CIP labels. What does work is specified versus unspecified: a
    // resolved stereocenter stops matching the same center with no label. The
    // stereocenter must carry four heavy substituents; a center with an implicit
    // hydrogen does not resolve through OEFP molecule preparation.
    std::vector<OEChem::OEGraphMol> stereo_mols(2);
    OEChem::OESmilesToMol(stereo_mols[0], "F[C@](Cl)(Br)I");
    OEChem::OESmilesToMol(stereo_mols[1], "FC(Cl)(Br)I");
    std::vector<OEChem::OEMolBase*> mols;
    for (auto& gm : stereo_mols) {
        mols.push_back(&static_cast<OEChem::OEMolBase&>(gm));
    }

    FingerprintOptions achiral;
    achiral.fp_type = "atom_pair";
    FingerprintComparison blind(mols, achiral);
    EXPECT_DOUBLE_EQ(blind.Compare(0, 1), 0.0);

    FingerprintOptions chiral = achiral;
    chiral.use_chirality = true;
    FingerprintComparison aware(mols, chiral);
    EXPECT_DOUBLE_EQ(aware.Compare(0, 1), 0.5714285714285714);
}

TEST_F(FingerprintComparisonTest, AtomPairDefaultsUseTheFullOEFPWindow) {
    // The fixture's benzene and phenol are too short; their maximum topological
    // distance saturates at 4, so any window from [1, 4] upward gives identical
    // distances and an off-by-one below 30 would survive. (Benzene vs octane was
    // abandoned for saturating at 1.0.) Instead, use a 31-carbon chain and an
    // oxygen + 30-carbon chain so the top of the window is reachable.
    std::vector<OEChem::OEGraphMol> long_mols(2);
    OEChem::OESmilesToMol(long_mols[0], std::string(31, 'C'));
    OEChem::OESmilesToMol(long_mols[1], "O" + std::string(30, 'C'));
    std::vector<OEChem::OEMolBase*> chains;
    for (auto& gm : long_mols) {
        chains.push_back(&static_cast<OEChem::OEMolBase&>(gm));
    }

    FingerprintOptions default_opts;
    default_opts.fp_type = "atom_pair";
    FingerprintComparison defaults(chains, default_opts);

    FingerprintOptions explicit_30 = default_opts;
    explicit_30.min_distance = 1;
    explicit_30.max_distance = 30;
    FingerprintComparison full_window(chains, explicit_30);

    // The default equals the explicit [1, 30] window, which is what pins the
    // default max_distance at exactly 30.
    EXPECT_DOUBLE_EQ(defaults.Compare(0, 1), full_window.Compare(0, 1));

    FingerprintOptions truncated = default_opts;
    truncated.min_distance = 1;
    truncated.max_distance = 29;
    FingerprintComparison narrow(chains, truncated);

    // And [1, 29] differs, so any off-by-one or clamp below 30 fails this.
    EXPECT_NE(defaults.Compare(0, 1), narrow.Compare(0, 1));
}

TEST_F(FingerprintComparisonTest, MinkowskiExponentIsHonored) {
    FingerprintOptions opts;
    opts.metric = "minkowski";
    opts.p = 0.5;
    FingerprintComparison comparison(mols_, opts);
    EXPECT_EQ(comparison.Facts().triangle, Capability::No);
}

TEST_F(FingerprintComparisonTest, AsymmetricTverskyIsRefusedByPDist) {
    FingerprintOptions opts;
    opts.metric = "tversky";
    opts.tversky_alpha = 0.2;
    opts.tversky_beta = 0.8;
    FingerprintComparison comparison(mols_, opts);

    DenseStorage storage(comparison.Size());
    PDistOptions pdist_opts;
    try {
        comparison.TryPDist(storage, pdist_opts);
        FAIL() << "expected ComparisonError";
    } catch (const ComparisonError& error) {
        // The OEFP kernel refuses this too, but with its own wording. Asserting
        // the OECluster phrasing is what distinguishes our guard from the
        // fallback conversion in the surrounding catch block.
        const std::string message(error.what());
        EXPECT_NE(message.find("cannot be used with pdist"), std::string::npos)
            << message;
        EXPECT_NE(message.find("tversky"), std::string::npos) << message;
    }
}

TEST_F(FingerprintComparisonTest, AsymmetricTverskyGuardRunsBeforeSizeCheck) {
    // The metric guard runs first, so asymmetric Tversky must fail with the
    // metric error even when storage size is wrong.
    FingerprintOptions opts;
    opts.metric = "tversky";
    opts.tversky_alpha = 0.2;
    opts.tversky_beta = 0.8;
    FingerprintComparison comparison(mols_, opts);

    DenseStorage storage(comparison.Size() + 1);
    PDistOptions pdist_opts;
    try {
        comparison.TryPDist(storage, pdist_opts);
        FAIL() << "expected ComparisonError";
    } catch (const ComparisonError& error) {
        const std::string message(error.what());
        EXPECT_NE(message.find("cannot be used with pdist"), std::string::npos)
            << message;
        EXPECT_EQ(message.find("storage size mismatch"), std::string::npos)
            << message;
    }
}

TEST_F(FingerprintComparisonTest, MahalanobisIsRefusedOnTheFingerprintSurface) {
    FingerprintOptions opts;
    opts.metric = "mahalanobis";
    EXPECT_THROW(FingerprintComparison(mols_, opts), ComparisonError);
}

TEST_F(FingerprintComparisonTest, AllFourFamiliesConstructAndScore) {
    // Exact benzene/phenol Jaccard distances at each family's defaults,
    // measured against OEFP 0.3.0. A range assertion here would pass even if
    // every family silently fell back to the same generator.
    const std::vector<std::pair<std::string, double>> expected{
        {"morgan", 0.72727272727272729},
        {"atom_pair", 0.57894736842105265},
        {"topological_atom_pair", 0.57894736842105265},
        {"topological_torsions", 0.77777777777777779},
    };

    for (const auto& [family, distance] : expected) {
        FingerprintOptions opts;
        opts.fp_type = family;
        FingerprintComparison comparison(mols_, opts);
        EXPECT_EQ(comparison.Size(), 3u) << family;
        EXPECT_DOUBLE_EQ(comparison.Compare(0, 1), distance) << family;
    }
}

TEST_F(FingerprintComparisonTest, TopologicalAtomPairIsAnAliasOfAtomPair) {
    FingerprintOptions plain;
    plain.fp_type = "atom_pair";
    FingerprintOptions aliased;
    aliased.fp_type = "topological_atom_pair";

    FingerprintComparison a(mols_, plain);
    FingerprintComparison b(mols_, aliased);
    EXPECT_DOUBLE_EQ(a.Compare(0, 1), b.Compare(0, 1));
    EXPECT_DOUBLE_EQ(a.Compare(0, 2), b.Compare(0, 2));
}

TEST_F(FingerprintComparisonTest, TorsionAtomCountChangesTheFingerprint) {
    // Benzene/phenol, not benzene/octane: octane shares no torsion with either
    // aromatic at any atom count, so that pair sits at distance 1.0 for every
    // setting and cannot observe the option at all.
    const std::vector<std::pair<unsigned int, double>> expected{
        {3u, 0.75},
        {4u, 0.77777777777777779},
        {5u, 0.90000000000000002},
    };

    for (const auto& [count, distance] : expected) {
        FingerprintOptions opts;
        opts.fp_type = "topological_torsions";
        opts.torsion_atom_count = count;
        FingerprintComparison comparison(mols_, opts);
        EXPECT_DOUBLE_EQ(comparison.Compare(0, 1), distance) << count;
    }

    // The default must be 4, not merely "one of the three".
    FingerprintOptions defaulted;
    defaulted.fp_type = "topological_torsions";
    FingerprintComparison implicit(mols_, defaulted);
    EXPECT_DOUBLE_EQ(implicit.Compare(0, 1), 0.77777777777777779);
}

TEST_F(FingerprintComparisonTest, DistanceAtomPairIsRejectedAsUnimplemented) {
    FingerprintOptions opts;
    opts.fp_type = "distance_atom_pair";
    try {
        FingerprintComparison comparison(mols_, opts);
        FAIL() << "expected ComparisonError";
    } catch (const ComparisonError& error) {
        // The dedicated branch exists only for this message. Asserting the
        // exception type alone cannot tell it from the generic unknown-type
        // error the name would otherwise fall through to.
        const std::string message(error.what());
        EXPECT_NE(message.find("not implemented by OEFP"), std::string::npos)
            << message;
        EXPECT_NE(message.find("'atom_pair'"), std::string::npos) << message;
        EXPECT_EQ(message.find("Unknown OEFP fingerprint type"), std::string::npos)
            << message;
    }
}

TEST_F(FingerprintComparisonTest, LegacyOpenEyeTypesStillNameTheReplacements) {
    for (const char* legacy : {"circular", "tree", "path", "maccs", "lingo"}) {
        FingerprintOptions opts;
        opts.fp_type = legacy;
        try {
            FingerprintComparison comparison(mols_, opts);
            FAIL() << "expected ComparisonError for " << legacy;
        } catch (const ComparisonError& error) {
            const std::string message(error.what());
            EXPECT_NE(message.find(legacy), std::string::npos) << message;
            EXPECT_NE(message.find("no longer supported"), std::string::npos)
                << message;
            // The point of the branch: every current family is offered.
            for (const char* replacement : {"'morgan'", "'atom_pair'",
                                            "'topological_atom_pair'",
                                            "'topological_torsions'"}) {
                EXPECT_NE(message.find(replacement), std::string::npos)
                    << legacy << ": " << message;
            }
        }
    }
}

TEST_F(FingerprintComparisonTest, AlternateFamilySpellingsResolveToTheSameFamily) {
    // normalize_family accepts more spellings than the four canonical names.
    // Each alternate must land on the same generator as its canonical form,
    // not merely construct without throwing.
    const std::vector<std::pair<std::string, std::string>> aliases{
        {"atompair", "atom_pair"},
        {"ATOM_PAIR", "atom_pair"},
        {"Morgan", "morgan"},
        {"topological_torsion", "topological_torsions"},
        {"TOPOLOGICAL_ATOM_PAIR", "atom_pair"},
    };

    for (const auto& [alternate, canonical] : aliases) {
        FingerprintOptions alternate_opts;
        alternate_opts.fp_type = alternate;
        FingerprintOptions canonical_opts;
        canonical_opts.fp_type = canonical;

        FingerprintComparison a(mols_, alternate_opts);
        FingerprintComparison b(mols_, canonical_opts);
        EXPECT_DOUBLE_EQ(a.Compare(0, 1), b.Compare(0, 1)) << alternate;
        EXPECT_DOUBLE_EQ(a.Compare(1, 2), b.Compare(1, 2)) << alternate;
    }
}

TEST_F(FingerprintComparisonTest, UnknownFamilyErrorNamesTheSupportedTypes) {
    FingerprintOptions opts;
    opts.fp_type = "ecfp4";
    try {
        FingerprintComparison comparison(mols_, opts);
        FAIL() << "expected ComparisonError";
    } catch (const ComparisonError& error) {
        const std::string message(error.what());
        // The unrecognized name must be echoed back, or the user cannot tell
        // which of several options was rejected.
        EXPECT_NE(message.find("ecfp4"), std::string::npos) << message;
        for (const char* supported : {"'morgan'", "'atom_pair'",
                                      "'topological_atom_pair'",
                                      "'topological_torsions'"}) {
            EXPECT_NE(message.find(supported), std::string::npos) << message;
        }
    }
}

TEST_F(FingerprintComparisonTest, NumbitsReachesEveryFamilyGenerator) {
    // The three-molecule fixture is too small to collide at any bit width, so
    // numbits is unobservable there. These two tripeptides differ enough to
    // fold differently at 64 bits and identically at neither.
    std::vector<OEChem::OEGraphMol> peptides(2);
    OEChem::OESmilesToMol(peptides[0], "CC(C)CC(N)C(=O)NC(Cc1ccccc1)C(=O)NC(CO)C(=O)O");
    OEChem::OESmilesToMol(peptides[1], "CC(C)CC(N)C(=O)NC(Cc1ccc(O)cc1)C(=O)NC(CS)C(=O)O");
    std::vector<OEChem::OEMolBase*> pair{&static_cast<OEChem::OEMolBase&>(peptides[0]),
                                         &static_cast<OEChem::OEMolBase&>(peptides[1])};

    // {family, distance at 64 bits, distance at the 2048-bit default}
    const std::vector<std::tuple<std::string, double, double>> expected{
        {"morgan", 0.30769230769230771, 0.34000000000000002},
        {"atom_pair", 0.0, 0.28614457831325302},
        {"topological_torsions", 0.096774193548387094, 0.29545454545454547},
    };

    for (const auto& [family, folded, unfolded] : expected) {
        FingerprintOptions narrow;
        narrow.fp_type = family;
        narrow.numbits = 64;
        FingerprintComparison narrow_comparison(pair, narrow);
        EXPECT_DOUBLE_EQ(narrow_comparison.Compare(0, 1), folded) << family << " at 64 bits";

        FingerprintOptions wide;
        wide.fp_type = family;
        FingerprintComparison wide_comparison(pair, wide);
        EXPECT_DOUBLE_EQ(wide_comparison.Compare(0, 1), unfolded) << family << " at 2048 bits";
    }
}

TEST_F(FingerprintComparisonTest, FifteenSupportedFamilyStorageCombinations) {
    const std::vector<std::string> families{"morgan", "atom_pair", "topological_atom_pair",
                                            "topological_torsions"};
    const std::vector<std::string> storages{"binary", "count", "sparse", "sparse_count"};

    size_t supported = 0;
    for (const std::string& family : families) {
        for (const std::string& storage : storages) {
            FingerprintOptions opts;
            opts.fp_type = family;
            opts.storage = storage;
            // Counts are real-valued, so the default boolean metric does not apply.
            opts.metric = (storage == "count" || storage == "sparse_count") ? "manhattan"
                                                                            : "tanimoto";
            const bool unsupported =
                (family == "topological_torsions" && storage == "sparse_count");
            if (unsupported) {
                EXPECT_THROW(FingerprintComparison(mols_, opts), ComparisonError)
                    << family << "/" << storage;
                continue;
            }
            FingerprintComparison comparison(mols_, opts);
            EXPECT_EQ(comparison.Size(), 3u) << family << "/" << storage;
            EXPECT_GE(comparison.Compare(0, 1), 0.0) << family << "/" << storage;
            ++supported;
        }
    }
    EXPECT_EQ(supported, 15u);
}

TEST_F(FingerprintComparisonTest, CountStorageDistinguishesRepeatedSubstructures) {
    std::vector<OEChem::OEGraphMol> alkanes(2);
    OEChem::OESmilesToMol(alkanes[0], "CCCCCCCCCCCCCCCC");
    OEChem::OESmilesToMol(alkanes[1], "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC");
    std::vector<OEChem::OEMolBase*> pair{&static_cast<OEChem::OEMolBase&>(alkanes[0]),
                                         &static_cast<OEChem::OEMolBase&>(alkanes[1])};

    FingerprintOptions binary;
    FingerprintComparison binary_comparison(pair, binary);
    EXPECT_DOUBLE_EQ(binary_comparison.Compare(0, 1), 0.0);

    // The metric here must read magnitudes, not just which bits are on. Both
    // alkanes occupy the same eight Morgan indices and differ only in their
    // counts (46 vs 88 total), so every set-based metric collapses: tanimoto
    // returns 1.0 and dice 0.0 on the counted fingerprints exactly as they do
    // on the binary ones. bray_curtis returns 0.3134. Do not substitute a
    // boolean metric here — the test would then assert the opposite of its name.
    FingerprintOptions counted;
    counted.storage = "count";
    counted.metric = "bray_curtis";
    FingerprintComparison count_comparison(pair, counted);
    // 42/134: sum|a-b| over sum(a+b) across the eight shared Morgan indices.
    EXPECT_NEAR(count_comparison.Compare(0, 1), 42.0 / 134.0, 1e-9);
}

TEST_F(FingerprintComparisonTest, BooleanMetricOnCountedStorageIsRejected) {
    for (const char* storage : {"count", "sparse_count"}) {
        FingerprintOptions explicit_metric;
        explicit_metric.storage = storage;
        explicit_metric.metric = "tanimoto";
        EXPECT_THROW(FingerprintComparison(mols_, explicit_metric), ComparisonError) << storage;

        // tanimoto is also the default, so a user who changes only the storage
        // meets the same rejection without ever naming a metric.
        FingerprintOptions default_metric;
        default_metric.storage = storage;
        EXPECT_THROW(FingerprintComparison(mols_, default_metric), ComparisonError) << storage;
    }
}

TEST_F(FingerprintComparisonTest, UnknownStorageIsRejected) {
    FingerprintOptions opts;
    opts.storage = "dense";
    EXPECT_THROW(FingerprintComparison(mols_, opts), ComparisonError);
}

TEST_F(FingerprintComparisonTest, PDistAgreesWithComparePairForEveryStorage) {
    for (const char* storage : {"binary", "count", "sparse", "sparse_count"}) {
        FingerprintOptions opts;
        opts.storage = storage;
        opts.metric = (std::string(storage).find("count") != std::string::npos) ? "manhattan"
                                                                                : "tanimoto";
        FingerprintComparison comparison(mols_, opts);

        DenseStorage storage_backend(comparison.Size());
        PDistOptions pdist_opts;
        pdist_opts.num_threads = 1;
        ASSERT_TRUE(comparison.TryPDist(storage_backend, pdist_opts)) << storage;
        EXPECT_NEAR(storage_backend.Get(0, 1), comparison.Compare(0, 1), 1e-12) << storage;
        EXPECT_NEAR(storage_backend.Get(1, 2), comparison.Compare(1, 2), 1e-12) << storage;
    }
}

TEST_F(FingerprintComparisonTest, CDistAgreesWithComparePairForEveryStorage) {
    for (const char* storage : {"binary", "count", "sparse", "sparse_count"}) {
        FingerprintOptions opts;
        opts.storage = storage;
        // manhattan is defined on all four storages and is not a bit-set metric,
        // so one metric covers the whole loop without tripping the counted-storage
        // rule.
        opts.metric = "manhattan";
        FingerprintComparison comparison(mols_, opts);

        std::vector<double> output(1 * 2, -1.0);
        CDistOptions cdist_opts;
        cdist_opts.num_threads = 1;
        ASSERT_TRUE(comparison.TryCDist(1, output.data(), cdist_opts)) << storage;
        EXPECT_NEAR(output[0], comparison.Compare(0, 1), 1e-12) << storage;
        EXPECT_NEAR(output[1], comparison.Compare(0, 2), 1e-12) << storage;
    }
}

TEST_F(FingerprintComparisonTest, TorsionAtomCountReachesCountStorage) {
    FingerprintOptions opts_4;
    opts_4.fp_type = "topological_torsions";
    opts_4.storage = "count";
    opts_4.metric = "manhattan";
    opts_4.torsion_atom_count = 4;
    FingerprintComparison comparison_4(mols_, opts_4);

    FingerprintOptions opts_5 = opts_4;
    opts_5.torsion_atom_count = 5;
    FingerprintComparison comparison_5(mols_, opts_5);

    EXPECT_NE(comparison_4.Compare(0, 1), comparison_5.Compare(0, 1));
}

TEST_F(FingerprintComparisonTest, UnsupportedCellIsNamedBeforeTheMetricRule) {
    // The default metric is tanimoto, so this request violates the counted-storage
    // metric rule as well. It must be told the thing it can act on: no choice of
    // metric makes this cell work.
    FingerprintOptions opts;
    opts.fp_type = "topological_torsions";
    opts.storage = "sparse_count";
    try {
        FingerprintComparison comparison(mols_, opts);
        FAIL() << "expected ComparisonError";
    } catch (const ComparisonError& error) {
        const std::string message(error.what());
        EXPECT_NE(message.find("does not support"), std::string::npos) << message;
        EXPECT_NE(message.find("sparse_count"), std::string::npos) << message;
        EXPECT_EQ(message.find("bit-set metric"), std::string::npos) << message;
    }
}
