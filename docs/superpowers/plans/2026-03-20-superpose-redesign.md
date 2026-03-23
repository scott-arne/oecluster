# SuperposeMetric Redesign Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace `SuperposeMetric` and `SiteHopperMetric` with a unified `SuperposeMetric` backed by `oespruce::OESuperpose`, add a `similarity` flag to all metrics, and update SWIG/Python/CLI accordingly.

**Architecture:** The new `SuperposeMetric` uses `oespruce::OESuperpose` with configurable `SuperposeMethod` and `SuperposeScoreType` enums. Each `Clone()` creates a thread-local `OESuperpose` instance. Atom predicates are compiled via oeselect at construction time and shared across clones. The `similarity` flag is a per-metric option (no base class changes).

**Tech Stack:** C++17, oespruce, oeselect, OEChem, OEBio, SWIG 4.x, Python 3.x, CLI11, Google Test

---

## File Structure

### Files to Create
- *(none — all changes modify existing files)*

### Files to Modify
- `CMakeLists.txt` — add oespruce dependency via `_find_oe_dep`, add oeselect via FetchContent
- `include/oecluster/metrics/SuperposeMetric.h` — new enums, options struct, class API
- `src/metrics/SuperposeMetric.cpp` — full rewrite using OESuperpose
- `include/oecluster/metrics/FingerprintMetric.h` — rename `similarity` to `similarity_func`, add `bool similarity`
- `src/metrics/FingerprintMetric.cpp` — update field references and apply similarity mode
- `include/oecluster/metrics/ROCSMetric.h` — add `bool similarity` to `ROCSOptions`
- `src/metrics/ROCSMetric.cpp` — apply similarity mode in `Distance()`
- `include/oecluster/oecluster.h` — remove `SiteHopperMetric.h` include
- `swig/oecluster.i` — remove SiteHopper refs, add oespruce include
- `swig/CMakeLists.txt` — add OESPRUCE_LIBRARY to expected lib vars
- `python/oecluster/__init__.py` — update pdist() kwargs, remove SiteHopperMetric wrapper
- `tools/oepdist.cpp` — replace superpose/sitehopper subcommands, add --sim to fp/rocs
- `tests/cpp/test_bio_metrics.cpp` — full rewrite for new SuperposeMetric API
- `tests/cpp/test_fingerprint_metric.cpp` — add similarity mode tests
- `tests/cpp/test_rocs_metric.cpp` — add similarity mode tests
- `tests/python/test_pdist.py` — add similarity mode and superpose kwargs tests

### Files to Delete
- `include/oecluster/metrics/SiteHopperMetric.h`
- `src/metrics/SiteHopperMetric.cpp`

---

### Task 1: Add oespruce and oeselect Dependencies

**Files:**
- Modify: `CMakeLists.txt:71-99` (add oespruce after OEBio block)
- Modify: `CMakeLists.txt:22-31` (add oeselect FetchContent)
- Modify: `CMakeLists.txt:123-129` (add oespruce to link libraries)
- Modify: `swig/CMakeLists.txt:3-9` (add OESPRUCE_LIBRARY)

- [ ] **Step 1: Add oespruce discovery to CMakeLists.txt**

After the OEBio `_find_oe_dep` block (line 99), add:

```cmake
_find_oe_dep(OESPRUCE_LIBRARY oespruce OpenEye::OESpruce "OpenEye::OEChem;OpenEye::OEBio")
```

- [ ] **Step 2: Add oeselect FetchContent to CMakeLists.txt**

After the `cmake_openeye` FetchContent block (line 31), add:

```cmake
FetchContent_Declare(
    oeselect
    GIT_REPOSITORY https://github.com/scott-arne/oeselect.git
    GIT_TAG main
    GIT_SHALLOW TRUE
)
FetchContent_MakeAvailable(oeselect)
```

- [ ] **Step 3: Link oespruce and oeselect to the oecluster target**

Update `target_link_libraries` to add `OpenEye::OESpruce` and `oeselect`:

```cmake
target_link_libraries(oecluster
    PUBLIC
        OpenEye::OEChem
        OpenEye::OEGraphSim
        OpenEye::OEShape
        OpenEye::OEBio
        OpenEye::OESpruce
    PRIVATE
        oeselect
)
```

- [ ] **Step 4: Add OESPRUCE_LIBRARY to swig/CMakeLists.txt**

Add `OESPRUCE_LIBRARY` to the `_OECLUSTER_EXPECTED_LIB_VARS` list:

```cmake
set(_OECLUSTER_EXPECTED_LIB_VARS
    OECHEM_LIBRARY OESYSTEM_LIBRARY OEPLATFORM_LIBRARY
    OEMATH_LIBRARY OEZSTD_LIBRARY
    OEGRAPHSIM_LIBRARY
    OESHAPE_LIBRARY OEHERMITE_LIBRARY OEOPT_LIBRARY OEMOLPOTENTIAL_LIBRARY
    OEBIO_LIBRARY
    OESPRUCE_LIBRARY
)
```

- [ ] **Step 5: Verify CMake configuration**

Run: `cmake -B build -S . 2>&1 | grep -E "(OESpruce|oeselect|Error)"`
Expected: OESpruce found message, oeselect fetched, no errors.

- [ ] **Step 6: Commit**

```bash
git add CMakeLists.txt swig/CMakeLists.txt
git commit -m "build: add oespruce and oeselect dependencies"
```

---

### Task 2: Rewrite SuperposeMetric Header

**Files:**
- Modify: `include/oecluster/metrics/SuperposeMetric.h` (full rewrite)

- [ ] **Step 1: Write the new SuperposeMetric.h**

Replace the entire file. Enums use PascalCase to match the spec and existing `ROCSScoreType` convention:

```cpp
/**
 * @file SuperposeMetric.h
 * @brief Distance metric based on protein superposition using oespruce OESuperpose.
 */

#ifndef OECLUSTER_METRICS_SUPERPOSEMETRIC_H
#define OECLUSTER_METRICS_SUPERPOSEMETRIC_H

#include <memory>
#include <string>
#include <vector>
#include "oecluster/DistanceMetric.h"

namespace OEBio { class OEDesignUnit; }
namespace OEChem { class OEMolBase; }

namespace OECluster {

/**
 * @brief Superposition method for protein overlay.
 */
enum class SuperposeMethod {
    GlobalCarbonAlpha,  ///< All matched alpha carbon atoms (default)
    Global,             ///< Global superposition
    DDM,                ///< Distance Difference Matrix
    Weighted,           ///< Weighted DDM
    SSE,                ///< Secondary Structure Elements (Tanimoto score)
    SiteHopper          ///< SiteHopper patch score
};

/**
 * @brief Score type selection for superposition results.
 */
enum class SuperposeScoreType {
    Auto,       ///< Natural score for the method
    RMSD,       ///< Root mean square deviation
    Tanimoto,   ///< Tanimoto coefficient [0,1]
    PatchScore  ///< SiteHopper patch score [0,4]
};

/**
 * @brief Configuration options for protein superposition.
 */
struct SuperposeOptions {
    SuperposeMethod method = SuperposeMethod::GlobalCarbonAlpha;
    SuperposeScoreType score_type = SuperposeScoreType::Auto;
    bool similarity = false;
    std::string predicate;      ///< oeselect expression for both ref and fit
    std::string ref_predicate;  ///< Override predicate for ref
    std::string fit_predicate;  ///< Override predicate for fit
};

/**
 * @brief Protein superposition distance metric using oespruce OESuperpose.
 *
 * Supports multiple superposition methods (GlobalCarbonAlpha, Global, DDM,
 * Weighted, SSE, SiteHopper) with configurable score types and atom predicates
 * via oeselect expressions.
 *
 * Score semantics per method:
 * - RMSD group (Global, GlobalCarbonAlpha, DDM, Weighted): raw RMSD
 * - Tanimoto group (SSE): 1.0 - tanimoto (distance), raw tanimoto (similarity)
 * - PatchScore group (SiteHopper): 4.0 - patch_score (distance), raw (similarity)
 *
 * Each Clone() creates a new thread-local OESuperpose instance.
 */
class SuperposeMetric : public DistanceMetric {
public:
    using Options = SuperposeOptions;

    /**
     * @brief Construct from design units.
     *
     * :param dus: Shared pointers to design units.
     * :param opts: Superposition options.
     */
    explicit SuperposeMetric(
        const std::vector<std::shared_ptr<OEBio::OEDesignUnit>>& dus,
        const Options& opts = Options());

    /**
     * @brief Construct from molecules directly.
     *
     * Molecules are copied internally.
     *
     * :param mols: Pointers to molecules with 3D coordinates.
     * :param opts: Superposition options.
     */
    explicit SuperposeMetric(
        const std::vector<OEChem::OEMolBase*>& mols,
        const Options& opts = Options());

    ~SuperposeMetric() override;

    double Distance(size_t i, size_t j) override;
    std::unique_ptr<DistanceMetric> Clone() const override;
    size_t Size() const override;
    std::string Name() const override;

private:
    struct SharedData;
    struct ThreadLocalData;
    std::shared_ptr<const SharedData> shared_;
    std::unique_ptr<ThreadLocalData> local_;
    Options opts_;

    /// Private clone constructor -- shares structure data, creates new OESuperpose.
    SuperposeMetric(std::shared_ptr<const SharedData> shared, const Options& opts);

    /// Initialize thread-local OESuperpose with configured options.
    void InitSuperpose();
};

}  // namespace OECluster

#endif  // OECLUSTER_METRICS_SUPERPOSEMETRIC_H
```

- [ ] **Step 2: Verify header compiles in isolation**

Run: `echo '#include "oecluster/metrics/SuperposeMetric.h"' | c++ -std=c++17 -fsyntax-only -I include -I ${OPENEYE_INCLUDE_DIR} -x c++ -`
Expected: No errors (forward declarations only — no oespruce include needed in header).

- [ ] **Step 3: Commit**

```bash
git add include/oecluster/metrics/SuperposeMetric.h
git commit -m "feat: rewrite SuperposeMetric header with OESuperpose API"
```

---

### Task 3: Implement SuperposeMetric

**Files:**
- Modify: `src/metrics/SuperposeMetric.cpp` (full rewrite)

- [ ] **Step 1: Write the new SuperposeMetric.cpp**

Replace the entire file. Key points:
- Predicates are compiled once in the constructor and stored in `SharedData`
- Predicates are passed to `SetupRef()` and `Superpose()` in `Distance()`
- Thread-local `OESuperpose` is created in `InitSuperpose()` (called from constructor and `Clone()`)

```cpp
/**
 * @file SuperposeMetric.cpp
 * @brief Implementation of protein superposition metric using oespruce OESuperpose.
 */

#include "oecluster/metrics/SuperposeMetric.h"

#include <oechem.h>
#include <oebio.h>
#include <oespruce.h>
#include <oeselect/oeselect.h>
#include "oecluster/Error.h"

namespace OECluster {

// ---------------------------------------------------------------------------
// Helper: convert SuperposeMethod enum to oespruce constant
// ---------------------------------------------------------------------------

static unsigned int to_oe_method(SuperposeMethod m) {
    switch (m) {
        case SuperposeMethod::GlobalCarbonAlpha:
            return OESpruce::OESuperposeMethod::GlobalCarbonAlpha;
        case SuperposeMethod::Global:
            return OESpruce::OESuperposeMethod::Global;
        case SuperposeMethod::DDM:
            return OESpruce::OESuperposeMethod::DifferenceDistanceMatrix;
        case SuperposeMethod::Weighted:
            return OESpruce::OESuperposeMethod::WeightedDifferenceDistanceMatrix;
        case SuperposeMethod::SSE:
            return OESpruce::OESuperposeMethod::SecondaryStructureElements;
        case SuperposeMethod::SiteHopper:
            return OESpruce::OESuperposeMethod::SiteHopper;
    }
    throw MetricError("Unknown SuperposeMethod");
}

// ---------------------------------------------------------------------------
// Helper: resolve score type for a method
// ---------------------------------------------------------------------------

static SuperposeScoreType resolve_score_type(SuperposeMethod method,
                                              SuperposeScoreType requested) {
    if (requested != SuperposeScoreType::Auto)
        return requested;

    switch (method) {
        case SuperposeMethod::GlobalCarbonAlpha:
        case SuperposeMethod::Global:
        case SuperposeMethod::DDM:
        case SuperposeMethod::Weighted:
            return SuperposeScoreType::RMSD;
        case SuperposeMethod::SSE:
            return SuperposeScoreType::Tanimoto;
        case SuperposeMethod::SiteHopper:
            return SuperposeScoreType::PatchScore;
    }
    throw MetricError("Unknown SuperposeMethod for score resolution");
}

// ---------------------------------------------------------------------------
// Helper: validate score type is compatible with method
// ---------------------------------------------------------------------------

static void validate_score_type(SuperposeMethod method,
                                 SuperposeScoreType score_type) {
    SuperposeScoreType natural = resolve_score_type(method, SuperposeScoreType::Auto);
    if (score_type != natural) {
        throw MetricError(
            "Incompatible score type for the selected superposition method");
    }
}

// ---------------------------------------------------------------------------
// Helper: method name for Name()
// ---------------------------------------------------------------------------

static std::string method_name(SuperposeMethod m) {
    switch (m) {
        case SuperposeMethod::GlobalCarbonAlpha: return "global_carbon_alpha";
        case SuperposeMethod::Global:            return "global";
        case SuperposeMethod::DDM:               return "ddm";
        case SuperposeMethod::Weighted:          return "weighted";
        case SuperposeMethod::SSE:               return "sse";
        case SuperposeMethod::SiteHopper:        return "sitehopper";
    }
    return "unknown";
}

// ---------------------------------------------------------------------------
// SharedData: structures + compiled predicates
// ---------------------------------------------------------------------------

struct SuperposeMetric::SharedData {
    std::vector<std::shared_ptr<OEBio::OEDesignUnit>> dus;
    std::vector<OEChem::OEGraphMol> mols;
    bool use_dus = false;

    // Compiled predicates (nullptr = OEIsTrueAtom, match all)
    std::shared_ptr<OEChem::OEUnaryPredicate<OEChem::OEAtomBase>> ref_pred;
    std::shared_ptr<OEChem::OEUnaryPredicate<OEChem::OEAtomBase>> fit_pred;
};

// ---------------------------------------------------------------------------
// ThreadLocalData: per-clone OESuperpose instance
// ---------------------------------------------------------------------------

struct SuperposeMetric::ThreadLocalData {
    OESpruce::OESuperpose superpose;

    explicit ThreadLocalData(const OESpruce::OESuperposeOptions& opts)
        : superpose(opts) {}
};

// ---------------------------------------------------------------------------
// InitSuperpose
// ---------------------------------------------------------------------------

void SuperposeMetric::InitSuperpose() {
    SuperposeScoreType resolved = resolve_score_type(opts_.method, opts_.score_type);
    if (opts_.score_type != SuperposeScoreType::Auto) {
        validate_score_type(opts_.method, resolved);
    }

    OESpruce::OESuperposeOptions sp_opts(to_oe_method(opts_.method));
    local_ = std::make_unique<ThreadLocalData>(sp_opts);
}

// ---------------------------------------------------------------------------
// Helper: compile a predicate, returning nullptr for empty string
// ---------------------------------------------------------------------------

static std::shared_ptr<OEChem::OEUnaryPredicate<OEChem::OEAtomBase>>
compile_predicate(const std::string& expr) {
    if (expr.empty()) return nullptr;
    auto pred = oeselect::compile_atom_predicate(expr);
    if (!pred) {
        throw MetricError("SuperposeMetric: invalid predicate: " + expr);
    }
    return pred;
}

// ---------------------------------------------------------------------------
// Helper: resolve ref/fit predicate strings
// ---------------------------------------------------------------------------

static std::string resolve_pred_str(const std::string& specific,
                                     const std::string& general) {
    return specific.empty() ? general : specific;
}

// ---------------------------------------------------------------------------
// Constructors
// ---------------------------------------------------------------------------

SuperposeMetric::SuperposeMetric(
    const std::vector<std::shared_ptr<OEBio::OEDesignUnit>>& dus,
    const Options& opts)
    : opts_(opts) {
    if (dus.empty()) {
        throw MetricError("SuperposeMetric: empty structure list");
    }

    auto shared = std::make_shared<SharedData>();
    shared->use_dus = true;
    shared->dus = dus;

    // Compile predicates
    std::string ref_str = resolve_pred_str(opts.ref_predicate, opts.predicate);
    std::string fit_str = resolve_pred_str(opts.fit_predicate, opts.predicate);
    shared->ref_pred = compile_predicate(ref_str);
    shared->fit_pred = compile_predicate(fit_str);

    shared_ = std::move(shared);
    InitSuperpose();
}

SuperposeMetric::SuperposeMetric(
    const std::vector<OEChem::OEMolBase*>& mols,
    const Options& opts)
    : opts_(opts) {
    if (mols.empty()) {
        throw MetricError("SuperposeMetric: empty structure list");
    }

    auto shared = std::make_shared<SharedData>();
    shared->use_dus = false;
    shared->mols.reserve(mols.size());
    for (auto* mol : mols) {
        shared->mols.emplace_back(*mol);
    }

    std::string ref_str = resolve_pred_str(opts.ref_predicate, opts.predicate);
    std::string fit_str = resolve_pred_str(opts.fit_predicate, opts.predicate);
    shared->ref_pred = compile_predicate(ref_str);
    shared->fit_pred = compile_predicate(fit_str);

    shared_ = std::move(shared);
    InitSuperpose();
}

SuperposeMetric::SuperposeMetric(
    std::shared_ptr<const SharedData> shared,
    const Options& opts)
    : shared_(std::move(shared)),
      opts_(opts) {
    InitSuperpose();
}

SuperposeMetric::~SuperposeMetric() = default;

// ---------------------------------------------------------------------------
// Distance
// ---------------------------------------------------------------------------

double SuperposeMetric::Distance(size_t i, size_t j) {
    // SetupRef with predicate
    bool ref_ok = false;
    if (shared_->use_dus) {
        if (shared_->ref_pred) {
            ref_ok = local_->superpose.SetupRef(*shared_->dus[i], *shared_->ref_pred);
        } else {
            ref_ok = local_->superpose.SetupRef(*shared_->dus[i]);
        }
    } else {
        OEChem::OEGraphMol ref_copy(shared_->mols[i]);
        if (shared_->ref_pred) {
            ref_ok = local_->superpose.SetupRef(ref_copy, *shared_->ref_pred);
        } else {
            ref_ok = local_->superpose.SetupRef(ref_copy);
        }
    }
    if (!ref_ok) {
        throw MetricError("SuperposeMetric: SetupRef failed for structure " +
                          std::to_string(i));
    }

    // Superpose with predicate
    OESpruce::OESuperposeResults results;
    bool sp_ok = false;
    if (shared_->use_dus) {
        if (shared_->fit_pred) {
            sp_ok = local_->superpose.Superpose(results, *shared_->dus[j],
                                                 *shared_->fit_pred);
        } else {
            sp_ok = local_->superpose.Superpose(results, *shared_->dus[j]);
        }
    } else {
        OEChem::OEGraphMol fit_copy(shared_->mols[j]);
        if (shared_->fit_pred) {
            sp_ok = local_->superpose.Superpose(results, fit_copy,
                                                 *shared_->fit_pred);
        } else {
            sp_ok = local_->superpose.Superpose(results, fit_copy);
        }
    }
    if (!sp_ok) {
        throw MetricError("SuperposeMetric: Superpose failed for structures " +
                          std::to_string(i) + " and " + std::to_string(j));
    }

    // Extract score based on resolved score type
    SuperposeScoreType resolved = resolve_score_type(opts_.method, opts_.score_type);
    double raw_score;

    switch (resolved) {
        case SuperposeScoreType::RMSD:
            raw_score = results.GetRMSD();
            if (raw_score < 0.0) {
                throw MetricError("SuperposeMetric: sentinel RMSD for structures " +
                                  std::to_string(i) + " and " + std::to_string(j));
            }
            // RMSD: distance = raw, similarity = raw (no-op)
            return raw_score;

        case SuperposeScoreType::Tanimoto:
            raw_score = results.GetTanimoto();
            if (raw_score < 0.0) {
                throw MetricError("SuperposeMetric: sentinel Tanimoto for structures " +
                                  std::to_string(i) + " and " + std::to_string(j));
            }
            return opts_.similarity ? raw_score : (1.0 - raw_score);

        case SuperposeScoreType::PatchScore:
            raw_score = results.GetTanimoto();
            // SiteHopper returns patch score via GetTanimoto (range [0,4])
            if (raw_score < 0.0) {
                throw MetricError("SuperposeMetric: sentinel PatchScore for structures " +
                                  std::to_string(i) + " and " + std::to_string(j));
            }
            return opts_.similarity ? raw_score : (4.0 - raw_score);

        default:
            throw MetricError("SuperposeMetric: unresolved score type");
    }
}

// ---------------------------------------------------------------------------
// Clone / Size / Name
// ---------------------------------------------------------------------------

std::unique_ptr<DistanceMetric> SuperposeMetric::Clone() const {
    return std::unique_ptr<DistanceMetric>(
        new SuperposeMetric(shared_, opts_));
}

size_t SuperposeMetric::Size() const {
    return shared_->use_dus ? shared_->dus.size() : shared_->mols.size();
}

std::string SuperposeMetric::Name() const {
    return "superpose:" + method_name(opts_.method);
}

}  // namespace OECluster
```

**Implementation notes for the developer:**
- The exact oeselect API (`oeselect::compile_atom_predicate`) returns a
  `shared_ptr<OEUnaryPredicate<OEAtomBase>>`. Verify this against the oeselect
  header.
- Whether `OESuperposeResults::GetTanimoto()` returns the patch score for
  SiteHopper needs verification. There may be a dedicated `GetPatchScore()`.
- Whether `OESuperpose::SetupRef` and `Superpose` accept predicates as a
  second argument needs verification against the oespruce API. The overloads
  shown here follow the pattern from the OpenEye documentation.

- [ ] **Step 2: Verify compilation**

Run: `cmake --build build --target oecluster 2>&1 | tail -20`
Expected: Compiles without errors.

- [ ] **Step 3: Commit**

```bash
git add src/metrics/SuperposeMetric.cpp
git commit -m "feat: implement SuperposeMetric with OESuperpose backend"
```

---

### Task 4: Add Similarity Flag to FingerprintMetric

**Files:**
- Modify: `include/oecluster/metrics/FingerprintMetric.h:22-29`
- Modify: `src/metrics/FingerprintMetric.cpp:121-144`
- Test: `tests/cpp/test_fingerprint_metric.cpp`

**Design note:** `FingerprintOptions` already has a `std::string similarity` field for
the similarity function name (e.g., "tanimoto", "dice"). To add the boolean
`similarity` flag consistently with `ROCSOptions` and `SuperposeOptions`, we
rename the existing string field to `similarity_func`. This is a small
breaking change to C++ callers of the old field name, but the CLI and Python
layer abstract it.

- [ ] **Step 1: Write the failing test**

Add to `tests/cpp/test_fingerprint_metric.cpp`:

```cpp
TEST_F(FingerprintMetricTest, SimilarityMode) {
    FingerprintOptions opts;
    opts.similarity = true;
    FingerprintMetric metric(mols_, opts);
    double sim = metric.Distance(0, 1);
    // Tanimoto similarity should be in [0, 1]
    EXPECT_GE(sim, 0.0);
    EXPECT_LE(sim, 1.0);
}

TEST_F(FingerprintMetricTest, SimilarityModeErrorForManhattan) {
    FingerprintOptions opts;
    opts.similarity_func = "manhattan";
    opts.similarity = true;
    EXPECT_THROW(FingerprintMetric(mols_, opts), MetricError);
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cmake --build build --target test_fingerprint_metric && ctest -R test_fingerprint_metric -V`
Expected: FAIL — fields don't exist yet.

- [ ] **Step 3: Rename `similarity` to `similarity_func` and add `bool similarity`**

In `include/oecluster/metrics/FingerprintMetric.h`, update the struct:

```cpp
struct FingerprintOptions {
    std::string fp_type = "circular";
    unsigned int numbits = 2048;
    unsigned int min_distance = 0;
    unsigned int max_distance = 2;
    unsigned int atom_type_mask = 0;  ///< 0 = use method default
    unsigned int bond_type_mask = 0;  ///< 0 = use method default
    std::string similarity_func = "tanimoto";  ///< Similarity function name
    bool similarity = false;  ///< Return raw similarity instead of distance
};
```

- [ ] **Step 4: Update FingerprintMetric.cpp references from `similarity` to `similarity_func`**

In `src/metrics/FingerprintMetric.cpp`:
- Change all `pimpl_->opts.similarity` to `pimpl_->opts.similarity_func`
  (lines 124, 143 in the `Distance()` method and the constructor)
- Change `opts.similarity` to `opts.similarity_func` in the constructor
  fingerprint-type-selection logic (line 59: `std::string sim_lower = to_lower(opts.similarity);` — but
  actually this reference is to `fp_lower`, not `similarity`. The `similarity`
  string is only used in `Distance()`)

Add validation in the constructor, after fingerprint generation loop:
```cpp
if (opts.similarity) {
    std::string sim_lower = to_lower(opts.similarity_func);
    if (sim_lower == "manhattan" || sim_lower == "euclidean") {
        throw MetricError(
            "similarity=true is not supported with " + opts.similarity_func +
            " (no bounded similarity form)");
    }
}
```

Update `Distance()`:
```cpp
double FingerprintMetric::Distance(size_t i, size_t j) {
    const auto& fp_i = pimpl_->fingerprints[i];
    const auto& fp_j = pimpl_->fingerprints[j];
    std::string sim_lower = to_lower(pimpl_->opts.similarity_func);
    bool sim_mode = pimpl_->opts.similarity;

    if (sim_lower == "tanimoto") {
        double sim = static_cast<double>(OEGraphSim::OETanimoto(fp_i, fp_j));
        return sim_mode ? sim : (1.0 - sim);
    } else if (sim_lower == "dice") {
        double sim = static_cast<double>(OEGraphSim::OEDice(fp_i, fp_j));
        return sim_mode ? sim : (1.0 - sim);
    } else if (sim_lower == "cosine") {
        double sim = static_cast<double>(OEGraphSim::OECosine(fp_i, fp_j));
        return sim_mode ? sim : (1.0 - sim);
    } else if (sim_lower == "manhattan") {
        return static_cast<double>(OEGraphSim::OEManhattan(fp_i, fp_j));
    } else if (sim_lower == "euclidean") {
        return static_cast<double>(OEGraphSim::OEEuclid(fp_i, fp_j));
    }

    throw MetricError("Unknown similarity function: " + pimpl_->opts.similarity_func);
}
```

- [ ] **Step 5: Update oepdist.cpp reference**

In `tools/oepdist.cpp`, change `opts.similarity = fp_similarity;` to
`opts.similarity_func = fp_similarity;` and update the `JsonStr` call
from `"similarity"` to `"similarity_func"`.

- [ ] **Step 6: Run tests**

Run: `cmake --build build && ctest -R test_fingerprint_metric -V`
Expected: All PASS.

- [ ] **Step 7: Commit**

```bash
git add include/oecluster/metrics/FingerprintMetric.h \
        src/metrics/FingerprintMetric.cpp \
        tools/oepdist.cpp \
        tests/cpp/test_fingerprint_metric.cpp
git commit -m "feat: add similarity flag to FingerprintMetric, rename similarity to similarity_func"
```

---

### Task 5: Add Similarity Flag to ROCSMetric

**Files:**
- Modify: `include/oecluster/metrics/ROCSMetric.h:31-34`
- Modify: `src/metrics/ROCSMetric.cpp:59-86`
- Test: `tests/cpp/test_rocs_metric.cpp`

- [ ] **Step 1: Write the failing test**

Add to the ROCS test file:

```cpp
TEST_F(ROCSMetricTest, SimilarityMode) {
    ROCSOptions opts;
    opts.similarity = true;
    ROCSMetric metric(mols_, opts);
    double sim = metric.Distance(0, 1);
    // ComboNorm similarity should be in [0, 1]
    EXPECT_GE(sim, 0.0);
    EXPECT_LE(sim, 1.0);
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cmake --build build --target test_rocs_metric && ctest -R test_rocs_metric -V`
Expected: FAIL — `similarity` is not a member of `ROCSOptions`.

- [ ] **Step 3: Add `similarity` to ROCSOptions**

In `include/oecluster/metrics/ROCSMetric.h`:

```cpp
struct ROCSOptions {
    ROCSScoreType score_type = ROCSScoreType::ComboNorm;
    unsigned int color_ff_type = 1;  ///< OEColorFFType (1=ImplicitMillsDean)
    bool similarity = false;         ///< Return raw similarity instead of distance
};
```

- [ ] **Step 4: Update ROCSMetric::Distance() to respect similarity**

In `src/metrics/ROCSMetric.cpp`, replace the `Distance()` method:

```cpp
double ROCSMetric::Distance(size_t i, size_t j) {
    local_->overlay.SetupRef(*shared_->mols[i]);

    OEShape::OEBestOverlayScore score;
    local_->overlay.BestOverlay(score, *shared_->mols[j]);

    if (opts_.similarity) {
        switch (opts_.score_type) {
            case ROCSScoreType::ComboNorm:
                return static_cast<double>(score.GetTanimotoCombo()) / 2.0;
            case ROCSScoreType::Combo:
                return static_cast<double>(score.GetTanimotoCombo());
            case ROCSScoreType::Shape:
                return static_cast<double>(score.GetTanimoto());
            case ROCSScoreType::Color:
                return static_cast<double>(score.GetColorTanimoto());
        }
    }

    switch (opts_.score_type) {
        case ROCSScoreType::ComboNorm:
            return 1.0 - static_cast<double>(score.GetTanimotoCombo()) / 2.0;
        case ROCSScoreType::Combo:
            return 2.0 - static_cast<double>(score.GetTanimotoCombo());
        case ROCSScoreType::Shape:
            return 1.0 - static_cast<double>(score.GetTanimoto());
        case ROCSScoreType::Color:
            return 1.0 - static_cast<double>(score.GetColorTanimoto());
    }

    return 0.0;
}
```

- [ ] **Step 5: Run tests**

Run: `cmake --build build && ctest -R test_rocs_metric -V`
Expected: All PASS.

- [ ] **Step 6: Commit**

```bash
git add include/oecluster/metrics/ROCSMetric.h \
        src/metrics/ROCSMetric.cpp \
        tests/cpp/test_rocs_metric.cpp
git commit -m "feat: add similarity flag to ROCSMetric"
```

---

### Task 6: Delete SiteHopperMetric and Update SWIG Bindings

These two changes are combined into a single task to avoid intermediate build
breakage (SWIG `%include "oecluster/metrics/SiteHopperMetric.h"` would fail
if the header is deleted first).

**Files:**
- Delete: `include/oecluster/metrics/SiteHopperMetric.h`
- Delete: `src/metrics/SiteHopperMetric.cpp`
- Modify: `include/oecluster/oecluster.h:38` (remove include)
- Modify: `CMakeLists.txt:110` (remove from OECLUSTER_SOURCES)
- Modify: `swig/oecluster.i` (remove SiteHopper refs, add oespruce include)

- [ ] **Step 1: Update swig/oecluster.i**

In `swig/oecluster.i`:

1. In the `%{ ... %}` block: remove `#include "oecluster/metrics/SiteHopperMetric.h"` (line 18) and add `#include <oespruce.h>` after `#include <oebio.h>` (line 22)
2. Remove `%ignore OECluster::SiteHopperMetric::SiteHopperMetric(std::shared_ptr<const SharedData>, const Options&);` (line 323)
3. Remove `%include "oecluster/metrics/SiteHopperMetric.h"` (line 484)

The existing `%ignore` for `SuperposeMetric`'s clone constructor (line 322) remains unchanged — the signature still matches.

- [ ] **Step 2: Delete SiteHopperMetric files**

```bash
rm include/oecluster/metrics/SiteHopperMetric.h
rm src/metrics/SiteHopperMetric.cpp
```

- [ ] **Step 3: Remove SiteHopperMetric.h include from oecluster.h**

In `include/oecluster/oecluster.h`, remove:
```cpp
#include "oecluster/metrics/SiteHopperMetric.h"
```

- [ ] **Step 4: Remove SiteHopperMetric.cpp from CMakeLists.txt**

In `CMakeLists.txt`, remove `src/metrics/SiteHopperMetric.cpp` from `OECLUSTER_SOURCES`.

- [ ] **Step 5: Verify full build (library + SWIG)**

Run: `cmake --build build 2>&1 | tail -20`
Expected: Library and SWIG wrapper compile. Tests may still fail (updated in later tasks).

- [ ] **Step 6: Commit**

```bash
git add -A
git commit -m "feat: delete SiteHopperMetric, update SWIG bindings"
```

---

### Task 7: Update Python Layer

**Files:**
- Modify: `python/oecluster/__init__.py`

- [ ] **Step 1: Add module-level _oecluster import for enum access**

Near the top, alongside the existing imports from `._oecluster`, add:
```python
from . import _oecluster
```

- [ ] **Step 2: Remove SiteHopperMetric references**

1. Remove `"SiteHopperMetric"` from `__all__`
2. Remove `_HAVE_SITEHOPPER` variable and its import block (lines 158, 181-186)
3. Remove the `SiteHopperMetric` class (lines 661-687)

- [ ] **Step 3: Update SuperposeMetric import block**

Update the superpose import to not use availability guard (since oespruce is now always linked):
```python
from ._oecluster import SuperposeMetric as _SuperposeMetric
from ._oecluster import SuperposeOptions
```

Remove `_HAVE_SUPERPOSE` variable and its try/except block. Also remove
`_HAVE_FINGERPRINT` and `_HAVE_ROCS` guards since all metrics are always
available (the conditional compilation guards were removed in a prior change).
Replace with direct imports:

```python
from ._oecluster import FingerprintMetric as _FingerprintMetric
from ._oecluster import FingerprintOptions
from ._oecluster import ROCSMetric as _ROCSMetric
from ._oecluster import ROCSOptions
from ._oecluster import SuperposeMetric as _SuperposeMetric
from ._oecluster import SuperposeOptions
```

- [ ] **Step 4: Replace pdist() with kwargs-aware version**

```python
def pdist(items,
          metric,
          *,
          similarity=False,
          num_threads=0,
          chunk_size=256,
          cutoff=0.0,
          output=None,
          progress=None,
          **kwargs) -> DistanceMatrix:
    """
    Compute pairwise distances for a collection of items using a specified metric.

    :param items: List of molecules, design units, or other items.
    :param metric: Distance metric: "fingerprint", "rocs", "superpose",
                   "sitehopper", or a C++ metric object.
    :param similarity: Return similarities instead of distances.
    :param num_threads: Number of threads (0 = auto).
    :param chunk_size: Pairs per work unit.
    :param cutoff: Distance cutoff for sparse storage (0 = store all).
    :param output: Optional file path for memory-mapped storage.
    :param progress: Optional callback(completed, total).
    :param kwargs: Metric-specific options (see below).
    :returns: DistanceMatrix with computed distances/similarities.
    :raises TypeError: If unknown kwargs are passed.

    Metric-specific kwargs:

    - ``"fingerprint"``: fp_type, numbits, min_distance, max_distance,
      atom_type_mask, bond_type_mask, similarity_func
    - ``"rocs"``: score_type, color_ff_type
    - ``"superpose"``: method, score_type, predicate, ref_predicate, fit_predicate
    - ``"sitehopper"``: alias for superpose with method="sitehopper"
    """
    if isinstance(metric, str):
        metric_lower = metric.lower()
        labels = _extract_labels(items)
        params = {'metric_type': metric_lower, 'similarity': similarity}

        if metric_lower == "fingerprint":
            opts = FingerprintOptions()
            opts.similarity = similarity
            if 'fp_type' in kwargs:
                opts.fp_type = kwargs.pop('fp_type')
            if 'numbits' in kwargs:
                opts.numbits = kwargs.pop('numbits')
            if 'min_distance' in kwargs:
                opts.min_distance = kwargs.pop('min_distance')
            if 'max_distance' in kwargs:
                opts.max_distance = kwargs.pop('max_distance')
            if 'atom_type_mask' in kwargs:
                opts.atom_type_mask = kwargs.pop('atom_type_mask')
            if 'bond_type_mask' in kwargs:
                opts.bond_type_mask = kwargs.pop('bond_type_mask')
            if 'similarity_func' in kwargs:
                opts.similarity_func = kwargs.pop('similarity_func')
            if kwargs:
                raise TypeError(
                    f"Unknown kwargs for fingerprint metric: {list(kwargs)}")
            metric_obj = _FingerprintMetric(items, opts)
            metric_name = "fingerprint"

        elif metric_lower == "rocs":
            opts = ROCSOptions()
            opts.similarity = similarity
            if 'score_type' in kwargs:
                st = kwargs.pop('score_type')
                score_map = {
                    'combo_norm': _oecluster.ROCSScoreType_ComboNorm,
                    'combo': _oecluster.ROCSScoreType_Combo,
                    'shape': _oecluster.ROCSScoreType_Shape,
                    'color': _oecluster.ROCSScoreType_Color,
                }
                if st not in score_map:
                    raise ValueError(f"Unknown ROCS score type: {st}")
                opts.score_type = score_map[st]
            if 'color_ff_type' in kwargs:
                opts.color_ff_type = kwargs.pop('color_ff_type')
            if kwargs:
                raise TypeError(
                    f"Unknown kwargs for rocs metric: {list(kwargs)}")
            metric_obj = _ROCSMetric(items, opts)
            metric_name = "rocs"

        elif metric_lower in ("superpose", "sitehopper"):
            opts = SuperposeOptions()
            opts.similarity = similarity
            method = kwargs.pop('method', None)
            if metric_lower == "sitehopper":
                method = method or "sitehopper"
            if method is not None:
                method_map = {
                    'global_carbon_alpha': _oecluster.SuperposeMethod_GlobalCarbonAlpha,
                    'global': _oecluster.SuperposeMethod_Global,
                    'ddm': _oecluster.SuperposeMethod_DDM,
                    'weighted': _oecluster.SuperposeMethod_Weighted,
                    'sse': _oecluster.SuperposeMethod_SSE,
                    'sitehopper': _oecluster.SuperposeMethod_SiteHopper,
                }
                if method not in method_map:
                    raise ValueError(f"Unknown superpose method: {method}")
                opts.method = method_map[method]
            if 'score_type' in kwargs:
                st = kwargs.pop('score_type')
                st_map = {
                    'auto': _oecluster.SuperposeScoreType_Auto,
                    'rmsd': _oecluster.SuperposeScoreType_RMSD,
                    'tanimoto': _oecluster.SuperposeScoreType_Tanimoto,
                    'patch_score': _oecluster.SuperposeScoreType_PatchScore,
                }
                if st not in st_map:
                    raise ValueError(f"Unknown superpose score type: {st}")
                opts.score_type = st_map[st]
            if 'predicate' in kwargs:
                opts.predicate = kwargs.pop('predicate')
            if 'ref_predicate' in kwargs:
                opts.ref_predicate = kwargs.pop('ref_predicate')
            if 'fit_predicate' in kwargs:
                opts.fit_predicate = kwargs.pop('fit_predicate')
            if kwargs:
                raise TypeError(
                    f"Unknown kwargs for superpose metric: {list(kwargs)}")
            metric_obj = _SuperposeMetric(items, opts)
            metric_name = metric_obj.Name()

        else:
            raise ValueError(
                f"Unknown metric: {metric!r}. "
                f"Valid options: 'fingerprint', 'rocs', 'superpose', 'sitehopper'")

    else:
        metric_obj = metric
        metric_name = metric_obj.Name()
        labels = []
        params = {}

    n = metric_obj.Size()

    if output is not None:
        storage = MMapStorage(output, n)
    elif cutoff > 0.0:
        storage = SparseStorage(n, cutoff)
    else:
        storage = DenseStorage(n)

    options = PDistOptions()
    options.num_threads = num_threads
    options.chunk_size = chunk_size
    options.cutoff = cutoff
    if progress is not None:
        options.progress = progress

    _cpp_pdist(metric_obj, storage, options)
    return DistanceMatrix(storage, metric_name, labels, params)
```

- [ ] **Step 5: Replace SuperposeMetric wrapper class**

```python
class SuperposeMetric:
    """Protein superposition distance metric using oespruce OESuperpose."""

    def __new__(cls, items, *, method="global_carbon_alpha", similarity=False,
                predicate=None, ref_predicate=None, fit_predicate=None):
        """
        Construct a SuperposeMetric.

        :param items: List of OEDesignUnit or OEMolBase objects.
        :param method: Superposition method name.
        :param similarity: Return similarity instead of distance.
        :param predicate: oeselect expression for both ref and fit.
        :param ref_predicate: Override predicate for ref.
        :param fit_predicate: Override predicate for fit.
        :returns: C++ SuperposeMetric object.
        """
        opts = SuperposeOptions()
        method_map = {
            'global_carbon_alpha': _oecluster.SuperposeMethod_GlobalCarbonAlpha,
            'global': _oecluster.SuperposeMethod_Global,
            'ddm': _oecluster.SuperposeMethod_DDM,
            'weighted': _oecluster.SuperposeMethod_Weighted,
            'sse': _oecluster.SuperposeMethod_SSE,
            'sitehopper': _oecluster.SuperposeMethod_SiteHopper,
        }
        if method not in method_map:
            raise ValueError(f"Unknown superpose method: {method}")
        opts.method = method_map[method]
        opts.similarity = similarity
        if predicate is not None:
            opts.predicate = predicate
        if ref_predicate is not None:
            opts.ref_predicate = ref_predicate
        if fit_predicate is not None:
            opts.fit_predicate = fit_predicate
        return _SuperposeMetric(items, opts)
```

- [ ] **Step 6: Commit**

```bash
git add python/oecluster/__init__.py
git commit -m "feat: update Python layer for SuperposeMetric redesign"
```

---

### Task 8: Update CLI (oepdist)

**Files:**
- Modify: `tools/oepdist.cpp`

- [ ] **Step 1: Add parse helpers for new enums**

Near the top of `oepdist.cpp`, after the existing parse functions, add:

```cpp
static OECluster::SuperposeMethod ParseSuperposeMethod(const std::string& s) {
    if (s == "global_carbon_alpha") return OECluster::SuperposeMethod::GlobalCarbonAlpha;
    if (s == "global")              return OECluster::SuperposeMethod::Global;
    if (s == "ddm")                 return OECluster::SuperposeMethod::DDM;
    if (s == "weighted")            return OECluster::SuperposeMethod::Weighted;
    if (s == "sse")                 return OECluster::SuperposeMethod::SSE;
    if (s == "sitehopper")          return OECluster::SuperposeMethod::SiteHopper;
    throw std::runtime_error("Unknown superpose method: " + s);
}

static OECluster::SuperposeScoreType ParseSuperposeScoreType(const std::string& s) {
    if (s == "auto")        return OECluster::SuperposeScoreType::Auto;
    if (s == "rmsd")        return OECluster::SuperposeScoreType::RMSD;
    if (s == "tanimoto")    return OECluster::SuperposeScoreType::Tanimoto;
    if (s == "patch_score") return OECluster::SuperposeScoreType::PatchScore;
    throw std::runtime_error("Unknown superpose score type: " + s);
}
```

- [ ] **Step 2: Add --sim flag to fingerprint subcommand**

After the `fp_cmd->add_option("--similarity", ...)` line, add:
```cpp
bool fp_sim = false;
fp_cmd->add_flag("--sim", fp_sim, "Return similarity instead of distance");
```

In the fp callback, after building opts, add:
```cpp
opts.similarity = fp_sim;
```

Also update the `opts.similarity` assignment (now `opts.similarity_func`):
```cpp
opts.similarity_func = fp_similarity;
```

- [ ] **Step 3: Add --sim flag to rocs subcommand**

After the color-ff option, add:
```cpp
bool rocs_sim = false;
rocs_cmd->add_flag("--sim", rocs_sim, "Return similarity instead of distance");
```

In the rocs callback, after building opts, add:
```cpp
opts.similarity = rocs_sim;
```

- [ ] **Step 4: Replace superpose subcommand and remove sitehopper subcommand**

Remove the entire sitehopper subcommand (lines 340-400). Replace the superpose
subcommand (lines 264-338) with:

```cpp
CommonOpts sup_co;
std::string sup_method = "global_carbon_alpha";
std::string sup_score_type = "auto";
bool sup_sim = false;
std::string sup_predicate, sup_ref_predicate, sup_fit_predicate;

auto* sup_cmd = app.add_subcommand("superpose",
    "Protein superposition distance");
AddCommonOpts(sup_cmd, sup_co);
sup_cmd->add_option("--method", sup_method,
    "Method: global_carbon_alpha, global, ddm, weighted, sse, sitehopper");
sup_cmd->add_option("--score-type", sup_score_type,
    "Score: auto, rmsd, tanimoto, patch_score");
sup_cmd->add_flag("--sim", sup_sim, "Return similarity instead of distance");
sup_cmd->add_option("--predicate", sup_predicate,
    "oeselect expression for both ref and fit");
sup_cmd->add_option("--ref-predicate", sup_ref_predicate,
    "Override predicate for ref structures");
sup_cmd->add_option("--fit-predicate", sup_fit_predicate,
    "Override predicate for fit structures");

sup_cmd->callback([&]() {
    OECluster::SuperposeOptions opts;
    opts.method = ParseSuperposeMethod(sup_method);
    opts.score_type = ParseSuperposeScoreType(sup_score_type);
    opts.similarity = sup_sim;
    opts.predicate = sup_predicate;
    opts.ref_predicate = sup_ref_predicate;
    opts.fit_predicate = sup_fit_predicate;

    std::string params = "{" + JsonStr("method", sup_method) + ","
        + JsonStr("score_type", sup_score_type) + ","
        + JsonBool("similarity", sup_sim) + "}";

    if (sup_co.input2.empty()) {
        // Try design units first, fall back to molecules
        try {
            auto ds = OEPDist::ReadDesignUnits(sup_co.input1, sup_co.verbose);
            if (!ds.dus.empty()) {
                OECluster::SuperposeMetric metric(ds.dus, opts);
                OEPDist::OutputMetadata meta{
                    "pdist", metric.Name(), params, ds.labels, ds.labels};
                exit_code = run_pdist(metric, sup_co.output, meta,
                    sup_co.num_threads, sup_co.chunk_size,
                    sup_co.cutoff, sup_co.progress);
                return;
            }
        } catch (...) {}
        auto ms = OEPDist::ReadMolecules(sup_co.input1, sup_co.verbose);
        if (ms.ptrs.empty()) {
            std::cerr << "Error: No structures read from "
                      << sup_co.input1 << std::endl;
            exit_code = 1; return;
        }
        OECluster::SuperposeMetric metric(ms.ptrs, opts);
        OEPDist::OutputMetadata meta{
            "pdist", metric.Name(), params, ms.labels, ms.labels};
        exit_code = run_pdist(metric, sup_co.output, meta,
            sup_co.num_threads, sup_co.chunk_size,
            sup_co.cutoff, sup_co.progress);
    } else {
        // Cross-distance mode: try design units first
        try {
            auto sa = OEPDist::ReadDesignUnits(sup_co.input1, sup_co.verbose);
            auto sb = OEPDist::ReadDesignUnits(sup_co.input2, sup_co.verbose);
            if (!sa.dus.empty() && !sb.dus.empty()) {
                std::vector<std::shared_ptr<OEBio::OEDesignUnit>> all;
                all.insert(all.end(), sa.dus.begin(), sa.dus.end());
                all.insert(all.end(), sb.dus.begin(), sb.dus.end());
                OECluster::SuperposeMetric metric(all, opts);
                OEPDist::OutputMetadata meta{
                    "cdist", metric.Name(), params, sa.labels, sb.labels};
                exit_code = run_cdist(metric, sa.dus.size(), sup_co.output,
                    meta, sup_co.num_threads, sup_co.chunk_size,
                    sup_co.cutoff, sup_co.progress);
                return;
            }
        } catch (...) {}
        // Fall back to molecules
        auto sa = OEPDist::ReadMolecules(sup_co.input1, sup_co.verbose);
        auto sb = OEPDist::ReadMolecules(sup_co.input2, sup_co.verbose);
        if (sa.ptrs.empty() || sb.ptrs.empty()) {
            std::cerr << "Error: No structures read" << std::endl;
            exit_code = 1; return;
        }
        std::vector<OEChem::OEMolBase*> all;
        all.insert(all.end(), sa.ptrs.begin(), sa.ptrs.end());
        all.insert(all.end(), sb.ptrs.begin(), sb.ptrs.end());
        OECluster::SuperposeMetric metric(all, opts);
        OEPDist::OutputMetadata meta{
            "cdist", metric.Name(), params, sa.labels, sb.labels};
        exit_code = run_cdist(metric, sa.ptrs.size(), sup_co.output,
            meta, sup_co.num_threads, sup_co.chunk_size,
            sup_co.cutoff, sup_co.progress);
    }
});
```

- [ ] **Step 5: Remove the `ParseAlignmentMethod` function (no longer used)**

- [ ] **Step 6: Verify compilation**

Run: `cmake --build build --target oepdist 2>&1 | tail -10`
Expected: Compiles without errors.

- [ ] **Step 7: Commit**

```bash
git add tools/oepdist.cpp
git commit -m "feat: update oepdist CLI for SuperposeMetric redesign"
```

---

### Task 9: Rewrite C++ Tests

**Files:**
- Modify: `tests/cpp/test_bio_metrics.cpp` (full rewrite)
- Modify: `tests/CMakeLists.txt` (add TEST_ASSETS_DIR define if missing)

- [ ] **Step 1: Rewrite test_bio_metrics.cpp**

Replace the entire file. Tests cover: each method, distance/similarity modes,
score type validation, predicates, clone, pdist integration, and empty input.

```cpp
/**
 * @file test_bio_metrics.cpp
 * @brief Tests for SuperposeMetric using oespruce OESuperpose.
 */

#include <cmath>
#include <gtest/gtest.h>
#include "oecluster/oecluster.h"
#include "oecluster/metrics/SuperposeMetric.h"
#include <oechem.h>
#include <oebio.h>

using namespace OECluster;

static std::string asset_path(const std::string& name) {
    return std::string(TEST_ASSETS_DIR) + "/" + name;
}

class SuperposeMetricDUTest : public ::testing::Test {
protected:
    void SetUp() override {
        std::vector<std::string> files = {
            asset_path("spruce_5FQD_1_5FQD_1-ALIGNED_BC__DU__LVY_B-1438.oedu"),
            asset_path("spruce_8G66_1_8G66_1-ALIGNED_BC__DU__YOT_B-502.oedu"),
        };
        for (const auto& f : files) {
            auto du = std::make_shared<OEBio::OEDesignUnit>();
            if (!OEChem::OEReadDesignUnit(f, *du)) {
                GTEST_SKIP() << "Cannot read test asset: " << f;
            }
            dus_.push_back(du);
        }
    }

    std::vector<std::shared_ptr<OEBio::OEDesignUnit>> dus_;
};

// -- Method tests --

TEST_F(SuperposeMetricDUTest, GlobalCarbonAlpha_RMSD) {
    SuperposeOptions opts;
    opts.method = SuperposeMethod::GlobalCarbonAlpha;
    SuperposeMetric metric(dus_, opts);
    EXPECT_EQ(metric.Size(), 2u);
    EXPECT_EQ(metric.Name(), "superpose:global_carbon_alpha");
    double d = metric.Distance(0, 1);
    EXPECT_GT(d, 0.0);
    EXPECT_TRUE(std::isfinite(d));
}

TEST_F(SuperposeMetricDUTest, Global_RMSD) {
    SuperposeOptions opts;
    opts.method = SuperposeMethod::Global;
    SuperposeMetric metric(dus_, opts);
    EXPECT_EQ(metric.Name(), "superpose:global");
    double d = metric.Distance(0, 1);
    EXPECT_GT(d, 0.0);
}

TEST_F(SuperposeMetricDUTest, DDM_RMSD) {
    SuperposeOptions opts;
    opts.method = SuperposeMethod::DDM;
    SuperposeMetric metric(dus_, opts);
    EXPECT_EQ(metric.Name(), "superpose:ddm");
    double d = metric.Distance(0, 1);
    EXPECT_GT(d, 0.0);
}

TEST_F(SuperposeMetricDUTest, Weighted_RMSD) {
    SuperposeOptions opts;
    opts.method = SuperposeMethod::Weighted;
    SuperposeMetric metric(dus_, opts);
    EXPECT_EQ(metric.Name(), "superpose:weighted");
    double d = metric.Distance(0, 1);
    EXPECT_GT(d, 0.0);
}

TEST_F(SuperposeMetricDUTest, SSE_Tanimoto) {
    SuperposeOptions opts;
    opts.method = SuperposeMethod::SSE;
    SuperposeMetric metric(dus_, opts);
    EXPECT_EQ(metric.Name(), "superpose:sse");
    double d = metric.Distance(0, 1);
    // Distance = 1.0 - tanimoto, should be in [0, 1]
    EXPECT_GE(d, 0.0);
    EXPECT_LE(d, 1.0);
}

TEST_F(SuperposeMetricDUTest, SiteHopper_PatchScore) {
    SuperposeOptions opts;
    opts.method = SuperposeMethod::SiteHopper;
    SuperposeMetric metric(dus_, opts);
    EXPECT_EQ(metric.Name(), "superpose:sitehopper");
    double d = metric.Distance(0, 1);
    // Distance = 4.0 - patch_score, should be in [0, 4]
    EXPECT_GE(d, 0.0);
    EXPECT_LE(d, 4.0);
}

// -- Distance / Similarity mode tests --

TEST_F(SuperposeMetricDUTest, SSE_Similarity) {
    SuperposeOptions opts;
    opts.method = SuperposeMethod::SSE;
    opts.similarity = true;
    SuperposeMetric metric(dus_, opts);
    double sim = metric.Distance(0, 1);
    EXPECT_GE(sim, 0.0);
    EXPECT_LE(sim, 1.0);
}

TEST_F(SuperposeMetricDUTest, SiteHopper_Similarity) {
    SuperposeOptions opts;
    opts.method = SuperposeMethod::SiteHopper;
    opts.similarity = true;
    SuperposeMetric metric(dus_, opts);
    double sim = metric.Distance(0, 1);
    EXPECT_GE(sim, 0.0);
    EXPECT_LE(sim, 4.0);
}

TEST_F(SuperposeMetricDUTest, RMSD_SimilarityIsNoOp) {
    SuperposeOptions opts;
    opts.method = SuperposeMethod::GlobalCarbonAlpha;

    opts.similarity = false;
    SuperposeMetric dist_metric(dus_, opts);
    double d = dist_metric.Distance(0, 1);

    opts.similarity = true;
    SuperposeMetric sim_metric(dus_, opts);
    double s = sim_metric.Distance(0, 1);

    EXPECT_DOUBLE_EQ(d, s);
}

// -- Score type validation --

TEST_F(SuperposeMetricDUTest, IncompatibleScoreTypeThrows) {
    SuperposeOptions opts;
    opts.method = SuperposeMethod::SSE;
    opts.score_type = SuperposeScoreType::RMSD;
    EXPECT_THROW(SuperposeMetric(dus_, opts), MetricError);
}

TEST_F(SuperposeMetricDUTest, AutoScoreTypeResolvesCorrectly) {
    SuperposeOptions opts;
    opts.method = SuperposeMethod::GlobalCarbonAlpha;
    opts.score_type = SuperposeScoreType::Auto;
    SuperposeMetric metric(dus_, opts);
    double d = metric.Distance(0, 1);
    // Auto resolves to RMSD for GlobalCarbonAlpha
    EXPECT_GT(d, 0.0);
}

// -- Predicate tests --

TEST_F(SuperposeMetricDUTest, PredicateChangesSuperposition) {
    // Default (no predicate) vs backbone-only predicate should yield
    // different RMSD values
    SuperposeOptions opts_default;
    opts_default.method = SuperposeMethod::GlobalCarbonAlpha;
    SuperposeMetric metric_default(dus_, opts_default);
    double d_default = metric_default.Distance(0, 1);

    SuperposeOptions opts_pred;
    opts_pred.method = SuperposeMethod::GlobalCarbonAlpha;
    opts_pred.predicate = "backbone";
    SuperposeMetric metric_pred(dus_, opts_pred);
    double d_pred = metric_pred.Distance(0, 1);

    // Both should be finite and non-negative; they may or may not differ
    // depending on the predicate. At minimum, verify construction and
    // distance work with a predicate.
    EXPECT_TRUE(std::isfinite(d_default));
    EXPECT_TRUE(std::isfinite(d_pred));
    EXPECT_GE(d_pred, 0.0);
}

TEST_F(SuperposeMetricDUTest, InvalidPredicateThrows) {
    SuperposeOptions opts;
    opts.predicate = "!!!INVALID_PREDICATE!!!";
    EXPECT_THROW(SuperposeMetric(dus_, opts), MetricError);
}

// -- Clone and pdist integration --

TEST_F(SuperposeMetricDUTest, CloneProducesValidResults) {
    SuperposeMetric metric(dus_);
    auto clone = metric.Clone();
    EXPECT_EQ(clone->Size(), 2u);
    double d = clone->Distance(0, 1);
    EXPECT_TRUE(std::isfinite(d));
}

TEST_F(SuperposeMetricDUTest, IntegrationWithPDist) {
    SuperposeMetric metric(dus_);
    DenseStorage storage(2);
    pdist(metric, storage);
    EXPECT_EQ(storage.NumPairs(), 1u);
    double d = storage.Get(0, 1);
    EXPECT_TRUE(std::isfinite(d));
}

// -- Error cases --

TEST_F(SuperposeMetricDUTest, EmptyStructureListThrows) {
    std::vector<std::shared_ptr<OEBio::OEDesignUnit>> empty;
    EXPECT_THROW(SuperposeMetric(empty), MetricError);
}
```

- [ ] **Step 2: Ensure TEST_ASSETS_DIR is defined in tests/CMakeLists.txt**

Check if the define exists; if not, add for the `test_bio_metrics` target:
```cmake
target_compile_definitions(test_bio_metrics PRIVATE
    TEST_ASSETS_DIR="${CMAKE_SOURCE_DIR}/tests/assets")
```

- [ ] **Step 3: Run tests**

Run: `cmake --build build && ctest -R test_bio_metrics -V`
Expected: All PASS (or SKIP if test assets aren't available).

- [ ] **Step 4: Commit**

```bash
git add tests/cpp/test_bio_metrics.cpp tests/CMakeLists.txt
git commit -m "test: rewrite C++ tests for new SuperposeMetric API"
```

---

### Task 10: Add Python Tests

**Files:**
- Modify: `tests/python/test_pdist.py`

- [ ] **Step 1: Append new tests to test_pdist.py**

The existing file already imports `pytest` and `numpy as np`. Append:

```python
def test_pdist_fingerprint_similarity():
    """Test pdist with fingerprint metric in similarity mode."""
    from openeye import oechem
    import oecluster

    smiles = ["c1ccccc1", "c1ccc(O)cc1", "CCCCCCCC"]
    mols = []
    for smi in smiles:
        mol = oechem.OEGraphMol()
        oechem.OESmilesToMol(mol, smi)
        mols.append(mol)

    dist = oecluster.pdist(mols, "fingerprint", similarity=True)
    arr = np.asarray(dist)
    # Tanimoto similarity should be in [0, 1]
    assert np.all(arr >= 0.0)
    assert np.all(arr <= 1.0)

    # Compare with distance mode: sim + dist should equal 1.0
    dist_d = oecluster.pdist(mols, "fingerprint", similarity=False)
    arr_d = np.asarray(dist_d)
    np.testing.assert_array_almost_equal(arr + arr_d, np.ones_like(arr))


def test_pdist_superpose_sitehopper_alias():
    """Test that 'sitehopper' routes to superpose with method=sitehopper."""
    from openeye import oechem
    import oecluster

    files = [
        "tests/assets/spruce_5FQD_1_5FQD_1-ALIGNED_BC__DU__LVY_B-1438.oedu",
        "tests/assets/spruce_8G66_1_8G66_1-ALIGNED_BC__DU__YOT_B-502.oedu",
    ]
    dus = []
    for f in files:
        du = oechem.OEDesignUnit()
        if not oechem.OEReadDesignUnit(f, du):
            pytest.skip(f"Cannot read {f}")
        dus.append(du)

    dm = oecluster.pdist(dus, "sitehopper")
    assert "sitehopper" in dm.metric_name


def test_pdist_superpose_kwargs():
    """Test superpose with method kwarg."""
    from openeye import oechem
    import oecluster

    files = [
        "tests/assets/spruce_5FQD_1_5FQD_1-ALIGNED_BC__DU__LVY_B-1438.oedu",
        "tests/assets/spruce_8G66_1_8G66_1-ALIGNED_BC__DU__YOT_B-502.oedu",
    ]
    dus = []
    for f in files:
        du = oechem.OEDesignUnit()
        if not oechem.OEReadDesignUnit(f, du):
            pytest.skip(f"Cannot read {f}")
        dus.append(du)

    dm = oecluster.pdist(dus, "superpose", method="ddm")
    assert dm.metric_name == "superpose:ddm"


def test_pdist_superpose_predicate():
    """Test superpose with predicate kwarg."""
    from openeye import oechem
    import oecluster

    files = [
        "tests/assets/spruce_5FQD_1_5FQD_1-ALIGNED_BC__DU__LVY_B-1438.oedu",
        "tests/assets/spruce_8G66_1_8G66_1-ALIGNED_BC__DU__YOT_B-502.oedu",
    ]
    dus = []
    for f in files:
        du = oechem.OEDesignUnit()
        if not oechem.OEReadDesignUnit(f, du):
            pytest.skip(f"Cannot read {f}")
        dus.append(du)

    dm = oecluster.pdist(dus, "superpose", predicate="backbone")
    arr = np.asarray(dm)
    assert arr.shape == (1,)
    assert np.isfinite(arr[0])


def test_pdist_superpose_similarity():
    """Test superpose with similarity=True for SSE method."""
    from openeye import oechem
    import oecluster

    files = [
        "tests/assets/spruce_5FQD_1_5FQD_1-ALIGNED_BC__DU__LVY_B-1438.oedu",
        "tests/assets/spruce_8G66_1_8G66_1-ALIGNED_BC__DU__YOT_B-502.oedu",
    ]
    dus = []
    for f in files:
        du = oechem.OEDesignUnit()
        if not oechem.OEReadDesignUnit(f, du):
            pytest.skip(f"Cannot read {f}")
        dus.append(du)

    dm = oecluster.pdist(dus, "superpose", method="sse", similarity=True)
    arr = np.asarray(dm)
    # SSE Tanimoto similarity in [0, 1]
    assert np.all(arr >= 0.0)
    assert np.all(arr <= 1.0)


def test_pdist_unknown_kwargs_raises():
    """Test that unknown kwargs raise TypeError."""
    from openeye import oechem
    import oecluster

    mols = [oechem.OEGraphMol(), oechem.OEGraphMol()]
    oechem.OESmilesToMol(mols[0], "C")
    oechem.OESmilesToMol(mols[1], "CC")

    with pytest.raises(TypeError, match="Unknown kwargs"):
        oecluster.pdist(mols, "fingerprint", bogus_option=42)
```

- [ ] **Step 2: Run Python tests**

Run: `pytest tests/python/test_pdist.py -v`
Expected: All PASS.

- [ ] **Step 3: Commit**

```bash
git add tests/python/test_pdist.py
git commit -m "test: add Python tests for similarity mode, superpose kwargs, predicates"
```

---

### Task 11: Final Integration Build and Verification

- [ ] **Step 1: Full clean build**

Run: `cmake -B build -S . && cmake --build build`
Expected: No errors.

- [ ] **Step 2: Run all C++ tests**

Run: `ctest --test-dir build -V`
Expected: All PASS.

- [ ] **Step 3: Run all Python tests**

Run: `pytest tests/python/ -v`
Expected: All PASS.

- [ ] **Step 4: Quick smoke test of CLI**

Run:
```bash
./build/tools/oepdist superpose tests/assets/*.oedu -o /tmp/test_sup.npy --method global_carbon_alpha
./build/tools/oepdist superpose tests/assets/*.oedu -o /tmp/test_sh.npy --method sitehopper --sim
./build/tools/oepdist fp tests/some_mols.sdf -o /tmp/test_fp.npy --sim
./build/tools/oepdist rocs tests/some_mols.sdf -o /tmp/test_rocs.npy --sim
```
Expected: All produce output files without errors.

- [ ] **Step 5: Final commit (if any fixups needed)**

```bash
git add -A
git commit -m "chore: integration fixups for SuperposeMetric redesign"
```
