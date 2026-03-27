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
// SharedData: structures + parsed selections
// ---------------------------------------------------------------------------

struct SuperposeMetric::SharedData {
    std::vector<std::shared_ptr<OEBio::OEDesignUnit>> dus;
    std::vector<OEChem::OEGraphMol> mols;
    bool use_dus = false;

    // Parsed oeselect selections (empty OESelection = match all atoms)
    OESel::OESelection ref_sele;
    OESel::OESelection fit_sele;
    bool has_ref_pred = false;
    bool has_fit_pred = false;
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
// Helper: parse a predicate string, returning a validated OESelection
// ---------------------------------------------------------------------------

static OESel::OESelection parse_predicate(const std::string& expr) {
    try {
        return OESel::OESelection::Parse(expr);
    } catch (const std::exception& e) {
        throw MetricError("SuperposeMetric: invalid predicate: " + expr
                          + " (" + e.what() + ")");
    }
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

static void init_predicates(OESel::OESelection& ref_sele, bool& has_ref,
                             OESel::OESelection& fit_sele, bool& has_fit,
                             const SuperposeOptions& opts) {
    std::string ref_str = resolve_pred_str(opts.ref_predicate, opts.predicate);
    std::string fit_str = resolve_pred_str(opts.fit_predicate, opts.predicate);
    if (!ref_str.empty()) {
        ref_sele = parse_predicate(ref_str);
        has_ref = true;
    }
    if (!fit_str.empty()) {
        fit_sele = parse_predicate(fit_str);
        has_fit = true;
    }
}

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
    init_predicates(shared->ref_sele, shared->has_ref_pred,
                    shared->fit_sele, shared->has_fit_pred, opts);

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
    init_predicates(shared->ref_sele, shared->has_ref_pred,
                    shared->fit_sele, shared->has_fit_pred, opts);

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
    // When predicates are used with design units, extract the protein and use
    // the molecule-based API so OESelect is bound to the same molecule instance
    // that SetupRef/Superpose operates on.
    bool ref_ok = false;
    if (shared_->use_dus) {
        if (shared_->has_ref_pred) {
            OEChem::OEGraphMol ref_prot;
            shared_->dus[i]->GetProtein(ref_prot);
            OESel::OESelect ref_pred(ref_prot, shared_->ref_sele);
            ref_ok = local_->superpose.SetupRef(ref_prot, ref_pred);
        } else {
            ref_ok = local_->superpose.SetupRef(*shared_->dus[i]);
        }
    } else {
        OEChem::OEGraphMol ref_copy(shared_->mols[i]);
        if (shared_->has_ref_pred) {
            OESel::OESelect ref_pred(ref_copy, shared_->ref_sele);
            ref_ok = local_->superpose.SetupRef(ref_copy, ref_pred);
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
        if (shared_->has_fit_pred) {
            OEChem::OEGraphMol fit_prot;
            shared_->dus[j]->GetProtein(fit_prot);
            OESel::OESelect fit_pred(fit_prot, shared_->fit_sele);
            sp_ok = local_->superpose.Superpose(results, fit_prot,
                                                 fit_pred);
        } else {
            sp_ok = local_->superpose.Superpose(results, *shared_->dus[j]);
        }
    } else {
        OEChem::OEGraphMol fit_copy(shared_->mols[j]);
        if (shared_->has_fit_pred) {
            OESel::OESelect fit_pred(fit_copy, shared_->fit_sele);
            sp_ok = local_->superpose.Superpose(results, fit_copy,
                                                 fit_pred);
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
