/**
 * @file DescriptorBuild.cpp
 * @brief Implementation of shared descriptor calculator and selection construction.
 */

#include "DescriptorBuild.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <exception>
#include <iterator>
#include <string>
#include <oechem.h>
#include "oecluster/Error.h"

namespace OECluster {

std::string to_lower(const std::string& value) {
    std::string result = value;
    std::transform(result.begin(), result.end(), result.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return result;
}

std::vector<const OEChem::OEMolBase*> checked_inputs(
        const std::vector<OEChem::OEMolBase*>& mols, const char* caller) {
    std::vector<const OEChem::OEMolBase*> inputs;
    inputs.reserve(mols.size());
    for (size_t i = 0; i < mols.size(); ++i) {
        if (mols[i] == nullptr) {
            throw ComparisonError(std::string(caller) +
                                  " received null molecule pointer at index " +
                                  std::to_string(i));
        }
        inputs.push_back(mols[i]);
    }
    return inputs;
}

namespace {

std::shared_ptr<const OEFP::DescriptorSource> make_source(const std::string& name) {
    const std::string key = to_lower(name);
    if (key == "openeye") {
        return std::make_shared<const OEFP::OpenEyePropertyDescriptorSource>();
    }
    if (key == "mordred") {
        return std::make_shared<const OEFP::MordredDescriptorSource>();
    }
    if (key == "rdkit") {
        return std::make_shared<const OEFP::RDKitDescriptorSource>();
    }
    throw ComparisonError("Unknown descriptor source: " + name +
                          ". Supported sources are 'openeye', 'mordred', 'rdkit'");
}

}  // namespace

std::shared_ptr<const OEFP::DescriptorCalculator> make_descriptor_calculator(
        const std::vector<std::string>& sources) {
    std::vector<std::string> names = sources;
    if (names.empty()) {
        names.push_back("openeye");
    }

    std::vector<OEFP::DescriptorSourceEntry> entries;
    entries.reserve(names.size());
    for (const std::string& name : names) {
        entries.emplace_back(make_source(name));
    }

    try {
        return std::make_shared<const OEFP::DescriptorCalculator>(std::move(entries));
    } catch (const std::exception& exc) {
        throw ComparisonError("Failed to build descriptor calculator: " +
                              std::string(exc.what()));
    }
}

std::vector<size_t> resolve_column_indices(const OEFP::DescriptorSchema& schema,
                                           const std::vector<std::string>& columns,
                                           const std::vector<std::string>& groups) {
    std::vector<size_t> indices;

    if (columns.empty() && groups.empty()) {
        indices.reserve(schema.Size());
        for (size_t i = 0; i < schema.Size(); ++i) {
            indices.push_back(i);
        }
        return indices;
    }

    for (const std::string& name : columns) {
        if (!schema.Contains(name)) {
            throw ComparisonError("Unknown descriptor column: " + name);
        }
        indices.push_back(schema.IndexOf(name));
    }

    for (const std::string& group : groups) {
        const std::vector<size_t> group_indices = schema.IndicesForGroup(group);
        if (group_indices.empty()) {
            throw ComparisonError("Unknown or empty descriptor group: " + group);
        }
        indices.insert(indices.end(), group_indices.begin(), group_indices.end());
    }

    std::sort(indices.begin(), indices.end());
    indices.erase(std::unique(indices.begin(), indices.end()), indices.end());
    return indices;
}

bool is_numeric_kind(OEFP::DescriptorValueKind kind) {
    return kind == OEFP::DescriptorValueKind::Bool || kind == OEFP::DescriptorValueKind::Int ||
           kind == OEFP::DescriptorValueKind::Float;
}

std::string variance_drop_reason(double variance, size_t present_count) {
    // OEFP leaves the variance at its NaN initialization when it never
    // computed one, which is every column present in fewer than two rows.
    // That is a different problem from a column that is genuinely constant:
    // the remedy is more molecules where the descriptor is defined, not
    // accepting that the column carries no signal.
    if (present_count < 2u) {
        return "insufficient-data";
    }
    // Present NaNs and an accumulation that overflowed both land here, and
    // the variance alone cannot tell them apart -- so the name says the
    // variance is not finite and claims nothing about the values.
    if (!std::isfinite(variance)) {
        return "non-finite-variance";
    }
    if (variance <= 0.0) {
        return "zero-variance";
    }
    return {};
}

std::string join_drop_reasons(const std::vector<std::string>& reasons) {
    static const std::vector<std::string> CANONICAL_ORDER = {
        "insufficient-data", "non-finite-variance", "zero-variance"};

    std::vector<std::string> distinct;
    for (const std::string& reason : reasons) {
        if (std::find(distinct.begin(), distinct.end(), reason) == distinct.end()) {
            distinct.push_back(reason);
        }
    }

    // An unknown reason lands at CANONICAL_ORDER.size(), which ranks it last,
    // and the sort is stable so several of them keep the order they were
    // collected in instead of being reordered arbitrarily against each other.
    const auto rank = [](const std::string& reason) {
        return std::distance(CANONICAL_ORDER.begin(),
                             std::find(CANONICAL_ORDER.begin(), CANONICAL_ORDER.end(), reason));
    };
    std::stable_sort(distinct.begin(), distinct.end(),
                     [&rank](const std::string& left, const std::string& right) {
                         return rank(left) < rank(right);
                     });

    std::string joined;
    for (const std::string& reason : distinct) {
        if (!joined.empty()) {
            joined += ", ";
        }
        joined += reason;
    }
    return joined;
}

}  // namespace OECluster
