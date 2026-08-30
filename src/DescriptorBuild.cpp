/**
 * @file DescriptorBuild.cpp
 * @brief Implementation of shared descriptor calculator and selection construction.
 */

#include "DescriptorBuild.h"

#include <algorithm>
#include <cctype>
#include <exception>
#include "oecluster/Error.h"

namespace OECluster {

namespace {

std::string to_lower(const std::string& value) {
    std::string result = value;
    std::transform(result.begin(), result.end(), result.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return result;
}

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

}  // namespace OECluster
