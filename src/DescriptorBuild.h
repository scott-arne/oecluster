/**
 * @file DescriptorBuild.h
 * @brief Shared construction of OEFP descriptor calculators and column selections.
 *
 * Private to the descriptor implementations: not installed, not exposed through SWIG.
 */

#ifndef OECLUSTER_DESCRIPTORBUILD_H
#define OECLUSTER_DESCRIPTORBUILD_H

#include <oefp/descriptor_calculator.h>
#include <oefp/descriptor_schema.h>
#include <oefp/descriptor_source.h>
#include <cstddef>
#include <memory>
#include <string>
#include <vector>

namespace OECluster {

/**
 * @brief Build a descriptor calculator over the named sources.
 *
 * :param sources: Source names; an empty list selects ``{"openeye"}``.
 * :returns: A calculator over the merged, deduplicated schema.
 * :raises ComparisonError: When a source name is not recognized or OEFP rejects
 *     the merged schema.
 */
std::shared_ptr<const OEFP::DescriptorCalculator> make_descriptor_calculator(
    const std::vector<std::string>& sources);

/**
 * @brief Resolve explicit column names and group names to schema indices.
 *
 * Names and groups are unioned, then deduplicated and returned in schema order,
 * so a column named both directly and through its group appears once. An empty
 * ``columns`` and empty ``groups`` selects every column in the schema.
 *
 * :param schema: The calculator's merged schema.
 * :param columns: Explicit column names.
 * :param groups: Group names.
 * :returns: Schema indices in ascending order.
 * :raises ComparisonError: When a name or group is not present in the schema.
 */
std::vector<size_t> resolve_column_indices(const OEFP::DescriptorSchema& schema,
                                           const std::vector<std::string>& columns,
                                           const std::vector<std::string>& groups);

/**
 * @brief Report whether a descriptor value kind is a numeric scalar.
 *
 * :param kind: Schema value kind.
 * :returns: ``true`` for Bool, Int, and Float.
 */
bool is_numeric_kind(OEFP::DescriptorValueKind kind);

}  // namespace OECluster

#endif  // OECLUSTER_DESCRIPTORBUILD_H
