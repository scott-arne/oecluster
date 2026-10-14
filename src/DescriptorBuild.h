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

namespace OEChem { class OEMolBase; }

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
 * Column and group names are matched exactly, unlike the source names
 * ``make_descriptor_calculator`` accepts, which are case-insensitive.
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

/**
 * @brief Name why a column's variance disqualifies it, or "" if it does not.
 *
 * The three causes need different responses from a caller, and one string
 * for all of them told a caller their column was constant when in truth it
 * had too few values to measure. Order is load-bearing: -inf is both
 * non-finite and negative, so present_count is asked about first and
 * finiteness second.
 *
 * :param variance: The column variance as OEFP reported it.
 * :param present_count: Rows in which the column had a value.
 * :returns: ``"insufficient-data"``, ``"non-finite-variance"``,
 *     ``"zero-variance"``, or ``""`` when the column is usable.
 */
std::string variance_drop_reason(double variance, size_t present_count);

/**
 * @brief Convert a string to lowercase for case-insensitive matching.
 *
 * :param value: Input string.
 * :returns: Lowercase copy of the input.
 */
std::string to_lower(const std::string& value);

/**
 * @brief Validate molecule pointers and convert to const pointers.
 *
 * :param mols: Input molecule pointers.
 * :param caller: Caller name for error messages.
 * :returns: Const-qualified molecule pointers.
 * :raises ComparisonError: When a null pointer is found.
 */
std::vector<const OEChem::OEMolBase*> checked_inputs(
    const std::vector<OEChem::OEMolBase*>& mols, const char* caller);

}  // namespace OECluster

#endif  // OECLUSTER_DESCRIPTORBUILD_H
