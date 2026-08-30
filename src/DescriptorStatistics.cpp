/**
 * @file DescriptorStatistics.cpp
 * @brief Implementation of fitted descriptor statistics.
 */

#include "oecluster/DescriptorStatistics.h"

#include <cmath>
#include <exception>
#include <oechem.h>
#include <oefp/descriptor_batch.h>
#include <oefp/descriptor_selection.h>
#include <oefp/descriptor_statistics.h>
#include "DescriptorBuild.h"
#include "oecluster/Error.h"

namespace OECluster {

DescriptorStatisticsResult descriptor_statistics(const std::vector<OEChem::OEMolBase*>& mols,
                                                 const DescriptorStatisticsOptions& options) {
    // Variance needs two observations, so a shorter input cannot produce a
    // single surviving column and would otherwise fail below as though every
    // descriptor were constant.
    if (mols.size() < 2) {
        throw ComparisonError(
            "descriptor_statistics requires at least two molecules to fit variances, got " +
            std::to_string(mols.size()));
    }

    std::vector<const OEChem::OEMolBase*> inputs;
    inputs.reserve(mols.size());
    for (size_t i = 0; i < mols.size(); ++i) {
        if (mols[i] == nullptr) {
            throw ComparisonError("descriptor_statistics received null molecule pointer at index " +
                                  std::to_string(i));
        }
        inputs.push_back(mols[i]);
    }

    const std::shared_ptr<const OEFP::DescriptorCalculator> calculator =
        make_descriptor_calculator(options.sources);
    const OEFP::DescriptorSchema& schema = calculator->Schema();

    const std::vector<size_t> requested =
        resolve_column_indices(schema, options.columns, options.groups);

    // Non-numeric columns cannot be summarized, so they are dropped before OEFP
    // is asked for statistics rather than causing it to throw. No source OEFP
    // 0.3.0 ships can reach this branch -- openeye, rdkit and mordred emit only
    // float, int and bool. A future string-valued column is handled here as
    // written: OEFP counts String as scalar, so the batch builds and the filter
    // drops it. A vector, matrix or fingerprint column is NOT handled here --
    // CalculateBatch builds a batch over the whole merged schema and rejects
    // any non-scalar column before this filter can matter, so such a source
    // would need the calculator narrowed to the numeric columns first.
    DescriptorStatisticsResult result;
    std::vector<size_t> numeric;
    for (size_t index : requested) {
        const OEFP::DescriptorDefinition& definition = schema.Definition(index);
        if (is_numeric_kind(definition.value_kind)) {
            numeric.push_back(index);
        } else {
            result.dropped_columns.push_back(definition.name);
            result.dropped_reasons.push_back("non-numeric");
        }
    }
    if (numeric.empty()) {
        throw ComparisonError("Descriptor selection resolved to zero numeric columns");
    }

    const OEFP::DescriptorBatch batch = [&]() {
        try {
            return calculator->CalculateBatch(inputs);
        } catch (const std::exception& exc) {
            throw ComparisonError("Failed to compute descriptors: " + std::string(exc.what()));
        }
    }();
    result.num_rows = inputs.size();

    const OEFP::DescriptorColumnStatistics raw = [&]() {
        try {
            return OEFP::ColumnStatistics(batch, OEFP::DescriptorSelection::Indices(numeric));
        } catch (const std::exception& exc) {
            throw ComparisonError("Failed to compute descriptor statistics: " +
                                  std::string(exc.what()));
        }
    }();

    // The survival loop below indexes raw's vectors by position within numeric.
    // That is safe because OEFP assigns names, mean, variance, minimum, maximum
    // and present_count to one common length in a single block, so checking one
    // of them checks all six. Pin the count rather than reading past the end if
    // that ever stops holding.
    if (raw.names.size() != numeric.size()) {
        throw ComparisonError("OEFP returned statistics for " +
                              std::to_string(raw.names.size()) + " columns, expected " +
                              std::to_string(numeric.size()));
    }

    std::vector<size_t> surviving;
    for (size_t k = 0; k < numeric.size(); ++k) {
        const std::string& name = schema.Definition(numeric[k]).name;
        const double variance = raw.variance[k];
        if (!std::isfinite(variance) || variance <= 0.0) {
            result.dropped_columns.push_back(name);
            result.dropped_reasons.push_back("zero-variance");
            continue;
        }
        surviving.push_back(numeric[k]);
        result.columns.push_back(name);
        result.mean.push_back(raw.mean[k]);
        result.variance.push_back(variance);
        result.minimum.push_back(raw.minimum[k]);
        result.maximum.push_back(raw.maximum[k]);
        result.present_count.push_back(static_cast<size_t>(raw.present_count[k]));
    }

    if (surviving.empty()) {
        throw ComparisonError(
            "Every selected descriptor column has zero variance over these molecules");
    }

    if (options.inverse_covariance) {
        try {
            const OEFP::DescriptorInverseCovariance inverse = OEFP::InverseCovarianceMatrix(
                batch, OEFP::DescriptorSelection::Indices(surviving));
            result.inverse_covariance = inverse.matrix;
            result.inverse_covariance_rank = inverse.rank;
            result.inverse_covariance_rows = static_cast<size_t>(inverse.row_count);
        } catch (const std::exception& exc) {
            throw ComparisonError("Failed to compute descriptor inverse covariance: " +
                                  std::string(exc.what()));
        }
    }

    return result;
}

}  // namespace OECluster
