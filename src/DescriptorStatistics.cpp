/**
 * @file DescriptorStatistics.cpp
 * @brief Implementation of fitted descriptor statistics.
 */

#include "oecluster/DescriptorStatistics.h"

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
    // single surviving column. Without this guard it would fail below as a
    // selection emptied by insufficient data, which is true of every column
    // at once and so points at the descriptors rather than at the one count
    // the caller can change.
    if (mols.size() < 2) {
        throw ComparisonError(
            "descriptor_statistics requires at least two molecules to fit variances, got " +
            std::to_string(mols.size()));
    }

    const std::vector<const OEChem::OEMolBase*> inputs =
        checked_inputs(mols, "descriptor_statistics");

    const std::shared_ptr<const OEFP::DescriptorCalculator> calculator =
        make_descriptor_calculator(options.sources);
    const OEFP::DescriptorSchema& schema = calculator->Schema();

    const std::vector<size_t> requested =
        resolve_column_indices(schema, options.columns, options.groups);

    // Non-numeric columns cannot be summarized, so they are dropped before OEFP
    // is asked for statistics rather than causing it to throw.
    //
    // The branch is reserved, not live, exactly as its counterpart in
    // numeric_selection is. oecluster counts Bool, Int and Float as numeric,
    // and every source make_descriptor_calculator can build emits only those
    // kinds, so no schema reaches the branch and no test can drive it -- it is
    // untested because it is unreachable, not because it was overlooked. It is
    // kept against OEFP adding a source that emits another kind, where
    // rebuilding the handling would cost more than carrying it.
    //
    // Which kinds it would then cover is worth recording, because the two
    // differ. A string-valued column is handled here as written: OEFP counts
    // String as scalar, so the batch builds and this filter drops it. A
    // vector, matrix or fingerprint column is NOT -- CalculateBatch builds a
    // batch over the whole merged schema and rejects any non-scalar column
    // before this filter can matter, so such a source would need the
    // calculator narrowed to the numeric columns first.
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
    // Reserved along with the filter above, and unreachable for the same
    // reason: emptying the selection here takes a column the schema declares
    // and this layer does not count as numeric.
    //
    // Kept apart from the all-dropped refusal below, which says the columns
    // were dropped "over these molecules": nothing has been computed over
    // them yet here, and the comparison surface words this same verdict the
    // same way from a validator that holds no molecules at all.
    if (numeric.empty()) {
        throw ComparisonError("Descriptor selection resolved to zero numeric columns");
    }

    const OEFP::DescriptorBatch batch = [&]() {
        try {
            return detail::calculate_descriptor_batch(*calculator, inputs);
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
        const size_t present = static_cast<size_t>(raw.present_count[k]);
        const std::string reason = variance_drop_reason(variance, present);
        if (!reason.empty()) {
            result.dropped_columns.push_back(name);
            result.dropped_reasons.push_back(reason);
            continue;
        }
        surviving.push_back(numeric[k]);
        result.columns.push_back(name);
        result.mean.push_back(raw.mean[k]);
        result.variance.push_back(variance);
        result.minimum.push_back(raw.minimum[k]);
        result.maximum.push_back(raw.maximum[k]);
        result.present_count.push_back(present);
    }

    if (surviving.empty()) {
        // The caller never receives the drop report in this case, so the
        // causes have to travel in the refusal itself. Naming one cause for
        // all of them told a caller their columns were constant when the
        // truth was that none had enough values to measure.
        //
        // Built from result.dropped_reasons, the one ledger both filters
        // write to, rather than from a second vector collected here. A
        // message that claims to say why every selected column was dropped
        // has to read every drop, and the non-numeric filter above records
        // its own: a selection half emptied there and half here would
        // otherwise be explained by half its causes.
        throw ComparisonError("Every selected descriptor column was dropped over these "
                              "molecules (" +
                              join_drop_reasons(result.dropped_reasons) + ")");
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
