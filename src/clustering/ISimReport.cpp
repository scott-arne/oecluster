/**
 * @file ISimReport.cpp
 * @brief iSIM set similarity and the approximate fingerprint-native cluster report.
 */

#include "oecluster/clustering/ISimReport.h"

#include <limits>
#include <stdexcept>
#include <string>

#include "ISimKernels.h"

namespace OECluster {

namespace {

double nan_value() {
    return std::numeric_limits<double>::quiet_NaN();
}

// Shared by both entry points so they refuse the same inputs in the same
// order: metric, then zero width (as bitbirch refuses it), then width, then
// batch size.
void validate_isim_input(const std::string& metric, const OEFP::OEFPBatch& fingerprints,
                         const char* caller) {
    if (metric != "tanimoto") {
        throw std::invalid_argument(std::string(caller) +
                                    ": metric must be \"tanimoto\", got \"" + metric + "\"");
    }
    if (fingerprints.Size() > 0 && fingerprints.SizeBits() == 0) {
        throw std::invalid_argument(std::string(caller) +
                                    " requires non-zero-width fingerprints");
    }
    detail::check_isim_width(fingerprints.SizeBits());
    detail::check_isim_batch_size(fingerprints.Size());
}

}  // namespace

double isim(const OEFP::OEFPBatch& fingerprints, const ISimOptions& options) {
    validate_isim_input(options.metric, fingerprints, "isim");
    const size_t num_fingerprints = fingerprints.Size();
    if (num_fingerprints < 2) {
        return nan_value();
    }
    const size_t size_bits = fingerprints.SizeBits();
    detail::BitCounts counts(size_bits, 0u);
    for (size_t row = 0; row < num_fingerprints; ++row) {
        detail::add_row_to_counts(counts, fingerprints.RowWords(row), size_bits);
    }
    const detail::ISimSums sums =
        detail::isim_sums(num_fingerprints, detail::count_moments(counts));
    return detail::isim_ratio(sums.intersections, sums.unions);
}

}  // namespace OECluster
