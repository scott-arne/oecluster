#ifndef OEPDIST_OUTPUTWRITER_H
#define OEPDIST_OUTPUTWRITER_H

#include <cstddef>
#include <string>
#include <vector>

namespace OEPDist {

struct OutputMetadata {
    std::string mode;
    std::string metric;
    std::string params_json;
    std::vector<std::string> row_labels;
    std::vector<std::string> col_labels;
};

void WritePDist(const std::string& output_path,
                const double* data,
                size_t n,
                const OutputMetadata& meta);

void WriteCDist(const std::string& output_path,
                const double* data,
                size_t n_rows, size_t n_cols,
                const OutputMetadata& meta);

}  // namespace OEPDist

#endif  // OEPDIST_OUTPUTWRITER_H
