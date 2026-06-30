/**
 * @file OutputWriter.h
 * @brief Distance matrix output writers (NumPy .npy, CSV, binary) for oepdist CLI.
 */

#ifndef OEPDIST_OUTPUTWRITER_H
#define OEPDIST_OUTPUTWRITER_H

#include <cstddef>
#include <string>
#include <vector>

namespace OEPDist {

/**
 * @brief Metadata for distance matrix output, including mode, method, parameters, and labels.
 */
struct OutputMetadata {
    std::string mode;
    std::string comparison;
    std::string params_json;
    std::vector<std::string> row_labels;
    std::vector<std::string> col_labels;
};

/**
 * @brief Writes a pairwise distance matrix in condensed form.
 *
 * :param output_path: Output file path (.npy, .csv, or .bin).
 * :param data: Condensed distance array (n*(n-1)/2 elements).
 * :param n: Number of items.
 * :param meta: Output metadata including labels and parameters.
 * :raises std::runtime_error: if file write fails or format is unsupported.
 */
void WritePDist(const std::string& output_path,
                const double* data,
                size_t n,
                const OutputMetadata& meta);

/**
 * @brief Writes a cross-distance matrix (rectangular).
 *
 * :param output_path: Output file path (.npy, .csv, or .bin).
 * :param data: Full distance matrix (n_rows * n_cols elements, row-major).
 * :param n_rows: Number of rows.
 * :param n_cols: Number of columns.
 * :param meta: Output metadata including labels and parameters.
 * :raises std::runtime_error: if file write fails or format is unsupported.
 */
void WriteCDist(const std::string& output_path,
                const double* data,
                size_t n_rows, size_t n_cols,
                const OutputMetadata& meta);

}  // namespace OEPDist

#endif  // OEPDIST_OUTPUTWRITER_H
