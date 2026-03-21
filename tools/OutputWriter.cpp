#include "OutputWriter.h"

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <stdexcept>

namespace OEPDist {

namespace {

std::string GetExtension(const std::string& path) {
    auto dot = path.rfind('.');
    if (dot == std::string::npos) return "";
    std::string ext = path.substr(dot);
    std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);
    return ext;
}

std::string SidecarPath(const std::string& path) {
    auto dot = path.rfind('.');
    if (dot == std::string::npos) return path + ".json";
    return path.substr(0, dot) + ".json";
}

std::string JsonEscape(const std::string& s) {
    std::string out;
    for (char c : s) {
        if (c == '"') out += "\\\"";
        else if (c == '\\') out += "\\\\";
        else out += c;
    }
    return out;
}

void WriteSidecar(const std::string& path,
                  const OutputMetadata& meta,
                  size_t n_rows, size_t n_cols) {
    std::ofstream f(path);
    if (!f) throw std::runtime_error("Failed to write: " + path);
    f << "{\n";
    f << "  \"mode\": \"" << meta.mode << "\",\n";
    f << "  \"metric\": \"" << meta.metric << "\",\n";
    f << "  \"params\": " << meta.params_json << ",\n";
    f << "  \"n_rows\": " << n_rows << ",\n";
    f << "  \"n_cols\": " << n_cols << ",\n";
    f << "  \"row_labels\": [";
    for (size_t i = 0; i < meta.row_labels.size(); ++i) {
        if (i > 0) f << ", ";
        f << "\"" << JsonEscape(meta.row_labels[i]) << "\"";
    }
    f << "],\n";
    f << "  \"col_labels\": [";
    for (size_t i = 0; i < meta.col_labels.size(); ++i) {
        if (i > 0) f << ", ";
        f << "\"" << JsonEscape(meta.col_labels[i]) << "\"";
    }
    f << "]\n}\n";
}

void WriteNpy(const std::string& path, const double* data,
              size_t n_elements, const std::string& shape_str) {
    std::ofstream f(path, std::ios::binary);
    if (!f) throw std::runtime_error("Failed to write: " + path);

    std::string header = "{'descr': '<f8', 'fortran_order': False, 'shape': "
                         + shape_str + ", }";

    size_t preamble = 10;
    size_t total = preamble + header.size() + 1;
    size_t padding = (64 - (total % 64)) % 64;
    header.append(padding, ' ');
    header += '\n';

    uint16_t header_len = static_cast<uint16_t>(header.size());

    const char magic[] = "\x93NUMPY";
    f.write(magic, 6);
    const char version[] = "\x01\x00";
    f.write(version, 2);
    f.write(reinterpret_cast<const char*>(&header_len), 2);
    f.write(header.data(), header.size());
    f.write(reinterpret_cast<const char*>(data),
            n_elements * sizeof(double));
}

void WritePDistCsv(const std::string& path, const double* data,
                   size_t n, const OutputMetadata& meta) {
    std::ofstream f(path);
    if (!f) throw std::runtime_error("Failed to write: " + path);

    f << "";
    for (size_t j = 0; j < n; ++j)
        f << "," << meta.col_labels[j];
    f << "\n";

    char buf[32];
    for (size_t i = 0; i < n; ++i) {
        f << meta.row_labels[i];
        for (size_t j = 0; j < n; ++j) {
            f << ",";
            if (i == j) {
                f << "0";
            } else {
                size_t a = std::min(i, j), b = std::max(i, j);
                size_t k = n * a - a * (a + 1) / 2 + b - a - 1;
                std::snprintf(buf, sizeof(buf), "%.8g", data[k]);
                f << buf;
            }
        }
        f << "\n";
    }
}

void WriteCDistCsv(const std::string& path, const double* data,
                   size_t n_rows, size_t n_cols,
                   const OutputMetadata& meta) {
    std::ofstream f(path);
    if (!f) throw std::runtime_error("Failed to write: " + path);

    f << "";
    for (size_t j = 0; j < n_cols; ++j)
        f << "," << meta.col_labels[j];
    f << "\n";

    char buf[32];
    for (size_t i = 0; i < n_rows; ++i) {
        f << meta.row_labels[i];
        for (size_t j = 0; j < n_cols; ++j) {
            std::snprintf(buf, sizeof(buf), "%.8g", data[i * n_cols + j]);
            f << "," << buf;
        }
        f << "\n";
    }
}

void WriteBin(const std::string& path, const double* data,
              size_t n_elements) {
    std::ofstream f(path, std::ios::binary);
    if (!f) throw std::runtime_error("Failed to write: " + path);
    f.write(reinterpret_cast<const char*>(data),
            n_elements * sizeof(double));
}

}  // namespace

void WritePDist(const std::string& output_path,
                const double* data, size_t n,
                const OutputMetadata& meta) {
    size_t n_pairs = n * (n - 1) / 2;
    std::string ext = GetExtension(output_path);

    if (ext == ".npy") {
        std::string shape = "(" + std::to_string(n_pairs) + ",)";
        WriteNpy(output_path, data, n_pairs, shape);
        WriteSidecar(SidecarPath(output_path), meta, n, n);
    } else if (ext == ".csv") {
        WritePDistCsv(output_path, data, n, meta);
    } else if (ext == ".bin") {
        WriteBin(output_path, data, n_pairs);
        WriteSidecar(SidecarPath(output_path), meta, n, n);
    } else {
        throw std::runtime_error(
            "Unknown output format: " + ext +
            " (use .npy, .csv, or .bin)");
    }
}

void WriteCDist(const std::string& output_path,
                const double* data,
                size_t n_rows, size_t n_cols,
                const OutputMetadata& meta) {
    std::string ext = GetExtension(output_path);

    if (ext == ".npy") {
        std::string shape = "(" + std::to_string(n_rows) + ", "
                          + std::to_string(n_cols) + ")";
        WriteNpy(output_path, data, n_rows * n_cols, shape);
        WriteSidecar(SidecarPath(output_path), meta, n_rows, n_cols);
    } else if (ext == ".csv") {
        WriteCDistCsv(output_path, data, n_rows, n_cols, meta);
    } else if (ext == ".bin") {
        WriteBin(output_path, data, n_rows * n_cols);
        WriteSidecar(SidecarPath(output_path), meta, n_rows, n_cols);
    } else {
        throw std::runtime_error(
            "Unknown output format: " + ext +
            " (use .npy, .csv, or .bin)");
    }
}

}  // namespace OEPDist
