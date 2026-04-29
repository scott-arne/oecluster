/**
 * @file StorageBackend.cpp
 * @brief Implementation of storage backend classes.
 */

#include "oecluster/StorageBackend.h"
#include <cassert>
#include <algorithm>
#include <mutex>
#include <stdexcept>
#include <cstring>
#ifdef _WIN32
// Win32 has no POSIX mmap; use CreateFileMapping + MapViewOfFile instead.
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#else
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <cerrno>
#endif

namespace OECluster {

// DenseStorage implementation

DenseStorage::DenseStorage(size_t n)
    : n_(n), data_(n * (n - 1) / 2, 0.0) {
}

void DenseStorage::Set(size_t i, size_t j, double value) {
    // Ensure i < j for upper triangle storage
    if (i > j) {
        std::swap(i, j);
    }
    assert(i < j && "Set requires i < j after swap");

    size_t index = CondensedIndex(n_, i, j);
    data_[index] = value;
}

double DenseStorage::Get(size_t i, size_t j) const {
    // Diagonal is always zero
    if (i == j) {
        return 0.0;
    }

    // Ensure i < j for upper triangle storage
    if (i > j) {
        std::swap(i, j);
    }

    size_t index = CondensedIndex(n_, i, j);
    return data_[index];
}

size_t DenseStorage::NumItems() const {
    return n_;
}

size_t DenseStorage::NumPairs() const {
    return data_.size();
}

double* DenseStorage::Data() {
    return data_.data();
}

const double* DenseStorage::Data() const {
    return data_.data();
}

size_t DenseStorage::CondensedIndex(size_t n, size_t i, size_t j) {
    // Formula: n * i - i * (i + 1) / 2 + j - i - 1
    // This matches scipy.spatial.distance.squareform condensed indexing
    assert(i < j && "CondensedIndex requires i < j");
    return n * i - i * (i + 1) / 2 + j - i - 1;
}

// MMapStorage implementation

#ifdef _WIN32

namespace {

// Format a human-readable message from GetLastError() for diagnostics.
std::string win32_error_message(DWORD code) {
    LPSTR buf = nullptr;
    DWORD len = FormatMessageA(
        FORMAT_MESSAGE_ALLOCATE_BUFFER | FORMAT_MESSAGE_FROM_SYSTEM |
            FORMAT_MESSAGE_IGNORE_INSERTS,
        nullptr, code, MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT),
        reinterpret_cast<LPSTR>(&buf), 0, nullptr);
    std::string msg = (len && buf) ? std::string(buf, len) : "Unknown error";
    if (buf) {
        LocalFree(buf);
    }
    while (!msg.empty() && (msg.back() == '\n' || msg.back() == '\r')) {
        msg.pop_back();
    }
    return msg + " (code " + std::to_string(code) + ")";
}

}  // namespace

MMapStorage::MMapStorage(const std::string& path, size_t n)
    : n_(n),
      num_pairs_(n * (n - 1) / 2),
      file_size_(num_pairs_ * sizeof(double)),
      path_(path),
      file_handle_(INVALID_HANDLE_VALUE),
      mapping_handle_(nullptr),
      data_(nullptr) {

    // CREATE_ALWAYS would truncate; OPEN_ALWAYS keeps existing contents,
    // matching the POSIX O_RDWR|O_CREAT behavior used on other platforms.
    HANDLE fh = CreateFileA(path.c_str(), GENERIC_READ | GENERIC_WRITE,
                            FILE_SHARE_READ, nullptr, OPEN_ALWAYS,
                            FILE_ATTRIBUTE_NORMAL, nullptr);
    if (fh == INVALID_HANDLE_VALUE) {
        throw std::runtime_error("MMapStorage: Failed to open file " + path +
                                 ": " + win32_error_message(GetLastError()));
    }
    file_handle_ = fh;

    LARGE_INTEGER existing_size;
    if (!GetFileSizeEx(fh, &existing_size)) {
        DWORD err = GetLastError();
        CloseHandle(fh);
        file_handle_ = INVALID_HANDLE_VALUE;
        throw std::runtime_error("MMapStorage: Failed to stat file " + path +
                                 ": " + win32_error_message(err));
    }

    const bool needs_extend =
        static_cast<size_t>(existing_size.QuadPart) < file_size_;

    if (static_cast<size_t>(existing_size.QuadPart) != file_size_) {
        LARGE_INTEGER new_size;
        new_size.QuadPart = static_cast<LONGLONG>(file_size_);
        if (!SetFilePointerEx(fh, new_size, nullptr, FILE_BEGIN) ||
            !SetEndOfFile(fh)) {
            DWORD err = GetLastError();
            CloseHandle(fh);
            file_handle_ = INVALID_HANDLE_VALUE;
            throw std::runtime_error("MMapStorage: Failed to resize file " +
                                     path + ": " + win32_error_message(err));
        }
    }

    const DWORD size_high =
        static_cast<DWORD>((static_cast<uint64_t>(file_size_) >> 32) & 0xFFFFFFFFu);
    const DWORD size_low =
        static_cast<DWORD>(static_cast<uint64_t>(file_size_) & 0xFFFFFFFFu);
    HANDLE mh = CreateFileMappingA(fh, nullptr, PAGE_READWRITE,
                                   size_high, size_low, nullptr);
    if (mh == nullptr) {
        DWORD err = GetLastError();
        CloseHandle(fh);
        file_handle_ = INVALID_HANDLE_VALUE;
        throw std::runtime_error("MMapStorage: Failed to create file mapping " +
                                 path + ": " + win32_error_message(err));
    }
    mapping_handle_ = mh;

    void* view = MapViewOfFile(mh, FILE_MAP_READ | FILE_MAP_WRITE, 0, 0,
                               file_size_);
    if (view == nullptr) {
        DWORD err = GetLastError();
        CloseHandle(mh);
        CloseHandle(fh);
        mapping_handle_ = nullptr;
        file_handle_ = INVALID_HANDLE_VALUE;
        throw std::runtime_error("MMapStorage: Failed to map view of file " +
                                 path + ": " + win32_error_message(err));
    }
    data_ = static_cast<double*>(view);

    // SetEndOfFile leaves new bytes undefined on some filesystems; match the
    // POSIX branch by zeroing any extended region.
    if (needs_extend) {
        std::memset(data_, 0, file_size_);
    }
}

MMapStorage::~MMapStorage() {
    if (data_ != nullptr) {
        UnmapViewOfFile(data_);
    }
    if (mapping_handle_ != nullptr) {
        CloseHandle(mapping_handle_);
    }
    if (file_handle_ != INVALID_HANDLE_VALUE && file_handle_ != nullptr) {
        CloseHandle(file_handle_);
    }
}

#else  // !_WIN32

MMapStorage::MMapStorage(const std::string& path, size_t n)
    : n_(n),
      num_pairs_(n * (n - 1) / 2),
      file_size_(num_pairs_ * sizeof(double)),
      path_(path),
      fd_(-1),
      data_(nullptr) {

    // Open or create the file
    fd_ = open(path.c_str(), O_RDWR | O_CREAT, 0644);
    if (fd_ == -1) {
        throw std::runtime_error("MMapStorage: Failed to open file " + path +
                                 ": " + std::strerror(errno));
    }

    // Check existing file size
    struct stat st;
    if (fstat(fd_, &st) == -1) {
        close(fd_);
        throw std::runtime_error("MMapStorage: Failed to stat file " + path +
                                 ": " + std::strerror(errno));
    }

    // Resize file if needed
    if (static_cast<size_t>(st.st_size) != file_size_) {
        if (ftruncate(fd_, file_size_) == -1) {
            close(fd_);
            throw std::runtime_error("MMapStorage: Failed to resize file " + path +
                                     ": " + std::strerror(errno));
        }
    }

    // Memory map the file
    data_ = static_cast<double*>(mmap(nullptr, file_size_, PROT_READ | PROT_WRITE,
                                      MAP_SHARED, fd_, 0));
    if (data_ == MAP_FAILED) {
        close(fd_);
        throw std::runtime_error("MMapStorage: Failed to mmap file " + path +
                                 ": " + std::strerror(errno));
    }

    // If we created a new file or extended it, zero-fill the new regions
    if (static_cast<size_t>(st.st_size) < file_size_) {
        std::memset(data_, 0, file_size_);
    }
}

MMapStorage::~MMapStorage() {
    if (data_ != nullptr && data_ != MAP_FAILED) {
        munmap(data_, file_size_);
    }
    if (fd_ != -1) {
        close(fd_);
    }
}

#endif  // _WIN32

void MMapStorage::Set(size_t i, size_t j, double value) {
    // Ensure i < j for upper triangle storage
    if (i > j) {
        std::swap(i, j);
    }
    assert(i < j && "Set requires i < j after swap");

    size_t index = CondensedIndex(n_, i, j);
    data_[index] = value;
}

double MMapStorage::Get(size_t i, size_t j) const {
    // Diagonal is always zero
    if (i == j) {
        return 0.0;
    }

    // Ensure i < j for upper triangle storage
    if (i > j) {
        std::swap(i, j);
    }

    size_t index = CondensedIndex(n_, i, j);
    return data_[index];
}

size_t MMapStorage::NumItems() const {
    return n_;
}

size_t MMapStorage::NumPairs() const {
    return num_pairs_;
}

double* MMapStorage::Data() {
    return data_;
}

const double* MMapStorage::Data() const {
    return data_;
}

size_t MMapStorage::CondensedIndex(size_t n, size_t i, size_t j) {
    // Formula: n * i - i * (i + 1) / 2 + j - i - 1
    // This matches scipy.spatial.distance.squareform condensed indexing
    assert(i < j && "CondensedIndex requires i < j");
    return n * i - i * (i + 1) / 2 + j - i - 1;
}

// SparseStorage implementation

SparseStorage::SparseStorage(size_t n, double cutoff)
    : n_(n), cutoff_(cutoff) {
}

void SparseStorage::Set(size_t i, size_t j, double value) {
    // Only store values at or below cutoff
    if (value > cutoff_) {
        return;
    }

    // Ensure i < j for consistent storage
    if (i > j) {
        std::swap(i, j);
    }
    assert(i < j && "Set requires i < j after swap");

    // Get the calling thread's buffer
    std::thread::id tid = std::this_thread::get_id();

    // Try with shared lock first (read access) to check if buffer exists
    bool buffer_exists = false;
    {
        std::shared_lock<std::shared_mutex> lock(buffers_mutex_);
        buffer_exists = thread_buffers_.find(tid) != thread_buffers_.end();
    }

    // Get unique lock for write access
    std::unique_lock<std::shared_mutex> lock(buffers_mutex_);
    thread_buffers_[tid].emplace_back(i, j, value);
}

double SparseStorage::Get(size_t i, size_t j) const {
    // Diagonal is always zero
    if (i == j) {
        return 0.0;
    }

    // Ensure i < j for consistent lookup
    if (i > j) {
        std::swap(i, j);
    }

    size_t index = CondensedIndex(n_, i, j);
    auto it = lookup_.find(index);
    if (it != lookup_.end()) {
        return it->second;
    }
    return 0.0;
}

size_t SparseStorage::NumItems() const {
    return n_;
}

size_t SparseStorage::NumPairs() const {
    return n_ * (n_ - 1) / 2;
}

double* SparseStorage::Data() {
    return nullptr;
}

const double* SparseStorage::Data() const {
    return nullptr;
}

void SparseStorage::Finalize() {
    // Merge all per-thread buffers into merged_entries_
    std::unique_lock<std::shared_mutex> lock(buffers_mutex_);

    for (auto& [tid, buffer] : thread_buffers_) {
        merged_entries_.insert(merged_entries_.end(), buffer.begin(), buffer.end());
    }

    // Build lookup map for fast Get() access
    for (const auto& [i, j, value] : merged_entries_) {
        size_t index = CondensedIndex(n_, i, j);
        lookup_[index] = value;
    }

    // Clear thread buffers to free memory
    thread_buffers_.clear();
}

const std::vector<std::tuple<size_t, size_t, double>>& SparseStorage::Entries() const {
    return merged_entries_;
}

size_t SparseStorage::CondensedIndex(size_t n, size_t i, size_t j) {
    // Formula: n * i - i * (i + 1) / 2 + j - i - 1
    // This matches scipy.spatial.distance.squareform condensed indexing
    assert(i < j && "CondensedIndex requires i < j");
    return n * i - i * (i + 1) / 2 + j - i - 1;
}

}  // namespace OECluster
