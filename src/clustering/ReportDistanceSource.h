/**
 * @file ReportDistanceSource.h
 * @brief Distance sources for the cluster_report engine.
 *
 * The engine states each pass as an ordered stream of rows -- one anchor
 * against a run of targets -- and reads the values back in that order. The
 * matrix source answers from storage as the engine always has. The comparison
 * source fills bounded blocks of rows in parallel through a clone pool, then
 * hands the values back serially, so the reducer's arithmetic sees the same
 * values in the same order on both paths.
 */

#ifndef OECLUSTER_SRC_CLUSTERING_REPORTDISTANCESOURCE_H
#define OECLUSTER_SRC_CLUSTERING_REPORTDISTANCESOURCE_H

#include <algorithm>
#include <cstddef>
#include <exception>
#include <functional>
#include <limits>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <thread>
#include <utility>
#include <vector>

#include "ChunkedComparisons.h"
#include "InternalIndices.h"
#include "oecluster/PairwiseComparison.h"
#include "oecluster/StorageBackend.h"
#include "oecluster/clustering/ClusterTypes.h"

namespace OECluster::detail {

/// The fill block cap, in distances: 8 MB of doubles. A row longer than this
/// is a block of its own.
constexpr size_t FILL_BLOCK_DISTANCES = size_t{1} << 20;

struct ReportRow {
    size_t anchor;
    const size_t* targets;
    size_t count;
};

/// Pairs (members[k][i], members[k][j > i]) for k in [first, last).
class IntraRows {
public:
    IntraRows(const Clusters& members, size_t first, size_t last)
        : members_(&members), k_(first), last_(last) {}

    bool Next(ReportRow& row) {
        while (k_ < last_) {
            const Cluster& cluster = (*members_)[k_];
            if (i_ + 1 < cluster.size()) {
                row = {cluster[i_], cluster.data() + i_ + 1, cluster.size() - i_ - 1};
                ++i_;
                return true;
            }
            ++k_;
            i_ = 0;
        }
        return false;
    }

private:
    const Clusters* members_;
    size_t k_;
    size_t last_;
    size_t i_ = 0;
};

/// Pairs (i in members[a], j in members[b]) for a < b.
class CrossRows {
public:
    explicit CrossRows(const Clusters& members) : members_(&members) {}

    bool Next(ReportRow& row) {
        const size_t k = members_->size();
        while (a_ + 1 < k) {
            if (b_ >= k) {
                ++a_;
                b_ = a_ + 1;
                i_ = 0;
                continue;
            }
            const Cluster& left = (*members_)[a_];
            if (i_ < left.size()) {
                const Cluster& right = (*members_)[b_];
                row = {left[i_], right.data(), right.size()};
                ++i_;
                return true;
            }
            ++b_;
            i_ = 0;
        }
        return false;
    }

private:
    const Clusters* members_;
    size_t a_ = 0;
    size_t b_ = 1;
    size_t i_ = 0;
};

/// Per cluster: first[k] against every member, then second[k] against every
/// member when second is given.
class MemberRows {
public:
    MemberRows(const Clusters& members, const std::vector<size_t>& first,
               const std::vector<size_t>* second)
        : members_(&members), first_(&first), second_(second) {}

    bool Next(ReportRow& row) {
        if (k_ >= members_->size()) {
            return false;
        }
        const Cluster& cluster = (*members_)[k_];
        if (!on_second_) {
            row = {(*first_)[k_], cluster.data(), cluster.size()};
            if (second_ != nullptr) {
                on_second_ = true;
            } else {
                ++k_;
            }
            return true;
        }
        row = {(*second_)[k_], cluster.data(), cluster.size()};
        on_second_ = false;
        ++k_;
        return true;
    }

private:
    const Clusters* members_;
    const std::vector<size_t>* first_;
    const std::vector<size_t>* second_;
    size_t k_ = 0;
    bool on_second_ = false;
};

/// Each point against every point in the list, itself included; the engine
/// reads and discards the self entry, which both sources answer as 0.
class PeerRows {
public:
    explicit PeerRows(const std::vector<size_t>& points) : points_(&points) {}

    bool Next(ReportRow& row) {
        if (a_ >= points_->size()) {
            return false;
        }
        row = {(*points_)[a_], points_->data(), points_->size()};
        ++a_;
        return true;
    }

private:
    const std::vector<size_t>* points_;
    size_t a_ = 0;
};

/// Each anchor against one fixed target.
class FixedTargetRows {
public:
    FixedTargetRows(const std::vector<size_t>& anchors, const size_t* target)
        : anchors_(&anchors), target_(target) {}

    bool Next(ReportRow& row) {
        if (a_ >= anchors_->size()) {
            return false;
        }
        row = {(*anchors_)[a_], target_, 1};
        ++a_;
        return true;
    }

private:
    const std::vector<size_t>* anchors_;
    const size_t* target_;
    size_t a_ = 0;
};

/// Every point in [0, num_points) against the target list.
class AllPointsRows {
public:
    AllPointsRows(size_t num_points, const std::vector<size_t>& targets)
        : num_points_(targets.empty() ? 0 : num_points), targets_(&targets) {}

    bool Next(ReportRow& row) {
        if (p_ >= num_points_) {
            return false;
        }
        row = {p_, targets_->data(), targets_->size()};
        ++p_;
        return true;
    }

private:
    size_t num_points_;
    const std::vector<size_t>* targets_;
    size_t p_ = 0;
};

/// Reads storage exactly as the engine always has.
class MatrixSource {
public:
    class Reader {
    public:
        explicit Reader(const StorageBackend& storage) : storage_(&storage) {}
        double Next(size_t anchor, size_t target) {
            return checked_distance(*storage_, anchor, target);
        }
        void Finish() {}

    private:
        const StorageBackend* storage_;
    };

    explicit MatrixSource(const StorageBackend& storage) : storage_(storage) {}

    size_t NumItems() const { return storage_.NumSamples(); }

    template <class Plan>
    Reader Open(Plan) {
        return Reader(storage_);
    }

private:
    const StorageBackend& storage_;
};

template <class Plan>
class ComparisonReader;

class ComparisonSource {
public:
    ComparisonSource(PairwiseComparison& comparison, size_t num_threads,
                     size_t chunk_size, size_t fill_block = FILL_BLOCK_DISTANCES,
                     std::function<void(size_t, size_t)> on_block = {})
        : comparison_(comparison),
          num_threads_(num_threads),
          chunk_size_(chunk_size),
          fill_block_(fill_block),
          on_block_(std::move(on_block)) {}

    size_t NumItems() const { return comparison_.Size(); }

    template <class Plan>
    ComparisonReader<Plan> Open(Plan plan);

    size_t ChunkSize() const { return chunk_size_; }
    size_t FillBlock() const { return fill_block_; }
    const std::function<void(size_t, size_t)>& OnBlock() const { return on_block_; }

    // One pool for the whole report: clones are made on first use and reused
    // by every later pass. A unit is one ParallelFor chunk. Zero is resolved
    // here rather than in the pool so the item cap applies to it too: units
    // can outnumber items, and an unresolved zero would let the hardware
    // concurrency, not n, bound the workers and clones.
    ChunkedComparisons& Pool() {
        if (!pool_) {
            const size_t threads = num_threads_ > 0
                ? num_threads_
                : std::max<size_t>(1, std::thread::hardware_concurrency());
            pool_ = std::make_unique<ChunkedComparisons>(
                comparison_, comparison_.Size(), threads, 1);
        }
        return *pool_;
    }

    /// Fills one row. Returns how many entries were written before Compare
    /// threw; on a throw the exception is stored and the rest is skipped.
    static size_t FillRow(PairwiseComparison& comparison, const ReportRow& row,
                          double* out, std::exception_ptr& error) {
        for (size_t k = 0; k < row.count; ++k) {
            const size_t t = row.targets[k];
            try {
                out[k] = row.anchor == t
                    ? 0.0
                    : comparison.Compare(std::min(row.anchor, t), std::max(row.anchor, t));
            } catch (...) {
                error = std::current_exception();
                return k;
            }
        }
        return row.count;
    }

private:
    PairwiseComparison& comparison_;
    size_t num_threads_;
    size_t chunk_size_;
    size_t fill_block_;
    std::function<void(size_t, size_t)> on_block_;
    std::unique_ptr<ChunkedComparisons> pool_;
};

template <class Plan>
class ComparisonReader {
public:
    ComparisonReader(ComparisonSource& source, Plan plan)
        : source_(&source), plan_(std::move(plan)) {}

    double Next(size_t anchor, size_t target) {
        if (row_index_ == rows_.size()) {
            Fill();
            if (rows_.empty()) {
                throw std::logic_error("cluster_report: read past the end of a row plan");
            }
        }
        const ReportRow& row = rows_[row_index_];
        if (row.anchor != anchor || row.targets[column_] != target) {
            throw std::logic_error("cluster_report: read out of plan order");
        }
        const size_t position = offsets_[row_index_] + column_;
        if (error_ && position == error_position_) {
            std::rethrow_exception(error_);
        }
        const double value = block_[position];
        if (++column_ == row.count) {
            column_ = 0;
            ++row_index_;
        }
        return finite_report_distance(value, anchor, target);
    }

    void Finish() {
        if (row_index_ != rows_.size() || has_pending_) {
            throw std::logic_error("cluster_report: a row plan was not read to the end");
        }
        ReportRow row{};
        if (plan_.Next(row)) {
            throw std::logic_error("cluster_report: a row plan was not read to the end");
        }
    }

private:
    // Gathers rows in plan order up to the cap, splits them into units of
    // about chunk_size distances, and fills the units in parallel. An error
    // is kept only at the earliest position, so the serial read meets the
    // same failure the matrix path would, whatever order units finished in.
    void Fill() {
        rows_.clear();
        offsets_.clear();
        row_index_ = 0;
        column_ = 0;
        error_ = nullptr;
        size_t total = 0;
        const size_t cap = source_->FillBlock();
        while (true) {
            ReportRow row{};
            if (has_pending_) {
                row = pending_;
                has_pending_ = false;
            } else if (!plan_.Next(row)) {
                break;
            }
            if (row.count == 0) {
                continue;
            }
            if (rows_.empty() || (total < cap && row.count <= cap - total)) {
                rows_.push_back(row);
                offsets_.push_back(total);
                total += row.count;
            } else {
                pending_ = row;
                has_pending_ = true;
                break;
            }
        }
        if (rows_.empty()) {
            return;
        }
        block_.resize(total);

        const size_t chunk = source_->ChunkSize();
        std::vector<size_t> unit_starts;
        size_t in_unit = 0;
        for (size_t r = 0; r < rows_.size(); ++r) {
            const size_t count = rows_[r].count;
            if (unit_starts.empty() || in_unit >= chunk || count > chunk - in_unit) {
                unit_starts.push_back(r);
                in_unit = 0;
            }
            in_unit += count;
        }
        unit_starts.push_back(rows_.size());

        std::mutex error_mutex;
        size_t earliest = std::numeric_limits<size_t>::max();
        std::exception_ptr earliest_error;
        source_->Pool().Run(unit_starts.size() - 1, [&](PairwiseComparison& clone,
                                                        size_t begin, size_t end) {
            for (size_t u = begin; u < end; ++u) {
                for (size_t r = unit_starts[u]; r < unit_starts[u + 1]; ++r) {
                    std::exception_ptr error;
                    const size_t filled = ComparisonSource::FillRow(
                        clone, rows_[r], block_.data() + offsets_[r], error);
                    if (error) {
                        const std::lock_guard<std::mutex> lock(error_mutex);
                        const size_t position = offsets_[r] + filled;
                        if (position < earliest) {
                            earliest = position;
                            earliest_error = error;
                        }
                        break;
                    }
                }
            }
        });
        error_ = earliest_error;
        error_position_ = earliest;
        if (source_->OnBlock()) {
            source_->OnBlock()(total, rows_.size());
        }
    }

    ComparisonSource* source_;
    Plan plan_;
    std::vector<ReportRow> rows_;
    std::vector<size_t> offsets_;
    std::vector<double> block_;
    size_t row_index_ = 0;
    size_t column_ = 0;
    ReportRow pending_{};
    bool has_pending_ = false;
    std::exception_ptr error_;
    size_t error_position_ = 0;
};

template <class Plan>
ComparisonReader<Plan> ComparisonSource::Open(Plan plan) {
    return ComparisonReader<Plan>(*this, std::move(plan));
}

}  // namespace OECluster::detail

#endif  // OECLUSTER_SRC_CLUSTERING_REPORTDISTANCESOURCE_H
