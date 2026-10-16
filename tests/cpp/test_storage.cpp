#include <gtest/gtest.h>
#include <cstddef>
#include <filesystem>
#include <limits>
#include <stdexcept>
#include <string>
#include <system_error>
#include "oecluster/StorageBackend.h"
#include "../../src/StorageOverflow.h"

using namespace OECluster;

TEST(DenseStorageTest, ConstructorSetsSize) {
    DenseStorage storage(5);
    EXPECT_EQ(storage.NumSamples(), 5);
    EXPECT_EQ(storage.NumPairs(), 10);  // 5*4/2
}

TEST(DenseStorageTest, SetAndGet) {
    DenseStorage storage(4);
    storage.Set(0, 1, 0.5);
    storage.Set(0, 2, 0.8);
    storage.Set(1, 3, 0.2);
    EXPECT_DOUBLE_EQ(storage.Get(0, 1), 0.5);
    EXPECT_DOUBLE_EQ(storage.Get(0, 2), 0.8);
    EXPECT_DOUBLE_EQ(storage.Get(1, 3), 0.2);
}

TEST(DenseStorageTest, GetSymmetric) {
    DenseStorage storage(4);
    storage.Set(0, 1, 0.5);
    EXPECT_DOUBLE_EQ(storage.Get(1, 0), 0.5);
}

TEST(DenseStorageTest, GetDiagonalIsZero) {
    DenseStorage storage(4);
    EXPECT_DOUBLE_EQ(storage.Get(0, 0), 0.0);
    EXPECT_DOUBLE_EQ(storage.Get(2, 2), 0.0);
}

TEST(DenseStorageTest, DataPointerNotNull) {
    DenseStorage storage(4);
    EXPECT_NE(storage.Data(), nullptr);
}

TEST(DenseStorageTest, DataPointerMatchesCondensedIndex) {
    DenseStorage storage(4);  // N=4, pairs=6
    storage.Set(0, 1, 1.0);
    storage.Set(0, 2, 2.0);
    storage.Set(0, 3, 3.0);
    storage.Set(1, 2, 4.0);
    storage.Set(1, 3, 5.0);
    storage.Set(2, 3, 6.0);
    const double* data = storage.Data();
    // Condensed order: (0,1), (0,2), (0,3), (1,2), (1,3), (2,3)
    EXPECT_DOUBLE_EQ(data[0], 1.0);
    EXPECT_DOUBLE_EQ(data[1], 2.0);
    EXPECT_DOUBLE_EQ(data[2], 3.0);
    EXPECT_DOUBLE_EQ(data[3], 4.0);
    EXPECT_DOUBLE_EQ(data[4], 5.0);
    EXPECT_DOUBLE_EQ(data[5], 6.0);
}

TEST(DenseStorageTest, InitializedToZero) {
    DenseStorage storage(4);
    for (size_t i = 0; i < storage.NumPairs(); ++i) {
        EXPECT_DOUBLE_EQ(storage.Data()[i], 0.0);
    }
}

TEST(DenseStorageTest, GetOutOfRangeThrows) {
    DenseStorage storage(2);
    EXPECT_THROW(storage.Get(0, 2), std::out_of_range);
    EXPECT_THROW(storage.Get(1, 1399), std::out_of_range);
    EXPECT_THROW(storage.Get(500, 1399), std::out_of_range);
}

TEST(DenseStorageTest, GetOutOfRangeDiagonalThrows) {
    // The diagonal shortcut used to answer 0.0 for an index naming no item.
    DenseStorage storage(2);
    EXPECT_THROW(storage.Get(1000, 1000), std::out_of_range);
}

TEST(DenseStorageTest, GetOutOfRangeMessageNamesIndexAndCount) {
    DenseStorage storage(2);
    try {
        storage.Get(1, 1399);
        FAIL() << "Expected std::out_of_range";
    } catch (const std::out_of_range& e) {
        EXPECT_STREQ(e.what(),
                     "DenseStorage index 1399 is outside the storage range of 2 samples");
    }
}

TEST(DenseStorageTest, SetOutOfRangeThrowsRatherThanAliasingAnotherPair) {
    // CondensedIndex(5, 2, 5) == 9 == CondensedIndex(5, 3, 4): just past the
    // end, the computed index lands on a real pair, so the write used to
    // succeed and replace a distance the caller never named.
    DenseStorage storage(5);
    storage.Set(3, 4, 0.25);
    EXPECT_THROW(storage.Set(2, 5, 0.99), std::out_of_range);
    EXPECT_DOUBLE_EQ(storage.Get(3, 4), 0.25);
}

TEST(DenseStorageTest, SetFarOutOfRangeThrowsRatherThanWritingPastTheBuffer) {
    // Far past the end there is no pair to collide with: the index leaves the
    // allocation and the write corrupts whatever follows it, undiagnosed.
    DenseStorage storage(5);
    EXPECT_THROW(storage.Set(0, 5000, 1.0), std::out_of_range);
}

TEST(DenseStorageTest, SetOutOfRangeMessageNamesIndexAndCount) {
    DenseStorage storage(2);
    try {
        storage.Set(1, 1399, 0.5);
        FAIL() << "Expected std::out_of_range";
    } catch (const std::out_of_range& e) {
        EXPECT_STREQ(e.what(),
                     "DenseStorage index 1399 is outside the storage range of 2 samples");
    }
}

TEST(DenseStorageTest, SetDiagonalThrowsRatherThanAliasingAnotherPair) {
    // CondensedIndex(5, 2, 2) == 6 == CondensedIndex(5, 1, 4). Get answers the
    // diagonal 0.0 without consulting storage, so the diagonal owns no slot to
    // write; the condensed formula hands the write a neighbour's slot instead.
    DenseStorage storage(5);
    storage.Set(1, 4, 0.75);
    EXPECT_THROW(storage.Set(2, 2, 0.99), std::invalid_argument);
    EXPECT_DOUBLE_EQ(storage.Get(1, 4), 0.75);
}

TEST(DenseStorageTest, SetDiagonalMessageNamesTheIndex) {
    DenseStorage storage(5);
    try {
        storage.Set(3, 3, 0.5);
        FAIL() << "Expected std::invalid_argument";
    } catch (const std::invalid_argument& e) {
        EXPECT_STREQ(e.what(),
                     "DenseStorage cannot store the diagonal pair (3, 3); "
                     "only distances between distinct items are stored");
    }
}

TEST(MMapStorageTest, CreateAndWrite) {
    auto path = std::filesystem::temp_directory_path() / "test_mmap.bin";
    {
        MMapStorage storage(path.string(), 4);
        storage.Set(0, 1, 0.5);
        storage.Set(2, 3, 0.8);
        EXPECT_DOUBLE_EQ(storage.Get(0, 1), 0.5);
        EXPECT_DOUBLE_EQ(storage.Get(2, 3), 0.8);
        EXPECT_NE(storage.Data(), nullptr);
    }
    EXPECT_TRUE(std::filesystem::exists(path));
    std::filesystem::remove(path);
}

TEST(MMapStorageTest, PersistsAfterDestruction) {
    auto path = std::filesystem::temp_directory_path() / "test_mmap_persist.bin";
    {
        MMapStorage storage(path.string(), 4);
        storage.Set(0, 1, 0.5);
        storage.Set(1, 2, 0.9);
    }
    {
        MMapStorage storage(path.string(), 4);
        EXPECT_DOUBLE_EQ(storage.Get(0, 1), 0.5);
        EXPECT_DOUBLE_EQ(storage.Get(1, 2), 0.9);
    }
    std::filesystem::remove(path);
}

TEST(MMapStorageTest, GetOutOfRangeThrows) {
    auto path = std::filesystem::temp_directory_path() / "test_mmap_range.bin";
    {
        MMapStorage storage(path.string(), 2);
        EXPECT_THROW(storage.Get(0, 2), std::out_of_range);
        EXPECT_THROW(storage.Get(1000, 1000), std::out_of_range);
    }
    std::filesystem::remove(path);
}

TEST(MMapStorageTest, SetOutOfRangeThrows) {
    // The mapping is sized for exactly NumPairs doubles, so an unchecked write
    // past the end is a write past the mapping, not merely past a vector.
    auto path = std::filesystem::temp_directory_path() / "test_mmap_set_range.bin";
    {
        MMapStorage storage(path.string(), 2);
        EXPECT_THROW(storage.Set(0, 2, 1.0), std::out_of_range);
        // Out of range on the diagonal is reported as a range error: the
        // indices name no item at all, which is the more basic complaint.
        EXPECT_THROW(storage.Set(1000, 1000, 1.0), std::out_of_range);
        EXPECT_THROW(storage.Set(1, 1, 1.0), std::invalid_argument);
    }
    std::filesystem::remove(path);
}

TEST(SparseStorageTest, StoresOnlyBelowCutoff) {
    SparseStorage storage(4, 0.5);  // cutoff = 0.5
    storage.Set(0, 1, 0.3);  // below cutoff, stored
    storage.Set(0, 2, 0.8);  // above cutoff, ignored
    storage.Set(1, 2, 0.5);  // at cutoff, stored
    storage.Finalize();
    auto& entries = storage.Entries();
    EXPECT_EQ(entries.size(), 2);
}

TEST(SparseStorageTest, DataReturnsNullptr) {
    SparseStorage storage(4, 0.5);
    EXPECT_EQ(storage.Data(), nullptr);
}

TEST(SparseStorageTest, GetReturnsStoredValues) {
    SparseStorage storage(4, 0.5);
    storage.Set(0, 1, 0.3);
    storage.Finalize();
    EXPECT_DOUBLE_EQ(storage.Get(0, 1), 0.3);
}

TEST(SparseStorageTest, GetReturnsZeroForUnstored) {
    SparseStorage storage(4, 0.5);
    storage.Finalize();
    EXPECT_DOUBLE_EQ(storage.Get(0, 1), 0.0);
}

TEST(SparseStorageTest, GetOutOfRangeThrowsRatherThanMissing) {
    // "Not stored" answers 0.0; "not an item" must not borrow that answer.
    SparseStorage storage(4, 0.5);
    storage.Finalize();
    EXPECT_THROW(storage.Get(0, 4), std::out_of_range);
    EXPECT_THROW(storage.Get(1000, 1000), std::out_of_range);
}

TEST(SparseStorageTest, SetOutOfRangeThrowsWhicheverSideOfTheCutoff) {
    // The cutoff shortcut returns before anything is stored. Checking the
    // indices behind it would make the diagnosis depend on the value the bad
    // call happened to carry: below the cutoff a refusal, above it silence.
    SparseStorage storage(4, 0.5);
    EXPECT_THROW(storage.Set(0, 4, 0.3), std::out_of_range);
    EXPECT_THROW(storage.Set(0, 4, 0.9), std::out_of_range);
    EXPECT_THROW(storage.Set(2, 2, 0.3), std::invalid_argument);
    EXPECT_THROW(storage.Set(2, 2, 0.9), std::invalid_argument);
}

TEST(SparseStorageTest, SetOutOfRangeStoresNothing) {
    // A refused Set must not leave a half-written entry behind for Finalize to
    // pick up -- an entry whose condensed index collides with a real pair.
    SparseStorage storage(4, 0.5);
    EXPECT_THROW(storage.Set(0, 4, 0.3), std::out_of_range);
    storage.Finalize();
    EXPECT_EQ(storage.Entries().size(), 0u);
}

// Sample counts whose derived buffer length wraps. DenseStorage and
// MMapStorage both size their storage from n * (n - 1) / 2 doubles, and a
// large enough n used to wrap one of those two multiplications and report a
// NumSamples() the buffer cannot back: DenseStorage(SIZE_MAX) constructed with
// NumPairs() == 1.
//
// The boundary values are tested against checked_pair_count rather than
// against a constructor, because the only run that matters is the one where
// the guard has regressed, and then most of these counts are not an assertion
// failure but a dead machine: 6074001001 wraps to 3327948884 pairs, which an
// unguarded DenseStorage zero-fills as 26.6 GB of memory and an unguarded
// MMapStorage writes to disk; 2147483649 has a byte count that wraps to
// exactly 8 GiB, which an unguarded MMapStorage writes. A test that is
// OOM-killed has reported nothing, and no destructor of ours runs to clean up
// after it. Constructing is reserved for the few counts that stay cheap even
// unguarded; those cases prove the wiring, and the helper cases prove the
// arithmetic.

namespace {

/// Runs `call` and returns the message of the std::overflow_error it threw, or
/// a sentinel no expectation below can match. gtest has no message matcher
/// without gmock, and EXPECT_THROW discards the exception it catches.
template <typename Call>
std::string overflow_message(Call call) {
    try {
        call();
    } catch (const std::overflow_error& error) {
        return error.what();
    } catch (...) {
        return "<threw some other exception type>";
    }
    return "<threw nothing>";
}

/// Path for an MMapStorage refusal probe, removed at both ends of the test.
/// Removing it first makes the existence check after the probe describe only
/// the call under test; removing it again on the way out clears any residual
/// artifact a regressed guard left, which for the counts constructed below is
/// an 8-byte file or an empty one.
class ProbePath {
public:
    explicit ProbePath(const char* name)
        : path_((std::filesystem::temp_directory_path() / name).string()) {
        std::filesystem::remove(path_);
    }

    ~ProbePath() {
        // A destructor must not throw, and a probe that left nothing behind is
        // the expected case rather than an error.
        std::error_code ignored;
        std::filesystem::remove(path_, ignored);
    }

    ProbePath(const ProbePath&) = delete;
    ProbePath& operator=(const ProbePath&) = delete;

    const std::string& Path() const { return path_; }

private:
    std::string path_;
};

/// Largest n accepted: 2305843008139952128 pairs, 18446744065119617024 bytes,
/// both inside size_t.
constexpr size_t LARGEST_ACCEPTED_N = 2147483648ULL;

/// Largest n whose pair count fits in size_t (18446744070963499500). Its byte
/// count does not, so it is refused -- but by the byte check, which is what
/// distinguishes an exact predicate from one that divides SIZE_MAX by n and
/// doubles, losing the remainder.
constexpr size_t LARGEST_REPRESENTABLE_PAIRS_N = 6074001000ULL;

/// Smallest n whose pair count exceeds size_t, one past the count above.
constexpr size_t SMALLEST_WRAPPING_PAIRS_N = 6074001001ULL;

/// Smallest n whose pair count fits in size_t while its byte count does not:
/// 2305843010287435776 pairs, whose byte count wraps to 8 GiB.
constexpr size_t SMALLEST_WRAPPING_BYTES_N = 2147483649ULL;

/// Smallest n whose byte count fits in size_t but exceeds the signed file
/// offset the resize calls take: 9223372043074251000 bytes.
constexpr size_t SMALLEST_UNMAPPABLE_N = 1518500251ULL;

}  // namespace

// The arithmetic. Nothing below this point until the wiring tests constructs
// a backend, so no sample count here can allocate anything.

TEST(CheckedPairCountTest, CountsNoPairsBelowTwoSamples) {
    // n - 1 underflows at n == 0, which the early return is there to avoid.
    EXPECT_EQ(checked_pair_count(0, "DenseStorage"), 0u);
    EXPECT_EQ(checked_pair_count(1, "DenseStorage"), 0u);
    EXPECT_EQ(checked_pair_count(2, "DenseStorage"), 1u);
    EXPECT_EQ(checked_pair_count(5, "DenseStorage"), 10u);
}

TEST(CheckedPairCountTest, AcceptsTheLargestSampleCountWhoseByteCountFits) {
    // The byte boundary itself: 2305843008139952128 pairs come to
    // 18446744065119617024 bytes, the last byte count size_t can hold. One
    // more sample crosses it, which the refusal test below pins. Nothing
    // constructs a backend for this count -- accepting it means asking for 18
    // exabytes.
    EXPECT_EQ(checked_pair_count(LARGEST_ACCEPTED_N, "DenseStorage"),
              2305843008139952128ULL);
}

TEST(CheckedPairCountTest, AcceptsTheSmallestSampleCountNoFileCanHold) {
    // Past the signed offset the file APIs size files with, but that is
    // MMapStorage's complaint, not this helper's: both of its checks pass.
    EXPECT_EQ(checked_pair_count(SMALLEST_UNMAPPABLE_N, "MMapStorage"),
              1152921505384281375ULL);
}

TEST(CheckedPairCountTest, RefusesTheSmallestSampleCountWhoseByteCountWraps) {
    EXPECT_EQ(overflow_message(
                  [] { checked_pair_count(SMALLEST_WRAPPING_BYTES_N, "DenseStorage"); }),
              "DenseStorage: sample count 2147483649 needs more than size_t "
              "bytes of storage");
}

TEST(CheckedPairCountTest, RefusesTheLargestRepresentablePairCountAtTheByteCheck) {
    // The pair count here is the largest size_t can hold, so a predicate that
    // discarded the remainder of its division would refuse this n for having
    // too many pairs. The refusal is right; the reason would not be.
    EXPECT_EQ(overflow_message(
                  [] { checked_pair_count(LARGEST_REPRESENTABLE_PAIRS_N, "DenseStorage"); }),
              "DenseStorage: sample count 6074001000 needs more than size_t "
              "bytes of storage");
}

TEST(CheckedPairCountTest, RefusesTheSmallestSampleCountWhosePairCountWraps) {
    EXPECT_EQ(overflow_message(
                  [] { checked_pair_count(SMALLEST_WRAPPING_PAIRS_N, "DenseStorage"); }),
              "DenseStorage: sample count 6074001001 has more pairs than "
              "size_t can hold");
}

TEST(CheckedPairCountTest, RefusesTheLargestSampleCountThereIs) {
    EXPECT_EQ(overflow_message([] {
                  checked_pair_count(std::numeric_limits<size_t>::max(),
                                     "DenseStorage");
              }),
              "DenseStorage: sample count 18446744073709551615 has more pairs "
              "than size_t can hold");
}

TEST(CheckedPairCountTest, NamesTheBackendItWasGiven) {
    EXPECT_EQ(overflow_message(
                  [] { checked_pair_count(SMALLEST_WRAPPING_PAIRS_N, "MMapStorage"); }),
              "MMapStorage: sample count 6074001001 has more pairs than "
              "size_t can hold");
}

// The wiring. Each case below constructs a backend, so each uses a sample
// count that stays cheap even if the guard it is testing is removed.

TEST(StorageOverflowTest, DenseRefusesASampleCountWhosePairCountWraps) {
    // Unguarded this allocates one double: SIZE_MAX * (SIZE_MAX - 1) wraps to
    // 2, so the buffer is a single pair while NumSamples() reports SIZE_MAX.
    // That is the cheapest possible regression and the reason this count, not
    // one of the large ones, is the one constructed here.
    const size_t n = std::numeric_limits<size_t>::max();
    EXPECT_EQ(overflow_message([n] { DenseStorage storage(n); }),
              "DenseStorage: sample count 18446744073709551615 has more pairs "
              "than size_t can hold");
}

TEST(StorageOverflowTest, MMapRefusesASampleCountWhosePairCountWraps) {
    // Unguarded this creates and maps an 8-byte file, for the same wrap as
    // above.
    const ProbePath probe("oecluster_overflow_pairs.bin");
    const std::string& path = probe.Path();
    const size_t n = std::numeric_limits<size_t>::max();
    EXPECT_EQ(overflow_message([&path, n] { MMapStorage storage(path, n); }),
              "MMapStorage: sample count 18446744073709551615 has more pairs "
              "than size_t can hold");
    EXPECT_FALSE(std::filesystem::exists(path));
}

TEST(StorageOverflowTest, MMapRefusesAByteCountPastTheFileOffsetLimit) {
    // This n passes both size_t checks. What it exceeds is the signed offset
    // the resize call takes, which would receive the byte count as a negative
    // size -- and be refused for it, so unguarded this costs an empty file
    // rather than the 9 exabytes the count asks for. It is the only test of
    // this third refusal, which has no arithmetic half: the limit it compares
    // against is the platform's.
    const ProbePath probe("oecluster_overflow_offset.bin");
    const std::string& path = probe.Path();
    const std::string message = overflow_message(
        [&path] { MMapStorage storage(path, SMALLEST_UNMAPPABLE_N); });
    EXPECT_NE(message.find("MMapStorage: sample count 1518500251 needs "
                           "9223372043074251000 bytes"),
              std::string::npos)
        << "message was: " << message;
    // The limit that completes the sentence is the platform's, so only the
    // wording around it is asserted here.
    EXPECT_NE(message.find(" byte file size limit"), std::string::npos)
        << "message was: " << message;
    EXPECT_FALSE(std::filesystem::exists(path));
}
