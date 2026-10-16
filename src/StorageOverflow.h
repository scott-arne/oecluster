/**
 * @file StorageOverflow.h
 * @brief Buffer-length arithmetic shared by the contiguous storage backends.
 *
 * Private to the storage implementation: not installed, not exposed through SWIG.
 *
 * Declaring the helper here gives it external linkage, which it did not have
 * while it sat in an anonymous namespace in StorageBackend.cpp. That is the
 * trade being made rather than an oversight: the arithmetic below decides
 * whether a sample count is representable, the counts that exercise it are
 * ones no test can afford to allocate, and reaching it from the test suite is
 * what makes them affordable. DescriptorBuild.h declares its private helpers
 * the same way, and one shared pattern is worth more than one fewer symbol.
 */

#ifndef OECLUSTER_STORAGEOVERFLOW_H
#define OECLUSTER_STORAGEOVERFLOW_H

#include <cstddef>

namespace OECluster {

/**
 * @brief Pair count for n items, refusing a count whose buffer cannot be sized.
 *
 * Both contiguous backends derive a buffer length from the sample count in two
 * steps, and either can wrap. A wrapped length is not a short allocation the
 * caller can notice: NumSamples() keeps reporting n while the buffer holds
 * whatever the wrapped product came to, so Set, Get and every reader that
 * walks Data() index past the end. Refusing here gives those readers the
 * precondition they already assume.
 *
 * Declared apart from the backends so that the arithmetic can be tested
 * without constructing one. The counts that matter here are the ones a
 * regressed guard would try to allocate, and a test that allocated them would
 * exhaust the machine rather than report the regression.
 *
 * :param n: Number of items the backend was asked to store.
 * :param backend: Storage class name, used in the message.
 * :returns: n * (n - 1) / 2, which is 0 for fewer than two items.
 * :raises std::overflow_error: If that pair count, or the byte count it
 *     implies, exceeds what size_t can represent.
 */
size_t checked_pair_count(size_t n, const char* backend);

}  // namespace OECluster

#endif  // OECLUSTER_STORAGEOVERFLOW_H
