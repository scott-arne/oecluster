/**
 * @file KernelOptions.h
 * @brief Translation of OECluster threading options into OEFP batch kernel options.
 *
 * Private to the comparison implementations: not installed, not exposed through SWIG.
 */

#ifndef OECLUSTER_COMPARISONS_KERNELOPTIONS_H
#define OECLUSTER_COMPARISONS_KERNELOPTIONS_H

#include <oefp/oefp.h>
#include <cstddef>

namespace OECluster {

/**
 * @brief Build OEFP batch kernel options from OECluster threading options.
 *
 * :param num_threads: Worker count; 0 lets OEFP choose.
 * :param chunk_size: Pairs per work unit; 0 selects the 256-pair default.
 * :returns: Populated kernel options.
 */
inline OEFP::BatchKernelOptions make_kernel_options(size_t num_threads, size_t chunk_size) {
    OEFP::BatchKernelOptions options;
    options.num_threads = num_threads;
    options.chunk_size = chunk_size > 0 ? chunk_size : 256;
    return options;
}

}  // namespace OECluster

#endif  // OECLUSTER_COMPARISONS_KERNELOPTIONS_H
