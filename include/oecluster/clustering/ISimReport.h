/**
 * @file ISimReport.h
 * @brief iSIM set similarity and the approximate fingerprint-native cluster report.
 *
 * iSIM values are union-weighted: summed pairwise intersections over summed
 * pairwise unions. They are exact for that quantity and are not the mean of
 * per-pair Tanimoto values.
 */

#ifndef OECLUSTER_CLUSTERING_ISIMREPORT_H
#define OECLUSTER_CLUSTERING_ISIMREPORT_H

#include <cstddef>
#include <limits>
#include <string>
#include <vector>

#include "oefp/batch.h"
#include "oecluster/clustering/ClusterReport.h"
#include "oecluster/clustering/ClusterTypes.h"

namespace OECluster {

/// Options for isim(). metric accepts only "tanimoto"; it is reserved so
/// metrics whose iSIM is the exact mean pair similarity can be added later.
struct ISimOptions {
    std::string metric = "tanimoto";
};

/**
 * @brief iSIM Tanimoto similarity of a set of binary fingerprints.
 *
 * NaN for fewer than two fingerprints; 1.0 when every fingerprint is all-zero.
 *
 * @throws std::invalid_argument for a metric other than "tanimoto", a
 *         non-empty batch of zero-width fingerprints, fingerprints of
 *         2^31 or more bits, or 2^32 or more fingerprints.
 */
double isim(const OEFP::OEFPBatch& fingerprints, const ISimOptions& options = ISimOptions());

}  // namespace OECluster

#endif  // OECLUSTER_CLUSTERING_ISIMREPORT_H
