/**
 * @file MCSComparisonSnapshotAccess.h
 * @brief Test-only read access to an MCSComparison's molecule snapshots.
 *
 * Private to the comparison implementations: not installed, not exposed through SWIG.
 */

#ifndef OECLUSTER_COMPARISONS_MCSCOMPARISONSNAPSHOTACCESS_H
#define OECLUSTER_COMPARISONS_MCSCOMPARISONSNAPSHOTACCESS_H

#include <vector>
#include <oechem.h>

#include "oecluster/comparisons/MCSComparison.h"

namespace OECluster {

/// Test-only accessor for the molecule snapshots a comparison holds.
///
/// ``SharedData`` is defined in the .cpp, so friendship alone does not let a
/// test dereference ``shared_``: the type is incomplete there. This declares
/// the single observation the test needs, and the definition sits beside
/// ``SharedData`` where it is complete.
///
/// The declaration lives in this private header rather than in
/// ``include/oecluster/comparisons/MCSComparison.h`` because that header is
/// installed. A declaration there is a callable, supported-looking entry point
/// for any native consumer of the shipped library, handing out raw pointers
/// into another object's private storage. Keeping it out of the installed tree
/// also retires two name-based guards that were protecting the same hook one
/// surface at a time: the ``#ifndef SWIG`` around the struct, since
/// ``swig/oecluster.i`` ``%include``s only headers from ``include/oecluster``
/// and never any from ``src``, and the ``EXCLUDE_SYMBOLS`` entry in
/// ``docs/conf.py``, since Doxygen's ``INPUT`` is ``include/oecluster``.
/// Neither can now be defeated by renaming the struct.
///
/// The mangled symbol still exists in the shipped archive; only a separate
/// test-only build configuration would remove it, and that would assert the
/// thread-safety guarantee against a binary we do not ship. What this removes
/// is any declaration a consumer could reach without hand-writing an extern.
struct MCSComparisonSnapshotAccess {
    /// Returns a pointer to every molecule snapshot, in storage order.
    ///
    /// The pointer is typed rather than ``const void*``, and that is the whole
    /// safeguard rather than a stylistic preference. A ``void*`` return needs an
    /// explicit cast in the definition, and that cast accepts
    /// ``static_cast<const void*>(&mol)`` -- the address of the ``shared_ptr``
    /// slot rather than of the molecule it owns -- as quietly as it accepts
    /// ``mol.get()``. The accessor would then report the snapshot vector's own
    /// element addresses, which differ between any two clones no matter what the
    /// snapshots point at, and the disjointness test would pass over an aliasing
    /// ``Clone()``. ``const OEMol*`` admits no conversion from the
    /// ``shared_ptr`` slot, so the same slip is a compile error. Do not widen
    /// this back to ``void*``.
    static std::vector<const OEChem::OEMol*> SnapshotAddresses(const MCSComparison& cmp);
};

}  // namespace OECluster

#endif  // OECLUSTER_COMPARISONS_MCSCOMPARISONSNAPSHOTACCESS_H
