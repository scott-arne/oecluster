/**
 * @file MurckoScaffold.h
 * @brief Bemis-Murcko scaffold assignment and scaffold-identity clustering.
 */

#ifndef OECLUSTER_CLUSTERING_MURCKOSCAFFOLD_H
#define OECLUSTER_CLUSTERING_MURCKOSCAFFOLD_H

#include <cstddef>
#include <string>
#include <utility>
#include <vector>

#include "oecluster/clustering/ClusterTypes.h"

namespace OEChem {
class OEMolBase;
}

namespace OECluster {

/**
 * @brief Scaffold extraction level.
 *
 * ``Framework`` is the classic Bemis-Murcko scaffold: ring systems plus the
 * linkers that connect them, with sidechains removed. ``Generic`` reduces that
 * framework to its topology by setting every heavy atom to carbon and every
 * bond to single, so that scaffolds differing only in element or bond order
 * collapse together.
 */
enum class ScaffoldType {
    Framework,
    Generic
};

/**
 * @brief Options for Murcko scaffold assignment.
 */
struct MurckoOptions {
    /** @brief Extraction level; Framework is the classic Murcko scaffold. */
    ScaffoldType scaffold = ScaffoldType::Framework;
    /**
     * @brief Worker threads; 0 auto-detects hardware concurrency.
     *
     * A positive value is clamped to the molecule count and to a multiple of
     * the hardware concurrency, so an over-large request is bounded rather
     * than refused. Extraction runs serially regardless of this value when the
     * process memory-pool mode is not thread-safe; see ``murcko_scaffolds``.
     */
    size_t num_threads = 0;
};

/**
 * @brief Murcko scaffold clustering result.
 *
 * Clusters are scaffold identity classes. Label ``i`` is the class whose
 * scaffold is ``ClusterScaffolds()[i]``; labels are assigned by sorting the
 * distinct non-empty scaffold strings lexicographically, so the labeling is
 * canonical and independent of input order. Molecules with no ring system have
 * an empty scaffold string, which is excluded from that set: they form no
 * cluster and carry ``NOISE_LABEL``.
 */
class MurckoResult : public ClusteringResult {
public:
    MurckoResult() = default;
    MurckoResult(std::vector<ClusterLabel> labels, Clusters members,
                 std::vector<std::string> scaffolds,
                 std::vector<std::string> cluster_scaffolds)
        : ClusteringResult(std::move(labels), std::move(members)),
          scaffolds_(std::move(scaffolds)),
          cluster_scaffolds_(std::move(cluster_scaffolds)) {}

    /** @brief Per-item scaffold SMILES; empty for an acyclic molecule. */
    const std::vector<std::string>& Scaffolds() const { return scaffolds_; }
    /** @brief Scaffold naming each cluster; entry i belongs to label i. */
    const std::vector<std::string>& ClusterScaffolds() const {
        return cluster_scaffolds_;
    }

    std::string Method() const override { return "murcko"; }

private:
    std::vector<std::string> scaffolds_;
    std::vector<std::string> cluster_scaffolds_;
};

/**
 * @brief Compute the Bemis-Murcko scaffold of each molecule.
 *
 * Returns one canonical SMILES per input molecule, in input order. A molecule
 * with no ring system yields an empty string, which is the same "missing
 * scaffold" convention ``scaffold_agreement`` already consumes. A molecule that
 * *has* rings but from which no framework can be extracted is an error, not an
 * empty string.
 *
 * Input molecules are never modified: each is copied before extraction.
 *
 * Extraction is parallelized across molecules only when the process
 * memory-pool mode reports a thread-safe setting; otherwise the call runs
 * serially and returns the same answer.
 *
 * :param mols: Molecules to process; must not contain null pointers.
 * :param options: Extraction level and threading.
 * :returns: Scaffold SMILES, one per input molecule, in input order.
 * :raises ComparisonError: If ``mols`` is empty, holds a null pointer, or holds
 *     a ring-containing molecule whose framework cannot be extracted.
 * :raises std::invalid_argument: If ``scaffold`` is not a known enumerator.
 */
std::vector<std::string> murcko_scaffolds(
    const std::vector<OEChem::OEMolBase*>& mols,
    const MurckoOptions& options = MurckoOptions());

/**
 * @brief Cluster molecules by Bemis-Murcko scaffold identity.
 *
 * Two molecules share a cluster exactly when their scaffolds canonicalize to
 * the same non-empty SMILES. An acyclic molecule canonicalizes to ``""``,
 * which names no cluster: acyclic molecules are noise rather than a shared
 * "no scaffold" cluster, however many of them the input holds.
 *
 * :param mols: Molecules to cluster; must not contain null pointers.
 * :param options: Extraction level and threading.
 * :returns: MurckoResult with labels, clusters, and scaffold strings.
 * :raises ComparisonError: If ``mols`` is empty, holds a null pointer, or holds
 *     a ring-containing molecule whose framework cannot be extracted.
 * :raises std::invalid_argument: If ``scaffold`` is not a known enumerator.
 */
MurckoResult murcko_cluster(const std::vector<OEChem::OEMolBase*>& mols,
                            const MurckoOptions& options = MurckoOptions());

}  // namespace OECluster

#endif  // OECLUSTER_CLUSTERING_MURCKOSCAFFOLD_H
