"""Cluster a small set of molecules from inline SMILES."""

from openeye import oechem

import oecluster


MOLECULES = (
    ("benzene", "c1ccccc1"),
    ("toluene", "Cc1ccccc1"),
    ("phenol", "Oc1ccccc1"),
    ("ethanol", "CCO"),
    ("propanol", "CCCO"),
    ("acetic_acid", "CC(=O)O"),
)


def molecules_from_smiles(records):
    """Build titled OpenEye molecules from ``(title, smiles)`` records."""
    molecules = []
    for title, smiles in records:
        mol = oechem.OEGraphMol()
        if not oechem.OESmilesToMol(mol, smiles):
            raise ValueError(f"Could not parse SMILES for {title}: {smiles}")
        mol.SetTitle(title)
        molecules.append(mol)
    return molecules


def main():
    """Run a minimal clustering workflow."""
    molecules = molecules_from_smiles(MOLECULES)
    titles = [mol.GetTitle() for mol in molecules]

    dm = oecluster.pdist(
        molecules,
        "fingerprint",
        fp_type="morgan",
        metric="tanimoto",
        num_threads=1,
    )
    result = oecluster.butina(dm, threshold=0.55)

    print(f"clusters: {len(result.clusters)}")
    print("representatives:")
    for cluster_id, cluster in enumerate(result.clusters):
        ranked = oecluster.rank_representatives(cluster, dm, method="medoid")
        representative = ranked[0]
        member_titles = ", ".join(titles[index] for index in cluster)
        print(
            f"  cluster {cluster_id}: "
            f"size={len(cluster)} "
            f"medoid={titles[representative.member]} "
            f"radius={representative.metrics.cluster_radius:.3f} "
            f"members=[{member_titles}]"
        )


if __name__ == "__main__":
    main()
