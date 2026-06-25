"""Rank and select multiple representatives from one cluster."""

from openeye import oechem

import oecluster


MOLECULES = (
    ("benzene", "c1ccccc1", 0.1, 0.3, "aryl"),
    ("toluene", "Cc1ccccc1", 0.0, 0.8, "aryl"),
    ("phenol", "Oc1ccccc1", 0.2, 0.6, "aryl"),
    ("anisole", "COc1ccccc1", 0.1, 0.7, "aryl"),
    ("ethanol", "CCO", 0.0, 0.2, "alcohol"),
    ("propanol", "CCCO", 0.0, 0.4, "alcohol"),
)


def molecules_from_records(records):
    """Build titled OpenEye molecules from example metadata records."""
    molecules = []
    for title, smiles, *_ in records:
        mol = oechem.OEGraphMol()
        if not oechem.OESmilesToMol(mol, smiles):
            raise ValueError(f"Could not parse SMILES for {title}: {smiles}")
        mol.SetTitle(title)
        molecules.append(mol)
    return molecules


def describe_selection(label, selected, titles):
    """Print selected representatives in a compact, scannable format."""
    members = ", ".join(
        f"{titles[item.member]}(rank={item.metrics.representative_rank})"
        for item in selected
    )
    print(f"{label}: {members}")


def main():
    """Rank representatives and compare score versus diversity selection."""
    molecules = molecules_from_records(MOLECULES)
    titles = [record[0] for record in MOLECULES]
    liabilities = [record[2] for record in MOLECULES]
    priorities = [record[3] for record in MOLECULES]
    scaffolds = [record[4] for record in MOLECULES]

    dm = oecluster.pdist(
        molecules,
        "fingerprint",
        fp_type="morgan",
        metric="tanimoto",
        num_threads=1,
    )
    cluster = tuple(range(len(molecules)))

    ranked = oecluster.rank_representatives(
        cluster,
        dm,
        method="weighted_medoid",
        threshold=0.55,
        alpha=1.0,
        beta=0.5,
        gamma=0.4,
        liability_penalties=liabilities,
        priority_scores=priorities,
        scaffold_labels=scaffolds,
    )
    print("ranked representatives:")
    for item in ranked[:3]:
        print(
            f"  {item.metrics.representative_rank}. {titles[item.member]} "
            f"score={item.score:.3f} "
            f"mean_distance={item.metrics.mean_distance_to_cluster:.3f} "
            f"scaffold_purity={item.metrics.scaffold_purity:.3f}"
        )

    score_selected = oecluster.select_representatives(
        cluster,
        dm,
        k=2,
        method="weighted_medoid",
        selection="score",
        threshold=0.55,
        alpha=1.0,
        beta=0.5,
        gamma=0.4,
        liability_penalties=liabilities,
        priority_scores=priorities,
        scaffold_labels=scaffolds,
    )
    diversity_selected = oecluster.select_representatives(
        cluster,
        dm,
        k=2,
        method="weighted_medoid",
        selection="diversity",
        threshold=0.55,
        alpha=1.0,
        beta=0.5,
        gamma=0.4,
        liability_penalties=liabilities,
        priority_scores=priorities,
        scaffold_labels=scaffolds,
    )

    describe_selection("score selection", score_selected, titles)
    describe_selection("diversity selection", diversity_selected, titles)


if __name__ == "__main__":
    main()
