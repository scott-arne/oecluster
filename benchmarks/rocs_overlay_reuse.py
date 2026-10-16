"""What it costs ROCS to build an OEOverlay per search instead of reusing one.

``ROCSComparison`` builds a fresh ``OEOverlay`` for every search. It did not
always: one overlay per clone was reused across pairs, which made a score
depend on what the clone had scored before and moved an accepted input's
``cdist`` cell between 0.546 and a saturated 1.0. The rebuild is a correctness
fix and is not justified by its cost, but the cost is quoted in
``docs/python-api.md`` and in ``ROCSComparison.h``, so it has to be a number
anyone can reproduce.

The two arms here are the two strategies, written against the OpenEye API with
the same calls and the same options ``ROCSComparison::InitOverlay`` uses, so
both can run in one process. That is the point: a paired measurement across
two *builds* cannot alternate, and an unpaired one cannot separate the effect
from the machine warming up, which is how two earlier attempts at this
disagreed by a factor of 1.4.

What it does not measure: anything outside the overlay. Both arms run the same
number of searches, so the doubling that symmetrization costs is in neither
number.

Run::

    python benchmarks/rocs_overlay_reuse.py
    python benchmarks/rocs_overlay_reuse.py --repetitions 15 --conformers 50
"""

from __future__ import annotations

import argparse
import itertools
import json
import random
import statistics
import time

from openeye import oechem, oeomega, oeshape

# Flexible, drug-like-ish, and varied enough in size that the conformer search
# is doing real work. Listed here rather than loaded from a file so that the
# benchmark is self-contained and the fixture is reviewable in the diff.
SMILES = (
    "CCCCCCCC",
    "CCCCCCCCO",
    "c1ccc(CCCCO)cc1",
    "CCCCCCCCN",
    "c1ccc(CCCC)cc1",
    "CCCCOCCCC",
    "c1ccc(OCCCC)cc1",
    "CCCCCCCCCC",
)

COLOR_FORCE_FIELD = 1  # ImplicitMillsDean, the ROCSOptions default.


def build_molecules(conformers: int) -> list[oechem.OEMol]:
    """Omega ensembles for :data:`SMILES`, prepared as ROCSComparison prepares them.

    :param conformers: Maximum conformers per molecule.
    :returns: The prepared molecules.
    """
    omega = oeomega.OEOmega()
    omega.SetMaxConfs(conformers)
    omega.SetStrictStereo(False)

    color_opts = oeshape.OEColorOptions()
    color_opts.SetColorForceField(COLOR_FORCE_FIELD)
    prep = oeshape.OEOverlapPrep()
    prep.SetAssignColor(True)
    prep.SetUseHydrogens(True)
    prep.Initialize(color_opts)

    mols = []
    for smi in SMILES:
        mol = oechem.OEMol()
        if not oechem.OESmilesToMol(mol, smi) or not omega(mol):
            raise RuntimeError(f"could not build conformers for {smi}")
        prep.Prep(mol)
        mols.append(mol)
    return mols


def overlay_options() -> oeshape.OEOverlayOptions:
    """The options ``ROCSComparison::InitOverlay`` configures for ComboNorm."""
    color_opts = oeshape.OEColorOptions()
    color_opts.SetColorForceField(COLOR_FORCE_FIELD)
    opts = oeshape.OEOverlayOptions()
    opts.SetColorOptions(color_opts)
    return opts


def score_shared(mols, opts, pairs) -> float:
    """Both orientations of every pair on one overlay, reused throughout."""
    overlay = oeshape.OEOverlay(opts)
    total = 0.0
    for i, j in pairs:
        for ref, fit in ((i, j), (j, i)):
            overlay.SetupRef(mols[ref])
            score = oeshape.OEBestOverlayScore()
            overlay.BestOverlay(score, mols[fit])
            total += score.GetTanimotoCombo()
    return total


def score_fresh(mols, opts, pairs) -> float:
    """Both orientations of every pair, each on an overlay of its own."""
    total = 0.0
    for i, j in pairs:
        for ref, fit in ((i, j), (j, i)):
            overlay = oeshape.OEOverlay(opts)
            overlay.SetupRef(mols[ref])
            score = oeshape.OEBestOverlayScore()
            overlay.BestOverlay(score, mols[fit])
            total += score.GetTanimotoCombo()
    return total


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--conformers", type=int, default=20)
    parser.add_argument("--repetitions", type=int, default=9)
    parser.add_argument("--seed", type=int, default=20261016,
                        help="Seeds the coin that orders the two arms.")
    parser.add_argument("--json", action="store_true")
    args = parser.parse_args()

    mols = build_molecules(args.conformers)
    opts = overlay_options()
    pairs = list(itertools.combinations(range(len(mols)), 2))

    # One untimed repetition of each arm. Without it the first measurement
    # carries the toolkit's one-off setup, and whichever arm ran first wears
    # it -- the bias that made an earlier unpaired benchmark unusable.
    score_shared(mols, opts, pairs)
    score_fresh(mols, opts, pairs)

    rng = random.Random(args.seed)
    rows = []
    for repetition in range(args.repetitions):
        # Randomised order within each pair of runs, so that any residual
        # drift over the benchmark falls on both arms equally rather than on
        # whichever one is always second.
        arms = [("shared", score_shared), ("fresh", score_fresh)]
        if rng.random() < 0.5:
            arms.reverse()
        timings = {}
        for name, fn in arms:
            start = time.perf_counter()
            fn(mols, opts, pairs)
            timings[name] = time.perf_counter() - start
        rows.append({
            "repetition": repetition,
            "first": arms[0][0],
            "shared_ms_per_pair": 1000.0 * timings["shared"] / len(pairs),
            "fresh_ms_per_pair": 1000.0 * timings["fresh"] / len(pairs),
            "ratio_fresh_over_shared": timings["fresh"] / timings["shared"],
        })

    ratios = [row["ratio_fresh_over_shared"] for row in rows]
    summary = {
        "molecules": len(mols),
        "conformers_requested": args.conformers,
        "conformers_built": [mol.NumConfs() for mol in mols],
        "pairs": len(pairs),
        "searches_per_pair": 2,
        "repetitions": args.repetitions,
        "ratio_min": min(ratios),
        "ratio_median": statistics.median(ratios),
        "ratio_max": max(ratios),
    }

    if args.json:
        print(json.dumps({"rows": rows, "summary": summary}, indent=2))
        return

    print(f"{len(mols)} molecules, {summary['conformers_built']} conformers, "
          f"{len(pairs)} pairs, 2 searches per pair")
    print(f"{'rep':>4} {'first':>7} {'shared ms/pair':>15} {'fresh ms/pair':>14} "
          f"{'fresh/shared':>13}")
    for row in rows:
        print(f"{row['repetition']:>4} {row['first']:>7} "
              f"{row['shared_ms_per_pair']:>15.3f} {row['fresh_ms_per_pair']:>14.3f} "
              f"{row['ratio_fresh_over_shared']:>13.3f}")
    print(f"ratio fresh/shared: min {summary['ratio_min']:.3f}, "
          f"median {summary['ratio_median']:.3f}, max {summary['ratio_max']:.3f}")


if __name__ == "__main__":
    main()
