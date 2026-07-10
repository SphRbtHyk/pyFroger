"""Run BioPython maximum-parsimony on the synthetic apparatus alignments.

Cluster (SLURM array) version: one array task processes one alignment,
selected by  IDX = SLURM_ARRAY_TASK_ID + BASE  against a deterministic,
sorted, scenario-filtered manifest. Falls back to a serial loop when run
outside an array (local testing).
"""
import glob
import os
import re
import sys

from Bio import AlignIO
from Bio.Phylo.TreeConstruction import (
    ParsimonyScorer,
    NNITreeSearcher,
    ParsimonyTreeConstructor,
)

# Paths from the environment so they change per cluster without editing code.
INPUT_FOLDER = os.environ.get("DATADIR", "../../datasets/synthetic/generated_trees")
OUTPUT_FOLDER = os.environ.get("OUTDIR", "../experiment_outputs/parsimony/")

SCENARIOS = {"scenario_1", "scenario_2", "scenario_3", "scenario_6"}

# BASE lets absolute indices exceed SLURM's MaxArraySize (index cap, default
# 1001): submit tranches with --export=ALL,BASE=N, keeping SLURM_ARRAY_TASK_ID
# in [0, 999].
BASE = int(os.environ.get("BASE", "0"))


def build_manifest():
    """Deterministic, sorted, scenario-filtered list of input alignments.

    Sorting matters: glob order is not stable, so an unsorted list would make
    the IDX -> file mapping differ between the dump step and the run step.
    """
    paths = sorted(glob.glob(os.path.join(INPUT_FOLDER, "*", "*.fasta")))
    return [p for p in paths
            if os.path.basename(os.path.dirname(p)) in SCENARIOS]


def load_manifest():
    """Prefer a pre-dumped MANIFEST file (exact index parity); else build it."""
    manifest_path = os.environ.get("MANIFEST")
    if manifest_path and os.path.exists(manifest_path):
        with open(manifest_path) as f:
            return [line.strip() for line in f if line.strip()]
    return build_manifest()


def run_one(apparatus_path):
    scenario = os.path.basename(os.path.dirname(apparatus_path))
    experiment_id = os.path.splitext(os.path.basename(apparatus_path))[0]
    out_path = os.path.join(OUTPUT_FOLDER, f"{scenario}_{experiment_id}.nwk")

    # Idempotent: skip work already done so a failed array can be re-run cheaply.
    if os.path.exists(out_path):
        print(f"[skip] {scenario}/{experiment_id} (exists)", flush=True)
        return

    print(f"[run]  {scenario}/{experiment_id}", flush=True)
    aln = AlignIO.read(apparatus_path, "fasta")

    constructor = ParsimonyTreeConstructor(NNITreeSearcher(ParsimonyScorer()))
    pars_tree = constructor.build_tree(aln)

    tree_pars = re.sub(r":[0-9.]+", "", pars_tree.format("newick"))
    with open(out_path, "w") as f:
        f.write(tree_pars)


def main():
    if "--dump-manifest" in sys.argv:
        for path in build_manifest():
            print(path)
        return

    os.makedirs(OUTPUT_FOLDER, exist_ok=True)
    manifest = load_manifest()

    task_id = os.environ.get("SLURM_ARRAY_TASK_ID")
    if task_id is None:
        # Local / serial mode: process everything.
        print(f"[serial] {len(manifest)} alignments", flush=True)
        for path in manifest:
            run_one(path)
        return

    idx = int(task_id) + BASE
    if idx >= len(manifest):
        print(f"[noop] IDX={idx} >= manifest size {len(manifest)}", flush=True)
        return
    run_one(manifest[idx])


if __name__ == "__main__":
    main()