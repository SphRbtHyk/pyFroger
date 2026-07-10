#!/usr/bin/env python3
"""Run IQ-TREE (ML), FastME (BME distance) and MrBayes (Bayesian) on a set of
simulated manuscript collations and emit reconstructed Newick Trees.
"""
from __future__ import annotations
import argparse
import csv
import glob
import os
import subprocess
import sys
import time
import dendropy
import stemma_input as sd


def run_tool(runtime, image, tool, args, workdir):
    workdir = os.path.abspath(workdir)
    if runtime == "docker":
        cmd = ["docker", "run", "--rm", "-u", f"{os.getuid()}:{os.getgid()}",
               "-v", f"{workdir}:/data", "-w", "/data", image] + args
    elif runtime in ("apptainer", "singularity"):
        cmd = [runtime, "exec", "--pwd", "/data", "-B", f"{workdir}:/data",
               image, tool] + args
    else:
        raise ValueError(runtime)
    p = subprocess.run(cmd, capture_output=True, text=True)
    return p


def reroot_relabel(in_nwk, out_nwk, rev, outgroup_id, schema="newick"):
    tns = dendropy.TaxonNamespace()
    if schema == "con":  # MrBayes .con.tre: pick the majority-rule tree
        tl = dendropy.TreeList.get(
            path=in_nwk, schema="nexus", taxon_namespace=tns)
        t = next(
            (x for x in tl if x.label and "majrule" in x.label.lower()), tl[0])
    else:
        t = dendropy.Tree.get(
            path=in_nwk, schema="newick", taxon_namespace=tns)
    node = t.find_node_with_taxon_label(outgroup_id)
    if node is not None and node.edge.head_node is not t.seed_node:
        t.reroot_at_edge(node.edge, update_bipartitions=True,
                         suppress_unifurcations=True)
    for lf in t.leaf_node_iter():
        if lf.taxon is not None:
            lf.taxon.label = rev.get(lf.taxon.label, lf.taxon.label)
    t.write(path=out_nwk, schema="newick", suppress_rooting=False,
            suppress_internal_node_labels=True)


def _strip_invariant(ids_seqs):
    """Drop constant columns so IQ-TREE's +ASC precondition holds.
    Accepts a dict {id: seq}; returns (filtered_dict, n_removed)."""
    items = list(ids_seqs.items())
    ids = [i for i, _ in items]
    seqs = [list(s) for _, s in items]
    keep = [c for c in range(len(seqs[0]))
            if len({row[c] for row in seqs}) > 1]
    filtered = {ids[r]: "".join(seqs[r][c] for c in keep)
                for r in range(len(ids))}
    return filtered, len(seqs[0]) - len(keep)


def method_iqtree(prep, work, cfg):
    ids_seqs = prep["ids_seqs"]
    if "+ASC" in cfg["iqtree_model"].upper():
        ids_seqs, n_removed = _strip_invariant(ids_seqs)
        if n_removed:
            print(f"[iqtree] stripped {n_removed} invariant sites for +ASC",
                  flush=True)
        if not ids_seqs or len(next(iter(ids_seqs.values()))) == 0:
            raise ValueError("no variable sites remain after +ASC filtering")
    sd.write_phylip_aln(os.path.join(work, "aln.phy"), ids_seqs)
    args = ["-s", "aln.phy", "-st", "MORPH", "-m", cfg["iqtree_model"],
            "-o", prep["outgroup_id"], "-pre", "iqtree",
            "-T", str(cfg["iqtree_threads"]),
            "--seed", str(cfg["seed"]), "-redo", "-quiet"]
    p = run_tool(cfg["runtime"], cfg["iqtree_image"], "iqtree", args, work)
    tree = os.path.join(work, "iqtree.treefile")
    return p, tree, "newick"


def method_fastme(prep, work, cfg):
    sd.write_phylip_matrix(os.path.join(work, "dist.phy"),
                           prep["id_order"], prep["D"])
    # BME+NNI default; -s adds SPR
    args = ["-i", "dist.phy", "-o", "fastme.nwk", "-s"]
    p = run_tool(cfg["runtime"], cfg["fastme_image"], "fastme", args, work)
    return p, os.path.join(work, "fastme.nwk"), "newick"


def method_mrbayes(prep, work, cfg):
    sd.write_nexus_standard(os.path.join(work, "run.nex"), prep["ids_seqs"],
                            prep["gmax"], prep["outgroup_id"],
                            cfg["ngen"], cfg["samplefreq"],
                            stopval=cfg.get("mb_stopval", 0.1),
                            diagnfreq=cfg.get("mb_diagnfreq", 10000))
    p = run_tool(cfg["runtime"], cfg["mrbayes_image"], "mb", ["run.nex"], work)
    return p, os.path.join(work, "run.nex.con.tre"), "con"


METHODS = {"iqtree": method_iqtree,
           "fastme": method_fastme,
           "mrbayes": method_mrbayes}


# ---- driver -----------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "inputs", help="glob for collation FASTAs, e.g. 'collations/*.fasta'")
    ap.add_argument("--out", default="bench")
    ap.add_argument("--methods", default="iqtree,fastme,mrbayes")
    ap.add_argument("--outgroup", default="w1")
    ap.add_argument("--runtime", choices=["docker", "apptainer", "singularity"],
                    default="docker")
    ap.add_argument("--iqtree-threads", default="AUTO",
                    help="pass $SLURM_CPUS_PER_TASK under SLURM")
    ap.add_argument("--metrics", default=None,
                    help="metrics CSV path")
    ap.add_argument("--iqtree-image", default="evolbioinfo/iqtree:v3.1.1")
    ap.add_argument("--fastme-image", default="evolbioinfo/fastme:v2.1.6.4")
    ap.add_argument("--mrbayes-image", default="evolbioinfo/mrbayes:v3.2.7")
    ap.add_argument("--iqtree-model", default="MK+ASC",
                    help="use MK if the matrix retains constant sites")
    ap.add_argument("--ngen", type=int, default=1000)  # default=1_000_000)
    ap.add_argument("--samplefreq", type=int, default=100)  # ,default=1000)
    ap.add_argument("--mb-stopval", type=float, default=0.01,
                    help="ASDSF stop-rule threshold; 0 disables (runs to ngen)")
    ap.add_argument("--mb-diagnfreq", type=int, default=10000,
                    help="MrBayes convergence-diagnostic frequency")

    ap.add_argument("--seed", type=int, default=12345)
    ap.add_argument("--force", action="store_true",
                    help="recompute existing outputs")
    a = ap.parse_args()

    methods = [m.strip() for m in a.methods.split(",") if m.strip()]
    cfg = {k: getattr(a, k) for k in
           ("runtime", "iqtree_model", "ngen", "samplefreq", "seed",
            "iqtree_threads", "mb_stopval", "mb_diagnfreq")}
    cfg.update({"iqtree_image": a.iqtree_image, "fastme_image": a.fastme_image,
                "mrbayes_image": a.mrbayes_image})

    files = sorted(glob.glob(a.inputs))
    if not files:
        sys.exit(f"no inputs matched {a.inputs!r}")
    for m in methods:
        os.makedirs(os.path.join(a.out, "trees", m), exist_ok=True)

    metrics_path = a.metrics or os.path.join(a.out, "metrics.csv")
    os.makedirs(os.path.dirname(os.path.abspath(metrics_path)), exist_ok=True)
    new_file = not os.path.exists(metrics_path)
    mfh = open(metrics_path, "a", newline="")
    writer = csv.writer(mfh)
    if new_file:
        writer.writerow(["method", "stem", "n_taxa", "n_chars",
                         "mean_hamming", "wall_s", "status"])

    for fasta in files:
        parent = os.path.basename(os.path.dirname(fasta))
        base = os.path.splitext(os.path.basename(fasta))[0]
        stem = f"{parent}_{base}"
        try:
            prep = sd.prepare(fasta, a.outgroup)
        except Exception as e:
            sys.stderr.write(f"[skip] {stem}: prepare failed: {e}\n")
            continue
        for m in methods:
            out_nwk = os.path.join(a.out, "trees", m, f"{stem}.nwk")
            if os.path.exists(out_nwk) and not a.force:
                continue
            work = os.path.join(a.out, "work", m, stem)
            os.makedirs(work, exist_ok=True)
            t0 = time.time()
            status = "ok"
            try:
                p, raw, schema = METHODS[m](prep, work, cfg)
                if p.returncode != 0 or not os.path.exists(raw):
                    status = f"toolerr(rc={p.returncode})"
                    sys.stderr.write(f"[{m}:{stem}] {p.stderr[-400:]}\n")
                else:
                    reroot_relabel(raw, out_nwk, prep["rev"],
                                   prep["outgroup_id"], schema)
            except Exception as e:
                status = f"exc:{type(e).__name__}"
                sys.stderr.write(f"[{m}:{stem}] {e}\n")
            writer.writerow([m, stem, prep["n_taxa"], prep["n_chars"],
                             f"{prep['mean_hamming']:.4f}",
                             f"{time.time()-t0:.1f}", status])
            mfh.flush()
            sys.stderr.write(
                f"[{m}:{stem}] {status} ({time.time()-t0:.1f}s)\n")
    mfh.close()


if __name__ == "__main__":
    main()
