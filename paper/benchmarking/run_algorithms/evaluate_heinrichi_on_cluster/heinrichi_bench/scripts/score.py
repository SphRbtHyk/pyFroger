#!/usr/bin/env python3
"""
score.py — Roos average sign distance of reconstructed trees vs the true stemma.

Reads the 65-node Heinrichi adjacency file and one or more reconstructed trees
(Newick or NEXUS), and reports the skeleton sign distance for each.

Two details that move the number more than the choice of method does:

  * Zero-length branches are collapsed to polytomies (PAUP*'s collapse=yes
    default). The metric counts unweighted edges, so an arbitrary bifurcation
    where the data support none destroys ties that the true stemma has.
  * Scoring is restricted to the witnesses present in BOTH the tree and the
    true stemma. The matrix has 37 taxa, the stemma 34 (Da, I, J are absent
    from it): they take part in the search, not in the scoring.

Baselines are printed for calibration:
    star        every triple ties -> the floor any method must beat
    true        the true stemma against itself -> 100.0 by construction

Usage
-----
    python score.py --true heinrichi.graph \
        --tree iqtree=results/iqtree/heinrichi.treefile \
        --tree mrbayes=results/mrbayes/heinrichi.mb.nex.con.tre \
        --tree fastme=results/fastme/heinrichi.nwk \
        --out results/scores.tsv
"""

from __future__ import annotations

import argparse
import csv
import sys
from itertools import permutations
from pathlib import Path

import dendropy
import networkx as nx


# ── the metric ────────────────────────────────────────────────────────────

def roos_sign(x: int, y: int) -> int:
    """Sign of the comparison of two path lengths; -1 encodes 'unreachable'."""
    if x == -1 and y == -1:
        return 0
    if x == -1:
        return -1
    if y == -1:
        return 1
    if x < y:
        return -1
    if x > y:
        return 1
    return 0


def average_sign_distance(gt, wt, gp, wp) -> tuple[float, int]:
    """Skeleton sign distance over the witnesses common to both graphs."""
    common = sorted(set(wt) & set(wp))
    dt = dict(nx.all_pairs_shortest_path_length(gt))
    dp = dict(nx.all_pairs_shortest_path_length(gp))

    def d(dd, w, a, b):
        return dd[w[a]].get(w[b], -1)

    s = tot = 0.0
    for i, j, k in permutations(common, 3):
        st = roos_sign(d(dt, wt, i, j), d(dt, wt, i, k))
        sp = roos_sign(d(dp, wp, i, j), d(dp, wp, i, k))
        s += 1 - abs(st - sp) / 2.0
        tot += 1
    return (round(100.0 * s / tot, 3) if tot else float("nan")), len(common)


# ── the true stemma ───────────────────────────────────────────────────────

def graph_from_heinrichi(path: Path) -> tuple[nx.Graph, dict[str, str]]:
    """Header N, then N labels, then an N x N adjacency matrix.

    The skeleton metric is undirected: we discard arc orientation here.
    """
    toks = Path(path).read_text().split()
    n = int(toks[0])
    labels = toks[1:1 + n]
    m = list(map(int, toks[1 + n:1 + n + n * n]))

    g = nx.Graph()
    g.add_nodes_from(labels)
    for i in range(n):
        for j in range(n):
            if m[i * n + j] == 1:
                g.add_edge(labels[i], labels[j])

    # Witnesses are the alphabetic labels; the numeric ones are hypothetical.
    wit = {lab: lab for lab in labels if lab[0].isalpha()}
    return g, wit


# ── reconstructed trees ───────────────────────────────────────────────────

def load_tree(path: Path, collapse_eps: float) -> dendropy.Tree:
    """Read a tree, whatever the tool decided to write.

    PAUP's `contree ... treefile=` emits a bare Newick string with no trailing
    semicolon, which DendroPy's Newick reader rejects. IQ-TREE writes proper
    Newick. MrBayes writes NEXUS. So: sniff, repair, and if it still fails,
    report every attempt rather than only the last one (reporting only the last
    makes a Newick failure masquerade as "Expecting '#NEXUS', but found '('").
    """
    raw = path.read_text(encoding="utf-8", errors="replace").strip()
    if not raw:
        raise SystemExit(f"error: {path} is empty")

    schema = "nexus" if raw.lstrip().upper().startswith("#NEXUS") else "newick"

    print(f"--- Detected schema: {schema}")

    if schema == "newick" and not raw.rstrip().endswith(";"):
        raw = raw.rstrip() + ";"

    errors = []
    t = None
    for sch in (schema, "nexus" if schema == "newick" else "newick"):
        try:
            # suppress_internal_node_taxa=True is essential: PAUP writes
            # majority-rule support values as internal node labels (100, 0,
            # ...), and with taxa suppressed=False DendroPy tries to make a
            # Taxon out of each, then chokes on the duplicates. The labels are
            # still readable as nd.label, so nothing is lost.
            t = dendropy.Tree.get(
                data=raw, schema=sch,
                preserve_underscores=True,
                suppress_internal_node_taxa=True,
            )
            break
        except Exception as exc:                     # noqa: BLE001
            errors.append(f"  as {sch}: {exc}")
            t = None

    if t is None:
        msg = "\n".join(errors)
        raise SystemExit(
            f"error: cannot parse {path}\n{msg}\n"
            f"  first 120 chars: {raw[:120]!r}"
        )

    if collapse_eps is not None and collapse_eps >= 0:
        for edge in list(t.postorder_edge_iter()):
            if edge.tail_node is None:
                continue                     # root edge
            if edge.head_node.is_leaf():
                continue                     # never collapse a terminal
            if edge.length is not None and edge.length <= collapse_eps:
                edge.collapse()
    return t


def tree_to_graph(tree: dendropy.Tree,
                  known: set[str] | None = None) -> tuple[nx.Graph, dict[str, str]]:
    """Graph of the tree, plus the map witness-label -> node id.

    Leaves are always witnesses. Internal nodes are witnesses ONLY if their
    label names one — which is how a stemmatic method (Froger, RHM, Semstem)
    represents a manuscript copied directly from another. For a phylogenetic
    tree the internal labels are support values (PAUP writes 100, 0, ...) and
    must not be mistaken for taxa; passing `known` filters them out.
    """
    g = nx.Graph()
    ids = {}
    for k, nd in enumerate(tree.preorder_node_iter()):
        ids[nd] = k
        g.add_node(k)
    for nd in tree.preorder_node_iter():
        if nd.parent_node is not None:
            g.add_edge(ids[nd.parent_node], ids[nd])

    def label_of(nd):
        lab = nd.taxon.label if nd.taxon is not None else nd.label
        return lab.strip().replace(" ", "_") if lab else None

    wit = {}
    for nd in tree.leaf_node_iter():
        if (lab := label_of(nd)):
            wit[lab] = ids[nd]

    if known:
        for nd in tree.preorder_internal_node_iter():
            lab = label_of(nd)
            if lab and lab in known and lab not in wit:
                wit[lab] = ids[nd]

    return g, wit


def star_graph(witnesses) -> tuple[nx.Graph, dict[str, str]]:
    g = nx.Graph()
    g.add_node("__hub__")
    wit = {}
    for w in witnesses:
        g.add_edge("__hub__", w)
        wit[w] = w
    return g, wit


# ── main ──────────────────────────────────────────────────────────────────

def main(argv=None) -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--true", required=True, type=Path, dest="true_graph")
    p.add_argument("--tree", action="append", default=[], metavar="NAME=PATH",
                   help="repeatable; e.g. --tree iqtree=results/x.treefile")
    p.add_argument("--collapse-eps", type=float, default=1e-6,
                   help="collapse internal branches <= this (default 1e-6); "
                        "pass -1 to disable")
    p.add_argument("--out", type=Path, default=None, help="TSV output")
    args = p.parse_args(argv)

    gt, wt = graph_from_heinrichi(args.true_graph)
    print(f"true stemma: {gt.number_of_nodes()} nodes, "
          f"{len(wt)} witnesses", file=sys.stderr)

    rows = []

    # baselines
    sc, n = average_sign_distance(gt, wt, gt, wt)
    rows.append({"method": "true_stemma", "skeleton": sc,
                 "n_scored": n, "n_leaves": len(wt), "note": "sanity check"})
    sg, sw = star_graph(sorted(wt))
    sc, n = average_sign_distance(gt, wt, sg, sw)
    rows.append({"method": "star", "skeleton": sc,
                 "n_scored": n, "n_leaves": len(sw),
                 "note": "floor: every triple ties"})

    for spec in args.tree:
        if "=" not in spec:
            raise SystemExit(f"error: --tree expects NAME=PATH, got {spec!r}")
        name, path = spec.split("=", 1)
        path = Path(path)
        if not path.exists():
            print(f"warning: {name}: {path} missing, skipped", file=sys.stderr)
            continue

        tree = load_tree(path, args.collapse_eps)
        gp, wp = tree_to_graph(tree, known=set(wt))

        missing = sorted(set(wt) - set(wp))
        extra = sorted(set(wp) - set(wt))
        sc, n = average_sign_distance(gt, wt, gp, wp)

        note = []
        if missing:
            note.append(f"absent from tree: {','.join(missing)}")
        if extra:
            note.append(f"not scored: {','.join(extra)}")
        rows.append({"method": name, "skeleton": sc, "n_scored": n,
                     "n_leaves": len(wp), "note": "; ".join(note)})

    fields = ["method", "skeleton", "n_scored", "n_leaves", "note"]
    w = csv.DictWriter(sys.stdout, fieldnames=fields, delimiter="\t")
    w.writeheader()
    w.writerows(rows)

    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        with args.out.open("w", newline="") as fh:
            w2 = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
            w2.writeheader()
            w2.writerows(rows)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())