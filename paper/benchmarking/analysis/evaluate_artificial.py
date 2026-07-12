"""Benchmark Heinrichi by running
1) Froger
2) PAUP (parsimony)
3) NJ

All outputs are made in Newick and then compared using RF.
"""
import re
import csv
from collections import deque
import dendropy
import networkx as nx

from Bio.Phylo.TreeConstruction import *
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio import Phylo
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment

from pyfroger.list_apparatus import ListApparatus
from pyfroger.group_method import UnorientedStemma

from pyfroger.evaluate_reconstructed_trees import evaluate, average_sign_distance, infer_witnesses


# --- Heinrichi NEXUS matrix reader (AlignIO nexus parser rejects this file's
#     `format ... labels=left;` line, so read the matrix block directly) ---
def read_nexus_matrix(path):
    body = re.search(r"matrix(.*?);", open(path).read(), re.S | re.I).group(1)
    recs = []
    for line in body.splitlines():
        parts = line.split()
        if len(parts) == 2:
            recs.append(SeqRecord(Seq(parts[1]), id=parts[0], name=parts[0]))
    return MultipleSeqAlignment(recs)


# --- true_graph.txt (directed contaminated DAG) -> Newick.
#     BFS spanning tree from the archetype; secondary (contamination) edges
#     dropped; hypothetical nodes unary-suppressed / dangling-pruned;
#     ancestor-witnesses rendered as a pendant tip beside their clade. ---
def graph_to_newick(path, root_label="1"):
    lines = [l.rstrip("\n") for l in open(path)]
    n = int(lines[0])
    labels = lines[1:1 + n]
    M = [[int(x) for x in r.split()] for r in lines[1 + n:1 + n + n]]
    idx = {labels[i]: i for i in range(n)}
    children = {i: [j for j in range(n) if M[i][j] == 1] for i in range(n)}
    root = idx[root_label]
    kids = {i: [] for i in range(n)}
    seen = {root}
    q = deque([root])
    while q:
        u = q.popleft()
        for v in children[u]:
            if v not in seen:
                seen.add(v)
                kids[u].append(v)
                q.append(v)

    def build(v):
        subs = [s for s in (build(c) for c in kids[v]) if s is not None]
        name = labels[v]
        if name.isdigit():
            if not subs:
                return None                     # prune hypo leaf
            if len(subs) == 1:
                return subs[0]                  # suppress unary
            return "(" + ",".join(subs) + ")"
        return name if not subs else "(" + ",".join(subs + [name]) + ")"
    return build(root) + ";"


def graph_from_heinrichi(path: str) -> nx.DiGraph:
    """Ground-truth Heinrichi : en-tête N, puis N labels, puis matrice
    d'adjacence orientée N×N (1 = arc i→j, 0 = diagonale, -1 = absence)."""
    toks = open(path).read().split()
    n = int(toks[0])
    labels = toks[1:1 + n]
    m = list(map(int, toks[1 + n:1 + n + n * n]))
    g = nx.DiGraph()
    g.add_nodes_from(labels)
    for i in range(n):
        base = i * n
        for j in range(n):
            if m[base + j] == 1:
                g.add_edge(labels[i], labels[j])
    return g


graph_heinrichi =\
    graph_from_heinrichi("../../datasets/artificial_scribes/true_graph.txt")

# --- Heinrichi collation (tab-separated, witness-per-column, `.txt` sigla,
#     no locus IDs) -> Froger format: `locus,<siglum>,...` with a locus id
#     column and clean sigla. csv.writer quotes the comma/`;`-bearing readings;
#     empty cells stay empty (witness absent at that locus). ---


def to_froger_csv(src, dst):
    with open(src, newline="") as f:
        rows = list(csv.reader(f, delimiter="\t"))
    sigla = [c[:-4] if c.endswith(".txt") else c for c in rows[0]]
    with open(dst, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["locus"] + sigla)
        for i, r in enumerate(rows[1:], 1):
            w.writerow(["locus_%d" % i] + r)


aln = read_nexus_matrix("../../datasets/artificial_scribes/heinrichi.nex")

with open("ground_truth_heinrichi.dnd", "w") as f:
    f.write(graph_to_newick("../../datasets/artificial_scribes/true_graph.txt"))
with open("ground_truth_heinrichi.dnd") as f:
    true_tree = f.read()


# print("---- TREE FROM NEIGHBOR JOINING")

# calculator = DistanceCalculator('identity')
# constructor = DistanceTreeConstructor(calculator, 'nj')
# dm = calculator.get_distance(aln)


# tree = constructor.build_tree(aln)

# Phylo.draw_ascii(tree.as_phyloxml())

# tree_nj = re.sub(r':[0-9.]+', '', tree.format('newick'))

# # print(tree)

# # print("---- TREE FROM PARSIMONY")

# # scorer = ParsimonyScorer()
# # searcher = NNITreeSearcher(scorer)
# # constructor = ParsimonyTreeConstructor(searcher)
# # pars_tree = constructor.build_tree(aln)


# # Phylo.draw_ascii(pars_tree.as_phyloxml())
# # tree_pars = re.sub(r':[0-9.]+', '', pars_tree.format('newick'))

# # Phylo.draw_ascii(pars_tree.as_phyloxml())


# # tree_PARS = dendropy.Tree.get(data=tree_pars, schema="newick",
# #                               suppress_internal_node_taxa=False)

# # print("--- Comparison true vs parsimony")
# # print(average_sign_distance(true_source=graph_heinrichi,
# #                             recon_source=tree_pars,
# #                             witnesses=infer_witnesses(tree_PARS)))

# print("--- Comparison true vs NJ")
# tree_NJ = dendropy.Tree.get(data=tree_nj, schema="newick",
#                             suppress_internal_node_taxa=False)

# print(average_sign_distance(true_source=graph_heinrichi,
#                             recon_source=tree_nj,
#                             witnesses=infer_witnesses(tree_NJ)))

print("--- Comparison Froger")

to_froger_csv("../../datasets/artificial_scribes/heinrichi.csv",
              "../../datasets/artificial_scribes/heinrichi_froger.csv")

apparatus = ListApparatus.from_csv(
    "../../datasets/artificial_scribes/heinrichi.csv",
    ref_witness="Bd",          # no extant archetype in Heinrichi; set per your
                              # orientation convention
    delimiter=",")


stemma = UnorientedStemma(list_apparatus=apparatus,
                          min_frequency=10,
                          verbose=False
                          )
stemma.build_tree(anomaly_method="raw_cutoff",
                  resolve_anomaly_threshold=1)
# pprint(stemma.children_of)
# print(stemma.to_newick())
# stemma.display_tree()
# print(evaluate(true_newick=true_tree,
#                recon_newick=stemma.to_newick()))

# Froger to dendropy for get witnesses

print(stemma.display_variant_groups())

# print(stemma.to_extended_newick())

# froger_tree = dendropy.Tree.get(data=stemma.to_newick(), schema="newick",
#                                 suppress_internal_node_taxa=False)

# print(average_sign_distance(true_source=graph_heinrichi,
#                             recon_source=stemma.to_newick(),
#                             witnesses=infer_witnesses(froger_tree)))
