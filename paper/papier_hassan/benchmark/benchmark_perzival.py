"""Benchmark Perzival by running
1) Froger
2) PAUP
3) NJ

All outputs are made in Newick and then compared using RF.
"""
import re
from Bio import AlignIO

from Bio.Phylo.TreeConstruction import *
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio import Phylo

from pyfroger.list_apparatus import ListApparatus
from pyfroger.group_method import UnorientedStemma

from pyfroger.scoring import evaluate, diagnose


aln = AlignIO.read(
    '../../datasets/artificial_scribes/parzival.fasta', 'fasta')


print("---- TREE FROM NEIGHBOR JOINING")

calculator = DistanceCalculator('identity')
constructor = DistanceTreeConstructor(calculator, 'nj')
dm = calculator.get_distance(aln)


tree = constructor.build_tree(aln)

Phylo.draw_ascii(tree.as_phyloxml())

tree_nj = re.sub(r':[0-9.]+', '', tree.format('newick'))

# print(tree)

print("---- TREE FROM PARSIMONY")

scorer = ParsimonyScorer()
searcher = NNITreeSearcher(scorer)
constructor = ParsimonyTreeConstructor(searcher)
pars_tree = constructor.build_tree(aln)

# Phylo.draw_ascii(pars_tree.as_phyloxml())
tree_pars = re.sub(r':[0-9.]+', '', pars_tree.format('newick'))

with open("ground_truth_perzival.dnd") as f:
    true_tree = f.read()

print("--- Comparison true vs parsimony")
print(evaluate(true_newick=true_tree, recon_newick=tree_pars))

print("--- Comparison true vs NJ")
print(evaluate(true_newick=true_tree, recon_newick=tree_nj))

print("--- Comparison Froger")

apparatus = ListApparatus.from_csv(
    "../../benchmarks/scripts/tmp/parzival.csv",
    ref_witness="p9",
    delimiter=";")


stemma = UnorientedStemma(list_apparatus=apparatus,
                          min_frequency=1,
                          verbose=False
                          )
stemma.build_tree(anomaly_method="raw_cutoff", resolve_anomaly_threshold=5)
# pprint(stemma.children_of)
# print(stemma.to_newick())
stemma.display_tree()
print(stemma.get_contamination_loci())
print(evaluate(true_newick=true_tree, recon_newick=stemma.to_newick()))
